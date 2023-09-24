from abc import ABC, abstractmethod
import pandas as pd
from constants import Col, Tier


class Variant(ABC):
    def __init__(self, df: pd.DataFrame):
        self._generate_data(df)
        self._assign_tier()
        self._sort()
    
    @property
    @abstractmethod
    def call_condition(self):
        pass

    @property
    @abstractmethod
    def nocall_condition(self):
        pass

    @property
    @abstractmethod
    def columns(self):
        pass


    def _generate_data(self, df):
        self.call = df.query(self.call_condition)[self.columns]
        nocall_columns = [x for x in self.columns if x != Col.TIER]
        self.nocall = df.query(self.nocall_condition)[nocall_columns]


    def _assign_tier(self):
        def populate_tier_default(row):
            clinical_significance = row[Col.CLINICAL_SIGNIFICANCE]
            hotspot = row[Col.HOTSPOT]

            tier = Tier.TIER_NA
            if pd.notna(clinical_significance):
                if 'benign' in clinical_significance.lower():
                    tier = Tier.TIER_4
                elif clinical_significance == 'not_provided'\
                    or clinical_significance == 'Uncertain_significance'\
                    or 'conflicting' in clinical_significance.lower():
                    tier = Tier.TIER_3_4
            if pd.notna(hotspot):
                if hotspot in ['Deleterious', 'Hotspot']:
                    tier = Tier.TIER_1_2

            return tier
        
        self.call[Col.TIER] = self.call.apply(populate_tier_default, axis=1)

    
    def _sort(self):
        Variant.sort_by_tier_column(self.call)
        

    def print_worksheet(self, writer: pd.ExcelWriter):
        self.call.to_excel(writer, sheet_name = self.__class__.__qualname__, index=False)
        self.call.to_excel(writer, sheet_name = f'{self.__class__.__qualname__}_nocall', index=False)


    @abstractmethod
    def generate_report_info(self):
        pass


    @staticmethod
    def filter_significant_tier(df: pd.DataFrame):
        return df.loc[df[Col.TIER] == Tier.TIER_1_2]
    

    @staticmethod
    def sort_by_tier_column(df: pd.DataFrame):
        df[Col.TIER] = pd.Categorical(df[Col.TIER], list(Tier), ordered=True)
        df.sort_values(by=Col.TIER, inplace=True)
    

    @staticmethod
    def fill_na_tier(*dfs):
        for df in dfs:
            df2 = df.loc[df[Col.TIER] == Tier.TIER_NA, Col.TIER] = Tier.TIER_3
            # TIER_NA가 NaN 처리되는 것으로 보임
            # df.loc[df[Col.TIER].isna(), Col.TIER] = Tier.TIER_3
    


class SNV(Variant):
    columns = [
        Col.GENE_NAME, Col.AA_CHANGE, Col.NUCLEOTIDE_CHANGE, Col.VAF, Col.TIER,
        Col.ROWTYPE, Col.COSM_ID, Col.TOTAL_DEPTH, Col.VARIANT_COUNT,
        Col.CLINICAL_SIGNIFICANCE, Col.MUTATION_TYPE, Col.LOCATION,
        Col.ONCOMINE_GENE_CLASS, Col.HOTSPOT, Col.POLYPHEN_SCORE, Col.SIFT_SCORE,
        Col.REFSEQ, Col.REFSNP_ID, Col.REFSNP_STAT, Col.FAIL_REASON
    ]

    call_condition = f'`{Col.CALL}` == "POS"'\
        f' and `{Col.LOCATION}` not in ["intronic", "utr_3", "utr_5"]'\
        f' and `{Col.ROWTYPE}` in ["snp", "del", "ins", "complex", "mnp", "RNAExonTiles"]'\
        f' and `{Col.MUTATION_TYPE}`.notna()'\
        f' and `{Col.MUTATION_TYPE}` != "synonymous"'
    
    nocall_condition = f'`{Col.CALL}` == "NOCALL"'\
        f' and `{Col.ROWTYPE}` not in ["CNV", "Fusion"]'
    

    def __init__(self, df: pd.DataFrame):
        super().__init__(df)
    

    def _assign_tier(self):
        super()._assign_tier()
        self.call.loc[(self.call[Col.GENE_NAME] == 'UGT1A1')
            & (self.call[Col.AA_CHANGE] == 'p.Gly71Arg'), Col.TIER] = Tier.TIER_3
        self.call.loc[self.call[Col.TOTAL_DEPTH] < 100, Col.TIER] = Tier.TIER_4
    

    def print_worksheet(self, writer: pd.ExcelWriter):
        super().print_worksheet(writer)
        call_sheet = writer.sheets[self.__class__.__qualname__]
        nocall_sheet = writer.sheets[f'{self.__class__.__qualname__}_nocall']

        vaf_loc = self.call.columns.get_loc(Col.VAF)
        percent_format = writer.book.add_format({'num_format': '0.0%'})
        call_sheet.set_column(vaf_loc, vaf_loc, None, percent_format)
        nocall_sheet.set_column(vaf_loc, vaf_loc, None, percent_format)
    
    
    def generate_report_info(self):
        mut = self.call[[Col.GENE_NAME, Col.AA_CHANGE, Col.NUCLEOTIDE_CHANGE,
                          Col.VAF, Col.TIER]]
        mut_sig_genes = Variant.filter_significant_tier(mut)[Col.GENE_NAME].tolist()
        mut.loc[:, Col.VAF] = mut[Col.VAF].map('{:.1%}'.format) # pylint: disable=consider-using-f-string
        mut.columns = ['Gene', 'Amino acid change', 'Nucleotide change',
                    'Variant allele frequency(%)', 'Tier']
        Variant.fill_na_tier(mut)
        return mut, mut_sig_genes
 


class CNV(Variant):
    columns = [
        Col.GENE_NAME, Col.COPY_NUMBER, Col.TIER, Col.ROWTYPE, Col.CALL,
        Col.CHROMOSOME, Col.POSITION, Col.ID, Col.QUALITY, Col.FILTER, Col.MAPD,
        Col.CI, Col.LENGTH, Col.END_POSITION, Col.ONCOMINE_GENE_CLASS,
        Col.HOTSPOT, Col.FAIL_REASON, Col.CLINICAL_SIGNIFICANCE
    ]

    call_condition = f'(`{Col.CALL}` in ["DEL", "AMP"]'\
        f' or `{Col.ROWTYPE}` == "LOH")'\
        f' and `{Col.GENE_NAME}`.notna()'
    
    nocall_condition = f'`{Col.CALL}` == "NOCALL"'\
        f' and `{Col.ROWTYPE}` in ["CNV", "LOH"]'
    

    def __init__(self, df: pd.DataFrame):
        super().__init__(df)


    def _assign_tier(self):
        tier_2_3_gene_names = [
            'AKT1', 'ALK', 'BRAF', 'CCND2', 'CCNE1', 'CD274', 'CDK4', 'CDK6', 
            'DDR1', 'DDR2', 'EGFR', 'EMSY', 'ERBB2', 'FGF23', 'FGF3', 'FGF4', 
            'FGF19', 'CCND1', 'FGF9', 'FGFR1', 'FGFR2', 'FGFR4', 'GNAS', 'KRAS',
            'MAP2K1', 'MCL1', 'MDM2', 'MET', 'MYC', 'NTRK1', 'PIK3CA', 'PTPN11',
            'RAF1', 'RICTOR', 'ROS1', 'SRC'
        ]
        super()._assign_tier()
        self.call.loc[(self.call[Col.CALL] == 'AMP')
                & (self.call[Col.GENE_NAME].isin(tier_2_3_gene_names)),
                Col.TIER] = Tier.TIER_1_2


    def print_worksheet(self, writer: pd.ExcelWriter):
        super().print_worksheet(writer)
        worksheet = writer.sheets[self.__class__.__qualname__]
        df = self.call.reset_index(drop=True)
        size = len(df.index)
        call_loc = df.columns.get_loc(Col.CALL)
        worksheet.autofilter(0, call_loc, size, call_loc)
        for index, row in df.iterrows():
            call = row[Col.CALL]
            if call == 'AMP':
                pass
            else:
                worksheet.set_row(index + 1, None, None, {'hidden': True})
    

    def generate_report_info(self):
        amp = self.call[[Col.GENE_NAME, Col.COPY_NUMBER, Col.TIER]]
        amp_sig_genes = Variant.filter_significant_tier(amp)[Col.GENE_NAME].tolist()
        amp.columns = ['Gene', 'Estimated copy number', 'Tier']
        Variant.fill_na_tier(amp)
        return amp, amp_sig_genes
    


class Fusion(Variant):
    columns = [
        Col.TIER, Col.FILTER, Col.CALL, Col.ROWTYPE, Col.ID, Col.CHROMOSOME,
        Col.GENE, Col.ALTERATION, Col.POSITION, Col.TOTAL_READ, Col.ANNOTATION,
        Col.EXON_NUMBER, Col.ONCOMINE_GENE_CLASS, Col.HOTSPOT, Col.FAIL_REASON,
        Col.CLINICAL_SIGNIFICANCE
    ]
    
    call_condition = f'`{Col.CALL}` == "POS"'\
        f' and `{Col.ROWTYPE}` in ["Fusion", "RNAExonVariant"]'
    
    nocall_condition = f'{Col.FILTER} == "FAIL"'\
        f' and `{Col.ROWTYPE}` in ["Fusion", "RNAExonVariant"]'
    

    def __init__(self, df: pd.DataFrame):
        super().__init__(df)


    def _assign_tier(self):
        tier_2_3_gene_names = (
            'ALK', 'BRAF', 'MET', 'ESR1', 'EGFR', 'ETV6', 'NTRK3', 'FLI1',
            'FGFR', 'FGFR3', 'NTRK2', 'NRG1', 'NTRK3', 'PAX8', 'RAF1', 'RELA',
            'RET', 'PIK3CA'
            )
        super()._assign_tier()
        self.call.loc[self.call[Col.GENE].str.startswith(tier_2_3_gene_names),
                      Col.TIER] = Tier.TIER_1_2
        
    
    def _sort(self):
        self.call[Col.TIER] = pd.Categorical(self.call[Col.TIER], list(Tier), ordered=True)
        self.call.sort_values(by=[Col.TOTAL_READ, Col.TIER], ascending=[False, True], inplace=True)
    

    def generate_report_info(self):
        fus = self.call.assign(ChBr = lambda x:
                        (x[Col.CHROMOSOME] + ':' + x[Col.POSITION].astype(str)))
        fus = fus[[Col.TOTAL_READ, Col.GENE, 'ChBr', Col.TIER]]

        fus = fus.groupby(Col.TOTAL_READ).agg({
            Col.GENE: list, 'ChBr': list, Col.TIER: 'min'})
        fus.sort_values(by=Col.TOTAL_READ, ascending=False, inplace=True)
        fus = fus.apply(pd.Series.explode, axis=1)
        fus.reset_index(inplace=True)
        assert len(fus.columns) == 6

        fus.columns = ['Total Read', 'GeneA', 'GeneB', 'ChBrA', 'ChBrB', 'Tier']
        fus = fus[['GeneA', 'ChBrA', 'GeneB', 'ChBrB', 'Total Read', 'Tier']]
        fus.rename(lambda x: x.replace('ChBr', 'Chromosome:Breakpoint'),
                axis='columns', inplace=True)
        
        fus_sig_genes = Variant.filter_significant_tier(fus).apply(
            lambda x: x['GeneA'] + '-' + x['GeneB'] + ' fusion', axis=1).tolist()
        
        Variant.sort_by_tier_column(fus)
        Variant.fill_na_tier(fus)
        return fus, fus_sig_genes
    