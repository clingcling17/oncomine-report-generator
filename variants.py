from abc import ABC, abstractmethod
import pandas as pd
from constants import Col, Tier
import re


def initialize_variant_blacklist(df: pd.DataFrame):
    Variant.blacklist = df


class Variant(ABC):
    blacklist = None

    def __init__(self, df: pd.DataFrame):
        self._generate_data(df)
        self._assign_tier()
        self._sort()
    
    @property
    @abstractmethod
    def call_condition(self) -> str:
        pass

    @property
    @abstractmethod
    def nocall_condition(self) -> str:
        pass

    @property
    @abstractmethod
    def columns(self) -> list[str]:
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
        Variant.sort_by_tier(self.call)
        

    def print_worksheet(self, writer: pd.ExcelWriter):
        self.call.to_excel(writer, sheet_name = self.__class__.__qualname__, index=False)
        self.call.to_excel(writer, sheet_name = f'{self.__class__.__qualname__}_nocall', index=False)


    @abstractmethod
    def generate_report_info(self):
        pass

    
    @staticmethod
    def sort_by_tier(df: pd.DataFrame):
        df[Col.TIER] = pd.Categorical(df[Col.TIER], list(Tier), ordered=True)
        df.sort_values(by=Col.TIER, inplace=True)
    

    @staticmethod
    def fill_na_tier(df, tier):
        df.loc[df[Col.TIER] == Tier.TIER_NA, Col.TIER] = tier


    @staticmethod
    def filter_out_blacklist(df):
        return df.drop(df[df[Col.TIER] == Tier.TIER_BLACKLIST].index)
    


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

    
    def _generate_data(self, df):
        super()._generate_data(df)
        self.call.loc[self.call[Col.AA_CHANGE].notna(),
                      Col.AA_CHANGE] = self.call[Col.AA_CHANGE].str.replace('Ter', '*')
    

    def _assign_tier(self):
        super()._assign_tier()
        self.call.loc[(self.call[Col.GENE_NAME] == 'UGT1A1')
            & (self.call[Col.AA_CHANGE] == 'p.Gly71Arg'), Col.TIER] = Tier.TIER_3

        self.call.loc[self.call[Col.TOTAL_DEPTH] < 500, Col.TIER] = Tier.TIER_4

        self.call.loc[(self.call[Col.GENE_NAME] == 'EGFR')
            & ((self.call[Col.AA_CHANGE].str.contains(r'p\.Glu746_.*del.*', regex=True))
                | (self.call[Col.AA_CHANGE].str.contains(r'p\.Leu747_.*del.*', regex=True))),
            Col.TIER] = Tier.TIER_1

        self.call.loc[(self.call[Col.GENE_NAME] == 'MAML3')
            & (self.call[Col.NUCLEOTIDE_CHANGE].str.contains(r'c\.1455_.*del.*', regex=True)
            & (self.call[Col.NUCLEOTIDE_CHANGE].str.len() >= 16)),
            Col.TIER] = Tier.TIER_BLACKLIST

        if Variant.blacklist is not None:
            blacklist = Variant.blacklist
            joined_df = self.call.merge(blacklist, how='left',
                                        indicator='indicator_info')
            joined_df.loc[joined_df['indicator_info'] == 'both', Col.TIER] = Tier.TIER_BLACKLIST
            joined_df.drop(columns='indicator_info', inplace=True)
            self.call = joined_df
    

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
        mut.loc[:, Col.VAF] = mut[Col.VAF].map('{:.1%}'.format) # pylint: disable=consider-using-f-string
        mut.columns = ['Gene', 'Amino acid change', 'Nucleotide change',
                    'Variant allele frequency(%)', 'Tier']
        return Variant.filter_out_blacklist(mut)
 


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

    
    def _generate_data(self, df):
        super()._generate_data(df)
        df = self.call
        df.drop(df[(df[Col.CALL] == 'AMP') & (df[Col.COPY_NUMBER] < 4)].index,
                inplace=True)


    def _assign_tier(self):
        tier_1_2_gene_names = [
            'AKT1', 'ALK', 'BRAF', 'CCND2', 'CCNE1', 'CD274', 'CDK4', 'CDK6', 
            'DDR1', 'DDR2', 'EGFR', 'EMSY', 'ERBB2', 'FGF23', 'FGF3', 'FGF4', 
            'FGF19', 'CCND1', 'FGF9', 'FGFR1', 'FGFR2', 'FGFR4', 'GNAS', 'KRAS',
            'MAP2K1', 'MCL1', 'MDM2', 'MET', 'MYC', 'NTRK1', 'PIK3CA', 'PTPN11',
            'RAF1', 'RICTOR', 'ROS1', 'SRC'
        ]
        super()._assign_tier()
        self.call.loc[(self.call[Col.CALL] == 'AMP')
                & (self.call[Col.GENE_NAME].isin(tier_1_2_gene_names)),
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
        amp = self.call.loc[self.call[Col.CALL] == 'AMP']
        amp = amp[[Col.GENE_NAME, Col.COPY_NUMBER, Col.TIER]]
        amp.columns = ['Gene', 'Estimated copy number', 'Tier']
        return amp
    


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
        if self.call.empty:
            return #빈 테이블일 때 .str 처리 에러
        tier_1_2_gene_names = (
            'ALK', 'BRAF', 'MET', 'ESR1', 'EGFR', 'ETV6', 'NTRK3', 'FLI1',
            'FGFR', 'FGFR3', 'NTRK2', 'NRG1', 'NTRK3', 'PAX8', 'RAF1', 'RELA',
            'RET', 'PIK3CA'
            )
        super()._assign_tier()
        self.call.loc[self.call[Col.GENE].str.startswith(tier_1_2_gene_names),
                      Col.TIER] = Tier.TIER_1_2
        
    
    def _sort(self):
        self.call[Col.TIER] = pd.Categorical(self.call[Col.TIER], list(Tier), ordered=True)
        self.call.sort_values(by=[Col.TOTAL_READ, Col.TIER], ascending=[False, True], inplace=True)
    

    def generate_report_info(self):
        if self.call.empty:
            return pd.DataFrame(columns=['GeneA', 'Chromosome:BreakpointA',
                                         'GeneB', 'Chromosome:BreakpointB',
                                         'Total Read', 'Tier'])
        fus = self.call.assign(ChBr = lambda x:
                        (x[Col.CHROMOSOME] + ':' + x[Col.POSITION].astype(str)))
        fus = fus[[Col.GENE, 'ChBr', Col.TOTAL_READ, Col.TIER]]

        fus = fus.groupby(Col.TOTAL_READ).agg({
            Col.GENE: list, 'ChBr': list, Col.TIER: 'min'})
        fus.sort_values(by=Col.TOTAL_READ, ascending=False, inplace=True)
        fus = fus.apply(pd.Series.explode, axis=1)
        fus.reset_index(inplace=True)
        assert len(fus.columns) == 6, f'column number is not 6, {len(fus.columns)}.'

        fus.columns = ['Total Read', 'GeneA', 'GeneB', 'ChBrA', 'ChBrB', 'Tier']
        fus = fus[['GeneA', 'ChBrA', 'GeneB', 'ChBrB', 'Total Read', 'Tier']]
        fus.rename(lambda x: x.replace('ChBr', 'Chromosome:Breakpoint'),
                axis='columns', inplace=True)
        Variant.sort_by_tier(fus)
        return fus
    