import pandas as pd
import unittest
from pathlib import Path
import constants
from constants import Tier

column_orig_names = [
    'FUNC1.gene', 'INFO.1.GENE_NAME', 'FUNC1.protein', 'FUNC1.coding', 
    'INFO.A.AF', 'INFO.1.FDP', 'INFO.A.FAO', 'FUNC1.transcript', 
    'FUNC1.function', 'FUNC1.oncomineGeneClass', 'FUNC1.oncomineVariantClass', 
    'FUNC1.location', 'rowtype', 'INFO...OID', 'FUNC1.CLNID1', 
    'FUNC1.CLNREVSTAT1', 'FUNC1.CLNSIG1', 'FUNC1.sift', 'FUNC1.polyphen', 
    'FUNC1.grantham', 'INFO.1.FAIL_REASON', 'CHROM', 'POS', 'INFO.1.END', 
    'FORMAT.1.CN', 'call', 'INFO...CI', 'ID', 'INFO...LEN', 'QUAL', 
    'INFO...CDF_MAPD', 'ALT', 'INFO...READ_COUNT', 'INFO.1.EXON_NUM', 
    'INFO.1.ANNOTATION', 'FILTER', 'Tier'
]
assert(len(constants.columns) == len(column_orig_names))

Col = constants.Col


def parse_oncomine_file(file: Path):
    df = pd.read_table(file, index_col='vcf.rownum', comment='#',
                       na_values=['.'], low_memory=False)
    df = df.reindex(columns=column_orig_names)
    df = df[[c for c in column_orig_names]]
    df.columns = constants.columns
    # df[Col.VAF] = df[Col.VAF].map('{:.1%}'.format)
    return df


def assign_default_tier(df: pd.DataFrame):
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
            if hotspot == 'Deleterious' or hotspot == 'Hotspot':
                tier = Tier.TIER_1_2

        return tier
    
    df[Col.TIER] = df.apply(populate_tier_default, axis=1)


def generate_intermediate_report(df: pd.DataFrame):
    
    snv_columns = [
        Col.GENE_NAME, Col.AA_CHANGE, Col.NUCLEOTIDE_CHANGE, Col.VAF, Col.TIER,
        Col.ROWTYPE, Col.COSM_ID, Col.TOTAL_DEPTH, Col.VARIANT_COUNT, 
        Col.REFSNP_ID, Col.CLINICAL_SIGNIFICANCE, Col.MUTATION_TYPE, 
        Col.GRANTHAM_SCORE, Col.LOCATION, Col.ONCOMINE_GENE_CLASS, Col.HOTSPOT,
        Col.POLYPHEN_SCORE, Col.SIFT_SCORE, Col.REFSEQ, Col.FAIL_REASON
    ]

    snv_nocall_columns = [
        Col.GENE_NAME, Col.AA_CHANGE, Col.REFSEQ, Col.NUCLEOTIDE_CHANGE,
        Col.VAF, Col.MUTATION_TYPE, Col.TOTAL_DEPTH, Col.VARIANT_COUNT,
        Col.ONCOMINE_GENE_CLASS, Col.HOTSPOT, Col.LOCATION, Col.ROWTYPE, 
        Col.COSM_ID, Col.REFSNP_ID, Col.REFSNP_STAT, Col.CLINICAL_SIGNIFICANCE,
        Col.FAIL_REASON
    ]

    cnv_columns = [
        Col.GENE_NAME, Col.COPY_NUMBER, Col.TIER, Col.ROWTYPE, Col.CALL, 
        Col.CHROMOSOME, Col.POSITION, Col.ID, Col.QUALITY, Col.FILTER, Col.MAPD, 
        Col.CI, Col.LENGTH, Col.END_POSITION, Col.ONCOMINE_GENE_CLASS,
        Col.HOTSPOT, Col.FAIL_REASON
    ]

    cnv_nocall_columns = [
        Col.GENE_NAME, Col.COPY_NUMBER, Col.CALL, Col.CI, Col.FILTER, 
        Col.ROWTYPE, Col.ID, Col.CHROMOSOME, Col.POSITION, Col.LENGTH,
        Col.END_POSITION, Col.ONCOMINE_GENE_CLASS, Col.HOTSPOT, Col.QUALITY,
        Col.MAPD, Col.FAIL_REASON
    ]

    fusion_columns = [
        Col.TIER, Col.FILTER, Col.CALL, Col.ROWTYPE, Col.ID, Col.CHROMOSOME,
        Col.GENE, Col.ALTERATION, Col.POSITION, Col.TOTAL_READ, Col.ANNOTATION,
        Col.EXON_NUMBER, Col.ONCOMINE_GENE_CLASS, Col.HOTSPOT, Col.FAIL_REASON
    ]

    fusion_nocall_columns = [
        Col.FILTER, Col.CALL, Col.ROWTYPE, Col.ID, Col.CHROMOSOME, Col.GENE, 
        Col.ALTERATION, Col.POSITION, Col.TOTAL_READ, Col.ANNOTATION, 
        Col.EXON_NUMBER, Col.ONCOMINE_GENE_CLASS, Col.FAIL_REASON
    ]

    snv_condition = f'`{Col.CALL}` == "POS"'\
        f' and `{Col.LOCATION}` not in ["intronic", "utr_3", "utr_5"]'\
        f' and `{Col.ROWTYPE}` in ["snp", "del", "ins", "complex", "mnp", "RNAExonTiles"]'\
        f' and `{Col.MUTATION_TYPE}`.notna()'\
        f' and `{Col.MUTATION_TYPE}` != "synonymous"'
    
    snv_nocall_condition = f'`{Col.CALL}` == "NOCALL"'\
        f' and `{Col.ROWTYPE}` not in ["CNV", "Fusion"]'
    
    cnv_condition = f'(`{Col.CALL}` in ["DEL", "AMP"]'\
        f' or `{Col.ROWTYPE}` == "LOH")'\
        f' and `{Col.GENE_NAME}`.notna()'
    
    cnv_nocall_condition = f'`{Col.CALL}` == "NOCALL"'\
        f' and `{Col.ROWTYPE}` in ["CNV", "LOH"]'
    
    fusion_condition = f'`{Col.CALL}` == "POS"'\
        f' and `{Col.ROWTYPE}` in ["Fusion", "RNAExonVariant"]'
    
    fusion_nocall_condition = f'{Col.FILTER} == "FAIL"'\
        f' and `{Col.ROWTYPE}` in ["Fusion", "RNAExonVariant"]'
    
    snv = df.query(snv_condition)[snv_columns]
    snv_nocall = df.query(snv_nocall_condition)[snv_nocall_columns]
    cnv = df.query(cnv_condition)[cnv_columns]
    cnv_nocall = df.query(cnv_nocall_condition)[cnv_nocall_columns]
    fusion = df.query(fusion_condition)[fusion_columns]
    fusion_nocall = df.query(fusion_nocall_condition)[fusion_nocall_columns]

    snv.loc[(snv[Col.GENE_NAME] == 'UGT1A1') 
            & (snv[Col.AA_CHANGE] == 'p.Gly71Arg'), Col.TIER] = Tier.TIER_3
    snv.loc[snv[Col.TOTAL_DEPTH] < 100, Col.TIER] = Tier.TIER_4
    snv = sort_by_tier_column(snv)

    cnv_tier_2_3_gene_names = [
            'AKT1', 'ALK', 'BRAF', 'CCND2', 'CCNE1', 'CD274', 'CDK4', 'CDK6', 
            'DDR1', 'DDR2', 'EGFR', 'EMSY', 'ERBB2', 'FGF23', 'FGF3', 'FGF4', 
            'FGF19', 'CCND1', 'FGF9', 'FGFR1', 'FGFR2', 'FGFR4', 'GNAS', 'KRAS',
            'MAP2K1', 'MCL1', 'MDM2', 'MET', 'MYC', 'NTRK1', 'PIK3CA', 'PTPN11',
            'RAF1', 'RICTOR', 'ROS1', 'SRC'
        ]
    cnv.loc[(cnv[Col.CALL] == 'AMP')
            & (cnv[Col.GENE_NAME].isin(cnv_tier_2_3_gene_names)),
            Col.TIER] = Tier.TIER_1_2
    
    fusion_tier_2_3_gene_names = (
            'ALK', 'BRAF', 'MET', 'ESR1', 'EGFR', 'ETV6', 'NTRK3', 'FLI1',
            'FGFR', 'FGFR3', 'NTRK2', 'NRG1', 'NTRK3', 'PAX8', 'RAF1', 'RELA',
            'RET', 'PIK3CA'
    )
    fusion.loc[fusion[Col.GENE].str.startswith(fusion_tier_2_3_gene_names),
               Col.TIER] = Tier.TIER_1_2
    cnv = sort_by_tier_column(cnv)
    
    return {
        'SNV': snv,
        'SNV_NOCALL': snv_nocall,
        'CNV': cnv,
        'CNV_NOCALL': cnv_nocall,
        'FUSION': fusion,
        'FUSION_NOCALL': fusion_nocall
    }


def write_dataframe_as_sheet(dataframes: dict, file: Path):
    with pd.ExcelWriter(file, engine='xlsxwriter') as writer:
        percent_format = writer.book.add_format({'num_format': '0.0%'})
        for key in dataframes:
            df = dataframes[key]
            df.to_excel(writer, sheet_name = key, index=False)
            worksheet = writer.sheets[key]
            if 'VAF' in df.columns:
                loc = df.columns.get_loc(Col.VAF)
                worksheet.set_column(loc, loc, None, percent_format)
            if key == 'CNV':
                df.reset_index(drop=True, inplace=True)
                size = len(df.index)
                loc = df.columns.get_loc(Col.CALL)
                worksheet.autofilter(0, loc, size, loc)
                for index, row in df.iterrows():
                    call = row[Col.CALL]
                    if call == 'AMP':
                        pass
                    else:
                        worksheet.set_row(index + 1, None, None, {'hidden': True})




def sort_by_tier_column(df: pd.DataFrame):
    df[Col.TIER] = pd.Categorical(df[Col.TIER], list(Tier), ordered=True)
    return df.sort_values(by=Col.TIER)


def sort_by_total_read_and_tier(df: pd.DataFrame):
    df[Col.TIER] = pd.Categorical(df[Col.TIER], list(Tier), ordered=True)
    return df.sort_values(by=[Col.TOTAL_READ, Col.TIER])


class ReadTests(unittest.TestCase):
    def test_read(self):
        file = Path('resources/M23-test-oncomine.tsv')
        print("Read dataframe.")
        df = parse_oncomine_file(file)
        self.assertIsNotNone(df)
        self.assertTrue('float' in str(df[Col.VAF].dtype).lower())
        assign_default_tier(df)


class IntermediateReportTests(unittest.TestCase):
    def test_report(self):
        source = Path('resources/M23-test-oncomine.tsv')
        dest = Path('test/M23-test.xlsx')
        df = parse_oncomine_file(source)
        assign_default_tier(df)
        reports = generate_intermediate_report(df)
        self.assertIsNotNone(reports)
        write_dataframe_as_sheet(reports, dest)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
