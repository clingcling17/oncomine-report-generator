from pathlib import Path
import unittest
import pandas as pd
import constants
import variants
from variants import Variant

Col = constants.Col


def parse_oncomine_file(file: Path):
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

    df = pd.read_table(file, index_col='vcf.rownum', comment='#',
                       na_values=['.'], low_memory=False)
    not_exist_columns = [x for x in column_orig_names if x not in df.columns.tolist()]
    not_exist_columns.remove('Tier')
    print(f'Columns not in current tsv file: {not_exist_columns}')
    df = df.reindex(columns=column_orig_names) #tsv 파일에 존재하지 않는 칼럼이 있을 경우 추가
    df = df[[c for c in column_orig_names]]
    df.columns = constants.columns
    df[Col.TIER] = df[Col.TIER].apply(str)
    return df


def generate_variants(df: pd.DataFrame):
    snv = variants.SNV(df)
    cnv = variants.CNV(df)
    fusion = variants.Fusion(df)

    return snv, cnv, fusion


def write_dataframe_as_sheet(file, snv: Variant, cnv: Variant, fusion: Variant):
    with pd.ExcelWriter(file, engine='xlsxwriter') as writer: # pylint: disable=abstract-class-instantiated
        snv.print_worksheet(writer)
        cnv.print_worksheet(writer)
        fusion.print_worksheet(writer)


def filter_significant_tier(df: pd.DataFrame):
    return df.loc[df[Col.TIER] == constants.Tier.TIER_1_2]
  

def generate_printable_gene_info(snv: Variant, cnv: Variant, fusion: Variant):
    mut = snv.generate_report_info()
    amp = cnv.generate_report_info()
    fus = fusion.generate_report_info()
    Variant.fill_na_tier(mut, constants.Tier.TIER_3)
    Variant.fill_na_tier(amp, constants.Tier.TIER_3)
    Variant.fill_na_tier(fus, constants.Tier.TIER_3)

    mut_sig_genes = filter_significant_tier(mut)['Gene'].tolist()
    amp_sig_genes = filter_significant_tier(amp)['Gene'].tolist()
    if fus.empty:
        fus_sig_genes = []
    else:
        fus_sig_genes = filter_significant_tier(fus).apply(
                lambda x: x['GeneA'] + '-' + x['GeneB'] + ' fusion', axis=1).tolist()
    sig_genes = list(set().union(mut_sig_genes, amp_sig_genes)) + list(set(fus_sig_genes)) #중복 제거

    return mut, amp, fus, sig_genes

# class ReadTests(unittest.TestCase):
#     def test_read(self):
#         file = Path('test/resources/M23-test-oncomine.tsv')
#         print("Read dataframe.")
#         df = parse_oncomine_file(file)
#         self.assertIsNotNone(df)
#         self.assertTrue('float' in str(df[Col.VAF].dtype).lower())
#         assign_default_tier(df)


# class IntermediateReportTests(unittest.TestCase):
#     def test_report(self):
#         source = Path('test/resources/M23-test-oncomine.tsv')
#         # dest = Path('test/M23-test.xlsx')
#         df = parse_oncomine_file(source)
#         assign_default_tier(df)
#         reports = generate_intermediate_tables(df)
#         self.assertIsNotNone(reports)
#         # write_dataframe_as_sheet(dest, **reports)
#         generate_printable_gene_info(reports['SNV'], reports['CNV'], reports['FUSION'])


# if __name__ == '__main__':
#     runner = unittest.TextTestRunner(verbosity=2)
#     unittest.main(testRunner=runner)
