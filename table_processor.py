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


def generate_printable_gene_info(snv: Variant, cnv: Variant, fusion: Variant):
    mut_info, mut_sig_genes = snv.generate_report_info()
    amp_info, amp_sig_genes = cnv.generate_report_info()
    fus_info, fus_sig_genes = fusion.generate_report_info()
    sig_genes = list(set().union(mut_sig_genes, amp_sig_genes)) + fus_sig_genes
    return mut_info, amp_info, fus_info, sig_genes

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
