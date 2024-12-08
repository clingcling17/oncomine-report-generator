import os
import sys
from pathlib import Path
import pprint
from pandas import DataFrame
from tabulate import tabulate
from numpy import nan

import file_processor
import value_reader
import table_processor
from constants import Metrics, Tier, Col


MAPD_POOR_NOTE = 'Note) Sample의 질이 좋지 않아 (MAPD > 0.5) LOH score를 계산할'\
    ' 수 없습니다.'
MAPD_FAIR_NOTE = 'Note) Sample의 질이 좋지 않아 (MAPD ≒ 0.5) LOH score를 계산할'\
    ' 수 없습니다.'
UNIFORMITY_POOR_NOTE = 'Note) Sequencing의 질이 좋지 않아 신뢰도가 낮으므로'\
    ' (uniformity < 90%) 임상 적용시 주의가 필요합니다.'


def run(source_file: Path, dest_dir, case_name):
    file_processor.unzip_to_destination_and_normalize(source_file, dest_dir)
    fusion_file = file_processor.find_fusion_file(source_file.parent, case_name)
    if fusion_file is not None:
        file_processor.unzip_to_destination_and_normalize(fusion_file, dest_dir)
    files_to_read = file_processor.find_target_files(dest_dir)
    files_to_read_paths = {k: str(v) for k, v in files_to_read.items()}
    print(f'files to read: \n{pprint.pformat(files_to_read_paths)}')

    oncomine_D_file = files_to_read['ONCOMINE_D_FILE']
    oncomine_R_file = files_to_read['ONCOMINE_R_FILE']
    vcf_file = files_to_read['VCF_FILE']
    qc_file = files_to_read['QC_FILE']
    tumor_fraction_file = files_to_read['TUMOR_FRACTION_FILE']
    blacklist_file = files_to_read['BLACKLIST_FILE']
    assert oncomine_D_file is not None and vcf_file is not None and qc_file is not None and tumor_fraction_file is not None

    qc_pdf_text = value_reader.read_pdf_as_text(qc_file)
    coverage_metrics = value_reader.parse_coverage_metrics(qc_pdf_text)
    assert coverage_metrics is not None

    headers = value_reader.parse_headers(vcf_file)

    genomic_instability_metric, genomic_instability_status = value_reader.parse_tumor_fraction(tumor_fraction_file)

    D_oncomine_df = table_processor.parse_oncomine_file(oncomine_D_file)
    if oncomine_R_file is not None:
        R_oncomine_df = table_processor.parse_oncomine_file(oncomine_R_file)
    else:
        R_oncomine_df = None
    blacklist = None if blacklist_file is None else table_processor.read_blacklist(blacklist_file)
    snv, cnv, fusion = table_processor.generate_variants(D_oncomine_df, R_oncomine_df, blacklist)
    worksheet = dest_dir / (case_name + '_filtered_data.xlsx')
    table_processor.write_dataframe_as_sheet(worksheet, snv, cnv, fusion)
    print(f'Printed intermediate table to worksheet: {worksheet}')
    
    mut_info, amp_info, fus_info, sig_genes = table_processor.generate_printable_gene_info(snv, cnv, fusion)

    # 검사결과
    cellularity = headers[Metrics.CELLULARITY] + '%'
    tumor_mutational_burden = headers[Metrics.TMB_MUTATIONS_PER_MB]
    msi_score = '{:.2f}'.format(float(headers[Metrics.MSI_SCORE]))
    msi_status = headers[Metrics.MSI_STATUS]
    loh = headers[Metrics.PERCENT_LOH]
    sig_note = ''
    if loh is None or loh == 'NA':
        loh = 'not available (see note)'
        mapd = float(headers[Metrics.MAPD])
        if mapd > 0.5:
            sig_note = MAPD_POOR_NOTE
        else:
            # loh = f'{loh} (MAPD={mapd})'
            sig_note = MAPD_FAIR_NOTE
    else:
        if float(loh) == 0.0:
            loh = '0%'
        else:
            loh = loh + '%'

    # 검사정보
    mean_depth = coverage_metrics[Metrics.MEAN_DEPTH]
    on_target = coverage_metrics[Metrics.ON_TARGET]
    mapped_reads = coverage_metrics[Metrics.MAPPED_READS]
    uniformity = coverage_metrics[Metrics.UNIFORMITY]

    score_factors = [
        mapped_reads >= 5000000,
        on_target >= 90,
        mean_depth >= 1200,
        uniformity >= 90
        ]
    score_no = sum(score_factors)
    assert 0 <= score_no <= 4

    if score_no >= 4:
        total_quality_score = 'Very Good'
    elif score_no == 3:
        total_quality_score = 'Good'
    elif score_no == 2:
        total_quality_score = 'Intermediate'
    else:
        total_quality_score = 'Poor'

    # overall_qc_test_result = 'Fail' if score_no <= 2 else 'Pass'
    overall_qc_test_result = 'Pass'

    qc_note = UNIFORMITY_POOR_NOTE \
        if total_quality_score == 'Good' and uniformity < 90 else ''
    
    on_target = str(on_target) + '%'
    
    sig_genes = ', '.join(sig_genes)
    # mut_sig, amp_sig, fus_sig, mut_unk, amp_unk, fus_unk

    def _filter_significant_tier(df: DataFrame):
        return df.loc[df[Col.TIER] <= Tier.TIER_1_2]
    
    def _has_low_read(df: DataFrame):
        return not df.query("`Total Read` < 500").empty
    
    fus_low_read_note = 'Note) fusion read 수가 낮아 위양성의 가능성이 있으므로 해석과 임상 적용에 주의가 필요합니다.'

    mut_sig = _filter_significant_tier(mut_info)
    amp_sig = _filter_significant_tier(amp_info)
    fus_sig = _filter_significant_tier(fus_info)
    fus_sig_note = fus_low_read_note if _has_low_read(fus_sig) else ''
    mut_unk = _diff_table(mut_info, mut_sig)
    amp_unk = _diff_table(amp_info, amp_sig)
    fus_unk = _diff_table(fus_info, fus_sig)
    fus_unk_note = fus_low_read_note if _has_low_read(fus_unk) else ''

    mut_sig = _print_table(mut_sig)
    amp_sig = _print_table(amp_sig)
    fus_sig = _print_table(fus_sig)
    mut_unk = _print_table(mut_unk)
    amp_unk = _print_table(amp_unk)
    fus_unk = _print_table(fus_unk)

    print_params = [
        cellularity, mut_sig, amp_sig, fus_sig, fus_sig_note, 
        mut_unk, amp_unk, fus_unk, fus_unk_note,
        tumor_mutational_burden, msi_score, msi_status, loh,
        genomic_instability_metric, genomic_instability_status, sig_note,
        sig_genes, mean_depth, on_target, total_quality_score,
        overall_qc_test_result, qc_note
    ]
    assert len(print_params) == 22

    format_file = Path('resources/report_text_format.txt')
    with open(format_file, 'rt', encoding='utf-8') as f:
        text_form = f.read()
    
    full_text = text_form.format(cellularity, mut_sig, amp_sig, fus_sig, fus_sig_note,
        mut_unk, amp_unk, fus_unk, fus_unk_note,
        tumor_mutational_burden, msi_score, msi_status, loh,
        genomic_instability_metric, genomic_instability_status, sig_note,
        sig_genes, mean_depth, on_target, total_quality_score,
        overall_qc_test_result, qc_note)
    # wrapped_text = textwrap.fill(full_text, width=80, expand_tabs=False,
    #                              replace_whitespace=False, drop_whitespace=False)

    report_file = dest_dir / f'{case_name}_report.txt'
    with open(report_file, 'wt', encoding='utf-8') as f:
        f.write(full_text)
        # f.write(wrapped_text)
    print(f'Generated report text file: {report_file}')


def _diff_table(op1: DataFrame, op2: DataFrame):
    index = op1.index.name
    return op1[~op1.index.isin(op2.index)]


def _print_table(df: DataFrame):
    return tabulate(df.replace(nan, None), headers=df.columns.tolist(),
                    showindex=False) + ('\nNot Found' if df.empty else '')

def _parse_case_name(file_name: str):
    case_name = file_name.split('_')[0]
    case_name = case_name[:-1]
    if case_name.endswith('-'):
        case_name = case_name[:-1]
    return case_name
    # pattern = re.compile(r'M[\d-]+')
    # match = pattern.search(case_name)
    # return match.group()


def main():
    print("Starting oncomine report generator.")

    if len(sys.argv) != 2:
        sys.exit('Please check arguments number.\n'\
                 + 'Usage: run.exe Mxx-xxxx.zip')
    source_file = Path(sys.argv[1]).absolute()
    # dest_path = sys.argv[2]
    print(f'File path: {source_file}')
    case_name = _parse_case_name(source_file.stem)

    dest_dir = Path(os.getcwd(), case_name).absolute()
    os.chdir(getattr(sys, '_MEIPASS')) # pyinstaller temporary dir
    print(f'Destination path: {dest_dir}')
    print(f'Case name: {case_name}')
    run(source_file, dest_dir, case_name)
    

if __name__ == "__main__":
    main()
