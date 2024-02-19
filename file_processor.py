import sys
from zipfile import ZipFile
from pathlib import Path


def find_fusion_file(dir: Path, case_name: Path):
    return next((x for x in dir.iterdir()
                if x.suffix == ('.zip') and
                x.name.startswith((case_name + 'R', case_name + '-R'))), None)


# 윈도우에서 인식하지 못하는 파일명 내 : 문자를 -로 변경
def unzip_to_destination_and_normalize(source_file: Path, dest_dir: Path):
    if not dest_dir.exists():
        dest_dir.mkdir(parents=True, exist_ok=True)

    with ZipFile(source_file, 'r') as zipdata:
        zipinfos = zipdata.infolist()
        for zipinfo in zipinfos:
            zipinfo.filename = zipinfo.filename.replace(':', '-')
            if not Path(dest_dir, zipinfo.filename).exists():
                zipdata.extract(zipinfo, dest_dir)
    
    print(f'Unzipped the file [{source_file}] to directory [{dest_dir}].')



def find_target_files(root: Path):
    case_name = root.name
    variants_dir = root / 'Variants'
    D_variants_case_dir = next(x for x in variants_dir.iterdir()
                    if x.is_dir and
                    x.name.startswith((case_name + 'D', case_name + '-D')))
    targets = [x for x in D_variants_case_dir.iterdir()
                if x.is_file and x.name.startswith(case_name)]
    D_oncomine_file = next(x for x in targets if x.name.endswith('-oncomine.tsv'))

    vcf_file = next(x for x in targets if x.suffix == '.vcf')
    qc_dir = root / 'QC'
    qc_file = next(x for x in qc_dir.iterdir()
                   if x.suffix == ('.pdf') and x.name.startswith(case_name))
    
    tumor_fraction_file = root / 'CnvActor' / 'TumorFraction' / 'tumor_fraction.json'

    R_variants_case_dir = next((x for x in variants_dir.iterdir()
                    if x.is_dir and
                    x.name.startswith((case_name + 'R', case_name + '-R'))), None)
    if R_variants_case_dir is not None:
        targets = [x for x in R_variants_case_dir.iterdir()
                if x.is_file and x.name.startswith(case_name)]
        R_oncomine_file = next(x for x in targets if x.name.endswith('-oncomine.tsv'))
    else:
        R_oncomine_file = None

    return {
        'ONCOMINE_D_FILE': D_oncomine_file,
        'ONCOMINE_R_FILE': R_oncomine_file,
        'VCF_FILE': vcf_file,
        'QC_FILE': qc_file,
        'TUMOR_FRACTION_FILE': tumor_fraction_file
    }
