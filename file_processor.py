import sys
from zipfile import ZipFile
from pathlib import Path


# 윈도우에서 인식하지 못하는 파일명 내 : 문자를 -로 변경
def unzip_to_destination_and_normalize(source_file: Path, dest_dir: Path):
    if dest_dir.exists():
        sys.exit(f'The specified directory {dest_dir} already exists.')

    dest_dir.mkdir(parents=True, exist_ok=True)

    with ZipFile(source_file, 'r') as zipdata:
        zipinfos = zipdata.infolist()
        for zipinfo in zipinfos:
            zipinfo.filename = zipinfo.filename.replace(':', '-')
            zipdata.extract(zipinfo, dest_dir)
    
    print(f'Unzipped the file [{source_file}] to directory [{dest_dir}].')


def find_target_files(root: Path):
    case_name = root.name
    variants_dir = root / 'Variants'
    variants_case_dir = next(x for x in variants_dir.iterdir()
                    if x.is_dir and x.name.startswith(case_name))
    targets = [x for x in variants_case_dir.iterdir()
                if x.is_file and x.name.startswith(case_name)]
    oncomine_file = next(x for x in targets if x.name.endswith('-oncomine.tsv'))
    vcf_file = next(x for x in targets if x.suffix == '.vcf')
    qc_dir = root / 'QC'
    qc_file = next(x for x in qc_dir.iterdir()
                   if x.suffix == ('.pdf') and x.name.startswith(case_name))
    
    return {
        'ONCOMINE_FILE': oncomine_file,
        'VCF_FILE': vcf_file,
        'QC_FILE': qc_file
    }
