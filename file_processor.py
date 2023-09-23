import sys
import zipfile
from pathlib import Path

def get_case_name(file: Path):
    return file.stem.split('_')[0]


# 윈도우에서 인식하지 못하는 파일명 내 : 문자를 -로 변경
def unzip_to_destination_and_normalize(file: Path, dest_dir: Path):
    if dest_dir.exists():
        sys.exit('The specified directory %s already exists.' % dest_dir)
    
    dest_dir.mkdir(parents=True, exist_ok=True)

    zipdata = zipfile.ZipFile(file)
    zipinfos = zipdata.infolist()
    
    for zipinfo in zipinfos:
        zipinfo.filename = zipinfo.filename.replace(':', '-')
        zipdata.extract(zipinfo, dest_dir)
    
    print(f'Unzipped the file {file} to directory {dest_dir}')
    
    return dest_dir


def find_target_files(dir: Path):
    case_name = dir.name
    variants_dir = dir / 'Variants'
    case_dir = next(x for x in variants_dir.iterdir() 
                    if x.is_dir and x.name.startswith(case_name))
    targets = [x for x in case_dir.iterdir() 
                if x.is_file and x.name.startswith(case_name)]    
    oncomine_file = next(x for x in targets if x.name.endswith('-oncomine.tsv'))
    vcf_file = next(x for x in targets if x.suffix == '.vcf')
    qc_dir = dir / 'QC'
    qc_file = next(x for x in qc_dir.iterdir()
                   if x.suffix == ('.pdf') and x.name.startswith(case_name))
    
    return {
        'oncomine_file': oncomine_file,
        'vcf_file': vcf_file,
        'qc_file': qc_file
    }
