import re
import unittest
from pathlib import Path
import pprint
from ast import literal_eval
import fitz
from constants import Metrics

def read_pdf_as_text(file: Path):
    with fitz.open(file) as doc:
        text = chr(12).join([page.get_text() for page in doc])
    return text


def parse_coverage_metrics(text: str):
    pattern = re.compile(r'\bCoverage metrics\b[\s\S]+\b'\
                        r'Sample Name\s+'\
                        r'BarCode\s+'\
                        r'Mapped Reads\s+'\
                        r'On Target\s+'\
                        r'Mean Depth\s+'\
                        r'Uniformity\s+'\
                        r'\S+\s+'\
                        r'\S+\s+'\
                        r'(\d+)\s+'
                        r'(\d+.\d+%)\s+'\
                        r'([0-9]+(?:\.[0-9]+)?)\s+'\
                        r'(\d+.\d+%)\s+')
    match = pattern.search(text)
    matched = match.groups()
    mapped_reads = int(matched[0].replace(",", ""))
    on_target = float(matched[1].strip("%"))
    mean_depth = literal_eval(matched[2].replace(",", ""))
    uniformity = float(matched[3].strip("%"))

    coverage_metrics = {
        Metrics.MAPPED_READS: mapped_reads,
        Metrics.ON_TARGET: on_target,
        Metrics.MEAN_DEPTH: mean_depth,
        Metrics.UNIFORMITY: uniformity
    }

    print(f'Coverage metrics: \n{pprint.pformat(coverage_metrics)}')

    return coverage_metrics


def parse_headers(file: Path):
    header_parameters = [
        Metrics.MSI_SCORE, Metrics.MSI_STATUS, Metrics.PERCENT_LOH,
        Metrics.TMB_MUTATIONS_PER_MB, Metrics.TOTAL_MAPPED_FUSION_PANEL_READS,
        Metrics.CELLULARITY, Metrics.MAPD
        ]
    
    result = {
        Metrics.PERCENT_LOH: None
    }
    
    with open(file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        for line in lines:
            if line == '':
                continue
            if not line.startswith('##'):
                break
            match = re.match(r'##(\S+)=(\S+)', line)
            if (match.group(1) in header_parameters):
                result[match.group(1)] = match.group(2)
            if len(result) == len(header_parameters):
                break
    
    assert all(item in result.keys() for item in header_parameters)

    return result

    
class ReadTests(unittest.TestCase):
    def test_read(self):
        file = Path('resources/M23-test_QC.pdf')
        text = read_pdf_as_text(file)
        self.assertIsNotNone(text)
        print('\nRead pdf.')
        metrics = parse_coverage_metrics(text)
        self.assertIsNotNone(metrics)
        print(metrics)

        file = Path('resources/M23-test.vcf')
        headers = parse_headers(file)
        self.assertIsNotNone(headers)
        print(headers)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
