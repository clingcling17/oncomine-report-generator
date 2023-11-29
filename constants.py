from enum import Enum

class Col(object):
    GENE_NAME = 'Gene_name'
    GENE = 'Gene'
    AA_CHANGE = 'AA_Change'
    NUCLEOTIDE_CHANGE = 'Nucleotide_change'
    VAF = 'VAF'
    TOTAL_DEPTH = 'Total_depth'
    VARIANT_COUNT = 'Variant_count'
    REFSEQ = 'RefSeq'
    MUTATION_TYPE = 'mutation_type'
    ONCOMINE_GENE_CLASS = 'Oncomine_gene_class' #Function
    HOTSPOT = 'Hotspot'
    LOCATION = 'Location'
    ROWTYPE = 'rowtype'
    COSM_ID = 'COSM_ID'
    REFSNP_ID = 'RefSNP_id'
    REFSNP_STAT = 'RefSNP_stat'
    CLINICAL_SIGNIFICANCE = 'Clinical_significance'
    SIFT_SCORE = 'SIFT_score'
    POLYPHEN_SCORE = 'PolyPhen_score'
    GRANTHAM_SCORE = 'Grantham_score'
    FAIL_REASON = 'Fail_reason'
    CHROMOSOME = 'Chromosome'
    POSITION = 'Position'
    END_POSITION = 'End_position'
    COPY_NUMBER = 'Copy_number'
    CALL = 'Call'
    CI = 'CI'
    ID = 'ID'
    LENGTH = 'Length'
    QUALITY = 'Quality'
    MAPD = 'MAPD'
    ALTERATION = 'Alteration'
    TOTAL_READ = 'Total_Read'
    EXON_NUMBER = 'Exon_number'
    ANNOTATION = 'Annotation'
    FILTER = 'Filter'
    TIER = 'Tier'

columns = [value for name, value in vars(Col).items() if not name.startswith('_')]


class Metrics:
    MAPPED_READS = 'Mapped Reads'
    ON_TARGET = 'On Target'
    MEAN_DEPTH = 'Mean Depth'
    UNIFORMITY = 'Uniformity'

    MSI_SCORE='MSIScore'
    MSI_STATUS='MSIStatus'
    PERCENT_LOH = 'percentLOH'
    TMB_MUTATIONS_PER_MB='TMBMutationsPerMb'
    CELLULARITY = 'manually_input_percent_tumor_cellularity'
    MAPD = 'mapd'
    

class Tier(Enum):
    TIER_1_2 = 'I/II'
    TIER_NA = 'N/A'
    TIER_3 = 'III'
    TIER_3_4 = 'III/IV'
    TIER_4 = 'IV'

    @property
    def index(self):
        return list(Tier).index(self)

    def __str__(self):
        return self.value

    def __lt__(self, other):
        return self.index < other.index

assert Tier.TIER_NA < Tier.TIER_3
