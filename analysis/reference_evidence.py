"""Stage 3C: provenance-tracked published reference *examples*.

No sequence alignment, annotation transfer, phage-host prediction, or
input sequence identification is carried out in this module.
"""
from dataclasses import dataclass, asdict

@dataclass(frozen=True)
class EvidenceRecord:
    reference_id: str
    organism_or_phage: str
    feature: str
    accession: str
    evidence_class: str
    finding: str
    interpretation_limit: str
    source: str
    source_url: str


EVIDENCE = (
    EvidenceRecord(
        'M102AD-GENOME', 'S. mutans phage M102AD', 'Complete phage genome',
        'DQ386162; RefSeq NC_028984.1', 'DATABASE_SEQUENCE',
        'Genome is 30,664 bp; 40 ORFs were reported in the 2012 study.',
        'Does not establish the function of a submitted sequence or its host range.',
        'Delisle et al. 2012; GenBank/RefSeq', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3302630/'),
    EvidenceRecord(
        'M102AD-ORF14', 'S. mutans phage M102AD', 'ORF14: putative receptor-binding protein',
        'Genome DQ386162; ORF14 nt 14973–17129 (paper Table 3)',
        'PUBLISHED_COMPUTATIONAL_ANNOTATION',
        'The authors assigned ORF14 a putative RBP function using comparisons to related phages.',
        'ORF14 binding to a specific purified bacterial molecule was not demonstrated; this is not an automatic match to uploaded ORFs.',
        'Delisle et al. 2012, Table 3 and discussion', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3302630/'),
    EvidenceRecord(
        'M102AD-VIRION', 'S. mutans phage M102AD', 'Structural virion proteins',
        'Genome DQ386162', 'EXPERIMENTAL_PROTEOMICS',
        'Authors identified virion structural proteins by mass spectrometry.',
        'Detection of proteins within virions does not experimentally establish ORF14 receptor-binding activity.',
        'Delisle et al. 2012', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3302630/'),
    EvidenceRecord(
        'M102AD-ADSORPTION', 'S. mutans phage M102AD', 'Serotype-dependent adsorption',
        'Genome DQ386162', 'EXPERIMENTAL_HOST_PHENOTYPE',
        'Adsorption and antiserum assays implicated serotype-specific surface carbohydrates; UA159 was not susceptible in their tested host-range panel.',
        'The exact individual receptor molecule and the responsible phage protein were not experimentally confirmed.',
        'Delisle et al. 2012, Results/Discussion', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3302630/'),
    EvidenceRecord(
        'SMHBZ8-GENOME', 'S. mutans phage SMHBZ8', 'Genome accession',
        'MT430910', 'DATABASE_SEQUENCE',
        'SMHBZ8 genome accession reported in the 2021 experimental study.',
        'An accession does not demonstrate any particular ORF function or protein-receptor contact.',
        'Wolfoviz-Zilberman et al. 2021', 'https://doi.org/10.3390/antibiotics10081015'),
    EvidenceRecord(
        'SMHBZ8-PHENOTYPE', 'S. mutans phage SMHBZ8', 'Phage activity and caries outcomes',
        'MT430910', 'EXPERIMENTAL_HOST_PHENOTYPE',
        'The study measured reduced bacterial burden and caries-preventive outcomes in its in vitro and mouse models.',
        'Does not identify a phage receptor or establish activity on every S. mutans strain.',
        'Wolfoviz-Zilberman et al. 2021, Results', 'https://doi.org/10.3390/antibiotics10081015'),
    EvidenceRecord(
        'SMHBZ8-LYSINS', 'S. mutans phage SMHBZ8', 'Reported lysin annotations',
        'QKE60409.1; QKE60410.1', 'PUBLISHED_COMPUTATIONAL_ANNOTATION',
        'Two lysins were identified as annotated candidate lysins in the referenced SMHBZ8 genome comparison.',
        'Not evidence of enzyme activity in the uploaded sample or a demonstrated infection-specific lysis pathway.',
        'Ben-Zaken et al. 2021 (cited by Wolfoviz-Zilberman et al.); Viruses 13:825', 'https://doi.org/10.3390/v13050825'),
    EvidenceRecord(
        'UA159-GENOME', 'S. mutans UA159', 'Reference bacterial genome',
        'NC_004350.2', 'DATABASE_SEQUENCE',
        'Public reference chromosome for S. mutans UA159.',
        'A bacterial FASTA upload is not verified as UA159 by this application.',
        'NCBI RefSeq; Ajdic et al. 2002', 'https://www.ncbi.nlm.nih.gov/nuccore/NC_004350.2'),
    EvidenceRecord(
        'UA159-SPAP', 'S. mutans UA159', 'SpaP / antigen I/II cell-surface adhesin',
        'AAN58348.1; gene spaP', 'DATABASE_PROTEIN_ANNOTATION',
        'NCBI annotates a UA159 surface antigen SpaP; research describes adhesion-related function.',
        'SpaP is not established here as the receptor for M102AD or SMHBZ8.',
        'NCBI Protein; Ajdic et al. 2002', 'https://www.ncbi.nlm.nih.gov/protein/AAN58348.1'),
    EvidenceRecord(
        'UA159-RGP', 'S. mutans UA159', 'Rhamnose-glucose cell-wall polysaccharides',
        'Genome NC_004350.2 (pathway feature; not a single receptor-protein accession)',
        'PUBLISHED_GENOME_ANNOTATION',
        'UA159 genome study describes an RGP biosynthesis operon; M102AD adsorption experiments implicated serotype-specific cell-wall carbohydrates.',
        'This does not establish the precise receptor in a submitted strain or confirm that SpaP mediates adsorption.',
        'Ajdic et al. 2002; Delisle et al. 2012', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC137901/'),
    EvidenceRecord(
        'ORAL-REVIEW', 'S. mutans (general)', 'Biofilm and EPS context',
        'DOI 10.3389/fmolb.2026.1766853', 'REVIEW_BACKGROUND',
        'Review discusses glucosyltransferases, extracellular polysaccharides and oral microbiome interventions.',
        'Secondary review: not direct experimental validation of a phage-RBP-receptor interaction.',
        'Radhamanalan 2026', 'https://doi.org/10.3389/fmolb.2026.1766853'),
)


def get_reference_records(category='All'):
    """Return fixed published-reference examples; never label user input with these."""
    if category not in ('All', 'Phage', 'Bacterial', 'Background'):
        raise ValueError('Unsupported reference category')
    records = EVIDENCE
    if category == 'Phage':
        records = tuple(x for x in records if x.organism_or_phage.startswith('S. mutans phage'))
    if category == 'Bacterial':
        records = tuple(x for x in records if x.organism_or_phage == 'S. mutans UA159')
    if category == 'Background':
        records = tuple(x for x in records if x.evidence_class == 'REVIEW_BACKGROUND')
    return [asdict(x) for x in records]


def evidence_category_counts():
    counts = {}
    for row in EVIDENCE:
        counts[row.evidence_class] = counts.get(row.evidence_class, 0) + 1
    return counts
