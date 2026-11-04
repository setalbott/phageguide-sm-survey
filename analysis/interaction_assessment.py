"""Stage 5: conservative evidence assessment, not a host-range predictor.

Sequence comparison is exact protein equality only. Reference FASTA provenance is
user-supplied, not verified by this application. Literature scenarios are curated
examples and never attached to submitted sequences automatically.
"""
from dataclasses import dataclass, asdict
from io import StringIO
from Bio import SeqIO
from analysis.sequence_utils import summarize_fasta
from analysis.orf_analysis import screen_phage_orfs

SOURCE = 'Delisle et al. (2012), Applied and Environmental Microbiology, DOI:10.1128/AEM.07726-11'
SOURCE_URL = 'https://doi.org/10.1128/AEM.07726-11'

@dataclass(frozen=True)
class Assessment:
    receptor_status: str
    lysis_status: str
    reason: str
    sequence_evidence: str
    experimental_evidence: str
    source: str
    source_url: str
    evidence_limit: str


def assess_uploaded_pair() -> Assessment:
    """No direct host receptor/lysis conclusions from the current inputs."""
    return Assessment(
        'Insufficient evidence', 'Insufficient evidence',
        'Neither FASTA validation, ORF detection, protein sequence identity, nor general literature establishes a receptor or lysis mechanism for this submitted phage–strain pair.',
        'No verified, strain-linked receptor interaction has been established for the submitted sequences.',
        'No strain-matched phage challenge or lysis experiment was provided or independently verified.',
        'Methods and limitations of PhageGuide-SM Stage 5', '',
        'There is no trained or calibrated prediction model; no confidence percentage is produced.',
    )


def reference_scenario(scenario: str) -> dict:
    """Published observations for named study cases, NOT the user-uploaded FASTA."""
    cases = {
        'M102AD / P42S': {
            'phage': 'M102AD', 'bacterial_strain': 'P42S (serotype c)',
            'observation': 'M102AD was grown on strain P42S and its adsorption to P42S was studied.',
            'evidence_class': 'EXPERIMENTAL_HOST_PHENOTYPE',
            'receptor_result': 'Insufficient evidence: serotype-specific surface carbohydrates implicated, exact molecular receptor unresolved.',
            'lysis_result': 'Published phage growth/host interaction for this named strain, not an identified molecular lysis pathway.',
            'limitation': 'Not inferred for any uploaded FASTA, and not proof of a specific protein–receptor contact.',
        },
        'M102AD / UA159': {
            'phage': 'M102AD', 'bacterial_strain': 'UA159',
            'observation': 'UA159 was in the tested host-range panel; only strain F of the additional 25 strains was susceptible.',
            'evidence_class': 'EXPERIMENTAL_HOST_PHENOTYPE',
            'receptor_result': 'Insufficient evidence: no specific molecular receptor established.',
            'lysis_result': 'UA159 was not susceptible in the reported assay; this is not a molecular lysis-pathway assignment.',
            'limitation': 'Observed under the reported experimental conditions; does not automatically identify uploaded sequences.',
        },
    }
    if scenario not in cases:
        raise ValueError('Unknown published scenario')
    return {**cases[scenario], 'source': SOURCE, 'source_url': SOURCE_URL}


def compare_reference_protein(phage_text: str, phage_kind: str, reference_fasta: str,
                              accession: str, provenance: str, minimum_aa: int = 30) -> list[dict]:
    """Exact AA equality only; no functional transfer or identity verification.

    For DNA, translated complete ORFs from Stage 3B are compared with supplied
    AA sequences. For phage protein input, whole FASTA records are compared.
    Output is an exact sequence comparison, not a functional/host prediction.
    """
    if not accession.strip() or not provenance.strip():
        raise ValueError('Reference accession and provenance are required.')
    if phage_kind not in ('dna', 'protein'):
        raise ValueError('Choose dna or protein for phage input.')
    summarize_fasta(reference_fasta, 'protein')
    summarize_fasta(phage_text, phage_kind)
    refs = [(rec.id, str(rec.seq).upper().rstrip('*')) for rec in SeqIO.parse(StringIO(reference_fasta.lstrip('\ufeff')), 'fasta')]
    if any('*' in sequence for _, sequence in refs):
        raise ValueError('Reference proteins must not contain internal stop symbols.')
    if phage_kind == 'protein':
        queries = [(rec.id, str(rec.seq).upper().rstrip('*'), 'uploaded protein')
                   for rec in SeqIO.parse(StringIO(phage_text.lstrip('\ufeff')), 'fasta')]
    else:
        queries = [(f'{orf.record_id}:{orf.strand}:{orf.start}-{orf.end}', orf.protein_sequence,
                    'Stage 3B ATG-to-stop ORF translation')
                   for orf in screen_phage_orfs(phage_text, min_amino_acids=minimum_aa)]
    output = []
    for q_id, q_seq, method in queries:
        for ref_id, ref_seq in refs:
            if q_seq == ref_seq and q_seq and '*' not in q_seq:
                output.append({
                    'query_id': q_id, 'reference_id': ref_id,
                    'reference_accession_user_supplied': accession.strip(),
                    'reference_provenance_user_supplied': provenance.strip(),
                    'comparison': 'EXACT_AMINO_ACID_SEQUENCE_MATCH',
                    'length_aa': len(q_seq), 'query_method': method,
                    'interpretation': 'Identical AA letters only; reference identity, protein function, receptor binding, and strain activity are NOT verified.',
                })
    return output
