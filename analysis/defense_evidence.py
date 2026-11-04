"""Stage 4: auditable defense literature and optional imported DefenseFinder calls.

A user's DNA FASTA is NEVER scanned for defense genes by this module.
Only rows explicitly imported from an external DefenseFinder systems TSV are
reported as predictions. Published cohort results are NOT mapped to user input.
"""
from dataclasses import dataclass, asdict
import csv
import io
import re

MAX_TSV_BYTES = 3_000_000
REQUIRED = ('sys_id', 'type', 'subtype', 'sys_beg', 'sys_end', 'protein_in_syst', 'genes_count', 'name_of_profiles_in_sys')


class DefenseImportError(ValueError):
    pass


@dataclass(frozen=True)
class PublishedDefenseFinding:
    organism: str
    scope: str
    finding: str
    evidence_class: str
    source: str
    url: str
    limitation: str


LITERATURE = (
    PublishedDefenseFinding(
        'Streptococcus mutans', '44 clinical isolates (NOT UA159 or an uploaded sequence)',
        'Computational surveys reported 6–12 anti-phage defense systems per isolate.',
        'PUBLISHED_COMPUTATIONAL_GENOME_SCREEN',
        'Claisse et al., Molecular Oral Microbiology (2026), DOI: 10.1111/omi.70014',
        'https://doi.org/10.1111/omi.70014',
        'Derived from DefenseFinder/PADLOC bioinformatics, not experimentally proven active in each isolate.'),
    PublishedDefenseFinding(
        'Streptococcus mutans', '44 clinical isolates (NOT UA159 or an uploaded sequence)',
        'Restriction-modification systems were most prevalent, followed by MazEF and CRISPR-Cas; all 44 carried at least one restriction-modification system and MazEF in the study.',
        'PUBLISHED_COMPUTATIONAL_GENOME_SCREEN',
        'Claisse et al. (2026), Results 3.1',
        'https://doi.org/10.1111/omi.70014',
        'Presence inferred from genomic detection; does not establish effectiveness against a particular phage.'),
    PublishedDefenseFinding(
        'Streptococcus mutans', '44 clinical isolates (NOT UA159 or an uploaded sequence)',
        'Bioinformatics identified CRISPR-Cas subtypes I-C, I-E, and II-A among the isolates.',
        'PUBLISHED_COMPUTATIONAL_GENOME_SCREEN',
        'Claisse et al. (2026), Results',
        'https://doi.org/10.1111/omi.70014',
        'Do not assign these subtypes to UA159 or the submitted genome without strain-specific evidence.'),
)


def published_defense_findings():
    return [asdict(row) for row in LITERATURE]


def decode_defense_tsv(blob: bytes) -> str:
    if not isinstance(blob, bytes):
        raise DefenseImportError('Upload must be file bytes.')
    if not blob:
        raise DefenseImportError('The DefenseFinder TSV is empty.')
    if len(blob) > MAX_TSV_BYTES:
        raise DefenseImportError('DefenseFinder TSV exceeds the 3 MB upload limit.')
    try:
        text = blob.decode('utf-8-sig')
    except UnicodeDecodeError as e:
        raise DefenseImportError('DefenseFinder TSV must be UTF-8 text.') from e
    if '\x00' in text:
        raise DefenseImportError('Binary content is not a valid TSV.')
    return text


def parse_defensefinder_systems(text: str, dataset_accession: str, tool_version: str, source_description: str):
    """Parse official defense_finder_systems.tsv format without interpreting activity.

    These are reported as *imported predicted systems*, not independently rerun.
    Dataset metadata is explicitly user-supplied and not independently verified.
    """
    for label, value in [('Dataset accession or sample ID', dataset_accession), ('DefenseFinder version', tool_version), ('Input/source description', source_description)]:
        if not isinstance(value, str) or not value.strip():
            raise DefenseImportError(f'{label} is required to retain provenance.')
        if len(value) > 300:
            raise DefenseImportError(f'{label} is too long.')
    if not isinstance(text, str) or not text.strip():
        raise DefenseImportError('The DefenseFinder TSV is empty.')
    if '\x00' in text or len(text.encode('utf-8')) > MAX_TSV_BYTES:
        raise DefenseImportError('Invalid or oversized TSV.')
    try:
        reader = csv.DictReader(io.StringIO(text, newline=''), delimiter='\t')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise DefenseImportError('Duplicate TSV column names.')
        missing = [key for key in REQUIRED if key not in fields]
        if missing:
            raise DefenseImportError('Missing required DefenseFinder columns: ' + ', '.join(missing))
        rows = []
        for index, item in enumerate(reader, 1):
            if index > 10000:
                raise DefenseImportError('More than 10,000 calls in TSV; import limited.')
            if None in item or any(k is None for k in item):
                raise DefenseImportError(f'Malformed TSV row {index}: extra columns.')
            if any(v is None for v in item.values()):
                raise DefenseImportError(f'Malformed TSV row {index}: missing fields.')
            if not item['sys_id'].strip() or not item['type'].strip():
                raise DefenseImportError(f'Malformed TSV row {index}: system identifier/type missing.')
            try:
                gene_count = int(item['genes_count'].strip())
            except ValueError as e:
                raise DefenseImportError(f'Invalid genes_count in row {index}.') from e
            if gene_count < 1:
                raise DefenseImportError(f'Invalid genes_count in row {index}.')
            rows.append({
                'system_id': item['sys_id'].strip(), 'system_type': item['type'].strip(),
                'subtype': item['subtype'].strip(), 'begin_protein': item['sys_beg'].strip(),
                'end_protein': item['sys_end'].strip(), 'proteins': item['protein_in_syst'].strip(),
                'genes_count': gene_count, 'profiles': item['name_of_profiles_in_sys'].strip(),
                'evidence_class': 'IMPORTED_DEFENSEFINDER_PREDICTION',
                'dataset_id_user_supplied': dataset_accession.strip(),
                'defensefinder_version_user_supplied': tool_version.strip(),
                'source_user_supplied': source_description.strip(),
                'interpretation': 'Predicted genomic system only; activity and phage resistance not established.'
            })
        return rows
    except csv.Error as e:
        raise DefenseImportError('Could not parse tab-delimited DefenseFinder results.') from e
