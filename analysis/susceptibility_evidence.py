"""Literature-derived *assay outcomes*, not predicted interactions.

Scope: Ben-Zaken et al. Viruses 2021 Table 3; Delisle et al. AEM 2012
host-range results; Sugai et al. iScience 2023 Table 1 aggregates.
Do not infer the identity of an uploaded genome or use records as labels
for untested combinations.
"""
from dataclasses import asdict, dataclass

SMHBZ8_DOI = 'https://doi.org/10.3390/v13050825'
M102AD_DOI = 'https://doi.org/10.1128/AEM.07726-11'
KSM96_DOI = 'https://doi.org/10.1016/j.isci.2023.108465'
NO_RECORD = 'No published result in this dataset.'

@dataclass(frozen=True)
class SusceptibilityRecord:
    phage: str
    strain: str
    assay: str
    observed_result: str
    outcome: str
    source: str
    source_url: str
    evidence_type: str = 'PUBLISHED_EXPERIMENTAL_HOST_RANGE'
    review_flag: str = ''
    notes: str = ''


def _smhbz8():
    positive = ['ATCC 700610', 'ATCC 27351'] + [f'ES{i}' for i in range(1, 11)] + ['MT8148']
    negative = ['MT703', 'LM7', 'OMZ175', 'M76219']
    records = []
    for strain in positive + negative:
        flag = ('IDENTIFIER CONFLICT: Table 3 prints M76219; Methods section prints MT6219. '
                'Review original strain identity before linking to genome accessions.') if strain == 'M76219' else ''
        records.append(SusceptibilityRecord(
            'SMHBZ8', strain, 'Qualitative phage spot/plaque assay (presence or absence of plaques)',
            'Plaques observed' if strain in positive else 'No plaques observed',
            'Sensitive (S)' if strain in positive else 'Resistant (R) under tested conditions',
            'Ben-Zaken et al. (2021), Viruses 13:825, Table 3', SMHBZ8_DOI,
            review_flag=flag,
            notes='Serotype c' if strain in positive else ('Serotype f' if strain in ('OMZ175','M76219') else 'Serotype e')
        ))
    return records


def _m102ad():
    # 25 *additional* strains in Delisle et al. 2012 Results; only F sensitive.
    tested = ('CM7', 'ATCC 202022', '7A', '123.1', 'Q', '59.1', 'UA96', '29B',
              'T9', 'T8', 'N', '87.1', '5A', 'ATCC 25175', '11A', '136.1',
              'JB', '3B', 'Ny257', '22A', 'F', 'NG003', 'T7B', 'JH 1140', 'UA159')
    records = [SusceptibilityRecord(
        'M102AD', 'P42S', 'Propagation host; adsorption kinetics and plaque-based phage growth',
        'Established propagation host; adsorption measured', 'Susceptible host (propagation)',
        'Delisle et al. (2012), Applied and Environmental Microbiology, Methods; Table 1', M102AD_DOI,
        notes='Separate established propagation host, not a member of the 25 additional-strain host-range panel.'
    )]
    for strain in tested:
        records.append(SusceptibilityRecord(
            'M102AD', strain, 'Host-range plaque assay (spotting serial phage dilutions)',
            'Plaques/sensitivity reported' if strain == 'F' else 'No susceptibility reported in host-range panel',
            'Sensitive' if strain == 'F' else 'Not susceptible under tested conditions',
            'Delisle et al. (2012), Applied and Environmental Microbiology, Host range Results', M102AD_DOI,
            notes='One of 25 additionally tested S. mutans strains; only F was sensitive.'
        ))
    return records


RECORDS = tuple(_smhbz8() + _m102ad())

# Aggregate outcomes; never treat the counts as 125 named individual records.
KSM96_TOTALS = (
    {'phage':'φKSM96 (4KSM96)', 'serotype':'c', 'tested':105, 'susceptible':102, 'insensitive':3},
    {'phage':'φKSM96 (4KSM96)', 'serotype':'e', 'tested':14, 'susceptible':11, 'insensitive':3},
    {'phage':'φKSM96 (4KSM96)', 'serotype':'k', 'tested':4, 'susceptible':0, 'insensitive':4},
    {'phage':'φKSM96 (4KSM96)', 'serotype':'f', 'tested':2, 'susceptible':0, 'insensitive':2},
)


def get_susceptibility_records(phage='All', strain='All'):
    """Filter exact published strain identifiers; no biological prediction or fuzzy matching."""
    if not isinstance(phage, str) or not isinstance(strain, str):
        raise ValueError('Phage and strain filters must be strings')
    if not phage.strip() or not strain.strip():
        raise ValueError('Phage and strain filters cannot be empty')
    return [asdict(r) for r in RECORDS
            if (phage == 'All' or r.phage.casefold() == phage.casefold())
            and (strain == 'All' or r.strain.casefold() == strain.casefold())]


def lookup_susceptibility(phage, strain):
    """Returns only exact named-pair records, or an explicit no-record status."""
    results = get_susceptibility_records(phage, strain)
    return {'status': 'Published result found' if results else NO_RECORD, 'records': results}


def get_ksm96_serotype_totals():
    return [dict(r, source='Sugai et al. (2023), iScience 26:108465, Table 1',
                 source_url=KSM96_DOI,
                 evidence_type='PUBLISHED_AGGREGATE_HOST_RANGE') for r in KSM96_TOTALS]
