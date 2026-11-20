"""Exact two-genome reference matching, never an infection predictor.

NCBI nucleotide reference sequences are fetched only when needed, cached by
Streamlit in the UI. Fail closed when a reference cannot be obtained or checked.
"""
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from Bio import SeqIO
from analysis.sequence_utils import summarize_fasta, SequenceInputError

PHAGE_ACCESSION = "KR153145.1"     # GenBank phage phiAPCM01 (Dalmasso 2015)
BACTERIUM_ACCESSION = "NZ_CP172847.1"  # RefSeq DPC6143 (Wall and Wall 2025)
REFS = {"phage": PHAGE_ACCESSION, "bacterium": BACTERIUM_ACCESSION}


@dataclass(frozen=True)
class ReferenceMatch:
    label: str
    accession: str
    status: str
    uploaded_length: int | None
    reference_length: int | None
    explanation: str


def parse_single_dna(text):
    """Require one full DNA FASTA record; never trust its header as identity."""
    summarize_fasta(text, "dna")
    records = list(SeqIO.parse(StringIO(text), "fasta"))
    if len(records) != 1:
        raise SequenceInputError("A single complete genome FASTA record is required for reference matching.")
    return str(records[0].seq).upper()


def fetch_ncbi_reference(accession):
    """Fetch accession.version as FASTA; do not infer matches from headers/length."""
    if accession not in REFS.values():
        raise ValueError("Accession is not a curated benchmark reference.")
    params = urlencode({"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"})
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
    request = Request(url, headers={"User-Agent": "PhageGuide-SM research prototype (reference verification)", "Accept": "text/plain"})
    with urlopen(request, timeout=18) as resp:
        text = resp.read(7_000_000).decode("utf-8")
    seq = parse_single_dna(text)
    record_id = next(SeqIO.parse(StringIO(text), "fasta")).id
    if record_id.upper() != accession.upper():
        raise ValueError(f"NCBI returned {record_id}, not requested {accession}.")
    return seq


def match_uploaded_fasta(label, uploaded_fasta, reference_sequence, accession):
    """Exact full-genome match only; a partial match never inherits study claims."""
    try:
        uploaded = parse_single_dna(uploaded_fasta)
    except (SequenceInputError, ValueError) as exc:
        return ReferenceMatch(label, accession, "invalid", None, len(reference_sequence) if reference_sequence else None, str(exc))
    if reference_sequence is None:
        return ReferenceMatch(label, accession, "reference_unavailable", len(uploaded), None,
            "The accession-specific reference could not be retrieved and verified. Check NCBI connectivity; no sequence identity was asserted.")
    reference = reference_sequence.upper()
    if uploaded == reference:
        return ReferenceMatch(label, accession, "matched", len(uploaded), len(reference),
            "Exact nucleotide-for-nucleotide match to the complete accession-version reference sequence.")
    if len(uploaded) != len(reference):
        detail = f"Length differs ({len(uploaded):,} vs {len(reference):,} nt)."
    else:
        differences = sum(a != b for a, b in zip(uploaded, reference))
        detail = f"Same length but {differences:,} nucleotide positions differ."
    return ReferenceMatch(label, accession, "not_matched", len(uploaded), len(reference),
        detail + " No study-specific findings can be assigned to this uploaded sequence.")


def assess_reference_pair(phage_fasta, bacteria_fasta, reference_phage, reference_bacterium):
    phage = match_uploaded_fasta("Phage φAPCM01", phage_fasta, reference_phage, PHAGE_ACCESSION)
    bacterium = match_uploaded_fasta("S. mutans DPC6143", bacteria_fasta, reference_bacterium, BACTERIUM_ACCESSION)
    both = phage.status == bacterium.status == "matched"
    return {"phage": phage, "bacterium": bacterium, "both_matched": both,
            "status": "verified_reference_pair" if both else "incomplete_reference_match"}
