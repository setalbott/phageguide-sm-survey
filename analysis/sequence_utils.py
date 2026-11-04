"""FASTA validation and descriptive summaries (no biological predictions)."""
from dataclasses import dataclass
from io import StringIO
from Bio import SeqIO

DNA_ALPHABET = set("ACGTRYSWKMBDHVN")  # IUPAC DNA ambiguity codes
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
MAX_INPUT_CHARACTERS = 10_000_000


class SequenceInputError(ValueError):
    """Invalid or unsupported sequence input."""


@dataclass(frozen=True)
class SequenceSummary:
    record_count: int
    total_length: int
    identifiers: tuple[str, ...]
    sequence_type: str
    gc_percent: float | None
    ambiguous_bases: int | None
    type_overlap_warning: bool = False


def decode_fasta_upload(content: bytes) -> str:
    """Decode an uploaded FASTA as UTF-8 text with a size limit."""
    if len(content) > MAX_INPUT_CHARACTERS:
        raise SequenceInputError("Uploaded FASTA is too large (10 MB maximum in this prototype).")
    try:
        decoded = content.decode("utf-8-sig")  # accepts normal UTF-8 or UTF-8 BOM
    except UnicodeDecodeError as exc:
        raise SequenceInputError("File must be plain-text UTF-8 FASTA, not binary or another encoding.") from exc
    if "\x00" in decoded:
        raise SequenceInputError("File contains NUL/binary bytes; upload a plain-text FASTA file.")
    return decoded


def summarize_fasta(fasta_text: str | None, sequence_type: str) -> SequenceSummary:
    """Validate FASTA syntax and symbols, then describe all sequence records.

    The user must choose 'dna' or 'protein'. Sequence letters cannot reliably
    establish biological origin or distinguish every DNA string from protein.
    GC% denominator is unambiguous A/C/G/T bases only.
    """
    if sequence_type not in {"dna", "protein"}:
        raise SequenceInputError("Choose DNA or protein as the input type.")
    if not isinstance(fasta_text, str) or not fasta_text.strip():
        raise SequenceInputError("Please provide a FASTA sequence.")
    if len(fasta_text) > MAX_INPUT_CHARACTERS:
        raise SequenceInputError("FASTA text is too large (10 million characters maximum).")
    lines = [line.strip() for line in fasta_text.lstrip("\ufeff").splitlines() if line.strip()]
    if not lines[0].startswith(">"):
        raise SequenceInputError("FASTA must begin with a header line starting with >.")
    alphabet = DNA_ALPHABET if sequence_type == "dna" else PROTEIN_ALPHABET
    identifiers: list[str] = []
    sequences: list[str] = []
    current: list[str] | None = None
    for line_number, line in enumerate(lines, start=1):
        if line.startswith(">"):
            header = line[1:].strip()
            if not header or header.startswith(">"):
                raise SequenceInputError(f"Missing/invalid FASTA identifier near nonblank line {line_number}.")
            if current is not None:
                sequences.append("".join(current))
            identifiers.append(header.split()[0])
            current = []
        else:
            if current is None:
                raise SequenceInputError("Sequence text appears before the first FASTA header.")
            if line.startswith(";"):
                raise SequenceInputError("FASTA comments are not supported in this version.")
            if any(char.isspace() for char in line):
                raise SequenceInputError(f"Whitespace inside sequence at nonblank line {line_number} is unsupported.")
            if ">" in line:
                raise SequenceInputError(f"Unexpected '>' in sequence at nonblank line {line_number}.")
            current.append(line)
    if current is not None:
        sequences.append("".join(current))
    for identifier, seq in zip(identifiers, sequences):
        if not seq:
            raise SequenceInputError(f"Record '{identifier}' has no sequence.")
        invalid = sorted(set(seq.upper()) - alphabet)
        if invalid:
            raise SequenceInputError(
                f"Record '{identifier}' contains invalid {sequence_type.upper()} characters: "
                + ", ".join(repr(c) for c in invalid)
            )

    # Parse through Biopython after stricter checks; do not depend on its
    # permissive handling of empty or malformed records.
    normalized = "\n".join(lines) + "\n"
    try:
        records = list(SeqIO.parse(StringIO(normalized), "fasta"))
    except (ValueError, TypeError) as exc:
        raise SequenceInputError(f"Could not parse FASTA: {exc}") from exc
    if len(records) != len(identifiers):
        raise SequenceInputError("FASTA record count did not match parsed headers.")
    total = sum(len(s) for s in sequences)
    gc_count = sum(s.upper().count("G") + s.upper().count("C") for s in sequences)
    canonical_count = sum(sum(s.upper().count(c) for c in "ACGT") for s in sequences)
    ambiguous = sum(sum(s.upper().count(c) for c in "RYSWKMBDHVN") for s in sequences)
    gc = round(100 * gc_count / canonical_count, 2) if sequence_type == "dna" and canonical_count else None
    overlap = sequence_type == "protein" and all(set(s.upper()) <= DNA_ALPHABET for s in sequences)
    return SequenceSummary(
        record_count=len(records), total_length=total, identifiers=tuple(identifiers),
        sequence_type=sequence_type, gc_percent=gc,
        ambiguous_bases=ambiguous if sequence_type == "dna" else None,
        type_overlap_warning=overlap,
    )
