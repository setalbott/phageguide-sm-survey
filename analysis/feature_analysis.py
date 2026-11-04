"""Stage 3A: descriptive sequence features only; no receptor predictions."""
from dataclasses import dataclass
from io import StringIO
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from analysis.sequence_utils import summarize_fasta


@dataclass(frozen=True)
class RecordFeature:
    identifier: str
    length: int
    molecule_type: str
    gc_percent: float | None
    ambiguous_positions: int | None
    molecular_weight_da: float | None
    weight_note: str


def describe_records(fasta_text: str, sequence_type: str) -> list[RecordFeature]:
    """Validate Stage 2 input and calculate reproducible properties per record.

    Molecular weight is computed only for proteins composed of the twenty
    standard amino acids (without stop/ambiguity markers).
    """
    summarize_fasta(fasta_text, sequence_type)  # reuse established validation
    records = list(SeqIO.parse(StringIO(fasta_text.lstrip('\ufeff')), 'fasta'))
    output = []
    for rec in records:
        seq = str(rec.seq).upper()
        if sequence_type == 'dna':
            canonical = sum(seq.count(base) for base in 'ACGT')
            gc = round(100 * (seq.count('G') + seq.count('C')) / canonical, 2) if canonical else None
            ambiguity = len(seq) - canonical
            weight, note = None, 'Not applicable for DNA input'
        else:
            gc, ambiguity = None, None
            if set(seq) <= set('ACDEFGHIKLMNPQRSTVWY'):
                weight = round(ProteinAnalysis(seq).molecular_weight(), 2)
                note = 'Estimated from standard amino-acid sequence'
            else:
                weight = None
                note = 'Not calculated: ambiguous/nonstandard residues or stop symbol present'
        output.append(RecordFeature(rec.id, len(seq), sequence_type, gc, ambiguity, weight, note))
    return output
