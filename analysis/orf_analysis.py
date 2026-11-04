"""Stage 3B: simple six-frame ATG-to-stop phage ORF screening.

Not gene calling, protein functional annotation, or RBP prediction.
Coordinates are 1-based inclusive relative to each supplied FASTA record.
"""
from dataclasses import dataclass
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from analysis.sequence_utils import summarize_fasta

STOP_CODONS = {'TAA', 'TAG', 'TGA'}


@dataclass(frozen=True)
class ORFCandidate:
    record_id: str
    strand: str
    frame: int
    start: int
    end: int
    nucleotide_length: int
    protein_length_aa: int
    protein_sequence: str
    stop_codon: str
    evidence: str = 'Sequence-based ATG-to-stop ORF; function not established'


def screen_phage_orfs(fasta_text: str, min_amino_acids: int = 30) -> list[ORFCandidate]:
    """Find first-start ATG-to-stop ORFs on both strands of each DNA record.

    In each reading frame, uses the first ATG since the previous stop, and
    emits an ORF only if it reaches an in-frame stop. Alternative start sites,
    overlapping ATGs in one stop-delimited segment, partial and circular ORFs
    are not assessed. Ambiguous codons do not count as starts/stops.
    """
    if type(min_amino_acids) is not int or min_amino_acids < 1:
        raise ValueError('Minimum protein length must be a positive integer.')
    summarize_fasta(fasta_text, 'dna')
    candidates = []
    for record in SeqIO.parse(StringIO(fasta_text.lstrip('\ufeff')), 'fasta'):
        forward = str(record.seq).upper()
        n = len(forward)
        for strand, dna in [('+', forward), ('-', str(Seq(forward).reverse_complement()))]:
            for offset in range(3):
                start_index = None
                for i in range(offset, len(dna) - 2, 3):
                    codon = dna[i:i+3]
                    if codon in STOP_CODONS:
                        if start_index is not None:
                            coding = dna[start_index:i]
                            protein = str(Seq(coding).translate(table=1))
                            if len(protein) >= min_amino_acids:
                                if strand == '+':
                                    begin, end = start_index + 1, i + 3
                                else:
                                    begin, end = n - (i + 3) + 1, n - start_index
                                candidates.append(ORFCandidate(
                                    record_id=record.id, strand=strand, frame=offset + 1,
                                    start=begin, end=end,
                                    nucleotide_length=i + 3 - start_index,
                                    protein_length_aa=len(protein),
                                    protein_sequence=protein, stop_codon=codon,
                                ))
                        start_index = None
                    elif codon == 'ATG' and start_index is None:
                        start_index = i
    return sorted(candidates, key=lambda x: (x.record_id, x.start, x.end, x.strand))
