"""Small protein-sequence global alignments for a teaching/survey interface.

Alignment is a comparison of letters, not evidence of protein function or binding.
The size and pair limits avoid quadratic cost for public Streamlit usage.
"""
from io import StringIO
from Bio import SeqIO, Align
from analysis.sequence_utils import summarize_fasta

MAX_SEQUENCE_AA = 2000
MAX_COMPARISONS = 40


def compare_protein_sequences(phage_text, phage_kind, reference_fasta, accession,
                              provenance, minimum_aa=30):
    """Globally align protein records or translated phage ORFs against reference proteins.

    Returns dictionaries with explicit aligned strings, symbols, and position counts.
    Positions are 1-based in the ungapped protein sequences; an alignment does not
    imply either sequence is experimentally verified or homologous.
    """
    from analysis.orf_analysis import screen_phage_orfs
    if not accession or not accession.strip() or not provenance or not provenance.strip():
        raise ValueError('Reference identifier and provenance are required.')
    if phage_kind not in ('dna', 'protein'):
        raise ValueError('Phage input type must be dna or protein.')
    summarize_fasta(phage_text, phage_kind)
    summarize_fasta(reference_fasta, 'protein')
    refs = [(r.id, str(r.seq).upper().rstrip('*')) for r in SeqIO.parse(StringIO(reference_fasta.lstrip('\ufeff')), 'fasta')]
    if phage_kind == 'protein':
        queries = [(r.id, str(r.seq).upper().rstrip('*')) for r in SeqIO.parse(StringIO(phage_text.lstrip('\ufeff')), 'fasta')]
    else:
        queries = [(f'{o.record_id}:{o.strand}:{o.start}-{o.end}', o.protein_sequence)
                   for o in screen_phage_orfs(phage_text, min_amino_acids=int(minimum_aa))]
    if not queries:
        return []
    if len(queries) * len(refs) > MAX_COMPARISONS:
        raise ValueError(f'Too many comparisons for this interactive demo (limit {MAX_COMPARISONS} protein pairs). Select fewer records.')
    for _, seq in queries + refs:
        if '*' in seq:
            raise ValueError('Internal stop symbols cannot be aligned as a standard protein sequence in this demo.')
        if len(seq) > MAX_SEQUENCE_AA:
            raise ValueError(f'Sequence exceeds the survey comparison limit of {MAX_SEQUENCE_AA} amino acids. Use a shorter protein or dedicated alignment software.')
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                                   open_gap_score=-3, extend_gap_score=-0.5)
    output = []
    for qid, query in queries:
        for rid, reference in refs:
            aln = aligner.align(query, reference)[0]
            qa, ra, bars = [], [], []
            matches = substitutions = gaps = 0
            for qi, ri in zip(aln.indices[0], aln.indices[1]):
                a = query[qi] if qi >= 0 else '-'
                b = reference[ri] if ri >= 0 else '-'
                qa.append(a); ra.append(b)
                if a == b and a != '-':
                    bars.append('|'); matches += 1
                elif a == '-' or b == '-':
                    bars.append(' '); gaps += 1
                else:
                    bars.append('.'); substitutions += 1
            output.append({
                'query_id': qid, 'reference_id': rid,
                'reference_accession_user_supplied': accession.strip(),
                'reference_provenance_user_supplied': provenance.strip(),
                'query_aligned': ''.join(qa), 'match_symbols': ''.join(bars),
                'reference_aligned': ''.join(ra), 'match_count': matches,
                'difference_count': substitutions + gaps, 'substitution_count': substitutions,
                'gap_positions': gaps, 'alignment_positions': len(qa),
                'percent_identity_over_alignment': 100*matches/len(qa) if qa else 0.,
                'exact_match': query == reference,
                'interpretation': 'Sequence-level comparison only: not a binding probability, receptor confirmation, functional annotation, or lysis prediction.',
            })
    return output


def format_alignment(result, width=60):
    """Return readable three-line alignment blocks with 1-based ungapped offsets."""
    if width < 1:
        raise ValueError('Block width must be positive.')
    q = result['query_aligned']; r = result['reference_aligned']; m = result['match_symbols']
    qpos = rpos = 0
    blocks = []
    for start in range(0, len(q), width):
        a = q[start:start+width]; b = r[start:start+width]; sym = m[start:start+width]
        qa0 = qpos + 1; ra0 = rpos + 1
        qpos += len(a.replace('-', '')); rpos += len(b.replace('-', ''))
        blocks.append(f'Input     {qa0:>5}  {a}  {qpos}\n'
                      f'                 {sym}\n'
                      f'Reference {ra0:>5}  {b}  {rpos}')
    return '\n\n'.join(blocks)
