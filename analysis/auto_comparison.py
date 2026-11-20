"""Optional reference protein comparison triggered by Analyze.

No biological function, binding probability, or susceptibility is inferred.
"""
from analysis.sequence_comparison import compare_protein_sequences


def automatic_comparison(phage_fasta, phage_type, minimum_aa, reference_fasta):
    """Return a structured alignment report or an informative unavailable state.

    Input FASTA validation for both phage and bacterium is performed by the UI before
    calling this function. Protein reference validation is delegated to the existing
    sequence comparison module.
    """
    if not reference_fasta or not reference_fasta.strip():
        return {"status": "no_reference", "rows": [],
                "message": "No reference protein provided. Add a protein FASTA to calculate alignment matches, substitutions, gaps and percent identity."}
    rows = compare_protein_sequences(
        phage_fasta, phage_type.lower(), reference_fasta,
        accession="USER-OR-SYNTHETIC-REFERENCE",
        provenance="Reference FASTA provided in Analyze (not independently verified)",
        minimum_aa=int(minimum_aa),
    )
    if not rows:
        return {"status": "no_orf", "rows": [],
                "message": "No qualifying complete ORF was found at the selected minimum length, so no translated phage protein can be compared. Lower the ORF minimum or supply phage protein FASTA."}
    return {"status": "ok", "rows": rows, "message": "Protein sequence alignment completed. This does not establish susceptibility."}
