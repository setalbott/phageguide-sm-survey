"""PhageGuide-SM Stage 3A: descriptive per-record features."""
import streamlit as st
from analysis.sequence_utils import SequenceInputError, summarize_fasta, decode_fasta_upload
from analysis.feature_analysis import describe_records
from analysis.orf_analysis import screen_phage_orfs
from analysis.reference_evidence import get_reference_records, evidence_category_counts
from analysis.defense_evidence import published_defense_findings, parse_defensefinder_systems, decode_defense_tsv, DefenseImportError
from analysis.interaction_assessment import assess_uploaded_pair, reference_scenario, compare_reference_protein

st.set_page_config(page_title="PhageGuide-SM", page_icon="🧬", layout="wide")

st.title("🧬 PhageGuide-SM")
st.caption("Survey prototype | Stage 5 evidence-assessment demo — educational research tool")
st.info(
    "This prototype assesses available evidence, but does NOT predict a specific molecular receptor or lysis pathway for uploaded genomes. Insufficient evidence is a valid result."
)

with st.sidebar:
    st.header("About this prototype")
    st.write("Based on Figure 1 of the BIOL 831 PhageGuide-SM proposal.")
    st.write("Stage 3C adds a separately labeled, literature-derived reference catalog. It does not match or assign functions to uploaded sequences.")
    st.write("Reference genomes are not bundled yet. Upload your own genome FASTA.")

st.subheader("Try the guided survey example")
st.write("This is a **synthetic, non-biological demonstration**. Use the sample to test the interface, not to assess whether a real phage infects a bacterium.")
if st.button("Load sample DNA and bacterial FASTA", key="survey_load_sample", type="secondary"):
    st.session_state["survey_phage_type"] = "DNA"
    st.session_state["survey_phage_source"] = "Paste FASTA"
    st.session_state["survey_bacteria_source"] = "Paste FASTA"
    st.session_state["survey_phage_fasta"] = ">synthetic_forward_orf\nATGAAATAG"
    st.session_state["survey_bacteria_fasta"] = ">synthetic_bacteria\nATGCATGC"
    st.session_state["survey_orf_min"] = 2
    st.success("Sample loaded. Scroll to the inputs and click Analyze sequences. The DNA example has a two-amino-acid ORF, MK.")
with st.expander("Survey instructions — 4 short tasks", expanded=True):
    st.markdown("""1. Click **Load sample DNA and bacterial FASTA**.
2. Click **Analyze sequences**. Find the phage ORF `MK` at positions 1–9.
3. Scroll to the **Published reference evidence** and **Stage 4 defense** sections. Notice that the cited findings do **not** describe your synthetic input.
4. In **Stage 5**, switch between **M102AD / P42S** and **M102AD / UA159**. Find **Insufficient evidence** for the uploaded pair and explain what the exact-match comparison can and cannot show.

After testing, complete the separate survey link provided by the researcher. Do not submit real patient, unpublished, or confidential genome data in this public demo.""")

st.subheader("Step 1 · Add your sequences")
phage_col, bacteria_col = st.columns(2)

with phage_col:
    st.markdown("**Phage genome or protein**")
    phage_type = st.radio("Phage sequence type", ["DNA", "Protein"], horizontal=True, key="survey_phage_type")
    phage_source = st.radio("Phage input method", ["Upload FASTA", "Paste FASTA"], horizontal=True, key="survey_phage_source")
    if phage_source == "Upload FASTA":
        phage_file = st.file_uploader("Phage FASTA file", type=["fasta", "fa", "faa", "fna", "txt"], key="phage")
        phage_text = None
        if phage_file is not None:
            try:
                phage_text = decode_fasta_upload(phage_file.getvalue())
            except SequenceInputError as exc:
                st.error(f"Phage: {exc}")
    else:
        phage_text = st.text_area("Paste phage FASTA", height=160, placeholder=">phage_example\nATGCGC", key="survey_phage_fasta")

with bacteria_col:
    st.markdown("**Streptococcus mutans strain genome**")
    st.caption("A UA159 reference selector is planned for a later stage; the genome is not bundled.")
    bacteria_source = st.radio("Bacterial input method", ["Upload FASTA", "Paste FASTA"], horizontal=True, key="survey_bacteria_source")
    if bacteria_source == "Upload FASTA":
        bacteria_file = st.file_uploader("Bacterial genome FASTA file", type=["fasta", "fa", "fna", "txt"], key="bacteria")
        bacteria_text = None
        if bacteria_file is not None:
            try:
                bacteria_text = decode_fasta_upload(bacteria_file.getvalue())
            except SequenceInputError as exc:
                st.error(f"Bacterium: {exc}")
    else:
        bacteria_text = st.text_area("Paste bacterial genome FASTA", height=160, placeholder=">bacteria_example\nATGCATGC", key="survey_bacteria_fasta")

st.subheader("Optional · Phage DNA ORF screening")
orf_minimum = st.number_input("Minimum translated protein length (amino acids)", min_value=1, max_value=10000, value=30, step=1, help="30 aa is the default exploratory threshold; 2 aa is useful for short synthetic tests.", key="survey_orf_min")

if st.button("Analyze sequences", type="primary", use_container_width=True):
    st.subheader("Step 2 · Results")
    summaries = {}
    for label, text, kind in [
        ("Phage", phage_text, phage_type.lower()),
        ("Bacterium", bacteria_text, "dna"),
    ]:
        try:
            summaries[label] = summarize_fasta(text, kind)
        except SequenceInputError as exc:
            st.error(f"{label}: {exc}")

    if len(summaries) == 2:
        st.success("Both inputs passed FASTA validation.")
        cols = st.columns(2)
        for col, (label, result) in zip(cols, summaries.items()):
            with col:
                st.markdown(f"### {label} input")
                st.metric("Records", result.record_count)
                st.metric("Total sequence length", f"{result.total_length:,}")
                st.write("**Sequence type:**", result.sequence_type.upper())
                st.write("**Record IDs:**", ", ".join(result.identifiers[:5]))
                if result.record_count > 5:
                    st.caption(f"Showing first 5 of {result.record_count} IDs")
                if result.type_overlap_warning:
                    st.warning("This sequence uses only DNA-compatible letters. Letters alone cannot verify whether it is a protein or DNA sequence. Check your input type and source annotation.")
                if result.sequence_type == "dna":
                    st.write("**GC content:**", f"{result.gc_percent:.2f}%" if result.gc_percent is not None else "N/A (no unambiguous DNA bases)")
                    st.write("**Ambiguous nucleotide positions:**", result.ambiguous_bases)
        st.divider()
        st.subheader("Stage 3A · Descriptive sequence features")
        st.caption("These features are computed from sequence text. They do not identify RBPs, bacterial receptors, defense systems, or lysis pathways.")
        phage_tab, bacterium_tab, evidence_tab = st.tabs(["Phage features", "Bacterial features", "Evidence and limitations"])
        with phage_tab:
            st.dataframe([vars(feature) for feature in describe_records(phage_text, phage_type.lower())], use_container_width=True)
            st.caption("Protein molecular weight is estimated only when all residues are among the 20 standard amino acids. DNA GC% excludes ambiguous bases from its denominator.")
            st.divider()
            st.markdown("#### Stage 3B · Candidate coding regions (phage DNA only)")
            if phage_type == "DNA":
                candidates = screen_phage_orfs(phage_text, min_amino_acids=int(orf_minimum))
                st.write(f"**ORF candidates meeting threshold:** {len(candidates)}")
                if candidates:
                    st.dataframe([vars(candidate) for candidate in candidates], use_container_width=True)
                else:
                    st.info("No complete ATG-to-stop ORFs met the chosen length threshold.")
                st.caption("Six-frame screen with genetic code 1; 1-based inclusive coordinates per record. Minus-strand coordinates reference the submitted plus-strand record. Uses the first ATG within each stop-delimited frame segment. Not a validated gene caller; no receptor-binding function inferred.")
            else:
                st.info("ORF screening requires phage DNA. Protein sequences retain their descriptive feature analysis.")
        with bacterium_tab:
            st.dataframe([vars(feature) for feature in describe_records(bacteria_text, "dna")], use_container_width=True)
            st.caption("The input has not been verified as S. mutans by the software. Sequence-derived properties do not identify a bacterial receptor.")
        with evidence_tab:
            st.write("**Evidence level:** Descriptive computation from user-supplied FASTA only.")
            st.write("**Method:** Biopython SeqIO; GC calculation over canonical A/C/G/T; ProteinAnalysis for standard amino acid molecular weight.")
            st.write("**Stage 3B ORF method:** Six-frame complete ATG-to-stop screening, standard genetic code. ORF candidates are not confirmed genes or RBPs; alternative starts, partial and circular genes are not called.")
            st.write("**Not available:** Receptor-binding annotation, host specificity, defense-system annotation, lysis prediction, or calibrated confidence.")
        st.divider()
        st.subheader("Receptor and lysis dashboard · evidence assessment")
        a, b, c = st.columns(3)
        a.metric("Receptor", "Insufficient evidence")
        b.metric("Lysis pathway", "Insufficient evidence")
        c.metric("Confidence", "Not calibrated")
        st.markdown("**Input-specific conclusion:** No experimentally established receptor or molecular lysis pathway can be assigned to these uploaded sequences. Reference examples and exact sequence comparisons below do not validate infection.")
        st.caption("GC% is computed across unambiguous A/C/G/T positions only; these summaries do not establish infection or receptor binding.")
    else:
        st.warning("Correct the input errors above and analyze again.")


st.divider()
st.subheader("Stage 3C · Published reference evidence (not matched to your input)")
st.warning("This is a fixed educational reference catalog. Entries are NOT computed from or matched to the uploaded sequences. No ORF is assigned a function or receptor by this panel.")
reference_filter = st.selectbox("Browse reference examples", ["All", "Phage", "Bacterial", "Background"])
reference_rows = get_reference_records(reference_filter)
st.dataframe(reference_rows, use_container_width=True, hide_index=True)
st.caption("Evidence labels distinguish deposited sequences, database annotations, published computational annotations, laboratory proteomics/phenotypes, and secondary review context. Open the cited URL in each entry to inspect the source. Published observations do not imply results for the supplied genome or strain.")
with st.expander("What each evidence label means"):
    st.markdown("""
    - **DATABASE_SEQUENCE**: an accession-linked deposited/reference genome, not a functional test.
    - **DATABASE_PROTEIN_ANNOTATION**: a database gene/protein description; not experimentally tested as a phage receptor here.
    - **PUBLISHED_COMPUTATIONAL_ANNOTATION / PUBLISHED_GENOME_ANNOTATION**: function suggested by genome comparison or gene annotation.
    - **EXPERIMENTAL_PROTEOMICS / EXPERIMENTAL_HOST_PHENOTYPE**: experiments support the specific observation stated, not necessarily a direct protein–receptor pair.
    - **REVIEW_BACKGROUND**: secondary context, not experimental confirmation of an interaction.
    """)
st.info("Stage 4: published defense-system examples and optional imported DefenseFinder calls (below). The app does not run DefenseFinder or infer defenses directly from your uploaded DNA. Stage 5 offers conservative evidence assessment below, but not a calibrated host-range or lysis prediction. There is no calibrated confidence score.")


st.divider()
st.subheader("Stage 4 · Bacterial anti-phage defense evidence")
st.warning("Published findings below describe a specific research cohort, NOT your uploaded genome. No DefenseFinder analysis runs inside this website, and predicted defense systems do not prove resistance to any phage.")
st.markdown("#### A. Published computational genome-screening findings")
st.dataframe(published_defense_findings(), use_container_width=True, hide_index=True)
st.caption("Claisse et al. (2026), DOI 10.1111/omi.70014: 44 S. mutans clinical isolates analyzed using DefenseFinder/PADLOC. These are bioinformatics findings, not experimentally validated defense activity, and cannot be assigned to UA159 or a user upload.")
st.markdown("#### B. Optional: import authentic DefenseFinder systems output")
st.write("If you have **actually run DefenseFinder** on a documented genome/protein dataset, you may import its `defense_finder_systems.tsv`. Otherwise leave this section empty. No simulated calls are bundled as evidence.")
imported_file = st.file_uploader("Upload defense_finder_systems.tsv", type=["tsv", "txt"], key="defensefinder")
if imported_file is not None:
    df_accession = st.text_input("Genome accession or sample ID used for the DefenseFinder run (you provide this)", key="df_accession")
    df_version = st.text_input("DefenseFinder version used (you provide this)", key="df_version")
    df_source = st.text_input("Input type and provenance (e.g., ordered proteins from specified accession)", key="df_source")
    if st.button("Inspect imported predicted defense systems", key="df_import"):
        try:
            tsv = decode_defense_tsv(imported_file.getvalue())
            systems = parse_defensefinder_systems(tsv, df_accession, df_version, df_source)
        except DefenseImportError as exc:
            st.error(f"DefenseFinder import: {exc}")
        else:
            st.success(f"Parsed {len(systems)} predicted systems from the supplied TSV (not independently verified).")
            if systems:
                st.dataframe(systems, use_container_width=True, hide_index=True)
            else:
                st.info("The file has valid headers but no detected systems. This does not demonstrate biological absence of defense.")
            st.warning("These predictions are tied ONLY to the provenance you entered. This application has NOT matched the TSV to the bacterial FASTA you uploaded, nor confirmed activity experimentally.")
            st.download_button("Download imported prediction table (TSV)", data=(__import__('pandas').DataFrame(systems).to_csv(sep="\t", index=False) if systems else "No predicted systems in the imported table\n"), file_name="phageguide_imported_defense_predictions.tsv", mime="text/tab-separated-values")
st.markdown("**Evidence categories:** `PUBLISHED_COMPUTATIONAL_GENOME_SCREEN` = literature reports from genomics tools; `IMPORTED_DEFENSEFINDER_PREDICTION` = externally generated predictions, with user-supplied provenance; **experimentally confirmed defense activity in your uploaded strain: none shown**.")
st.caption("Official DefenseFinder output schema: https://github.com/mdmparis/defense-finder ; research source: https://doi.org/10.1111/omi.70014")


st.divider()
st.subheader("Stage 5 · Receptor and lysis evidence assessment")
st.warning("These are transparent evidence checks, NOT a validated host-range predictor. Uploaded bacterial FASTA is not independently assigned a strain identity. No percentage confidence scores are generated.")
assessment = assess_uploaded_pair()
st.markdown("#### A. Assessment of the uploaded phage–bacterium pair")
ac, bc = st.columns(2)
ac.metric("Molecular receptor", assessment.receptor_status)
bc.metric("Molecular lysis pathway", assessment.lysis_status)
st.write("**Reason:**", assessment.reason)
st.write("**Sequence-based evidence:**", assessment.sequence_evidence)
st.write("**Experiment-based evidence:**", assessment.experimental_evidence)
st.caption(assessment.evidence_limit)

st.markdown("#### B. Published strain-specific experimental examples (reference only)")
st.caption("These findings concern the named study strains. Choosing an example does NOT verify the identity or biological response of either FASTA uploaded above.")
example_choice = st.selectbox("Choose a published phage–strain example", ["M102AD / P42S", "M102AD / UA159"], key="stage5_case")
case = reference_scenario(example_choice)
st.write("**Evidence category:**", case['evidence_class'])
st.write("**Published observation:**", case['observation'])
st.write("**Receptor interpretation:**", case['receptor_result'])
st.write("**Lysis interpretation:**", case['lysis_result'])
st.write("**Source:**", case['source'])
st.write("**Publication:**", case['source_url'])
st.caption(case['limitation'])

st.markdown("#### C. Optional exact protein-sequence comparison")
st.write("Supply a **separate protein reference FASTA** and its provenance. The application compares it with phage protein input or translated Stage 3B ORFs. Identical amino-acid letters are only a sequence match, not proof of function, adsorption, susceptibility, or lysis. No reference protein FASTA is bundled as biological evidence.")
reference_file = st.file_uploader("Reference protein FASTA (optional)", type=["fasta", "faa", "fa", "txt"], key="stage5_reference")
reference_id = st.text_input("Reference accession / identifier (source-verified by you)", key="stage5_reference_id")
reference_provenance = st.text_input("Reference sequence source and version / provenance", key="stage5_reference_provenance")
if st.button("Compare phage sequences to supplied reference", key="stage5_compare"):
    if reference_file is None:
        st.error("Choose a reference protein FASTA file first.")
    elif not phage_text:
        st.error("Paste or upload a phage sequence in Step 1 first.")
    else:
        try:
            supplied_reference = decode_fasta_upload(reference_file.getvalue())
            rows = compare_reference_protein(phage_text, phage_type.lower(), supplied_reference,
                                             reference_id, reference_provenance, int(orf_minimum))
        except (ValueError, SequenceInputError) as exc:
            st.error(f"Comparison: {exc}")
        else:
            if rows:
                st.success(f"{len(rows)} exact protein-sequence matches found. No biological function or host interaction inferred.")
                st.dataframe(rows, hide_index=True, use_container_width=True)
            else:
                st.info("No exact protein-sequence matches found under these rules. This does not establish lack of homology or phage activity.")
            st.warning("Reference identifier/provenance was entered by the user, not independently verified; an exact match does NOT confirm receptor binding or lysis.")
