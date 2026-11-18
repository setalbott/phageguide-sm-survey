"""PhageGuide-SM survey prototype — science-inspired responsive UI.

The analysis functions and evidence data remain in their original modules.
Only the presentation and persistent display state are reorganized here.
"""
import streamlit as st
from analysis.sequence_utils import SequenceInputError, summarize_fasta, decode_fasta_upload
from analysis.feature_analysis import describe_records
from analysis.orf_analysis import screen_phage_orfs
from analysis.reference_evidence import get_reference_records
from analysis.defense_evidence import published_defense_findings, parse_defensefinder_systems, decode_defense_tsv, DefenseImportError
from analysis.interaction_assessment import assess_uploaded_pair, reference_scenario
from analysis.sequence_comparison import compare_protein_sequences, format_alignment
from analysis.susceptibility_evidence import get_susceptibility_records, get_ksm96_serotype_totals, lookup_susceptibility, NO_RECORD, RECORDS

st.set_page_config(page_title="PhageGuide-SM | Bioinformatics Workbench", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

# Global style is additive, not a replacement for Streamlit's accessible controls.
from ui_theme import inject_theme, alignment_html
inject_theme()

SAMPLE_PHAGE = ">synthetic_forward_orf\nATGAAATAG"
SAMPLE_BACTERIA = ">synthetic_bacteria\nATGCATGC"

def load_sample():
    st.session_state["survey_phage_type"] = "DNA"
    st.session_state["survey_phage_source"] = "Paste FASTA"
    st.session_state["survey_bacteria_source"] = "Paste FASTA"
    st.session_state["survey_phage_fasta"] = SAMPLE_PHAGE
    st.session_state["survey_bacteria_fasta"] = SAMPLE_BACTERIA
    st.session_state["survey_orf_min"] = 2
    st.session_state.pop("analyzed_inputs", None)  # Never show stale results after loading.


def input_widget(label, source_key, text_key, uploader_key, allowed_types, example):
    method = st.radio(f"{label} input", ["Paste FASTA", "Upload FASTA"], horizontal=True, key=source_key, label_visibility="collapsed")
    if method == "Paste FASTA":
        return st.text_area("FASTA sequence", height=145, placeholder=example, key=text_key, label_visibility="collapsed")
    uploaded = st.file_uploader("Upload FASTA", type=allowed_types, key=uploader_key)
    if uploaded is None:
        return None
    try:
        return decode_fasta_upload(uploaded.getvalue())
    except SequenceInputError as exc:
        st.error(f"{label}: {exc}")
        return None


def show_sequence_summary(label, summary):
    st.markdown(f"**{label}**")
    st.metric("Sequence length", f"{summary.total_length:,}")
    st.caption(f"{summary.record_count} record(s) · {summary.sequence_type.upper()} · {', '.join(summary.identifiers[:3])}")
    if summary.sequence_type == "dna":
        st.metric("GC content", f"{summary.gc_percent:.2f}%" if summary.gc_percent is not None else "N/A")
        if summary.ambiguous_bases:
            st.caption(f"{summary.ambiguous_bases} ambiguous nucleotide positions")
    if summary.type_overlap_warning:
        st.caption("Letters alone cannot verify whether this is DNA or protein.")


def render_results(phage_text, bacteria_text, phage_type, orf_minimum):
    summaries = {}
    for label, text, kind in (("Phage", phage_text, phage_type.lower()), ("Bacterium", bacteria_text, "dna")):
        try:
            summaries[label] = summarize_fasta(text, kind)
        except SequenceInputError as exc:
            st.error(f"{label}: {exc}")
    if len(summaries) != 2:
        st.warning("Correct the FASTA inputs and click Analyze again.")
        return

    st.success("Validated · Both FASTA inputs passed")
    a, b = st.columns(2)
    with a:
        show_sequence_summary("Phage", summaries["Phage"])
    with b:
        show_sequence_summary("Bacterium", summaries["Bacterium"])

    features_tab, assessment_tab, comparison_tab = st.tabs(["Features", "Receptor / lysis", "Compare"])
    with features_tab:
        st.markdown("**Phage sequence features**")
        st.dataframe([vars(f) for f in describe_records(phage_text, phage_type.lower())], use_container_width=True, hide_index=True)
        if phage_type == "DNA":
            candidates = screen_phage_orfs(phage_text, min_amino_acids=int(orf_minimum))
            st.markdown(f"**Candidate ORFs · {len(candidates)}**")
            if candidates:
                st.dataframe([vars(c) for c in candidates], use_container_width=True, hide_index=True)
            else:
                st.info("No complete ATG-to-stop ORFs meet this minimum length.")
        else:
            st.info("ORF screening requires phage DNA. Protein properties are still available above.")
        with st.expander("Bacterial features and sequence-method limitations"):
            st.dataframe([vars(f) for f in describe_records(bacteria_text, "dna")], use_container_width=True, hide_index=True)
            st.caption("GC% excludes ambiguous bases. Six-frame ORFs are candidates, not confirmed genes. DNA uploaded as S. mutans is not strain-verified.")

    with assessment_tab:
        render_assessment()
    with comparison_tab:
        render_comparison(phage_text, phage_type, orf_minimum)


def render_assessment():
    assessment = assess_uploaded_pair()
    st.caption("Uploaded pair · Sequence-level evidence only")
    st.markdown(f'<div class="pg-status"><strong>Receptor:</strong> {assessment.receptor_status}<br><strong>Lysis pathway:</strong> {assessment.lysis_status}<br><small>Missing verified strain identity, binding data, and infection/lysis measurements.</small></div>', unsafe_allow_html=True)
    with st.expander("Why this is insufficient evidence"):
        st.write(assessment.reason)
        st.write("**Sequence evidence:**", assessment.sequence_evidence)
        st.write("**Experimental evidence:**", assessment.experimental_evidence)
        st.write("GC content and amino-acid identity are not binding probabilities. A probability would require validated binding measurements and an independently evaluated calibrated model.")


def render_comparison(phage_text, phage_type, orf_minimum):
    st.caption("Align translated ORFs or uploaded phage proteins to a reference. Similarity does not prove binding.")
    comparison_mode = st.radio("Reference sequence source", ["Built-in synthetic example", "Upload reference protein FASTA"], horizontal=True, key="comparison_mode")
    reference_text = None
    if comparison_mode == "Built-in synthetic example":
        st.caption("Try MK vs MK (100%) or MK vs MT (50%). Synthetic only.")
        synthetic_choice = st.selectbox("Synthetic reference to compare to phage ORF MK", ["MK — identical (2 matches)", "MT — one different amino acid (1 match)"], key="reference_sample_choice")
        reference_seq = "MK" if synthetic_choice.startswith("MK") else "MT"
        reference_text = f">synthetic_reference_{reference_seq}\n{reference_seq}"
        reference_id = "SYNTHETIC-DEMO-ONLY"
        reference_provenance = "Built-in artificial example; no real protein or experiment"
        st.code(reference_text, language="text")
    else:
        reference_file = st.file_uploader("Reference protein FASTA (optional)", type=["fasta", "faa", "fa", "txt"], key="stage5_reference")
        reference_id = st.text_input("Reference accession / identifier (entered by you; not independently checked)", key="stage5_reference_id")
        reference_provenance = st.text_input("Reference sequence source / provenance", key="stage5_reference_provenance")
        if reference_file is not None:
            try:
                reference_text = decode_fasta_upload(reference_file.getvalue())
            except SequenceInputError as exc:
                st.error(f"Reference: {exc}")

    if st.button("Show aligned protein sequence comparison", key="stage5_compare_alignment"):
        if not reference_text:
            st.error("Provide a reference protein FASTA or choose the synthetic example.")
        elif not phage_text:
            st.error("Supply phage DNA/protein FASTA in Step 1 first. The guided sample works here.")
        else:
            try:
                comparison_rows = compare_protein_sequences(
                    phage_text, phage_type.lower(), reference_text,
                    reference_id, reference_provenance, int(orf_minimum))
            except (ValueError, SequenceInputError) as exc:
                st.error(f"Sequence comparison: {exc}")
            else:
                if not comparison_rows:
                    st.info("No qualifying translated ORFs or protein records to compare. For the synthetic DNA test, set minimum ORF length to 2.")
                else:
                    st.success(f"Compared {len(comparison_rows)} protein pair(s). These are sequence comparisons, not biological interaction predictions.")
                    for i, item in enumerate(comparison_rows, start=1):
                        st.markdown(f"**Comparison {i}: {item['query_id']} vs {item['reference_id']}**")
                        c1, c2, c3 = st.columns(3)
                        c1.metric("Matches", item['match_count'])
                        c2.metric("Differences", item['difference_count'])
                        c3.metric("Sequence identity", f"{item['percent_identity_over_alignment']:.1f}%")
                        st.markdown('<div class="pg-legend"><span><b style="color:#83e0bd">|</b> match</span><span><b style="color:#f3bd73">·</b> substitution</span><span><b style="color:#f08585">–</b> gap</span></div>', unsafe_allow_html=True)
                        st.markdown(alignment_html(item), unsafe_allow_html=True)
                        with st.expander("Plain-text alignment with exact coordinates"):
                            st.code(format_alignment(item), language="text")
                        st.caption(f"{item['substitution_count']} substitutions; {item['gap_positions']} gap positions. This is a global amino-acid alignment under fixed scoring rules, not a validated homology threshold. Source: {item['reference_provenance_user_supplied']}. No functional inference.")
                    st.warning("**Interpretation for your own upload:** GC% describes DNA composition; matches and differences compare amino-acid sequences. Neither is a measure of binding probability. Experimental binding/host-range data plus a validated calibrated predictor would be required for that estimate.")


def analysis_page():
    st.markdown("""<div class="pg-brand">
      <div class="pg-brand-kicker">PHAGEGUIDE / GENOMICS WORKBENCH</div>
      <div class="pg-brand-title">Sequence intelligence. <em>Evidence first.</em></div>
      <div class="pg-brand-sub">Compare phage and bacterial FASTA · Explore features · Review published studies</div>
    </div>""", unsafe_allow_html=True)
    left, right = st.columns([0.95, 1.05], gap="large")
    with left:
        st.markdown('<div class="pg-section-kicker">01 / INPUTS</div>', unsafe_allow_html=True)
        st.markdown("### Analyze sequences")
        st.button("Load Sample", key="survey_load_sample", on_click=load_sample, use_container_width=True)
        st.markdown("**Phage · DNA or protein**")
        phage_type = st.radio("Sequence type", ["DNA", "Protein"], horizontal=True, key="survey_phage_type")
        phage_text = input_widget("Phage", "survey_phage_source", "survey_phage_fasta", "phage", ["fasta", "fa", "faa", "fna", "txt"], ">phage_example\nATGAAATAG")
        st.markdown("**S. mutans · bacterial DNA**")
        bacteria_text = input_widget("Bacterium", "survey_bacteria_source", "survey_bacteria_fasta", "bacteria", ["fasta", "fa", "fna", "txt"], ">bacterial_example\nATGCATGC")
        with st.expander("ORF options"):
            orf_minimum = st.number_input("Minimum translated protein length (aa)", min_value=1, max_value=10000, value=30, step=1, key="survey_orf_min", help="The built-in sample uses 2 amino acids.")
        if st.button("Analyze", type="primary", use_container_width=True, key="run_analysis"):
            st.session_state["analyzed_inputs"] = (phage_text, bacteria_text, phage_type, int(orf_minimum))
        st.caption("Use public or synthetic data only. No confidential or patient sequences.")
    with right:
        st.markdown('<div class="pg-section-kicker">02 / RESULTS</div>', unsafe_allow_html=True)
        st.markdown("### Analysis results")
        if "analyzed_inputs" not in st.session_state:
            st.markdown('<div class="pg-empty"><div class="pg-empty-icon">⌬</div><strong>READY FOR ANALYSIS</strong><p>Enter two FASTA sequences or choose Load Sample. Results appear here after Analyze.</p><span>INPUT → FEATURES → EVIDENCE</span></div>', unsafe_allow_html=True)
        else:
            render_results(*st.session_state["analyzed_inputs"])
    st.caption("Sequence characteristics and published observations are not validated binding predictions.")


def evidence_page():
    st.markdown('<div class="pg-kicker">REFERENCE LIBRARY · VERIFIED PROVENANCE</div>', unsafe_allow_html=True)
    st.title("Published evidence")
    st.caption("Studies, database annotations, and computational defense screens are classified separately from uploaded data.")
    ref_tab, host_tab, susceptibility_tab, defense_tab = st.tabs(["Reference catalog", "Phage–strain studies", "Susceptibility Results", "Defense systems"])
    with ref_tab:
        st.caption("Published examples only · Not matched to your upload · No function assigned to uploaded ORFs")
        reference_filter = st.selectbox("Browse reference examples", ["All", "Phage", "Bacterial", "Background"])
        reference_rows = get_reference_records(reference_filter)
        st.dataframe(reference_rows, use_container_width=True, hide_index=True)
        st.caption("See source links in each record. Published results are not input-specific predictions.")
        with st.expander("What each evidence label means"):
            st.markdown("""
            - **DATABASE_SEQUENCE**: an accession-linked deposited/reference genome, not a functional test.
            - **DATABASE_PROTEIN_ANNOTATION**: a database gene/protein description; not experimentally tested as a phage receptor here.
            - **PUBLISHED_COMPUTATIONAL_ANNOTATION / PUBLISHED_GENOME_ANNOTATION**: function suggested by genome comparison or gene annotation.
            - **EXPERIMENTAL_PROTEOMICS / EXPERIMENTAL_HOST_PHENOTYPE**: experiments support the specific observation stated, not necessarily a direct protein–receptor pair.
            - **REVIEW_BACKGROUND**: secondary context, not experimental confirmation of an interaction.
            """)
    with host_tab:
        st.caption("Study-specific observations, not predictions for the uploaded sample.")
        st.markdown("**Published phage–strain examples**")
        st.caption("These findings concern the named study strains. Choosing an example does NOT verify the identity or biological response of either FASTA uploaded above.")
        example_choice = st.selectbox("Choose a published phage–strain example", ["M102AD / P42S", "M102AD / UA159"], key="stage5_case")
        case = reference_scenario(example_choice)
        st.write("**Phage:**", case['phage'])
        st.write("**S. mutans study strain:**", case['bacterial_strain'])
        st.write("**Study method:**", case['experiment_type'])
        st.write("**Evidence category:**", case['evidence_class'])
        st.write("**Published observation:**", case['observation'])
        st.write("**Receptor status:**", case['receptor_evidence_status'])
        st.write("**Receptor candidate / interpretation:**", case['receptor_candidate'])
        st.write("**Receptor interpretation:**", case['receptor_result'])
        st.write("**Still missing:**", case['missing_evidence'])
        st.write("**Lysis interpretation:**", case['lysis_result'])
        st.write("**Source:**", case['source'])
        st.write("**Publication:**", case['source_url'])
        st.caption(case['limitation'])
    with susceptibility_tab:
        st.markdown("**Experimentally tested susceptibility — named phage–strain pairs**")
        st.caption("Published plaque or propagation observations only. Nothing on this page identifies your uploaded bacterial strain or predicts an untested pair.")
        phage_filter = st.selectbox("Filter by phage", ["All", "SMHBZ8", "M102AD"], key="susceptibility_phage")
        strain_filter = st.text_input("Filter by exact strain identifier (or leave blank for all)", key="susceptibility_strain", placeholder="e.g. UA159, LM7, ES1, F")
        if strain_filter.strip():
            result = lookup_susceptibility(phage_filter, strain_filter.strip())
            rows = result["records"]
        else:
            rows = get_susceptibility_records(phage_filter, "All")
        if rows:
            st.dataframe(rows, use_container_width=True, hide_index=True,
                         column_order=["phage", "strain", "assay", "observed_result", "outcome", "review_flag", "source", "source_url", "notes"])
            flags = [r for r in rows if r["review_flag"]]
            if flags:
                st.warning("Identifier review required: Table 3 labels one SMHBZ8 strain M76219, while the Methods call it MT6219. Do not automatically merge these identifiers or attach a genome accession.")
        else:
            st.info(NO_RECORD)
        with st.expander("What these susceptibility labels mean"):
            st.write("SMHBZ8: sensitive/resistant refers to plaque presence in Ben-Zaken et al. (2021), Table 3. M102AD: host-range results and a separately established propagation host from Delisle et al. (2012). Neither reports a validated susceptibility probability for uploaded sequences.")
            st.write("In the SMHBZ8 paper, its Table 3 supports 13 S and 4 R S. mutans records, despite a contradictory sentence in the surrounding narrative. M102AD's 25-strain additional panel contains one sensitive strain (F); P42S is listed separately as its established propagation host.")
        st.divider()
        st.markdown("**φKSM96: study totals by bacterial serotype**")
        st.dataframe(get_ksm96_serotype_totals(), use_container_width=True, hide_index=True)
        st.caption("113/125 susceptible overall, from Sugai et al. (2023), Table 1. These are aggregated counts, **not** individual isolate lookup records; the exact isolate-level labels/genomes remain unverified here.")
        with st.expander("Why not make an untested-pair prediction?"):
            st.write("This is a literature lookup. A trained predictor needs reliably linked phage and strain genomes, consistent experimental endpoint labels, and independently held-out testing data. Attachment, plaque susceptibility, growth inhibition, and receptor identification are not interchangeable.")
    with defense_tab:
        st.markdown("**Anti-phage defense evidence**")
        st.caption("Published predicted defense systems are NOT assigned to the uploaded genome.")
        st.markdown("**Published genome screens**")
        st.dataframe(published_defense_findings(), use_container_width=True, hide_index=True)
        st.caption("Claisse et al. (2026), DOI 10.1111/omi.70014: 44 S. mutans clinical isolates analyzed using DefenseFinder/PADLOC. These are bioinformatics findings, not experimentally validated defense activity, and cannot be assigned to UA159 or a user upload.")
        st.markdown("**Import a DefenseFinder result (optional)**")
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


def help_page():
    st.markdown('<div class="pg-kicker">LAB NOTES · START HERE</div>', unsafe_allow_html=True)
    st.title("How to use PhageGuide-SM")
    st.markdown("1. Go to **Analyze** and choose **Load Sample**.\n2. Click **Analyze** and look for ORF **MK** at bases **1–9**.\n3. In **Compare**, try **MK vs MK** and **MK vs MT**.\n4. Use **Published evidence** to view studies, receptor status, and predicted defenses.")
    st.info("**Insufficient evidence** means the upload lacks validated binding, infection, or lysis measurements. Sequence identity and GC% do not estimate binding probability.")
    with st.expander("What do matches and GC content mean?"):
        st.write("GC% is the fraction of unambiguous DNA bases that are G or C. Alignment matches count identical amino acids; mismatches and gaps are different positions. These measurements do not demonstrate a receptor interaction.")
    with st.expander("Scientific limitations and methods"):
        st.write("This software does not establish the identity of an uploaded S. mutans strain, run DefenseFinder on submitted DNA, or verify molecular receptors experimentally.")
        st.write("The ORF screen searches both strands for complete ATG-to-stop candidates under genetic code 1; it is not a validated gene caller.")
        st.write("The sequence comparator uses a global amino-acid alignment; no binding probability or calibrated confidence is produced.")
    with st.expander("Research context and participant privacy"):
        st.write("BIOL 831 research prototype. Reference studies are available on the Published evidence page. Do not submit patient, unpublished, or confidential data to the public demonstration website.")

# Native top navigation: supported by Streamlit >=1.36.
nav = st.navigation([
    st.Page(analysis_page, title="Analyze", icon="🧬", default=True),
    st.Page(evidence_page, title="Published evidence", icon="📚"),
    st.Page(help_page, title="Help", icon="❔"),
], position="top")
nav.run()
