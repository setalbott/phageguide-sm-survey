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
from analysis.two_upload_reference import assess_reference_pair, fetch_ncbi_reference, PHAGE_ACCESSION, BACTERIUM_ACCESSION
from analysis.susceptibility_evidence import get_susceptibility_records, get_ksm96_serotype_totals, lookup_susceptibility, NO_RECORD, RECORDS

st.set_page_config(page_title="PhageGuide-SM | Bioinformatics Workbench", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

# Global style is additive, not a replacement for Streamlit's accessible controls.
from ui_theme import inject_theme, alignment_html
inject_theme()

SAMPLE_PHAGE = ">synthetic_forward_orf\nATGAAATAG"
SAMPLE_BACTERIA = ">synthetic_bacteria\nATGCATGC"
SAMPLE_REFERENCE = ">synthetic_reference_MT\nMT"

def load_sample():
    st.session_state["survey_phage_type"] = "DNA"
    st.session_state["survey_phage_source"] = "Paste FASTA"
    st.session_state["survey_bacteria_source"] = "Paste FASTA"
    st.session_state["survey_phage_fasta"] = SAMPLE_PHAGE
    st.session_state["survey_bacteria_fasta"] = SAMPLE_BACTERIA
    st.session_state["survey_orf_min"] = 2
    st.session_state["analysis_reference_source"] = "Paste FASTA"
    st.session_state["analysis_reference_fasta"] = SAMPLE_REFERENCE
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


@st.cache_data(ttl=86400, show_spinner=False)
def cached_reference(accession):
    return fetch_ncbi_reference(accession)


def render_two_upload_benchmark(phage_text, bacteria_text, phage_type):
    st.markdown("#### Two-upload benchmark · φAPCM01 / DPC6143")
    st.caption("**Reference match — not a prediction.** Requires two complete DNA genomes identical to the versioned published reference sequences. FASTA headers and lengths alone are not proof of identity.")
    if phage_type != "DNA":
        st.info("Phage protein input cannot verify the φAPCM01 DNA genome. Upload phage DNA for the two-genome benchmark; other analysis features remain available.")
        return
    references = {}
    for label, accession in (("phage", PHAGE_ACCESSION), ("bacterium", BACTERIUM_ACCESSION)):
        try:
            references[label] = cached_reference(accession)
        except Exception:
            references[label] = None
    outcome = assess_reference_pair(phage_text, bacteria_text, references["phage"], references["bacterium"])
    for key in ("phage", "bacterium"):
        item = outcome[key]
        summary = f"**{item.label} — {item.accession}**: {item.explanation}"
        if item.status == "matched":
            st.success(summary)
        elif item.status == "reference_unavailable":
            st.warning(summary)
        else:
            st.info(summary)
    if outcome["both_matched"]:
        st.success("**Verified reference match: φAPCM01 / DPC6143.** Both uploaded DNA genomes match the accession-specific published sequences exactly. The results below describe the documented laboratory pair, not a newly predicted outcome.")
        st.markdown("**Attachment and receptor evidence**")
        st.write("Wall & Wall (2025), Figure 1: fluorescently labeled φAPCM01 attached to DPC6143 and injected its DNA. Resistance mutations in rhamnose–glucose polysaccharide (RGP) biosynthesis genes impaired adsorption (Figures 3–4); RGP surface glycans are implicated in receptor function. This does not resolve every molecular contact or attachment protein.")
        st.markdown("[Wall & Wall (2025), Figure 1 and receptor experiments](https://doi.org/10.1099/mic.0.001575)")
        st.markdown("**Infection and growth observations**")
        st.write("Dalmasso et al. (2015): φAPCM01 infected DPC6143 in the reported host-range panel (only 1 of 17 strains). In broth and artificial saliva the phage inhibited DPC6143 growth (Figure 3); in preformed biofilms metabolic activity and viable counts declined (Figure 4). These are findings for the tested strains and conditions.")
        st.markdown("[Dalmasso et al. (2015), Figures 2–4](https://doi.org/10.1371/journal.pone.0138651)")
        st.warning("**Molecular lysis pathway: undetermined.** Putative endolysin/lysis annotations do not experimentally establish the exact pathway or timing in this specific pair.")
        with st.expander("Evidence interpretation and limitations"):
            st.write("Both exact reference sequences match the known published pair. This is an identity-based evidence retrieval workflow; it does not derive infection from novel sequence combinations. Attachment, entry, growth inhibition and lysis-mechanism identification are distinct endpoints. Wall & Wall also observed DNA entry into strain UA159 without productive susceptibility, so binding does not imply infection.")
    else:
        st.warning("**Insufficient verification for the two-genome benchmark.** The uploaded pair is not confirmed as the published φAPCM01 / DPC6143 reference combination; study-specific receptor, growth and lysis findings are not transferred to an unverified upload.")
    st.caption("Predicting susceptibility for unfamiliar phage–strain combinations remains an open research requirement.")


def render_results(phage_text, bacteria_text, phage_type, orf_minimum, reference_text=None):
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
    render_two_upload_benchmark(phage_text, bacteria_text, phage_type)
    st.markdown("#### Protein alignment")
    from analysis.auto_comparison import automatic_comparison
    try:
        report = automatic_comparison(phage_text, phage_type, int(orf_minimum), reference_text)
    except (ValueError, SequenceInputError) as exc:
        st.error(f"Reference comparison: {exc}")
    else:
        if report["status"] != "ok":
            st.info(report["message"])
        else:
            st.caption(f"{len(report['rows'])} aligned protein pair(s) · Sequence similarity only")
            for i, item in enumerate(report["rows"], start=1):
                st.markdown(f"**Alignment {i}: {item['query_id']} vs {item['reference_id']}**")
                m1, m2, m3 = st.columns(3)
                m1.metric("Matches", item["match_count"])
                m2.metric("Substitutions", item["substitution_count"])
                m3.metric("Identity", f"{item['percent_identity_over_alignment']:.1f}%")
                st.caption(f"Gaps: {item['gap_positions']} · Total differences: {item['difference_count']}")
                st.markdown(alignment_html(item), unsafe_allow_html=True)
                with st.expander("View alignment coordinates and scoring notes"):
                    st.code(format_alignment(item), language="text")
                    st.caption("Global protein alignment, not a verified functional match.")
    st.warning("**Susceptibility for an untested phage–strain pair remains undetermined.** Amino-acid identity does not measure binding or infection probability; validated host-range or binding data would be required.")
    with st.expander("Sequence details · Length, GC content and identifiers"):
        a, b = st.columns(2)
        with a:
            show_sequence_summary("Phage", summaries["Phage"])
        with b:
            show_sequence_summary("Bacterium", summaries["Bacterium"])

    features_tab, assessment_tab, comparison_tab = st.tabs(["Features", "Receptor / lysis", "Compare more"])
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


def pair_results_page():
    """Landing page for published, named phage–strain records only.

    Delisle's P42S propagation/adsorption experiment is presented separately
    from genome-derived lysis annotations; nothing here predicts a new pair.
    """
    st.markdown("""<div class="pg-brand">
      <div class="pg-brand-kicker">PHAGEGUIDE / GENOMICS WORKBENCH</div>
      <div class="pg-brand-title">Pair results. <em>Evidence first.</em></div>
      <div class="pg-brand-sub">Select a published phage–strain pair to inspect its reported findings.</div>
    </div>""", unsafe_allow_html=True)
    st.markdown("### Published findings for a selected pair")
    st.caption("These are study-specific results, not predictions or analyses of your uploaded FASTA.")

    available_phages = sorted({record.phage for record in RECORDS})
    # M102AD is the introductory example, even if the reference set changes.
    phage_choices = ["M102AD"] + [name for name in available_phages if name != "M102AD"]
    left, right = st.columns([0.95, 1.05], gap="large")
    with left:
        st.markdown('<div class="pg-section-kicker">01 / SELECT PAIR</div>', unsafe_allow_html=True)
        phage = st.selectbox("Phage", phage_choices, index=0, key="home_phage")
        strains = sorted({record.strain for record in RECORDS if record.phage == phage})
        if phage == "M102AD":
            strains = ["P42S"] + [name for name in strains if name != "P42S"]
        options = strains + ["Other strain (enter identifier)"]
        strain_choice = st.selectbox("S. mutans strain", options, index=0, key="home_strain_choice")
        if strain_choice == "Other strain (enter identifier)":
            strain = st.text_input("Strain identifier", placeholder="Enter an exact strain ID", key="home_custom_strain").strip()
        else:
            strain = strain_choice
        st.caption("Pairs are matched to curated records by exact identifier. Unlisted does not mean resistant.")
        st.info("For sequence uploads and automated protein alignments, use **Analyze** in the top navigation.")

    with right:
        st.markdown('<div class="pg-section-kicker">02 / PUBLISHED RESULTS</div>', unsafe_allow_html=True)
        if not strain:
            st.info(NO_RECORD)
            return
        matches = lookup_susceptibility(phage, strain)["records"]
        if not matches:
            st.info(NO_RECORD)
            st.caption("No assay outcome is available for this exact pair in the curated reference set. This is not a susceptibility prediction.")
            return

        record = matches[0]
        st.markdown(f"#### {phage} / {strain}")
        st.success("Published finding for the selected pair")
        if phage == "M102AD" and strain == "P42S":
            st.markdown("**What researchers observed**")
            st.write("Delisle et al. (2012) grew phage M102AD on *S. mutans* P42S (serotype c) and measured phage attachment (adsorption) to this established propagation host.")
            st.caption("Adsorption rate reported in Table 1: 1.5 × 10⁻¹⁰ mL/min. This is an attachment measurement, not a receptor identification or an identified lysis pathway.")
            st.markdown("**Receptor evidence**")
            st.info("**Exact receptor unresolved.** Adsorption/antiserum experiments implicated serotype-associated cell-wall carbohydrates, but did not identify and validate the precise P42S receptor molecule.")
            st.markdown("**Lysis pathway**")
            st.warning("**Undetermined for this pair.** The M102AD genome has *annotated* ORF18 (putative holin), ORF19 (putative endolysin), and ORF20 (putative endolysin). Those genome annotations are not an experimentally demonstrated P42S-specific lysis pathway.")
            with st.expander("Assay details, interpretation and evidence limits"):
                st.write("Propagation and plaque-based phage growth: P42S was used to propagate M102AD. Adsorption kinetics were measured separately. Neither experiment identifies a receptor molecule or demonstrates the roles of the annotated lysis proteins in this specific pair.")
                st.write("Genome annotations: Delisle et al. (2012), Table 3, lists ORF18 as holin, ORF19 as endolysin/Glyco 25/Acm, and ORF20 as endolysin/CHAP. The study describes these as predicted/putative protein functions.")
        else:
            st.markdown("**What researchers observed**")
            st.write(record["observed_result"])
            st.caption(f"Assay: {record['assay']} · Reported classification: {record['outcome']}")
            st.markdown("**Receptor evidence**")
            st.info("**Insufficient pair-specific evidence.** This record reports an assay outcome; it does not establish the exact bacterial receptor.")
            st.markdown("**Lysis pathway**")
            st.warning("**Undetermined for this pair.** A plaque or host-range result does not identify the molecular lysis pathway.")
            if record.get("review_flag"):
                st.warning(f"**Identifier needs review:** {record['review_flag']}")
            with st.expander("Study note and evidence limits"):
                st.write(record.get("notes") or "No additional verified strain-specific mechanism in this record.")
                st.write("The published outcome applies only to the named phage and strain under the reported assay conditions.")

        st.markdown(f"**Source:** [{record['source']}]({record['source_url']})")
        st.caption("**Evidence class:** Published experimental observation. Attachment, susceptibility, receptor identity, and lysis mechanism are separate endpoints.")


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
        st.caption("Two-upload benchmark: use full φAPCM01 phage DNA and DPC6143 bacterial chromosome FASTA. Download versioned records from NCBI: [KR153145.1](https://www.ncbi.nlm.nih.gov/nuccore/KR153145.1?report=fasta) and [NZ_CP172847.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP172847.1?report=fasta). The built-in synthetic sample is for ORF/alignment testing, not the reference match.")
        st.button("Load Sample", key="survey_load_sample", on_click=load_sample, use_container_width=True)
        st.markdown("**Phage · DNA or protein**")
        phage_type = st.radio("Sequence type", ["DNA", "Protein"], horizontal=True, key="survey_phage_type")
        phage_text = input_widget("Phage", "survey_phage_source", "survey_phage_fasta", "phage", ["fasta", "fa", "faa", "fna", "txt"], ">phage_example\nATGAAATAG")
        st.markdown("**S. mutans · bacterial DNA**")
        bacteria_text = input_widget("Bacterium", "survey_bacteria_source", "survey_bacteria_fasta", "bacteria", ["fasta", "fa", "fna", "txt"], ">bacterial_example\nATGCATGC")
        st.markdown("**Reference protein · optional**")
        st.caption("Paste or upload a protein FASTA to include an alignment in the report.")
        reference_text = input_widget("Reference protein", "analysis_reference_source", "analysis_reference_fasta", "analysis_reference_upload", ["fasta", "faa", "fa", "txt"], ">reference_protein\nMT")
        with st.expander("ORF options"):
            orf_minimum = st.number_input("Minimum translated protein length (aa)", min_value=1, max_value=10000, value=30, step=1, key="survey_orf_min", help="The built-in sample uses 2 amino acids.")
        if st.button("Analyze", type="primary", use_container_width=True, key="run_analysis"):
            st.session_state["analyzed_inputs"] = (phage_text, bacteria_text, phage_type, int(orf_minimum), reference_text)
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
    with st.expander("Related literature: M102 is a different phage"):
        st.write("Shibata, Yamashita & van der Ploeg (2009) investigated **M102**, not φAPCM01. M102 adsorption depended on the serotype-c rhamnose–glucose polysaccharide (RGP) glucose side chain: RGP-glucose-deficient mutants failed to adsorb M102. This mechanistic result must not be attributed directly to φAPCM01 or M102AD.")
        st.markdown("[Shibata et al. (2009), FEMS Microbiology Letters](https://doi.org/10.1111/j.1574-6968.2009.01546.x)")
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
    st.markdown("1. Go to **Analyze** and choose **Load Sample**.\n2. Click **Analyze** and look for ORF **MK** at bases **1–9**.\n3. The built-in sample automatically aligns translated **MK** against reference **MT** in the main report.\n4. Use **Published evidence** to view studies, receptor status, and predicted defenses.")
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
    st.Page(pair_results_page, title="Pair Results", icon="🔬", default=True),
    st.Page(analysis_page, title="Analyze", icon="🧬"),
    st.Page(evidence_page, title="Published evidence", icon="📚"),
    st.Page(help_page, title="Help", icon="❔"),
], position="top")
nav.run()
