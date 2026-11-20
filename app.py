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
