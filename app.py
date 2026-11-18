import io
import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from crispr_designer import analyze_sequence, parse_fasta, validate_dna

st.set_page_config(page_title="CRISPR gRNA Designer", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.block-container {max-width: 1180px; padding-top: 2rem;}
[data-testid="stMetric"] {background:#f5f8fb; border:1px solid #e4ebf2; padding:14px; border-radius:12px;}
h1, h2, h3 {color:#18324b;}
.small-muted {color:#64748b; font-size:0.9rem;}
</style>
""", unsafe_allow_html=True)

st.title("CRISPR gRNA Designer")
st.caption("A sequence-based SpCas9 guide RNA candidate finder • Python + Biopython")
st.info("Educational prototype: candidate ranking is based on simple sequence heuristics. It does not predict genome-wide off-targets or guarantee editing efficiency.")

with st.sidebar:
    st.header("Design settings")
    pam = st.selectbox("PAM pattern", ["NGG"], index=0, help="NGG is the canonical PAM for SpCas9.")
    guide_length = st.number_input("Guide length (nt)", min_value=18, max_value=22, value=20, step=1)
    min_gc = st.slider("Minimum GC (%)", 0, 100, 40)
    max_gc = st.slider("Maximum GC (%)", 0, 100, 60)
    exclude_poly_t = st.checkbox("Exclude guides with 4 or more consecutive T bases", value=True)
    st.markdown("---")
    st.caption("SpCas9 default: 20-nt protospacer followed by NGG PAM.")

tab1, tab2, tab3 = st.tabs(["01  Sequence input", "02  Candidate guides", "03  How it works"])

with tab1:
    st.subheader("Enter a target DNA sequence")
    uploaded = st.file_uploader("Upload a FASTA file", type=["fa", "fasta", "fna"])
    default_seq = "ATGGCTAGCTACGATCGGATCCGATGCTAGCTAGGCTAACCGGTTAACCGGATCGATCGGCTAACGATCGATCGG"
    sequence_text = st.text_area("Or paste DNA sequence (5′ → 3′)", value=default_seq, height=170,
                                 help="Whitespace and FASTA header lines are removed. Use A, C, G and T only.")
    analyze = st.button("Find candidate gRNAs", type="primary", use_container_width=True)

    if uploaded is not None:
        raw = uploaded.getvalue().decode("utf-8", errors="replace")
        try:
            seq = parse_fasta(raw)
            st.success(f"Loaded FASTA sequence: {len(seq):,} bp")
        except ValueError as exc:
            st.error(str(exc))
            seq = None
    else:
        seq = None

    if analyze:
        try:
            if uploaded is None:
                seq = validate_dna(sequence_text)
            if len(seq) < int(guide_length) + 3:
                raise ValueError(f"Sequence must be at least {int(guide_length) + 3} bases long.")
            results, stats = analyze_sequence(
                seq, guide_length=int(guide_length), min_gc=min_gc, max_gc=max_gc,
                exclude_poly_t=exclude_poly_t, pam_pattern=pam
            )
            st.session_state["results"] = results
            st.session_state["stats"] = stats
            st.session_state["sequence"] = seq
            st.success("Analysis completed.")
        except ValueError as exc:
            st.error(str(exc))

with tab2:
    if "results" not in st.session_state:
        st.subheader("Your results will appear here")
        st.write("Paste a sequence or upload a FASTA file, then select **Find candidate gRNAs**.")
    else:
        results = st.session_state["results"]
        stats = st.session_state["stats"]
        seq = st.session_state["sequence"]
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Sequence length", f"{stats['length']:,} bp")
        c2.metric("GC content", f"{stats['sequence_gc']:.1f}%")
        c3.metric("PAM sites scanned", stats["pam_sites"])
        c4.metric("Candidates passing filters", len(results))
        st.caption("Coordinates are 1-based and inclusive on the submitted sequence. Reverse-strand guides are reported in their 5′ → 3′ guide orientation.")
        if results:
            df = pd.DataFrame(results)
            st.dataframe(df, use_container_width=True, hide_index=True)
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("Download candidate guides (CSV)", csv, "crispr_grna_candidates.csv", "text/csv")
            best = results[0]
            st.markdown("#### Top-ranked sequence-based candidate")
            st.code(best["Guide (5'→3')"], language=None)
            st.write(f"Position: {best['Start']}-{best['End']} • Strand: {best['Strand']} • PAM: {best['PAM']} • GC: {best['GC %']:.1f}%")
        else:
            st.warning("No candidates passed the selected filters. Try widening the GC range or checking the input sequence.")
        st.markdown("**Scoring note:** Candidates are sorted by closeness to 50% GC. This is a transparent teaching heuristic, not a validated activity or specificity score.")

with tab3:
    st.subheader("Workflow")
    st.markdown("""
1. **Input:** Read a DNA sequence pasted by the user or uploaded in FASTA format.
2. **Validation:** Normalize to uppercase and check that the sequence contains only A, C, G and T.
3. **PAM scan:** Search both strands for SpCas9 `NGG` PAM sites.
4. **Guide extraction:** Extract the adjacent protospacer (20 nt by default), accounting for strand orientation.
5. **Filtering:** Apply GC-content limits and optionally exclude four-or-more consecutive T bases.
6. **Ranking:** Sort remaining candidates by closeness to 50% GC and export a CSV report.
    """)
    st.warning("A real experimental design pipeline should additionally use a genome reference, align candidate guides to the relevant organism, calculate off-target risk, and validate candidates experimentally.")
