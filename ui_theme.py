"""PhageGuide-SM presentation only: navy/gold responsive genomics styling.

No analytical results or scientific evidence are generated in this module.
"""
import html


def alignment_html(result, width=48):
    """Accessible, color-coded HTML alignment; position ticks are 1-based.

    Match, substitution, and gap information always comes from the original
    Biopython-derived alignment result, never from presentation logic.
    """
    qa = result['query_aligned']
    ra = result['reference_aligned']
    symbols = result['match_symbols']
    qpos = rpos = 0
    panels = []
    for start in range(0, len(qa), width):
        qchunk, rchunk, marks = qa[start:start + width], ra[start:start + width], symbols[start:start + width]
        qstart, rstart = qpos + 1, rpos + 1
        qpos += len(qchunk.replace('-', ''))
        rpos += len(rchunk.replace('-', ''))
        def colored(sequence):
            out = []
            for letter, mark in zip(sequence, marks):
                cls = 'match' if mark == '|' else ('gap' if letter == '-' else 'difference')
                out.append(f'<span class="pg-aa-{cls}">{html.escape(letter)}</span>')
            return ''.join(out)
        rendered_marks = ''.join('<span class="pg-aa-' + ('match' if m == '|' else ('gap' if m == ' ' else 'difference')) + '">' + ('|' if m == '|' else ('·' if m == '.' else ' ')) + '</span>' for m in marks)
        panels.append('<div class="pg-alignment-block" role="group" aria-label="Aligned protein sequences">'
                      f'<div><span class="pg-alignment-label">INPUT {qstart:>4}</span> {colored(qchunk)} <span class="pg-position">{qpos}</span></div>'
                      f'<div><span class="pg-alignment-label">           </span> {rendered_marks}</div>'
                      f'<div><span class="pg-alignment-label">REF   {rstart:>4}</span> {colored(rchunk)} <span class="pg-position">{rpos}</span></div>'
                      '</div>')
    return '<div class="pg-alignment" aria-label="Position-by-position amino acid alignment">' + ''.join(panels) + '</div>'


def inject_theme():
    import streamlit as st
    st.markdown("""
<style>
:root{--pg-navy:#091626;--pg-mid:#10253a;--pg-gold:#e8bd64;--pg-gold-soft:#f8d995;--pg-text:#eff5fb;--pg-muted:#aec4d5}
.stApp{background:radial-gradient(ellipse at 83% 0%,rgba(64,98,127,.19),transparent 44%),linear-gradient(160deg,#091626 0%,#0b1c2f 63%,#0b1729 100%);color:var(--pg-text)}
[data-testid="stHeader"]{background:rgba(9,22,38,.96)}
[data-testid="stToolbar"]{color:#f6e3b5}
.block-container{padding-top:1rem;max-width:1600px;padding-bottom:3.5rem}
[data-testid="stMain"]{min-width:0}
.pg-brand{position:relative;isolation:isolate;overflow:hidden;border:1px solid rgba(232,189,100,.4);border-radius:18px;background:linear-gradient(114deg,#0a192c,#112b42 70%,#152a3f);padding:20px 25px;margin-bottom:13px;box-shadow:0 14px 35px rgba(0,0,0,.16)}
.pg-brand:after{content:"";position:absolute;right:-4px;top:-18px;width:320px;height:170px;opacity:.27;pointer-events:none;z-index:-1;background-image:repeating-linear-gradient(150deg,transparent 0 15px,rgba(232,189,100,.45) 16px 17px,transparent 18px 31px),radial-gradient(ellipse at 25% 15%,#deb86c 0 3px,transparent 4px),radial-gradient(ellipse at 75% 75%,#deb86c 0 3px,transparent 4px);background-size:100% 100%,28px 30px,28px 30px;transform:skewX(-15deg)}
.pg-brand-kicker{color:var(--pg-gold);font-weight:750;letter-spacing:.22em;font-size:.68rem}
.pg-brand-title{font-size:clamp(1.55rem,2.9vw,2.4rem);font-weight:800;letter-spacing:-.035em;color:#f6f7fa;margin:.25rem 0}
.pg-brand-title em{font-style:normal;color:var(--pg-gold-soft)}
.pg-brand-sub{font-size:.87rem;color:#c5d8e7;margin-top:4px}
.pg-section-kicker{color:var(--pg-gold);font-weight:800;font-size:.7rem;letter-spacing:.14em;margin-bottom:.1rem}
.pg-empty{min-height:305px;padding:2rem;border:1px dashed rgba(232,189,100,.45);border-radius:14px;background:linear-gradient(135deg,rgba(27,56,79,.38),rgba(16,37,58,.4)),repeating-linear-gradient(0deg,transparent 0,transparent 31px,rgba(232,189,100,.055) 32px),repeating-linear-gradient(90deg,transparent 0,transparent 31px,rgba(232,189,100,.055) 32px);display:flex;flex-direction:column;justify-content:center;align-items:center;text-align:center;gap:.65rem}
.pg-empty-icon{font-size:2.3rem;color:var(--pg-gold-soft)}.pg-empty strong{letter-spacing:.12em;color:#f7e2b0}.pg-empty p{max-width:340px;color:#bad0df;margin:0}.pg-empty span{font-size:.7rem;color:#e5c07b;letter-spacing:.1em}
.pg-status{border-left:3px solid var(--pg-gold);background:rgba(232,189,100,.09);padding:.72rem .95rem;border-radius:0 9px 9px 0;color:#f8e4ba;font-size:.91rem;margin:9px 0}
.pg-alignment{overflow-x:auto;border:1px solid rgba(232,189,100,.25);border-radius:12px;background:#081422;padding:12px;max-width:100%;}
.pg-alignment-block{white-space:pre;min-width:max-content;font:600 .88rem/1.7 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;margin-bottom:8px}
.pg-alignment-label{color:#9eb8ca}.pg-position{color:#8198aa;font-weight:400}
.pg-aa-match{color:#83e0bd;background:rgba(53,139,107,.10)}.pg-aa-difference{color:#f3bd73;background:rgba(217,141,57,.18)}.pg-aa-gap{color:#f08585;background:rgba(211,86,86,.13)}
.pg-legend{display:flex;flex-wrap:wrap;gap:12px;font-size:.75rem;color:#b9ccdb;margin-bottom:7px}
.pg-legend b{font-family:ui-monospace,Menlo,monospace}
h1,h2,h3,h4{letter-spacing:-.02em}h1{color:#f6f7fb!important}h3{font-size:1.12rem!important;color:#f7e5b8!important}
[data-testid="stMetric"]{background:rgba(20,49,72,.75);border:1px solid rgba(232,189,100,.22);padding:11px 13px;border-radius:12px}
[data-testid="stMetricLabel"]{color:#c7d7e5}[data-testid="stMetricValue"]{color:#fff2ce}
.stButton>button[kind="primary"]{background:#e8bd64;color:#142033;border:1px solid #f3d89e;font-weight:800}
.stButton>button[kind="primary"]:hover{background:#ffdb90;color:#071321}
[data-testid="stExpander"]{border-color:rgba(232,189,100,.25)}
[data-testid="stMarkdownContainer"] a{color:#f3ca78}
[data-testid="stTabs"] [data-baseweb="tab"]{font-weight:650}
[data-testid="stTextArea"] textarea, [data-testid="stCodeBlock"] code, [data-testid="stCode"]{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace!important}
/* Native Streamlit navigation remains keyboard-accessible and responsive. */
@media(max-width:760px){.block-container{padding:.75rem .7rem 3rem}.pg-brand{padding:15px 15px}.pg-brand:after{width:160px;opacity:.14}.pg-brand-sub{font-size:.8rem}.pg-empty{min-height:140px;padding:1.1rem}.pg-alignment-block{font-size:.76rem}[data-testid="stMetric"]{padding:8px}}
</style>
""", unsafe_allow_html=True)
