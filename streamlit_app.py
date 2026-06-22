"""Gaze — SPPS impurity attribution for LC-MS. Research use only."""

import os, sys, tempfile
from datetime import datetime, timezone

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

st.set_page_config(page_title="GAZE", layout="wide")


st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@300;400;500;600&display=swap');


* {
  font-family: 'IBM Plex Mono', monospace !important;
}

[data-testid="stIconMaterial"],
[data-testid="stTooltipIcon"] {
  font-family: 'Material Symbols Rounded', 'Material Symbols Outlined' !important;
  text-shadow: none !important;
  -webkit-text-stroke: 0 !important;
}

html, body { background: #0b0f1e !important; }
.stApp { background: transparent !important; }

.stApp *:not([data-testid="stIconMaterial"]):not([data-testid="stTooltipIcon"]) {
  -webkit-text-stroke: 1px #ff3d9a;
  paint-order: stroke fill;
}

[data-testid="stVegaLiteChart"], [data-testid="stVegaLiteChart"] *,
[data-testid="stArrowVegaLiteChart"], [data-testid="stArrowVegaLiteChart"] *,
svg text {
  -webkit-text-stroke: 0 !important;
  text-shadow: none !important;
  font-family: "Source Sans Pro", -apple-system, sans-serif !important;
}

[data-testid="stVegaLiteChart"], [data-testid="stArrowVegaLiteChart"] {
  background: #0b1020 !important;
  border: 1px solid rgba(143,224,212,0.12) !important;
  border-radius: 10px !important;
  padding: 8px 16px 4px !important;
}

[data-testid="stVegaLiteChart"] svg,
[data-testid="stArrowVegaLiteChart"] svg { overflow: visible !important; }


section[data-testid="stSidebar"] {
  position: fixed !important;
  bottom: 0 !important; left: 0 !important; right: 0 !important; top: auto !important;
  width: 100vw !important;
  min-width: 100vw !important;
  max-width: 100vw !important;
  height: 232px !important;
  min-height: 232px !important;
  max-height: 232px !important;
  background: #070b16 !important;
  border-right: none !important;
  border-top: 1px solid rgba(6,182,212,0.18) !important;
  z-index: 100 !important;
  transform: none !important;
}

[data-testid="stSidebarResizer"],
[data-testid="stSidebarCollapseButton"] { display: none !important; }

section[data-testid="stSidebar"] > div:first-child {
  overflow-x: auto !important;
  overflow-y: auto !important;
  height: 232px !important;
  max-height: 232px !important;
  padding: 12px 24px 16px !important;
}

section[data-testid="stSidebar"] [data-testid="stVerticalBlock"] {
  display: flex !important;
  flex-direction: row !important;
  flex-wrap: wrap !important;
  align-items: flex-start !important;
  gap: 20px !important;
}
section[data-testid="stSidebar"] [data-testid="stVerticalBlock"] > div {
  flex: 0 0 auto !important;
  min-width: 170px !important;
  max-width: 260px !important;
}

section[data-testid="stSidebar"] [data-testid="stVerticalBlock"] [data-testid="stVerticalBlock"] {
  display: flex !important;
  flex-direction: column !important;
  flex-wrap: nowrap !important;
  gap: 0.6rem !important;
  min-width: unset !important;
  max-width: unset !important;
  flex: unset !important;
}
section[data-testid="stSidebar"] [data-testid="stVerticalBlock"] [data-testid="stVerticalBlock"] > div {
  min-width: unset !important;
  max-width: unset !important;
  flex: unset !important;
}

.main, [data-testid="stAppViewContainer"] > .main {
  margin-left: 0 !important;
}

[data-testid="stMain"] {
  height: calc(100vh - 236px) !important;
  max-height: calc(100vh - 236px) !important;
}
[data-testid="stMainBlockContainer"],
.block-container {
  padding-bottom: 48px !important;
}

section[data-testid="stSidebar"] p,
section[data-testid="stSidebar"] label,
section[data-testid="stSidebar"] .stMarkdown p,
section[data-testid="stSidebar"] .stWidgetLabel span { color: #94a3b8 !important; }


h1, h2, h3, h4 {
  color: #e2e8f0 !important;
  letter-spacing: -0.02em !important;
}


.gz-brand-wrap {
  padding: 1rem 0 0.75rem;
  border-bottom: 1px solid rgba(6,182,212,0.12);
  margin-bottom: 1.4rem;
}
.gz-brand-name {
  font-size: 2.9rem;
  font-weight: 700;
  letter-spacing: 0.04em;
  color: #8fe0d4;
  line-height: 1;
  display: inline-block;
  position: relative;
}

.gz-brand-mini {
  font-size: 1.25rem;
  font-weight: 700;
  letter-spacing: 0.06em;
  color: #8fe0d4 !important;
  line-height: 1;
}
.gz-brand-pill {
  display: inline-block;
  background: rgba(159,193,242,0.10);
  border: 1px solid rgba(159,193,242,0.30);
  color: #9fc1f2 !important;
  font-size: 0.65rem;
  font-weight: 600;
  letter-spacing: 0.06em;
  padding: 2px 8px;
  border-radius: 20px;
  margin-left: 10px;
  vertical-align: middle;
  position: relative;
  top: -8px;
}
.gz-brand-sub {
  font-size: 0.72rem;
  color: #475569;
  letter-spacing: 0.1em;
  text-transform: uppercase;
  font-weight: 400;
  margin-top: 5px;
}


.gz-section {
  font-size: 0.65rem !important;
  font-weight: 600 !important;
  text-transform: uppercase !important;
  letter-spacing: 0.14em !important;
  color: #8fe0d4 !important;
  margin: 1.8rem 0 0.85rem !important;
  padding-bottom: 0.45rem !important;
  border-bottom: 1px solid rgba(143,224,212,0.18) !important;
}


.gz-metric {
  background: rgba(255,255,255,0.025);
  border: 1px solid rgba(255,255,255,0.07);
  border-radius: 10px;
  padding: 16px 18px;
  position: relative;
  overflow: hidden;
  transition: border-color 0.2s ease, background 0.2s ease;
}
.gz-metric::after {
  content: '';
  position: absolute;
  top: 0; left: 0; right: 0;
  height: 2px;
  background: #06b6d4;
  opacity: 0.4;
}
.gz-metric:hover {
  border-color: rgba(6,182,212,0.3);
  background: rgba(6,182,212,0.04);
}
.gz-metric-v {
  font-size: 1.55rem;
  font-weight: 600;
  color: #e2e8f0;
  line-height: 1.1;
  letter-spacing: -0.02em;
}
.gz-metric-l {
  font-size: 0.65rem;
  color: #475569;
  margin-top: 6px;
  text-transform: uppercase;
  letter-spacing: 0.09em;
  font-weight: 500;
}


.gz-pass {
  background: rgba(16,185,129,0.07);
  border: 1px solid rgba(16,185,129,0.2);
  border-left: 3px solid #10b981;
  padding: 10px 14px;
  border-radius: 6px;
  color: #6ee7b7;
  font-size: 0.87rem;
  margin: 6px 0;
}
.gz-fail {
  background: rgba(239,68,68,0.07);
  border: 1px solid rgba(239,68,68,0.2);
  border-left: 3px solid #ef4444;
  padding: 10px 14px;
  border-radius: 6px;
  color: #fca5a5;
  font-size: 0.87rem;
  margin: 6px 0;
}
.gz-warn {
  background: rgba(245,158,11,0.07);
  border: 1px solid rgba(245,158,11,0.2);
  border-left: 3px solid #f59e0b;
  padding: 10px 14px;
  border-radius: 6px;
  color: #fcd34d;
  font-size: 0.87rem;
  margin: 6px 0;
}


.stButton > button {
  border-radius: 8px !important;
  font-weight: 500 !important;
  letter-spacing: 0.01em !important;
  transition: all 0.18s ease !important;
}
.stButton > button[kind="primary"] {
  background: #7fcfc4 !important;
  border: none !important;
  color: #0b1220 !important;
  font-weight: 600 !important;
}
.stButton > button[kind="primary"]:hover {
  background: #9be0d6 !important;
  transform: translateY(-1px) !important;
}


.stTextInput input,
.stNumberInput input {
  background: rgba(255,255,255,0.04) !important;
  border: 1px solid rgba(255,255,255,0.1) !important;
  border-radius: 8px !important;
  color: #e2e8f0 !important;
  transition: border-color 0.15s !important;
}
.stTextInput input:focus,
.stNumberInput input:focus {
  border-color: rgba(6,182,212,0.5) !important;
  box-shadow: 0 0 0 3px rgba(6,182,212,0.08) !important;
}
.stTextInput input::placeholder { color: #334155 !important; }


[data-testid="stFileUploaderDropzone"] {
  background: rgba(6,182,212,0.02) !important;
  border: 1.5px dashed rgba(6,182,212,0.22) !important;
  border-radius: 10px !important;
  transition: all 0.2s ease !important;
}
[data-testid="stFileUploaderDropzone"]:hover {
  border-color: rgba(6,182,212,0.45) !important;
  background: rgba(6,182,212,0.05) !important;
}


[data-testid="stDataFrame"] > div {
  border: 1px solid rgba(255,255,255,0.07) !important;
  border-radius: 8px !important;
  overflow: hidden !important;
}


[data-testid="stExpander"] {
  background: rgba(255,255,255,0.02) !important;
  border: 1px solid rgba(255,255,255,0.07) !important;
  border-radius: 8px !important;
}


hr {
  border: none !important;
  border-top: 1px solid rgba(255,255,255,0.06) !important;
}


::-webkit-scrollbar { width: 5px; height: 5px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb { background: rgba(6,182,212,0.2); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: rgba(6,182,212,0.4); }


[data-testid="stAppViewContainer"],
[data-testid="stMain"],
.block-container { background: transparent !important; }

.gz-ambient {
  position: fixed;
  inset: 0;
  z-index: -10;
  pointer-events: none;
  overflow: visible;
  will-change: transform;
}
.gz-dot-bg {
  position: absolute;
  border-radius: 50%;
  box-sizing: border-box;
  filter: blur(3px);
}


[data-baseweb="popover"] [role="listbox"],
ul[role="listbox"],
[data-baseweb="menu"] {
  background: #0d1426 !important;
  border: 1px solid rgba(6,182,212,0.35) !important;
  border-radius: 10px !important;
  box-shadow: 0 14px 44px rgba(0,0,0,0.55), 0 0 0 1px rgba(6,182,212,0.08) !important;
  overflow: hidden !important;
}
li[role="option"] { transition: background 0.12s ease !important; }
li[role="option"]:hover {
  background: rgba(6,182,212,0.14) !important;
}
li[role="option"][aria-selected="true"] {
  background: rgba(6,182,212,0.20) !important;
  box-shadow: inset 2px 0 0 #06b6d4 !important;
}


.gz-dots { display: inline-flex; gap: 9px; margin-left: 16px; vertical-align: middle; position: relative; top: -10px; }
.gz-dot {
  width: 11px; height: 11px; border-radius: 50%; display: inline-block;
  box-sizing: border-box;
  border: 2px solid #ffffff;
}
.gz-dot.teal  { background: #2dd4bf; }
.gz-dot.blue  { background: #4f9bff; }
.gz-dot.green { background: #34d399; }
.gz-dot.pink  { background: #ff5cae; }


.gz-section { position: relative; padding-left: 14px !important; }
.gz-section::before {
  content: ''; position: absolute; left: 0; top: 50%;
  width: 5px; height: 5px; margin-top: -4px; border-radius: 1px;
  background: #8fe0d4;
}


#MainMenu, .stDeployButton, footer,
[data-testid="stToolbar"], [data-testid="stDecoration"] { display: none !important; }
header[data-testid="stHeader"] { background: transparent !important; height: 0 !important; }
</style>
""", unsafe_allow_html=True)


@st.cache_resource
def _core():
    from core.mzml_parser     import (
        parse_mzml_apex,
        parse_mzml_xic_peak,
        coadd_centroid_spectra,
    )
    from core.peak_detector   import detect_peaks
    from core.peak_matcher    import match_peaks
    from core.impurity_engine import enumerate_impurities
    from core.mass_calculator import (
        calculate_monoisotopic_mass, calculate_average_mass,
        MEDBZ_LINKER_MONO, MENBZ_LINKER_MONO,
    )
    from core.batch_intelligence import get_root_cause
    return locals()


TAG_PRESETS = {
    "None / standard peptide":    {"suffix": "",        "fixed_da": 0.0,    "amide_default": False},
    "ArgTag (R6, no linker)":     {"suffix": "RRRRRR",  "fixed_da": 0.0,    "amide_default": True},
    "SynTag (MeDbz + R6)":        {"suffix": "RRRRRR",  "fixed_da": None,   "amide_default": True,
                                   "linker": "MEDBZ_LINKER_MONO"},
    "SynTag-MeNbz (cyclized)":    {"suffix": "RRRRRR",  "fixed_da": None,   "amide_default": True,
                                   "linker": "MENBZ_LINKER_MONO"},
    "Custom":                     {"suffix": "",        "fixed_da": 0.0,    "amide_default": False},
}


def _run_analysis(mzml_bytes, fname, sequence, c_amide,
                  fixed_mod_da, tolerance_ppm, coupling_reagent,
                  scan_window=3,
                  auto_window=True,
                  min_envelope_score=0.30, min_confidence=3.0,
                  triage_mode=False):
    c = _core()
    with tempfile.NamedTemporaryFile(suffix=".mzML", delete=False) as tf:
        tf.write(mzml_bytes)
        tmp = tf.name


    parent_mono = c["calculate_monoisotopic_mass"](
        sequence, c_terminal_amide=c_amide, fixed_modification_da=fixed_mod_da
    )
    parent_avg = c["calculate_average_mass"](
        sequence, c_terminal_amide=c_amide, fixed_modification_da=fixed_mod_da
    )

    window_spectra = []
    xic_info       = None
    try:
        if scan_window and scan_window > 0:


            (apex, window_spectra, tic_vals, tic_rts, n_ms1,
             parse_warnings, xic_info) = c["parse_mzml_xic_peak"](
                tmp, parent_mass=parent_mono, window_scans=scan_window,
                auto_window=auto_window,
                triage_mode=triage_mode,
            )
            if apex is None:
                return None, ["No usable MS1 spectrum found in the file."]
            spec = c["coadd_centroid_spectra"](window_spectra) if len(window_spectra) > 1 else apex
            extra_warns = []
            if len(window_spectra) > 1:
                _eff_win = xic_info.get('effective_window', scan_window)
                _fwhm_note = (f' [FWHM≈{xic_info["fwhm_scans"]} scans, auto]'
                              if xic_info.get('fwhm_scans') else '')
                extra_warns.append(
                    f'Multi-scan integration: co-added {len(window_spectra)} scans '
                    f'(apex ±{_eff_win} scans{_fwhm_note}, RT ~{apex.rt_min:.3f} min).'
                )
            if xic_info.get('fallback_to_tic'):
                extra_warns.append(
                    'XIC apex picker found no parent peak — fell back to TIC apex. '
                    'The reported apex may not contain the target peptide.'
                )
            elif xic_info.get('tolerance_used') is not None:
                extra_warns.append(
                    f'XIC apex picker located parent at RT {apex.rt_min:.3f} min '
                    f'using ±{xic_info["tolerance_used"]:.0f} ppm window.'
                )
            parse_warnings = list(parse_warnings) + extra_warns
        else:
            spec, tic_vals, tic_rts, n_ms1, parse_warnings = c["parse_mzml_apex"](tmp)
            apex = spec
            if spec is None:
                return None, ["No usable MS1 spectrum found in the file."]
    finally:
        os.unlink(tmp)

    pl = c["detect_peaks"](spec, threshold_fraction=0.003, triage_mode=triage_mode)


    cal_anchors = xic_info.get('apex_hits') if xic_info else None
    mr = c["match_peaks"](
        pl, sequence,
        tolerance_ppm=tolerance_ppm,
        c_terminal_amide=c_amide,
        fixed_modification_da=fixed_mod_da,
        coupling_reagent=coupling_reagent,
        min_envelope_score=min_envelope_score,
        min_confidence=min_confidence,
        source_scans=window_spectra if len(window_spectra) > 1 else None,
        cal_anchors=cal_anchors,
        triage_mode=triage_mode,
    )


    enriched = []
    for m in mr.matches:
        name, why, fix = c["get_root_cause"](m.impurity_type)
        enriched.append({
            "tier":            m.triage_tier,
            "impurity":        m.impurity_name,
            "type":            m.impurity_type,
            "obs_mz":          m.observed_mz,
            "ppm":             m.ppm_error,
            "z":               m.charge_state,
            "intensity":       m.observed_intensity,
            "rel_int_pct":     round(100 * m.observed_intensity / max(pl.base_peak_int, 1), 2),
            "iso_confirmed":   m.isotope_confirmed,
            "envelope_score":  m.envelope_score,
            "delta_da":        m.delta_mono,
            "predicted_mass":  m.predicted_mass,
            "risk":            m.risk_level,
            "confidence":      m.confidence,
            "position":        m.position,
            "n_degenerate":    m.n_degenerate,
            "scan_persistence": m.scan_persistence,
            "iso_warning":     m.iso_warning,
            "why":             why,
            "fix":             fix,
        })

    return {
        "sequence":     sequence,
        "c_amide":      c_amide,
        "fixed_mod_da": fixed_mod_da,
        "fname":        fname,
        "tolerance_ppm": tolerance_ppm,
        "coupling":     coupling_reagent,
        "triage_mode":  triage_mode,
        "parent_mono":  parent_mono,
        "parent_avg":   parent_avg,
        "scan_number":  apex.scan_number,
        "rt_min":       apex.rt_min,
        "tic_at_apex":  apex.total_ion_current,
        "n_analytical_scans":   n_ms1,
        "scan_window_used":     scan_window,
        "n_scans_coadded":      len(window_spectra) if window_spectra else 1,

        "xic_apex_used":             bool(xic_info and not xic_info.get('fallback_to_tic')),
        "xic_tolerance_ppm":         (xic_info.get('tolerance_used') if xic_info else None),
        "xic_fallback_to_tic":       bool(xic_info and xic_info.get('fallback_to_tic')),
        "xic_info_effective_window": (xic_info.get('effective_window') if xic_info else None),
        "n_peaks":      len(pl.peaks),
        "base_peak_mz": pl.base_peak_mz,
        "base_peak_int": pl.base_peak_int,
        "tic_rts":      tic_rts,
        "tic_vals":     tic_vals,
        "matches":      enriched,
        "n_predicted":  mr.n_impurities_predicted,
        "n_matched":    mr.n_peaks_matched,
        "bear_case":    mr.bear_case_applied,
        "warnings":     list(parse_warnings),
        "spectrum_mz":  list(spec.mz_array),
        "spectrum_int": list(spec.intensity_array),

        "parent_envelope_pct":     round(mr.parent_envelope_fraction * 100, 1),
        "matched_impurity_pct":    round(mr.matched_impurity_fraction * 100, 2),
        "unexplained_pct":         round(mr.unexplained_fraction * 100, 1),
        "n_unexplained_above_1pct": mr.n_unexplained_peaks_above_1pct,
        "top_unexplained":         mr.top_unexplained_peaks,

        "iemm_applied":            mr.iemm_applied,
        "native_parent_pct":       round(mr.native_parent_fraction * 100, 1),
        "deamidated_parent_pct":   round(mr.deamidated_parent_fraction * 100, 1),
    }, []


def _make_pdf(result):
    from app.report import build_pdf
    pdf_bytes = build_pdf(result)
    return pdf_bytes


def _tic_chart(rts, tics, apex_rt):
    df = pd.DataFrame({"RT (min)": rts, "TIC": tics})
    return df


def _spectrum_chart(mz_arr, int_arr, matches, parent_mono, parent_avg):
    return pd.DataFrame({"m/z": mz_arr, "intensity": int_arr})


def _sidebar():
    st.sidebar.markdown("""
<div style="padding:10px 4px 8px;">
  <div class="gz-brand-mini">GAZE</div>
  <div style="font-size:0.6rem;color:#334155;letter-spacing:0.1em;text-transform:uppercase;
              font-weight:500;margin-top:3px;">SPPS Impurity Attribution · v5.4</div>
</div>
<hr style="border:none;border-top:1px solid rgba(6,182,212,0.12);margin:0 0 6px;">
""", unsafe_allow_html=True)

    tol_ppm = st.sidebar.select_slider(
        "Mass tolerance",
        options=[2, 5, 10, 15, 20, 25, 30, 50],
        value=15,
        format_func=lambda v: f"{v} ppm",
    )
    st.sidebar.caption("5–10 ppm for Orbitrap, 15–30 ppm for Q-TOF / Waters.")

    analyst = st.sidebar.text_input("Analyst", value="", placeholder="optional")

    with st.sidebar.expander("Coupling reagent", expanded=False):
        coupling = st.radio(
            "Coupling reagent",
            ["HATU", "HBTU", "DIC", "PyBOP", "HATU/HOAt"],
            index=0,
            label_visibility="collapsed",
            help="Adds reagent-specific predicted impurities.",
        )

    with st.sidebar.expander("Advanced", expanded=False):
        triage_mode = st.checkbox(
            "Triage mode (return all candidates)",
            value=False,
            help="When ON, Gaze returns every mass-coincidence candidate above "
                 "the spectrum noise floor with its evidence (env_score, "
                 "confidence, persistence, isotope_confirmed) and a triage "
                 "tier (A/B/C). Use this when surveying for impurities — "
                 "you'll get 3-5× more rows but with everything visible so "
                 "you can dismiss noise at sight. Tier A is high-confidence; "
                 "Tier C is mostly noise. "
                 "When OFF (default), only high-confidence matches survive "
                 "the quality gates."
        )
        window_mode = st.selectbox(
            "Integration window",
            ["Auto (FWHM-based)", "Fixed 7 scans", "Fixed 15 scans", "Custom"],
            index=0,
            help="Auto: co-add window = FWHM ± 50% of the parent XIC peak. "
                 "Fixed: always co-add the specified number of scans on each side of the apex.",
        )
        if window_mode == "Auto (FWHM-based)":
            scan_window = 3
            auto_window = True
        elif window_mode == "Fixed 7 scans":
            scan_window = 3
            auto_window = False
        elif window_mode == "Fixed 15 scans":
            scan_window = 7
            auto_window = False
        else:
            scan_window = st.slider("Custom window (scans each side)", 0, 10, 3, 1)
            auto_window = False
        if triage_mode:
            st.caption(
                "ℹ️ Triage mode is ON — the envelope and confidence sliders "
                "below are ignored (all gates dropped, every candidate kept)."
            )
        min_env = st.slider(
            "Min envelope-fit score",
            min_value=0.0, max_value=1.0, value=0.30, step=0.05,
            disabled=triage_mode,
            help="Reject candidate matches whose isotope envelope shape doesn't "
                 "fit the theoretical averagine prediction at least this well. "
                 "0 = disabled (v2 behavior). 0.3 is a permissive default; raise "
                 "to 0.5+ for stricter filtering.",
        )
        min_conf = st.slider(
            "Min confidence score",
            min_value=0.0, max_value=10.0, value=3.0, step=0.5,
            disabled=triage_mode,
            help="Reject matches with confidence below this threshold. Confidence "
                 "combines risk weight, ppm error, intensity, and envelope quality.",
        )

    return {
        "tolerance_ppm":      int(tol_ppm),
        "coupling":           coupling,
        "analyst":            analyst,
        "scan_window":        int(scan_window),
        "auto_window":        bool(auto_window),
        "min_envelope_score": float(min_env),
        "min_confidence":     float(min_conf),
        "triage_mode":        bool(triage_mode),
    }


def _input_form():
    st.markdown('<p class="gz-section">1 · Sequence &amp; Construct</p>', unsafe_allow_html=True)

    cA, cB = st.columns([3, 2])

    with cA:
        sequence = st.text_input(
            "Peptide sequence (single-letter, no spaces)",
            value="",
            placeholder="e.g. HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
            help="Use single-letter amino acid codes. The C-terminal tag is added "
                 "automatically based on your tag selection below.",
        ).strip().upper()

    with cB:
        tag_label = st.selectbox(
            "C-terminal tag",
            list(TAG_PRESETS.keys()),
            index=0,
            help="ArgTag = +RRRRRR (hexaarginine tag). "
                 "SynTag = +MeDbz linker + RRRRRR. "
                 "Cyclized form has -H2O from the MeDbz → MeNbz cyclization.",
        )

    preset = TAG_PRESETS[tag_label]
    suffix = preset["suffix"]
    if "linker" in preset:
        c = _core()
        fixed_mod_da = c[preset["linker"]]
    elif preset["fixed_da"] is None:
        fixed_mod_da = 0.0
    else:
        fixed_mod_da = preset["fixed_da"]
    amide_default = preset["amide_default"]

    cC, cD = st.columns([1, 1])
    with cC:
        c_amide = st.checkbox(
            "C-terminal amide (Rink amide resin)",
            value=amide_default,
            help="Check if the resin was Rink amide. Uncheck for free acid (Wang resin etc.).",
        )
    with cD:
        if tag_label == "Custom":
            fixed_mod_da = st.number_input(
                "Custom C-terminal modification (Da)",
                value=0.0, step=0.001, format="%.4f",
                help="Mass to add to predicted parent mono mass. Use for any covalent "
                     "C-terminal modification not covered by the preset list.",
            )
        else:
            st.write("")

    full_sequence = sequence + suffix
    sequence_summary = full_sequence
    if not full_sequence:
        sequence_summary = "(enter a sequence)"

    st.markdown('<p class="gz-section">2 · mzML File</p>', unsafe_allow_html=True)
    upload = st.file_uploader(
        "MS1 spectrum file (.mzML)",
        type=["mzML"],
        accept_multiple_files=False,
        help="Convert vendor data with MSConvert (centroided, peak-picked).",
    )


    if full_sequence:
        c = _core()
        try:
            mono = c["calculate_monoisotopic_mass"](
                full_sequence, c_terminal_amide=c_amide,
                fixed_modification_da=fixed_mod_da,
            )
            avg = c["calculate_average_mass"](
                full_sequence, c_terminal_amide=c_amide,
                fixed_modification_da=fixed_mod_da,
            )
            st.caption(
                f"Sequence used for analysis: **{full_sequence}** ({len(full_sequence)} residues)  · "
                f"Predicted parent **{mono:.3f} Da** monoisotopic, {avg:.2f} Da average  · "
                f"C-amide: {'yes' if c_amide else 'no'}  · "
                f"Linker mass: +{fixed_mod_da:.4f} Da" if fixed_mod_da else
                f"Sequence used for analysis: **{full_sequence}** ({len(full_sequence)} residues)  · "
                f"Predicted parent **{mono:.3f} Da** monoisotopic, {avg:.2f} Da average  · "
                f"C-amide: {'yes' if c_amide else 'no'}"
            )
        except Exception as e:
            st.warning(f"Cannot compute parent mass: {e}")

    return {
        "sequence":     full_sequence,
        "c_amide":      c_amide,
        "fixed_mod_da": fixed_mod_da,
        "tag_label":    tag_label,
        "upload":       upload,
    }


def _results_page(result):
    sb = st.session_state.get("sidebar_state", {})
    analyst = sb.get("analyst", "")


    parent_pct = result.get("parent_envelope_pct", 0)
    if parent_pct < 10.0:
        st.warning(
            f"**Parent envelope only {parent_pct:.1f}% of base-peak intensity.** "
            f"The apex picker may have locked onto a noise apex that doesn't "
            f"contain the target peptide. The impurity table below should be "
            f"treated as a list of mass-coincident peaks at this retention "
            f"time, not as confident impurity assignments. Cross-check with "
            f"UV chromatogram or known retention time before drawing conclusions."
        )
    elif result.get("xic_fallback_to_tic"):
        st.warning(
            "XIC apex picker found no parent envelope at the predicted m/z; "
            "fell back to TIC apex. The reported apex may not contain the "
            "target peptide."
        )


    st.markdown('<p class="gz-section">Analysis Summary</p>', unsafe_allow_html=True)
    cols = st.columns(4)
    cols[0].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result["parent_mono"]:.3f}</div>'
        f'<div class="gz-metric-l">parent monoisotopic (Da)</div></div>',
        unsafe_allow_html=True)
    cols[1].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result["n_predicted"]}</div>'
        f'<div class="gz-metric-l">impurities predicted</div></div>',
        unsafe_allow_html=True)
    cols[2].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result["n_matched"]}</div>'
        f'<div class="gz-metric-l">matched in spectrum</div></div>',
        unsafe_allow_html=True)
    cols[3].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result["scan_number"]}</div>'
        f'<div class="gz-metric-l">apex scan (RT {result["rt_min"]:.2f} min)</div></div>',
        unsafe_allow_html=True)


    st.markdown("&nbsp;", unsafe_allow_html=True)
    cov_cols = st.columns(4)
    cov_cols[0].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result.get("parent_envelope_pct", 0):.1f}%</div>'
        f'<div class="gz-metric-l">parent envelope</div></div>',
        unsafe_allow_html=True)
    cov_cols[1].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result.get("matched_impurity_pct", 0):.2f}%</div>'
        f'<div class="gz-metric-l">matched impurities</div></div>',
        unsafe_allow_html=True)
    cov_cols[2].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result.get("unexplained_pct", 0):.1f}%</div>'
        f'<div class="gz-metric-l">unexplained intensity</div></div>',
        unsafe_allow_html=True)
    n_coadded = result.get("n_scans_coadded", 1)
    _xic_eff_win = result.get("xic_info_effective_window")
    if n_coadded > 1 and _xic_eff_win is not None:
        coadd_label = f'{n_coadded}-scan coadd (±{_xic_eff_win})'
    elif n_coadded > 1:
        coadd_label = f'{n_coadded}-scan coadd'
    else:
        coadd_label = 'single apex scan'
    cov_cols[3].markdown(
        f'<div class="gz-metric"><div class="gz-metric-v">{result.get("n_unexplained_above_1pct", 0)}</div>'
        f'<div class="gz-metric-l">unexplained peaks &gt;1% of base · {coadd_label}</div></div>',
        unsafe_allow_html=True)


    for w in result["warnings"]:
        st.markdown(f'<div class="gz-warn">⚠ {w}</div>', unsafe_allow_html=True)


    high_intensity_matches = [m for m in result["matches"]
                              if m["rel_int_pct"] > 5.0 and m["risk"] == "HIGH"]
    if high_intensity_matches:
        st.markdown(
            '<div class="gz-fail">Major impurities detected at >5% of base peak. '
            'Synthesis quality is likely compromised. See table below.</div>',
            unsafe_allow_html=True)
    else:
        st.markdown(
            '<div class="gz-pass">No major impurities at &gt;5% of base peak. '
            'See table for trace-level findings.</div>',
            unsafe_allow_html=True)


    unexplained_pct = result.get("unexplained_pct", 0)
    if unexplained_pct > 20:
        st.markdown(
            f'<div class="gz-warn">{unexplained_pct:.0f}% of detected intensity is not '
            f'attributed to the parent envelope or any predicted impurity. '
            f'See "Top unexplained peaks" below for the largest unaccounted-for signals — '
            f'they may indicate impurities not in Gaze\'s catalogue, dimers/multimers, '
            f'or co-eluting species.</div>',
            unsafe_allow_html=True)

    st.divider()


    st.markdown('<p class="gz-section">Matched Impurities</p>', unsafe_allow_html=True)
    if not result["matches"]:
        st.info("No predicted impurities matched any peak above the intensity threshold. "
                "Either the synthesis is very clean, or your sequence/linker setting is wrong "
                "(check the predicted parent mass against the dominant peak in the spectrum).")
    else:
        df = pd.DataFrame(result["matches"])


        triage_on = result.get("triage_mode", False)
        display_cols = {}
        if triage_on:
            display_cols["tier"] = "Tier"
        display_cols.update({
            "impurity":         "Impurity",
            "position":         "Pos",
            "z":                "z",
            "obs_mz":           "Obs m/z",
            "ppm":               "ppm",
            "rel_int_pct":      "% base",
            "envelope_score":   "Env",
            "iso_confirmed":    "Iso?",
            "scan_persistence": "Pers",
            "iso_warning":      "⚠ Iso?",
            "confidence":       "Conf",
            "delta_da":         "Δ Da",
            "risk":             "Risk",
            "why":              "Likely cause",
            "fix":              "Recommended action",
        })

        df_cols = [c for c in display_cols if c in df.columns]
        df_show = df[df_cols].rename(columns={c: display_cols[c] for c in df_cols})


        if triage_on and "Tier" in df_show.columns:
            tier_rank = {"A": 0, "B": 1, "C": 2, "": 3}
            df_show["_r"] = df_show["Tier"].map(lambda t: tier_rank.get(t, 9))
            df_show = df_show.sort_values(["_r", "% base"],
                                          ascending=[True, False]).drop(columns="_r")
        elif "Conf" in df_show.columns:
            df_show = df_show.sort_values(["Conf", "% base"], ascending=[False, False])
        else:
            risk_rank = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
            df_show["_r"] = df_show["Risk"].map(lambda r: risk_rank.get(r, 9))
            df_show = df_show.sort_values(["_r", "% base"],
                                          ascending=[True, False]).drop(columns="_r")
        if triage_on:
            st.caption(
                "**Triage mode is on.** Tier A is high-confidence; "
                "Tier B is plausible; Tier C is mostly noise — scan and dismiss."
            )
        st.dataframe(df_show, use_container_width=True, hide_index=True)


    top_unexp = result.get("top_unexplained", [])
    if top_unexp:
        with st.expander(f"Top unexplained peaks ({len(top_unexp)})", expanded=False):
            st.caption(
                "Peaks not assigned to the parent envelope and not matched to any "
                "predicted impurity. Sorted by descending intensity. Inspect these "
                "manually to determine if they represent impurities outside Gaze's "
                "catalogue, parent multimers, or co-eluting species."
            )
            unexp_df = pd.DataFrame(top_unexp, columns=["m/z", "intensity", "% base"])
            st.dataframe(unexp_df, use_container_width=True, hide_index=True)

    st.divider()


    st.markdown('<p class="gz-section">Apex Scan Preview</p>', unsafe_allow_html=True)
    st.caption(f"Scan {result['scan_number']} at RT {result['rt_min']:.3f} min · "
               f"{result['n_peaks']} peaks above 0.3% of base · "
               f"base peak at m/z {result['base_peak_mz']:.4f}")


    pairs = sorted(zip(result["spectrum_mz"], result["spectrum_int"]),
                   key=lambda x: -x[1])[:200]
    df_peaks = pd.DataFrame(pairs, columns=["m/z", "intensity"])
    st.bar_chart(df_peaks.set_index("m/z"), height=240)

    st.divider()


    st.markdown('<p class="gz-section">Download Report</p>', unsafe_allow_html=True)
    cD1, cD2 = st.columns(2)
    with cD1:
        try:
            pdf = _make_pdf({**result, "analyst": analyst})
            ts  = datetime.now().strftime("%Y%m%d-%H%M%S")
            fname = (result["fname"] or "gaze").rsplit(".", 1)[0]
            st.download_button(
                "PDF report",
                data=pdf,
                file_name=f"gaze_{fname}_{ts}.pdf",
                mime="application/pdf",
                use_container_width=True,
            )
        except Exception as e:
            st.error(f"PDF generation failed: {e}")

    with cD2:
        if result["matches"]:
            csv = pd.DataFrame(result["matches"]).to_csv(index=False).encode("utf-8")
            st.download_button(
                "Matches as CSV",
                data=csv,
                file_name=f"gaze_matches.csv",
                mime="text/csv",
                use_container_width=True,
            )


def main():
    st.markdown("""
<div class="gz-ambient">
  <span class="gz-dot-bg" style="top:6vh;    left:82vw; width:130px; height:130px; background:rgba(45,212,191,0.11);"></span>
  <span class="gz-dot-bg" style="top:22vh;   left:11vw; width:90px;  height:90px;  background:rgba(96,165,250,0.10);"></span>
  <span class="gz-dot-bg" style="top:48vh;   left:91vw; width:110px; height:110px; background:rgba(255,92,174,0.10);"></span>
  <span class="gz-dot-bg" style="top:72vh;   left:30vw; width:85px;  height:85px;  background:rgba(52,211,153,0.10);"></span>
  <span class="gz-dot-bg" style="top:90vh;   left:64vw; width:100px; height:100px; background:rgba(45,212,191,0.10);"></span>
  <span class="gz-dot-bg" style="top:-18vh;  left:48vw; width:95px;  height:95px;  background:rgba(255,92,174,0.11);"></span>
  <span class="gz-dot-bg" style="top:-48vh;  left:8vw;  width:120px; height:120px; background:rgba(96,165,250,0.11);"></span>
  <span class="gz-dot-bg" style="top:-78vh;  left:74vw; width:100px; height:100px; background:rgba(52,211,153,0.11);"></span>
  <span class="gz-dot-bg" style="top:-112vh; left:35vw; width:110px; height:110px; background:rgba(45,212,191,0.11);"></span>
  <span class="gz-dot-bg" style="top:-148vh; left:88vw; width:90px;  height:90px;  background:rgba(255,92,174,0.10);"></span>
  <span class="gz-dot-bg" style="top:-182vh; left:18vw; width:105px; height:105px; background:rgba(96,165,250,0.11);"></span>
  <span class="gz-dot-bg" style="top:-216vh; left:58vw; width:95px;  height:95px;  background:rgba(52,211,153,0.11);"></span>
  <span class="gz-dot-bg" style="top:-248vh; left:5vw;  width:110px; height:110px; background:rgba(45,212,191,0.11);"></span>
</div>
<div class="gz-brand-wrap">
  <div>
    <span class="gz-brand-name">GAZE</span>
    <span class="gz-brand-pill">v5.4</span>
    <span class="gz-dots">
      <span class="gz-dot teal"></span>
      <span class="gz-dot blue"></span>
      <span class="gz-dot green"></span>
      <span class="gz-dot pink"></span>
    </span>
  </div>
  <div class="gz-brand-sub">Post-synthesis SPPS impurity attribution · LC-MS · MS1 · Research use only — confirm with RT / MS/MS</div>
</div>
""", unsafe_allow_html=True)

    components.html("""
<script>
(function(){
  try {
    var doc = window.parent.document, FACTOR = 1.40;
    function apply(top){
      var amb = doc.querySelector('.gz-ambient');
      if (amb) amb.style.transform = 'translate3d(0,' + (top * FACTOR) + 'px,0)';
    }
    function attach(el, getTop){
      if (!el || el.__gzPar) return;
      el.__gzPar = true;
      el.addEventListener('scroll', function(){ apply(getTop()); }, {passive:true});
    }
    function init(){
      var main = doc.querySelector('[data-testid="stMain"]');
      if (main) attach(main, function(){ return main.scrollTop; });
      attach(window.parent, function(){
        return window.parent.scrollY || (doc.scrollingElement ? doc.scrollingElement.scrollTop : 0);
      });
    }
    init(); setTimeout(init, 400); setTimeout(init, 1200);
  } catch (e) {}
})();
</script>
""", height=0)

    sidebar_state = _sidebar()
    st.session_state["sidebar_state"] = sidebar_state

    inputs = _input_form()

    can_run = bool(inputs["sequence"]) and inputs["upload"] is not None
    if not can_run:
        if not inputs["sequence"]:
            st.info("Enter a peptide sequence above, then upload an mzML file to run.")
        else:
            st.info("Upload an mzML file to run the analysis.")
        return

    if st.button("Run analysis", type="primary", use_container_width=True):
        with st.spinner("Parsing mzML, detecting peaks, matching impurities…"):
            data = inputs["upload"].getvalue()
            result, errs = _run_analysis(
                data, inputs["upload"].name,
                inputs["sequence"],
                inputs["c_amide"],
                inputs["fixed_mod_da"],
                sidebar_state["tolerance_ppm"],
                sidebar_state["coupling"],
                scan_window       = sidebar_state["scan_window"],
                auto_window       = sidebar_state.get("auto_window", True),
                min_envelope_score= sidebar_state["min_envelope_score"],
                min_confidence    = sidebar_state["min_confidence"],
                triage_mode       = sidebar_state["triage_mode"],
            )
        if errs:
            for e in errs:
                st.error(e)
            return
        st.session_state["last_result"] = result

    if "last_result" in st.session_state:
        st.divider()
        _results_page(st.session_state["last_result"])


if __name__ == "__main__":
    main()
