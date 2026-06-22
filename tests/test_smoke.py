"""Smoke tests for the Gaze pipeline against the SynTag mzML fixture."""
import os
import sys

try:
    import pytest
    _HAS_PYTEST = True
except ImportError:
    _HAS_PYTEST = False
    class _PytestStub:
        @staticmethod
        def approx(value, abs=None, rel=None):
            class _Approx:
                def __eq__(self, other):
                    if abs is not None: return __builtins__.abs(other - value) <= abs
                    if rel is not None: return __builtins__.abs(other - value) <= __builtins__.abs(value) * rel
                    return other == value
            return _Approx()
        @staticmethod
        def mark_skipif(cond, reason=""):
            def deco(fn):
                fn._skip_if_true = (cond, reason)
                return fn
            return deco
    pytest = _PytestStub()

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

SYNTAG_PATHS = [
    os.path.join(_ROOT, "tests", "fixtures", "Galaxy2-_Syntag_mzML_.mzml"),
    "/mnt/user-data/uploads/Galaxy2-_Syntag_mzML_.mzml",
]
SYNTAG_FILE = next((p for p in SYNTAG_PATHS if os.path.exists(p)), None)
SKIP_REASON = (
    "SynTag fixture not found. Place Galaxy2-_Syntag_mzML_.mzml in "
    "tests/fixtures/ or update SYNTAG_PATHS in test_smoke.py."
)
if _HAS_PYTEST:
    needs_fixture = pytest.mark.skipif(SYNTAG_FILE is None, reason=SKIP_REASON)
else:
    def needs_fixture(fn):
        fn._skip_if_no_fixture = SYNTAG_FILE is None
        return fn

SEQUENCE        = "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRRRRRRR"
EXPECTED_PARENT = 4232.269
PARENT_TOL_DA   = 0.5
EXPECTED_RT_MIN = 4.80
RT_TOL_MIN      = 0.10


@needs_fixture
def test_parser_filters_waters_multifunction():
    """parse_mzml_apex keeps only function=1 and warns on multi-function files."""
    from core.mzml_parser import parse_mzml_apex

    spec, tic_vals, tic_rts, n_ms1, warnings = parse_mzml_apex(SYNTAG_FILE)
    assert spec is not None, "no apex spectrum returned"
    assert n_ms1 == 538, f"expected 538 analytical MS1 scans, got {n_ms1}"
    assert spec.function_id == 1, f"apex should be from function=1, got {spec.function_id}"
    assert any("multi-function" in w.lower() for w in warnings), \
        "expected a Waters multi-function warning"


@needs_fixture
def test_apex_scan_picked():
    """Apex should be at ~RT 4.80, scan 258 in this file."""
    from core.mzml_parser import parse_mzml_apex

    spec, *_ = parse_mzml_apex(SYNTAG_FILE)
    assert abs(spec.rt_min - EXPECTED_RT_MIN) < RT_TOL_MIN, \
        f"apex RT {spec.rt_min:.3f} not near {EXPECTED_RT_MIN}"
    assert 845 < spec.mz_array.max() > 1000, "apex spectrum has implausibly small m/z range"
    assert spec.total_ion_current > 1e6, "apex TIC is suspiciously small"


@needs_fixture
def test_parent_mass_matches():
    """GLP-1(7-36) + ArgTag with C-amide should give 4232.27 Da."""
    from core.mass_calculator import calculate_monoisotopic_mass
    M = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)
    assert abs(M - EXPECTED_PARENT) < 0.01, \
        f"predicted parent {M:.4f} differs from expected {EXPECTED_PARENT} by >0.01 Da"


@needs_fixture
def test_self_calibration_reduces_error():
    """Linear self-calibration should keep residuals under 20 ppm at every parent anchor."""
    import numpy as np
    from core.mzml_parser import parse_mzml_apex
    from core.peak_detector import detect_peaks
    from core.peak_matcher import _self_calibrate
    from core.mass_calculator import calculate_mz, calculate_monoisotopic_mass

    spec, *_ = parse_mzml_apex(SYNTAG_FILE)
    pl = detect_peaks(spec, threshold_fraction=0.003)
    M  = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)
    correct = _self_calibrate(pl, SEQUENCE, c_terminal_amide=True)

    for z, raw_obs in [(3, 1411.864), (4, 1059.126), (5, 847.484), (6, 706.391)]:
        theo = calculate_mz(M, z)
        corrected = correct(raw_obs)
        residual_ppm = (corrected - theo) / theo * 1e6
        assert abs(residual_ppm) < 20, \
            f"z={z}: post-calibration residual {residual_ppm:+.1f} ppm > 20 ppm"


@needs_fixture
def test_impurity_match_clean_spectrum():
    """On the clean SynTag file the matcher reports few matches, none above 10% of base peak."""
    from core.mzml_parser import parse_mzml_apex
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    spec, *_ = parse_mzml_apex(SYNTAG_FILE)
    pl = detect_peaks(spec, threshold_fraction=0.003)
    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0,
                     c_terminal_amide=True, coupling_reagent="HATU")

    assert abs(mr.parent_mass - EXPECTED_PARENT) < PARENT_TOL_DA, \
        f"parent {mr.parent_mass:.4f} differs from expected {EXPECTED_PARENT} by >{PARENT_TOL_DA}"
    assert mr.n_peaks_matched <= 12, \
        f"matched {mr.n_peaks_matched} impurities; this is supposed to be a clean spectrum"
    assert mr.bear_case_applied, \
        "bear case should be applied for a 4.2 kDa peptide"

    base = pl.base_peak_int
    for m in mr.matches:
        rel = 100 * m.observed_intensity / base
        assert rel < 10.0, (
            f"impurity {m.impurity_name} is at {rel:.1f}% of base peak. "
            f"Either the spectrum is not clean, or a parent-isotope match "
            f"slipped through the collision guard."
        )


@needs_fixture
def test_chromatographic_peak_window():
    """parse_mzml_chromatographic_peak returns a window of scans including the apex."""
    from core.mzml_parser import parse_mzml_chromatographic_peak

    apex, window, tic_v, tic_r, n_ms1, warns = parse_mzml_chromatographic_peak(
        SYNTAG_FILE, window_scans=3)
    assert apex is not None
    assert apex.scan_number == 258, f"expected apex scan 258, got {apex.scan_number}"
    assert len(window) == 7, f"expected 7 scans (apex ±3), got {len(window)}"
    apex_in_window = any(s.scan_number == apex.scan_number for s in window)
    assert apex_in_window, "apex spectrum must be inside the window"


@needs_fixture
def test_coadd_centroid_spectra_increases_signal():
    """Co-adding scans should grow the base-peak intensity and surface more peaks."""
    from core.mzml_parser import parse_mzml_chromatographic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks

    apex, window, *_ = parse_mzml_chromatographic_peak(SYNTAG_FILE, window_scans=3)
    coadd = coadd_centroid_spectra(window)
    pl_apex  = detect_peaks(apex,  threshold_fraction=0.003)
    pl_coadd = detect_peaks(coadd, threshold_fraction=0.003)
    assert pl_coadd.n_peaks >= pl_apex.n_peaks, (
        f"coadd has {pl_coadd.n_peaks} peaks vs {pl_apex.n_peaks} in single apex"
    )
    ratio = pl_coadd.base_peak_int / pl_apex.base_peak_int
    assert ratio > 2.5, f"coadd base peak {pl_coadd.base_peak_int:.0f} not "\
                         f"sufficiently larger than single-apex {pl_apex.base_peak_int:.0f}"


@needs_fixture
def test_spectrum_coverage_metrics():
    """MatchResult coverage fractions should sum to ~100%."""
    from core.mzml_parser import parse_mzml_apex
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    spec, *_ = parse_mzml_apex(SYNTAG_FILE)
    pl = detect_peaks(spec, threshold_fraction=0.003)
    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0,
                     c_terminal_amide=True, coupling_reagent="HATU")

    assert mr.total_intensity > 0, "total_intensity not set"
    assert 0 <= mr.parent_envelope_fraction <= 1.0
    assert 0 <= mr.matched_impurity_fraction <= 1.0
    assert 0 <= mr.unexplained_fraction <= 1.0

    total_frac = (mr.parent_envelope_fraction +
                  mr.matched_impurity_fraction +
                  mr.unexplained_fraction)
    assert 0.99 <= total_frac <= 1.01, \
        f"coverage fractions sum to {total_frac:.3f}, expected ~1.0"

    assert mr.parent_envelope_fraction > 0.5, \
        f"parent envelope only {mr.parent_envelope_fraction*100:.1f}% of intensity"


@needs_fixture
def test_envelope_score_on_matches():
    """Every match should carry an envelope_score in [0, 1]."""
    from core.mzml_parser import parse_mzml_apex
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    spec, *_ = parse_mzml_apex(SYNTAG_FILE)
    pl = detect_peaks(spec, threshold_fraction=0.003)
    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0,
                     c_terminal_amide=True, coupling_reagent="HATU")

    for m in mr.matches:
        assert hasattr(m, 'envelope_score')
        assert 0.0 <= m.envelope_score <= 1.0, \
            f"envelope_score {m.envelope_score} out of [0,1] for {m.impurity_name}"


@needs_fixture
def test_envelope_filter_reduces_false_positives():
    """Compared to envelope-disabled, envelope filtering should reduce match count."""
    from core.mzml_parser import parse_mzml_apex
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    spec, *_ = parse_mzml_apex(SYNTAG_FILE)
    pl = detect_peaks(spec, threshold_fraction=0.003)
    mr_off = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                         coupling_reagent="HATU",
                         min_envelope_score=0.0, min_confidence=0.0)
    mr_on  = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                         coupling_reagent="HATU",
                         min_envelope_score=0.30, min_confidence=3.0)
    assert mr_on.n_peaks_matched <= mr_off.n_peaks_matched, \
        f"envelope filter increased matches ({mr_off.n_peaks_matched} -> {mr_on.n_peaks_matched})"


@needs_fixture
def test_mass_degeneracy_collapsed():
    """The 7 single-Arg deletions should collapse to one row with n_degenerate=7."""
    from core.impurity_engine import enumerate_impurities, collapse_mass_degenerate

    df = enumerate_impurities(SEQUENCE, c_terminal_amide=True)
    df = collapse_mass_degenerate(df)
    assert 'n_degenerate' in df.columns

    r_dels = df[df['impurity_type'] == 'deletion_single'].copy()
    r_collapsed = r_dels[r_dels['n_degenerate'] >= 7]
    assert len(r_collapsed) >= 1, \
        f"R-deletion mass-degeneracy not collapsed; got {dict(r_dels['n_degenerate'].value_counts())}"
    pos = str(r_collapsed.iloc[0]['position'])
    assert pos.count('/') >= 6, f"position string '{pos}' should list 7 R positions"


@needs_fixture
def test_persistence_filter_drops_coadd_artifacts():
    """With source_scans, matches absent from every individual scan are dropped."""
    from core.mzml_parser import parse_mzml_chromatographic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    apex, window, *_ = parse_mzml_chromatographic_peak(SYNTAG_FILE, window_scans=3)
    coadd = coadd_centroid_spectra(window)
    pl = detect_peaks(coadd, threshold_fraction=0.003)

    mr_no_filter = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0,
                               c_terminal_amide=True, coupling_reagent="HATU")
    mr_with_filter = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0,
                                 c_terminal_amide=True, coupling_reagent="HATU",
                                 source_scans=window)
    assert mr_with_filter.n_peaks_matched < mr_no_filter.n_peaks_matched, \
        "persistence filter should drop coadd-only artifacts"
    for m in mr_with_filter.matches:
        assert m.scan_persistence >= 0.30, \
            f"match {m.impurity_name} survived with persistence={m.scan_persistence}"


@needs_fixture
def test_persistence_filter_disabled_when_no_window():
    """With source_scans=None the filter is a no-op and every match has persistence=1.0."""
    from core.mzml_parser import parse_mzml_apex
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    spec, *_ = parse_mzml_apex(SYNTAG_FILE)
    pl = detect_peaks(spec, threshold_fraction=0.003)
    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                     coupling_reagent="HATU")
    for m in mr.matches:
        assert m.scan_persistence == 1.0, \
            f"single-apex match has persistence != 1.0: {m.scan_persistence}"


@needs_fixture
def test_v4_clean_synthesis_reports_few_matches():
    """End-to-end: the clean SynTag synthesis should report between 1 and 5 matches."""
    from core.mzml_parser import parse_mzml_chromatographic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    apex, window, *_ = parse_mzml_chromatographic_peak(SYNTAG_FILE, window_scans=3)
    coadd = coadd_centroid_spectra(window)
    pl = detect_peaks(coadd, threshold_fraction=0.003)
    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                     coupling_reagent="HATU", source_scans=window)
    assert mr.n_peaks_matched <= 5, \
        f"clean SynTag synthesis reports {mr.n_peaks_matched} matches; expected ≤5"
    assert mr.n_peaks_matched >= 1, \
        f"expected at least 1 real match (ammonia_loss), got 0"


@needs_fixture
def test_triage_mode_returns_more_than_default():
    """Triage mode should surface more candidates than default mode (≥10) on a clean file."""
    from core.mzml_parser import parse_mzml_chromatographic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    apex, window, *_ = parse_mzml_chromatographic_peak(SYNTAG_FILE, window_scans=3)
    coadd = coadd_centroid_spectra(window)
    pl_default = detect_peaks(coadd, threshold_fraction=0.003)
    pl_triage  = detect_peaks(coadd, triage_mode=True)

    mr_default = match_peaks(pl_default, SEQUENCE, tolerance_ppm=15.0,
                             c_terminal_amide=True, coupling_reagent="HATU",
                             source_scans=window)
    mr_triage  = match_peaks(pl_triage, SEQUENCE, tolerance_ppm=15.0,
                             c_terminal_amide=True, coupling_reagent="HATU",
                             source_scans=window, triage_mode=True)

    assert mr_triage.n_peaks_matched > mr_default.n_peaks_matched, (
        f"triage_mode produced {mr_triage.n_peaks_matched} matches; "
        f"default produced {mr_default.n_peaks_matched}. "
        f"triage_mode should return strictly more on a clean synthesis."
    )
    assert mr_triage.n_peaks_matched >= 10, (
        f"triage_mode should surface ≥10 candidates on a clean synthesis; "
        f"got {mr_triage.n_peaks_matched}"
    )


@needs_fixture
def test_triage_mode_tags_every_match_with_tier():
    """Every match returned by triage_mode must have triage_tier in {A,B,C}."""
    from core.mzml_parser import parse_mzml_chromatographic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    apex, window, *_ = parse_mzml_chromatographic_peak(SYNTAG_FILE, window_scans=3)
    coadd = coadd_centroid_spectra(window)
    pl    = detect_peaks(coadd, triage_mode=True)
    mr    = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                        coupling_reagent="HATU", source_scans=window,
                        triage_mode=True)

    valid = {'A', 'B', 'C'}
    for m in mr.matches:
        assert m.triage_tier in valid, (
            f"match {m.impurity_name} z={m.charge_state} has invalid "
            f"triage_tier={m.triage_tier!r}"
        )

    seen_b = False
    seen_c = False
    for m in mr.matches:
        if m.triage_tier == 'B':
            seen_b = True
        if m.triage_tier == 'C':
            seen_c = True
        if seen_b and m.triage_tier == 'A':
            assert False, "found tier-A match after tier-B match (sort broken)"
        if seen_c and m.triage_tier in ('A', 'B'):
            assert False, f"found tier-{m.triage_tier} match after tier-C (sort broken)"


@needs_fixture
def test_triage_mode_default_off_preserves_old_behavior():
    """triage_mode is off by default; existing callers must see no behavior change."""
    from core.mzml_parser import parse_mzml_chromatographic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks

    apex, window, *_ = parse_mzml_chromatographic_peak(SYNTAG_FILE, window_scans=3)
    coadd = coadd_centroid_spectra(window)
    pl    = detect_peaks(coadd, threshold_fraction=0.003)
    mr    = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                        coupling_reagent="HATU", source_scans=window)

    for m in mr.matches:
        assert m.triage_tier == '', (
            f"default-mode match has unexpected triage_tier={m.triage_tier!r}; "
            f"triage_tier should only be set when triage_mode=True"
        )


@needs_fixture
def test_xic_apex_finds_parent_on_clean_file():
    """On a clean file the XIC apex should agree with the TIC apex within 0.1 min."""
    from core.mzml_parser import (
        parse_mzml_xic_peak, parse_mzml_chromatographic_peak,
    )
    from core.mass_calculator import calculate_monoisotopic_mass

    M = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)

    tic_apex, _, _, _, _, _ = parse_mzml_chromatographic_peak(SYNTAG_FILE, window_scans=3)
    xic_apex, _, _, _, _, _, xic_info = parse_mzml_xic_peak(
        SYNTAG_FILE, parent_mass=M, window_scans=3,
    )
    assert tic_apex is not None and xic_apex is not None
    assert not xic_info['fallback_to_tic'], \
        "XIC apex picker should find the parent on a clean SynTag file, not fall back"
    assert abs(tic_apex.rt_min - xic_apex.rt_min) < 0.1, \
        f"XIC apex RT {xic_apex.rt_min} differs from TIC apex {tic_apex.rt_min} by >0.1 min"


@needs_fixture
def test_xic_apex_handles_wrong_mass_gracefully():
    """A parent mass absent from the spectrum should fall back to TIC apex, not crash."""
    from core.mzml_parser import parse_mzml_xic_peak

    apex, window, _, _, _, _, xic_info = parse_mzml_xic_peak(
        SYNTAG_FILE, parent_mass=99999.99, window_scans=3,
    )
    assert apex is not None, "should have fallen back to TIC apex"
    assert xic_info['fallback_to_tic'], "should report fallback"
    assert xic_info['apex_intensity'] is None


@needs_fixture
def test_xic_isotope_envelope_required():
    """The XIC builder requires M0 + M+1 at the matched charge, rejecting lone noise spikes."""
    from core.mzml_parser import (
        parse_mzml_xic_peak, _build_xic_from_scan_slices,
    )
    import numpy as np

    fake_slice = (np.array([606.328]), np.array([1000000.0]))
    target_mzs = [(2, 606.328)]
    xic, hits = _build_xic_from_scan_slices(
        [fake_slice], target_mzs, tolerance_ppm=25.0,
    )
    assert xic[0] == 0.0, \
        f"single-peak (no isotope partner) should give 0 XIC; got {xic[0]}"

    fake_slice2 = (np.array([606.328, 606.830]), np.array([1000000.0, 500000.0]))
    xic2, hits2 = _build_xic_from_scan_slices(
        [fake_slice2], target_mzs, tolerance_ppm=25.0,
    )
    assert xic2[0] > 0, \
        f"M0 + M+1 peaks should give non-zero XIC; got {xic2[0]}"


@needs_fixture
def test_calibrate_from_anchors_handles_q_tof_drift():
    """XIC anchors fed via cal_anchors should recover the parent envelope under Q-TOF drift."""
    from core.mzml_parser import parse_mzml_xic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks
    from core.mass_calculator import calculate_monoisotopic_mass

    M = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)
    apex, window, _, _, _, _, xic_info = parse_mzml_xic_peak(
        SYNTAG_FILE, parent_mass=M, window_scans=3,
    )
    spec = coadd_centroid_spectra(window) if len(window) > 1 else apex
    pl = detect_peaks(spec, threshold_fraction=0.003)

    anchors = xic_info.get('apex_hits') or []
    assert len(anchors) >= 2, f"expected ≥2 cal anchors, got {len(anchors)}"

    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                     coupling_reagent='HATU', source_scans=window,
                     cal_anchors=anchors)
    assert mr.parent_envelope_fraction > 0.5, \
        f"with cal_anchors, expected parent envelope >50%, got {mr.parent_envelope_fraction*100:.1f}%"


@needs_fixture
def test_iemm_split_clean_synthesis_attributes_to_native():
    """On a clean synthesis IEMM should attribute >85% of the parent envelope to native form."""
    from core.mzml_parser import parse_mzml_xic_peak, coadd_centroid_spectra
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks
    from core.mass_calculator import calculate_monoisotopic_mass

    M = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)
    apex, window, _, _, _, _, xic_info = parse_mzml_xic_peak(
        SYNTAG_FILE, parent_mass=M, window_scans=3,
    )
    spec = coadd_centroid_spectra(window) if len(window) > 1 else apex
    pl = detect_peaks(spec, threshold_fraction=0.003)

    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                     coupling_reagent='HATU', source_scans=window,
                     cal_anchors=xic_info.get('apex_hits'))

    assert mr.iemm_applied, "IEMM should have been applied"
    if mr.parent_envelope_fraction > 0:
        native_share = mr.native_parent_fraction / mr.parent_envelope_fraction
        assert native_share > 0.85, \
            f"on clean synthesis, native share of parent should be >85%, got {native_share*100:.0f}%"


def test_iemm_decomposes_synthetic_mixture():
    """_iemm_deconvolve_charge should recover a synthetic 70:30 native:deamidated ratio."""
    from core.peak_matcher import _iemm_deconvolve_charge

    Tn = [1.0, 1.5, 1.4, 0.9, 0.4, 0.15, 0.05]
    Td = [0.0] + Tn[:-1]

    a_true, b_true = 70.0, 30.0
    O = [a_true * tn + b_true * td for tn, td in zip(Tn, Td)]

    a, b = _iemm_deconvolve_charge(O, Tn, Td)
    assert abs(a - a_true) < 1.0, f"native scale: expected ~{a_true}, got {a}"
    assert abs(b - b_true) < 1.0, f"deam scale: expected ~{b_true}, got {b}"


def test_iemm_handles_pure_native():
    """IEMM should report b=0 for a pure native envelope (no deamidation)."""
    from core.peak_matcher import _iemm_deconvolve_charge
    Tn = [1.0, 1.5, 1.4, 0.9, 0.4, 0.15, 0.05]
    Td = [0.0] + Tn[:-1]
    O = [50.0 * x for x in Tn]
    a, b = _iemm_deconvolve_charge(O, Tn, Td)
    assert abs(a - 50.0) < 0.5, f"native: expected ~50, got {a}"
    assert b < 0.1, f"deam should be ~0, got {b}"


def test_iemm_handles_pure_deamidated():
    """IEMM should report a=0 for a pure deamidated envelope."""
    from core.peak_matcher import _iemm_deconvolve_charge
    Tn = [1.0, 1.5, 1.4, 0.9, 0.4, 0.15, 0.05]
    Td = [0.0] + Tn[:-1]
    O = [40.0 * x for x in Td]
    a, b = _iemm_deconvolve_charge(O, Tn, Td)
    assert a < 0.1, f"native should be ~0, got {a}"
    assert abs(b - 40.0) < 0.5, f"deam: expected ~40, got {b}"


def test_peptide_elemental_formula_basic():
    """peptide_elemental_formula gives correct atom counts."""
    from core.mass_calculator import peptide_elemental_formula
    C, H, N, O, S = peptide_elemental_formula('G')
    assert (C, H, N, O, S) == (2, 5, 1, 2, 0), \
        f"glycine formula expected C2H5NO2, got C{C}H{H}N{N}O{O}S{S}"

    C, H, N, O, S = peptide_elemental_formula('A')
    assert (C, H, N, O, S) == (3, 7, 1, 2, 0), \
        f"alanine formula expected C3H7NO2, got C{C}H{H}N{N}O{O}S{S}"

    seq = 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'
    C, H, N, O, S = peptide_elemental_formula(seq, c_terminal_amide=True)
    assert 140 < C < 160 and 220 < H < 240, \
        f"GLP-1 amide: got C{C}H{H}N{N}O{O}S{S}"


def test_peptide_isotope_envelope_matches_averagine_for_typical_peptide():
    """Exact and averagine envelopes should agree within ~15% for a typical peptide."""
    from core.mass_calculator import peptide_isotope_envelope
    from core.peak_matcher import theoretical_isotope_envelope

    seq = 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'
    exact = peptide_isotope_envelope(seq, c_terminal_amide=True, n_iso=4)

    M = 3295.66
    averagine = theoretical_isotope_envelope(M, n_iso=4)

    assert abs(exact[1] / averagine[1] - 1.0) < 0.15, \
        f"exact M+1={exact[1]:.3f} differs >15% from averagine M+1={averagine[1]:.3f}"


@needs_fixture
def test_iterparse_does_not_leak_memory():
    """Repeated parses must not grow RSS unboundedly (iterparse parent-element leak)."""
    import resource
    from core.mzml_parser import parse_mzml_xic_peak
    from core.mass_calculator import calculate_monoisotopic_mass

    M = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)

    def rss_kb():
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    parse_mzml_xic_peak(SYNTAG_FILE, parent_mass=M, window_scans=3)
    rss_warm = rss_kb()

    for _ in range(3):
        parse_mzml_xic_peak(SYNTAG_FILE, parent_mass=M, window_scans=3)
    rss_after = rss_kb()

    growth_mb = (rss_after - rss_warm) / 1024
    assert growth_mb < 50, \
        f"3 parses grew RSS by {growth_mb:.1f} MB — possible memory leak in iterparse"


def test_strip_ns_cache_works():
    """_strip_ns should cache repeated lookups for the same tag."""
    from core.mzml_parser import _strip_ns, _NS_CACHE
    _NS_CACHE.clear()
    tag = '{http://psi.hupo.org/ms/mzml}spectrum'
    assert _strip_ns(tag) == 'spectrum'
    assert tag in _NS_CACHE
    assert _strip_ns(tag) == 'spectrum'


@needs_fixture
def test_pdf_builds():
    """The full pipeline plus PDF report should run without exceptions."""
    from core.mzml_parser import parse_mzml_apex
    from core.peak_detector import detect_peaks
    from core.peak_matcher import match_peaks
    from core.mass_calculator import (
        calculate_monoisotopic_mass, calculate_average_mass,
    )
    from core.batch_intelligence import get_root_cause
    from app.report import build_pdf

    spec, tic_vals, tic_rts, n_ms1, parse_warnings = parse_mzml_apex(SYNTAG_FILE)
    pl = detect_peaks(spec, threshold_fraction=0.003)
    mr = match_peaks(pl, SEQUENCE, tolerance_ppm=15.0, c_terminal_amide=True,
                     coupling_reagent="HATU")

    parent_mono = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)
    parent_avg  = calculate_average_mass(SEQUENCE, c_terminal_amide=True)

    enriched = []
    for m in mr.matches:
        _, why, fix = get_root_cause(m.impurity_type)
        enriched.append({
            "impurity":       m.impurity_name,
            "type":           m.impurity_type,
            "obs_mz":         m.observed_mz,
            "ppm":            m.ppm_error,
            "z":              m.charge_state,
            "intensity":      m.observed_intensity,
            "rel_int_pct":    round(100 * m.observed_intensity / max(pl.base_peak_int, 1), 2),
            "iso_confirmed":  m.isotope_confirmed,
            "envelope_score": m.envelope_score,
            "delta_da":       m.delta_mono,
            "predicted_mass": m.predicted_mass,
            "risk":           m.risk_level,
            "confidence":     m.confidence,
            "why":            why,
            "fix":            fix,
        })

    result = {
        "sequence":            SEQUENCE,
        "c_amide":             True,
        "fixed_mod_da":        0.0,
        "fname":               os.path.basename(SYNTAG_FILE),
        "analyst":             "smoke test",
        "tolerance_ppm":       15,
        "coupling":            "HATU",
        "parent_mono":         parent_mono,
        "parent_avg":          parent_avg,
        "scan_number":         spec.scan_number,
        "rt_min":              spec.rt_min,
        "tic_at_apex":         spec.total_ion_current,
        "n_analytical_scans":  n_ms1,
        "scan_window_used":    0,
        "n_scans_coadded":     1,
        "n_peaks":             len(pl.peaks),
        "base_peak_mz":        pl.base_peak_mz,
        "base_peak_int":       pl.base_peak_int,
        "matches":             enriched,
        "n_predicted":         mr.n_impurities_predicted,
        "n_matched":           mr.n_peaks_matched,
        "bear_case":           mr.bear_case_applied,
        "warnings":            list(parse_warnings),
        "parent_envelope_pct":      round(mr.parent_envelope_fraction * 100, 1),
        "matched_impurity_pct":     round(mr.matched_impurity_fraction * 100, 2),
        "unexplained_pct":          round(mr.unexplained_fraction * 100, 1),
        "n_unexplained_above_1pct": mr.n_unexplained_peaks_above_1pct,
        "top_unexplained":          mr.top_unexplained_peaks,
    }

    pdf = build_pdf(result)
    assert isinstance(pdf, (bytes, bytearray)), "build_pdf must return bytes"
    assert len(pdf) > 1000, f"PDF suspiciously small ({len(pdf)} bytes)"
    assert pdf[:4] == b"%PDF", "output is not a PDF"


if __name__ == "__main__":
    if SYNTAG_FILE is None:
        print(f"SKIP: {SKIP_REASON}")
        sys.exit(0)
    print("Running smoke tests directly (without pytest)…\n")
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"  PASS  {name}")
            except AssertionError as e:
                print(f"  FAIL  {name}: {e}")
            except Exception as e:
                print(f"  ERROR {name}: {type(e).__name__}: {e}")
    print("\nDone.")
