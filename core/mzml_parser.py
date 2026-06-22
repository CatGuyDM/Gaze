import os
import zlib
import base64
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Optional
import numpy as np

CV_MS_LEVEL         = 'MS:1000511'
CV_SCAN_START_TIME  = 'MS:1000016'
CV_MZ_ARRAY         = 'MS:1000514'
CV_INTENSITY_ARRAY  = 'MS:1000515'
CV_ZLIB_COMPRESSION = 'MS:1000574'
CV_NO_COMPRESSION   = 'MS:1000576'
CV_64BIT_FLOAT      = 'MS:1000523'
CV_32BIT_FLOAT      = 'MS:1000521'
CV_POSITIVE_SCAN    = 'MS:1000130'
CV_NEGATIVE_SCAN    = 'MS:1000129'
CV_TIC              = 'MS:1000285'
CV_BASE_PEAK_INT    = 'MS:1000505'


CV_CENTROID_SPECTRUM = 'MS:1000127'
CV_PROFILE_SPECTRUM  = 'MS:1000128'

@dataclass
class Spectrum:
    scan_index:       int
    scan_number:      int
    ms_level:         int
    rt_min:           Optional[float]
    mz_array:         np.ndarray
    intensity_array:  np.ndarray
    n_peaks:          int
    is_empty:         bool
    malformed:        bool = False
    malformed_reason: str  = ''
    polarity:         str  = 'positive'
    centroided:       bool = False
    function_id:      Optional[int] = None

    @property
    def base_peak_intensity(self):
        return float(self.intensity_array.max()) if not self.is_empty and len(self.intensity_array) else 0.0

    @property
    def total_ion_current(self):
        return float(self.intensity_array.sum()) if not self.is_empty and len(self.intensity_array) else 0.0


_NS_CACHE = {}

def _strip_ns(tag):
    cached = _NS_CACHE.get(tag)
    if cached is not None:
        return cached
    stripped = tag.split('}', 1)[1] if '}' in tag else tag
    _NS_CACHE[tag] = stripped
    return stripped


def _iterparse_spectra(file_path):
    context = ET.iterparse(file_path, events=('start', 'end'))


    root = None
    iterator = iter(context)
    for event, elem in iterator:
        if event == 'start':
            root = elem
            break


    for event, elem in iterator:
        if event != 'end':
            continue
        if _strip_ns(elem.tag) == 'spectrum':
            yield elem


            if root is not None:
                root.clear()


def _decode_binary(b64_str, compressed, bit_width):
    if not b64_str:
        return np.array([], dtype=np.float64)

    raw = base64.b64decode(b64_str)
    if compressed:
        raw = zlib.decompress(raw)
    if not raw:
        return np.array([], dtype=np.float64)


    if bit_width == 32:
        arr = np.frombuffer(raw, dtype=np.float32).astype(np.float64)
    else:
        arr = np.frombuffer(raw, dtype=np.float64).copy()
    return arr


def _parse_spectrum_elem(elem, scan_index):
    attrib = elem.attrib
    scan_id = attrib.get('id', f'scan={scan_index+1}')
    scan_number = scan_index + 1
    function_id = None
    for part in scan_id.split():
        if part.startswith('scan='):
            try:
                scan_number = int(part.split('=')[1])
            except ValueError:
                pass
        elif part.startswith('function='):
            try:
                function_id = int(part.split('=')[1])
            except ValueError:
                pass

    default_length = int(attrib.get('defaultArrayLength', 0))
    ms_level = 1
    rt_min = None
    polarity = 'positive'
    malformed = False
    malformed_reason = ''
    centroided  = False
    mz_arr = int_arr = None

    for child in elem.iter():
        tag = _strip_ns(child.tag)
        acc = child.attrib.get('accession', '')
        val = child.attrib.get('value', '')

        if tag == 'cvParam':
            if acc == CV_MS_LEVEL:
                try: ms_level = int(val)
                except ValueError: pass
            elif acc == CV_SCAN_START_TIME:
                try: rt_min = float(val)
                except ValueError: pass
            elif acc == CV_POSITIVE_SCAN:
                polarity = 'positive'
            elif acc == CV_NEGATIVE_SCAN:
                polarity = 'negative'
            elif acc == CV_CENTROID_SPECTRUM:
                centroided = True

        elif tag == 'binaryDataArray':
            is_mz = is_int = compressed = False
            bit_width = 64
            binary_text = None
            for cv in child:
                if _strip_ns(cv.tag) != 'cvParam':
                    continue
                a = cv.attrib.get('accession', '')
                if a == CV_MZ_ARRAY:       is_mz = True
                elif a == CV_INTENSITY_ARRAY: is_int = True
                elif a == CV_ZLIB_COMPRESSION: compressed = True
                elif a == CV_32BIT_FLOAT:  bit_width = 32
            for sub in child:
                if _strip_ns(sub.tag) == 'binary' and sub.text:
                    binary_text = sub.text
                    break
            if not (is_mz or is_int):
                continue
            if binary_text is None:
                if is_mz:
                    malformed, malformed_reason = True, f'No m/z binary in scan {scan_number}'
                continue
            try:
                arr = _decode_binary(binary_text, compressed, bit_width)
            except Exception as e:
                malformed, malformed_reason = True, f'Decode failed scan {scan_number}: {e}'
                arr = np.array([])
            if is_mz:   mz_arr = arr
            elif is_int: int_arr = arr

    if default_length == 0:
        return Spectrum(scan_index=scan_index, scan_number=scan_number,
                        ms_level=ms_level, rt_min=rt_min,
                        mz_array=np.array([]), intensity_array=np.array([]),
                        n_peaks=0, is_empty=True, polarity=polarity,
                        centroided=centroided, function_id=function_id)

    if mz_arr is None:
        mz_arr = np.array([])
        malformed = True
        malformed_reason = malformed_reason or f'No m/z array in scan {scan_number}'
    if int_arr is None:
        int_arr = np.zeros(len(mz_arr))
    if len(mz_arr) != len(int_arr):
        min_len = min(len(mz_arr), len(int_arr))
        mz_arr, int_arr = mz_arr[:min_len], int_arr[:min_len]
        malformed = True
        malformed_reason = f'Array length mismatch in scan {scan_number}'

    return Spectrum(scan_index=scan_index, scan_number=scan_number,
                    ms_level=ms_level, rt_min=rt_min,
                    mz_array=mz_arr, intensity_array=int_arr,
                    n_peaks=len(mz_arr), is_empty=len(mz_arr) == 0,
                    malformed=malformed, malformed_reason=malformed_reason,
                    polarity=polarity, centroided=centroided,
                    function_id=function_id)


def parse_mzml_apex(file_path, ms_level=1, analytical_function=1):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'File not found: {file_path}')

    best_spec   = None
    best_tic    = -1.0
    tic_list    = []
    rt_list     = []
    warnings    = []
    n_ms1       = 0
    n_skipped_function = 0
    scan_index  = 0
    have_seen_function_tag = False

    for elem in _iterparse_spectra(file_path):
        spec = _parse_spectrum_elem(elem, scan_index)
        elem.clear()
        scan_index += 1

        if spec.ms_level != ms_level:
            continue


        if spec.function_id is not None:
            have_seen_function_tag = True
            if spec.function_id != analytical_function:
                n_skipped_function += 1
                continue

        n_ms1 += 1
        tic   = spec.total_ion_current
        rt    = spec.rt_min if spec.rt_min is not None else float(n_ms1)
        tic_list.append(tic)
        rt_list.append(rt)

        if spec.malformed:
            warnings.append(f'Scan {spec.scan_number}: {spec.malformed_reason}')


        if not spec.is_empty and len(spec.mz_array) > 0:
            mz_max = float(spec.mz_array.max())
        else:
            mz_max = 0.0

        if tic > best_tic and not spec.is_empty and mz_max > 500.0:
            best_tic  = tic
            best_spec = spec

    if have_seen_function_tag and n_skipped_function > 0:
        warnings.append(
            f'Waters multi-function file: skipped {n_skipped_function} reference '
            f'channel scans (function != {analytical_function}). Used {n_ms1} '
            f'analytical scans for analysis.'
        )

    return best_spec, tic_list, rt_list, n_ms1, warnings


def parse_mzml_chromatographic_peak(
    file_path,
    ms_level=1,
    analytical_function=1,
    window_scans=3,
    rt_window_min=None,
):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'File not found: {file_path}')


    apex_scan_idx_in_pass = -1
    apex_rt              = None
    apex_tic             = -1.0
    analytical_scan_indices = []
    analytical_rts          = []
    n_skipped_function = 0
    have_seen_function_tag = False
    n_ms1 = 0
    scan_index = 0

    for elem in _iterparse_spectra(file_path):
        spec = _parse_spectrum_elem(elem, scan_index)
        elem.clear()
        scan_index += 1
        if spec.ms_level != ms_level:
            continue
        if spec.function_id is not None:
            have_seen_function_tag = True
            if spec.function_id != analytical_function:
                n_skipped_function += 1
                continue
        n_ms1 += 1
        analytical_scan_indices.append(spec.scan_index)
        analytical_rts.append(spec.rt_min if spec.rt_min is not None else float(n_ms1))
        if spec.is_empty or len(spec.mz_array) == 0:
            continue
        mz_max = float(spec.mz_array.max())
        tic    = spec.total_ion_current
        if tic > apex_tic and mz_max > 500.0:
            apex_tic = tic
            apex_scan_idx_in_pass = len(analytical_scan_indices) - 1
            apex_rt = spec.rt_min

    warnings = []
    if have_seen_function_tag and n_skipped_function > 0:
        warnings.append(
            f'Waters multi-function file: skipped {n_skipped_function} reference '
            f'channel scans (function != {analytical_function}). Used {n_ms1} '
            f'analytical scans for analysis.'
        )

    if apex_scan_idx_in_pass < 0:
        return None, [], [], analytical_rts, n_ms1, warnings


    if rt_window_min is not None and apex_rt is not None:
        keep_idx = [i for i, rt in enumerate(analytical_rts)
                    if abs(rt - apex_rt) <= rt_window_min]
    else:
        lo = max(0, apex_scan_idx_in_pass - window_scans)
        hi = min(len(analytical_scan_indices) - 1,
                 apex_scan_idx_in_pass + window_scans)
        keep_idx = list(range(lo, hi + 1))

    keep_file_indices = set(analytical_scan_indices[i] for i in keep_idx)


    window_spectra = []
    apex_spectrum  = None
    tic_list   = []
    scan_index = 0
    for elem in _iterparse_spectra(file_path):
        spec = _parse_spectrum_elem(elem, scan_index)
        elem.clear()
        scan_index += 1
        if spec.ms_level != ms_level:
            continue
        if spec.function_id is not None and spec.function_id != analytical_function:
            continue

        tic_list.append(spec.total_ion_current)
        if spec.scan_index in keep_file_indices:
            window_spectra.append(spec)

            apex_file_idx = analytical_scan_indices[apex_scan_idx_in_pass]
            if spec.scan_index == apex_file_idx:
                apex_spectrum = spec

    return (apex_spectrum, window_spectra, tic_list, analytical_rts,
            n_ms1, warnings)


def _build_xic_from_scan_slices(slices_per_scan, target_mzs, tolerance_ppm,
                                 charges_for_isotopes=None, n_isotopes=3):
    n_scans = len(slices_per_scan)
    xic     = np.zeros(n_scans, dtype=np.float64)
    hits    = [[] for _ in range(n_scans)]

    targets_have_charge = (
        len(target_mzs) > 0
        and isinstance(target_mzs[0], (tuple, list))
        and len(target_mzs[0]) == 2
    )

    C13 = 1.003355
    LOG2 = math.log(2.0)

    for i, (mzs, intens) in enumerate(slices_per_scan):
        if mzs is None or len(mzs) == 0:
            continue


        if len(mzs) > 1 and mzs[0] > mzs[-1]:
            ord_idx = np.argsort(mzs)
            mzs = mzs[ord_idx]; intens = intens[ord_idx]

        for tgt in target_mzs:
            if targets_have_charge:
                z, m0_theo = tgt
            else:
                z, m0_theo = None, tgt


            window_m0 = m0_theo * tolerance_ppm * 1e-6
            lo = np.searchsorted(mzs, m0_theo - window_m0, side='left')
            hi = np.searchsorted(mzs, m0_theo + window_m0, side='right')
            if lo >= hi:
                continue
            cand_idx_m0 = np.arange(lo, hi)

            if z is None:

                local = cand_idx_m0[np.argmax(intens[cand_idx_m0])]
                xic[i] += float(intens[local])
                hits[i].append((None, float(mzs[local]), float(intens[local]), 1))
                continue

            spacing = C13 / z


            neutral_mass_est    = m0_theo * z - z * 1.00728
            expected_m1_over_m0 = neutral_mass_est / 1800.0 if neutral_mass_est > 0 else 1.0

            best_m0_idx       = -1
            best_envelope_sum = 0.0
            best_n_iso        = 0
            best_score        = float('inf')

            for cand in cand_idx_m0:
                cand_mz       = float(mzs[cand])
                cand_int      = float(intens[cand])
                envelope_sum  = cand_int
                n_iso_seen    = 1
                m1_intensity  = 0.0
                for k in range(1, n_isotopes + 1):
                    target_mk = cand_mz + k * spacing
                    win_mk    = target_mk * tolerance_ppm * 1e-6
                    lo_k = np.searchsorted(mzs, target_mk - win_mk, side='left')
                    hi_k = np.searchsorted(mzs, target_mk + win_mk, side='right')
                    if lo_k < hi_k:

                        slice_mz = mzs[lo_k:hi_k]
                        best_local = lo_k + int(np.argmin(np.abs(slice_mz - target_mk)))
                        envelope_sum += float(intens[best_local])
                        n_iso_seen   += 1
                        if k == 1:
                            m1_intensity = float(intens[best_local])
                    else:
                        break
                if n_iso_seen < 2:
                    continue


                if cand_int > 0 and m1_intensity > 0 and expected_m1_over_m0 > 0:
                    obs_ratio   = m1_intensity / cand_int
                    ratio_score = abs(math.log(obs_ratio / expected_m1_over_m0) / LOG2)
                else:
                    ratio_score = 10.0

                if (n_iso_seen > best_n_iso) or (
                    n_iso_seen == best_n_iso and ratio_score < best_score
                ):
                    best_m0_idx       = int(cand)
                    best_envelope_sum = envelope_sum
                    best_n_iso        = n_iso_seen
                    best_score        = ratio_score

            if best_m0_idx < 0:
                continue

            obs_m0 = float(mzs[best_m0_idx])
            m0_int = float(intens[best_m0_idx])
            xic[i] += best_envelope_sum
            hits[i].append((z, obs_m0, m0_int, best_n_iso))

    return xic, hits


def _smooth_savgol_simple(y, window=5):
    if len(y) < window:
        return y.copy()
    half = window // 2
    out  = np.zeros_like(y, dtype=np.float64)
    for i in range(len(y)):
        lo = max(0, i - half)
        hi = min(len(y), i + half + 1)
        out[i] = np.mean(y[lo:hi])
    return out


def _xic_apex_index(xic, min_fwhm_scans=3, snr_threshold=3.0):
    if len(xic) < min_fwhm_scans:
        return None

    sm = _smooth_savgol_simple(xic, window=5)
    apex_idx = int(np.argmax(sm))
    apex_int = sm[apex_idx]
    if apex_int <= 0:
        return None


    sorted_sm = np.sort(sm)
    n_low     = max(5, len(sorted_sm) // 2)
    baseline  = float(np.median(sorted_sm[:n_low]))
    if baseline <= 0:


        baseline = float(np.percentile(sm[sm > 0], 50)) if (sm > 0).any() else 1.0

    if apex_int < snr_threshold * baseline:
        return None


    half_max = 0.5 * apex_int
    lo = apex_idx
    while lo > 0 and sm[lo - 1] >= half_max:
        lo -= 1
    hi = apex_idx
    while hi < len(sm) - 1 and sm[hi + 1] >= half_max:
        hi += 1
    fwhm_span = hi - lo + 1
    if fwhm_span < min_fwhm_scans:
        return None

    return apex_idx


def _xic_fwhm_scans(xic, apex_idx):
    if len(xic) < 3 or apex_idx is None:
        return 0
    sm = _smooth_savgol_simple(xic, window=5)
    apex_int = float(sm[apex_idx])
    if apex_int <= 0:
        return 0
    half_max = 0.5 * apex_int
    lo = apex_idx
    while lo > 0 and sm[lo - 1] >= half_max:
        lo -= 1
    hi = apex_idx
    while hi < len(sm) - 1 and sm[hi + 1] >= half_max:
        hi += 1
    return hi - lo + 1


def parse_mzml_xic_peak(
    file_path,
    parent_mass,
    charges=(2, 3, 4, 5, 6, 7, 8),
    ms_level=1,
    analytical_function=1,
    window_scans=3,
    rt_window_min=None,
    auto_window=True,
    progressive_ppm=(25.0, 50.0, 100.0, 200.0, 400.0),
    triage_mode=False,
    debug=False,
):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'File not found: {file_path}')

    target_mzs = [(z, (parent_mass + z * 1.00728) / z) for z in charges]
    max_ppm    = max(progressive_ppm)


    C13_iso = 1.003355
    slice_anchors = []
    for z, m0 in target_mzs:
        for k in range(0, 4):
            slice_anchors.append(m0 + k * C13_iso / z)


    analytical_scan_indices = []
    analytical_rts          = []
    analytical_tics         = []
    slices_per_scan         = []
    n_skipped_function      = 0
    have_seen_function_tag  = False
    n_ms1                   = 0
    scan_index              = 0

    for elem in _iterparse_spectra(file_path):
        spec = _parse_spectrum_elem(elem, scan_index)
        elem.clear()
        scan_index += 1
        if spec.ms_level != ms_level:
            continue
        if spec.function_id is not None:
            have_seen_function_tag = True
            if spec.function_id != analytical_function:
                n_skipped_function += 1
                continue
        n_ms1 += 1
        analytical_scan_indices.append(spec.scan_index)
        analytical_rts.append(spec.rt_min if spec.rt_min is not None else float(n_ms1))
        analytical_tics.append(spec.total_ion_current)


        if spec.is_empty or len(spec.mz_array) == 0:
            slices_per_scan.append((np.array([]), np.array([])))
            continue
        mzs  = spec.mz_array
        ints = spec.intensity_array
        keep_mask = np.zeros(len(mzs), dtype=bool)
        for anchor_mz in slice_anchors:
            window = anchor_mz * max_ppm * 1e-6
            lo = np.searchsorted(mzs, anchor_mz - window, side='left')
            hi = np.searchsorted(mzs, anchor_mz + window, side='right')
            if lo < hi:
                keep_mask[lo:hi] = True
        slices_per_scan.append((mzs[keep_mask], ints[keep_mask]))

    warnings = []
    if have_seen_function_tag and n_skipped_function > 0:
        warnings.append(
            f'Waters multi-function file: skipped {n_skipped_function} reference '
            f'channel scans (function != {analytical_function}). Used {n_ms1} '
            f'analytical scans for analysis.'
        )

    if n_ms1 == 0:
        return (None, [], [], analytical_rts, n_ms1, warnings,
                {'apex_rt': None, 'apex_intensity': None,
                 'tolerance_used': None, 'fallback_to_tic': False,
                 'xic_trace': np.array([]), 'target_mzs': target_mzs})


    apex_idx_in_pass = None
    tolerance_used   = None
    xic_trace        = None
    apex_hits        = None
    _best_hits       = []


    apex_snr  = 1.5 if triage_mode else 3.0
    apex_fwhm = 2   if triage_mode else 3

    per_charge_apexes = {}

    for tol in sorted(progressive_ppm):
        for z, m0_theo in target_mzs:
            if z in per_charge_apexes:
                continue
            xic_z, hits_z = _build_xic_from_scan_slices(
                slices_per_scan, [(z, m0_theo)], tolerance_ppm=tol,
            )
            cand = _xic_apex_index(xic_z, min_fwhm_scans=apex_fwhm, snr_threshold=apex_snr)
            if cand is not None:
                apex_int = float(xic_z[cand])
                n_iso = max((h[3] for h in hits_z[cand]), default=0) if cand < len(hits_z) and hits_z[cand] else 0
                per_charge_apexes[z] = {
                    'apex_idx': cand, 'score': apex_int * max(n_iso, 1),
                    'xic': xic_z, 'hits': hits_z, 'tol': tol,
                }
        if per_charge_apexes:
            break

    if per_charge_apexes:
        best_z = max(per_charge_apexes, key=lambda z: per_charge_apexes[z]['score'])
        best = per_charge_apexes[best_z]
        apex_idx_in_pass = best['apex_idx']
        xic_trace = best['xic']
        tolerance_used = best['tol']
        _best_hits = best['hits']
    else:
        apex_idx_in_pass = None
        xic_trace = None
        tolerance_used = None
        _best_hits = []

    if debug:
        import sys
        apex_rt_dbg = analytical_rts[apex_idx_in_pass] if apex_idx_in_pass is not None else None
        print(f"[debug-apex] best_z={best_z if per_charge_apexes else None}, "
              f"apex_idx={apex_idx_in_pass}, apex_rt={apex_rt_dbg}, tol={tolerance_used}", file=sys.stderr)
        if apex_idx_in_pass is not None:
            print(f"[debug-apex] target_mzs={target_mzs}", file=sys.stderr)
            win_lo = max(0, apex_idx_in_pass - 3)
            win_hi = min(len(slices_per_scan) - 1, apex_idx_in_pass + 3)
            for si in range(win_lo, win_hi + 1):
                mzs, ints = slices_per_scan[si]
                if mzs is not None and len(mzs) > 0:
                    top_n = min(20, len(mzs))
                    top_idx = np.argpartition(-ints, top_n - 1)[:top_n]
                    top_idx = top_idx[np.argsort(-ints[top_idx])]
                    print(f"[debug-apex]  scan_offset={si-apex_idx_in_pass}: "
                          f"{list(zip(mzs[top_idx].round(4), ints[top_idx].round(0)))}", file=sys.stderr)


    if apex_idx_in_pass is not None:
        widest_tol = max(progressive_ppm)
        _, all_hits = _build_xic_from_scan_slices(
            slices_per_scan, target_mzs, tolerance_ppm=widest_tol,
        )
        window_lo = max(0, apex_idx_in_pass - window_scans)
        window_hi = min(len(all_hits) - 1, apex_idx_in_pass + window_scans)
        apex_hits = []
        for i in range(window_lo, window_hi + 1):
            apex_hits.extend(all_hits[i])

    fallback = False
    if apex_idx_in_pass is None:


        fallback = True
        warnings.append(
            f'XIC apex not found at any tolerance up to {max(progressive_ppm)} ppm. '
            f'Falling back to TIC apex.'
        )


        if not analytical_tics:
            return (None, [], [], analytical_rts, n_ms1, warnings,
                    {'apex_rt': None, 'apex_intensity': None,
                     'tolerance_used': None, 'fallback_to_tic': True,
                     'xic_trace': np.array([]), 'target_mzs': target_mzs})
        apex_idx_in_pass = int(np.argmax(analytical_tics))

    apex_rt        = analytical_rts[apex_idx_in_pass]
    apex_intensity = (None if fallback else float(xic_trace[apex_idx_in_pass]))

    fwhm = None
    if auto_window and xic_trace is not None and apex_idx_in_pass is not None:
        fwhm = _xic_fwhm_scans(xic_trace, apex_idx_in_pass)
        if fwhm >= 3:
            _half = max(2, min(12, int(round(fwhm * 0.5))))
        else:
            _half = window_scans
        effective_window = _half
    else:
        effective_window = window_scans

    if rt_window_min is not None and apex_rt is not None:
        keep_idx = [i for i, rt in enumerate(analytical_rts)
                    if abs(rt - apex_rt) <= rt_window_min]
    else:
        lo = max(0, apex_idx_in_pass - effective_window)
        hi = min(len(analytical_scan_indices) - 1,
                 apex_idx_in_pass + effective_window)
        keep_idx = list(range(lo, hi + 1))

    keep_file_indices = set(analytical_scan_indices[i] for i in keep_idx)
    apex_file_idx     = analytical_scan_indices[apex_idx_in_pass]


    window_spectra = []
    apex_spectrum  = None
    tic_list       = []
    scan_index     = 0
    for elem in _iterparse_spectra(file_path):
        spec = _parse_spectrum_elem(elem, scan_index)
        elem.clear()
        scan_index += 1
        if spec.ms_level != ms_level:
            continue
        if spec.function_id is not None and spec.function_id != analytical_function:
            continue
        tic_list.append(spec.total_ion_current)
        if spec.scan_index in keep_file_indices:
            window_spectra.append(spec)
            if spec.scan_index == apex_file_idx:
                apex_spectrum = spec

    xic_info = {
        'apex_rt':         apex_rt,
        'apex_intensity':  apex_intensity,
        'tolerance_used':  tolerance_used,
        'fallback_to_tic': fallback,
        'xic_trace':       xic_trace if xic_trace is not None else np.array([]),
        'target_mzs':      target_mzs,
        'apex_hits':       apex_hits if apex_hits is not None else [],
        'auto_window':     auto_window,
        'fwhm_scans':      fwhm if (auto_window and xic_trace is not None) else None,
        'effective_window': effective_window,
    }

    return (apex_spectrum, window_spectra, tic_list, analytical_rts,
            n_ms1, warnings, xic_info)


def coadd_centroid_spectra(spectra, mz_tolerance_ppm=20.0):
    if not spectra:
        return None
    if len(spectra) == 1:
        return spectra[0]


    all_mz  = np.concatenate([s.mz_array        for s in spectra if not s.is_empty])
    all_int = np.concatenate([s.intensity_array for s in spectra if not s.is_empty])

    if len(all_mz) == 0:

        return spectra[0]


    order   = np.argsort(all_mz)
    all_mz  = all_mz[order]
    all_int = all_int[order]


    out_mz  = []
    out_int = []
    cluster_mz_sum  = all_mz[0]  * all_int[0]
    cluster_int_sum = all_int[0]
    cluster_anchor  = all_mz[0]

    for i in range(1, len(all_mz)):
        ppm_gap = (all_mz[i] - cluster_anchor) / cluster_anchor * 1e6
        if ppm_gap < mz_tolerance_ppm:
            cluster_mz_sum  += all_mz[i] * all_int[i]
            cluster_int_sum += all_int[i]


        else:
            if cluster_int_sum > 0:
                out_mz.append(cluster_mz_sum / cluster_int_sum)
                out_int.append(cluster_int_sum)
            cluster_mz_sum  = all_mz[i] * all_int[i]
            cluster_int_sum = all_int[i]
            cluster_anchor  = all_mz[i]

    if cluster_int_sum > 0:
        out_mz.append(cluster_mz_sum / cluster_int_sum)
        out_int.append(cluster_int_sum)

    out_mz  = np.array(out_mz,  dtype=float)
    out_int = np.array(out_int, dtype=float)


    apex = max(spectra, key=lambda s: s.total_ion_current)
    return Spectrum(
        scan_index       = apex.scan_index,
        scan_number      = apex.scan_number,
        ms_level         = apex.ms_level,
        rt_min           = apex.rt_min,
        mz_array         = out_mz,
        intensity_array  = out_int,
        n_peaks          = len(out_mz),
        is_empty         = len(out_mz) == 0,
        malformed        = False,
        malformed_reason = f'Co-added from {len(spectra)} scans (multi-scan integration)',
        polarity         = apex.polarity,
        centroided       = apex.centroided,
        function_id      = apex.function_id,
    )

