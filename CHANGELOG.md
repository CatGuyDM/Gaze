# Gaze CHANGELOG

## [1.2.0]
Chemistry derivations moved to `LOGIC.md`. No functional changes. I am now ready to make this a public repository :D (I want to start sending it out for testing) to be totally honest I think there is still a lot to do before this tool will be usefull, BUT I think I have made enough progress to start asking for feedback. I am not even in college yet, and thus I wont be a part of this tools target audience for another few years.

Next update will contain adjustments (or more thorough overhaul)
---

## [1.1.0]

bug fix lol, its late and I am too tired to talk about it. (this is a privite repo rn so its lowkey a diary)

---

## [1.0.0]

1.0 yippeeeee (i have no idea how I am supposed to number these :D )
 (`triage_mode=True` on `match_peaks`, `detect_peaks`,
`parse_mzml_xic_peak`).
- **`match_peaks(triage_mode=True)`**
- **`detect_peaks(triage_mode=True)`**: 
- **`parse_mzml_xic_peak(triage_mode=True)`**: 
- **`PeakMatch.triage_tier`** 

---

## [0.9.0]

Fixed some minor bugs i love you alllll <3 <3 <3 

- `compute_spectrum_coverage`: O(n²) Python loop → O(n log n) via
  `np.searchsorted` + `np.argpartition`
- `detect_peaks`: ~2× faster via `bisect` over Python lists
- `_decode_binary`: ~33% faster via `np.frombuffer`
- `_build_xic_from_scan_slices`: vectorized with `searchsorted`
- Memory leak fix in `_iterparse_spectra` 
- `compute_spectrum_coverage` was hardcoding `c_terminal_amide=False`;
  now threaded through from `match_peaks`

---

## [0.8.0]

IT IS 2 AM AND I AM TIRED BOSS D:

- **`parse_mzml_xic_peak()`**
- **`_calibrate_from_anchors()`**
- **`_two_stage_calibrate()`**
- **IEMM deamidation deconvolution** (`_iemm_deconvolve_charge`)
  `MatchResult` gains `iemm_applied`, `native_parent_fraction`,
  `deamidated_parent_fraction`.

---

## [0.7.0]

Filters to reduce overcalling.

- **Mass-degeneracy collapse** (`collapse_mass_degenerate` in
  `impurity_engine.py`): I merged candidates with identical monoisotopic mass into single rows with a position list
- **Scan persistence filter**: I decided matches must appear in ≥30% of input scans
  (`min_persistence_fraction`). `PeakMatch.scan_persistence` field added.
  Disabled when `triage_mode=True`.
- **`PeakMatch.position`** and **`PeakMatch.n_degenerate`** 
  every match

---

## [0.6.0]

6th commit, I did a few important things

- **Multi-scan integration**: `parse_mzml_chromatographic_peak()` and
  `coadd_centroid_spectra()` The idea behind this is that it co-adds 2N+1 scans around apex
- **Spectrum coverage metric**: `compute_spectrum_coverage()` 
- **Isotope envelope fitting**: `envelope_fit_score()` scores match quality
  0.0–1.0 against averagine theoretical envelope. `PeakMatch.envelope_score`
  added. I also maded it so low-scoring matches are penalised or rejected.

---

## [0.5.0]

- I did a bunch of little stupid shi, I found a few actual mzML files to test it on, so now i know that it "works" (i have no idea if its output is actually any good, I js know that there is A output) 
- I added a bunch of little settings stuff.

---

## [0.4.0]

- This is my fourth git commit. With the basic structure I am going to start the actual "work"

---

## [0.3.0]

- `core/impurity_engine.py`: enumeration engine yipeee :D
- Average mass support; m/z for z=1..6
- My friend said I should name this, and we decided on "Gaze" cuz your like "gazing" at the impurities (ik its stupid but im tired gng)

---

## [0.2.0]

- `core/mass_calculator.py`: 2nd commit, this is for monoisotopic mass, m/z, impurity deltas

---

## [0.1.0]

- its 1am and I am highkey tweaking rn, my readme/logic reference materials will be much better so dont worry, this commit is mostly js structure.
