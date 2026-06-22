# Gaze

Gaze matches peaks in an LC-MS spectrum to a list of predicted peptide
synthesis impurities. This file explains the logic and lists the impurities
it knows about.

## The steps

1. Read the mzML file and build an extracted-ion chromatogram (XIC) for the
   peptide across charge states z=2 to 8.
2. Pick the apex (strongest scan) and average a few neighboring scans.
3. Detect the peaks in that averaged spectrum and guess their charge states.
4. List the impurities expected for the input sequence.
5. Correct small m/z errors using the peptide's own isotope peaks.
6. Match impurities to peaks by mass and isotope pattern.
7. In default mode, only keep matches that show up in several scans.
8. Add up parent vs impurity vs unexplained signal and write a PDF.

## Masses

Monoisotopic atomic masses used:

```
C = 12.0000000
H = 1.0078250
N = 14.0030740
O = 15.9949146
S = 31.9720707
P = 30.9737616
```

A peptide's mass is the sum of its residue masses plus one water. A C-terminal
amide subtracts 0.984 Da. Each disulfide bond subtracts 2 H.

## Impurities

Each impurity has a fixed mass shift from the real peptide. Gaze checks for
these:

| Impurity | Mass shift (Da) | What it is |
|---|---|---|
| aspartimide | -18.011 | Asp cyclizes and loses water |
| deamidation | +0.984 | Asn/Gln turns into Asp/Glu |
| water loss | -18.011 | loss of H2O |
| ammonia loss | -17.027 | loss of NH3 |
| pyroglutamate | -17.027 | N-terminal Gln cyclizes |
| Met oxidation | +15.995 | oxygen added to methionine |
| Trp oxidation | +3.995 / +31.990 | mono / di oxidation |
| Cys oxidation | +15.995 / +31.990 / +47.985 | sulfenic / sulfinic / cysteic acid |
| phosphorylation | +79.966 | +HPO3 |
| hydroxylation | +15.995 | extra oxygen |
| N-methylation | +14.016 | extra CH3 |
| N-formylation | +27.995 | extra HCO |
| N-acetylation | +42.011 | extra COCH3 |
| Fmoc residual | +222.068 | leftover Fmoc protecting group |
| tBu residual | +56.063 | leftover tBu group |
| Pbf residual | +252.082 | leftover Pbf group (on Arg) |
| Trt residual | +242.110 | leftover Trt group |
| Boc residual | +100.052 | leftover Boc group (on Lys) |
| TFA adduct | +113.992 | TFA stuck to the peptide |
| DMF adduct | +73.053 | DMF added |
| piperidine adduct | +84.081 | piperidine added |
| disulfide | -2.016 | one S-S bond forms (loses 2 H) |

## Charge states

Bigger peptides pick up more charges (about one per 500-1000 Da), so Gaze
scales the number of charge states it checks with peptide size. When two charge
states could explain the same peaks, it picks the one whose isotope pattern
matches the expected shape.

## Overlapping peaks

The real peptide and its deamidated form differ by only 0.984 Da, so their
isotope patterns overlap. Gaze separates them by fitting the observed pattern
as a mix of the two expected patterns.
