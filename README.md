# Gaze

Gaze finds SPPS (solid-phase peptide synthesis) impurities in LC-MS data.

You give it a peptide sequence and an mzML file. It looks for common synthesis
byproducts by matching the masses of the peaks in the spectrum against a list
of predicted impurities, then shows how much of the signal is the real peptide
versus impurities.

## Process

1. Find the peptide in the data across its charge states and pick the best scan
2. Detect the peaks in that scan
3. Build a list of likely impurities for the sequence
4. Match impurities to peaks by mass and isotope pattern
5. Add up how much signal is parent vs impurity vs unexplained
6. Make a PDF report

There are two modes: a default mode that only keeps confident matches, and a
triage mode that returns everything above the noise floor (sorted into A/B/C
tiers) for when you don't want to miss anything.

## Install

```bash
git clone https://github.com/CatGuyDM/gaze.git
cd gaze
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

Needs Python 3.10+.

## Run

```bash
streamlit run streamlit_app.py
```

Open http://localhost:8501, enter a sequence, pick a C-terminal tag, upload an
mzML file, and click Run.

## Input

mzML files with centroided MS1 data (use MSConvert to peak-pick).

See LOGIC.md for the impurity list and how the matching works.
