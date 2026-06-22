"""
Coupling-reagent–induced impurities.

The full impurity catalogue lives in mass_calculator.IMPURITY_DELTAS and is
enumerated by impurity_engine.enumerate_impurities. This module adds the
small set of impurities that depend on which activation reagent was used
(HATU, HBTU, DIC/Oxyma, PyBOP).

References for deltas:
  HATU N-terminal guanylation: El-Faham & Albericio, Chem Rev 2011
  HBTU uronium adduct:         El-Faham & Albericio, Chem Rev 2011
  DIC/Oxyma adduct:            Subiros-Funosas et al. Chem Eur J 2009
  PyBOP phosphonium adduct:    Coste et al. Tetrahedron Lett 1990
"""
import pandas as pd


COUPLING_REAGENT_IMPURITIES = {


    'HATU': [{'name': 'guanylation_nterm', 'type': 'guanylation',
              'delta_mono': 42.02180, 'delta_avg': 42.0399, 'risk': 'MEDIUM',
              'notes': 'N-terminal guanylation from HATU. El-Faham 2011.'}],
    'HATU/HOAt': [{'name': 'guanylation_nterm', 'type': 'guanylation',
                   'delta_mono': 42.02180, 'delta_avg': 42.0399, 'risk': 'MEDIUM',
                   'notes': 'N-terminal guanylation from HATU/HOAt. El-Faham 2011.'}],

    'HBTU': [{'name': 'uronium_adduct', 'type': 'uronium_adduct',
              'delta_mono': 57.02146, 'delta_avg': 57.0524, 'risk': 'MEDIUM',
              'notes': 'Uronium adduct from HBTU. El-Faham 2011.'}],

    'DIC':  [{'name': 'oxyma_adduct', 'type': 'oxyma_adduct',
              'delta_mono': 169.03821, 'delta_avg': 169.1351, 'risk': 'LOW',
              'notes': 'Oxyma adduct from DIC/Oxyma activation. Subiros-Funosas 2009.'}],

    'PyBOP': [{'name': 'phosphonium_adduct', 'type': 'phosphonium_adduct',
               'delta_mono': 77.02146, 'delta_avg': 77.0594, 'risk': 'LOW',
               'notes': 'Phosphonium adduct from PyBOP. Coste 1990.'}],
}


def get_coupling_reagent_impurities(parent_mono, parent_avg, coupling_reagent='HATU'):
    from core.impurity_engine import DF_COLUMNS, MAX_CHARGE
    from core.mass_calculator import calculate_mz
    if coupling_reagent not in COUPLING_REAGENT_IMPURITIES:
        return pd.DataFrame(columns=DF_COLUMNS)
    rows = []
    for i, imp in enumerate(COUPLING_REAGENT_IMPURITIES[coupling_reagent]):
        mono = round(parent_mono + imp['delta_mono'], 5)
        avg  = round(parent_avg  + imp['delta_avg'],  4)
        row  = {
            'impurity_id':    f'{imp["type"]}_cr_{i+1:04d}',
            'impurity_name':  imp['name'],
            'impurity_type':  imp['type'],
            'position':       '1',
            'sequence_variant': '',
            'delta_mono':     round(imp['delta_mono'], 5),
            'delta_avg':      round(imp['delta_avg'],  4),
            'mono_mass':      mono,
            'avg_mass':       avg,
            'risk_level':     imp['risk'],
            'confidence':     'HIGH',
            'notes':          imp['notes'],
        }
        for z in range(1, MAX_CHARGE + 1):
            row[f'mz_z{z}'] = calculate_mz(mono, z)
        rows.append(row)
    return pd.DataFrame(rows, columns=DF_COLUMNS)
