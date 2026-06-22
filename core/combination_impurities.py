"""Combination impurities — second-order events (e.g. aspartimide + piperidide)
that occur on the same peptide molecule. Each rule is gated by a sequence
predicate so combinations only fire when the chemistry can occur."""
import pandas as pd
from core.mass_calculator import (
    IMPURITY_DELTAS, IMPURITY_DELTAS_AVG,
    calculate_monoisotopic_mass, calculate_average_mass,
    calculate_mz, NON_STANDARD_RESIDUES,
)

MAX_CHARGE = 6


def _has_asp_risk(seq): return any(seq[i]=='D' and seq[i+1] in 'GSAT' for i in range(len(seq)-1))
def _has_met(seq):      return 'M' in seq
def _has_trp(seq):      return 'W' in seq
def _has_asn_gln(seq):  return 'N' in seq or 'Q' in seq
def _has_arg(seq):      return 'R' in seq
def _has_2cys(seq):     return seq.count('C') >= 2
def _has_tbu(seq):      return any(aa in 'DEST' for aa in seq)


COMBINATION_RULES = [
    {'name': 'aspartimide_then_piperidine', 'type': 'combination_aspartimide_piperidine',
     'mono_components': ['aspartimide', 'piperidine_aspartimide'],
     'risk': 'HIGH',
     'condition': _has_asp_risk,
     'notes': 'Aspartimide + piperidine ring-opening. Sequential events. Palasek 2007.'},
    {'name': 'met_oxidation_plus_deamidation', 'type': 'combination_oxidation_deamidation',
     'mono_components': ['met_oxidation', 'deamidation'],
     'risk': 'MEDIUM',
     'condition': lambda s: _has_met(s) and _has_asn_gln(s),
     'notes': 'Met oxidation + deamidation under oxidative conditions.'},
    {'name': 'trp_oxidation_plus_met_oxidation', 'type': 'combination_double_oxidation',
     'mono_components': ['trp_oxidation_mono', 'met_oxidation'],
     'risk': 'MEDIUM',
     'condition': lambda s: _has_met(s) and _has_trp(s),
     'notes': 'Trp kynurenine + Met sulfoxide in same molecule.'},
    {'name': 'disulfide_plus_deamidation', 'type': 'combination_disulfide_deamidation',
     'mono_components': ['disulfide', 'deamidation'],
     'risk': 'MEDIUM',
     'condition': lambda s: _has_2cys(s) and _has_asn_gln(s),
     'notes': 'SS bond + deamidation on same molecule.'},
    {'name': 'tbu_plus_met_oxidation', 'type': 'combination_tbu_oxidation',
     'mono_components': ['tbu_residual', 'met_oxidation'],
     'risk': 'HIGH',
     'condition': lambda s: _has_met(s) and _has_tbu(s),
     'notes': 'Residual tBu + Met oxidation.'},
    {'name': 'pbf_plus_deamidation', 'type': 'combination_pbf_deamidation',
     'mono_components': ['pbf_residual', 'deamidation'],
     'risk': 'HIGH',
     'condition': lambda s: _has_arg(s) and _has_asn_gln(s),
     'notes': 'Residual Pbf on Arg + deamidation. Both occur during long TFA exposure.'},
    {'name': 'acetylation_plus_met_oxidation', 'type': 'combination_acetylation_oxidation',
     'mono_components': ['n_acetylation', 'met_oxidation'],
     'risk': 'MEDIUM',
     'condition': _has_met,
     'notes': 'N-terminal acetylation + Met oxidation.'},
]


def enumerate_combination_impurities(
    sequence, c_terminal_amide=False, disulfide_bonds=0,
    fixed_modification_da=0.0, extra_residues=None,
    extra_residues_avg=None,
):
    if extra_residues is None:
        extra_residues = NON_STANDARD_RESIDUES
    if extra_residues_avg is None:
        extra_residues_avg = extra_residues

    seq = sequence.upper().strip()
    parent_mono = calculate_monoisotopic_mass(
        seq, c_terminal_amide=c_terminal_amide,
        disulfide_bonds=disulfide_bonds,
        fixed_modification_da=fixed_modification_da,
        extra_residues=extra_residues,
    )
    parent_avg = calculate_average_mass(
        seq, c_terminal_amide=c_terminal_amide,
        disulfide_bonds=disulfide_bonds,
        fixed_modification_da=fixed_modification_da,
        extra_residues=extra_residues_avg,
    )


    from core.impurity_engine import DF_COLUMNS

    rows = []
    for i, rule in enumerate(COMBINATION_RULES):
        if not rule['condition'](seq):
            continue
        dm = sum(IMPURITY_DELTAS[c] for c in rule['mono_components'])
        da = sum(IMPURITY_DELTAS_AVG.get(c, IMPURITY_DELTAS[c]) for c in rule['mono_components'])
        mono = round(parent_mono + dm, 5)
        avg  = round(parent_avg  + da, 4)
        row  = {'impurity_id': f'{rule["type"]}_{i+1:04d}',
                'impurity_name': rule['name'], 'impurity_type': rule['type'],
                'position': '', 'sequence_variant': '',
                'delta_mono': round(dm, 5), 'delta_avg': round(da, 4),
                'mono_mass': mono, 'avg_mass': avg,
                'risk_level': rule['risk'], 'confidence': 'HIGH',
                'notes': f'[COMBINATION] {rule["notes"]} '
                         f'Components: {rule["mono_components"]}'}
        for z in range(1, MAX_CHARGE + 1):
            row[f'mz_z{z}'] = calculate_mz(mono, z)
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=DF_COLUMNS)
    return pd.DataFrame(rows, columns=DF_COLUMNS)
