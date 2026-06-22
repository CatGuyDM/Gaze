import warnings
import pandas as pd

from core.mass_calculator import (
    RESIDUE_MASSES, RESIDUE_MASSES_AVG, NON_STANDARD_RESIDUES,
    NON_STANDARD_RESIDUES_AVG, WATER, WATER_AVG, PROTON,
    IMPURITY_DELTAS, IMPURITY_DELTAS_AVG, ASPARTIMIDE_HIGH_RISK,
    ASPARTIMIDE_POSSIBLE, RACEMIZATION_HIGH_RISK, TBU_PROTECTED,
    TRT_PROTECTED, calculate_monoisotopic_mass, calculate_average_mass,
    calculate_mz,
)
from core.synthesis_conditions import SynthesisConditions

MAX_CHARGE = 6

DF_COLUMNS = [
    'impurity_id', 'impurity_name', 'impurity_type', 'position',
    'sequence_variant', 'delta_mono', 'delta_avg', 'mono_mass', 'avg_mass',
] + [f'mz_z{z}' for z in range(1, MAX_CHARGE + 1)] + [
    'risk_level', 'confidence', 'notes',
]


def _row(uid, name, itype, pos, variant, dmono, davg, parent_mono, parent_avg,
         risk, conf, notes):
    mono = round(parent_mono + dmono, 5)
    avg  = round(parent_avg  + davg,  4)
    row  = {
        'impurity_id': uid, 'impurity_name': name, 'impurity_type': itype,
        'position': str(pos), 'sequence_variant': variant,
        'delta_mono': round(dmono, 5), 'delta_avg': round(davg, 4),
        'mono_mass': mono, 'avg_mass': avg,
        'risk_level': risk, 'confidence': conf, 'notes': notes,
    }
    for z in range(1, MAX_CHARGE + 1):
        row[f'mz_z{z}'] = calculate_mz(mono, z)
    return row


def enumerate_impurities(
    sequence,
    c_terminal_amide=False,
    disulfide_bonds=0,
    fixed_modification_da=0.0,
    extra_residues=None,
    extra_residues_avg=None,
    max_single_deletions=True,
    max_double_deletions=True,
    max_truncation_length=3,
    include_protecting_groups=True,
    include_adducts=True,
    include_oxidation=True,
    include_deamidation=True,
    include_modified_residues=True,
    synthesis_conditions=None,
    tfa_min_intensity_fraction=0.005,
):
    if extra_residues is None:
        extra_residues = NON_STANDARD_RESIDUES
    if extra_residues_avg is None:
        extra_residues_avg = {}
        for code, mono_mass in extra_residues.items():
            if code in NON_STANDARD_RESIDUES_AVG:
                extra_residues_avg[code] = NON_STANDARD_RESIDUES_AVG[code]
            else:
                extra_residues_avg[code] = round(mono_mass * 1.0006, 4)

    seq = sequence.upper().strip()
    if not seq:
        raise ValueError('Sequence cannot be empty.')

    _combined_check = {**RESIDUE_MASSES, **(extra_residues or NON_STANDARD_RESIDUES)}
    _unknown = sorted(set(aa for aa in seq if aa not in _combined_check))
    if _unknown:
        warnings.warn(
            f'Unrecognized residue(s) in "{seq}": {_unknown}. '
            f'Treated as Gly. Pass extra_residues={{code: mass_da}}.',
            UserWarning, stacklevel=2
        )
        for _unk in _unknown:
            seq = seq.replace(_unk, 'G')

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

    rows = []
    uid  = [0]

    def add(name, itype, pos, variant, dmono, davg, risk, conf, notes):
        uid[0] += 1
        rows.append(_row(
            f'{itype}_{uid[0]:04d}', name, itype, pos, variant,
            dmono, davg, parent_mono, parent_avg, risk, conf, notes
        ))

    combined_mono = {**RESIDUE_MASSES, **extra_residues}
    combined_avg  = {**RESIDUE_MASSES_AVG, **extra_residues_avg}


    if max_single_deletions:
        for i, aa in enumerate(seq):
            if aa not in combined_mono:
                continue
            variant = seq[:i] + seq[i+1:]
            if not variant:
                continue
            try:
                dmono = -(combined_mono[aa] - 0)
                variant_mass = calculate_monoisotopic_mass(
                    variant, c_terminal_amide=c_terminal_amide,
                    disulfide_bonds=disulfide_bonds,
                    fixed_modification_da=fixed_modification_da,
                    extra_residues=extra_residues)
                dmono = variant_mass - parent_mono
                davg  = -(combined_avg.get(aa, combined_mono[aa]))
            except Exception:
                continue
            add(f'del_{aa}{i+1}', 'deletion_single', str(i+1), variant,
                dmono, davg, 'HIGH', 'HIGH',
                f'Single residue deletion at position {i+1} ({aa}).')


    _DKP_SECONDARY = {'P', 'X'}
    _DKP_MOTIFS    = {('A','P'),('G','P'),('V','P'),('P','P'),('S','P'),
                      ('R','P'),('K','P'),('L','P'),('I','P'),('F','P')}

    if max_double_deletions:
        for i in range(len(seq) - 1):
            aa1, aa2 = seq[i], seq[i+1]
            if aa1 not in combined_mono or aa2 not in combined_mono:
                continue
            variant = seq[:i] + seq[i+2:]
            if not variant:
                continue
            try:
                variant_mass = calculate_monoisotopic_mass(
                    variant, c_terminal_amide=c_terminal_amide,
                    disulfide_bonds=disulfide_bonds,
                    fixed_modification_da=fixed_modification_da,
                    extra_residues=extra_residues)
                dmono = variant_mass - parent_mono
                davg  = -(combined_avg.get(aa1, 0) + combined_avg.get(aa2, 0))
            except Exception:
                continue

            is_dkp = aa2 in _DKP_SECONDARY
            itype  = 'dkp_deletion' if is_dkp else 'deletion_double'
            note   = (f'DKP (diketopiperazine) double deletion at {i+1}-{i+2} ({aa1}{aa2}). '
                      f'X-Pro motif: Fmoc deprotection of the {aa1}{aa2}-resin intermediate '
                      f'triggers intramolecular cyclization releasing the dipeptide as a DKP. '
                      f'Wang et al. ACS Omega 2022.' if is_dkp else
                      f'Double deletion at {i+1}-{i+2} ({aa1}{aa2}).')
            add(f'del2_{aa1}{i+1}_{aa2}{i+2}', itype,
                f'{i+1}-{i+2}', variant, dmono, davg, 'HIGH', 'HIGH', note)


    for k in range(1, min(max_truncation_length + 1, len(seq))):
        variant = seq[k:]
        if not variant:
            continue
        try:
            variant_mass = calculate_monoisotopic_mass(
                variant, c_terminal_amide=c_terminal_amide,
                disulfide_bonds=0,
                fixed_modification_da=fixed_modification_da,
                extra_residues=extra_residues)
            dmono = variant_mass - parent_mono
            davg  = sum(-combined_avg.get(seq[j], combined_mono.get(seq[j], 0))
                        for j in range(k))
        except Exception:
            continue
        add(f'trunc_n_{k}aa', 'truncation_n', f'1-{k}', variant,
            dmono, davg, 'HIGH', 'HIGH', f'N-terminal truncation, {k} residue(s) removed.')


    for k in range(1, min(max_truncation_length + 1, len(seq))):
        variant = seq[:-k]
        if not variant:
            continue
        try:
            variant_mass = calculate_monoisotopic_mass(
                variant, c_terminal_amide=False,
                disulfide_bonds=0,
                fixed_modification_da=fixed_modification_da,
                extra_residues=extra_residues)
            dmono = variant_mass - parent_mono
            davg  = sum(-combined_avg.get(seq[-(j+1)], combined_mono.get(seq[-(j+1)], 0))
                        for j in range(k))
        except Exception:
            continue
        add(f'trunc_c_{k}aa', 'truncation_c', f'{len(seq)-k+1}-{len(seq)}', variant,
            dmono, davg, 'HIGH', 'HIGH', f'C-terminal truncation, {k} residue(s) removed.')


    for i in range(len(seq) - 1):
        if seq[i] == 'D':
            next_aa = seq[i + 1]
            risk = 'HIGH' if next_aa in ASPARTIMIDE_HIGH_RISK else (
                   'POSSIBLE' if next_aa in ASPARTIMIDE_POSSIBLE else None)
            if risk:
                dm = IMPURITY_DELTAS['aspartimide']
                da = IMPURITY_DELTAS_AVG['aspartimide']

                risk_level = risk if risk == 'HIGH' else 'MEDIUM'
                add(f'aspartimide_D{i+1}{next_aa}{i+2}', 'aspartimide',
                    f'{i+1}-{i+2}', '', dm, da, risk_level, risk,
                    f'Aspartimide at D{i+1}-{next_aa}{i+2}. Palasek 2007.')
                add(f'pip_aspartimide_D{i+1}{next_aa}{i+2}', 'piperidine_aspartimide',
                    f'{i+1}-{i+2}', '',
                    IMPURITY_DELTAS['piperidine_aspartimide'],
                    IMPURITY_DELTAS_AVG['piperidine_aspartimide'],
                    risk_level, risk,
                    f'Piperidine ring-opening of aspartimide at D{i+1}-{next_aa}{i+2}.')


    if include_protecting_groups:
        n_tbu = sum(1 for aa in seq if aa in TBU_PROTECTED)
        for n in range(1, min(n_tbu, 3) + 1):
            add(f'tbu_residual_x{n}', 'tbu_residual', '', '',
                IMPURITY_DELTAS['tbu_residual'] * n,
                IMPURITY_DELTAS_AVG['tbu_residual'] * n,
                'HIGH', 'HIGH', f'Residual tBu on {n} Asp/Glu/Ser/Thr/Tyr residue(s).')
        if 'R' in seq:
            n_arg = seq.count('R')
            for n in range(1, min(n_arg, 2) + 1):
                add(f'pbf_residual_x{n}', 'pbf_residual', '', '',
                    IMPURITY_DELTAS['pbf_residual'] * n,
                    IMPURITY_DELTAS_AVG['pbf_residual'] * n,
                    'HIGH', 'HIGH', f'Residual Pbf on {n} Arg residue(s).')
        trt_residues = [aa for aa in seq if aa in TRT_PROTECTED]
        if trt_residues:
            for n in range(1, min(len(trt_residues), 3) + 1):
                add(f'trt_residual_x{n}', 'trt_residual', '', '',
                    IMPURITY_DELTAS['trt_residual'] * n,
                    IMPURITY_DELTAS_AVG['trt_residual'] * n,
                    'HIGH', 'HIGH', f'Residual Trt on {n} Cys/His/Asn/Gln residue(s).')
        if 'K' in seq:
            add('boc_lys_residual', 'boc_lys_residual', '', '',
                IMPURITY_DELTAS['boc_lys_residual'],
                IMPURITY_DELTAS_AVG['boc_lys_residual'],
                'HIGH', 'HIGH', 'Residual Boc on Lys.')
        add('fmoc_incomplete_deprotection', 'fmoc_residual', '', '',
            IMPURITY_DELTAS['fmoc_residual'],
            IMPURITY_DELTAS_AVG['fmoc_residual'],
            'HIGH', 'HIGH', 'Incomplete Fmoc removal.')


    if include_adducts:
        add('piperidine_adduct_nterm', 'piperidine_adduct', '1', '',
            IMPURITY_DELTAS['piperidine_adduct'],
            IMPURITY_DELTAS_AVG['piperidine_adduct'],
            'LOW', 'MEDIUM',
            'N-terminal piperidine adduct (+84.081 Da). Rare in fresh crude; '
            'dominant piperidine-derived adducts are aspartimide-piperidide '
            '(+67) and DBF-piperidine (+178).')
        add('dmf_adduct_nterm', 'dmf_adduct', '1', '',
            IMPURITY_DELTAS['dmf_adduct'],
            IMPURITY_DELTAS_AVG['dmf_adduct'],
            'LOW', 'HIGH', 'N-terminal DMF adduct (+73.053 Da).')


    if include_oxidation:
        if 'M' in seq:
            add('met_oxidation', 'met_oxidation', '', '',
                IMPURITY_DELTAS['met_oxidation'],
                IMPURITY_DELTAS_AVG['met_oxidation'],
                'MEDIUM', 'HIGH', 'Met sulfoxide (+15.995 Da).')
        if 'W' in seq:
            add('trp_oxidation_mono', 'trp_oxidation', '', '',
                IMPURITY_DELTAS['trp_oxidation_mono'],
                IMPURITY_DELTAS_AVG['trp_oxidation_mono'],
                'MEDIUM', 'HIGH', 'Trp kynurenine oxidation (+4.000 Da).')
            add('trp_oxidation_di', 'trp_oxidation', '', '',
                IMPURITY_DELTAS['trp_oxidation_di'],
                IMPURITY_DELTAS_AVG['trp_oxidation_di'],
                'LOW', 'HIGH', 'Trp di-oxidation (+32.000 Da).')
        if 'C' in seq:


            add('cys_oxidation_mono', 'cys_oxidation', '', '',
                +15.99491, +15.9994,
                'MEDIUM', 'HIGH', 'Cys sulfenic acid (+15.995 Da).')
            add('cys_oxidation_di', 'cys_oxidation', '', '',
                +31.98983, +31.9988,
                'MEDIUM', 'HIGH', 'Cys sulfinic acid (+31.990 Da).')
            add('cys_oxidation_tri', 'cys_oxidation', '', '',
                +47.98474, +47.9982,
                'LOW', 'HIGH',
                'Cys cysteic acid (+47.985 Da). Often co-occurs with disulfide '
                'reduction during workup with reducing agents like TCEP.')


    if include_deamidation:
        if 'N' in seq or 'Q' in seq:
            add('deamidation', 'deamidation', '', '',
                IMPURITY_DELTAS['deamidation'],
                IMPURITY_DELTAS_AVG['deamidation'],
                'MEDIUM', 'HIGH', 'Deamidation Asn→Asp or Gln→Glu (+0.984 Da).')
        if seq[0] == 'Q':
            add('pyroglutamate', 'pyroglutamate', '1', '',
                IMPURITY_DELTAS['pyroglutamate'],
                IMPURITY_DELTAS_AVG['pyroglutamate'],
                'HIGH', 'HIGH', 'N-terminal pyroglutamate from Gln cyclization.')

    add('n_acetylation', 'n_acetylation', '1', '',
        IMPURITY_DELTAS['n_acetylation'],
        IMPURITY_DELTAS_AVG['n_acetylation'],
        'MEDIUM', 'HIGH', 'N-terminal acetylation (+42.011 Da).')


    rac = sorted({aa for aa in seq if aa in RACEMIZATION_HIGH_RISK})
    if rac:
        conds_known = synthesis_conditions is not None and synthesis_conditions.conditions_known
        note = (f'Racemization risk at {rac}. No mass change — HPLC chiral column required.'
                + ('' if conds_known else
                   ' sequence-based estimate only — add synthesis conditions for full assessment.'))
        add('racemization_risk', 'racemization', '', '', 0.0, 0.0, 'HIGH', 'HIGH', note)


    add('tfa_adduct_nterm', 'tfa_adduct', '1', '',
        IMPURITY_DELTAS['tfa_adduct'],
        IMPURITY_DELTAS_AVG['tfa_adduct'],
        'LOW', 'HIGH',
        f'TFA adduct (+113.998 Da). Min intensity gate: {tfa_min_intensity_fraction*100:.1f}%.')


    if seq.count('C') >= 4:
        add('disulfide_scrambling', 'disulfide_scrambling', '', '',
            IMPURITY_DELTAS['disulfide_scrambling'],
            IMPURITY_DELTAS_AVG.get('disulfide_scrambling', 0.0),
            'MEDIUM', 'HIGH',
            f'Disulfide scrambling ({seq.count("C")} Cys, ≥2 SS bonds). '
            f'No net mass change — distinguish from parent by RT and MS/MS only.')


    if include_modified_residues:
        phos = [aa for aa in seq if aa in 'STY']
        if phos:
            for n in range(1, min(len(phos), 3) + 1):
                add(f'phosphorylation_x{n}', 'phosphorylation', '', '',
                    IMPURITY_DELTAS['phosphorylation'] * n,
                    IMPURITY_DELTAS_AVG.get('phosphorylation', 79.979) * n,
                    'LOW', 'HIGH', f'Phosphorylation on {n} Ser/Thr/Tyr (+{79.966*n:.3f} Da).')
        if 'P' in seq:
            add('hydroxylation_pro', 'hydroxylation', '', '',
                IMPURITY_DELTAS['hydroxylation'],
                IMPURITY_DELTAS_AVG.get('hydroxylation', 15.999),
                'LOW', 'HIGH', 'Pro hydroxylation → Hyp (+15.995 Da).')
        add('n_methylation', 'n_methylation', '1', '',
            IMPURITY_DELTAS['n_methylation'],
            IMPURITY_DELTAS_AVG.get('n_methylation', 14.027),
            'LOW', 'HIGH', 'N-methylation (+14.016 Da). Can arise from formaldehyde in DMF.')
        if any(aa in seq for aa in 'KRQN'):
            add('ammonia_loss', 'ammonia_loss', '', '',
                IMPURITY_DELTAS['ammonia_loss'],
                IMPURITY_DELTAS_AVG.get('ammonia_loss', -17.027),
                'LOW', 'HIGH', 'Ammonia loss (-17.027 Da) from Lys/Arg/Gln/Asn or N-terminus.')

    if not rows:
        return pd.DataFrame(columns=DF_COLUMNS)

    df = pd.DataFrame(rows, columns=DF_COLUMNS)
    risk_order = {'HIGH': 0, 'POSSIBLE': 1, 'MEDIUM': 2, 'LOW': 3}
    df['_r'] = df['risk_level'].map(lambda x: risk_order.get(x, 9))
    df = df.sort_values(['_r', 'delta_mono']).drop(columns='_r').reset_index(drop=True)
    return df


def filter_by_type(df, impurity_type):
    return df[df['impurity_type'] == impurity_type].copy()


def filter_high_risk(df):
    return df[df['risk_level'] == 'HIGH'].copy()


def impurity_summary(df):
    counts = df['risk_level'].value_counts()
    parts = [f'{counts.get(r, 0)} {r}' for r in ['HIGH', 'MEDIUM', 'LOW']]
    return f'Total: {len(df)} impurities predicted: ' + ', '.join(parts)


def collapse_mass_degenerate(df, mass_tolerance_da=0.001):
    if len(df) == 0:
        return df.assign(n_degenerate=pd.Series(dtype=int))


    df = df.copy().reset_index(drop=True)
    df['_mass_bin'] = (df['mono_mass'] / mass_tolerance_da).round().astype(int)

    out_rows = []
    for (itype, mass_bin), group in df.groupby(['impurity_type', '_mass_bin'], sort=False):
        if len(group) == 1:
            row = group.iloc[0].to_dict()
            row['n_degenerate'] = 1
            out_rows.append(row)
            continue


        rep = group.iloc[0].to_dict()
        positions = [str(r['position']) for _, r in group.iterrows() if str(r.get('position', '')).strip()]

        names = [r['impurity_name'] for _, r in group.iterrows()]

        rep['impurity_name'] = f"{names[0]} (+{len(names)-1} mass-equiv)"

        try:
            sorted_pos = sorted(set(positions), key=lambda s: int(''.join(c for c in s if c.isdigit()) or 0))
        except Exception:
            sorted_pos = sorted(set(positions))
        rep['position'] = '/'.join(sorted_pos) if sorted_pos else rep.get('position', '')

        original_note = str(rep.get('notes', ''))
        alt_summary = ', '.join(names[:6]) + ('...' if len(names) > 6 else '')
        rep['notes'] = (
            f'{len(names)} mass-equivalent candidates: {alt_summary}. '
            f'Mass alone cannot distinguish; needs MS/MS or RT to localize. '
            f'{original_note}'
        )
        rep['n_degenerate'] = len(names)
        out_rows.append(rep)

    out = pd.DataFrame(out_rows).drop(columns=['_mass_bin'], errors='ignore')
    return out
