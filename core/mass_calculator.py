"""Mass calculations for peptides and SPPS impurities."""


PROTON = 1.00728
WATER     = 18.01056
WATER_AVG = 18.01528


RESIDUE_MASSES = {
    'A': 71.03711,  'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}

RESIDUE_MASSES_AVG = {
    'A': 71.0788,  'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
    'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
    'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
    'M': 131.1926, 'F': 147.1766, 'P': 97.1167,  'S': 87.0782,
    'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326,
}


RESIDUE_FORMULAS = {
    'A': (3,  5, 1, 1, 0), 'R': (6, 12, 4, 1, 0),
    'N': (4,  6, 2, 2, 0), 'D': (4,  5, 1, 3, 0),
    'C': (3,  5, 1, 1, 1), 'E': (5,  7, 1, 3, 0),
    'Q': (5,  8, 2, 2, 0), 'G': (2,  3, 1, 1, 0),
    'H': (6,  7, 3, 1, 0), 'I': (6, 11, 1, 1, 0),
    'L': (6, 11, 1, 1, 0), 'K': (6, 12, 2, 1, 0),
    'M': (5,  9, 1, 1, 1), 'F': (9,  9, 1, 1, 0),
    'P': (5,  7, 1, 1, 0), 'S': (3,  5, 1, 2, 0),
    'T': (4,  7, 1, 2, 0), 'W': (11, 10, 2, 1, 0),
    'Y': (9,  9, 1, 2, 0), 'V': (5,  9, 1, 1, 0),
}


def peptide_elemental_formula(sequence, c_terminal_amide=False, disulfide_bonds=0):
    NONSTD_TO_STD = {'Aib': 'A', 'Nle': 'L', 'Sar': 'G',
                      'Hyp': 'P', 'Orn': 'K', 'X':   'A'}
    C, H, N, O, S = 0, 0, 0, 0, 0
    seq = sequence.upper().strip()
    i = 0
    while i < len(seq):
        matched = False
        for ns, std in NONSTD_TO_STD.items():
            if seq[i:i+3].title() == ns and (i+3 == len(seq) or not seq[i+3].isalpha()):
                c, h, n, o, s = RESIDUE_FORMULAS[std]
                C += c; H += h; N += n; O += o; S += s
                i += 3
                matched = True
                break
        if matched:
            continue
        r = seq[i]
        if r in RESIDUE_FORMULAS:
            c, h, n, o, s = RESIDUE_FORMULAS[r]
            C += c; H += h; N += n; O += o; S += s
        i += 1

    H += 2; O += 1
    if c_terminal_amide:
        H += 1; N += 1; O -= 1
    if disulfide_bonds:
        H -= 2 * disulfide_bonds
    return C, H, N, O, S


def peptide_isotope_envelope(sequence, c_terminal_amide=False,
                              disulfide_bonds=0, n_iso=6):
    AB = {
        'C': [0.98930, 0.01070, 0.0, 0.0, 0.0],
        'H': [0.99989, 0.00011, 0.0, 0.0, 0.0],
        'N': [0.99636, 0.00364, 0.0, 0.0, 0.0],
        'O': [0.99757, 0.00038, 0.00205, 0.0, 0.0],
        'S': [0.94990, 0.00750, 0.04250, 0.0, 0.01],
    }
    counts = peptide_elemental_formula(
        sequence, c_terminal_amide=c_terminal_amide,
        disulfide_bonds=disulfide_bonds,
    )
    elements = ['C', 'H', 'N', 'O', 'S']
    envelope = [1.0] + [0.0] * n_iso

    for elem, n_atoms in zip(elements, counts):
        if n_atoms == 0:
            continue
        elem_dist = AB[elem][:n_iso + 1]
        elem_total = [1.0] + [0.0] * n_iso
        base = list(elem_dist) + [0.0] * (n_iso + 1 - len(elem_dist))
        n = n_atoms
        while n > 0:
            if n & 1:
                new = [0.0] * (n_iso + 1)
                for i in range(n_iso + 1):
                    if elem_total[i] == 0:
                        continue
                    for j in range(n_iso + 1 - i):
                        new[i + j] += elem_total[i] * base[j]
                elem_total = new
            n >>= 1
            if n > 0:
                new = [0.0] * (n_iso + 1)
                for i in range(n_iso + 1):
                    if base[i] == 0:
                        continue
                    for j in range(n_iso + 1 - i):
                        new[i + j] += base[i] * base[j]
                base = new

        new_env = [0.0] * (n_iso + 1)
        for i in range(n_iso + 1):
            if envelope[i] == 0:
                continue
            for j in range(n_iso + 1 - i):
                new_env[i + j] += envelope[i] * elem_total[j]
        envelope = new_env

    if envelope[0] > 0:
        envelope = [v / envelope[0] for v in envelope]
    return envelope


NON_STANDARD_RESIDUES = {
    'Aib': 85.05276, 'Sar': 57.02146, 'Nle': 113.08406,
    'Hyp': 113.04768, 'Orn': 114.07931,
    'X':   85.05276,
}

NON_STANDARD_RESIDUES_AVG = {
    'Aib': 85.1060, 'Sar': 57.0519, 'Nle': 113.1594,
    'Hyp': 113.1158, 'Orn': 114.1472,
    'X':   85.1060,
}


AMIDE_CORRECTION     = -0.98402
AMIDE_CORRECTION_AVG = -0.98400
DISULFIDE_H2         = -2.01565
DISULFIDE_H2_AVG     = -2.01600


IMPURITY_DELTAS = {
    'aspartimide':            -18.01056,
    'piperidine_aspartimide': +67.07859,
    'piperazine_aspartimide': +68.07384,

    'fmoc_residual':     +222.06808,
    'tbu_residual':       +56.06260,
    'pbf_residual':      +252.08197,
    'trt_residual':      +242.10955,
    'boc_lys_residual':  +100.05243,

    'piperidine_adduct': +84.08133,
    'dmf_adduct':        +73.05276,
    'n_formylation':     +27.99491,
    'n_acetylation':     +42.01057,

    'met_oxidation':       +15.99491,
    'trp_oxidation_mono':  +3.99491,
    'trp_oxidation_di':    +31.98983,

    'deamidation':         +0.98402,
    'pyroglutamate':       -17.02655,
    'water_loss':          -18.01056,
    'ammonia_loss':        -17.02655,

    'tfa_adduct':          +113.99239,

    'disulfide':           -2.01565,
    'disulfide_scrambling':  0.0,

    'phosphorylation':     +79.96633,
    'hydroxylation':       +15.99491,
    'n_methylation':       +14.01565,
}

IMPURITY_DELTAS_AVG = {
    'aspartimide':            -18.0153,
    'piperidine_aspartimide': +67.0902,
    'piperazine_aspartimide': +68.1206,
    'fmoc_residual':         +222.2400,
    'tbu_residual':           +56.1063,
    'pbf_residual':          +252.3168,
    'trt_residual':          +242.3163,
    'boc_lys_residual':      +100.1158,
    'piperidine_adduct':      +84.1404,
    'dmf_adduct':             +73.0938,
    'n_formylation':          +28.0101,
    'n_acetylation':          +42.0367,
    'met_oxidation':          +15.9994,
    'trp_oxidation_mono':      +3.9988,
    'trp_oxidation_di':       +31.9988,
    'deamidation':             +0.9847,
    'pyroglutamate':          -17.0305,
    'water_loss':             -18.0153,
    'ammonia_loss':           -17.0305,
    'tfa_adduct':            +114.0233,
    'disulfide':              -2.0160,
    'disulfide_scrambling':    0.0,
    'phosphorylation':        +79.9799,
    'hydroxylation':          +15.9994,
    'n_methylation':          +14.0269,
}


ASPARTIMIDE_HIGH_RISK   = set('GSAT')
ASPARTIMIDE_POSSIBLE    = set('NQHC')
RACEMIZATION_HIGH_RISK  = set('HCD')
TBU_PROTECTED           = set('DESTY')
TRT_PROTECTED           = {'C', 'H', 'N', 'Q'}


MEDBZ_LINKER_MONO = 148.06367
MEDBZ_LINKER_AVG  = 148.16
MENBZ_LINKER_MONO = MEDBZ_LINKER_MONO - WATER
MENBZ_LINKER_AVG  = MEDBZ_LINKER_AVG  - WATER_AVG
SYNTAG_C_AMIDE_MONO = MEDBZ_LINKER_MONO + 6 * RESIDUE_MASSES['R']
SYNTAG_C_AMIDE_AVG  = MEDBZ_LINKER_AVG  + 6 * RESIDUE_MASSES_AVG['R']
SYNTAG_MENBZ_MONO   = MENBZ_LINKER_MONO + 6 * RESIDUE_MASSES['R']
SYNTAG_MENBZ_AVG    = MENBZ_LINKER_AVG  + 6 * RESIDUE_MASSES_AVG['R']
ARGTAG_MONO = 6 * RESIDUE_MASSES['R']
ARGTAG_AVG  = 6 * RESIDUE_MASSES_AVG['R']


SEMAGLUTIDE_SEQUENCE    = 'HXEGTFTSDVSSYLEGQAAKEFIAWLVRGR'
SEMAGLUTIDE_LYS_MOD_DA  = 773.431
SEMAGLUTIDE_AVG_MW      = 4113.58
SEMAGLUTIDE_FULL_MONO   = 4111.115

OXYTOCIN_SEQUENCE       = 'CYIQNCPLG'
OXYTOCIN_MONO_PUBLISHED = 1006.4365
BRADYKININ_SEQUENCE     = 'RPPGFSPFR'


def calculate_monoisotopic_mass(
    sequence,
    c_terminal_amide=False,
    disulfide_bonds=0,
    fixed_modification_da=0.0,
    extra_residues=None,
):
    res = {**RESIDUE_MASSES,
           **(extra_residues if extra_residues is not None else NON_STANDARD_RESIDUES)}
    seq = sequence.upper().strip()
    if not seq:
        raise ValueError('Sequence cannot be empty.')
    unknown = sorted(set(aa for aa in seq if aa not in res))
    if unknown:
        raise ValueError(
            f'Unknown residue(s): {unknown}. '
            f'Pass extra_residues={{code: mass_da}} for non-standard amino acids.'
        )
    total = sum(res[aa] for aa in seq) + WATER
    if c_terminal_amide:
        total += AMIDE_CORRECTION
    if disulfide_bonds:
        total += disulfide_bonds * DISULFIDE_H2
    total += fixed_modification_da
    return round(total, 5)


def calculate_average_mass(
    sequence,
    c_terminal_amide=False,
    disulfide_bonds=0,
    fixed_modification_da=0.0,
    extra_residues=None,
):
    res = {**RESIDUE_MASSES_AVG,
           **(extra_residues if extra_residues is not None else NON_STANDARD_RESIDUES_AVG)}
    seq = sequence.upper().strip()
    if not seq:
        raise ValueError('Sequence cannot be empty.')
    unknown = sorted(set(aa for aa in seq if aa not in res))
    if unknown:
        raise ValueError(f'Unknown residue(s): {unknown}.')
    total = sum(res[aa] for aa in seq) + WATER_AVG
    if c_terminal_amide:
        total += AMIDE_CORRECTION_AVG
    if disulfide_bonds:
        total += disulfide_bonds * DISULFIDE_H2_AVG
    total += fixed_modification_da
    return round(total, 4)


def calculate_mz(neutral_mass, charge):
    if charge == 0:
        raise ValueError('Charge state cannot be 0.')
    return round((neutral_mass + charge * PROTON) / charge, 5)


def calculate_ppm_error(theoretical, observed):
    if theoretical == 0:
        return 999.0
    return abs((observed - theoretical) / theoretical) * 1e6


def calculate_all_charge_states(
    sequence, max_charge=6, c_terminal_amide=False,
    disulfide_bonds=0, fixed_modification_da=0.0, extra_residues=None,
):
    mass = calculate_monoisotopic_mass(
        sequence, c_terminal_amide=c_terminal_amide,
        disulfide_bonds=disulfide_bonds,
        fixed_modification_da=fixed_modification_da,
        extra_residues=extra_residues,
    )
    return {z: calculate_mz(mass, z) for z in range(1, max_charge + 1)}
