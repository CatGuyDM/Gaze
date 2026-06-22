"""Root-cause map for SPPS impurity types.

Given an impurity_type string (from impurity_engine), returns a tuple
(plain_english_name, why_it_forms, what_to_do). Used to populate the
report's "Likely cause" and "Recommended action" columns.
"""


_ROOT_CAUSE = {
    'deletion_single': (
        'Deletion peptide',
        'One amino acid coupling step failed or was skipped. The growing chain '
        'continued without incorporating that residue. Common causes: incomplete '
        'activation, resin swelling, or coupling reagent failure at a difficult '
        'position.',
        'Review coupling monitoring data at this position. Consider double-coupling '
        'or extended coupling time. Check for Pro or β-branched residues nearby.'
    ),
    'deletion_double': (
        'Double deletion peptide',
        'Two consecutive coupling steps failed. Indicates a systematic problem '
        'at that region of the sequence — likely aggregation or poor solvation.',
        'This region of the sequence may be aggregation-prone. Add pseudoproline '
        'dipeptides if Ser or Thr are in the vicinity. Consider chaotropic additives.'
    ),
    'truncation_n': (
        'N-terminal truncation',
        'The peptide chain terminated early, producing a shortened version missing '
        'residues from the N-terminus. Often caused by incomplete Fmoc deprotection '
        'at a specific cycle.',
        'Check Fmoc deprotection efficiency at early coupling cycles. Monitor UV '
        'absorbance at 301 nm during deprotection.'
    ),
    'truncation_c': (
        'C-terminal truncation',
        'Early chain termination at the C-terminus. Can result from incomplete '
        'loading of the first residue onto the resin.',
        'Verify resin loading and first residue attachment efficiency. Use HATU '
        'or double coupling for initial loading.'
    ),
    'aspartimide': (
        'Aspartimide',
        'Cyclization of an Asp residue with its backbone NH to form a succinimide '
        'ring (-18 Da). A classic Fmoc SPPS side reaction at Asp-Gly, Asp-Ser, '
        'Asp-Ala, and Asp-Thr sequences. Catalyzed by repeated piperidine treatment.',
        'Add 0.1 M HOBt to all piperidine deprotection steps. Consider using '
        'Dmb-protected Asp or pseudoproline dipeptides at risk positions. '
        'Minimize deprotection time and temperature.'
    ),
    'piperidine_aspartimide': (
        'Piperidine adduct on aspartimide',
        'After aspartimide forms, the piperidine used for Fmoc removal '
        'ring-opens the succinimide (+67 Da net). Both the aspartimide and '
        'this secondary product are often present together.',
        'Same intervention as aspartimide: 0.1 M HOBt in deprotection. '
        'If this product is detected, aspartimide was already forming before ring-opening.'
    ),
    'piperidine_adduct': (
        'N-terminal piperidine adduct',
        'Piperidine reacted with the free N-terminal amine (+84 Da). Rare '
        'in fresh crude; the dominant piperidine-derived adducts are '
        'aspartimide-piperidide (+67 Da, listed separately) and DBF-piperidine '
        '(+178 Da). When this candidate hits, verify carefully.',
        'Increase DMF wash steps after Fmoc deprotection. Verify wash volumes '
        'are sufficient. Check wash conductivity if monitoring inline.'
    ),
    'fmoc_residual': (
        'Incomplete Fmoc removal',
        'The Fmoc protecting group (+222 Da) was not fully removed from one '
        'position. That residue\'s nitrogen was blocked, halting chain elongation '
        'at the next coupling.',
        'Extend deprotection time or increase piperidine concentration. '
        'Double deprotection recommended for sequences following Pro residues.'
    ),
    'tbu_residual': (
        'Residual tBu protecting group',
        'The acid-labile tBu group (+56 Da per occurrence) on Asp, Glu, Ser, '
        'Thr, or Tyr was not fully removed during TFA cleavage. Often caused by '
        'insufficient TFA exposure time or scavenger depletion.',
        'Extend TFA cleavage time. Increase scavenger concentration (TIS, water, '
        'or thioanisole). For long or hydrophobic sequences, consider two-step cleavage.'
    ),
    'pbf_residual': (
        'Residual Pbf on Arg',
        'The Pbf protecting group (+254 Da) on Arg was not removed. Pbf is '
        'the slowest side-chain protecting group to cleave and requires extended '
        'TFA treatment.',
        'Extend TFA cleavage to at least 3 hours for Arg-containing sequences. '
        'Pbf removal is the rate-limiting step in cleavage for most Arg peptides. '
        'Ensure scavengers are fresh.'
    ),
    'trt_residual': (
        'Residual Trt protecting group',
        'Trt (+242 Da) on Cys, His, Asn, or Gln was not fully removed. '
        'Less common than tBu or Pbf retention but occurs with short TFA treatment.',
        'Standard TFA conditions (2-4h with scavengers) should be sufficient. '
        'If persistent, add 2.5% water and 2.5% TIS as scavengers.'
    ),
    'boc_lys_residual': (
        'Residual Boc on Lys',
        'The Boc group (+100 Da) on Lys ε-amine was not removed. Boc is '
        'acid-stable and requires strong TFA conditions to cleave from Lys.',
        'Ensure TFA concentration is sufficient (>90%). Extend cleavage time. '
        'Check that scavengers are not quenching TFA prematurely.'
    ),
    'met_oxidation': (
        'Met oxidation',
        'Met was oxidized to its sulfoxide (+16 Da). Caused by dissolved oxygen, '
        'peroxides in solvents/reagents, or oxidative conditions during synthesis '
        'or workup.',
        'Degas all solvents before synthesis. Add methionine or free-radical '
        'scavengers to cleavage cocktail. Store peptide under inert atmosphere. '
        'Use fresh reagents — peroxides accumulate in aged DMF and NMP.'
    ),
    'trp_oxidation': (
        'Trp oxidation',
        'Trp was oxidized to kynurenine (+4 Da) or di-oxidation product (+32 Da). '
        'Trp is the most oxidation-sensitive amino acid in peptides.',
        'Add 2% thioanisole to cleavage cocktail. Minimize UV exposure during '
        'purification. Store peptide protected from light.'
    ),
    'cys_oxidation': (
        'Cys oxidation',
        'Cys thiol was oxidized to sulfenic (+16), sulfinic (+32), or '
        'cysteic acid (+48). Often arises post-cleavage when the peptide '
        'is exposed to air after a reducing workup; can also result from '
        'transition-metal contamination.',
        'Reduce free Cys with TCEP just before LC-MS injection. Store '
        'Cys-containing peptides under N2 with EDTA in the storage buffer. '
        'Use degassed solvents during purification.'
    ),
    'deamidation': (
        'Deamidation',
        'Asn or Gln lost an amide group (-1 Da, net +0.984 Da). '
        'Asn-Gly is the fastest-deamidating sequence (t½ ~24h at pH 7.4). '
        'Also occurs under acidic or basic conditions during synthesis.',
        'Minimize exposure to neutral or basic pH. Lyophilize immediately after '
        'cleavage. For Asn-Gly sequences, consider using Gln instead of Asn.'
    ),
    'racemization': (
        'Racemization risk',
        'His, Cys, or Asp residues are prone to epimerization during activation '
        'and coupling. Cannot be detected by mass alone — requires chiral HPLC '
        'or NMR. Flagged here because the sequence contains high-risk residues.',
        'Use collidine as base instead of DIPEA for coupling steps involving His, '
        'Cys, or Asp. Minimize pre-activation time. Keep coupling temperature '
        'at or below 50°C.'
    ),
    'tfa_adduct': (
        'TFA covalent adduct',
        'TFA (+114 Da) formed a covalent adduct, most commonly at the N-terminus. '
        'Rare but observed at high TFA concentrations or with extended exposure.',
        'Ensure adequate post-cleavage workup. Precipitate peptide in cold ether '
        'to remove TFA. Multiple lyophilization cycles from water reduce TFA salts.'
    ),
    'disulfide_scrambling': (
        'Disulfide scrambling',
        'Two Cys residues formed a disulfide bond in an incorrect pairing. '
        'Identical mass to correct SS formation — distinguishable only by HPLC '
        'retention time or MS2.',
        'Use orthogonal Cys protection during synthesis (e.g. Cys(StBu) + Cys(Trt)). '
        'Perform directed disulfide formation post-deprotection under dilute, '
        'controlled oxidation conditions.'
    ),
    'phosphorylation': (
        'Phosphorylation',
        'Ser, Thr, or Tyr acquired a phosphate group (+80 Da). May indicate '
        'phosphorylating contaminants in reagents or unintended phosphate transfer.',
        'Check reagent purity. If unintended, review all buffers and reagents for '
        'phosphate contamination.'
    ),
    'hydroxylation': (
        'Pro hydroxylation',
        'Pro was converted to hydroxyproline (+16 Da) under oxidative conditions.',
        'Degas solvents. Reduce oxidant exposure during synthesis and workup.'
    ),
    'n_methylation': (
        'N-methylation',
        'N-methylation (+14 Da) of the backbone or side chain. Can arise from '
        'formaldehyde contamination in degraded DMF.',
        'Use fresh, high-purity DMF. Avoid using DMF that has been open >1 month. '
        'Consider switching to NMP which does not generate formaldehyde on degradation.'
    ),
    'ammonia_loss': (
        'Ammonia loss',
        'Loss of NH3 (-17 Da) from Lys, Arg, Gln, Asn, or N-terminus. '
        'Can indicate thermal degradation or harsh acidic conditions.',
        'Review cleavage and workup conditions. Avoid prolonged heating.'
    ),
    'dmf_adduct': (
        'DMF adduct',
        'N,N-dimethylformamide (+73 Da) reacted with the N-terminal amine. '
        'Caused by formylation from degraded DMF or insufficient washing.',
        'Use fresh DMF. Verify wash steps after coupling are complete. '
        'Consider switching to NMP for sensitive sequences.'
    ),
    'n_acetylation': (
        'N-terminal acetylation',
        'The N-terminus was acetylated (+42 Da). Can arise from acetonitrile '
        'in mobile phase, acetic acid in workup, or acetic anhydride contamination.',
        'Review all mobile phase components and workup reagents for acetic acid '
        'or acetonitrile exposure. Avoid acetonitrile contact during synthesis.'
    ),
    'guanylation': (
        'N-terminal guanylation',
        'HATU reacted with the free N-terminal amine to form a guanidine adduct '
        '(+42 Da). Occurs when excess HATU is present before the amino acid is '
        'added — pre-activation gone wrong.',
        'Add HATU to the amino acid solution first, then add base, then add '
        'to the resin. Do not pre-activate resin with HATU alone.'
    ),
    'dkp_deletion': (
        'DKP double deletion (diketopiperazine)',
        'During Fmoc deprotection of an X-Pro resin intermediate, the newly freed '
        'N-terminal amine attacks the adjacent amide bond intramolecularly, forming '
        'a 6-membered diketopiperazine ring. Both residues are lost from the chain '
        'simultaneously. This is NOT a coupling failure — it happens during deprotection. '
        'The desXaaXaa impurity has the same mass as a random double deletion but a '
        'completely different mechanism. Documented in tirzepatide manufacturing at Eli Lilly. '
        'Source: Wang et al. ACS Omega 2022;7(50):46809.',
        'Shorten Fmoc deprotection time. Couple the next residue immediately after '
        'deprotection — do not allow the deprotected X-Pro intermediate to age. '
        'Switch from 20% piperidine/DMF to 2% DBU/5% piperazine/NMP to reduce DKP. '
        'Consider Bsmoc protection instead of Fmoc at Pro residues. '
        'Keep deprotection temperature low. Use CTC resin if Pro is C-terminal.'
    ),
    'combination_aspartimide_piperidine': (
        'Aspartimide + piperidine ring-opening (combination)',
        'Two events occurred on the same molecule: aspartimide formed, then '
        'piperidine opened the ring. Both failures happened sequentially to '
        'the same peptide chain.',
        'This combination indicates a severe aspartimide problem. Use 0.1 M HOBt '
        'throughout. Consider switching to Asp(OtBu) pseudoamide protection.'
    ),
    'combination_oxidation_deamidation': (
        'Oxidation + deamidation (combination)',
        'Both Met oxidation and Asn/Gln deamidation occurred on the same molecule. '
        'Suggests oxidative stress conditions during synthesis or storage.',
        'Degas all solvents. Store under N2. Review storage conditions and timeline.'
    ),
    'combination_tbu_oxidation': (
        'Residual tBu + oxidation (combination)',
        'Incomplete deprotection and Met oxidation both occurred. Two separate '
        'failure modes on the same chain.',
        'Address both root causes: extend TFA cleavage and degas solvents.'
    ),
}

_DEFAULT_CAUSE = (
    'Synthesis impurity',
    'This impurity arises from an incomplete or side reaction during Fmoc SPPS.',
    'Review synthesis conditions at the affected position.'
)


def get_root_cause(impurity_type: str):
    for key in _ROOT_CAUSE:
        if key in impurity_type:
            return _ROOT_CAUSE[key]
    return _DEFAULT_CAUSE

