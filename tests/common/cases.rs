#![allow(non_upper_case_globals)]

mod graph;
mod molecules;

pub use graph::TestMolecule;
use pauling::{AtomId, BondId, ResonanceSystem};

pub struct ResonanceSystemExpectation {
    pub atoms: fn() -> Vec<AtomId>,
    pub bonds: fn() -> Vec<BondId>,
}

impl ResonanceSystemExpectation {
    pub fn to_system(&self) -> ResonanceSystem {
        let atoms = (self.atoms)();
        let bonds = (self.bonds)();
        ResonanceSystem::new(atoms, bonds)
    }
}

pub struct ResonanceCase {
    pub slug: &'static str,
    pub title: &'static str,
    pub build: fn() -> TestMolecule,
    pub expected: &'static [ResonanceSystemExpectation],
}

impl ResonanceCase {
    pub fn expected_systems(&self) -> Vec<ResonanceSystem> {
        self.expected
            .iter()
            .map(ResonanceSystemExpectation::to_system)
            .collect()
    }
}

macro_rules! resonance_cases {
    ( $( $name:ident => {
            title: $title:expr,
            molecule: $builder:path,
            systems: [ $( { atoms: $atoms:expr, bonds: $bonds:expr } ),* $(,)? ],
        }
    ),* $(,)? ) => {
        $(
            pub const $name: ResonanceCase = ResonanceCase {
                slug: stringify!($name),
                title: $title,
                build: $builder,
                expected: &[
                    $(
                        ResonanceSystemExpectation {
                            atoms: $atoms,
                            bonds: $bonds,
                        }
                    ),*
                ],
            };
        )*

        pub const ALL_CASES: &[ResonanceCase] = &[
            $( $name ),*
        ];

        macro_rules! for_each_resonance_case {
            ($macro:ident) => {
                $(
                    $macro!($name);
                )*
            };
        }
    };
}

resonance_cases! {
    // Amino acid zwitterions
    alanine_zwitterion => {
        title: "Alanine Zwitterion",
        molecule: molecules::build_alanine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    asparagine_zwitterion => {
        title: "Asparagine Zwitterion",
        molecule: molecules::build_asparagine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8], bonds: || vec![6, 7] },
        ],
    },
    aspartate_zwitterion => {
        title: "Aspartate Zwitterion",
        molecule: molecules::build_aspartate_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8], bonds: || vec![6, 7] },
        ],
    },
    cysteine_zwitterion => {
        title: "Cysteine Zwitterion",
        molecule: molecules::build_cysteine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    glutamate_zwitterion => {
        title: "Glutamate Zwitterion",
        molecule: molecules::build_glutamate_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![7, 8, 9], bonds: || vec![7, 8] },
        ],
    },
    glutamine_zwitterion => {
        title: "Glutamine Zwitterion",
        molecule: molecules::build_glutamine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![7, 8, 9], bonds: || vec![7, 8] },
        ],
    },
    glycine_zwitterion => {
        title: "Glycine Zwitterion",
        molecule: molecules::build_glycine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    isoleucine_zwitterion => {
        title: "Isoleucine Zwitterion",
        molecule: molecules::build_isoleucine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    leucine_zwitterion => {
        title: "Leucine Zwitterion",
        molecule: molecules::build_leucine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    lysine_zwitterion => {
        title: "Lysine Zwitterion",
        molecule: molecules::build_lysine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    methionine_zwitterion => {
        title: "Methionine Zwitterion",
        molecule: molecules::build_methionine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    proline_zwitterion => {
        title: "Proline Zwitterion",
        molecule: molecules::build_proline_zwitterion,
        systems: [
            { atoms: || vec![5, 6, 7], bonds: || vec![6, 7] },
        ],
    },
    serine_zwitterion => {
        title: "Serine Zwitterion",
        molecule: molecules::build_serine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    threonine_zwitterion => {
        title: "Threonine Zwitterion",
        molecule: molecules::build_threonine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    valine_zwitterion => {
        title: "Valine Zwitterion",
        molecule: molecules::build_valine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },

    // Aromatic amino acid zwitterions
    histidine_zwitterion_aromatic => {
        title: "Histidine Zwitterion (Aromatic Representation)",
        molecule: molecules::build_histidine_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10], bonds: || vec![6, 7, 8, 9, 10] },
        ],
    },
    histidine_zwitterion_kekule => {
        title: "Histidine Zwitterion (Kekulé Representation)",
        molecule: molecules::build_histidine_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10], bonds: || vec![6, 7, 8, 9, 10] },
        ],
    },
    phenylalanine_zwitterion_aromatic => {
        title: "Phenylalanine Zwitterion (Aromatic Representation)",
        molecule: molecules::build_phenylalanine_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11], bonds: || vec![6, 7, 8, 9, 10, 11] },
        ],
    },
    phenylalanine_zwitterion_kekule => {
        title: "Phenylalanine Zwitterion (Kekulé Representation)",
        molecule: molecules::build_phenylalanine_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11], bonds: || vec![6, 7, 8, 9, 10, 11] },
        ],
    },
    tryptophan_zwitterion_aromatic => {
        title: "Tryptophan Zwitterion (Aromatic Representation)",
        molecule: molecules::build_tryptophan_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14], bonds: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14, 15] },
        ],
    },
    tryptophan_zwitterion_kekule => {
        title: "Tryptophan Zwitterion (Kekulé Representation)",
        molecule: molecules::build_tryptophan_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14], bonds: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14, 15] },
        ],
    },
    tyrosine_zwitterion_aromatic => {
        title: "Tyrosine Zwitterion (Aromatic Representation)",
        molecule: molecules::build_tyrosine_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12], bonds: || vec![6, 7, 8, 9, 10, 11, 12] },
        ],
    },
    tyrosine_zwitterion_kekule => {
        title: "Tyrosine Zwitterion (Kekulé Representation)",
        molecule: molecules::build_tyrosine_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12], bonds: || vec![6, 7, 8, 9, 10, 11, 12] },
        ],
    },

    // Nucleobases
    adenine_aromatic => {
        title: "Adenine (Aromatic Representation)",
        molecule: molecules::build_adenine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    adenine_kekule => {
        title: "Adenine (Kekulé Representation)",
        molecule: molecules::build_adenine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    cytosine_aromatic => {
        title: "Cytosine (Aromatic Representation)",
        molecule: molecules::build_cytosine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    cytosine_kekule => {
        title: "Cytosine (Kekulé Representation)",
        molecule: molecules::build_cytosine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    guanine_aromatic => {
        title: "Guanine (Aromatic Representation)",
        molecule: molecules::build_guanine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] },
        ],
    },
    guanine_kekule => {
        title: "Guanine (Kekulé Representation)",
        molecule: molecules::build_guanine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] },
        ],
    },
    thymine_aromatic => {
        title: "Thymine (Aromatic Representation)",
        molecule: molecules::build_thymine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    thymine_kekule => {
        title: "Thymine (Kekulé Representation)",
        molecule: molecules::build_thymine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    uracil_aromatic => {
        title: "Uracil (Aromatic Representation)",
        molecule: molecules::build_uracil_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    uracil_kekule => {
        title: "Uracil (Kekulé Representation)",
        molecule: molecules::build_uracil_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },

    // Nucleosides and phosphate fragments
    adenosine_aromatic => {
        title: "Adenosine (Aromatic Representation)",
        molecule: molecules::build_adenosine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    adenosine_kekule => {
        title: "Adenosine (Kekulé Representation)",
        molecule: molecules::build_adenosine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    deoxyadenosine_aromatic => {
        title: "Deoxyadenosine (Aromatic Representation)",
        molecule: molecules::build_deoxyadenosine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    deoxyadenosine_kekule => {
        title: "Deoxyadenosine (Kekulé Representation)",
        molecule: molecules::build_deoxyadenosine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    dimethyl_phosphate_anion => {
        title: "Dimethyl Phosphate Anion",
        molecule: molecules::build_dimethyl_phosphate_anion,
        systems: [
            { atoms: || vec![0, 1, 2], bonds: || vec![0, 1] },
        ],
    },
    dinucleotide_backbone => {
        title: "Dinucleotide Backbone Fragment",
        molecule: molecules::build_dinucleotide_backbone,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![4, 5] },
        ],
    },

    // Aromatic heterocycles and polycycles
    abaxes_dithiacyclohexene => {
        title: "ABAXES - 1,4-dithiacyclohex-2-ene",
        molecule: molecules::build_abaxes,
        systems: [
            { atoms: || vec![0, 1, 2, 3], bonds: || vec![0, 1, 2] },
        ],
    },
    acridine_aromatic => {
        title: "Acridine (Aromatic Representation)",
        molecule: molecules::build_acridine_aromatic,
        systems: [
            { atoms: || (0..14).collect(), bonds: || (0..16).collect() },
        ],
    },
    acridine_kekule => {
        title: "Acridine (Kekulé Representation)",
        molecule: molecules::build_acridine_kekule,
        systems: [
            { atoms: || (0..14).collect(), bonds: || (0..16).collect() },
        ],
    },
    methyl_pyridine_aromatic => {
        title: "4-Methylpyridine (Aromatic Representation)",
        molecule: molecules::build_4_methylpyridine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5], bonds: || vec![0, 1, 2, 3, 4, 5] },
        ],
    },
    methyl_pyridine_kekule => {
        title: "4-Methylpyridine (Kekulé Representation)",
        molecule: molecules::build_4_methylpyridine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5], bonds: || vec![0, 1, 2, 3, 4, 5] },
        ],
    },
    trinitrobenzene_aromatic => {
        title: "1,3,5-Trinitrobenzene (Aromatic Representation)",
        molecule: molecules::build_trinitrobenzene_aromatic,
        systems: [
            { atoms: || (0..15).collect(), bonds: || (0..15).collect() },
        ],
    },
    trinitrobenzene_kekule => {
        title: "1,3,5-Trinitrobenzene (Kekulé Representation)",
        molecule: molecules::build_trinitrobenzene_kekule,
        systems: [
            { atoms: || (0..15).collect(), bonds: || (0..15).collect() },
        ],
    },

    // Carbonyl and sulfur resonance benchmarks
    acbuol_isolated_lactone => {
        title: "ACBUOL - Isolated Lactone Resonance",
        molecule: molecules::build_acbuol,
        systems: [
            { atoms: || vec![4, 5, 6], bonds: || vec![8, 9] },
        ],
    },
    chloroacetyl_chloride => {
        title: "Chloroacetyl Chloride",
        molecule: molecules::build_chloroacetyl_chloride,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    methanesulfonamide_cross_conjugation => {
        title: "Methanesulfonamide (Cross-Conjugation)",
        molecule: molecules::build_methanesulfonamide,
        systems: [
            { atoms: || vec![0, 1, 2], bonds: || vec![0, 1] },
        ],
    },
    tetramethylthiourea_resonance => {
        title: "Tetramethylthiourea Resonance",
        molecule: molecules::build_tetramethylthiourea,
        systems: [
            { atoms: || vec![0, 1, 2, 3], bonds: || vec![0, 1, 2] },
        ],
    },

    // Negative controls (minimal resonance)
    adamantane_saturated_polycycle => {
        title: "Adamantane (Saturated Polycyclic Alkane)",
        molecule: molecules::build_adamantane,
        systems: [],
    },
    afurpo10_cage_ether => {
        title: "AFURPO10 - Saturated Cage Ether",
        molecule: molecules::build_afurpo10,
        systems: [],
    },
    choline_cation_saturated_salt => {
        title: "Choline Cation (Saturated Quaternary Ammonium Salt)",
        molecule: molecules::build_choline_cation,
        systems: [],
    },
    decalin_saturated_fused_bicycle => {
        title: "Decalin (Saturated Fused Bicyclic Alkane)",
        molecule: molecules::build_decalin,
        systems: [],
    },
    dimethyl_sulfoxide_polar_bond => {
        title: "Dimethyl Sulfoxide (Polar S-O bond)",
        molecule: molecules::build_dimethyl_sulfoxide,
        systems: [],
    },
    phosphinane_saturated_heterocycle => {
        title: "Phosphinane (Saturated Heterocycle)",
        molecule: molecules::build_phosphinane,
        systems: [],
    },

    // Inorganic hypervalent systems
    perchlorate_anion_delocalization => {
        title: "Perchlorate Anion Delocalization",
        molecule: molecules::build_perchlorate_anion,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4], bonds: || vec![0, 1, 2, 3] },
        ],
    },
    perchloric_acid_hypervalent_center => {
        title: "Perchloric Acid",
        molecule: molecules::build_perchloric_acid,
        systems: [
            { atoms: || vec![0, 1, 2, 3], bonds: || vec![0, 1, 2] },
        ],
    },

    // Extended conjugation benchmarks
    fullerene_c60_conjugated_sphere => {
        title: "Fullerene C60 (Kekulé Representation)",
        molecule: molecules::build_fullerene_c60,
        systems: [
            { atoms: || (0..60).collect(), bonds: || (0..90).collect() },
        ],
    },
}
