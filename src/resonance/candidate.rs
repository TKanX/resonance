use crate::core::atom::Element;
use crate::perception::ChemicalPerception;
use crate::perception::Hybridization;

pub fn determine(perception: &mut ChemicalPerception) {
    for atom in &mut perception.atoms {
        atom.is_conjugation_candidate =
            matches!(atom.hybridization, Hybridization::SP | Hybridization::SP2);
    }

    for atom_idx in 0..perception.atoms.len() {
        if perception.atoms[atom_idx].is_conjugation_candidate {
            continue;
        }

        if is_lone_pair_candidate(perception, atom_idx)
            || is_charged_carbon_candidate(perception, atom_idx)
        {
            perception.atoms[atom_idx].is_conjugation_candidate = true;
        }
    }
}

fn is_lone_pair_candidate(perception: &ChemicalPerception, atom_idx: usize) -> bool {
    let atom = &perception.atoms[atom_idx];

    if atom.lone_pairs > 0 {
        for &(neighbor_idx, _) in &perception.adjacency[atom_idx] {
            if perception.atoms[neighbor_idx].is_conjugation_candidate {
                return true;
            }
        }
    }

    false
}

fn is_charged_carbon_candidate(perception: &ChemicalPerception, atom_idx: usize) -> bool {
    let atom = &perception.atoms[atom_idx];

    if atom.element != Element::C {
        return false;
    }

    match atom.formal_charge {
        1 if atom.total_degree == 3 => true,
        -1 => true,
        _ => false,
    }
}
