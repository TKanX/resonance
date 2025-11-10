use crate::core::atom::Element;
use crate::core::bond::BondOrder;
use crate::perception::ChemicalPerception;
use std::collections::{HashMap, HashSet};

pub fn perceive(perception: &mut ChemicalPerception) {
    apply_explicit_aromaticity(perception);
    apply_topological_aromaticity(perception);
}

/// Phase 1: Handles bonds explicitly marked as `BondOrder::Aromatic`.
/// This step "trusts the input" and seeds the aromaticity perception.
fn apply_explicit_aromaticity(perception: &mut ChemicalPerception) {
    let aromatic_bonds: Vec<usize> = perception
        .bonds
        .iter()
        .enumerate()
        .filter(|(_, bond)| bond.order == BondOrder::Aromatic)
        .map(|(idx, _)| idx)
        .collect();

    for bond_idx in aromatic_bonds {
        let bond = &mut perception.bonds[bond_idx];
        if bond.is_aromatic {
            continue;
        } // Already processed
        bond.is_aromatic = true;

        if let Some(start_idx) = perception.atom_id_to_index.get(&bond.start_atom_id) {
            perception.atoms[*start_idx].is_aromatic = true;
        }
        if let Some(end_idx) = perception.atom_id_to_index.get(&bond.end_atom_id) {
            perception.atoms[*end_idx].is_aromatic = true;
        }
    }
}
