use crate::core::atom::AtomId;
use crate::core::bond::BondOrder;
use crate::errors::PerceptionError;
use crate::perception::ChemicalPerception;
use std::collections::{HashMap, VecDeque};

const KEKULIZATION_ATTEMPT_LIMIT: usize = 1000;

pub fn kekulize(perception: &mut ChemicalPerception) -> Result<(), PerceptionError> {
    let mut visited_bonds = vec![false; perception.bonds.len()];
    let mut total_attempts = 0;

    for bond_idx in 0..perception.bonds.len() {
        if perception.bonds[bond_idx].is_aromatic && !visited_bonds[bond_idx] {
            let component_bond_indices =
                collect_aromatic_component(perception, bond_idx, &mut visited_bonds);

            assign_kekule_orders(perception, &component_bond_indices, &mut total_attempts)?;
        }
    }
    Ok(())
}

fn collect_aromatic_component(
    perception: &ChemicalPerception,
    start_bond_idx: usize,
    visited_bonds: &mut [bool],
) -> Vec<usize> {
    let mut component_indices = Vec::new();
    let mut queue = VecDeque::new();

    visited_bonds[start_bond_idx] = true;
    queue.push_back(start_bond_idx);

    while let Some(bond_idx) = queue.pop_front() {
        component_indices.push(bond_idx);
        let bond = &perception.bonds[bond_idx];

        for atom_id in [bond.start_atom_id, bond.end_atom_id] {
            if let Some(&atom_idx) = perception.atom_id_to_index.get(&atom_id) {
                for &(_, neighbor_bond_id) in &perception.adjacency[atom_idx] {
                    if let Some(&neighbor_bond_idx) =
                        perception.bond_id_to_index.get(&neighbor_bond_id)
                    {
                        if perception.bonds[neighbor_bond_idx].is_aromatic
                            && !visited_bonds[neighbor_bond_idx]
                        {
                            visited_bonds[neighbor_bond_idx] = true;
                            queue.push_back(neighbor_bond_idx);
                        }
                    }
                }
            }
        }
    }
    component_indices
}

fn assign_kekule_orders(
    perception: &mut ChemicalPerception,
    component_bond_indices: &[usize],
    total_attempts: &mut usize,
) -> Result<(), PerceptionError> {
    let unassigned_bond_indices: Vec<usize> = component_bond_indices
        .iter()
        .copied()
        .filter(|&idx| perception.bonds[idx].kekule_order.is_none())
        .collect();

    if unassigned_bond_indices.is_empty() {
        return Ok(());
    }

    let mut atom_double_bond_counts = HashMap::new();

    if !kekule_backtrack(
        perception,
        0,
        &unassigned_bond_indices,
        &mut atom_double_bond_counts,
        total_attempts,
    ) {
        return Err(PerceptionError::KekulizationFailed(*total_attempts));
    }

    for &bond_idx in component_bond_indices {
        if perception.bonds[bond_idx].kekule_order.is_none() {
            perception.bonds[bond_idx].kekule_order = Some(BondOrder::Single);
        }
    }

    Ok(())
}

fn kekule_backtrack(
    perception: &mut ChemicalPerception,
    position: usize,
    unassigned_bonds: &[usize],
    atom_counts: &mut HashMap<AtomId, u8>,
    attempts: &mut usize,
) -> bool {
    if *attempts >= KEKULIZATION_ATTEMPT_LIMIT {
        return false;
    }

    if position == unassigned_bonds.len() {
        return true;
    }

    *attempts += 1;
    let bond_idx = unassigned_bonds[position];
    let (start_id, end_id) = {
        let bond = &perception.bonds[bond_idx];
        (bond.start_atom_id, bond.end_atom_id)
    };

    let can_assign_double = atom_counts.get(&start_id).copied().unwrap_or(0) == 0
        && atom_counts.get(&end_id).copied().unwrap_or(0) == 0;

    if can_assign_double {
        perception.bonds[bond_idx].kekule_order = Some(BondOrder::Double);
        *atom_counts.entry(start_id).or_insert(0) += 1;
        *atom_counts.entry(end_id).or_insert(0) += 1;

        if kekule_backtrack(
            perception,
            position + 1,
            unassigned_bonds,
            atom_counts,
            attempts,
        ) {
            return true;
        }

        *atom_counts.get_mut(&start_id).unwrap() -= 1;
        *atom_counts.get_mut(&end_id).unwrap() -= 1;
        perception.bonds[bond_idx].kekule_order = None;
    }

    perception.bonds[bond_idx].kekule_order = Some(BondOrder::Single);

    if kekule_backtrack(
        perception,
        position + 1,
        unassigned_bonds,
        atom_counts,
        attempts,
    ) {
        return true;
    }

    perception.bonds[bond_idx].kekule_order = None;

    false
}
