use super::graph::TestMolecule;
use pauling::{AtomId, BondId, BondOrder, Element};

pub fn build_glycine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::H, 0);
    mol.add_atom(6, Element::H, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 0, 5, BondOrder::Single);
    mol.add_bond(5, 0, 6, BondOrder::Single);
    mol.add_bond(6, 0, 7, BondOrder::Single);
    mol.add_bond(7, 1, 8, BondOrder::Single);
    mol.add_bond(8, 1, 9, BondOrder::Single);

    mol
}
