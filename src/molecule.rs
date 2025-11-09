use crate::core::atom::{AtomId, Element};
use crate::core::bond::{BondId, BondOrder};
use crate::graph::traits::{AtomView, BondView};

#[derive(Clone, Debug)]
pub struct Atom {
    id: AtomId,
    element: Element,
    formal_charge: i8,
}

impl AtomView for Atom {
    fn id(&self) -> AtomId {
        self.id
    }
    fn element(&self) -> Element {
        self.element
    }
    fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
}

#[derive(Clone, Debug)]
pub struct Bond {
    id: BondId,
    order: BondOrder,
    start: AtomId,
    end: AtomId,
}

impl BondView for Bond {
    fn id(&self) -> BondId {
        self.id
    }
    fn order(&self) -> BondOrder {
        self.order
    }
    fn start_atom_id(&self) -> AtomId {
        self.start
    }
    fn end_atom_id(&self) -> AtomId {
        self.end
    }
}
