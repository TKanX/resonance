//! Resonance detection logic and associated data structures.
//!
//! The perception pipeline delegates the final grouping step to this module
//! once atoms and bonds have been annotated with conjugation metadata.

pub mod candidate;
mod find;
mod system;

/// Identifies conjugated components and constructs [`ResonanceSystem`] values.
pub use find::find_systems;
/// Canonical representation of a resonance system.
pub use system::ResonanceSystem;
