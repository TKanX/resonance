#[macro_use]
#[path = "common/cases.rs"]
mod cases;

use cases::ResonanceCase;
use pauling::{ResonanceSystem, find_resonance_systems};
use std::cmp::Ordering;
use std::collections::HashSet;

fn system_cmp(a: &ResonanceSystem, b: &ResonanceSystem) -> Ordering {
    a.atoms.cmp(&b.atoms).then_with(|| a.bonds.cmp(&b.bonds))
}

fn run_resonance_case(case: &ResonanceCase) {
    let molecule = (case.build)();
    let mut actual: Vec<_> = find_resonance_systems(&molecule)
        .expect("perception should succeed")
        .into_iter()
        .map(|system| ResonanceSystem::new(system.atoms, system.bonds))
        .collect();
    actual.sort_by(system_cmp);

    let mut expected = case.expected_systems();
    expected.sort_by(system_cmp);

    assert_eq!(
        actual.len(),
        expected.len(),
        "case {} ({}): expected {} systems but found {}",
        case.slug,
        case.title,
        expected.len(),
        actual.len()
    );

    for (idx, (actual_system, expected_system)) in actual.iter().zip(expected.iter()).enumerate() {
        assert_eq!(
            actual_system, expected_system,
            "case {} ({}), system {} mismatch",
            case.slug, case.title, idx
        );
    }
}

macro_rules! declare_case_test {
    ($case:ident) => {
        #[test]
        fn $case() {
            run_resonance_case(&cases::$case);
        }
    };
}

for_each_resonance_case!(declare_case_test);

#[test]
fn catalog_has_unique_slugs() {
    let mut slugs = HashSet::new();
    for case in cases::ALL_CASES {
        assert!(
            slugs.insert(case.slug),
            "duplicate slug discovered: {}",
            case.slug
        );
        assert!(
            !case.title.is_empty(),
            "case {} should expose a non-empty title",
            case.slug
        );
    }
}
