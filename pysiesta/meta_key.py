"""Metadata fields and display labels used by siesta_getmeta.py.

Edit this file to add, remove, or reorder fields in the generated output.
Each entry is (JSON key, table label).
"""

META_KEYS = {
    "main": [
        ("molecular_formula", "Molecular formula"),
        ("molecule_name", "Molecule name"),
        ("data_id", "Data ID"),
        ("atom_count", "Atom count"),
        ("charge_e", "Charge (e)"),
        ("spin_multiplicity", "Spin multiplicity"),
        ("homo_lumo_gap_ev", "HOMO–LUMO gap (eV)"),
        ("exchange_correlation", "Ex-Corr"),
        ("basis_size", "Basis set size"),
        ("pseudopotential_species", "Pseudopotential species"),
        ("program", "Package"),
        ("database_source", "DB source"),
    ],
    "structure": {
        "group": "1_molecular_structure",
        "title": "Molecular structure",
        "keys": [("structure_3d", "3D structure")],
    },
    "calculation": {
        "group": "2_calculation_conditions",
        "title": "Calculation conditions",
        "keys": [
            ("program_version", "Program version"),
            ("dispersion_correction", "Dispersion correction"),
            ("spin_mode", "Spin mode"),
            ("relativistic_effects", "Relativistic effects"),
            ("basis_details", "Basis details"),
        ],
    },
    "results": {
        "group": "3_electronic_structure_results",
        "title": "Electronic structure results",
        "keys": [
            ("total_energy_ev", "Total energy (eV)"),
            ("homo_energy_ev", "HOMO energy (eV)"),
            ("lumo_energy_ev", "LUMO energy (eV)"),
            ("charge_density", "Charge density"),
            ("dos_pdos", "DOS / PDOS"),
            ("molecular_orbitals", "Molecular orbitals"),
            ("dipole_moment_debye", "Dipole moment (Debye)"),
            ("magnetic_moment", "Magnetic moment"),
        ],
    },
    "data_source": {
        "group": "4_provenance",
        "title": "Provenance",
        "keys": [
            ("original_data_id", "Original data ID"),
            ("references", "References"),
            ("dates", "Dates"),
            ("provider", "Provider"),
            ("original_calculation_files", "Original calculation files"),
        ],
    },
}

META_LABELS = {key: label for key, label in META_KEYS["main"]}
for section in ("structure", "calculation", "results", "data_source"):
    META_LABELS.update(dict(META_KEYS[section]["keys"]))

GROUP_LABELS = {
    detail["group"]: detail["title"]
    for section in ("structure", "calculation", "results", "data_source")
    for detail in (META_KEYS[section],)
}
