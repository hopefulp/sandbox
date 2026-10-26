#!/usr/bin/env python3
"""Extract molecule metadata from a SIESTA calculation or dataset directory.

Print tables to stdout and save metadata.json / molecule_names.tsv in the current
working directory; use -o DIR to choose a different output directory.
Single and group modes use the same records array, one record per calculation.
Missing/unverified values are null. Large scientific files are referenced, not read.
HOMO/LUMO extraction is limited to integer, closed-shell, single-k-point runs.
"""
import argparse
import csv
from collections import Counter
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import subprocess
import textwrap

from meta_key import GROUP_LABELS, META_KEYS, META_LABELS

NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
ELEMENTS = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()


def number(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


def key(s):
    return re.sub(r'[_.-]', '', s).lower()


def read_fdf(path, stack=()):
    """Read scalar/block FDF values, includes and redirected blocks; first wins."""
    path = path.resolve()
    if path in stack:
        raise ValueError('Cyclic FDF include: ' + str(path))
    values, blocks = {}, {}
    lines = path.read_text(errors='replace').splitlines()
    i = 0
    while i < len(lines):
        line = re.split(r'[#!]', lines[i], maxsplit=1)[0].strip()
        i += 1
        if not line:
            continue
        parts = line.split()
        tag = parts[0].lower()
        if tag == '%include':
            v, b = read_fdf(path.parent / parts[1], (*stack, path))
            for k, x in v.items():
                values.setdefault(k, x)
            for k, x in b.items():
                blocks.setdefault(k, x)
        elif tag == '%block':
            rows = []
            if '<' in parts:
                rows = (path.parent / parts[parts.index('<') + 1]).read_text().splitlines()
            else:
                while i < len(lines):
                    row = re.split(r'[#!]', lines[i], maxsplit=1)[0].strip()
                    i += 1
                    if row.lower().startswith('%endblock'):
                        break
                    if row:
                        rows.append(row)
            blocks.setdefault(key(parts[1]), rows)
        elif len(parts) > 1 and not tag.startswith('%'):
            if parts[1] == '<':
                v, _ = read_fdf(path.parent / parts[2], (*stack, path))
                if key(parts[0]) in v:
                    values.setdefault(key(parts[0]), v[key(parts[0])])
            else:
                values.setdefault(key(parts[0]), ' '.join(parts[1:]))
    return values, blocks


def structure(folder, label, values, blocks):
    p = folder / (label + '.STRUCT_OUT')
    if p.exists():
        rows = p.read_text().splitlines()
        cell = [[number(x) for x in row.split()] for row in rows[:3]]
        n = int(rows[3])
        atoms = []
        for row in rows[4:4+n]:
            f = row.split()
            z = int(f[1])
            xyz = [sum(number(f[j+2]) * cell[j][k] for j in range(3)) for k in range(3)]
            atoms.append({'element': ELEMENTS[z], 'position_angstrom': xyz})
        if len(atoms) != n:
            raise ValueError('Truncated STRUCT_OUT')
        return atoms, cell, str(p)
    species = {int(r.split()[0]): ELEMENTS[int(r.split()[1])] for r in blocks.get('chemicalspecieslabel', [])}
    fmt = key(values.get('atomiccoordinatesformat', ''))
    if fmt not in ('ang', 'notscaledcartesianang', 'bohr', 'notscaledcartesianbohr'):
        raise ValueError('No STRUCT_OUT; unsupported/missing input coordinate format: ' + fmt)
    scale = 0.529177210903 if 'bohr' in fmt else 1.0
    atoms = [{'element': species[int(f[3])], 'position_angstrom': [number(x)*scale for x in f[:3]]}
             for f in (r.split() for r in blocks.get('atomiccoordinatesandatomicspecies', []))]
    if not atoms:
        raise ValueError('No atomic structure found')
    return atoms, None, 'input FDF AtomicCoordinatesAndAtomicSpecies'


def load_names(source):
    """Read TABS molecule names and coordinates from PDF or text."""
    path = Path(source).expanduser().resolve()
    if path.suffix.lower() == '.pdf':
        if not shutil.which('pdftotext'):
            raise ValueError('Reading PDF names requires pdftotext (Poppler)')
        result = subprocess.run(['pdftotext', '-layout', str(path), '-'],
                                capture_output=True, text=True)
        if result.returncode:
            raise ValueError(f'Cannot read name PDF {path}: {result.stderr.strip()}')
        text = result.stdout
    else:
        text = path.read_text(encoding='utf-8')
    names = []
    ids = set()
    for entry in re.split(r'Molecule#\s*', text)[1:]:
        match = re.match(r'(\d+)\s*\n(.*?)PointGroup:', entry, re.S)
        if not match:
            raise ValueError(f'Malformed TABS entry in {path}')
        ident, name = int(match[1]), ' '.join(match[2].split())
        if not name or ident in ids:
            raise ValueError(f'Empty or duplicate molecule name for ID {ident} in {path}')
        ids.add(ident)
        atoms = [{'element': m[1], 'position_angstrom': [number(m[i]) for i in (2, 3, 4)]}
                 for m in re.finditer(r'^\s*([A-Z][a-z]?)\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + r')\s*$', entry, re.M)]
        count = re.search(r'^\s*(\d+)\s*$', entry[match.end():], re.M)
        if not atoms or count is None or len(atoms) != int(count[1]):
            raise ValueError(f'Incomplete coordinates for molecule {ident} in {path}')
        names.append({'id': ident, 'name': name, 'signature': structure_signature(atoms), 'graph': bond_graph(atoms)})
    if not names:
        raise ValueError(f'No TABS Molecule# / PointGroup entries in {path}')
    return names


def structure_signature(atoms):
    """Element-pair distances are invariant to translation, rotation and ordering."""
    distances = {}
    for i, a in enumerate(atoms):
        for b in atoms[i+1:]:
            pair = tuple(sorted((a['element'], b['element'])))
            distances.setdefault(pair, []).append(math.dist(a['position_angstrom'], b['position_angstrom']))
    return (Counter(a['element'] for a in atoms),
            {pair: sorted(values) for pair, values in distances.items()})


def bond_graph(atoms):
    """Approximate connectivity for the eight elements present in TABS."""
    radii = dict(H=0.31, C=0.76, N=0.71, O=0.66, F=0.57, S=1.05, Cl=1.02, Br=1.20)
    if any(a['element'] not in radii for a in atoms):
        return None
    neighbors = [set() for _ in atoms]
    for i, a in enumerate(atoms):
        for j in range(i):
            b = atoms[j]
            distance = math.dist(a['position_angstrom'], b['position_angstrom'])
            if 0.4 < distance < 1.2 * (radii[a['element']] + radii[b['element']]):
                neighbors[i].add(j)
                neighbors[j].add(i)
    labels = [(a['element'], tuple(sorted(atoms[j]['element'] for j in neighbors[i])))
              for i, a in enumerate(atoms)]
    return labels, neighbors


def same_connectivity(left, right):
    if left is None or right is None:
        return False
    labels, edges = left
    other_labels, other_edges = right
    if Counter(labels) != Counter(other_labels):
        return False
    candidates = [[j for j, label in enumerate(other_labels) if label == value] for value in labels]
    order = sorted(range(len(labels)), key=lambda i: (len(candidates[i]), -len(edges[i])))
    mapping, used = {}, set()
    def search(depth):
        if depth == len(order):
            return True
        i = order[depth]
        for j in candidates[i]:
            if j in used or any((k in edges[i]) != (v in other_edges[j]) for k, v in mapping.items()):
                continue
            mapping[i] = j
            used.add(j)
            if search(depth + 1):
                return True
            used.remove(j)
            del mapping[i]
        return False
    return search(0)


def match_name(atoms, names):
    counts, distances = structure_signature(atoms)
    matches = []
    for entry in names:
        other_counts, other_distances = entry['signature']
        if counts != other_counts:
            continue
        if all(abs(a-b) <= 0.005 for pair, values in distances.items()
               for a, b in zip(values, other_distances[pair])):
            matches.append(entry)
    if matches:
        return matches, 'interatomic distances (0.005 angstrom tolerance)'
    graph = bond_graph(atoms)
    matches = [entry for entry in names if counts == entry['signature'][0]
               and same_connectivity(graph, entry['graph'])]
    return matches, 'unique inferred bond connectivity; geometry may differ'


def describe_name_match(matches, method):
    # Certainty is the matching rule's result, not a stereochemical validation.
    certain = len(matches) == 1
    display_name = ' / '.join(entry['name'] for entry in matches) or 'Unknown'
    if not certain:
        display_name += ' (?)'
    return {'display_name': display_name,
            'candidates': [{'id': entry['id'], 'name': entry['name']} for entry in matches],
            'method': method, 'certainty': '확실' if certain else '불확실'}


def extract(args):
    folder = args.folder.expanduser().resolve()
    if not folder.is_dir():
        raise ValueError('Not a directory: ' + str(folder))
    inp = folder / args.input
    values, blocks = read_fdf(inp)
    logs = sorted(folder.glob('fdf-*.log'))
    if args.fdf_log:
        logs = [folder / args.fdf_log]
    if len(logs) > 1:
        raise ValueError('Multiple fdf logs; select one with --fdf-log')
    if logs:
        effective, _ = read_fdf(logs[0])
        values.update(effective)
    get = lambda name: values.get(key(name))
    out = folder / args.stdout
    text = out.read_text(errors='replace') if out.exists() else ''
    def last(pattern, numeric=False):
        hits = re.findall(pattern, text, re.M | re.I)
        return (number(hits[-1]) if numeric else hits[-1].strip()) if hits else None
    def scalar(name):
        v = get(name)
        return number(v.split()[0]) if v is not None else None
    def files(*patterns):
        return sorted({str(p) for pat in patterns for p in folder.glob(pat) if p.is_file()})
    label = get('SystemLabel') or 'siesta'
    atoms, cell, structure_source = structure(folder, label, values, blocks)
    counts = Counter(a['element'] for a in atoms)
    order = (['C'] + (['H'] if 'H' in counts else []) + sorted(set(counts)-{'C','H'})) if 'C' in counts else sorted(counts)
    formula = ''.join(s + (str(counts[s]) if counts[s] != 1 else '') for s in order)
    molecule_name = args.name
    matched_id = None
    name_match = {'candidates': [], 'method': 'user supplied' if args.name else 'unavailable',
                  'certainty': '불확실'}
    names = getattr(args, 'name_map', None)
    if names is not None:
        matches, match_method = match_name(atoms, names)
        name_match = describe_name_match(matches, match_method)
        molecule_name = name_match['display_name']
        matched_id = matches[0]['id'] if len(matches) == 1 else None
    spin = get('Spin')
    if spin is None and get('SpinPolarized') is not None:
        spin = 'polarized' if get('SpinPolarized').lower() in ('t', 'true', '.true.') else 'none'
    electrons = last(r'Total number of electrons:\s*(' + NUM + ')', True)
    homo = lumo = gap = fermi = None
    eig = folder / (label + '.EIG')
    if eig.exists():
        tokens = eig.read_text().split()
        fermi = number(tokens[0])
        nb, ns, nk = map(int, tokens[1:4])
        if len(tokens) != 4 + nk*(1 + nb*ns):
            raise ValueError('Unexpected or truncated EIG contents')
        if ns == 1 and nk == 1 and spin in ('none', 'unpolarized') and electrons is not None and abs(electrons/2-round(electrons/2)) < 1e-8:
            occupied = round(electrons/2)
            es = [number(x) for x in tokens[5:]]
            if 0 < occupied < nb and es == sorted(es) and es[occupied-1] < fermi < es[occupied]:
                homo, lumo = es[occupied-1:occupied+1]
                gap = lumo-homo
    pseudos = []
    used = {r.split()[2] for r in blocks.get('chemicalspecieslabel', [])}
    for species in sorted(used):
        for suffix in ('.psf', '.psml', '.vps'):
            p = folder / (species + suffix)
            if p.exists():
                item = {'species': species, 'header': None}
                if suffix == '.psf':
                    with p.open(errors='replace') as handle:
                        item['header'] = [handle.readline().strip() for _ in range(3)]
                pseudos.append(item)
    dip = last(r'Electric dipole \(Debye\)\s*=\s*([^\n]+)')
    dip = [number(v) for v in dip.split()] if dip else None
    main = {
        'molecular_formula': formula, 'molecule_name': molecule_name, 'data_id': folder.name,
        'atom_count': len(atoms), 'charge_e': scalar('NetCharge'), 'spin_multiplicity': args.multiplicity,
        'homo_lumo_gap_ev': gap, 'exchange_correlation': {'functional': get('XC.functional'), 'authors': get('XC.authors')},
        'basis_size': get('PAO.BasisSize'),
        'pseudopotential_species': sorted(used),
        'program': 'SIESTA', 'database_source': str(folder)}
    groups = {
        '1_molecular_structure': {
            'structure_3d': {'source': structure_source, 'preview_files': files(label+'.png')},
            'atoms': atoms,
            'structure_files': files('*.fdf', '*.STRUCT_OUT', '*.XV', '*.xyz', '*.xsf', '*.cif', 'POSCAR', 'CONTCAR'),
            'identifiers': {'smiles': None, 'inchi': None, 'cas': None}},
        '2_calculation_conditions': {
            'program_version': last(r'^Siesta Version\s*:\s*(.+)$'),
            'dispersion_correction': {k:v for k,v in values.items() if 'vdw' in k or 'dispersion' in k} or None,
            'spin_mode': spin,
            'relativistic_effects': {'spin_orbit': get('SpinOrbit'), 'pseudopotential_headers': [p['header'] for p in pseudos]},
            'scf_convergence': {k:v for k,v in values.items() if ('tolerance' in k and (k.startswith('scf') or k.startswith('dm'))) or k == 'maxscfiterations'},
            'basis_details': {k:v for k,v in values.items() if k.startswith('pao') or k == 'meshcutoff'},
            'pseudopotential_details': pseudos},
        '3_electronic_structure_results': {
            'total_energy_ev': last(r'^siesta:\s+Total\s*=\s*(' + NUM + ')', True),
            'homo_energy_ev': homo, 'lumo_energy_ev': lumo,
            'charge_density': {'files': files(label+'.RHO', label+'.RHO.*', '*.cube'), 'preview_files': files('rho*.png')},
            'dos_pdos': {'files': files('*.DOS', '*.PDOS', '*.dos', '*.pdos')},
            'molecular_orbitals': {'wavefunction_files': files('*.WFSX', '*.WFS'), 'eigenvalue_files': files('*.EIG'), 'orbital_index_files': files('*.ORB_INDX')},
            'dipole_moment_debye': {'vector': dip, 'magnitude': math.sqrt(sum(x*x for x in dip)) if dip else None},
            'magnetic_moment': {'total_mu_b': None, 'local_mu_b': None}},
        '4_provenance': {
            'original_data_id': {'database': args.database, 'id': folder.name},
            'references': {'title': args.reference_title, 'doi': args.doi, 'bibliography_files': files('*.bib')},
            'dates': {'calculation_start': last(r'>> Start of run:\s*(.+)'),
                      'calculation_end': last(r'>> End of run:\s*(.+)'), 'registration': None},
            'provider': args.provider,
            'original_calculation_files': [p for p in files('*') if Path(p).suffix.lower() not in ('.psf', '.psml', '.vps')]}}
    return compact_metadata({'main_metadata': main, 'groups': groups,
                             'file_paths': sorted(str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file()),
                             '_name_match': name_match})



def compact_metadata(record):
    """One representative scalar per field; the structure image is supplied later."""
    main = record['main_metadata']
    xc = main.get('exchange_correlation')
    if isinstance(xc, dict):
        main['exchange_correlation'] = xc.get('authors') or xc.get('functional')
    species = main.get('pseudopotential_species')
    if isinstance(species, list):
        main['pseudopotential_species'] = ', '.join(species) or None
    groups = record['groups']
    calc = groups['2_calculation_conditions']
    electronic = groups['3_electronic_structure_results']
    provenance = groups['4_provenance']

    source = Path(main['database_source'])

    def file_reference(value):
        """Keep a valid path relative to this calculation's DB source."""
        if not value or value == '없음':
            return '없음'
        path = Path(value)
        if not path.is_absolute():
            path = source / path
        try:
            relative = path.relative_to(source)
        except ValueError:
            return '없음'
        if '..' in relative.parts or not path.is_file():
            return '없음'
        return str(relative)

    def first_file(value, keys):
        if not isinstance(value, dict):
            return file_reference(value)
        for key in keys:
            for candidate in value.get(key) or []:
                reference = file_reference(candidate)
                if reference != '없음':
                    return reference
        return '없음'

    dispersion = calc.get('dispersion_correction')
    if isinstance(dispersion, dict):
        dispersion = next((dispersion[k] for k in ('vdwfunctional', 'dispersioncorrection')
                           if dispersion.get(k)), None)
    relativistic = calc.get('relativistic_effects')
    if isinstance(relativistic, dict):
        soc = str(relativistic.get('spin_orbit') or '').lower()
        relativistic = 'SOC' if soc in ('t', 'true', '.true.', 'yes') else None
    basis = calc.get('basis_details')
    if isinstance(basis, dict):
        basis = basis.get('paobasissize')
    dipole = electronic.get('dipole_moment_debye')
    if isinstance(dipole, dict):
        dipole = dipole.get('magnitude')
    magnetic = electronic.get('magnetic_moment')
    if isinstance(magnetic, dict):
        magnetic = magnetic.get('total_mu_b')
    original_id = provenance.get('original_data_id')
    if isinstance(original_id, dict):
        original_id = ' / '.join(str(original_id[k]) for k in ('database', 'id')
                                 if original_id.get(k) is not None)
    reference = provenance.get('references')
    if isinstance(reference, dict):
        reference = reference.get('doi') or reference.get('title')
    dates = provenance.get('dates')
    if isinstance(dates, dict):
        dates = dates.get('calculation_end') or dates.get('calculation_start')
    originals = provenance.get('original_calculation_files')
    if isinstance(originals, list):
        inputs = [p for p in originals if Path(p).suffix.lower() == '.fdf']
        preferred = next((p for p in inputs if Path(p).name == 'RUN.fdf'), None)
        originals = preferred or (inputs[0] if inputs else None)
    originals = file_reference(originals)
    raw_groups = {
        '1_molecular_structure': {'structure_3d': None},
        '2_calculation_conditions': {
            'program_version': calc.get('program_version'),
            'dispersion_correction': dispersion,
            'spin_mode': calc.get('spin_mode'),
            'relativistic_effects': relativistic,
            'basis_details': basis},
        '3_electronic_structure_results': {
            'total_energy_ev': electronic.get('total_energy_ev'),
            'homo_energy_ev': electronic.get('homo_energy_ev'),
            'lumo_energy_ev': electronic.get('lumo_energy_ev'),
            'charge_density': first_file(electronic.get('charge_density'), ('files',)),
            'dos_pdos': first_file(electronic.get('dos_pdos'), ('files',)),
            'molecular_orbitals': first_file(electronic.get('molecular_orbitals'), ('wavefunction_files',)),
            'dipole_moment_debye': dipole,
            'magnetic_moment': magnetic},
        '4_provenance': {
            'original_data_id': original_id,
            'references': reference,
            'dates': dates,
            'provider': provenance.get('provider'),
            'original_calculation_files': originals}}
    record['main_metadata'] = {
        key: main.get(key) for key, _label in META_KEYS['main']}
    record['groups'] = {
        META_KEYS[section]['group']: {
            key: raw_groups[META_KEYS[section]['group']].get(key)
            for key, _label in META_KEYS[section]['keys']}
        for section in ('structure', 'calculation', 'results', 'data_source')}
    return record


# Labels and selected keys are maintained in meta_key.py.
DISPLAY_LABELS = {**META_LABELS, **GROUP_LABELS}


def display_label(name):
    """Map a metadata key to its display label, preserving unknown FDF keys."""
    if name in DISPLAY_LABELS:
        return DISPLAY_LABELS[name]
    if '_' in name:
        return name.replace('_', ' ').capitalize()
    if name.startswith(('pao', 'scf', 'dm', 'vdw')):
        return name
    return name[:1].upper() + name[1:]


def table_rows(value, prefix=''):
    """Flatten nested metadata without dropping list entries or empty values."""
    if isinstance(value, dict) and value:
        for name, child in value.items():
            label = display_label(name)
            if name == 'homo_lumo_gap_ev' and isinstance(child, (int, float)):
                yield f'{prefix} / {label}' if prefix else label, f'{child:.3f}'
                continue
            if name == 'exchange_correlation' and isinstance(child, dict):
                rendered = ' / '.join(str(child.get(k) or 'null') for k in ('functional', 'authors'))
                yield f'{prefix} / {label}' if prefix else label, rendered
                continue
            yield from table_rows(child, f'{prefix} / {label}' if prefix else label)
    elif isinstance(value, list) and value:
        for i, child in enumerate(value, 1):
            yield from table_rows(child, f'{prefix}[{i}]')
    else:
        rendered = '' if value is None else (value if isinstance(value, str) else json.dumps(value, ensure_ascii=False, allow_nan=False))
        yield prefix, rendered


def print_table(title, payload):
    rows = list(table_rows(payload))
    width = max(40, min(120, shutil.get_terminal_size((100, 24)).columns))
    field_width = min(max(len(field) for field, _ in rows), 48, (width - 7) // 2)
    field_width = max(field_width, len('Field'))
    value_width = width - field_width - 7
    border = '+' + '-' * (field_width + 2) + '+' + '-' * (value_width + 2) + '+'

    def row(field, value):
        left = textwrap.wrap(field, field_width) or ['']
        right = [part for line in value.splitlines() for part in (textwrap.wrap(line, value_width) or [''])] or ['']
        for i in range(max(len(left), len(right))):
            a = left[i] if i < len(left) else ''
            b = right[i] if i < len(right) else ''
            print(f'| {a:<{field_width}} | {b:<{value_width}} |')

    print(title)
    print(border)
    row('Field', 'Value')
    print(border)
    for field, value in rows:
        row(field, value)
    print(border)


def calculation_folders(root, stdout, limit=None):
    """Find calculations in natural traversal order, stopping once limit is reached."""
    if not root.is_dir():
        raise ValueError('Not a directory: ' + str(root))
    if Path(stdout).is_absolute() or len(Path(stdout).parts) != 1:
        raise ValueError('--stdout must be a filename within each calculation directory')
    def fail(error):
        raise error
    def natural_name(name):
        return [(1, int(part)) if part.isdigit() else (0, part.lower())
                for part in re.split(r'(\d+)', name)]
    found = []
    for directory, children, filenames in os.walk(root, onerror=fail):
        children.sort(key=lambda name: (natural_name(name), name))
        if stdout in filenames:
            found.append(Path(directory))
            children[:] = []
            if limit is not None and len(found) >= limit:
                break
    def natural_order(path):
        return [(1, int(part)) if part.isdigit() else (0, part.lower())
                for part in re.split(r'(\d+)', str(path.relative_to(root)))]
    return sorted(found, key=natural_order)


def print_dataset_table(records, section):
    """One metadata set per row, with web display labels as column headers."""
    rows = []
    for record in records:
        data = record['main_metadata']
        if section not in ('main', 'all'):
            data = next(v for k, v in record['groups'].items() if k.startswith(section + '_'))
        row = {'Data ID': record['main_metadata']['data_id']}
        row.update(table_rows({k: v for k, v in data.items() if k != 'database_source_directory'}))
        rows.append(row)
    headers = list(dict.fromkeys(key for row in rows for key in row))
    # Keep each record on one line, including strings containing newlines/tabs.
    values = [[' '.join(str(row.get(header, 'null')).split()) for header in headers] for row in rows]
    widths = [max(len(header), *(len(row[i]) for row in values)) for i, header in enumerate(headers)]
    border = '+' + '+'.join('-' * (width + 2) for width in widths) + '+'
    def line(cells):
        print('| ' + ' | '.join(f'{cell:<{width}}' for cell, width in zip(cells, widths)) + ' |')
    print(border)
    line(headers)
    print(border)
    for row in values:
        line(row)
    print(border)
    if section == 'all':
        for number in ('1', '2', '3', '4'):
            print()
            print_dataset_table(records, number)


def write_json(path, payload):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def print_name_table(records, errors, handle=None):
    """Print four TSV columns to stdout or a supplied file handle."""
    rows = []
    for record in records:
        main = record['main_metadata']
        match = record['_name_match']
        name = main['molecule_name'] or 'Unknown'
        if match['certainty'] == '불확실' and not name.endswith(' (?)'):
            name += ' (?)'
        rows.append([main['data_id'], name,
                     ' / '.join(f"{entry['id']:04d}" for entry in match['candidates']) or '-',
                     match['certainty']])
    for error in errors:
        rows.append([Path(error['source_directory']).name, 'Extraction error (?)', '-', '불확실'])
    rows.sort(key=lambda row: [(1, int(part)) if part.isdigit() else (0, part.lower())
                              for part in re.split(r'(\d+)', row[0])])
    writer = csv.writer(sys.stdout if handle is None else handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['Data ID', 'Molecule name', '두번째 파일 id', '확실/불확실'])
    writer.writerows(rows)


def save_name_table(path, records, errors):
    with path.open('w', encoding='utf-8', newline='') as handle:
        print_name_table(records, errors, handle)


def save_records(output_dir, root, records, errors):
    """Use one JSON envelope for both single and group extraction."""
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_json(output_dir / 'metadata.json', {
            'root_directory': str(root), 'count': len(records),
            'records': [{k: record[k] for k in ('main_metadata', 'groups', 'file_paths')} for record in records],
            'errors': errors})
        save_name_table(output_dir / 'molecule_names.tsv', records, errors)
        print(f'Saved: {(output_dir / "metadata.json").resolve()}', file=sys.stderr)
        print(f'Saved: {(output_dir / "molecule_names.tsv").resolve()}', file=sys.stderr)


def extract_dataset(args, root, folders):
    records, errors = [], []
    for folder in folders:
        options = argparse.Namespace(**vars(args))
        options.folder = folder
        try:
            data = extract(options)
            records.append({'source_directory': str(folder), **data})
        except (OSError, ValueError, KeyError, IndexError) as exc:
            errors.append({'source_directory': str(folder), 'error': str(exc)})
            print(f'error: {folder}: {exc}', file=sys.stderr)
    save_records(args.output_dir, root, records, errors)
    if args.section == 'names':
        print_name_table(records, errors)
    elif records:
        print_dataset_table(records, args.section)
    print(f'{len(records)} calculation(s) extracted; {len(errors)} failed.', file=sys.stderr)
    return 1 if errors else 0


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''주요 옵션
-n, --name FILE      PDF/text 구조와 비교해 분자 이름과 PDF ID를 찾습니다.
                     파일 대신 분자 이름을 직접 지정할 수도 있습니다.
-nt, --ntest N       -g에서 계산 폴더 N개를 찾으면 탐색을 멈춥니다.
                     생략하면 전체 계산 폴더를 처리합니다.
-o, --output-dir DIR metadata.json과 molecule_names.tsv의 저장 폴더입니다.
                     기본값은 명령을 실행한 현재 폴더입니다.
-s, --section       main 기본 표, all 전체, 1 구조, 2 계산 조건,
                     3 전자구조, 4 출처, names 이름 비교표를 선택합니다.
-g, --group         하위 계산 폴더를 함께 처리합니다.

실행 예제

1. 하위 계산 Data folder: Metadata main
siesta_getmeta.py DATA_FOLDER -n TABS.pdf
한 계산 폴더의 기본 표를 출력하고 PDF의 이름과 ID를 찾습니다.

1-1. 하위 계산 Data folder: Metadata detail
siesta_getmeta.py DATA_FOLDER -n TABS.pdf -s all
한 데이터의 기본 표와 구조, 계산 조건, 전자구조, 출처 정보를 출력합니다.

2. 계산 Data Group
siesta_getmeta.py DATA_ROOT -g -nt 5 -n TABS.pdf
하위 계산 폴더 5개를 처리하고 폴더별 메타데이터 표를 출력합니다.

3. 출력 파일과 저장 위치

현재 폴더에 저장합니다. -o DIR로 위치를 바꿀 수 있으며, 재실행하면 덮어씁니다.

metadata.json: 기본·상세 메타데이터, 파일 경로, 추출 오류를 저장합니다.
파일 경로는 DB source 기준 상대경로입니다. DB source와 합쳐 실제 파일을 찾습니다.

molecule_names.tsv: -n 파일 지정 시 Data ID, 이름, PDF ID, 확실 여부를 저장합니다.
파일이 없으면 "없음"으로 표시합니다.
''')
    parser.add_argument('-u', '--usage', action='help', help='show usage, options, and CLI examples, then exit')
    parser.add_argument('folder', type=Path, help='SIESTA calculation directory or parent dataset directory')
    parser.add_argument('-g', '--group', action='store_true', help='recursively extract calculations as a horizontal table and JSON records')
    parser.add_argument('-nt', '--ntest', type=int, metavar='N', help='with -g, stop searching after N calculation folders (positive integer; default: all)')
    parser.add_argument('-o', '--output-dir', type=Path, default=Path('.'), help='output directory for metadata.json and molecule_names.tsv (default: current working directory)')
    parser.add_argument('-s', '--section', choices=('names', 'all', 'main', '1', '2', '3', '4'), default='main', help='stdout section (default: main; all shows full metadata)')
    parser.add_argument('--input', default='RUN.fdf', help='input FDF filename (default: RUN.fdf)')
    parser.add_argument('--stdout', default='stdout.txt', help='SIESTA output filename (default: stdout.txt)')
    parser.add_argument('--fdf-log', help='effective-parameter log; required if multiple logs exist')
    parser.add_argument('--database', default='TABS', help='source DB (user-supplied, default: TABS)')
    parser.add_argument('-n', '--name', help='TABS PDF/text filename for structure-matched molecule names, or a literal name')
    parser.add_argument('--multiplicity', type=int, help='known spin multiplicity; not inferred')
    parser.add_argument('--provider')
    parser.add_argument('--reference-title')
    parser.add_argument('--doi')
    args = parser.parse_args()
    if args.multiplicity is not None and args.multiplicity < 1:
        parser.error('--multiplicity must be positive')
    if args.ntest is not None:
        if args.ntest < 1:
            parser.error('--ntest must be positive')
        if not args.group:
            parser.error('--ntest requires -g / --group')
    try:
        args.name_map = None
        args.name_source = None
        if args.name:
            source = Path(args.name).expanduser()
            if source.is_file() or source.suffix.lower() in ('.pdf', '.txt'):
                args.name_map = load_names(source)
                args.name_source = str(source.resolve())
        root = args.folder.expanduser().resolve()
        if args.group:
            folders = calculation_folders(root, args.stdout, args.ntest)
            if not folders:
                raise ValueError(f'No calculation directories containing {args.stdout} found under {root}')
            return extract_dataset(args, root, folders)
        data = extract(args)
        save_records(args.output_dir, root, [{'source_directory': str(root), **data}], [])
        if args.section == 'names':
            print_name_table([{'source_directory': str(root), **data}], [])
        elif args.section == 'main':
            print_table('Main metadata', data['main_metadata'])
        elif args.section != 'all':
            detail = META_KEYS[('structure', 'calculation', 'results', 'data_source')[int(args.section) - 1]]
            selected = data['groups'][detail['group']]
            print_table(detail['title'], selected)
        else:
            print_table('Main metadata', data['main_metadata'])
            for title, payload in data['groups'].items():
                print()
                print_table(GROUP_LABELS.get(title, display_label(title)), payload)
    except (OSError, ValueError, KeyError, IndexError) as exc:
        parser.exit(2, f'error: {exc}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
