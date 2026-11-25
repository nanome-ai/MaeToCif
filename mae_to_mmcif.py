#!/usr/bin/env python3
"""Convert Maestro CT entries to one mmCIF data block with one model per entry.

Run python3 mae_to_mmcif.py --help for the command-line options.
Uses the Python 3.10+ standard library. Coordinates stay in their source frames.
Partial CTs inherit missing fields from the preceding full CT. Components with
incompatible bond templates get distinct identifiers.
_geom_bond records every bond per model. _maestro_m_bond retains its source
order. _maestro_secondary_structure records assignments for individual models.
Standard helix and strand ranges use model-specific chain IDs. _maestro_chain
maps those IDs back to the original chain names.
The other _maestro_* tables retain the original CT properties and tables.
Ligand properties also use standard chemical component fields and
_pdbx_chem_comp_model_feature for additional source values. The model number
and label_asym_id identify the ligand. Atom features use atom_site.<id>.<key>.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from itertools import chain, groupby
from math import dist, isfinite
from pathlib import Path
import re
from tempfile import NamedTemporaryFile

from mae_parser import MAEParser, MaeBlock, MaeTable


NA = object()
ELEMENTS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
            "Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn "
            "Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W "
            "Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf "
            "Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()
AMINO_ACIDS = dict(zip("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(),
                      "ACDEFGHIKLMNPQRSTVWY"))
PARENT = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "ASH": "ASP", "GLH": "GLU",
          "LYN": "LYS", "CYX": "CYS", "CYM": "CYS", "SEP": "SER", "TPO": "THR",
          "PTR": "TYR", "MSE": "MET"}
DNA = {"DA", "DC", "DG", "DT", "DI", "DU"}
RNA = {"A", "C", "G", "U", "I"}
WATER = {"D2O", "DOD", "H2O", "HOH", "SOL", "SPC", "TIP", "TIP3", "TIP4", "WAT", "WATER"}
ORDERS = {"1": "sing", "2": "doub", "3": "trip", "4": "quad"}
SECONDARY = {"-1": "unassigned", "0": "coil", "1": "alpha_helix", "2": "beta_strand"}
COORDINATES = ("r_m_x_coord", "r_m_y_coord", "r_m_z_coord")
COMPONENT_PROPERTIES = {
    "name": ("chemicalname", "compoundname", "iupacname", "name", "title"),
    "formula": ("chemicalformula", "molecularformula", "molformula", "formula"),
    "formula_weight": ("formulaweight", "molecularweight", "molweight", "molmw", "mw"),
    "pdbx_formal_charge": ("formalcharge", "netformalcharge"),
}
DESCRIPTOR_PROPERTIES = {
    "SMILES": ("smiles", "isomericsmiles"),
    "SMILES_CANONICAL": ("canonicalsmiles", "smilescanonical"),
    "InChI": ("inchi",),
    "InChIKey": ("inchikey",),
}
ATOM_SITE_PROPERTIES = {
    *COORDINATES, "i_m_atomic_number", "i_m_residue_number", "i_m_formal_charge",
    "s_m_chain_name", "s_m_insertion_code", "s_m_pdb_residue_name", "s_m_pdb_atom_name",
    "s_m_atom_name", "s_pdb_altloc_chars", "s_m_pdb_altloc", "s_m_pdb_alt_loc",
    "r_m_pdb_occupancy", "r_m_pdb_tfactor",
}


def value(row, name):
    return (row.get(name) or "").strip()


def parent(component):
    return PARENT.get(component, component)


@dataclass
class Residue:
    chain: str
    number: int
    insertion: str
    component: str
    atoms: list[Atom] = field(default_factory=list)
    kind: str = ""
    author_chain: str = ""
    asym: str = ""
    entity: int = 0
    seq: int | None = None
    variant_component: str | None = None

    @property
    def component_id(self):
        return self.variant_component or self.component

    def label(self):
        return self.component_id, self.asym, self.seq if self.seq is not None else NA

    def author(self):
        return self.component_id, self.author_chain or NA, self.number

    def secondary(self):
        codes = [(atom.label, value(atom.raw, "i_m_secondary_structure")) for atom in self.atoms]
        alpha_carbon = {code for label, code in codes if label == "CA" and code}
        assigned = alpha_carbon or {code for _, code in codes if code}
        return next(iter(assigned)) if len(assigned) == 1 else None


@dataclass
class Atom:
    site: int
    residue: Residue
    label: str
    author: str
    alt: str
    element: str | None
    xyz: tuple[float, float, float]
    raw: dict[str, str | None]


@dataclass
class Bond:
    index: int
    first: Atom
    second: Atom
    raw: dict[str, str | None]

    @property
    def order(self):
        return self.raw.get("i_m_order")


@dataclass
class Model:
    number: int
    source: MaeBlock
    seqres: dict[str, list[str]]
    residues: list[Residue]
    atoms: dict[int, Atom]
    bonds: list[Bond]
    scope_chains: bool = False

    def chain_id(self, source_chain):
        return f"{source_chain}_{self.number}" if self.scope_chains else source_chain


def resolve_partial_ct(full, partial):
    if full is None:
        raise ValueError("A partial MAE CT must follow a full CT")

    tables = full.tables | partial.tables
    for name, changes in partial.tables.items():
        base = full.tables.get(name)
        if base is None:
            continue
        columns = list(dict.fromkeys(chain(base.columns, changes.columns)))
        rows = base.rows.copy()
        for index, row in changes.rows.items():
            rows[index] = base.rows.get(index, {}) | row
        tables[name] = MaeTable(columns, rows)
    return MaeBlock(partial.name, full.properties | partial.properties, tables)


def make_model(number, block, offset):
    seqres = {}
    if table := block.tables.get("m_PDB_SEQRES"):
        for row in table.rows.values():
            seqres.setdefault(value(row, "s_pdb_chain_id"), []).extend(value(row, "s_pdb_SEQRES").split())
    residues = {}
    atoms = {}
    names = defaultdict(set)
    table = block.tables.get("m_atom")
    if table is None or not table.rows:
        raise ValueError(f"Model {number} has no atoms")
    for index, row in table.rows.items():
        key = (value(row, "s_m_chain_name"), int(value(row, "i_m_residue_number") or "1"),
               value(row, "s_m_insertion_code"), value(row, "s_m_pdb_residue_name") or "UNK")
        residue = residues.get(key)
        if residue is None:
            residue = residues[key] = Residue(*key)
        atomic_number = int(value(row, "i_m_atomic_number") or 0)
        element = ELEMENTS[atomic_number] if 0 < atomic_number < len(ELEMENTS) else None
        author = value(row, "s_m_pdb_atom_name") or value(row, "s_m_atom_name") or f"{element or 'X'}{index}"
        alt = (value(row, "s_pdb_altloc_chars") or value(row, "s_m_pdb_altloc")
               or value(row, "s_m_pdb_alt_loc"))
        label = author
        suffix = 2
        while (label, alt) in names[key]:
            label = f"{author}_{suffix}"
            suffix += 1
        names[key].add((label, alt))
        xyz = tuple(float(value(row, name)) for name in COORDINATES)
        if not all(isfinite(coordinate) for coordinate in xyz):
            raise ValueError(f"Model {number}, atom {index}: non-finite coordinate")
        atom = Atom(offset + len(atoms) + 1, residue, label, author, alt, element, xyz, row)
        atoms[index] = atom
        residue.atoms.append(atom)

    bonds = []
    if table := block.tables.get("m_bond"):
        for index, row in table.rows.items():
            first = int(value(row, "i_m_from"))
            second = int(value(row, "i_m_to"))
            if first not in atoms or second not in atoms:
                raise ValueError(f"Model {number}, bond {index}: unknown atom {first} or {second}")
            bonds.append(Bond(index, atoms[first], atoms[second], row))
    return Model(number, block, seqres, list(residues.values()), atoms, bonds)


def residue_kind(residue, sequence):
    component = residue.component
    if component in WATER:
        return "water"
    if parent(component) in AMINO_ACIDS or component in DNA | RNA:
        return "polymer"
    for atom in residue.atoms:
        index = int(value(atom.raw, "i_pdb_seqres_index") or 0)
        is_hetero = value(atom.raw, "i_ppw_het") == "1"
        if 0 < index <= len(sequence) and sequence[index - 1] == component and not is_hetero:
            return "polymer"
    return "non-polymer"


def resolve_component_bonds(models):
    components = defaultdict(list)
    reserved = {
        component.casefold() for model in models for sequence in model.seqres.values()
        for component in sequence
    }
    for model in models:
        local_bonds = defaultdict(lambda: defaultdict(set))
        for bond in model.bonds:
            first, second = bond.first, bond.second
            if first.residue is second.residue and first.label != second.label:
                pair = tuple(sorted((first.label, second.label)))
                local_bonds[id(first.residue)][pair].add(bond.order)
        for residue in model.residues:
            names = Counter(atom.label for atom in residue.atoms)
            components[residue.component].append(
                (residue, names, local_bonds[id(residue)], "m_bond" in model.source.tables)
            )
            reserved.add(residue.component.casefold())

    templates = []
    for component, instances in sorted(components.items()):
        pairs = {pair for _, _, bonds, _ in instances for pair in bonds}
        variants = []
        for residue, names, bonds, has_bond_table in instances:
            observed = {}
            ambiguous = set()
            for pair in pairs:
                first, second = pair
                if not has_bond_table or first not in names or second not in names:
                    continue
                orders = bonds.get(pair, ())
                order = next(iter(orders), "0")
                if (names[first] != 1 or names[second] != 1 or len(orders) > 1
                        or (order not in ORDERS and order != "0")):
                    ambiguous.add(pair)
                    continue
                observed[pair] = order

            for identifier, known, blocked in variants:
                if all(pair not in known or known[pair] == order for pair, order in observed.items()):
                    break
            else:
                identifier = component
                if variants:
                    suffix = len(variants) + 1
                    identifier = f"{component}_{suffix}"
                    while identifier.casefold() in reserved:
                        suffix += 1
                        identifier = f"{component}_{suffix}"
                    reserved.add(identifier.casefold())
                known, blocked = {}, set()
                variants.append((identifier, known, blocked))
            known.update(observed)
            blocked.update(ambiguous)
            if identifier != component:
                residue.variant_component = identifier

        for identifier, known, blocked in variants:
            for (first, second), order in known.items():
                if order in ORDERS and (first, second) not in blocked:
                    templates.append((identifier, first, second, ORDERS[order], None))
    return sorted(templates)


def assign_entities(models):
    sequences = defaultdict(lambda: defaultdict(dict))
    polymers = defaultdict(list)
    for model in models:
        for chain, sequence in model.seqres.items():
            for position, component in enumerate(sequence, 1):
                sequences[chain][position][component] = None
        for residue in model.residues:
            sequence = model.seqres.get(residue.chain, [])
            residue.kind = residue_kind(residue, sequence)
            if residue.kind == "polymer":
                polymers[residue.chain].append(residue)

    for chain, residues in polymers.items():
        positions = sequences[chain]
        observed_keys = sorted({(r.number, r.insertion) for r in residues})
        if not positions:
            for position, key in enumerate(observed_keys, 1):
                positions[position] = {}
            fallback = dict(zip(observed_keys, positions))
        else:
            fallback = {}
        for residue in residues:
            indices = {int(index) for atom in residue.atoms
                       if (index := value(atom.raw, "i_pdb_seqres_index")) and int(index) > 0}
            if len(indices) > 1:
                raise ValueError(f"Conflicting SEQRES indices in chain {chain} residue {residue.number}")
            if fallback:
                residue.seq = fallback[residue.number, residue.insertion]
            elif indices:
                residue.seq = next(iter(indices))

        known = defaultdict(set)
        for residue in residues:
            if residue.seq is not None:
                known[residue.number, residue.insertion].add(residue.seq)
        for residue in residues:
            if residue.seq is None:
                candidates = known[residue.number, residue.insertion]
                if len(candidates) == 1:
                    residue.seq = next(iter(candidates))
                elif not residue.insertion and residue.component in positions.get(residue.number, {}):
                    residue.seq = residue.number
            if residue.seq is not None:
                positions[residue.seq][residue.component_id] = None

    entities = {}
    asyms = {}
    for chain_name in sequences:
        entities["polymer", chain_name] = len(entities) + 1
    for model in models:
        model.scope_chains = len(models) > 1
        for residue in model.residues:
            residue.author_chain = model.chain_id(residue.chain)
            if residue.kind == "polymer":
                key = (residue.chain,)
                entity_key = ("polymer", residue.chain)
            else:
                entity_key = (residue.kind, residue.component_id)
                key = (residue.chain, residue.number, residue.insertion, residue.component_id)
            if entity_key not in entities:
                entities[entity_key] = len(entities) + 1
            asym_key = (model.number, residue.kind, key)
            if asym_key not in asyms:
                asyms[asym_key] = (asym_id(len(asyms)), entities[entity_key])
            residue.asym, residue.entity = asyms[asym_key]
    return entities, asyms, sequences


def asym_id(index):
    identifier = ""
    while index >= 0:
        index, letter = divmod(index, 26)
        identifier = chr(ord("A") + letter) + identifier
        index -= 1
    return identifier


def cif_value(value):
    if value is None:
        return "?"
    if value is NA:
        return "."
    value = str(value)
    if (value and value not in {".", "?"} and not any(char.isspace() for char in value)
            and value[0] not in "_#$;'\"[]"
            and not value.lower().startswith(("data_", "save_", "loop_", "stop_", "global_"))):
        return value
    if "\n" not in value and "\r" not in value:
        for quote in ("'", '"'):
            if quote not in value:
                return quote + value + quote
    if value.startswith(";") or "\n;" in value:
        raise ValueError("A Maestro string contains an mmCIF text-field delimiter")
    return f"\n;{value}\n;\n"


def write_loop(output, category, columns, rows):
    rows = iter(rows)
    first = next(rows, None)
    if first is None:
        return
    output.write("loop_\n")
    output.writelines(f"_{category}.{column}\n" for column in columns)
    for row in chain((first,), rows):
        output.write(" ".join(map(cif_value, row)) + "\n")
    output.write("#\n")


def common_value(values):
    return next(iter(values)) if len(values) == 1 else None


def numeric_value(raw):
    if raw is None:
        return None
    try:
        number = Decimal(raw)
    except InvalidOperation:
        return None
    return number if number.is_finite() else None


def formal_charge(raw):
    return int(raw) if raw is not None and re.fullmatch(r"[+-]?[0-9]+", raw.strip()) else None


def ligand_properties(models):
    for model in models:
        ligands = [residue for residue in model.residues if residue.kind == "non-polymer"]
        standalone = len(model.residues) == 1
        properties = {}
        if len(ligands) == 1:
            for key, raw in model.source.properties.items():
                name = key.split("_", 2)[-1].replace("_", "").casefold()
                if raw and raw.strip() and (standalone or name.startswith("ligand")):
                    properties[key] = (name.removeprefix("ligand"), raw)
        for ligand in ligands:
            yield model.number, ligand, properties


def component_metadata(ligands):
    values = defaultdict(lambda: defaultdict(set))
    descriptors = defaultdict(set)
    for _, ligand, properties in ligands:
        fields = defaultdict(set)
        for name, raw in properties.values():
            if raw.strip() not in {".", "?"}:
                fields[name].add(raw.strip())

        for field, aliases in COMPONENT_PROPERTIES.items():
            for alias in aliases:
                candidates = set()
                for raw in fields.get(alias, ()):
                    if field == "pdbx_formal_charge":
                        result = formal_charge(raw)
                    elif field == "formula_weight":
                        mass = numeric_value(raw)
                        result = mass if mass is not None and mass >= 1 else None
                    else:
                        result = raw
                    if result is not None:
                        candidates.add(result)
                if candidates:
                    values[ligand.component_id][field].update(candidates)
                    break

        for kind, aliases in DESCRIPTOR_PROPERTIES.items():
            for alias in aliases:
                if not (source_values := fields.get(alias)):
                    continue
                program = (common_value(fields.get(alias + "program", ())) or
                           common_value(fields.get("descriptorprogram", ())) or "OTHER")
                version = (common_value(fields.get(alias + "programversion", ())) or
                           common_value(fields.get("descriptorprogramversion", ())) or "unknown")
                descriptors[ligand.component_id, kind, program, version].update(source_values)

    return (
        {component: {field: common_value(candidates) for field, candidates in fields.items()}
         for component, fields in values.items()},
        [(*key, common_value(candidates)) for key, candidates in sorted(descriptors.items())
         if len(candidates) == 1],
    )


def write_ligand_features(output, ligands):
    def rows():
        for model_number, ligand, properties in ligands:
            identifier = f"M_{model_number}_{ligand.asym}"
            yield identifier, "pdbx_PDB_model_num", model_number
            yield identifier, "label_asym_id", ligand.asym
            if ligand.variant_component:
                yield identifier, "source_comp_id", ligand.component
            for key, (_, raw) in properties.items():
                yield identifier, key, raw
            for atom in ligand.atoms:
                for key, raw in atom.raw.items():
                    if key not in ATOM_SITE_PROPERTIES and raw is not None and raw.strip():
                        yield identifier, f"atom_site.{atom.site}.{key}", raw

    write_loop(output, "pdbx_chem_comp_model", ("id", "comp_id"),
               ((f"M_{model_number}_{ligand.asym}", ligand.component_id)
                for model_number, ligand, _ in ligands))
    write_loop(output, "pdbx_chem_comp_model_feature", ("model_id", "feature_name", "feature_value"), rows())


def write_entities(output, entities, asyms, sequences, models, metadata):
    polymer_rows, sequence_rows, nonpoly_rows = [], [], []
    component_types = {}
    polymer_sources = {}
    chain_ids = defaultdict(list)
    for model in models:
        for name in dict.fromkeys(residue.chain for residue in model.residues if residue.kind == "polymer"):
            chain_ids[name].append(model.chain_id(name))
        polymer_sources.update(
            (residue.variant_component, residue.component) for residue in model.residues
            if residue.kind == "polymer" and residue.variant_component
        )
    for (kind, name), entity in entities.items():
        if kind != "polymer":
            component_name = metadata.get(name, {}).get("name") or name
            nonpoly_rows.append((entity, "water" if kind == "water" else component_name, name))
            component_types.setdefault(name, "non-polymer")
            continue
        positions = sequences[name]
        primary = [next(iter(options), "UNK") for _, options in sorted(positions.items())]
        primary = [polymer_sources.get(component, component) for component in primary]
        components = {monomer for options in positions.values() for monomer in options}
        parents = {parent(polymer_sources.get(component, component)) for component in components}
        if parents & AMINO_ACIDS.keys():
            polymer_type = "polypeptide(L)"
            component_type = "L-peptide linking"
            canonical = "".join(AMINO_ACIDS.get(parent(component), "X") for component in primary)
        elif parents <= DNA | RNA:
            polymer_type = ("polydeoxyribonucleotide" if parents <= DNA else "polyribonucleotide"
                            if parents <= RNA else "polydeoxyribonucleotide/polyribonucleotide hybrid")
            component_type = None
            canonical = "".join(component[-1] for component in primary)
        else:
            polymer_type, component_type, canonical = "other", "other", None

        for component in components:
            original = polymer_sources.get(component, component)
            component_types[component] = (
                "peptide linking" if original == "GLY" else component_type or
                ("DNA linking" if original in DNA else "RNA linking")
            )
        strand_ids = ",".join(chain_ids[name]) or NA
        polymer_rows.append((entity, polymer_type, strand_ids, canonical))
        for position, options in sorted(positions.items()):
            sequence_rows.extend(
                (entity, position, component, "y" if len(options) > 1 else "n")
                for component in options
            )

    for variant, original in polymer_sources.items():
        component_types.setdefault(variant, component_types.get(original, "other"))
    parent_components = PARENT | {variant: parent(original) for variant, original in polymer_sources.items()}
    write_loop(output, "entity", ("id", "type", "pdbx_description"),
               ((entity, kind, f"Chain {name or '(blank)'}" if kind == "polymer" else name)
                for (kind, name), entity in entities.items()))
    write_loop(output, "struct_asym", ("id", "entity_id"), asyms.values())
    write_loop(output, "entity_poly", ("entity_id", "type", "pdbx_strand_id", "pdbx_seq_one_letter_code_can"),
               polymer_rows)
    write_loop(output, "entity_poly_seq", ("entity_id", "num", "mon_id", "hetero"), sequence_rows)
    write_loop(output, "pdbx_entity_nonpoly", ("entity_id", "name", "comp_id"), nonpoly_rows)
    columns = tuple(field for field in COMPONENT_PROPERTIES
                    if any(values.get(field) is not None for values in metadata.values()))
    write_loop(output, "chem_comp", ("id", "type", "mon_nstd_parent_comp_id", *columns),
               ((component, kind, parent_components.get(component, NA),
                 *(metadata.get(component, {}).get(field) for field in columns))
                for component, kind in sorted(component_types.items())))


def ordered_residues(model):
    return sorted(model.residues, key=lambda residue: (
        residue.asym, residue.number, residue.insertion, residue.component
    ))


def atom_sites(models):
    # Keep each residue contiguous for readers that apply ranges in file order.
    for model in models:
        for residue in ordered_residues(model):
            for atom in residue.atoms:
                row = atom.raw
                group = "ATOM" if residue.kind == "polymer" and value(row, "i_ppw_het") != "1" else "HETATM"
                occupancy = row.get("r_m_pdb_occupancy")
                yield (group, atom.site, atom.element, atom.label, atom.alt or NA, residue.component_id,
                       residue.asym, residue.entity, residue.seq if residue.seq is not None else NA,
                       residue.insertion or None, *(row[axis] for axis in COORDINATES),
                       occupancy if occupancy is not None else -1, row.get("r_m_pdb_tfactor"),
                       row.get("i_m_formal_charge"), residue.number, residue.component_id,
                       residue.author_chain or NA, atom.author, model.number)


ATOM_COLUMNS = (
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
    "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code",
    "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
    "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num",
)


def write_bonds(output, models, templates):
    component_atoms = defaultdict(set)
    formal_charges = defaultdict(set)
    partial_charges = defaultdict(set)
    for model in models:
        for residue in model.residues:
            for atom in residue.atoms:
                key = (residue.component_id, atom.label)
                component_atoms[key].add(atom.element)
                if residue.kind == "non-polymer":
                    if (charge := formal_charge(atom.raw.get("i_m_formal_charge"))) is not None:
                        formal_charges[key].add(charge)
                    if (partial := numeric_value(atom.raw.get("r_m_charge1"))) is not None:
                        partial_charges[key].add(partial)
    charge_columns = tuple((field, charges) for field, charges in
                           (("charge", formal_charges), ("partial_charge", partial_charges))
                           if any(len(values) == 1 for values in charges.values()))
    write_loop(output, "chem_comp_atom", ("comp_id", "atom_id", "type_symbol",
                                         *(field for field, _ in charge_columns)),
               ((*key, common_value(elements), *(common_value(charges.get(key, ()))
                                                 for _, charges in charge_columns))
                for key, elements in sorted(component_atoms.items())))

    write_loop(output, "chem_comp_bond", ("comp_id", "atom_id_1", "atom_id_2", "value_order",
                                        "pdbx_aromatic_flag"), templates)

    write_loop(output, "geom_bond", ("pdbx_PDB_model_num", "atom_site_id_1", "atom_site_id_2",
                                    "dist", "site_symmetry_1", "site_symmetry_2"),
               ((model.number, bond.first.site, bond.second.site, f"{dist(bond.first.xyz, bond.second.xyz):.6f}",
                 "1_555", "1_555") for model in models for bond in model.bonds))
    write_connections(output, models)


def is_backbone_link(first, second):
    if first.residue.asym != second.residue.asym or first.residue.seq is None or second.residue.seq is None:
        return False
    if first.residue.seq > second.residue.seq:
        first, second = second, first
    return (second.residue.seq == first.residue.seq + 1
            and (first.label, second.label) in {("C", "N"), ("O3'", "P")})


def write_connections(output, models):
    def endpoint(atom):
        comp, asym, seq = atom.residue.label()
        return (asym, comp, seq, atom.label, atom.alt or NA, atom.residue.insertion or None,
                atom.residue.author_chain or NA, atom.residue.number, comp, "1_555")

    def rows():
        for model in models:
            for bond in model.bonds:
                a, b = bond.first, bond.second
                if a.residue is b.residue or (bond.order == "1" and is_backbone_link(a, b)):
                    continue
                kind = "disulf" if (a.residue.component, a.label, b.residue.component, b.label) == (
                    "CYS", "SG", "CYS", "SG") else "covale"
                yield (f"conn{model.number}_{bond.index}", kind, *endpoint(a), *endpoint(b), ORDERS.get(bond.order))

    fields = ("label_asym_id", "label_comp_id", "label_seq_id", "label_atom_id")
    columns = ["id", "conn_type_id"]
    for i in (1, 2):
        columns.extend(f"ptnr{i}_{name}" for name in fields)
        columns.extend((f"pdbx_ptnr{i}_label_alt_id", f"pdbx_ptnr{i}_PDB_ins_code",
                        f"ptnr{i}_auth_asym_id", f"ptnr{i}_auth_seq_id", f"ptnr{i}_auth_comp_id",
                        f"ptnr{i}_symmetry"))
    columns.append("pdbx_value_order")
    connections = list(rows())
    write_loop(output, "struct_conn_type", ("id",), sorted({(row[1],) for row in connections}))
    write_loop(output, "struct_conn", columns, connections)


def secondary_segments(annotations):
    for (_, code), annotated in groupby(annotations, key=lambda pair: (pair[0].asym, pair[1])):
        if code not in {"1", "2"}:
            continue
        segment = []
        for residue, _ in annotated:
            if segment:
                previous = segment[-1]
                seq_gap = (previous.seq is not None and residue.seq is not None
                           and residue.seq != previous.seq + 1)
                if seq_gap or residue.number - previous.number not in (0, 1):
                    yield code, segment
                    segment = []
            segment.append(residue)
        if segment:
            yield code, segment


def write_secondary_structure(output, models):
    annotations = [
        [(residue, residue.secondary()) for residue in ordered_residues(model) if residue.kind == "polymer"]
        for model in models
    ]
    write_loop(output, "maestro_secondary_structure",
               ("pdbx_PDB_model_num", "label_comp_id", "label_asym_id", "label_seq_id",
                "auth_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "code", "type"),
               ((model.number, *residue.label(), *residue.author(), residue.insertion or None,
                 code, SECONDARY.get(code))
                for model, residues in zip(models, annotations) for residue, code in residues))

    helices, sheets, strands = [], [], []
    for residues in annotations:
        for code, segment in secondary_segments(residues):
            first, last = segment[0], segment[-1]
            endpoints = (*first.label(), first.insertion or None, *last.label(), last.insertion or None,
                         *first.author(), *last.author())
            if code == "1":
                identifier = f"H{len(helices) + 1}"
                helices.append(("HELX_P", identifier, identifier, *endpoints, 1, len(segment)))
            else:
                identifier = f"S{len(sheets) + 1}"
                sheets.append((identifier, 1))
                strands.append((identifier, 1, *endpoints))
    columns = (
        "beg_label_comp_id", "beg_label_asym_id", "beg_label_seq_id", "pdbx_beg_PDB_ins_code",
        "end_label_comp_id", "end_label_asym_id", "end_label_seq_id", "pdbx_end_PDB_ins_code",
        "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
        "end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id",
    )
    if helices:
        output.write("_struct_conf_type.id HELX_P\n#\n")
    write_loop(output, "struct_conf", ("conf_type_id", "id", "pdbx_PDB_helix_id", *columns,
                                      "pdbx_PDB_helix_class", "pdbx_PDB_helix_length"), helices)
    write_loop(output, "struct_sheet", ("id", "number_strands"), sheets)
    write_loop(output, "struct_sheet_range", ("sheet_id", "id", *columns), strands)


def write_properties(output, headers, models):
    properties = dict.fromkeys(name for model in models for name in model.source.properties)
    write_loop(output, "maestro_model", ("pdbx_PDB_model_num", "ct_block", *properties),
               ((model.number, model.source.name, *(model.source.properties.get(name) for name in properties))
                for model in models))
    write_loop(output, "maestro_chain", ("pdbx_PDB_model_num", "auth_asym_id", "mae_chain_name"),
               ((model.number, model.chain_id(chain_name) or NA, chain_name)
                for model in models for chain_name in dict.fromkeys(r.chain for r in model.residues)))
    header_fields = dict.fromkeys(name for block in headers for name in block.properties)
    write_loop(output, "maestro_header", ("id", *header_fields),
               ((index, *(block.properties.get(name) for name in header_fields))
                for index, block in enumerate(headers, 1)))

    tables = defaultdict(list)
    for model in models:
        for name, table in model.source.tables.items():
            tables[name].append((model, table))
    for name, members in tables.items():
        fields = dict.fromkeys(column for _, table in members for column in table.columns)
        links = {
            "m_atom": ("atom_site_id",),
            "m_bond": ("atom_site_id_1", "atom_site_id_2"),
        }.get(name, ())

        def rows():
            for model, table in members:
                for index, row in table.rows.items():
                    if name == "m_atom":
                        sites = (model.atoms[index].site,)
                    elif name == "m_bond":
                        sites = (model.atoms[int(row["i_m_from"])].site, model.atoms[int(row["i_m_to"])].site)
                    else:
                        sites = ()
                    yield model.number, index, *sites, *(row.get(column) for column in fields)

        write_loop(output, f"maestro_{name}", ("pdbx_PDB_model_num", "mae_index", *links, *fields), rows())


def convert(source, destination):
    if source.resolve() == destination.resolve():
        raise ValueError("The input and output must be different files")
    blocks = MAEParser().parse(source)
    headers, models = [], []
    last_full = None
    atom_count = 0
    for block in blocks:
        if not block.name:
            headers.append(block)
            continue
        if block.name == "f_m_ct":
            last_full = block
        else:
            block = resolve_partial_ct(last_full, block)
        model = make_model(len(models) + 1, block, atom_count)
        models.append(model)
        atom_count += len(model.atoms)
    if not models:
        raise ValueError("No Maestro CT entries found")

    templates = resolve_component_bonds(models)
    entities, asyms, sequences = assign_entities(models)
    ligands = list(ligand_properties(models))
    metadata, descriptors = component_metadata(ligands)
    entry_id = re.sub(r"[^A-Za-z0-9_-]+", "_", source.stem).strip("_") or "maestro"
    with NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n", dir=destination.parent,
                            prefix=f".{destination.name}.", delete=False) as output:
        temporary = Path(output.name)
        try:
            output.write(f"data_{entry_id}\n_entry.id {cif_value(entry_id)}\n#\n")
            experiment = {model.source.properties.get("s_pdb_PDB_EXPDTA") for model in models}
            if len(experiment) == 1 and (method := next(iter(experiment))):
                write_loop(output, "exptl", ("entry_id", "method"), [(entry_id, method)])
            write_entities(output, entities, asyms, sequences, models, metadata)
            write_loop(output, "pdbx_chem_comp_descriptor",
                       ("comp_id", "type", "program", "program_version", "descriptor"), descriptors)
            write_loop(output, "atom_site", ATOM_COLUMNS, atom_sites(models))
            write_bonds(output, models, templates)
            write_ligand_features(output, ligands)
            write_secondary_structure(output, models)
            write_properties(output, headers, models)
            output.close()
            temporary.replace(destination)
        finally:
            temporary.unlink(missing_ok=True)
    return len(models), atom_count, sum(len(model.bonds) for model in models)


def main(argv=None):
    parser = argparse.ArgumentParser(
        usage="%(prog)s [-h] (SOURCE | --input SOURCE) [-o DESTINATION]",
        description=(
            "Convert a Maestro file to one mmCIF with one model per CT entry.\n"
            "The export preserves coordinates, bonds and secondary structures.\n"
            "Supply the input as SOURCE or with --input."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples\n"
            "  python3 %(prog)s ligand.mae\n"
            "  python3 %(prog)s --input poses.maegz --output ensemble.cif"
        ),
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("source", nargs="?", type=Path, metavar="SOURCE",
                        help="input .mae, .maegz or .mae.gz file")
    inputs.add_argument("-i", "--input", type=Path, metavar="SOURCE",
                        help="input .mae, .maegz or .mae.gz file")
    parser.add_argument("-o", "--output", type=Path, metavar="DESTINATION",
                        help="output path (default: replace the input's last extension with .cif)")
    args = parser.parse_args(argv)
    source = args.input or args.source
    destination = args.output or source.with_suffix(".cif")
    try:
        models, atoms, bonds = convert(source, destination)
    except (KeyError, OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")
    print(f"{destination}: {models} models, {atoms} atoms, {bonds} bonds")


if __name__ == "__main__":
    main()
