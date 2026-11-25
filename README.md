# Maestro to mmCIF

Convert a Schrödinger Maestro file (`.mae`, `.maegz` or `.mae.gz`) into one mmCIF data block with one model per Maestro entry.

The converter uses the Python 3.10+ standard library and the included `mae_parser.py`.

## Usage

```bash
python mae_to_mmcif.py input.mae
python mae_to_mmcif.py --input poses.maegz --output ensemble.cif
python mae_to_mmcif.py --help
```

Provide the input path as a positional argument or with `-i` / `--input`. Without `-o` / `--output`, the script writes beside the input and replaces its final extension with `.cif`.

## Export details

- Preserves coordinates, bond orders and per-model secondary structure. Partial Maestro entries (`p_m_ct`) inherit fields from the preceding full entry (`f_m_ct`).
- Uses model-specific chain IDs for secondary structure and distinct component IDs for incompatible bond templates. The mmCIF retains the original Maestro identifiers.
- Stores ligand properties in standard chemical component fields and `_pdbx_chem_comp_model_feature`. Parsed Maestro properties are also retained in `_maestro_*` categories.
- Writes `-1` for missing occupancies for Nanome compatibility.
