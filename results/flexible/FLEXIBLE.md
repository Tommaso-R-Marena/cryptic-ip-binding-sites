# Study J: does a flexible receptor recover study A's failures?

Pre-registered in `docs/FLEXIBLE_PLAN.md`, amended in `docs/FLEXIBLE_PLAN_AMENDMENT_1.md`, `_2.md` and `_3.md`. Seed 20261003, 2000 resamples of strict homology groups, exhaustiveness 32.
Flexible side chains: within 4.0 A of the ligand, at most 8.

Movable side chains per copy: mean 6.99, max 8, and 2 copies had none.

## J1 - top-pose success at 2 A (primary)

Decision: **worse**.
flex 0.026 against rigid 0.095; paired difference -0.078 [-0.140, -0.031] over 207 copies in 29 groups.

## J2 - best-of-list ceiling (secondary)

Decision: **worse**.
flex 0.177 against rigid 0.364; paired difference -0.169 [-0.243, -0.098] over 207 copies in 29 groups.

## The receptor pipeline, audited

The rigid arm scores 0.095 under this study's receptor, against study F's 0.110 under the project-prepared one, over 207 copies. A material gap is a finding about the pipeline, not about side-chain freedom, and is not folded into J1.

## J3 - by burial class (exploratory, not Holm-corrected)

| burial class | decision | difference | copies | groups |
| --- | --- | --- | --- | --- |
| cryptic | not evaluable | -0.155 [-0.500, +0.179] | 10 | 4 |
| semi_cryptic | worse | -0.137 [-0.226, -0.053] | 44 | 11 |
| surface | worse | -0.022 [-0.044, -0.005] | 153 | 22 |

## Receptor residues deleted to satisfy Meeko

93 residue(s) across 30 copies, each more than 8.0 A beyond the docking box, where Vina's own interaction cutoff puts them out of reach of any pose. Both arms of a copy see the same deletions or the copy fails. Per docs/FLEXIBLE_PLAN_AMENDMENT_2.md.

## Copies not docked, by cause

| cause | copies | examples |
| --- | --- | --- |
| ReceptorError | 19 | 5ED1:A:801, 5HP2:A:801, 5Y88:A:3000 |
| timed out after 18000 s | 17 | 5DGH:A:402, 6KWY:c:1901, 7MOF:A:501 |
| RuntimeError | 12 | 1FHW:B:1002, 5J16:A:303, 9HZG:A:701 |
| PolymerCreationError | 9 | 5ICN:A:401, 9D5K:A:801, 9D5K:B:801 |
| AtomValenceException | 2 | 6J6G:A:3000, 2P1N:B:601 |
| TypeError | 2 | 6QW6:5A:2401, 9DTR:A:2500 |
| ValueError | 2 | 7XZI:A:2001, 7XZJ:A:2001 |

### ReceptorError, at 5ED1:A:801

```
flexres)
                             ^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/scripts/flexible.py", line 574, in dock_copy
    ctx.receptor()
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/scripts/redocking.py", line 340, in receptor
    prepared = prepare_receptor(self.arrays, self.work, keep_metals=metals, protonated=protonated, name=name)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/cryptic_ip/docking/receptor.py", line 299, in prepare_receptor
    protonated = protonate(workdir / "polymer.pdb", workdir, ph=ph)
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/cryptic_ip/docking/receptor.py", line 174, in protonate
    raise ReceptorError(f"pdb2pqr failed: {' | '.join(tail)}"[:400])
cryptic_ip.docking.receptor.ReceptorError: pdb2pqr failed:     newcoords = hatom.coords |                 ^^^^^^^^^^^^ | AttributeError: 'NoneType' object has no attribute 'coords'
```

### RuntimeError, at 1FHW:B:1002

```
PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:
Residue A:265 matched with template 'PRO' has H discrepancy: 0 missing, 2 excess. 
These discrepancies may compromise the validity of the charge assignment from PQR, making the charges inapplicable to the processed receptor. 



The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/scripts/flexible.py", line 623, in cmd_dock_one
    runs, applied, dropped = dock_copy(ctx, flexres)
                             ^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/scripts/flexible.py", line 584, in dock_copy
    rigid_only, none_flex, _, dropped_rigid = prepare_pair(protonated, ctx.work / "arm_rigid",
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/scripts/flexible.py", line 289, in prepare_pair
    raise RuntimeError(
RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:265
```

### PolymerCreationError, at 5ICN:A:401

```
ptic-ip-binding-sites/scripts/flexible.py", line 281, in prepare_pair
    polymer = Polymer.from_pqr_string(text, mk_prep=mk, **options)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/meeko/polymer.py", line 1829, in from_pqr_string
    handle_parsing_situations(
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/meeko/polymer.py", line 908, in handle_parsing_situations
    raise PolymerCreationError(err, recs)
meeko.polymer.PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res C:18 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match templates, and --default_altloc to set
a default altloc variant. Use these at your own risk.

2. (processing individual structure) Inspecting and fixing the input structure is recommended.
Use --wanted_altloc to set variants for specific residues.


```

### AtomValenceException, at 6J6G:A:3000

```
_pair(protonated, ctx.work / "arm_rigid",
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/scripts/flexible.py", line 281, in prepare_pair
    polymer = Polymer.from_pqr_string(text, mk_prep=mk, **options)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/meeko/polymer.py", line 1769, in from_pqr_string
    tmp_raw_input_mols = cls._pqr_to_residue_mols(
                         ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/meeko/polymer.py", line 2695, in _pqr_to_residue_mols
    pdbmol, _, missed_altloc, needed_altloc = _aux_altloc_mol_build(
                                              ^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/meeko/utils/rdkitutils.py", line 415, in _aux_altloc_mol_build
    _ = Chem.SanitizeMol(pdbmol)
        ^^^^^^^^^^^^^^^^^^^^^^^^
rdkit.Chem.rdchem.AtomValenceException: Explicit valence for atom # 1 C, 5, is greater than permitted
```

### TypeError, at 6QW6:5A:2401

```
.py", line 601, in dock_copy
    result = dock_arm(rigid, flex, to_pdbqt(pose.mol), ctx.centre, size, seed=seed)
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/runner/work/cryptic-ip-binding-sites/cryptic-ip-binding-sites/scripts/flexible.py", line 348, in dock_arm
    v.set_receptor(rigid_pdbqt_filename=str(rigid))
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/vina/vina.py", line 157, in set_receptor
    self._vina.set_receptor(rigid_pdbqt_filename)
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/vina/vina_wrapper.py", line 704, in set_receptor
    return _vina_wrapper.Vina_set_receptor(self, *args)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: 

PDBQT parsing error: Coordinate "5 198.33" is not valid.
 > ATOM  100000  N   LYS j 291     169.943 235.945 198.335  1.00  0.00    -0.344 N 

Additional information:
Wrong number or type of arguments for overloaded function 'Vina_set_receptor'.
  Possible C/C++ prototypes are:
    Vina::set_receptor(std::string const &,std::string const &)
    Vina::set_receptor(std::string const &)
    Vina::set_receptor()

```

### ValueError, at 7XZI:A:2001

```
eko/polymer.py", line 1769, in from_pqr_string
    tmp_raw_input_mols = cls._pqr_to_residue_mols(
                         ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/meeko/polymer.py", line 2685, in _pqr_to_residue_mols
    raise ValueError(msg)
ValueError: each residue key must have exactly 1 resname
but got violations={':3': {'LEU', 'ALA', 'LYN', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLN', 'GLU', 'GLY', 'TRP', 'LYS', 'THR', 'ASH', 'ILE', 'TYR', 'ASN', 'ARG', 'MET', 'GLH', 'ASP', 'PRO'}, ':4': {'LEU', 'ALA', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLU', 'GLN', 'GLY', 'TRP', 'LYS', 'THR', 'CYS', 'ILE', 'TYR', 'ASN', 'ARG', 'MET', 'ASP', 'PRO'}, ':5': {'LEU', 'ALA', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLY', 'GLN', 'GLU', 'TRP', 'LYS', 'THR', 'ILE', 'TYR', 'ASN', 'ARG', 'GLH', 'ASP', 'PRO'}, ':7': {'LEU', 'ALA', 'LYN', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLU', 'GLY', 'GLN', 'TRP', 'LYS', 'THR', 'CYS', 'ILE', 'TYR', 'ASN', 'ARG', 'MET', 'GLH', 'ASP', 'PRO'}, ':9': {'LEU', 'ALA', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLU', 'GLN', 'GLY', 'HIP', 'TRP', 'LYS', 'THR', 'CYS', 'ASH', 'ILE', 'TYR', 'ASN', 'ARG', 'MET', 'GLH', 'ASP', 'PRO'}}
```

## Copies not docked

- AtomValenceException: Explicit valence for atom # 1 C, 5, is greater than permitted
- AtomValenceException: Explicit valence for atom # 4 C, 5, is greater than permitted
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res A:130 is within radius of box edge (6.3553 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match templates, and --default_altloc to set
a default altloc variant. Use these at your own risk.

2. (pro
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res A:531 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res A:534 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match temp
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res B:146 is within radius of box edge (0.4384 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res D:146 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match temp
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res B:146 is within radius of box edge (0.8325 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res D:146 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match temp
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res B:146 is within radius of box edge (1.9982 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res D:146 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match temp
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res B:320 is within radius of box edge (6.4847 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res B:341 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res B:345 is within radius of box edge (0.3257 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Traceback (most recent c
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res B:36 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res B:106 is within radius of box edge (4.2191 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match templ
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res B:656 is within radius of box edge (0.4456 <= delete_bad_res_from_box_radius=8.0000 A).

Bad res B:698 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match temp
- PolymerCreationError: 
Error: Creation of data structure for receptor failed.

Details:

Bad res C:18 is within radius of box edge (0.0000 <= delete_bad_res_from_box_radius=8.0000 A).
These residues can be ignored with option --delete_bad_res or --delete_bad_res_from_box_radius.

Recommendations:
1. (for batch processing) Use -a/--allow_bad_res to automatically remove residues
that do not match templates, and --default_altloc to set
a default altloc variant. Use these at your own risk.

2. (proc
- ReceptorError: PDB2PQR outputs disagree: 133912 atoms, 133532 charges
- ReceptorError: PDB2PQR outputs disagree: 136188 atoms, 135808 charges
- ReceptorError: PDB2PQR outputs disagree: 137887 atoms, 137507 charges
- ReceptorError: PDB2PQR outputs disagree: 141708 atoms, 141148 charges
- ReceptorError: PDB2PQR outputs disagree: 141887 atoms, 141367 charges
- ReceptorError: pdb2pqr failed:     coords = [bondatom.coords, nextatom.coords] |                                ^^^^^^^^^^^^^^^ | AttributeError: 'NoneType' object has no attribute 'coords'
- ReceptorError: pdb2pqr failed:     newcoords = hatom.coords |                 ^^^^^^^^^^^^ | AttributeError: 'NoneType' object has no attribute 'coords'
- ReceptorError: pdb2pqr failed:   File "/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/pdb2pqr/main.py", line 802, in main_driver |     raise RuntimeError from err | RuntimeError
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:174
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:188
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:188, A:82, B:625
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:265
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:366
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:378
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:465
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:4666
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at A:581
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at B:265
- RuntimeError: hydrogens disagree with Meeko's template within 8.0 A of the box, at B:96
- TypeError: 

PDBQT parsing error: Coordinate "1 266.70" is not valid.
 > ATOM  100000  N   LYS u 107     329.345 153.631 266.701  1.00  0.00    -0.344 N 

Additional information:
Wrong number or type of arguments for overloaded function 'Vina_set_receptor'.
  Possible C/C++ prototypes are:
    Vina::set_receptor(std::string const &,std::string const &)
    Vina::set_receptor(std::string const &)
    Vina::set_receptor()

- TypeError: 

PDBQT parsing error: Coordinate "5 198.33" is not valid.
 > ATOM  100000  N   LYS j 291     169.943 235.945 198.335  1.00  0.00    -0.344 N 

Additional information:
Wrong number or type of arguments for overloaded function 'Vina_set_receptor'.
  Possible C/C++ prototypes are:
    Vina::set_receptor(std::string const &,std::string const &)
    Vina::set_receptor(std::string const &)
    Vina::set_receptor()

- ValueError: each residue key must have exactly 1 resname
but got violations={':3': {'LEU', 'ALA', 'LYN', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLN', 'GLU', 'GLY', 'TRP', 'LYS', 'THR', 'ASH', 'ILE', 'TYR', 'ASN', 'ARG', 'MET', 'GLH', 'ASP', 'PRO'}, ':4': {'LEU', 'ALA', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLU', 'GLN', 'GLY', 'TRP', 'LYS', 'THR', 'CYS', 'ILE', 'TYR', 'ASN', 'ARG', 'MET', 'ASP', 'PRO'}, ':5': {'LEU', 'ALA', 'PHE', 'SER', 'HIE', 'HID', 'VAL', 'GLY', 'GLN', 'GLU', 'TRP', 'LYS', 'THR', '
- ValueError: each residue key must have exactly 1 resname
but got violations={':3': {'VAL', 'TYR', 'GLY', 'HID', 'PHE', 'ASP', 'ARG', 'LEU', 'THR', 'SER', 'LYS', 'GLN', 'GLU', 'ALA', 'MET', 'PRO'}, ':4': {'CYS', 'HID', 'TRP', 'LYS', 'TYR', 'GLY', 'LYN', 'ARG', 'LEU', 'ASP', 'PRO', 'HIE', 'SER', 'GLU', 'THR', 'MET', 'ALA', 'ILE', 'VAL', 'PHE', 'ASN', 'GLN'}, ':7': {'CYS', 'HID', 'TRP', 'LYS', 'HIP', 'TYR', 'GLY', 'LYN', 'ARG', 'LEU', 'PRO', 'ASP', 'ASH', 'GLH', 'HIE', 'SER', 'GLU', 'THR', 'MET', '
- timed out after 18000 s
