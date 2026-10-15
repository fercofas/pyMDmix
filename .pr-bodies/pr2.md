Depends on #__DEP1__ (merge that one first; until then this diff also contains its commits).

**Summary**

Makes the Amber backend work with AmberTools 22 to 26 and with water models that include virtual sites.

**Changes**

- Topology and restart files use `.parm7` and `.rst7` (min, eq, md, `check_com.sh`, `getSimVolume`).
- `check_com.sh` and the system PDB keep virtual sites (`include_ep`). The system PDB is written with cpptraj so that it has the same atoms as the topology and trajectories.
- Residue and atom numbers are renumbered when `ambpdb` writes hybrid-36 numbering (more than 9999 residues), which Biskit cannot parse.
- Neutralization with NaCl loads the ion parameters of the water model of the solvent box used to solvate (`leaprc.water.tip3p`, `opc`, `opc3`, ...). Before, `neutralizeWNaCl` took the model from the replica, which does not exist when a system is solvated, and fell back to `DEF_AMBER_WATBOX`.
- `Leaper.flush` never ended when tleap exited (`readline()` returns an empty string forever), running at 100 % CPU. It now stops at EOF. `test_tleap` passed a string as `extraff` and used `leaprc.ff99SB`, which no longer exists.
- `InvalidPath` reports which `objectfile` or `frcmod` was not found.

**Compatibility**

Projects created before this change keep working; new projects write `.parm7` and `.rst7`.

**Validation**

- Ethanol/OPC project on 1YER: 14415 WAT with 14415 EPW, net charge 0.0001, topology and PDB consistent.
- Legacy project (TIP3P boxes, `leaprc.protein.ff14SB`, `leaprc.gaff`, `leaprc.water.tip3p`) created with this branch: `leap.log` shows `solvateOct` with the TIP3P box followed by `leaprc.water.tip3p` and `addIons`, with no OPC file loaded after solvation; no EPW in the topology. The dynamics of that system was not run.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_016M7mpPRFVJcSvhJWhCpc6p
