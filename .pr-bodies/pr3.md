Depends on #__DEP2__ (merge that one first; until then this diff also contains its commits).

**Summary**

Updates the package defaults to the combination recommended for mixed-solvent MD with current Amber. Only the defaults change: the previous combination (ff14SB with TIP3P) remains fully supported.

**Changes**

| Setting | Before | After | Reason |
|---|---|---|---|
| Water model | TIP3P | OPC | Designed together with ff19SB; better description of water-cosolvent mixtures |
| Protein force field | `leaprc.ff14SB` | `leaprc.protein.ff19SB` | Current force field. `leaprc.ff14SB` is not found by tleap in AmberTools 26; `leaprc.protein.ff14SB` is the current name |
| Small molecules | `leaprc.gaff` | `leaprc.gaff2` | Current general force field |
| Solvation buffer | 13 A and 14 A | 10 A | Smaller systems; adequate with a 9 A cutoff |
| Temperature | 300 K | 310 K | Physiological temperature |
| Minimization | `pmemd.cuda` | `mpirun -np 10 pmemd.MPI` | `pmemd.cuda` tends to blow up in the minimization of a freshly solvated system; the CPU code is more robust. Configurable with `AMBER_MIN_EXE` |

**Compatibility**

Nothing is removed. To keep the previous behavior a project sets `SOLVENTS` with TIP3P boxes, `EXTRAFF` with `leaprc.protein.ff14SB, leaprc.gaff, leaprc.water.tip3p`, and the previous buffer and temperature. Legacy boxes must declare `watermodel = TIP3P` in their `.config`.

**Validation**

- 1YER with ETAO and the new defaults: min, eq and 20 ns of production with 4 fs and HMR on three replicas.
- Legacy project (see PR 2) created with the combination above: ff14SB, gaff and TIP3P are loaded after the new defaults, OPC is not reloaded after solvation, temperature of the production input is 300 K. The dynamics of the legacy system was not run.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_016M7mpPRFVJcSvhJWhCpc6p
