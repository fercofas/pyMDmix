Depends on #__DEP3__ (merge that one first; until then this diff also contains its commits).

**Summary**

Tools to build, validate and register mixed-solvent boxes with OPC water, and the OPC versions of the solvents of the library.

**Contents**

- `examples/solvent_boxes/build_box.sh`: box builder in 8 stages (packmol, tleap, heating, 10 ns NPT, density and drift analysis, representative frame, `.off`, `.config` and aggregation RDF). Profiles for ETA, ANT, MAM, ISO, ISO5, MOH, PYR1 and ION; new cosolvents through variables.
- `boxcalc.py` (replaces `molecule_calculator.ods`) and `rdf_recalc.sh`.
- `opc_boxes/`: ETAO, ANTO, MOHO, ISOO, ISO5O, MAMO, IONO and PYR1O with their `.config` files and `register_opc_boxes.sh`.
- `TUTORIAL_CAJAS.md`.
- Updated `SOLVENTS.db` **[to be added after the boxes are final]**.

**Validation** (10 ns NPT at 300 K; density, drift between halves, RDF tails)

| Box | Composition | Density (g/cm3) | Drift |
|---|---|---|---|
| ETAO | 136 ETA + 1672 WAT | 0.9678 | none |
| ANTO | 136 ANT + 1672 WAT | 0.9565 | < 0.001 |
| MOHO | 448 MOH + 3760 WAT | 0.9682 | < 0.001 |
| ISOO | 112 ISO + 1672 WAT | 0.9696 | < 0.001 |
| ISO5O | 32 ISO + 2480 WAT | 0.9890 | < 0.001 |
| IONO | 96 CN3 + 96 COO + 2520 WAT | 1.0472 | < 0.001 |
| PYR1O | 8 PYR + 3504 WAT | 0.9957 | 0.0000 |

**[NOT VALIDATED]** MAMO and comparison with experimental densities of the mixtures.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_016M7mpPRFVJcSvhJWhCpc6p
