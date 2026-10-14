# Building mixed-solvent boxes with OPC water, step by step

This document describes how a mixed-solvent box (organic cosolvent in OPC water) is built, validated and registered in pyMDMix. It applies to the cosolvents of the library (ETA, ANT, MAM, ISO, ISO5, MOH, PYR1 and ION) and to any new cosolvent whose parameters are already stored in an `.off` file. The files are in `examples/solvent_boxes/`.

> Validation status. The ethanol box (ETAO) was built and validated on a real installation: density 0.9678 ± 0.0052 g/cm³ (approximate experimental value for 17 % w/w at 300 K: 0.97 g/cm³), no drift between the halves of the trajectory and 1672 EPW virtual sites. The other cosolvent profiles follow the same procedure and are validated one by one. Every box must pass the validation of section 5 before it is used.

## 1. Rationale

1. **Cosolvents are not reparametrized.** The RESP charges and atom types of the cosolvent were derived in the gas phase or with GAFF and do not depend on the water model. What changes is the solvent-water balance, which is why the box is rebuilt, re-equilibrated and validated (density and aggregation).
2. **OPC is a four-point model.** The water molecule includes a virtual site, `EPW`. Therefore `leaprc.water.opc` must be loaded, PDB files must be written with `include_ep`, and the number of `EPW` sites must be checked against the number of waters.
3. **ff19SB and OPC were designed together.** The script loads `leaprc.protein.ff19SB` and `leaprc.gaff2` (organic cosolvents use GAFF2 atom types).
4. **Box size.** The original MDMix boxes measure about 20 Å. With a 9 Å cutoff, `pmemd` requires the radius of the largest sphere that fits in the cell to be larger than cutoff + 2 Å, and `pmemd.cuda` rejects small boxes ("Small box detected"). The script replicates the original box `REP` times per side (by default, the smallest integer that gives at least 36 Å). The mixture fraction does not change.

## 2. Requirements

```bash
conda activate pymdmix3
source $AMBERHOME/amber.sh
which tleap cpptraj packmol pmemd pmemd.cuda
ls $AMBERHOME/dat/leap/cmd | grep -E "^leaprc\.(protein.ff19SB|gaff2|water.opc)$"
```

## 3. Choosing the composition with `boxcalc.py`

`boxcalc.py` replaces the `molecule_calculator.ods` spreadsheet. It computes the number of molecules from the volumes of the pure components (V = M / (ρ N_A)). It does not consider the excess volume of the mixture, which is why the final density is obtained from the NPT run.

```bash
cd examples/solvent_boxes
python boxcalc.py table                                   # available cosolvents
python boxcalc.py library                                 # actual concentration of the bundled TIP3P boxes
python boxcalc.py mix ETA --vv 20 --side 20               # molecules for 20 % v/v in 20 Å
python boxcalc.py check ETA --nco 136 --nwat 1672 --side 39.67
python boxcalc.py pure ETA --side 23                      # pure cosolvent box
python boxcalc.py molar UREA --molarity 7 --volume 0.1    # solution of a given molarity
```

Actual concentration of the TIP3P boxes of the library, computed with `boxcalc.py library`:

| Box | Cosolvent + water | % v/v | % w/w | Mole fraction | Molarity |
|---|---|---|---|---|---|
| ETAWAT20 | 17 + 209 | 20.9 | 17.2 | 0.0752 | 3.5 M |
| ANTWAT20 | 17 + 209 | 19.1 | 15.6 | 0.0752 | 3.6 M |
| MAMWAT20 | 18 + 209 | 19.6 | 22.0 | 0.0793 | 3.8 M |
| ISOWAT20 | 14 + 209 | 22.1 | 18.3 | 0.0628 | 2.9 M |
| ISOWAT5 | 4 + 310 | 5.2 | 4.1 | 0.0127 | 0.7 M |
| MOHWAT20 | 56 + 470 | 21.1 | 17.5 | 0.1065 | 5.2 M |
| PYRWAT1 | 1 + 438 | 1.0 | 1.0 | 0.0023 | 0.08 M (see note) |

Remarks on the original spreadsheet:

- In the `ethanol_20%` sheet, the molecular volume of ethanol (`ETHOH_vol`, cell B16) is computed with the molar mass of urea (60.06 g/mol, cell B6) instead of that of ethanol (46.07 g/mol). The correct result is 97.0 Å³ per molecule, not 126.4 Å³. With the correct value, 17 ETA + 209 WAT corresponds to 20.9 % v/v, consistent with the name "20".
- The densities of ANT, MAM, ISO, MOH and PYR used by `boxcalc.py` were not in the spreadsheet. They are literature values (marked `lit` in `boxcalc.py table`); confirm them with your supplier before publishing.
- The PYRWAT1 box stored in the library has 438 waters in a side of about 28 Å (starting density close to 0.60 g/mL), that is, it contains a large void. `build_box.sh` does not use that side: for PYR1 it starts from 23.65 Å, which corresponds to the actual composition.

## 4. Building the box

One empty working folder per solvent:

```bash
mkdir -p ~/boxes/ANT && cd ~/boxes/ANT
nohup bash /path/to/pyMDmix/examples/solvent_boxes/build_box.sh ANT > build.log 2>&1 &
tail -f build.log
```

Bundled profiles: `ETA`, `ANT`, `MAM`, `ISO`, `ISO5`, `MOH`, `PYR1`, `ION`. Each profile takes the TIP3P `.off` from `pyMDMix/data/solventlib/`, the number of molecules and the side of the original box.

Stages:

| Stage | What it does | Main output |
|---|---|---|
| 1 | Extracts the cosolvent units from the TIP3P `.off` | `RES.pdb`, `RES.lib` |
| 2 | Obtains an OPC water molecule and checks the `EPW` atom | `wat_opc.pdb` |
| 3 | Packs with packmol and builds the topology with tleap | `mix.parm7`, `mix.rst7` |
| 4 | Minimizes and heats (50 ps) | `heat.rst7` |
| 5 | 10 ns NPT on GPU | `npt.nc`, `npt.out` |
| 6 | Computes density and box size; selects the representative frame | `6_analysis.log`, `best_frame.txt` |
| 7 | Writes `best.pdb` (with `include_ep`) and the final `.off` | `UNIT.off` |
| 8 | Generates the solvent `.config` and the aggregation RDFs | `SOLV.config`, `rdf_*.dat` |

To resume from a given stage: `bash build_box.sh ANT 5`. Variables that can be overridden from the environment: `REP`, `NPT_NS`, `DISCARD_NS`, `TEMP`, `FFS`, `HEAT_CMD`, `NPT_CMD`, `SRC_OFF`, `SOLV`, `UNIT`, `IONS`.

Examples:

```bash
REP=3 bash build_box.sh MOH            # larger box
NPT_NS=20 bash build_box.sh ISO5       # longer NPT
REP=1 HEAT_CMD=pmemd NPT_CMD=pmemd bash build_box.sh PYR1   # CPU, only if the box is already large
```

Note: the messages printed by the scripts at run time are still in Spanish (for example `LISTO`, `AVISO`, `FALLO`, meaning done, warning and failure).

### 4.1 New cosolvent

If you already have an `.off` with the parametrized cosolvent unit, define the profile through variables:

```bash
python boxcalc.py mix TFE --vv 20 --side 21 --rep 2    # gives NCO and NWAT
COSOLV=TFE SRC_OFF=/path/TFE.off COMP="TFE:NN" NW0=NNN SIDE_EQ=21 \
SOLV=TFEO UNIT=TFEOPC20 CFG_SRC=/path/tfe.config RDFA=":TFE@O1" RDFB=":TFE@O1" \
bash build_box.sh TFE
```

Here `COMP` and `NW0` describe the original, unreplicated box (the script applies `REP`). The `.off` must contain one unit named after each residue in `COMP`. If the cosolvent is not parametrized, charges and parameters must be derived first (antechamber, RESP), which this procedure does not do (see section 9).

## 5. Validating the box

Minimum conditions before registering:

1. `6_analysis.log`: the density drift between the first and the second half must be close to zero (|drift| smaller than the standard deviation, typically < 0.002 g/cm³), and the mean density must be within a few percent of the experimental value of the mixture at that weight fraction and temperature.
2. Stage 7 must end with "sitios virtuales EPW en best.pdb" (EPW virtual sites in best.pdb) equal to the expected value, and the `.off` must list the cosolvent units, `WAT` and the box unit.
3. Stage 8 prints the height of the first peak and the tail value of g(r) between 9 and 12 Å (it must be close to 1). A cosolvent-cosolvent g(r) with a high first maximum (> 3) below 5 Å indicates aggregation. Compare with the TIP3P box if in doubt. Boxes with few cosolvent molecules (ISO5, PYR1) give noisy RDFs.
4. Inspect the last trajectory visually (no molecules crossing the box in an anomalous way).

To obtain the experimental reference density of an ethanol-water mixture, for example, use the density tables by % w/w and temperature (see the "% w/w" column of `boxcalc.py check`).

## 6. Registering the solvent in pyMDMix

```bash
mkdir ~/mdmix_boxes && cd ~/mdmix_boxes
cp ~/boxes/ANT/ANTOPC20.off ~/boxes/ANT/ANTO.config .
grep -n "name\|objectfile\|boxunit\|watermodel" ANTO.config    # check
mdmix create solvents -f ANTO.config
mdmix info solvents | grep ANTO
```

The `.config` generated by stage 8 reuses the `[PROBES]` and `[TYPES]` sections of the original `.config`. If the command fails with `InvalidPath`, the message states which file was not found (normally `objectfile`).

Solvent names used by this construction: ETAO, ANTO, MAMO, ISOO, ISO5O, MOHO, PYR1O and IONO. They keep the TIP3P boxes (ETA, ANT, ...) separate from the OPC ones. The validated boxes are stored in `examples/solvent_boxes/opc_boxes/` and can be registered all at once with `register_opc_boxes.sh`.

## 7. Testing the box with a system

```ini
[SYSTEM]
NAME = PEP
OFF = pep.off
# With OPC water its force field must be loaded explicitly (or rely on the package default)
EXTRAFF = leaprc.water.opc

[MDSETTINGS]
SOLVENTS = ETAO
NREPL = 1
NANOS = 1
RESTR = FREE
```

```bash
mdmix create project -n PEP_OPC -f pep_opc.cfg
```

In the created replica, `tleap` must load `leaprc.water.opc` and the topology must keep the `EPW` sites.

## 8. Common problems

| Symptom | Cause and solution |
|---|---|
| `max pairlist cutoff must be less than unit cell max sphere radius` or `Small box detected` | The box is too small. Increase `REP`. |
| tleap cannot find an atom type or parameter of the cosolvent | The force field of that type is not loaded (for example `leaprc.gaff2`) or an `frcmod` is missing. Adjust `FFS`. |
| `saveOff` rejects the unit list | Check `7_off.log`. The `saveOff { ... }` syntax was verified with the ethanol box (ETAOPC20.off contains ETA, WAT and the box unit). |
| The number of `EPW` sites does not match | The PDB was written without `include_ep`, or the water residue is not `WAT`. |
| `tleap` warning with "Errors = N" different from 0 | Check `3_build.log` before continuing. |
| The density does not converge | Lengthen the NPT with `NPT_NS=20` and increase `DISCARD_NS`. |

## 9. New probes: parametrization with ORCA and Multiwfn

This section parametrizes a cosolvent that is not in the library (benzene, isobutane, imidazole, acetaldehyde, guanidinium, methyl phosphate or others). It produces the Amber unit with RESP charges and the missing GAFF parameters, and the box is then built as in section 4.

### 9.1 Required programs

| Program | Use | Source |
|---|---|---|
| ORCA 6.x | HF/6-31G* optimization and ESP | download with an academic account from the ORCA forum; `tar -xJf orca_*.tar.xz` |
| Multiwfn | two-stage RESP fit | download from the Multiwfn website; `unzip` |
| AmberTools | `antechamber`, `parmchk2`, `tleap` | already installed |
| RDKit (or Open Babel) | SMILES to 3D structure | `python -m pip install rdkit` in the `pymdmix3` environment |

Installation notes:

- On Ubuntu `/usr/bin/orca` exists, which is the GNOME screen reader and not the quantum chemistry program. Do not put ORCA in the PATH under that name. Define a variable and use the full path:

```bash
export ORCA_HOME=/path/to/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA_HOME:$LD_LIBRARY_PATH
$ORCA_HOME/orca input.inp > output.out
```

- Multiwfn needs an unlimited stack: `ulimit -s unlimited` and `export KMP_STACKSIZE=200M` (the script does this).
- ORCA check: a water molecule at HF/6-31G* must give -76.0090 Eh and finish with `ORCA TERMINATED NORMALLY`.

### 9.2 A single probe

```bash
conda activate pymdmix3
source $AMBERHOME/amber.sh
python -m pip install rdkit            # once
bash examples/solvent_boxes/make_probe.sh BNZ "c1ccccc1"
```

Arguments: residue name (1 to 4 uppercase characters), SMILES or an `.xyz` file, net charge and multiplicity (optional; 0 and 1 by default). The script:

1. generates the 3D structure (RDKit, MMFF);
2. optimizes at HF/6-31G* with ORCA and computes the wavefunction;
3. performs the two-stage RESP fit with Multiwfn and averages the charges of topologically equivalent atoms (resonance forms included);
4. assigns GAFF types with `antechamber` reading those charges, generates the `frcmod` with `parmchk2` and warns if `parmchk2` flagged parameters with `ATTN`;
5. saves `NAME.off` with the unit and the new parameters (object `NAMEparms`).

With `QM_ONLY=1` it stops after the charges, and if `NAME_opt.xyz` and `NAME.charges` already exist the quantum step is skipped.

### 9.3 Several probes

```bash
mkdir -p ~/boxes/probes && cd ~/boxes/probes
nohup bash /path/to/examples/solvent_boxes/run_probes.sh > run_probes.log 2>&1 &
```

`run_probes.sh` processes in series benzene (BNZ), isobutane (IBN), imidazole (IMD), acetaldehyde (ACD), guanidinium (GDM, charge +1) and methyl phosphate (MPH, monoanion). A subset is chosen with `PROBES="BNZ IBN"`. At the end it prints a summary with `OK` or `FALLO` (failure) per probe.

### 9.4 Validation of the protocol

The protocol was checked with ethanol (HF/6-31G*, trans conformer, ORCA 6.1.1 and Multiwfn 2026.10.1) against the charges of the ETA unit in `ETAWAT20.off`:

| Atom | Library | This protocol |
|---|---|---|
| C1 (CH3) | -0.250 | -0.271 |
| C2 (CH2) | 0.519 | 0.537 |
| O1 | -0.711 | -0.706 |
| H1 (OH) | 0.415 | 0.411 |
| H2, H3 (CH2) | -0.082 | -0.089 |
| H4 to H6 (CH3) | 0.064 | 0.069 |

The differences are 0.02 e at most. The gauche conformer (0.09 kcal/mol higher) gives charges that differ more (up to 0.07 e on the methyl carbon), so the conformer dependence is larger than the difference with the library.

### 9.5 Limitations

- A single conformer. For molecules with rotors or OH groups it is advisable to fit several conformers with weights; this first version does not do so.
- Gas-phase RESP at HF/6-31G*, which is the Amber convention. For anions (methyl phosphate) and cations the gas phase is an approximation; dianions (phosphate at pH 7) are problematic without diffuse functions and are not included.
- The protonation state is a user decision. Methyl phosphate is parametrized as a monoanion.
- Benzene and isobutane are not miscible with water at 20 %. They are used at low concentration with `boxcalc.py dilute`:

```bash
python boxcalc.py dilute BENZ --molar 0.25 --side 40
```

### 9.6 Concentration of the probe boxes

The 20 % v/v convention used for the boxes of miscible cosolvents (ethanol, acetonitrile, methanol, isopropanol) follows from c = 0.20 · ρ · 1000 / M and cannot be applied to every probe. The criteria used in `build_probe_boxes.sh` are:

| Probe | Concentration | Molecules in 40 Å | Criterion |
|---|---|---|---|
| BNZ (benzene), IBN (isobutane) | 0.25 M | 10 | Apolar probes with very low solubility in water (benzene about 0.02 M); a 20 % v/v box would separate into two phases. 0.25 M is a supersaturated concentration that allows sampling without massive aggregation. |
| IMD (imidazole), ACD (acetaldehyde) | 0.5 M | 19 | Polar, soluble probes. |
| GDM (guanidinium, +1), MPH (methyl phosphate, -1) | 0.5 M | 19 | Ionic probes: there is no pure liquid of reference, so % v/v is undefined. A 20 % v/v equivalent (about 2.5 to 3 M of salt) would give an extreme ionic strength and ion pairing with the counterion that would distort the RDF. |

Ionic probes are built with the same number of counterions (Cl- for GDM, Na+ for MPH, using the OPC ion parameters of `leaprc.water.opc`), so that the box is neutral. They are specified with the `IONS` variable of `build_box.sh` (for example `IONS="Cl-:19"`). With 19 molecules the sampling is limited (as in the pyridine box, with only 8), so the first peaks of the cosolvent-cosolvent RDFs must be interpreted with caution.

### 9.7 Building the probe box

Use `NAME.off`, the object names and the counts from `boxcalc.py` in the variables of `build_box.sh` (section 4.1). For example, for benzene at 0.25 M in 40 Å:

```bash
mkdir -p ~/boxes/BNZ && cd ~/boxes/BNZ
COSOLV=BNZ SRC_OFF=~/boxes/probes/BNZ/BNZ.off PARMS=BNZparms \
COMP="BNZ:10" NW0=2091 SIDE_EQ=40 REP=1 \
SOLV=BNZO UNIT=BNZOPC025 CFG_SRC= RDFA=":BNZ@C1" RDFB=":BNZ@C1" \
bash /path/to/build_box.sh BNZ
```

For an ionic probe add `IONS="Cl-:19"` (or `Na+:19`) to the variables above. Since there is no source `.config`, the script does not generate `BNZO.config`; it has to be written by hand (section 6 and `examples/solvent_boxes/eta_opc.config`), with probes specific to the molecule (for example, `BNZ=BNZ@C1,C2,C3,C4,C5,C6` and type `Arom`). The ready-made `.config` files of the probes are in `examples/solvent_boxes/probes/`.
