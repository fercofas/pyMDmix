#!/bin/bash
# =============================================================================
# Parametrization of a probe (cosolvent) for pyMDMix: RESP charges with ORCA and Multiwfn
#
# Usage:
#   bash make_probe.sh NAME INPUT [CHARGE] [MULT]
#     NAME    residue name, 3 uppercase letters (e.g. BNZ)
#     INPUT   .xyz file with an initial 3D structure, or a quoted SMILES
#             (the SMILES is converted with RDKit or, if missing, with obabel)
#     CHARGE  net charge (default 0)           MULT  spin multiplicity (default 1)
#
# Output (in the NAME/ folder):
#   NAME.off     NAME unit + new parameters (NAMEparms object), ready for build_box.sh
#   NAME.mol2, NAME.frcmod, NAME.charges, NAME_opt.xyz
#
# Protocol: HF/6-31G* optimization, ESP at the same level, two-stage RESP fit
# (Multiwfn, equivalence by connectivity), GAFF2 atom types (antechamber) and frcmod (parmchk2).
# Single conformer: in flexible molecules (OH, rotors) the charges depend on the conformer;
# see TUTORIAL_CAJAS.md, section "Sondas nuevas".
#
# Variables (defaults are those of ragnarok):
#   ORCA_HOME  ORCA 6.x folder              MULTIWFN_HOME  folder with the Multiwfn executable
#   ATYPE=gaff2 atom type set                  NPROCS=1      cores for ORCA (uses the full path)
#   QM_ONLY=1  stop after computing the charges (to run the quantum part on another machine)
# Resuming: if NAME/NAME_opt.xyz and NAME/NAME.charges already exist, the quantum part is skipped
# and only the Amber part is run.
# =============================================================================
set -euo pipefail

NAME=${1:-}; INPUT=${2:-}; CHG=${3:-0}; MULT=${4:-1}
[ -n "$NAME" ] && [ -n "$INPUT" ] || { sed -n '2,20p' "$0"; exit 1; }
[[ "$NAME" =~ ^[A-Z0-9]{1,4}$ ]] || { echo "ERROR: NOMBRE debe ser 1 a 4 caracteres en mayusculas/digitos"; exit 1; }

ORCA_HOME=${ORCA_HOME:-/mnt/data_ubuntu/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg}
MULTIWFN_HOME=${MULTIWFN_HOME:-/mnt/data_ubuntu/software/Multiwfn/Multiwfn_2026.10.1_bin_Linux}
ATYPE=${ATYPE:-gaff2}
NPROCS=${NPROCS:-1}

[ -x "$ORCA_HOME/orca" ] || { echo "ERROR: no existe $ORCA_HOME/orca (define ORCA_HOME)"; exit 1; }
[ -x "$MULTIWFN_HOME/Multiwfn" ] || { echo "ERROR: no existe $MULTIWFN_HOME/Multiwfn (define MULTIWFN_HOME)"; exit 1; }
for c in antechamber parmchk2 tleap python; do command -v $c > /dev/null || { echo "ERROR: falta '$c' en el PATH (source \$AMBERHOME/amber.sh; conda activate)"; exit 1; }; done
[ -n "${AMBERHOME:-}" ] || { echo "ERROR: AMBERHOME no definido"; exit 1; }
[ -f "$AMBERHOME/dat/leap/cmd/leaprc.$ATYPE" ] || { echo "ERROR: no existe leaprc.$ATYPE"; exit 1; }

export LD_LIBRARY_PATH=$ORCA_HOME:${LD_LIBRARY_PATH:-}
export Multiwfnpath=$MULTIWFN_HOME
export PATH=$MULTIWFN_HOME:$PATH
ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=200M

mkdir -p "$NAME" && cd "$NAME"

# ---------------------------- 1. initial structure ----------------------------
if [ -s "${NAME}_opt.xyz" ] && [ -s "$NAME.charges" ]; then
  echo ">> $NAME: reanudando con ${NAME}_opt.xyz y $NAME.charges ya calculados"
  cp "${NAME}_opt.xyz" opt.xyz
  NAT=$(head -1 opt.xyz | tr -d ' \r')
else
if [ -f "../$INPUT" ]; then
  cp "../$INPUT" start.xyz
elif [ -f "$INPUT" ]; then
  cp "$INPUT" start.xyz
else
  python - "$INPUT" "$NAME" << 'PYEOF' || true
import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    sys.exit(0)
smi, name = sys.argv[1:3]
m = Chem.MolFromSmiles(smi)
if m is None:
    sys.exit('SMILES invalido: ' + smi)
m = Chem.AddHs(m)
if AllChem.EmbedMolecule(m, randomSeed=7) != 0:
    sys.exit('RDKit no pudo generar la estructura 3D')
AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
conf = m.GetConformer()
with open('start.xyz', 'w') as f:
    f.write('%d\n%s\n' % (m.GetNumAtoms(), smi))
    for a in m.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        f.write('%s %.6f %.6f %.6f\n' % (a.GetSymbol(), p.x, p.y, p.z))
# topological symmetry classes (equivalent atoms) to average the charges
g = Chem.RWMol(m)                      # generic graph: single bonds, no charges or aromaticity,
for b in g.GetBonds():                 # so that resonance forms (guanidinium, carboxylate,
    b.SetBondType(Chem.BondType.SINGLE)   # phosphate) come out equivalent
    b.SetIsAromatic(False)
for a in g.GetAtoms():
    a.SetFormalCharge(0)
    a.SetIsAromatic(False)
    a.SetNoImplicit(True)
g.UpdatePropertyCache(strict=False)
ranks = list(Chem.CanonicalRankAtoms(g, breakTies=False))
open(name + '.classes', 'w').write(' '.join(str(r) for r in ranks) + '\n')
PYEOF
  if [ ! -s start.xyz ]; then
    command -v obabel > /dev/null || { echo "ERROR: '$INPUT' no es un archivo; no hay RDKit ni obabel para interpretarlo como SMILES"; exit 1; }
    obabel -:"$INPUT" -oxyz --gen3d -O start.xyz > /dev/null 2>&1
  fi
  [ -s start.xyz ] || { echo "ERROR: no se genero la estructura para el SMILES '$INPUT'"; exit 1; }
fi
NAT=$(head -1 start.xyz | tr -d ' \r')
echo ">> $NAME | $NAT atomos | carga $CHG | multiplicidad $MULT"

# ---------------------------- 2. ORCA: HF/6-31G* optimization -------------------
echo "== Optimizacion HF/6-31G* (ORCA)"
{
  echo "! HF 6-31G* Opt TightSCF"
  [ "$NPROCS" -gt 1 ] && echo "%pal nprocs $NPROCS end"
  echo "* xyzfile $CHG $MULT start.xyz"
} > opt.inp
"$ORCA_HOME/orca" opt.inp > opt.out 2>&1 || true
grep -q "HURRAY" opt.out || { echo "ERROR: la optimizacion no convergio (ver $NAME/opt.out)"; exit 1; }
[ -f opt.xyz ] || { echo "ERROR: ORCA no escribio opt.xyz"; exit 1; }
cp opt.xyz "${NAME}_opt.xyz"
EOPT=$(grep "FINAL SINGLE POINT ENERGY" opt.out | tail -1 | awk '{print $5}')
echo "   energia final: $EOPT Eh"

# ---------------------------- 3. RESP with Multiwfn ---------------------------
echo "== Ajuste RESP de dos etapas (Multiwfn)"
"$ORCA_HOME/orca_2mkl" opt -molden > /dev/null 2>&1
[ -f opt.molden.input ] || { echo "ERROR: orca_2mkl no genero opt.molden.input"; exit 1; }
printf '7\n18\n1\ny\n0\nq\n' | Multiwfn opt.molden.input > resp.log 2>&1 || true
[ -f opt.chg ] || { echo "ERROR: Multiwfn no escribio opt.chg (ver $NAME/resp.log)"; exit 1; }
awk '{print $NF}' opt.chg > "$NAME.charges"
[ "$(wc -l < "$NAME.charges")" -eq "$NAT" ] || { echo "ERROR: numero de cargas distinto del de atomos"; exit 1; }
grep -E "Sum of charges|RMSE" resp.log | tail -2
# topologically equivalent atoms (SMILES input only): their charges are averaged
if [ -s "$NAME.classes" ]; then
python - "$NAME.charges" "$NAME.classes" << 'PYEOF'
import sys
from collections import defaultdict
q = [float(x) for x in open(sys.argv[1]).read().split()]
c = open(sys.argv[2]).read().split()
if len(c) == len(q):
    g = defaultdict(list)
    for i, k in enumerate(c):
        g[k].append(i)
    for idx in g.values():
        m = sum(q[i] for i in idx) / len(idx)
        for i in idx:
            q[i] = m
    open(sys.argv[1], 'w').write('\n'.join('%.6f' % x for x in q) + '\n')
    print('   cargas promediadas sobre atomos equivalentes; suma = %.6f' % sum(q))
PYEOF
fi
rm -f opt.*.tmp* opt_trj.xyz opt.bas* opt.densities* opt.hostnames opt.carthess opt.opt opt.property.txt
fi

if [ "${QM_ONLY:-0}" = "1" ]; then
  echo "QM_ONLY: cargas en $PWD/$NAME.charges y geometria en $PWD/${NAME}_opt.xyz. Falta la parte de Amber."
  exit 0
fi

# ---------------------------- 4. PDB with atom names, mol2 and frcmod -----------
python - "$NAME" "${NAME}_opt.xyz" "$NAME.charges" << 'PYEOF'
import sys
from collections import Counter
name, xyz, chg = sys.argv[1:4]
lines = open(xyz).read().split('\n')
n = int(lines[0])
atoms = [l.split() for l in lines[2:2 + n]]
q = [float(x) for x in open(chg).read().split()]
cnt, out = Counter(), []
for i, (el, x, y, z) in enumerate(atoms, 1):
    el = el.capitalize()
    cnt[el] += 1
    out.append('HETATM%5d %-4s %-3s A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s' %
               (i, el.upper() + str(cnt[el]), name, float(x), float(y), float(z), el.upper()))
open(name + '.pdb', 'w').write('\n'.join(out) + '\nEND\n')
PYEOF
antechamber -i "$NAME.pdb" -fi pdb -o "$NAME.mol2" -fo mol2 -c rc -cf "$NAME.charges" -at "$ATYPE" -rn "$NAME" -nc "$CHG" -pf y > antechamber.log 2>&1 \
  || { echo "ERROR: antechamber fallo (ver $NAME/antechamber.log)"; tail -5 antechamber.log; exit 1; }
parmchk2 -i "$NAME.mol2" -f mol2 -o "$NAME.frcmod" -s "$ATYPE" > parmchk2.log 2>&1 || { echo "ERROR: parmchk2 fallo"; exit 1; }
if grep -q "ATTN" "$NAME.frcmod"; then
  echo "   AVISO: parmchk2 marco parametros con ATTN en $NAME.frcmod; revisalos antes de usar la sonda"
fi

# ---------------------------- 5. Amber unit (.off) ------------------------------
cat > make_off.leap << EOF
source leaprc.$ATYPE
$NAME = loadMol2 $NAME.mol2
${NAME}parms = loadAmberParams $NAME.frcmod
check $NAME
saveOff { $NAME ${NAME}parms } $NAME.off
quit
EOF
tleap -f make_off.leap > make_off.log 2>&1 || { echo "ERROR: tleap fallo (ver $NAME/make_off.log)"; tail -5 make_off.log; exit 1; }
[ -s "$NAME.off" ] || { echo "ERROR: tleap no escribio $NAME.off"; exit 1; }

echo "== Resultado"
echo "   archivo: $PWD/$NAME.off  (unidades: $(grep -E '^ "' $NAME.off | tr -d '"' | tr '\n' ' '))"
echo "   cargas RESP (nombre, tipo, carga):"
awk '/@<TRIPOS>ATOM/{f=1;next} /@<TRIPOS>BOND/{f=0} f{printf "      %-5s %-4s %9.4f\n",$2,$6,$9}' "$NAME.mol2"
echo
echo "SIGUIENTE (caja OPC):"
echo "  COSOLV=$NAME SRC_OFF=$PWD/$NAME.off PARMS=${NAME}parms COMP=\"$NAME:NN\" NW0=NNN SIDE_EQ=LL \\"
echo "  SOLV=${NAME}O UNIT=${NAME}OPC20 CFG_SRC= RDFA=\":$NAME@C1\" RDFB=\":$NAME@C1\" bash build_box.sh $NAME"
echo "  (NN, NNN y LL salen de: python boxcalc.py mix SONDA --vv X --side L)"
