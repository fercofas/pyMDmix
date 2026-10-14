#!/bin/bash
# =============================================================================
# Construccion de la caja ETANOL / AGUA OPC para pyMDMix
#   17 etanoles + 209 aguas OPC (misma composicion que ETAWAT20, TIP3P)
#
# Uso (desde la carpeta de trabajo, con pymdmix3 y AMBER activados):
#   bash build_eta_opc.sh            # corre todas las etapas
#   bash build_eta_opc.sh 3          # reanuda desde la etapa 3
#
# Etapas:
#   1  Extraer etanol (coordenadas y parametros) de ETAWAT20.off
#   2  Obtener una molecula de agua OPC (con sitio virtual EPW)
#   3  Empaquetar con packmol y generar parm7/rst7 con tleap
#   4  Minimizar y calentar (pmemd)
#   5  Equilibrar en NPT (GPU)
#   6  Analizar densidad y caja, elegir el cuadro representativo
#   7  Construir el .off final con tleap
# =============================================================================
set -euo pipefail

# ---------------------------- parametros -------------------------------------
# La caja original de MDMix (17 ETA + 209 WAT, ~20 A) es demasiado pequena para el corte de 9 A:
# pmemd aborta con "max pairlist cutoff must be less than unit cell max sphere radius" y pmemd.cuda
# la rechaza. Por eso la caja se replica REP veces por lado (REP=2: 8x moleculas, ~40 A).
# Las fracciones de mezcla no cambian.
REP=${REP:-2}
NETA=${NETA:-$(( 17 * REP * REP * REP ))}   # moleculas de etanol
NWAT=${NWAT:-$(( 209 * REP * REP * REP ))}  # moleculas de agua
SIDE0=${SIDE0:-$(awk "BEGIN{print 21.0*$REP}")}  # lado inicial de la caja cubica (A)
TEMP=${TEMP:-300.0}               # K (valor por defecto de MDMix)
NPT_NS=${NPT_NS:-10}              # duracion de la fase NPT (ns)
DISCARD_NS=${DISCARD_NS:-2}       # ns de NPT que se descartan en el analisis
FF=${FF:-leaprc.protein.ff19SB}   # campo de fuerzas del soluto (ff19SB se diseno para usarse con OPC)
UNIT=${UNIT:-ETAOPC20}            # nombre de la unidad caja en el .off final
SRC_OFF=${SRC_OFF:-/mnt/data_ubuntu/mdmix/pymdmix/pyMDMix/data/solventlib/ETAWAT20.off}
PMEMD_CPU=${PMEMD_CPU:-pmemd}
# Calentamiento y NPT en GPU. Si la caja es pequena (REP=1) use CPU: NPT_CMD=pmemd HEAT_CMD=pmemd
PMEMD_GPU=${PMEMD_GPU:-pmemd.cuda}
HEAT_CMD=${HEAT_CMD:-$PMEMD_GPU}
NPT_CMD=${NPT_CMD:-$PMEMD_GPU}
START=${1:-1}

NSTEPS=$(( NPT_NS * 500000 ))     # dt = 2 fs
NFRAMES=$(( NPT_NS * 100 ))       # un cuadro cada 10 ps
DISCARD_FR=$(( DISCARD_NS * 100 ))

# ---------------------------- verificaciones ---------------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: no se encuentra '$1' en el PATH"; exit 1; }; }
need tleap; need cpptraj; need packmol; need python
command -v "$PMEMD_CPU" >/dev/null 2>&1 || PMEMD_CPU=sander
need "$PMEMD_CPU"; need "${HEAT_CMD%% *}"; need "${NPT_CMD%% *}"
[ -f "$SRC_OFF" ] || { echo "ERROR: no existe SRC_OFF=$SRC_OFF"; exit 1; }
python -c "import netCDF4, numpy" || { echo "ERROR: faltan netCDF4/numpy en este entorno"; exit 1; }
[ -n "${AMBERHOME:-}" ] || { echo "ERROR: AMBERHOME no esta definido (source \$AMBERHOME/amber.sh)"; exit 1; }
for f in "$FF" leaprc.water.opc; do
  [ -f "$AMBERHOME/dat/leap/cmd/$f" ] || { echo "ERROR: no existe $AMBERHOME/dat/leap/cmd/$f (ajusta FF=... o revisa tu instalacion)"; echo "Disponibles:"; ls "$AMBERHOME/dat/leap/cmd" | grep -E "^leaprc\.(protein|ff|water)" ; exit 1; }
done
echo ">> Etapa inicial: $START | $NETA ETA + $NWAT WAT | caja inicial $SIDE0 A | NPT ${NPT_NS} ns | UNIT=$UNIT"

run_stage() { [ "$START" -le "$1" ]; }

# ------------------------------ ETAPA 1 --------------------------------------
if run_stage 1; then
echo "== Etapa 1: extraer etanol de $SRC_OFF"
cat > 1_eta.leap << EOF
source $FF
loadOff $SRC_OFF
savePdb ETA eta.pdb
saveOff ETA eta.lib
quit
EOF
tleap -f 1_eta.leap > 1_eta.log || { echo "ERROR: tleap fallo con 1_eta.leap, revisa 1_eta.log"; tail -5 1_eta.log; exit 1; }
[ -s eta.pdb ] && [ -s eta.lib ] || { echo "ERROR etapa 1: revisa 1_eta.log"; exit 1; }
fi

# ------------------------------ ETAPA 2 --------------------------------------
if run_stage 2; then
echo "== Etapa 2: molecula de agua OPC"
cat > 2_wat.leap << EOF
source $FF
source leaprc.water.opc
savePdb WAT wat_opc.pdb
quit
EOF
tleap -f 2_wat.leap > 2_wat.log || { echo "ERROR: tleap fallo con 2_wat.leap, revisa 2_wat.log"; tail -5 2_wat.log; exit 1; }
[ -s wat_opc.pdb ] || { echo "ERROR etapa 2: revisa 2_wat.log"; exit 1; }
WRES=$(awk '/^ATOM|^HETATM/{print substr($0,18,3); exit}' wat_opc.pdb | tr -d ' ')
WATOMS=$(awk '/^ATOM|^HETATM/{printf "%s ", substr($0,13,4)}' wat_opc.pdb)
echo "   residuo del agua: '$WRES' | atomos: $WATOMS"
[ "$WRES" = "WAT" ] || { echo "ERROR: el residuo de agua OPC se llama '$WRES', no 'WAT'. Avisame antes de seguir."; exit 1; }
echo "$WATOMS" | grep -q "EPW" || { echo "ERROR: no aparece el sitio virtual EPW en wat_opc.pdb"; exit 1; }
fi

# ------------------------------ ETAPA 3 --------------------------------------
if run_stage 3; then
echo "== Etapa 3: empaquetar y construir topologia"
cat > mix.inp << EOF
tolerance 2.0
filetype pdb
output mix_packed.pdb
add_amber_ter
seed 12345

structure eta.pdb
  number $NETA
  inside box 0. 0. 0. $SIDE0 $SIDE0 $SIDE0
end structure

structure wat_opc.pdb
  number $NWAT
  inside box 0. 0. 0. $SIDE0 $SIDE0 $SIDE0
end structure
EOF
packmol < mix.inp > 3_packmol.log
grep -q "Success" 3_packmol.log || { echo "ERROR etapa 3: packmol no convergio (ver 3_packmol.log)"; exit 1; }

cat > 3_build.leap << EOF
source $FF
source leaprc.water.opc
loadOff eta.lib
mix = loadPdb mix_packed.pdb
set mix box { $SIDE0 $SIDE0 $SIDE0 }
check mix
charge mix
saveAmberParm mix mix.parm7 mix.rst7
quit
EOF
tleap -f 3_build.leap > 3_build.log || { echo "ERROR: tleap fallo con 3_build.leap, revisa 3_build.log"; tail -5 3_build.log; exit 1; }
[ -s mix.parm7 ] && [ -s mix.rst7 ] || { echo "ERROR etapa 3: revisa 3_build.log"; exit 1; }
echo "   carga total reportada por tleap:"; grep -i "total unperturbed charge\|Total perturbed charge" 3_build.log || true
grep -E "Errors = [1-9]|Fatal Error|FATAL" 3_build.log && { echo "AVISO: tleap reporto errores, revisa 3_build.log"; } || true
fi

# ------------------------------ ETAPA 4 --------------------------------------
if run_stage 4; then
echo "== Etapa 4: minimizacion y calentamiento"
cat > min.in << EOF
minimizacion
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=9.0, ntpr=500,
 /
EOF
cat > heat.in << EOF
calentamiento NVT
 &cntrl
  imin=0, irest=0, ntx=1,
  nstlim=25000, dt=0.002,
  ntc=2, ntf=2, cut=9.0,
  ntb=1,
  ntt=3, gamma_ln=2.0, ig=-1,
  tempi=100.0, temp0=$TEMP,
  ntpr=1000, ntwx=0, ntwr=25000,
 /
EOF
$PMEMD_CPU -O -i min.in  -o min.out  -p mix.parm7 -c mix.rst7  -r min.rst7
$HEAT_CMD -O -i heat.in -o heat.out -p mix.parm7 -c min.rst7  -r heat.rst7
fi

# ------------------------------ ETAPA 5 --------------------------------------
if run_stage 5; then
echo "== Etapa 5: NPT de ${NPT_NS} ns con: $NPT_CMD"
cat > npt.in << EOF
NPT
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=$NSTEPS, dt=0.002,
  ntc=2, ntf=2, cut=9.0,
  ntb=2, ntp=1, pres0=1.0, barostat=2,
  ntt=3, gamma_ln=2.0, temp0=$TEMP, ig=-1,
  iwrap=1,
  ntpr=5000, ntwx=5000, ntwr=500000, ioutfm=1,
 /
EOF
$NPT_CMD -O -i npt.in -o npt.out -p mix.parm7 -c heat.rst7 -r npt.rst7 -x npt.nc
fi

# ------------------------------ ETAPA 6 --------------------------------------
if run_stage 6; then
echo "== Etapa 6: analisis de densidad y caja"
cat > analyze_density.py << 'EOF'
import sys
import numpy as np
from netCDF4 import Dataset

parm, nc, discard = sys.argv[1], sys.argv[2], int(sys.argv[3])

# masa total (uma) desde la seccion MASS del parm7 (EPW tiene masa 0)
lines = open(parm).read().split('\n')
i = next(k for k, l in enumerate(lines) if l.startswith('%FLAG MASS'))
masses = []
for l in lines[i + 2:]:
    if l.startswith('%FLAG'):
        break
    masses += [float(x) for x in l.split()]
mass = sum(masses)

ds = Dataset(nc)
cell = np.array(ds.variables['cell_lengths'][:], dtype=float)
angles = np.array(ds.variables['cell_angles'][:], dtype=float)
nf = len(cell)
if nf <= discard + 10:
    sys.exit('ERROR: muy pocos cuadros (%d) para descartar %d' % (nf, discard))
if np.abs(angles - 90.0).max() > 1e-3:
    sys.exit('ERROR: la caja no es ortogonal, este script asume caja ortogonal')

vol = cell.prod(axis=1)
rho = mass * 1.66053907 / vol          # g/cm3
sl = slice(discard, nf)
v, r, c = vol[sl], rho[sl], cell[sl]
half = len(r) // 2
print('cuadros totales: %d | usados: %d' % (nf, len(r)))
print('masa total (uma): %.3f' % mass)
print('densidad media: %.4f +/- %.4f g/cm3' % (r.mean(), r.std()))
print('  primera mitad: %.4f | segunda mitad: %.4f (deriva: %.4f)' %
      (r[:half].mean(), r[half:].mean(), r[half:].mean() - r[:half].mean()))
print('lado medio de la caja (A): %.3f %.3f %.3f' % tuple(c.mean(axis=0)))
print('volumen medio (A^3): %.1f' % v.mean())
best = int(np.argmin(np.abs(v - v.mean())))
idx = discard + best
print('cuadro representativo (1-based): %d | caja: %.4f %.4f %.4f | densidad %.4f' %
      (idx + 1, cell[idx][0], cell[idx][1], cell[idx][2], rho[idx]))
open('best_frame.txt', 'w').write('%d %.4f %.4f %.4f\n' % (idx + 1, cell[idx][0], cell[idx][1], cell[idx][2]))
np.savetxt('density_vs_frame.dat', np.c_[np.arange(1, nf + 1), rho, cell], fmt='%d %.5f %.4f %.4f %.4f',
           header='frame density(g/cm3) a b c')
EOF
python analyze_density.py mix.parm7 npt.nc "$DISCARD_FR" | tee 6_analysis.log
fi

# ------------------------------ ETAPA 7 --------------------------------------
if run_stage 7; then
echo "== Etapa 7: construir el .off final ($UNIT)"
read FRAME BA BB BC < best_frame.txt
cpptraj -p mix.parm7 > 7_cpptraj.log << EOF
trajin npt.nc $FRAME $FRAME
trajout best.pdb pdb include_ep
EOF
[ -s best.pdb ] || { echo "ERROR etapa 7: revisa 7_cpptraj.log"; exit 1; }
NEP=$(grep -c " EPW" best.pdb || true)
echo "   sitios virtuales EPW en best.pdb: $NEP (esperados: $NWAT)"
[ "$NEP" -eq "$NWAT" ] || { echo "ERROR: el numero de EPW no coincide"; exit 1; }

cat > 7_off.leap << EOF
source $FF
source leaprc.water.opc
loadOff eta.lib
$UNIT = loadPdb best.pdb
set $UNIT box { $BA $BB $BC }
check $UNIT
saveOff { ETA WAT $UNIT } $UNIT.off
quit
EOF
tleap -f 7_off.leap > 7_off.log || { echo "ERROR: tleap fallo con 7_off.leap, revisa 7_off.log"; tail -5 7_off.log; exit 1; }
[ -s $UNIT.off ] || { echo "ERROR etapa 7: tleap no escribio $UNIT.off (revisa 7_off.log; si rechazo la lista de unidades, avisame)"; exit 1; }
echo "   unidades guardadas en $UNIT.off:"
sed -n '2,5p' $UNIT.off
echo
echo "LISTO: $UNIT.off. Caja: $BA x $BB x $BC A."
fi
