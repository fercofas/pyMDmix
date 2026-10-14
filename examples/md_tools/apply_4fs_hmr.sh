#!/bin/bash
# Opcional: pasa un proyecto de pyMDMix a pasos de 4 fs con reparticion de masa del hidrogeno (HMR).
#
# Uso (con el entorno pymdmix3 activo, que incluye parmed):
#   cd MI_PROYECTO/MD
#   bash /ruta/a/pymdmix/examples/md_tools/apply_4fs_hmr.sh
#
# 1. Edita todos los *.in de */*/ (min, eq y md de cada replica):
#      dt 0.002 -> 0.004;  ntwx y ntpr 500 -> 250;  nstlim 500000 -> 250000 (md.in y eq5.in);
#      nstlim 100000 -> 50000 (eq1 a eq4);  iwrap 1 -> 0 en eq1 a eq5 y md.
#    Cada ps de simulacion produce el mismo numero de fotogramas y cada archivo md conserva 1 ns.
# 2. Reparte la masa del hidrogeno (parmed HMassRepartition, 3.024 Da por defecto) en la topologia de la
#    primera replica y copia ese resultado a las demas, con el nombre de topologia de cada una.
#    Las replicas de pyMDMix parten del mismo sistema solvatado (mismo .rst7), por eso basta una.
#    Se guarda una copia *.parm7.sin_hmr de cada topologia original.
# Es idempotente: si se vuelve a ejecutar no cambia nada (marca HMR_4FS_HECHO).
# Variables: HMASS=3.024 (masa del hidrogeno tras la reparticion).

set -e
cd "${1:-.}"
HMASS=${HMASS:-3.024}
command -v python > /dev/null || { echo "ERROR: python no esta en el PATH (conda activate pymdmix3)"; exit 1; }
python -c "import parmed" 2> /dev/null || { echo "ERROR: falta parmed (conda activate pymdmix3)"; exit 1; }
ls */*/*.in > /dev/null 2>&1 || { echo "ERROR: no hay */*/*.in aqui; ejecuta el script desde la carpeta MD del proyecto"; exit 1; }

echo "== Archivos de entrada: dt 4 fs"
for f in */*/*.in; do
  sed -i -e 's/dt=0.002/dt=0.004/' -e 's/ntwx=500/ntwx=250/' -e 's/ntpr=500/ntpr=250/' \
         -e 's/nstlim= 100000/nstlim= 50000/' "$f"
done
for f in */*/md.in */*/eq5.in; do
  [ -f "$f" ] && sed -i 's/nstlim=500000/nstlim=250000/' "$f"
done
for f in */*/eq[1-5].in */*/md.in; do
  [ -f "$f" ] && sed -i 's/iwrap=1/iwrap=0/' "$f"
done
grep -h "dt=" */*/*.in | sort | uniq -c

echo "== Reparticion de masa del hidrogeno (${HMASS} Da)"
if [ -f HMR_4FS_HECHO ]; then
  echo "   ya aplicada (HMR_4FS_HECHO); no se repite"
else
  REPS=$(ls -d */ | tr -d '/' | while read d; do ls "$d"/*.parm7 > /dev/null 2>&1 && echo "$d"; done)
  FIRST=$(echo "$REPS" | head -1)
  TOP1=$(ls "$FIRST"/*.parm7 | head -1)
  python - "$TOP1" "$HMASS" << 'PYEOF'
import sys, shutil, parmed
from parmed.tools import HMassRepartition
top, mass = sys.argv[1], float(sys.argv[2])
shutil.copy(top, top + '.sin_hmr')
p = parmed.load_file(top)
HMassRepartition(p, str(mass)).execute()
p.save(top, overwrite=True)
hs = [a.mass for a in p.atoms if a.atomic_number == 1 and a.residue.name not in ('WAT', 'HOH')]
print('   %s: masa minima de H del soluto/cosolvente = %.3f, maxima = %.3f' % (top, min(hs), max(hs)))
PYEOF
  for r in $REPS; do
    [ "$r" = "$FIRST" ] && continue
    t=$(ls "$r"/*.parm7 | head -1)
    cp "$t" "$t.sin_hmr"
    cp "$TOP1" "$t"
    echo "   $t <- $TOP1"
  done
  touch HMR_4FS_HECHO
fi
echo "LISTO. Revisa un .in (por ejemplo cat */md/md.in | head -20) y lanza con run_replicas.sh."
