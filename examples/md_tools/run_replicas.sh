#!/bin/bash
# Automatiza la corrida de todas las replicas de un proyecto de pyMDMix con una sola equilibracion.
#
# Uso (con el entorno pymdmix3 y Amber cargados):
#   cd MI_PROYECTO/MD
#   bash /ruta/a/pymdmix/examples/md_tools/run_replicas.sh
#
# 1. Corre COMMANDS.sh completo (minimizacion, eq1 a eq5 y md) en la primera replica.
# 2. Copia su eq/eq5.rst7 a eq/ de las demas replicas (mismo sistema de partida).
# 3. En las demas replicas ejecuta solo la parte md de COMMANDS.sh (sin min ni eq).
#    Las replicas 2..N se lanzan de forma secuencial, o en paralelo repartidas entre GPUs.
#
# Variables:
#   GPUS="0 1"    GPUs disponibles; con varias, las replicas 2..N corren en paralelo por tandas (default "0")
#   FIRST=ETAO_1  carpeta de la replica que equilibra (default: la primera en orden alfabetico)
#   SOLO_PRIMERA=1  corre solo el paso 1 y se detiene
#   CORRER_PRIMERA=0  omite el paso 1 (la primera replica ya corrio)

set -e
cd "${1:-.}"
GPUS=${GPUS:-0}
REPS=$(ls -d */ | tr -d '/' | while read d; do [ -f "$d/COMMANDS.sh" ] && echo "$d"; done)
[ -n "$REPS" ] || { echo "ERROR: no hay carpetas con COMMANDS.sh; ejecuta desde la carpeta MD del proyecto"; exit 1; }
FIRST=${FIRST:-$(echo "$REPS" | head -1)}
OTHERS=$(echo "$REPS" | grep -vx "$FIRST" || true)
echo "Replica con equilibracion: $FIRST | resto: $(echo $OTHERS) | GPUs: $GPUS"

if [ "${CORRER_PRIMERA:-1}" = "1" ]; then
  echo "== $FIRST: minimizacion, equilibracion y produccion"
  (cd "$FIRST" && CUDA_VISIBLE_DEVICES=$(echo $GPUS | awk '{print $1}') sh COMMANDS.sh) || echo "AVISO: COMMANDS.sh de $FIRST termino con codigo distinto de cero"
fi
[ "${SOLO_PRIMERA:-0}" = "1" ] && exit 0
[ -s "$FIRST/eq/eq5.rst7" ] || { echo "ERROR: falta $FIRST/eq/eq5.rst7 (la equilibracion no termino)"; exit 1; }

# Parte md de COMMANDS.sh: desde la linea 'cd ../md' en adelante, ejecutada desde la carpeta de la replica
for r in $OTHERS; do
  cp "$FIRST/eq/eq5.rst7" "$r/eq/"
  sed -n '/^cd \.\.\/md$/,$p' "$r/COMMANDS.sh" | sed '1s#.*#cd md#' > "$r/COMMANDS_md.sh"
done

NG=$(echo $GPUS | wc -w)
i=0
for r in $OTHERS; do
  g=$(echo $GPUS | awk -v k=$(( i % NG + 1 )) '{print $k}')
  echo "== $r: produccion en GPU $g"
  ( cd "$r" && CUDA_VISIBLE_DEVICES=$g sh COMMANDS_md.sh || echo "AVISO: $r termino con codigo distinto de cero" ) &
  i=$(( i + 1 ))
  # tanda completa: espera a que terminen antes de lanzar mas
  [ $(( i % NG )) -eq 0 ] && wait
done
wait
echo "LISTO: todas las replicas terminaron."
