#!/bin/bash
# Registra en la base de solventes de pyMDMix las cajas OPC validadas de opc_boxes/:
# ETAO, ANTO, MOHO, ISOO, MAMO, ISO5O e IONO.
# Uso (con el entorno pymdmix3 activo):  bash register_opc_boxes.sh
HERE=$(cd "$(dirname "$0")" && pwd)/opc_boxes
command -v mdmix > /dev/null || { echo "ERROR: mdmix no esta en el PATH (conda activate pymdmix3)"; exit 1; }
cd "$HERE" || exit 1
for cfg in ETAO ANTO MOHO ISOO MAMO ISO5O IONO; do
  echo "##### $cfg"
  mdmix create solvents -f $cfg.config || { echo "ERROR registrando $cfg"; exit 1; }
done
echo
mdmix info solvents | grep -E "ETAO|ANTO|MOHO|ISOO|MAMO|ISO5O|IONO"
