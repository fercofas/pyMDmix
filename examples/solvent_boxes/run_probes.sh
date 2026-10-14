#!/bin/bash
# =============================================================================
# Parametrizes the new probes in series with make_probe.sh (ORCA + Multiwfn + Amber).
#
# Usage (on ragnarok, with AmberTools active and python with RDKit or obabel available):
#   mkdir -p /mnt/data_ubuntu/new_boxes/probes && cd /mnt/data_ubuntu/new_boxes/probes
#   nohup bash /path/to/run_probes.sh > run_probes.log 2>&1 &
#
# Each probe ends up in its own NAME/ folder with NAME.off ready for build_box.sh.
# To choose a subset:  PROBES="BNZ IBN" bash run_probes.sh
# A failure in one probe does not stop the others; a summary is printed at the end.
# =============================================================================
MP=$(dirname "$(readlink -f "$0")")/make_probe.sh
[ -f "$MP" ] || { echo "ERROR: no encuentro make_probe.sh junto a este script"; exit 1; }

# NAME  SMILES  CHARGE  MULT
read -r -d '' TABLE << 'TAB'
BNZ|c1ccccc1|0|1
IBN|CC(C)C|0|1
IMD|c1c[nH]cn1|0|1
ACD|CC=O|0|1
GDM|NC(N)=[NH2+]|1|1
MPH|COP(=O)([O-])O|-1|1
TAB

PROBES=${PROBES:-$(echo "$TABLE" | cut -d'|' -f1 | tr '\n' ' ')}
SUMMARY=""
for p in $PROBES; do
  line=$(echo "$TABLE" | grep "^$p|") || { echo "$p: no esta en la tabla"; SUMMARY="$SUMMARY\n$p NO-DEFINIDA"; continue; }
  IFS='|' read -r name smi chg mult <<< "$line"
  echo "################ $name  ($smi, carga $chg)"
  if bash "$MP" "$name" "$smi" "$chg" "$mult"; then
    SUMMARY="$SUMMARY\n$name OK"
  else
    SUMMARY="$SUMMARY\n$name FALLO (ver $PWD/$name)"
  fi
done
echo
echo "================ RESUMEN"
echo -e "$SUMMARY"
