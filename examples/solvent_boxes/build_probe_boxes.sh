#!/bin/bash
# =============================================================================
# Builds in series the OPC boxes of the probes parametrized with run_probes.sh.
# Concentrations (rationale in TUTORIAL_CAJAS.md, section 9.6):
#   BNZ and IBN 0.25 M: apolar probes with very low solubility in water (benzene ~0.02 M), a 20 % v/v
#                       box would phase-separate; 0.25 M is a supersaturated compromise for sampling.
#   IMD and ACD 0.5 M:  polar probes.
#   GDM and MPH 0.5 M:  ionic probes, no pure liquid exists so % v/v is undefined. A 20 % v/v equivalent
#                       (~2.5-3 M salt) would give extreme ionic strength and ion pairing; 0.5 M gives ~19
#                       molecules in 40 A plus the same number of counterions (Cl- for GDM, Na+ for MPH, OPC ion
#                       parameters from leaprc.water.opc), which keeps the box neutral.
# 40 A boxes (numbers from: python boxcalc.py dilute COSOLV --molar M --side 40; waters for ionic
# probes are slightly reduced to account for the counterion volume, the NPT stage relaxes the box).
#
# Usage (ragnarok, pymdmix3 environment and AmberTools active):
#   cd /mnt/data_ubuntu/new_boxes && mkdir -p probe_boxes && cd probe_boxes
#   nohup bash /path/to/build_probe_boxes.sh > build_probe_boxes.log 2>&1 &
# Variables: PROBES_DIR (run_probes.sh folder), PROBES="BNZ IBN" for a subset.
# =============================================================================
HERE=$(dirname "$(readlink -f "$0")")
B=$HERE/build_box.sh
PROBES_DIR=${PROBES_DIR:-/mnt/data_ubuntu/new_boxes/probes}

# NAME|COMP|NW0|SIDE|SOLVENT|RDF_MASK|IONS
read -r -d '' TABLE << 'TAB'
BNZ|BNZ:10|2091|40|BNZO|:BNZ@C1|
IBN|IBN:10|2082|40|IBNO|:IBN@C1|
IMD|IMD:19|2080|40|IMDO|:IMD@N2|
ACD|ACD:19|2079|40|ACDO|:ACD@O1|
GDM|GDM:19|2070|40|GDMO|:GDM@C1|Cl-:19
MPH|MPH:19|2065|40|MPHO|:MPH@P1|Na+:19
TAB

PROBES=${PROBES:-$(echo "$TABLE" | cut -d'|' -f1 | tr '\n' ' ')}
SUMMARY=""
for p in $PROBES; do
  line=$(echo "$TABLE" | grep "^$p|") || { SUMMARY="$SUMMARY\n$p NO-DEFINIDA"; continue; }
  IFS='|' read -r name comp nw0 side solv rdf ions <<< "$line"
  [ -f "$PROBES_DIR/$name/$name.off" ] || { SUMMARY="$SUMMARY\n$name FALTA $PROBES_DIR/$name/$name.off"; continue; }
  echo "################ $name  ($comp + $nw0 WAT, caja $side A)"
  mkdir -p "$name"
  if ( cd "$name" && COSOLV=$name SRC_OFF="$PROBES_DIR/$name/$name.off" PARMS=${name}parms \
       COMP="$comp" IONS="$ions" NW0=$nw0 SIDE_EQ=$side REP=1 SOLV=$solv UNIT=${name}OPC \
       CFG_SRC="$HERE/probes/$name.config" RDFA="$rdf" RDFB="$rdf" bash "$B" "$name" > build.log 2>&1 ); then
    SUMMARY="$SUMMARY\n$name OK"
  else
    SUMMARY="$SUMMARY\n$name FALLO (ver $PWD/$name/build.log)"
  fi
done
echo
echo "================ RESUMEN"
echo -e "$SUMMARY"
