#!/usr/bin/env bash
# Open the six pyMDMix 0.3.1 pull requests in CBDD/pyMDmix, in dependency order.
# Usage: bash open_prs.sh            (creates the pull requests)
#        DRY=1 bash open_prs.sh      (only prints what it would run)
set -euo pipefail
UP=CBDD/pyMDmix
ME=fercofas
cd "$(dirname "$0")"

declare -A NUM
create() {   # n branch dep title draft
  local n=$1 branch=$2 dep=$3 title=$4 draft=$5
  local body; body=$(mktemp)
  if [ -n "$dep" ]; then sed "s/__DEP${dep}__/${NUM[$dep]:-$dep}/" "pr${n}.md" > "$body"; else cp "pr${n}.md" "$body"; fi
  local cmd=(gh pr create --repo "$UP" --base master --head "$ME:$branch" --title "$title" --body-file "$body")
  [ "$draft" = "draft" ] && cmd+=(--draft)
  if [ "${DRY:-0}" = "1" ]; then printf '%q ' "${cmd[@]}"; echo; NUM[$n]="$n"; return; fi
  local url; url=$("${cmd[@]}")
  echo "PR $n: $url"
  NUM[$n]=${url##*/}
}

create 1 pr1-python3          ""  "Port pyMDMix to Python 3 (3.11 to 3.13)" final
create 2 pr2-amber-opc-compat 1   "Support current AmberTools and four-point water models with virtual sites (OPC)" draft
create 3 pr3-new-defaults     2   "Change defaults to OPC water, ff19SB, GAFF2, 10 A buffer, 310 K" draft
create 4 pr4-runtime-fixes    2   "Fix stale lock files (#6), '&>' under dash, replica concurrency and RAM of the COM density" draft
create 5 pr5-box-toolkit      3   "Add scripts to build OPC mixed-solvent boxes and OPC boxes for the library solvents" draft
create 6 pr6-docs-tools       4   "Add tutorial, complete project template, MD helper scripts; version 0.3.1" draft
