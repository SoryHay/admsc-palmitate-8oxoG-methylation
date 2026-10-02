#!/usr/bin/env bash
# Step runner of the 5mC/5hmC pipeline.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -u
S=/data/8oxo_project/AdMSC_paper/sessions/2026-09-03
NAME=${1:?step name}; SCRIPT=${2:?script}
L=$S/p1_5hmC/logs; mkdir -p "$L"; rm -f "$L/$NAME.DONE"
echo "step $NAME start $(date -Is)" >> "$S/p1_5hmC/RUN_LOG.md"
setsid nohup bash -c "bash '$SCRIPT'; ec=\$?; echo exit=\$ec > '$L/$NAME.DONE'; echo \"step $NAME end \$(date -Is) exit=\$ec\" >> '$S/p1_5hmC/RUN_LOG.md'" \
  > "$L/$NAME.log" 2>&1 < /dev/null &
echo "launched $NAME (pid $!) -> $L/$NAME.log ; DONE marker $L/$NAME.DONE"
