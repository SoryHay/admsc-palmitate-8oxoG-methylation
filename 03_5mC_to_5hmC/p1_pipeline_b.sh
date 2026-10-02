#!/usr/bin/env bash
# Second part of the 5mC/5hmC pipeline.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail; source "$(dirname "$0")/p1_common.sh"; SC=$S/scripts; L=$D/logs
run(){ name=$1; shift; echo "step $name start $(date -Is)" >> "$LOG"; "$@" > "$L/$name.log" 2>&1; ec=$?
       echo "exit=$ec" > "$L/$name.DONE"; echo "step $name end $(date -Is) exit=$ec" >> "$LOG"
       [ $ec -eq 0 ] || { echo "ABORT at $name exit=$ec $(date -Is)" >> "$LOG"; exit $ec; }; }
echo "pipeline_b: step 7 relaunched with --modified-bases 5mC 5hmC (modkit 0.6.4 requires it; first attempt exit=2 was a usage error, no data touched)" >> "$LOG"
run step7_pileup_NT bash $SC/p1_step7_pileup.sh NT
run step7_pileup_PA bash $SC/p1_step7_pileup.sh PA
run step8_delta     python3 $SC/p1_step8_delta.py "$D/NT_cpg_thr07.bed" "$D/PA_cpg_thr07.bed" "$S/results/RESULT_P1_5hmC_2026-09-03.md" "P-1, 10 % subset per arm"
echo "PIPELINE_B COMPLETE $(date -Is)" >> "$LOG"
