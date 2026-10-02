#!/usr/bin/env bash
# 5mC/5hmC pipeline: subset, fast5, pod5, benchmark, basecall, align, pileup, both arms.
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
set -uo pipefail; source "$(dirname "$0")/p1_common.sh"; SC=$S/scripts; L=$D/logs
run(){ name=$1; shift; echo "step $name start $(date -Is)" >> "$LOG"; "$@" > "$L/$name.log" 2>&1; ec=$?
       echo "exit=$ec" > "$L/$name.DONE"; echo "step $name end $(date -Is) exit=$ec" >> "$LOG"
       [ $ec -eq 0 ] || { echo "ABORT at $name exit=$ec $(date -Is)" >> "$LOG"; exit $ec; }; }
run step2_fast5_NT    bash $SC/p1_step2_fast5.sh NT
run step2_fast5_PA    bash $SC/p1_step2_fast5.sh PA
run step3_pod5_NT     bash $SC/p1_step3_pod5.sh NT
run step3_pod5_PA     bash $SC/p1_step3_pod5.sh PA
run step4_bench       bash $SC/p1_step4_bench.sh NT
run step5_basecall_NT bash $SC/p1_step5_basecall.sh NT
run step5_basecall_PA bash $SC/p1_step5_basecall.sh PA
run step6_align_NT    bash $SC/p1_step6_align.sh NT
run step6_align_PA    bash $SC/p1_step6_align.sh PA
run step7_pileup_NT   bash $SC/p1_step7_pileup.sh NT
run step7_pileup_PA   bash $SC/p1_step7_pileup.sh PA
echo "PIPELINE COMPLETE $(date -Is)" >> "$LOG"
