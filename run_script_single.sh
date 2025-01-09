#!/bin/bash
START=$(date +%s)

./HH_Recogn_2023_03_15_4x5.out 1 1 0.17 0.30 0.2496 0.00 1.05 0.60 0.60 0.5040 0  0  0  0
#                                  g_syn_RN              I_app_RN             E_syn_RN
#                                       g_syn_SN              I_app_SN           E_syn_SN
#                                            g_syn_IN              I_app_IN         E_syn_IN
#                                                   g_syn_CN            I_app_CN       E_syn_CN
matlab -nodisplay -nosplash -nodesktop -r "run('./Patterns_Viewer.m'); exit;"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "$(($DIFF / 60)) min $(($DIFF % 60)) sec"
