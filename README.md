1. Compilation under OS Linux:

UNN 19:05:55 makovkin_s@master ~$ lsb_release -a
LSB Version:    :core-4.1-amd64:core-4.1-noarch:cxx-4.1-amd64:cxx-4.1-noarch:desktop-4.1-amd64:desktop-4.1-noarch:languages-4.1-amd64:languages-4.1-noarch:printing-4.1-amd64:printing-4.1-noarch
Distributor ID: CentOS
Description:    CentOS Linux release 7.2.1511 (Core)
Release:        7.2.1511
Codename:       Core

icc -USE_OMP -qopenmp HH_Recogn_2023_03_15_4x5.cpp -o HH_Recogn_2023_03_15_4x5.out

2. Run the script with arguments (files Patterns.txt and SN_input_array.txt need to be exist in the running directory):
/run_script_single.sh

#!/bin/bash
START=$(date +%s)

./HH_Recogn_2023_03_15_4x5.out 1 1 0.17 0.30 0.2496 0.00 1.05 0.60 0.60 0.5040 0  0  0  0
#                                  g_syn_RN              I_app_RN              E_syn_RN
#                                       g_syn_SN              I_app_SN            E_syn_SN
#                                            g_syn_IN              I_app_IN          E_syn_IN
#                                                   g_syn_CN            I_app_CN        E_syn_CN
matlab -nodisplay -nosplash -nodesktop -r "run('./Patterns_Viewer.m'); exit;"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "$(($DIFF / 60)) min $(($DIFF % 60)) sec"

3. Draw results as images:
matlab -nodisplay -nosplash -nodesktop -r "run('./Plot_I_app_CN_g_syn_IN_Compare.m'); exit;"
