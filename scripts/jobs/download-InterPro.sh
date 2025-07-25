#!/bin/bash
# Job name:
#SBATCH --job-name=interpro_download
#
# Account:
#SBATCH --account=ac_stallerlab
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
## Command(s) to run:

python3 script-InterPro-parallelized.py IPR013087 IPR001356 IPR011598 IPR004827 IPR001766 IPR000536 IPR009071 IPR000418 IPR046360 IPR017956 IPR000327 IPR024752 IPR006612 IPR006600 IPR036388 IPR003656 IPR000679 IPR011539 IPR001132 IPR000770 IPR002857 IPR001346 IPR001739 IPR003150 IPR000232 IPR001275 IPR003350 IPR013801 IPR001606 IPR057520 IPR002100 IPR011129 IPR013854 IPR000818 IPR003690 IPR032200 IPR001523 IPR004212 IPR041686 IPR006578 IPR012295 IPR011615 IPR013524 IPR015351 IPR007889 IPR023082 Q12986 IPR024061 IPR003902 IPR005559 IPR033467 IPR018586 IPR003958 IPR007588