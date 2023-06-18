#!/bin/bash -l
#SBATCH --time=20:00:00
#SBATCH --ntasks=9
#SBATCH --mem=2g
#SBATCH --mail-type=FAIL
#SBATCH -A katagirf
cd /home/katagirf/liux4215/Documents/fit_MCM_Rpt2
module load R/3.6.3
Rscript --vanilla fit_mcm_Rpt2_morecompart.R 1 101 log_normal_morecompart1.Rdata fitting_morecompart1.Rdata
wait
