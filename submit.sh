#! /bin/sh

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-0:15:0
#SBATCH -p shared
#SBATCH --mem=0
#SBATCH --mail-user=hkim171@jhu.edu
#SBATCH --mail-type=end

module load gcc/5.5.0
module load R/3.6.1

ml

Rscript MARCC_test_script.R

