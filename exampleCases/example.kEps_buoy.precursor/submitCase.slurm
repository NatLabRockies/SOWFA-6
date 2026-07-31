#!/bin/bash
#SBATCH --job-name="prec_NtoS"
#SBATCH --output="case.log"
#SBATCH --partition=abc
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=3960MB
##SBATCH --account=your_hpc_billing_account
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=abc@efg.com


source ~/OpenFOAM/etc/OF-6-env
OpenFOAM-6-env

solver=superDeliciousVanilla
#cores=$NTASKS

blockMesh > log.blockMesh 2>&1


# Decompose large cases
#if [ ! -d "processor0" ]; then
#        echo "*=* No decomposition found. Running decomposePar. *=*"
#        decomposePar > log.decomposePar 2>&1
#else
#        echo "*=* Decomposition already exists. Skipping decomposePar. *=*"
#fi

#running in parallel
#srun -n 128 setFieldsABL -parallel > log.setFieldsABL 2>&1
#srun -n 128 $solver -parallel > log.$solver 2>&1


#running in series on a core
srun setFieldsABL > log.setFieldsABL 2>&1
srun $solver > log.$solver 2>&1
