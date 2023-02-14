module load osca

T=$( echo $1 | cut -d"." -f1)
A=$( echo $T | cut -d"/" -f5)
C=$( echo $A | cut -d"_" -f5)

osca \
--vqtl \
--bfile /home/rhillary/imputed/chr22 \
--pheno $1 \
--vqtl-mtd 2 \
--thread-num $SLURM_CPUS_PER_TASK \
--out /scratch/rhillary/residualised_vqtls/chr22_${A}
