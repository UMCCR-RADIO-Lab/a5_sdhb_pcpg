while read -r d;
do 

basedir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/"

echo "
#!/bin/bash
#PBS -P pq08
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=96GB
#PBS -l ncpus=24
#PBS -l software=umccrise
#PBS -l wd
#PBS -l storage=gdata/pq08
#PBS -l jobfs=150GB

# Load umccrise
source /g/data/pq08/software/umccrise/load_umccrise.sh
# Run umccriseon this node using the local reference and 24 cores
cd ${basedir}/${d}
umccrise final --genomes /g/data/pq08/reference/umccrise/genomes/
" > run_umccrise_${d}.qsub

done < <(ls ${basedir} | grep "^E")
