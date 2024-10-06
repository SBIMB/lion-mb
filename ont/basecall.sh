#PBS -N lion
#PBS -q gpu_1
#PBS -l select=1:ncpus=1:ngpus=1
#PBS -P XXXXXXXX
#PBS -l walltime=12:00:00





cd $lustre
dorado basecaller sup  ${sample}/  > output/${sample}.bam
echo "Finished before time-out" > output/${sample}.finish
dorado summary  output/${sample}.bam > output/${sample}.tsv