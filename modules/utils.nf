

process to_fq {
    cpus 2
    maxForks=3
    input:
	path bam
    output:
	path fq
    storeDir params.fq_dir
    script:
    fq="${bam.simpleName}.fq.gz"
    """
       samtools fastq -@2 $bam | bgzip > $fq
    """
}


process map_sr_to_bam {
    cpus 10
    input:
	tuple val(base), path(sr), path(lr_assembly)
    output:
	tuple val(base), tuple(bam), tuple(bai) tuple(lr_assembly)
    script:
	bam="${base}.bam"
        bai="${bam}.bai"
	"""
        hostname
        bwa-mem2 mem -t 10  ${lr_assembly}/base | samtools view -bT ${lr_assembly}/assemble.fa
        """
}

