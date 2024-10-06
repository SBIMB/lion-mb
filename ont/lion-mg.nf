



include {nano_qc_plot as raw_plot } from "./modules/nano-qc.nf" \
    addParams(qc_dir: "${params.qc_dir}/raw",status:"raw")
include {qc_summary; nano_qc_plot as qc_plot } from "./modules/nano-qc.nf" \
    addParams(qc_dir: "${params.qc_dir}/filtered",status:"filtered")


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

process filter_lr_poor_quality {
    maxForks 6
    cpus 4
    input:
       path(fq)
    output:
      path("${base}.qc.fq.gz"), emit: seqs
      path("lr_qc.txt"), topic: 'pipeline'
    publishDir "${params.qc_dir}/lr-reads/"
    script:
       base=fq.simpleName
       qc_cmd = "chopper -l ${params.min_lr_len} -t 2 -q ${params.min_lr_qual} -i "
       """
        #!/bin/bash
        set -euxo pipefail 
        hostname
        $qc_cmd  $fq | bgzip > ${base}.qc.fq.gz
        echo $qc_cmd > lr_qc.txt
        chopper --version >> lr_qc.txt
       """
}

process flye_assemble {
    cpus 10
    memory '120GB'
    input:
	tuple val(base), path(fqs)
    output:
	tuple val(base), path(outdir)
        path("flye.txt"), topic: 'pipeline'
    script:
      outdir="${base}-nano-assemble"
      fly_cmd = "flye --meta  -g 100m -t 10  --out-dir $outdir --nano-raw"
      """
      #!/bin/bash
      mkdir -p $outdir
      $fly_cmd  $fqs 
      echo $fly_cmd > flye.txt
      flye --version >> flye.txt
      """
}

process sr_qc {
    input:
	tuple val(base), path(fq)
    output:
	tuple val(base), path("*.qc.*")
        path("fastp.txt"), topic ('pipeline')
    script:
        fastp_cmd = "fastp  --average_qual ${params.sr_qual} "
	"""
        $fastp_cmd -i  ${fq[0]}  -I ${fq[1]}  \
           -o ${base}_R1.qc.fq.gz -O  ${base}_R2.qc.fq.gz -h ${base}.qc.html -j ${base}.qc.json
        echo $fastp_cmd > fastp.txt
        fastp --version >> fastp.txt
        """
}

process assembly_index {
    input:
	tuple val(basec), path(input_dir)
    output:
        tuple val(base), path(input_dir)
    script:
	m=input_dir =~ /(.*)-nano-assemble/
	base=m[0][1]
        """
	bwa-mem2 index -p ${input_dir}/$base ${input_dir}/assembly.fasta
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


process unicycler {
    input:
	tuple val(base), path(lr), path(sr1), path(sr2)
    output:
    file("${base}.txt")
    script:
	 """
         touch ${base}.txt
         """
}



def flat_list (it) {
    return 
}

workflow {
    lr_ch = Channel.fromPath(params.bams)
    sr_ch = Channel.fromFilePairs(params.sr_dir+"/*{R1,R2}*") \
	{ it -> it.name.replaceAll("_.*","") }
    main:
	to_fq(lr_ch) | 	(raw_plot & filter_lr_poor_quality )
        filter_lr_poor_quality.out.seqs | qc_plot | toList | qc_summary
        // combined long and short reads on base name using "join"
        //  -- probably need to put SRs through QC
        joined = filter_lr_poor_quality.out.seqs.map { it -> [it.simpleName, it]}.join(sr_ch)
        // now flatten the lists from [base, lr, [sr0, sr1]] to [base, lr, sr0, sr1]
        unicycler(joined.map { [it[0], it[1], it[2][0], it[2][1]] })
}
