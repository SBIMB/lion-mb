
nextflow.enable.moduleBinaries = true

include {nano_qc_plot as raw_plot } from "./modules/nano-qc.nf" addParams(qc_dir: "${params.qc_dir}/raw",status:"raw")
include {qc_summary; nano_qc_plot as qc_plot } from "./modules/nano-qc.nf" \
                                    addParams(qc_dir: "${params.qc_dir}/filtered",status:"filtered")

process filter_lr_poor_quality {
    maxForks 6
    cpus 2
    input:
       path(fq)
    output:
      path("${base}.qc.fq.gz"), emit: seqs
      path("lr_qc.txt")
    publishDir "${params.qc_dir}/lr-reads/"
    script:
       base=fq.simpleName
       qc_cmd = "fastplong -i ${fq} -o ${base}.qc.fq.gz --length_required ${params.min_lr_len} --qualified_quality_phred ${params.min_lr_qual} --json ${base}.qc.json --html ${base}.qc.html --disable_adapter_trimming"
       """
        #!/bin/bash
        set -euxo pipefail 
        hostname
        $qc_cmd 
        echo $qc_cmd > lr_qc.txt
        fastplong --version >> lr_qc.txt
       """
}


process sr_qc {
    input:
	tuple val(base), path(fq)
    output:
	tuple val(base), path("*.qc.*")
        path("fastp.txt")
    script:
        fastp_cmd = "fastp  --average_qual ${params.sr_qual} "
	"""
        $fastp_cmd -i  ${fq[0]}  -I ${fq[1]}  \
           -o ${base}_R1.qc.fq.gz -O  ${base}_R2.qc.fq.gz -h ${base}.qc.html -j ${base}.qc.json
        echo $fastp_cmd > fastp.txt
        fastp --version >> fastp.txt
        """
}


workflow {
    lr_ch = Channel.fromPath(params.lr_raw)
    sr_ch = Channel.fromFilePairs(params.sr_raw+"/*{R1,R2}*") \
	{ it -> it.name.replaceAll("_.*","") }
    main:
	if (params.ont_input_type == "bam") {
    	   ont_ch  = to_fq(lr_ch)
        } else {
  	   ont_ch = lr_ch
        }
        raw_plot(lr_ch)
        filter_lr_poor_quality(ont_ch)
        filtered_lr = filter_lr_poor_quality.out.seqs
	qc_plot(filtered_lr)  | toList | qc_summary


}
