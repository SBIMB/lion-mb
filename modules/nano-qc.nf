
params.qc_dir=false

println "Dir is ${params.qc_dir}"

process nano_qc_plot {
    input:
	path(fq)
    output:
        path(qcdir), topic: 'report'	
	if (params.qc_dir) {
	   publishDir params.qc_dir
        }
    script:
        base=fq.simpleName
	qcdir="${base}-qc"
        """
        NanoPlot --fastq $fq --outdir $qcdir -p ${base}-${params.status} 
        """
}

process qc_summary {
    input:
	path f
    output:
	path "lr_qc.tex"
	if (params.qc_dir) {
	   publishDir params.qc_dir
        }
    script:
	"""
        summarise_qc.py */*NanoStats.txt
        """
}



