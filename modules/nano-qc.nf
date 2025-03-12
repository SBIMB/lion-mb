

process nano_qc_plot {
    input:
	path(fq)
    output:
        path(qcdir)
	publishDir params.qc_dir, enabled: params.qc_dir != ""
    script:
        base=fq.simpleName
	qcdir="${base}-qc"
	status=params.status
        """
        NanoPlot --fastq $fq --outdir $qcdir -p ${base}-${status} 
        """
}

process qc_summary {
    input:
	path f
    output:
	path "lr_qc.tex"
    script:
	"""
        summarise_qc.py */*NanoStats.txt
        """
}



