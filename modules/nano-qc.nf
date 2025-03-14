


process nano_qc_plot {
    input:
        val(keywords)
	path(fq)
    output:
        path(qc_out)
    publishDir "output", saveAs: { qc_dir == "" ? null : "$qc_dir/$tag/$it" }
    script:
        base   =  fq.simpleName
	qc_out =  "${base}-qc"
	tag    =  keywords['tag']
	qc_dir =  keywords['qc_dir']
        """
        NanoPlot --fastq $fq --outdir ${qc_out} -p ${base}-${tag} 
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



