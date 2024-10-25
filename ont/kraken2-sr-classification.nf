
process sr_kraken2 {
	cpus 30
	
	input:
		tuple val(base), path(fq)
	
	output:
		tuple val(base), path("${base}_out.txt"), path("${base}_report.txt")
	
	publishDir params.kraken2_output_dir
	
	....
	script: 
		sample_out = "${base}_out.txt"
		sample_report = "${base}_report.txt"
	"""
	#run Kraken2 without the --report option
	kraken2 --db ${params.db} --memory-mapping --use-names --paired --threads 30 --output ${sample_out} ${fq[0]} ${fq[1]}

	#optional: keep a seperate summary report if needed
	kraken2 --db ${params.db} --memory-mapping --use-names --paired --threads 30 \
	--report ${params.report_out_dir}/${base}_report.txt \
	"${fq[0]}" "${fq[1]}"
	"""
}

process extract_eukaryotic_reads {
	cpus 4
	input:
		tuple val(base), path(kraken_out)

	output:
		tuple val(base), path("${base}_eukaryotic_reads.txt")

	script:
		eukaryotic_reads = "${base}_eukaryotic_reads.txt"
    """
    # Extracting reads classified as eukaryotic (using taxonomy ID 2759 for Eukaryota)
    awk '$3 == "2759" { print $2 }' ${kraken_out} > ${eukaryotic_reads}
    """
}

process filter_srs {
    cpus 4

    input:
        tuple val(base), path(fq), path(read_ids)

    output:
        path("${base}_filtered.fastq")

    script:
        """
        # Filter the original FASTQ file based on the eukaryotic read IDs
        seqtk subseq ${fq} ${read_ids} > ${base}_filtered.fastq
        """
}

process MAG_alignment {
	cpus 10
	memory '36GB'

	input:
		tuple val(base), path(filtered_fastq),
		path mag_files

	output:
		path("${base}_coverm_output.txt")

	script: 
		base = euk_reads.simpleName
	"""
	# Run CoverM to align eukaryotic reads to MAGs
	coverm genome --genome-files $mag_files --reads ${filtered_fastq} --out-file ${base}_coverm_output.txt --threads 10
	"""
}

workflow {
	sr_ch = Channel.fromFilePairs(params.sr_dir+"/*{1,2}*") \
		{it -> it.name.replaceAll("_.*","") }
	mag_ch = Channel.fromPath(params.mags_dir+"/*.fa")

	main:
	sr_kraken2(sr_ch) | extract_eukaryotic_reads | filter_srs | MAG-alignment(mag_ch)
}
	
	
