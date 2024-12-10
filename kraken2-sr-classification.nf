


process sr_kraken2 {
    maxForks 8
    cpus params.kraken_threads
    memory params.kraken_mem
    input:
	    tuple val(base), path(fq)
    output:
	    tuple val(base), path("${base}.txt"), path("${base}_report.txt")
    publishDir params.kraken2_output_dir  
    script: 
       sample_out = "${base}_out.txt"
       sample_report = "${base}_report.txt"
       """
       hostname
       #run Kraken2 without the --report option
       kraken2 --db ${params.db} --memory-mapping --use-names --paired --threads $task.cpus \
	       --output ${base}.txt  --report ${base}_report.txt   ${fq[0]} ${fq[1]}
       """
}

process extract_eukaryotic_reads {
    cpus 2
    input:
       tuple val(base), path(kraken_out), path(kraken_report)
    output:
       tuple val(base), path("${base}_eukaryotic_reads.txt")
    script:
       eukaryotic_reads = "${base}_eukaryotic_reads.txt"
       """
           awk '\$3 == "(taxid 2759)" { print \$2 }' ${kraken_out} > ${base}_eukaryotic_reads.txt
       """
}

process filter_srs {
    cpus 2
    errorStrategy 'finish'
    input:
        tuple val(base), path(read_ids), path(fq1), path(fq2)
    output:
        tuple val(base), path("${base}_filtered_{1,2}.fq.gz")
    script:
        """
        # Filter the original FASTQ file based on the eukaryotic read IDs
        seqtk subseq $fq1 ${read_ids}  |bgzip  > ${base}_filtered_1.fq.gz
	seqtk subseq $fq2 ${read_ids}   |bgzip  > ${base}_filtered_2.fq.gz
        """
}

process MAG_alignment {
	cpus 2
	memory '36GB'
	input:
  	   tuple val(base), path(filtered_fastq), path(mag_files)
	output:
		path("${base}_coverm_output.txt")
	publishDir params.coverm_output_dir

	script: 
	"""
	# Run CoverM to align eukaryotic reads to MAGs
	# filtered_fastq has two parts for forward and reverse ${filtered_fastq[0]} and ${filtered_fastq[1]}
	coverm genome --genome-fasta-directory $mag_files \
                -1 ${filtered_fastq[0]} -2 ${filtered_fastq[1]} \
                -x .fa \
		--mapper bwa-mem2 \
                --output-file ${base}_coverm_output.txt --threads 2
	"""
}

workflow {
    sr_ch = Channel.fromFilePairs(params.sr+"/*{1,2}*")  {it -> it.simpleName.replaceAll("_.*","") }
    mag_ch = Channel.fromPath(params.mags_dir + "/*", type: 'dir')\
	    .map { dir -> 
		def sample = dir.name.replaceAll("dastool_", "")
		[sample, dir.toString()]
		}

  
sr_kraken2(sr_ch)
           |  extract_eukaryotic_reads
           |  cross  (sr_ch)    // produces list of [[base, report], [base, fq1, fq2]]
           |  map { it -> it[0]+it[1][1] } //synactic sugar to look better: make list of [base, report, fq1,fq2]
           |  filter_srs
           |  join(mag_ch)
		| view()
           |  MAG_alignment
}
	
	
