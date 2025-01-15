
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
		tuple val(base), path(kraken_out), path(kraken_report), path(fq_forward), path(fq_reverse)

	output:
		tuple val(base), path("${base}_eukaryotic_reads.txt")

	publishDir params.eukaryotic_reads_output_dir
	
	script:
		eukaryotic_reads = "${base}_eukaryotic_reads.txt"
    """
    # Extracting reads classified as eukaryotic (using taxonomy ID 2759 for Eukaryota)
    extract_kraken_reads.py -k ${kraken_out} -s1 ${fq_forward} -s2 ${fq_reverse} \
    -o ${base}_euk_reads1.fq \
    -o2 ${base}_euk_reads2.fq \
    -t 2759 \
    --include-children
    """
}

workflow {
    sr_ch = Channel.fromFilePairs(params.sr+"/*{1,2}*")  {it -> it.simpleName.replaceAll("_.*","") }
                
//run Kraken2    
kraken2_ch = sr_kraken2(sr_ch)
//view
kraken2_ch.view { it }

//add srs to channel
kraken2_combined_reads_raw_ch = kraken2_ch.cross(sr_ch)    // produces list of [[base, report], [base, fq1, fq2]]
//flatten
kraken2_combined_reads_ch = kraken2_combined_reads_raw_ch.map { it -> it[0]+it[1][1] }
//view
kraken2_combined_reads_ch.view {it}

//extract eukaryotics reads
eukaryotic_reads_ch = extract_eukaryotic_reads(kraken2_combined_reads_ch)
//view
eukaryotic_reads_ch.view { it }

//extract microbial reads 
// microbial_reads_ch = extract_microbial_reads(kraken2_combined_reads_ch)
//view
//microbial_reads_ch.view { it }
}

