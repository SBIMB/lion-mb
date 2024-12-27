



process versions {
    output:
	path("*.txt")
    script:
	"""
        OPERA-MS.pl | head -n2 > mag-versions.txt
        flye --version >> mag-versions.txt
        unicyler --version >> mag-versions.txt
    """
}


process flye_assemble {
    cpus 10
    memory '120GB'
    input:
	tuple val(base), path(fqs)
    output:
	tuple val(base), path("${outdir}/*")
    publishDir "${params.mags_dir}/metaflye"    
    script:
      outdir="${base}-nano-assemble"
      """
      #!/bin/bash
      mkdir -p $outdir
      flye --meta  -g 100m -t 10  --out-dir $outdir --nano-hq  $fqs 
      echo flye --meta  -g 100m -t 10  --out-dir $outdir --nano-hq > flye.txt
      flye --version >> flye.txt
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


process opera_ms_hybrid {
    cpus 10
    memory 32.GB
    afterScript "/bin/rm lr.fq"
    input:
	tuple val(base), path(lr), path(sr1), path(sr2)
    output:
	tuple val(base), file("$base/*")
    publishDir "${params.mags_dir}/opera-ms"        
    script:
    """	
      hostname
      zcat $lr > lr.fq
      /usr/bin/time -f "%e %M" \
           OPERA-MS.pl --short-read1 $sr1  --short-read2 $sr2 --long-read lr.fq \
                  --no-ref-clustering  --no-strain-clustering --no-polishing \
                  --num-processors 10 --out-dir $base
    """
}

process unicycler_hybrid {
    cpus 10
    memory '180.GB'
    errorStrategy 'finish'
    input:
	tuple val(base), path(lr), path(sr1), path(sr2)
    output:
        path("$base/*")
    publishDir "${params.mags_dir}/unicyler"
    script:
	 """
         hostname
         mkdir -p $base
         /usr/bin/time -f "%e %M" unicycler -t 10 --mode normal --min_fasta_length 1000 \
                    -1 $sr1 -2 $sr2 -l $lr -o $base
         """
}

process kraken_raw {
   cpus 16
   memory "240G"
   clusterOptions " -x n01,n02,n03,n04,n08,n10,n11,n14,n15,n24,n29,n30,n31,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n45"
   errorStrategy 'finish'
   input:
      tuple val(base), path(input)
   output:
      tuple val(base), path('contigs.fasta'), path(kraken), path(report)
   script:
      kraken="${base}.kraken"
      report="${base}-kraken.txt"
      """
        hostname
	kraken2 --confidence 0.4 --db ${params.kraken_db} --threads 16 contigs.fasta --output ${kraken} --report ${report} --use-names
      """

}




process clean_kraken {
   errorStrategy 'finish'
   input:
       tuple val(base), path(input), path(kraken), path(report)
   output:
      tuple val(base),  path("${base}-virus.fa"), emit: virus
      tuple val(base), path(prok_fa), emit: prokaryote
   script:
      virus_fa = "${base}-virus.fa"
      prok_fa  = "${base}-prok.fa"      
   """
	extract_kraken_reads.py -k ${kraken} -r ${report} -s $input -o vraw.fa --include-children -t 10239
        seqkit seq vraw.fa -m 2000 -o vraw2.fa
        annotate_fasta_header.py vraw2.fa $kraken ${virus_fa}
	extract_kraken_reads.py -k ${kraken} -r ${report} -s $input -o braw.fa --include-children -t 2
        seqkit seq braw.fa -m 100000 -o braw2.fa
        annotate_fasta_header.py braw2.fa $kraken bacteria.fa
	extract_kraken_reads.py -k ${kraken} -r ${report} -s $input -o araw.fa --include-children -t 2157
        seqkit seq araw.fa -m 100000 -o araw2.fa
        annotate_fasta_header.py araw2.fa $kraken archaea.fa
        cat bacteria.fa archaea.fa > $prok_fa
   """
}


process checkv {
   cpus 10
   input:
	tuple val(base), path(mags)
   output:
        path("${base}-virus.fa")
   errorStrategy 'finish'
   """
     export PATH=${params.diamond_path}:\$PATH
     checkv end_to_end -t 10 $mags quality
     if [ -e quality/quality_summary.tsv ]; then 
       extract_qual_virus.py quality/quality_summary.tsv
       seqkit grep -f quality_viruses.txt $mags -o ${base}-virus.fa
     else
       touch ${base}-virus.fa
     fi
   """

}


process checkm2 {
   cpus 10
   input:
	tuple val(base), path(mags)
   output:
        path("${base}-bact.fa")
   errorStrategy 'finish'
   script:
      report = "${base}.txt"
      """
      mkdir bact temp
      seqkit split -i $mags  -O bact 
      for s in bact/*; do 
          checkm2 predict --threads 10 --input \$s   --remove_intermediates -o temp 
          grep -v GC_Content temp/quality_report.tsv >> report.txt
          /bin/rm -rf temp/*
          gunc run --threads 10 --input_fasta \$s --out_dir temp
          tail -n1 temp/*tsv >> gunc.txt
          /bin/rm -rf temp/*
      done
      extract_seqs_cc.py bact report.txt  gunc.txt  ${params.mag_completeness} ${params.mag_contamination} seq_ids
      seqkit grep -f seq_ids $mags -o ${base}-bact.fa
      """
}
      
workflow split_kingdoms {
    take: samples
    main:
       kraken_raw(samples)  | clean_kraken
       checkm2(clean_kraken.out.prokaryote)
       checkv(clean_kraken.out.virus)
}




workflow {
    lr_ch = Channel.fromPath(params.lr_qc).map { it -> [it.simpleName, it]}
    sr_ch = Channel.fromFilePairs(params.sr_qc)	{  it.name.replaceAll("_.*","") }
    main:
    lr_ch.join(sr_ch)
         .map { [it[0], it[1], it[2][0], it[2][1]] } // flatten list \
         |(unicycler_hybrid & opera_ms_hybrid)
    flye_assemble(lr_ch)
    split_kingdoms (opera_ms_hybrid.out)
}
