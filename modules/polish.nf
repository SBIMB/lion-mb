

process index_map_samples {
  cpus 16
  input:
     tuple val(base), path(sr1), path(sr2), path(contig)
  output:
     tuple val(base), path(contig), path(sam1), path(sam2)
  script:
    sam1= "${base}_1.sam"
    sam2= "${base}_2.sam"
    """
      bwa-mem2 index $contig
      bwa-mem2 mem -t 16 -a $contig $sr1 -Obam -o $sam1
      bwa-mem2 mem -t 16 -a $contig $sr2 -Obam -o $sam2

    """
}

process polypolish {
  maxForks 10
  input:
      tuple val(base), path(contig), path(sam1), path (sam2)
  output:
      tuple val(base), path("${cname}-polished.fa")
  script:
    cname=contig.baseName
    """
    polypolish filter --in1 $sam1 --in2 $sam2 --out1 filt_1.sam --out2 filt_2.sam
    polypolish polish $contig filt_1.sam filt_2.sam > ${cname}-polished.fa
    """
}


process split_samples {
   input:
      tuple val(base), path(contigs)
   output:
      path("indivs/*")
   script:
      """
      mkdir indivs
      seqkit split -s 1 $contigs -O indivs
       rnr -f "assembly.part_(.*).fasta" '${base}-\${1}.fa' indivs/*
      """
}

workflow short_read_polish {
   take:
      samples
      sr
   main:
       contigs =split_samples (samples)
          | flatten
          | map { it -> [it.simpleName.replaceAll(".part.*",""),it] }
	sr.cross(contigs) 
        |  map { [it[0][0], it[0][1][0], it[0][1][1], it[1][1] ] }
	| index_map_samples
	| polypolish
   emit:
      polypolish.out





}


process medaka {
   // incomplete as change to unicycler obviates
   cpus 16
   input:
      tuple val(base), path(contig), path(lr)
   script:
     """
     mkdir indivs
     mkdir medaka
     seqkit split -s 1 $contig -O indivs
     for seq in indivs/*fasta; do
         medaka_consensus -m  r941_e81_sup_g514 --bacteria -i $lr -d \$seq -o ${OUTDIR} \
                           -t 16 -o medaka
     """
}
