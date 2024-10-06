





def get_ids () {
    samples=channel.empty()
    codef  = new File(params.samples)
    codef.each { line ->
	def parts=line.split("\t")
	data=Channel.of(parts[0]).combine(Channel.fromPath(params.input_path+parts[1]))
	samples=samples.mix(data)
    }
    return samples
}



process fast5_to_pod {
    cpus 2
    maxForks params.max_forks
    input:
	tuple val(base), path(fast5)
    output:
	path(pod)
    publishDir "{params.output}/$base"
    script:
	core=fast5.simpleName
        pod ="${core}.pod5"
        """
          pod5 convert  fast5 -t 2 $fast5 -o $pod
        """
}


workflow {
    all_data=get_ids()
    main:
	fast5_to_pod(all_data)
}
