def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()

rule cp_fastq_pe:
    input: 
        get_fastq
    output: 
        fastq1=temp("fastq/{sample}-{lane}.1.fastq.gz"),
        fastq2=temp("fastq/{sample}-{lane}.2.fastq.gz")
    message: 
        "Copying fastq files {input}"
    run:
        file1,file2 = input
        shell("""
            ln -s {file1} {output.fastq1}
            ln -s {file2} {output.fastq2}
            """)


rule cp_fastq_se:
    input: 
        get_fastq
    output: 
        temp("fastq/{sample}-{lane}.fastq.gz"),
    message: 
        "Copying fastq files {input}"
    run:
        shell("""
            ln -s {input} {output}
            """)


rule fastp_pe:
	input:
		fw = lambda w: expand("fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=units.loc[w.sample].itertuples()),
		rv = lambda w: expand("fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=units.loc[w.sample].itertuples())
	output: 
		fastq1 = "fastq/{sample}.1.fastq.gz",
		fastq2 = "fastq/{sample}.2.fastq.gz"
	log: 
		"00log/fastp/{sample}.log"
	threads: 
		CLUSTER["fastp"]["cpu"]
	params: 
		fastp_params = config["params"]["fastp-pe"],
		tmp_fw       = "{sample}.1.fastq.tmp.gz",
		tmp_rv       = "{sample}.2.fastq.tmp.gz"
	message: 
		"Processing fastq files from {input}"
	shadow: "minimal"
	benchmark:
		".benchmarks/{sample}.merge_fastqs.benchmark.txt"
	run:
		shell("""
			cat {input.fw} > fastq/{params.tmp_fw}
			cat {input.rv} > fastq/{params.tmp_rv}
			/hpcnfs/data/DP/software/fastp -A -z 4 \
			-i fastq/{params.tmp_fw} \
			-I fastq/{params.tmp_rv} \
			-o {output.fastq1} \
			-O {output.fastq2} \
			-w {threads} \
			{params.fastp_params} 2> {log}
		    """)


rule fastp_se:
	input:
		lambda w: expand("fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=units.loc[w.sample].itertuples()),
	output: 
		"fastq/{sample}.se.fastq.gz"
	log: 
		"00log/fastp/{sample}.log"
	threads: 
		CLUSTER["fastp"]["cpu"]
	params: 
		fastp_params = config["params"]["fastp-se"],
		tmp_fw       = "{sample}.se.fastq.tmp.gz"
	message: 
		"Processing fastq files from {input}"
	shadow: "minimal"
	benchmark:
		".benchmarks/{sample}.merge_fastqs.benchmark.txt"
	run:
		shell("""
			cat {input} > fastq/{params.tmp_fw}
			/hpcnfs/data/DP/software/fastp -A -z 4 \
			-i fastq/{params.tmp_fw} \
			-o {output} \
			-w {threads} {params.fastp_params} 2> {log}
		    """)