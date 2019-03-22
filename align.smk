def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards)
    # single end sample
    return "fastq/{sample}.se.fastq.gz".format(**wildcards)


rule star:
    input:
        get_trimmed
    output:
        "01alignments/{sample}/{sample}.bam"
    log:
        "00log/alignments/{sample}.log"
    params:
        tmp_bam = "01alignments/{sample}/Aligned.sortedByCoord.out.bam",
        out_dir = "01alignments/{sample}/",
        # path to STAR reference genome index
        index   = config["ref"]["index"]
    threads: 10
    shadow: "minimal"
    run: 
        shell("""
        STAR --genomeDir {params.index} \
        --runThreadN {threads} \
        --readFilesIn {input} \
        --outFileNamePrefix {params.out_dir} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard \
        --outFilterMultimapNmax 1 \
        --bamRemoveDuplicatesType UniqueIdentical \
        --readFilesCommand zcat \
        --outStd Log > {log} 2>&1

        # Remove duplicates marked by STAR with the flag 0x400
        samtools view -@ {threads} -b -F 0x400 {params.tmp_bam} > {output}
        samtools index {output}
        """)


rule featureCounts:
    input:
         "01alignments/{sample}/{sample}.bam"
    output:
        featureCounts = "02featureCounts/{sample}/{sample}.featureCounts",
        counts        = "02featureCounts/{sample}/{sample}.counts",
        rpkm          = "02featureCounts/{sample}/{sample}.rpkm"
        # tpm           = "02featureCounts/{sample}/{sample}.tpm"
    log:
        "00log/featureCounts/{sample}.log"
    params:
        tmp      = "02featureCounts/{sample}/{sample}.seqDepth",
        gtf      = config["ref"]["annotation"],
        options  = config["params"]["featureCounts"],
        # Add -p option for pair-end data if it's the case
        pair_end = lambda w: "-p" if not is_single_end(w.sample) else str()
    threads: 10
    run:
        shell("""
        featureCounts -T {threads} \
        -a {params.gtf} \
        -o {output.featureCounts} \
        {params.pair_end} {params.options} {input} > {log} 2>&1

        # Create a file with the counts, another for RPKM and another for TPM
        cut -f 1,7 {output.featureCounts} | tail -n +2 > {output.counts}
        seqdepth=$(awk '{{sum+=$2}}END{{print sum}}' {output.counts})
        awk -v s=${{seqdepth}} 'BEGIN{{OFS="\t"}}NR>2{{print $1,(($7)*1000000*1000)/($6*s)}}' {output.featureCounts} > {output.rpkm}
        """)