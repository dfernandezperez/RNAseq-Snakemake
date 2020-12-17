# the commands of the index rule were obtained from here: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
rule salmon_index:
    input:
        primary_assembly = config["ref"]["assembly"],
        transcripts      = config["ref"]["transcriptome"],
    output:
        index = directory("results/02_salmon/salmon_index"),
    log:
        "results/00log/salmonIndex/log"
    threads:
        CLUSTER["salmon_index"]["cpu"]
    shell:
        """
        grep "^>" <(gunzip -c {input.primary_assembly}) \
        | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt
        cat {input.transcripts} {input.primary_assembly} > gentrome.fa.gz
        salmon index -t gentrome.fa.gz -d decoys.txt -p {threads} -i {output} --gencode
        rm decoys.txt gentrome.fa.gz
        """


def set_reads(wildcards, input):
        n = len(input.fastq)
        if n == 1:
            reads = "-r {}".format(*input.fastq)
            return reads
        else:
            reads = "-1 {} -2 {}".format(*input.fastq)
            return reads

rule salmon_quant:
    input:
        index = "results/02_salmon/salmon_index",
        fastq = get_fq,
    output:
        quant       = "results/02_salmon/{sample}/quant.sf",
        quant_genes = "results/02_salmon/{sample}/quant.genes.sf",
        sam         = "results/02_salmon/{sample}/{sample}.sam",
    log: 
        "results/00log/salmonQuant/{sample}.log"
    params:
        gtf           = config["ref"]["annotation"],
        options       = config["params"]["salmon"],
        library       = config["params"]["salmon_library"],
        out_fold      = "results/02_salmon/{sample}",
        reads  	      = set_reads,
    threads:    
        CLUSTER["salmon_quant"]["cpu"]
    shell:
        """
        salmon quant \
        -i {input.index} \
        -p {threads} \
        -g {params.gtf} \
        -l {params.library}\
        {params.reads} \
        -o {params.out_fold} \
        {params.options} \
        2> {log} \
        > {output.sam}
        """


# # Most star parameters taken from https://www.biorxiv.org/content/biorxiv/early/2019/10/31/657874.full.pdf
rule star:
    input:
        get_fq
    output:
        bam   = "results/02alignments/{sample}/{sample}.bam",
        log   = "results/02alignments/{sample}/Log.final.out"
    log:
        align   = "results/00log/alignments/{sample}.log",
        rm_dups = "results/00log/alignments/rm_dup/{sample}.log",
    params:
        out_dir      = "results/02alignments/{sample}/",
        star_params  = config["params"]["star"],
        # path to STAR reference genome index
        index        = config["ref"]["index"],
        samtools_mem = config["params"]["samtools_mem"]
    threads:
        CLUSTER["star"]["cpu"]
    shadow: 
        "minimal"
    shell: 
        """
        STAR --genomeDir {params.index} \
        --runThreadN {threads} \
        --readFilesIn {input} \
        --outFileNamePrefix {params.out_dir} \
        --outSAMtype None \
        --outStd Log \
        --quantMode TranscriptomeSAM \
        --outSAMunmapped Within \
        --quantTranscriptomeBan Singleend \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --outFilterMultimapNmax 20 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        {params.star_params} > {log.align}

        samtools view -h {params.out_dir}Aligned.toTranscriptome.out.bam \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb - > {output.bam} 2>> {log.align}
        """

# rule bam2bigwig:
#     input:
#         rules.featureCounts.output.annot_sam
#     output:
#         "results/05bigwig/{sample}.bw"
#     params:
#         tmp = "results/05bigwig/{sample}_tmp"
#     log:
#         "results/00log/bam2bigwig/{sample}.log"
#     threads:
#         CLUSTER["bam2bigwig"]["cpu"]
#     shell:
#         """
#         fgrep -v "Unassigned_" {input} | samtools view -@ {threads} -Sb - | samtools sort -@ {threads} -o {params.tmp}
#     samtools index {params.tmp}
#     bamCoverage --normalizeUsing CPM -p {threads} -bs 1 -b {params.tmp} -o {output} 2> {log}
#         rm {params.tmp} {params.tmp}.bai
#         """
