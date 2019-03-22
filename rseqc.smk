## RSEQC

rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"]
    output:
        bed="qc/rseqc/annotation.bed",
        db=temp("qc/rseqc/annotation.db")
    log:
        "00log/rseqc_gtf2bed.log"
    conda: "/hpcnfs/home/ieo4462/.conda/envs/DPbase"
    script:
        "scripts/gtf2bed.py"
       
  
rule rseqc_junction_annotation:
    input:
        bam="01alignments/{sample}/{sample}.bam",
        bed="qc/rseqc/annotation.bed"
    output:
        "qc/rseqc/{sample}.junctionanno.junction.bed"
    priority: 1
    log:
        "00log/rseqc/rseqc_junction_annotation/{sample}.log"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="qc/rseqc/{sample}.junctionanno"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"

        
rule rseqc_junction_saturation:
    input:
        bam="01alignments/{sample}/{sample}.bam",
        bed="qc/rseqc/annotation.bed"
    output:
        "qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf"
    priority: 1
    log:
        "00log/rseqc/rseqc_junction_saturation/{sample}.log"
    params:
        extra=r"-q 255", 
        prefix="qc/rseqc/{sample}.junctionsat"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "01alignments/{sample}/{sample}.bam",
    output:
        "qc/rseqc/{sample}.stats.txt"
    priority: 1
    log:
        "00log/rseqc/rseqc_stat/{sample}.log"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer:
    input:
        bam="01alignments/{sample}/{sample}.bam",
        bed="qc/rseqc/annotation.bed"
    output:
        "qc/rseqc/{sample}.infer_experiment.txt"
    priority: 1
    log:
        "00log/rseqc/rseqc_infer/{sample}.log"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

        
rule rseqc_innerdis:
    input:
        bam="01alignments/{sample}/{sample}.bam",
        bed="qc/rseqc/annotation.bed"
    output:
        "qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        "00log/rseqc/rseqc_innerdis/{sample}.log"
    params:
        prefix="qc/rseqc/{sample}.inner_distance_freq"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="01alignments/{sample}/{sample}.bam",
        bed="qc/rseqc/annotation.bed"
    output:
        "qc/rseqc/{sample}.readdistribution.txt"
    priority: 1
    log:
        "00log/rseqc/rseqc_readdis/{sample}.log"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "01alignments/{sample}/{sample}.bam"
    output:
        "qc/rseqc/{sample}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        "00log/rseqc/rseqc_readdup/{sample}.log"
    params:
        prefix="qc/rseqc/{sample}.readdup"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"

        
rule rseqc_readgc:
    input:
        "01alignments/{sample}/{sample}.bam"
    output:
        "qc/rseqc/{sample}.readgc.GC_plot.pdf"
    priority: 1
    log:
        "00log/rseqc/rseqc_readgc/{sample}.log"
    params:
        prefix="qc/rseqc/{sample}.readgc"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"
        

rule multiqc:
    input:
        expand("01alignments/{sample}/{sample}.bam", sample = SAMPLES),
        expand("qc/rseqc/{sample}.junctionanno.junction.bed", sample = SAMPLES),
        expand("qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf", sample = SAMPLES),
        expand("qc/rseqc/{sample}.infer_experiment.txt", sample = SAMPLES),
        expand("qc/rseqc/{sample}.stats.txt", sample = SAMPLES),
        expand("qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt", sample = SAMPLES),
        expand("qc/rseqc/{sample}.readdistribution.txt", sample = SAMPLES),
        expand("qc/rseqc/{sample}.readdup.DupRate_plot.pdf", sample = SAMPLES),
        expand("qc/rseqc/{sample}.readgc.GC_plot.pdf", sample = SAMPLES),
        expand("00log/rseqc/rseqc_junction_annotation/{sample}.log", sample = SAMPLES)
    output:
        "qc/multiqc_report.html"
    log:
        "00log/multiqc.log"
    params:
        log_name = "multiQC_log"
    shell:
        """
        multiqc {input} -o 10multiQC -f -v -n {params.log_name} 2> {log}
        """
