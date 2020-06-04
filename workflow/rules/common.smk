def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])


# Get raw or trimmed reads based on trimming configuration
def get_fq(wildcards):
    if config["trimming"]:
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/trimmed{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)


# Get raw or trimmed reads based on trimming configuration. Used for fastqc
def get_fq_forward(wildcards):
    if config["trimming"]:
        if not is_single_end(**wildcards):
            # paired-end sample
            return "{tmp}/fastq/trimmed{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return "{tmp}/fastq/{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)