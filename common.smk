def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards)
    # single end sample
    return "fastq/{sample}.se.fastq.gz".format(**wildcards)