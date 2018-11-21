rule call_narrow_peaks:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{sample}_peaks.narrowPeak"
    message:
        "Calling narrowPeak for {wildcards.sample}"
    params:
        name        = "{sample}",        #this option will give the output name, has to be similar to the output
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{sample}_peaks.narrowPeak.log"
    conda:
        "../envs/macs2_env.yaml"
    shell:
        """
        macs2 callpeak -t {input} {params.format} {params.genomesize} --name {params.name} --bdg -q {params.qvalue} --nomodel --extsize 147 --outdir results/bed/ &>{log}
        """

rule call_broad_peaks:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{sample}_peaks.broadPeak"
    message:
        "Calling broadPeak for {wildcards.sample}"
    params:
        name        = "{sample}",
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{sample}_peaks.broadPeak.log"
    conda:
        "../envs/macs2_env.yaml"
    shell:
        """
        macs2 callpeak -t {input} {params.format} --broad --broad-cutoff 0.1 {params.genomesize} --name {params.name} --bdg -q {params.qvalue} --nomodel --extsize 147 --outdir results/bed/ &>{log}
        """
