#! /usr/bin/env python

configfile: "config.yaml"
metadata = [r.strip('\n').split('\t') for r in open(config['run_table'], 'r')]
header = metadata.pop(0)
dpath = {r[header.index('Run')]:r[header.index('download_path')] for r in metadata}
run2name={r[header.index('Run')]:r[header.index('SampleName')] for r in metadata}
sp_samp= {r[header.index('SampleName')]:r[header.index('Sp')] for r in metadata}
out_f= {r[header.index('SampleName')]:r[header.index('outFolder')] for r in metadata}
fold2sp={r[header.index('outFolder')]:r[header.index('Sp')] for r in metadata}
vstout_ext=["micX","info","exskX","eej2","MULTI3X"]
uniq_outFold=set(out_f.values())

concat =dict()
for k,v in run2name.items():
    if v in concat:
        concat[v].append(k)
    else:
        concat[v] = [k]

samples_fold =dict()
for k,v in out_f.items():
    if v in samples_fold:
        samples_fold[v].append(k)
    else:
        samples_fold[v] = [k]

num_samp=dict()
for k,v in samples_fold.items():
    num_samp[k]=str(len(v))

sp_name=dict()
for k,v in fold2sp.items():
    sp_name[k]=config['translateSP'][v]

rule all:
    input:
        expand("{r}/Tidyvasttools-minFr_0.8-minSD_5-noVLOW-min_ALT_use25-Tidy.tab", r=set(uniq_outFold)),
        expand("{r}/Salmon_Gene_table_counts.csv", r=set(uniq_outFold)),
        expand("{r}/Salmon_Transcript_table_counts.csv", r=set(uniq_outFold)),
        expand("{r}/Salmon_Gene_table_TPM.csv", r=set(uniq_outFold)),
        expand("{r}/tidy_stats.log", r=set(uniq_outFold))

rule copySRA:
    output:
        sra = "SRA/{runid}.sra"
    params:
        url = lambda x: dpath[x.runid]
    shell:
        "cp {params.url} {output.sra} "

rule fastq_dump:
    input:
        "SRA/{runid}.sra"
    output:
        r1 = "SRA/{runid}_1.fastq.gz",
        r2 = "SRA/{runid}_2.fastq.gz"
    params:
        dir = "SRA/",
        name = lambda x: x.runid
    threads: 2
    shell:
        "parallel-fastq-dump"
        " -t {threads}"
        " --tmpdir ./"
        " --gzip"
        " -O {params.dir}"
        " --split-3"
        " -F -B -Q 33 "
        " --skip-technical"
        " --defline-qual '+'"
        " -s {input} 2 &> {params.dir}/{params.name}_fastqdump.log && rm {input}"

rule concatenate:
    input:
        fastqs1 = lambda wildcards:
            expand("SRA/{runid}_1.fastq.gz", runid=concat[wildcards.name]),
        fastqs2 = lambda wildcards:
            expand("SRA/{runid}_2.fastq.gz", runid=concat[wildcards.name])
    output:
         ofast1= "concat/{name}_1.fastq.gz",
         ofast2= "concat/{name}_2.fastq.gz"
    params:
        dir = "concat/{name}"
    shell:
        "cat {input.fastqs1} > {output.ofast1} && cat {input.fastqs2} > {output.ofast2} && rm {input.fastqs1} {input.fastqs2}"

rule vast_tools_align:
    input:
        r1= "concat/{name}_1.fastq.gz",
        r2= "concat/{name}_2.fastq.gz"
    output:
        "{outFolder}/to_combine/{name}.micX",
        "{outFolder}/to_combine/{name}.info",
        "{outFolder}/to_combine/{name}.exskX",
        "{outFolder}/to_combine/{name}.eej2",
        "{outFolder}/to_combine/{name}.MULTI3X"
    params:
        outFolder=lambda x: out_f[x.name],
        sp=lambda x: sp_samp[x.name],
        name=lambda x: x.name
    threads: workflow.cores * 0.5
    shell:
        "/users/mirimia/projects/vast-tools/vast-tools align"
        " {input.r1} {input.r2}"
        " -c 1"
        " --name {params.name}"
        " --output {params.outFolder}"
        " --sp {params.sp} --noIR &> {params.outFolder}/{params.name}_VASTTOLS.log"

rule vast_tools_combine:
    input:
        lambda wildcards:
             expand(wildcards.outFolder+"/to_combine/{name}.{exts}", zip, name=samples_fold[wildcards.outFolder], exts=vstout_ext),
    output:
        e="{outFolder}/INCLUSION_LEVELS_FULL-{sp2name}-{nsamp}.tab",
    params:
        outFolder="{outFolder}/",
        sp=lambda x: fold2sp[x.outFolder]
    shell:
        "/users/mirimia/projects/vast-tools/vast-tools combine"
        " --cores 1"
        " -o {params.outFolder}"
        " --sp {params.sp}"
        " --IR_version 2 -no_expr --onlyEX &> {params.outFolder}/COMBINE.log"

rule vast_tools_tidy:
    input:
        lambda wildcards: wildcards.outFolder+"/INCLUSION_LEVELS_FULL-"+ sp_name[wildcards.outFolder] + "-"+ num_samp[wildcards.outFolder]+".tab"
    output:
        tidy="{outFolder}/Tidyvasttools-minFr_0.8-minSD_5-noVLOW-min_ALT_use25-Tidy.tab",
        stats="{outFolder}/tidy_stats.log"
    params:
        outFolder=lambda x: set(uniq_outFold)
    shell:
        "/users/mirimia/projects/vast-tools/vast-tools tidy"
        " {input}"
        " -min_Fr 0.8"
        " -outFile {output.tidy}"
        " --noVLOW --log 1> {output.stats}"

rule downloadANN:
    output:
        "annotation_{spec}.gtf.gz"
    params:
        url2=lambda wildcards: config['annotation'][wildcards.spec]
    shell:
        "wget -O {output} {params.url2}"

rule downloadTR:
    output:
        "{spec}_transcriptome.fa"
    params:
        url1=lambda wildcards: config['transcriptome'][wildcards.spec]
    shell:
        "wget -O {output}.gz {params.url1} &&"
        " gunzip {output}.gz"

rule salmon_index:
    input:
        tr="{spec}_transcriptome.fa",
    output:
        config['directories']['references']+"Salmon/{spec}/pos.bin"
    params:
        d=config['directories']['references']+"Salmon/{spec}/"
    shell:
        "salmon index -t {input.tr} -i {params.d} -k 29 &&"
	    " rm {input.tr}"

rule transcript2gene:
    input:
        "annotation_{spec}.gtf.gz"
    output:
        config['directories']['references']+"Salmon/{spec}_transcript2gen.txt"
    params:
        refdir= config['directories']['references']+"Salmon/"
    shell:
        "mkdir -p {params.refdir} &&"
        " Rscript --vanilla /users/mirimia/liniguez/scripts/make_tx2gene.R {input} {output}&&"
        " rm {input} "

rule salmonMapping:
    input:
        r1= "concat/{name}_1.fastq.gz",
        r2= "concat/{name}_2.fastq.gz",
        idx = lambda x: config['directories']['references']+"Salmon/"+fold2sp[x.outFolder]+"/pos.bin"
    output:
        tr="{outFolder}/Salmon/{name}/quant.sf"
    params:
        idxTRUE = lambda x: config['directories']['references']+"Salmon/"+fold2sp[x.outFolder],
        dir="{outFolder}/Salmon/{name}"
    threads: workflow.cores * 0.5
    shell:
        "salmon quant -l A -i {params.idxTRUE} -1 {input.r1} -2 {input.r2} -p {threads}"
        " -o {params.dir} --seqBias --gcBias &> {params.dir}/SALMON.log "

rule Salmon_tables:
    input:
        lambda wildcards:
             expand(wildcards.outFolder+"/Salmon/{sampleName}/quant.sf", sampleName=samples_fold[wildcards.outFolder]),
        experiments = config['run_table'],
        t2g = lambda x: config['directories']['references']+"Salmon/"+fold2sp[x.outFolder]+"_transcript2gen.txt"
    output:
        gc = "{outFolder}/Salmon_Gene_table_counts.csv",
        tc = "{outFolder}/Salmon_Transcript_table_counts.csv",
        gtpm = "{outFolder}/Salmon_Gene_table_TPM.csv"
    params:
        fold= "{outFolder}"
    shell:
        "Rscript --vanilla /users/mirimia/liniguez/scripts/merge_tables.R {input.experiments} {params.fold} {input.t2g}"

