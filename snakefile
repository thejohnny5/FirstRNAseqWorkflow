import pandas as pd

df = pd.read_csv("config/samples.tsv", sep='\t').set_index("experiments", drop=False)

def split_urls(url: str) -> str:
	#Split Url and return file name at end
    return url.split('/')[-1]

def make_genome_name(genome_fasta: str) -> str:
    #Return genome name from file name in the form GCF_00000000.1
    g_split = genome_fasta.split('_')
    return "_".join([g_split[0], g_split[1]])

def fix_unzip_name(genome_fasta: str) -> str:
    #Remove .gz from end of names
    if genome_fasta.endswith('.gz'):   
        name_list = genome_fasta.split('.')    
        return ".".join(name_list[:-1])
    else:
        return genome_fasta

def make_samples_list(experiment: str) -> list[str]:
    samples = df["sra"][experiment]
    samples = samples.strip("[]")
    samples_list = samples.split(",")
    return samples_list

def make_final_count_rule(df: pd.DataFrame) -> list[str]:
    out_list = []
    for experiment in df.index:
        samples_list = make_samples_list(experiment)
        for sample in samples_list:
            path_out = "workflow/results/experiments/{experiment}/counts/{sample}.counts.txt".format(experiment=experiment, sample=sample.strip())
            out_list.append(path_out)
    return out_list

def make_exp_name(df: pd.DataFrame) -> str:
    exp_names = df.index.values.tolist()
    exp_names = [name[0:5] for name in exp_names]
    return ".".join(exp_names)
#print(make_final_count_rule(df))

def get_accession_links(file_shrunk: str, fasta_gtf: pd.DataFrame) -> str:
    for link in df[fasta_gtf]:
        if file_shrunk in link:
            return link



rule all_deseq2:
	input:
		#Expect Final Count files
		#make_final_count_rule(df)
		#expand("workflow/results/experiments/{experiment}/counts/{sample}.counts.txt", experiment="NMR_iPS", sample=["SRR5985204","SRR5985215"])
		expand("workflow/results/DESeq2/{experiments}/normalized_counts/{experiments}.norm_counts.tsv", experiments=make_exp_name(df))
		#expand("{experiment}/hisat2_ref/{genome_fasta}", experiment="NMR_iPS", genome_fasta="GCF_000247695.1_HetGla_female_1.0_genomic.fna")
	
rule all_counts:
	input:
		make_final_count_rule(df)

rule DESeq2:
	#Run DESEq2 on multiple experiments by calling all the files in the folder
	input:
		#Require Count Matrix with geneids and sample names generated from rule
		raw_counts = "workflow/results/DESeq2/{experiments}/raw_counts/{experiments}.raw_counts.tsv"
	output:
		#Outputs normalized counts and log2 normalized counts for each sample and the means for the experiments
		norm_counts = "workflow/results/DESeq2/{experiments}/normalized_counts/{experiments}.norm_counts.tsv",
		norm_counts_log2 = "workflow/results/DESeq2/{experiments}/normalized_counts/{experiments}.norm_counts_log2.tsv",
		norm_counts_mean = "workflow/results/DESeq2/{experiments}/normalized_counts/{experiments}.mean_norm_counts.tsv",
		norm_counts_mean_log2 = "workflow/results/DESeq2/{experiments}/normalized_counts/{experiments}.mean_norm_counts_log2.tsv",
		#Outputs pairwise comparisons for experiments
		pairwise_dir = directory("workflow/results/DESeq2/{experiments}/pairwise/")
	conda:
		"workflow/envs/deseq2_renv.yaml"
	script:
		#Rscript that runs the analysis
		"workflow/scripts/deseq2.R"

rule make_raw_counts:
	input:
		counts = make_final_count_rule(df)
	output:
		merged_counts = "workflow/results/DESeq2/{experiments}/raw_counts/{experiments}.raw_counts.tsv"
	script:
		"workflow/scripts/join_counts.py"


rule featureCounts:
	#Call counts using the bam file and gtf files
	input:
		#GTF files
		#Calls rule get_gtf
		gtf = lambda wildcards: "workflow/results/ref/gtf/{genome_gtf}".format(genome_gtf=split_urls(df["gtf_url"][wildcards.experiment])),
		#BAM files
		#Calls rule generate_and_index_bams
		bam = "workflow/results/experiments/{experiment}/bam/{sample}.bam"
	output:
		#Final Count file with header removed
		#Tempfile before removing header
		#tempfile = temp("{experiment}/counts/{sample}.temp.txt"),
		counts="workflow/results/experiments/{experiment}/counts/{sample}.counts.txt",
		summary="workflow/results/experiments/{experiment}/counts/{sample}.counts.txt.summary"
	conda:
		"workflow/envs/workflow.yaml"
	shell:
		#featureCounts
		#Remove first line
		#Put in final file
		"""
		featureCounts -p --extraAttributes product -a {input.gtf} -o {output.counts} {input.bam}
		tail -n +2 {output.counts} > {wildcards.experiment}_{wildcards.sample}.temp && mv {wildcards.experiment}_{wildcards.sample}.temp {output.counts}
		"""

rule get_gtf:
	output:
		#GTF file from provided tsv file with header gtf_url
		gtf = "workflow/results/ref/gtf/{genome_gtf}"
	params:
		#Use wildcard parameter to grab specific link from tsv file
		url = lambda wildcards: "{}".format(get_accession_links(wildcards.genome_gtf, "gtf_url")),
		#Url_name fixes the file name so that the base file name is without the .gz extension
		url_name = lambda wildcards: split_urls("{}".format(get_accession_links(wildcards.genome_gtf, "gtf_url")))
	shell:
		#Download url if not already downloaded to path ref/gtf
		#Create symbolic link to output
		"""
		wget -nc -P workflow/results/ref/gtf {params.url}
		zgrep -v 'gene_id ""' workflow/results/ref/gtf/{params.url_name} > {output.gtf}.temp.txt
		gzip {output.gtf}.temp.txt
		rm workflow/results/ref/gtf/{params.url_name}
		mv {output.gtf}.temp.txt.gz workflow/results/ref/gtf/{params.url_name}
		"""

rule generate_and_index_bams:
	#Transform sam file to sorted bam file and index the bam
	input:
		#sam file
		#Depends on rule hisat2_map
		"workflow/results/experiments/{experiment}/sam/{sample}.sam"
	output:
		#Bam files
		bam = "workflow/results/experiments/{experiment}/bam/{sample}.bam",
		bai = "workflow/results/experiments/{experiment}/bam/{sample}.bam.bai"
	conda:
		"workflow/envs/workflow.yaml"
	shell:
		#sort and index
		"""
		samtools sort {input} > {output.bam}
		samtools index {output.bam}
		"""

rule hisat2_map:
	#Map trimmed fastq files to reference genome
	input:
		#Fastq files require rule trim_fastqs
		fastq1 = "workflow/results/experiments/{experiment}/fastqs_trimmed/{sample}_1_val_1.fq.gz",
		fastq2 = "workflow/results/experiments/{experiment}/fastqs_trimmed/{sample}_2_val_2.fq.gz",
		#Reference genome requires rule hisat2_build
		#ref_genome = "{experiment}/hisat2_ref/{genome_fasta}/"
		ref_genome = lambda wildcards: "workflow/results/ref/hisat2_ref/{genome_fasta}".format(experiment=wildcards.experiment,
		genome_fasta=fix_unzip_name(split_urls(df["fasta_url"][wildcards.experiment])))
	params:
		genome_fasta = lambda wildcards: "{genome}".format(genome=fix_unzip_name(split_urls(df["fasta_url"][wildcards.experiment])))
	conda:
		"workflow/envs/workflow.yaml"
	output:
		#Sam file
		temp("workflow/results/experiments/{experiment}/sam/{sample}.sam")
	shell:
		#Map fastqs to indexed hisat2 genome
		"""
		hisat2 -x {input.ref_genome}/{params.genome_fasta} \
		-1 {input.fastq1} \
		-2 {input.fastq2} \
		-S {output}
		"""


rule hisat2_build:
	input:
		#FASTA file
		#requires rule download_fasta
		ref_file = "workflow/results/ref/fasta/{genome_fasta}"
	output:
		#directory in name of genome with hisat2 indexed genome
		directory("workflow/results/ref/hisat2_ref/{genome_fasta}")
	conda:
		"workflow/envs/workflow.yaml"
	shell:
		#Check if directory exists for genome, if not build genome with
		#hisat2-build
		#After genome is built or if it already exists, create symbolic link
		"""
		mkdir -p workflow/results/ref/hisat2_ref/{wildcards.genome_fasta}
		hisat2-build {input.ref_file} workflow/results/ref/hisat2_ref/{wildcards.genome_fasta}/{wildcards.genome_fasta}




		"""

rule download_fasta:
	#Download genomic sequence
	output:
		#Genome in fasta format unzipped
		fasta = "workflow/results/ref/fasta/{genome_fasta}"
	params:
		#Url is from tsv file for the experiment
		url = lambda wildcards: "{}".format(get_accession_links(wildcards.genome_fasta, "fasta_url")),
		#Url_name fixes the file name so that the base file name is without the .gz extension
		url_name = lambda wildcards: fix_unzip_name(split_urls("{}".format(get_accession_links(wildcards.genome_fasta, "fasta_url"))))
	shell:
		#Download genome to specific folder, unzip, and create symbolic link
		"""
		wget -nc -P workflow/results/ref/fasta {params.url}
		[ ! -f "workflow/results/ref/fasta/{params.url_name}" ] && gzip -dc workflow/results/ref/fasta/{params.url_name}.gz > workflow/results/ref/fasta/{params.url_name}
		"""

rule trim_fastqs:
	#trim adapters off fastq files
	input:
		#Raw reads
		#require rule download_and_compress
		fq1 = "workflow/results/experiments/{experiment}/fastq/{sample}_1.fastq.gz",
		fq2 = "workflow/results/experiments/{experiment}/fastq/{sample}_2.fastq.gz"
	output:
		#Trimmed fastq files
		fastq1 = "workflow/results/experiments/{experiment}/fastqs_trimmed/{sample}_1_val_1.fq.gz",
		fastq2 = "workflow/results/experiments/{experiment}/fastqs_trimmed/{sample}_2_val_2.fq.gz"
	params:
		#Output directory for fastqs
		outdir="workflow/results/experiments/{experiment}/fastqs_trimmed"
	conda:
		"workflow/envs/workflow.yaml"
	shell:
		"trim_galore --paired {input.fq1} {input.fq2} --output_dir {params.outdir}"

rule download_and_compress:
	output:
		#Raw reads from SRA
		fq1 = "workflow/results/experiments/{experiment}/fastq/{sample}_1.fastq.gz",
		fq2 = "workflow/results/experiments/{experiment}/fastq/{sample}_2.fastq.gz",
		temp = temp("workflow/results/experiments/{experiment}/sra/{sample}.sra")
	params:
		outdir = "workflow/results/experiments/{experiment}/fastq"
	conda:
		"workflow/envs/workflow.yaml"
	shell:
		#Sample wildcard should be SRA accession number which will download
		#and split fastq files into directory {experiment/fastq}
		"""
		prefetch -o {output.temp} {wildcards.sample}
		fasterq-dump --split-files {output.temp} --outdir {params.outdir}
		gzip {params.outdir}/{wildcards.sample}_1.fastq {params.outdir}/{wildcards.sample}_2.fastq
		"""
