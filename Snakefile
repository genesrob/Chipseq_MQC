import glob
sample_list = []
for file in glob.iglob('Novogene/Chipseq_1.12.20/Fastq/*.fq.gz'):
     split_file = file.split('/')
     file_name = split_file[-1]
     sample_name = file_name.split('.')[0]
     sample_list.append(sample_name)

sample_list_BAM =[]
for file in sample_list:
    sample_name = file.split('_')[0:-1]
    sample_name_str = '_'.join(sample_name)
    sample_list_BAM.append(sample_name_str)
    sample_list_BAM = list(set(sample_list_BAM))

QC_html= expand("Novogene/Chipseq_1.12.20/QC/fastqc/{sample}_fastqc.zip",sample= sample_list)
Align_bam= expand("Novogene/Chipseq_1.12.20/QC/bowtie2/{sample}.log",sample= sample_list_BAM)
MultiQC_listC= QC_html+Align_bam

rule all:
    input:
        expand("Novogene/Chipseq_1.12.20/QC/fastqc/{sample}_fastqc.html",sample=sample_list),
        expand("Novogene/Chipseq_1.12.20/Fastq/trimmed/{sample}_1_val_1.fq.gz",sample=sample_list_BAM),
        expand("Novogene/Chipseq_1.12.20/mapped/{sample}.bam",sample= sample_list_BAM),
        "Novogene/Chipseq_1.12.20/QC/multiqc/multiqc.html"

rule fastq:
    input:
        "Novogene/Chipseq_1.12.20/Fastq/{sample}.fq.gz"
    output:
        html="Novogene/Chipseq_1.12.20/QC/fastqc/{sample}_fastqc.html",
        zip="Novogene/Chipseq_1.12.20/QC/fastqc/{sample}_fastqc.zip"
    params:
        ""
    log:
        "Novogene/Chipseq_1.12.20/QC/fastqc/{sample}.log"
    threads:
        2
    wrapper:
        "0.49.0/bio/fastqc"

rule trim_galore:
    input:
        ["Novogene/Chipseq_1.12.20/Fastq/{sample}_1.fq.gz", "Novogene/Chipseq_1.12.20/Fastq/{sample}_2.fq.gz"]
    output:
        fastq1="Novogene/Chipseq_1.12.20/Fastq/trimmed/{sample}_1_val_1.fq.gz",
        fastq2="Novogene/Chipseq_1.12.20/Fastq/trimmed/{sample}_2_val_2.fq.gz",
        fastqc_dir=directory("Novogene/Chipseq_1.12.20/QC/fastqc/{sample}")
    params:
        extra="--length 1 -q 20"
    log:
         "Novogene/Chipseq_1.12.20/QC/cutadapt/{sample}.log"
    threads: 4 
    run:
        # Scale down the number of threads used. Since trim_galore uses ~ 4 times more threads then specified in pairend mode.
        trim_galore_threads = int(threads/4)               # Make fastqc out dir for later multiqc (as rule will complain if it doesn't exist)
        if not os.path.exists(output.fastqc_dir):
            os.mkdir(output.fastqc_dir)              # Figure out the output directory from the fastq1 file (as trim galore takes dir as a param)
        _outdir = os.path.dirname(output.fastq1)
        # Make sure that outdir for both fastqs is the same!
        assert _outdir == os.path.dirname(output.fastq2)        # Run trim_galore with less cores and specified fastqc out directory
        shell('trim_galore {params.extra} -j {threads} --gzip --paired --fastqc_args "--outdir {output.fastqc_dir}" -o {_outdir} {input} 2> {log}')
    #wrapper:
    #    "0.68.0/bio/cutadapt/pe"
    
    
rule bowtie2:
    input:
        sample=["Novogene/Chipseq_1.12.20/Fastq/trimmed/{sample}_1_val_1.fq.gz", "Novogene/Chipseq_1.12.20/Fastq/trimmed/{sample}_2_val_2.fq.gz"]
        #the real input would be inside a {sample} probably a problem
    output:
        temp("Novogene/Chipseq_1.12.20/mapped/{sample}.bam")
    log:
        "Novogene/Chipseq_1.12.20/QC/bowtie2/{sample}.log"
    params:
        index="/mnt/ife-ssd-home/shared-data/bowtie-indices/human/hg38",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "0.36.0/bio/bowtie2/align"  # "0.72.0/bio/bowtie2/align"


rule samtools_sort:
    input:
        "Novogene/Chipseq_1.12.20/mapped/{sample}.bam"
    output:
        "Novogene/Chipseq_1.12.20/mapped/{sample}.sorted.bam"
    params:
        "-m 4G"
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@.
    wrapper:
        "0.67.0/bio/samtools/sort"
    
rule samtools_index:
    input:
        "Novogene/Chipseq_1.12.20/mapped/{sample}.sorted.bam"
    output:
        "Novogene/Chipseq_1.12.20/mapped/{sample}.sorted.bam.bai"
    params:
   # optional params string
    wrapper:
        "0.67.0/bio/samtools/index"

rule multiqc:
    input:
        MultiQC_listC
    output:
        "Novogene/Chipseq_1.12.20/QC/multiqc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "Novogene/Chipseq_1.12.20/QC/multiqc/multiqc.log"
    wrapper:
        "0.67.0/bio/multiqc"

        

        # el nombre en bowtie debe cambiarse
        #Hecho cambio en bowtie Writing validated paired-end read 1 reads to A_11_FKDL202615954-1a_HJY3HDSXY_L2_1_val_1.fq.gz
        #Writing validated paired-end read 2 reads to A_11_FKDL202615954-1a_HJY3HDSXY_L2_2_val_2.fq.gz
        #4.3.21 10:39