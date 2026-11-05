import os
import pandas as pd
from itertools import combinations
import re
from os.path import join
import subprocess

workdir: config['Output_dir']
config_path = config["config_path"]
Output_dir = config['Output_dir']
scMIST = config['scMIST']
fq_path = config['fq_path']
ncbi_blast = config['ncbi_blast']
Kraken2Uniq_path = config['Kraken2Uniq_path']
kraken_database_path = config['kraken_library_path']
sample = config['sample'].strip().split(',')
umi = config['umi']
cb = config['cb']
r1=config['R1']
r2=config['R2']
bcumi_len = cb + umi


rule all:
    input:
##        expand('result/cellranger/{sample}.done',sample=sample),
##        expand("result/cellranger/{sample}/outs/possorted_genome_bam.bam",sample=sample),
        expand("result/unmap_human_fq/{sample}_unmap_human.bam",sample=sample),      
        expand('result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R1_001.fastq.gz',sample=sample),
        expand('result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R2_001.fastq.gz',sample=sample),
        expand('result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R1_001_filter.fastq.gz',sample=sample),
        expand('result/step1_re_unmap_score_6/{sample}.kraken.report.txt', sample=sample),
        expand('result/step1_re_unmap_score_6/{sample}.kraken.report.mpa.txt', sample=sample),
        expand('result/step2_re_unmap_score_6/{sample}_1.fa', sample=sample),
        expand('result/step2_re_unmap_score_6/{sample}_2.fa', sample=sample),
        expand('result/step3_re_unmap_score_6/{sample}.microbiome.output.txt', sample=sample),
        expand('result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_reads_map_lineage',sample=sample),
        expand('result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_species_1.fq',sample=sample),
        expand('result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_species_2.fq',sample=sample),
        expand('result/bwa_MHC_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',sample=sample),
        expand('result/bwa_MHC_0.95_gs_newlib_unmap_score_6/{sample}_2.fq',sample=sample),
        expand('result/bwa_rRNA_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',sample=sample),
        expand('result/bwa_rRNA_0.95_gs_newlib_unmap_score_6/{sample}_2.fq',sample=sample),   
        expand('result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt', sample=sample),
        expand('result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.mpa.txt', sample=sample),
        expand('result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_1.fa', sample=sample),
        expand('result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_2.fa', sample=sample),
        expand('result/step3_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.microbiome.output.txt', sample=sample),
        expand('result/bwa_black_0.95_gs_newlib_unmap_score_6/{sample}_reads_map_lineage', sample=sample),
        expand('result/bwa_black_0.95_gs_newlib_unmap_score_6/{sample}_bwa_black_species_id', sample=sample),
        expand('result/step4_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.sckmer.txt', sample=sample),
        expand('result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.barcode.kmer.hits.tsv', sample=sample),
        expand('result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.rpmm.tsv', sample=sample),
        'result/sample_denoising_v5_species_bwa_0.95_gs_newlib_unmap_score_6/sample_denoising.txt',
        expand('result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.cell_line_quantile_hits.tsv', sample=sample),
        expand('result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.cell_line_quantile_hits_taxa.tsv', sample=sample),
        expand('result/step7_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.counts.txt', sample=sample),
        'result/decontamn_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/decontamn.txt',
        expand('result/scMIST_cellranger_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_scMIST_cellranger.txt',sample=sample),


rule cellranger:
        input: join(fq_path,'{sample}_S1_L001_R1_001.fastq.gz'),
            join(fq_path,'{sample}_S1_L001_R2_001.fastq.gz')
        output: touch('result/cellranger/{sample}.done'),
                "result/cellranger/{sample}/outs/possorted_genome_bam.bam"
        params: sam='{sample}',
                outdir='result/cellranger/{sample}/'
        threads:  2
        shell:
                '''
#rm -r  {params.outdir}
/users/project_name/projects/MPXV/software/cellranger-7.2.0/cellranger count --id={params.sam} --transcriptome=/users/project_name/reference/Ensembl-GRCh38-2020A/refdata-gex-GRCh38-2020-A/ \
--fastqs={fq_path}/ \
--sample={params.sam} \
--output-dir={params.outdir}
                '''

rule unmap_bam2fq:
    input: "result/cellranger/{sample}/outs/possorted_genome_bam.bam"
    output:
        "result/unmap_human_fq/{sample}_unmap_human.bam",        
        "result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R1_001.fastq.gz",
        "result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R2_001.fastq.gz"
        #"result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/bamtofastq_S1_L000_R1_001.fastq.gz",
        #"result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/bamtofastq_S1_L000_R2_001.fastq.gz",
    params:
        outdir="result/unmap_human_fq/{sample}",
        sam='{sample}'
    threads: 2
    shell:
        """
        set -euo pipefail
        samtools view -b -f 4 {input} > {output[0]}
        rm -rf {params.outdir}
        /users/project_name/software/bamtofastq --nthreads=10 {output[0]} {params.outdir} --traceback
        src_dir=$(ls -d "{params.outdir}/{params.sam}_0_1_"* 2>/dev/null | head -n 1)
        dst_dir="{params.outdir}/{params.sam}_0_1_unknow_flowcell"

        if [[ -z "$src_dir" ]]; then
            echo "ERROR: bamtofastq produced no {params.sam}_0_1_* directory under {params.outdir}" >&2
            exit 1
        fi
        if [[ "$src_dir" != "$dst_dir" ]]; then
            mv "$src_dir" "$dst_dir"
        fi
        cat  {params.outdir}/{params.sam}_0_1_unknow_flowcell/bamtofastq_S1_L00*_R1_00*.fastq.gz> {output[1]}
        cat  {params.outdir}/{params.sam}_0_1_unknow_flowcell/bamtofastq_S1_L00*_R2_00*.fastq.gz> {output[2]}
        rm {params.outdir}/{params.sam}_0_1_unknow_flowcell/bamtofastq_S1_L00*_R1_00*.fastq.gz
        rm {params.outdir}/{params.sam}_0_1_unknow_flowcell/bamtofastq_S1_L00*_R2_00*.fastq.gz
        """

rule bcumi:
    input: "result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R1_001.fastq.gz"
    output: "result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R1_001_filter.fastq.gz"
    threads: 2
    params: bcumi_len=bcumi_len
    shell:
        '''
    /users/project_name/miniforge/envs/homer/bin/seqkit subseq -r  1:{params.bcumi_len} {input[0]}| pigz > {output[0]}
        '''

rule run_kraken:
    input: 
#        join(fq_path, '{sample}_{r1}.fastq.gz'),
#        join(fq_path, '{sample}_{r2}.fastq.gz')
#        fq1 = lambda wildcards: join(fq_path, f"{wildcards.sample}_{r1}.fastq.gz"),
#        fq2 = lambda wildcards: join(fq_path, f"{wildcards.sample}_{r2}.fastq.gz")
         fq1='result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R1_001_filter.fastq.gz',
         fq2='result/unmap_human_fq/{sample}/{sample}_0_1_unknow_flowcell/merge_bamtofastq_S1_L000_R2_001.fastq.gz'
    output: 
        'result/step1_re_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_re_unmap_score_6/{sample}.kraken.report.mpa.txt',
        'result/step1_re_unmap_score_6/{sample}.kraken.output.txt',
        'result/step1_re_unmap_score_6/{sample}_1.fq',
        'result/step1_re_unmap_score_6/{sample}_2.fq'
    params:
        sam='{sample}',
        outdir='result/step1_re_unmap_score_6/'
    threads: 2
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript {scMIST}/functions/run_kraken.r \
        --sample {params.sam} \
        --fq1 {input[0]} \
        --fq2 {input[1]} \
        --out_path {params.outdir} \
        --ncbi_blast_path {ncbi_blast} \
        --Kraken2Uniq_path {Kraken2Uniq_path} \
        --kraken_database_path {kraken_database_path} \
        --kreport2mpa_path {scMIST}/functions/kreport2mpa.py \
        --paired T
        '''

rule extract:
    input: 
        'result/step1_re_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_re_unmap_score_6/{sample}.kraken.report.mpa.txt',
        'result/step1_re_unmap_score_6/{sample}_1.fq',
        'result/step1_re_unmap_score_6/{sample}_2.fq'
    output: 
        'result/step2_re_unmap_score_6/{sample}_1.fa',
        'result/step2_re_unmap_score_6/{sample}_2.fa'
    params:
        sam='{sample}',
        outdir='result/step2_re_unmap_score_6/'
    threads: 2
    shell:
        '''
        Rscript {scMIST}/functions/extract_microbiome_reads.r \
        --sample_name {params.sam}_1 \
        --fq {input[2]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}


        Rscript {scMIST}/functions/extract_microbiome_reads.r \
        --sample_name {params.sam}_2 \
        --fq {input[3]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}
        '''

rule extract_microbiome_output:
    input: 
        'result/step1_re_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_re_unmap_score_6/{sample}.kraken.report.mpa.txt',
        'result/step1_re_unmap_score_6/{sample}.kraken.output.txt'
    output:
        'result/step3_re_unmap_score_6/{sample}.microbiome.output.txt'
    params:
        sam='{sample}',
        outdir='result/step3_re_unmap_score_6/'
    threads: 2
    shell:
        '''
        Rscript {scMIST}/functions/extract_microbiome_output.r --sample_name {params.sam} \
        --output_file {input[2]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}
        '''

#####################################
rule bwa:
    input:
        'result/step3_re_unmap_score_6/{sample}.microbiome.output.txt',
        'result/step1_re_unmap_score_6/{sample}.kraken.report.txt',
        'result/step2_re_unmap_score_6/{sample}_1.fa',
        'result/step2_re_unmap_score_6/{sample}_2.fa'
    output:
         'result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_reads_map_lineage',
         'result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_species_1.fq',
         'result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_species_2.fq'
    params:
        sam='{sample}',
        outdir='result/bwa_bwa_0.95_gs_newlib_unmap_score_6'
    threads: 2
    shell:
        '''
# /users/project_name/reference/kraken_database/database_ncbi_T2T_split_taxid_score_6/asshuman/bwa_library
# /users/project_name/micro/reference/kraken_database/database_ncbi_T2T_split_taxid/asshuman_score_6/bwa_library/

# revise ≥score 
Rscript {scMIST}/functions/reads_map_lineage_redatabase_1.r {input[0]} /users/project_name/reference/taxdum_matrix/taxonkit_lineage_matrix.txt \
         /users/project_name/reference/kraken_database/redatabase/kraken_overlap_databases.xlsx  {output[0]} 6
Rscript {scMIST}/functions/species_bwa_fa_v3.r --sample {params.sam} \
     --fa1 {input[2]} \
     --fa2 {input[3]} \
     --read_linea  {output[0]}  \
    --reference_fa  /users/project_name/reference/kraken_database/database_ncbi_T2T_split_taxid_score_6/asshuman/bwa_library \
    --condaenv /users/project_name/miniforge/envs/homer/bin \
    --work_path  {Output_dir} \
    --out_file {params.outdir} \
    --sfq1 {output[1]} \
    --sfq2 {output[2]} \
    --bwa_filter {scMIST}/functions/bwa_filter_v2.sh \
    --bc {cb} \
    --umi {umi} 
        '''

rule MHC:
    input:
         'result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_species_1.fq',
         'result/bwa_bwa_0.95_gs_newlib_unmap_score_6/{sample}_species_2.fq'
    output:
         'result/bwa_MHC_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',
         'result/bwa_MHC_0.95_gs_newlib_unmap_score_6/{sample}_2.fq'
    params:
        sam='{sample}',
        outdir='result/bwa_MHC_0.95_gs_newlib_unmap_score_6'
    threads: 2
    shell:
        '''
Rscript {scMIST}/functions/bwa_mhc_rrna.r --sample {params.sam} \
     --fq1 {input[0]} \
     --fq2 {input[1]} \
    --reference_fa  /users/project_name/reference/HLA/hla_gen.fasta \
    --condaenv /users/project_name/miniforge/envs/homer/bin \
    --work_path  {Output_dir} \
    --out_file {params.outdir} \
    --sfq1 {output[0]} \
    --sfq2 {output[1]}
        '''

rule rRNA:
    input:
         'result/bwa_MHC_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',
         'result/bwa_MHC_0.95_gs_newlib_unmap_score_6/{sample}_2.fq'
    output:
         'result/bwa_rRNA_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',
         'result/bwa_rRNA_0.95_gs_newlib_unmap_score_6/{sample}_2.fq'
    params:
        sam='{sample}',
        outdir='result/bwa_rRNA_0.95_gs_newlib_unmap_score_6'
    threads: 2
    shell:
        '''
Rscript {scMIST}/functions/bwa_mhc_rrna.r --sample {params.sam} \
     --fq1 {input[0]} \
     --fq2 {input[1]} \
    --reference_fa  /users/project_name/reference/rRNA/ncbi_human_rRNA/ncbi_human_18_28S_MT_rrna.fasta \
    --condaenv /users/project_name/miniforge/envs/homer/bin \
    --work_path  {Output_dir} \
    --out_file {params.outdir} \
    --sfq1 {output[0]} \
    --sfq2 {output[1]}
        '''

rule run_kraken_species:
    input: 
        'result/bwa_rRNA_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',
        'result/bwa_rRNA_0.95_gs_newlib_unmap_score_6/{sample}_2.fq'
    output: 
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.mpa.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.output.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_2.fq'
    params:
        sam='{sample}',
        outdir='result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/'
    threads: 2
    shell:
        '''
       # mkdir -p {params.outdir}
        Rscript {scMIST}/functions/run_kraken.r \
        --sample {params.sam} \
        --fq1 {input[0]} \
        --fq2 {input[1]} \
        --out_path {params.outdir} \
        --ncbi_blast_path {ncbi_blast} \
        --Kraken2Uniq_path {Kraken2Uniq_path} \
        --kraken_database_path {kraken_database_path} \
        --kreport2mpa_path {scMIST}/functions/kreport2mpa.py \
        --paired T
        '''
rule extract_species:
    input: 
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.mpa.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_1.fq',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_2.fq'
    output: 
        'result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_1.fa',
        'result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_2.fa'
    params:
        sam='{sample}',
        outdir='result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/'
    threads: 2
    shell:
        '''
        Rscript {scMIST}/functions/extract_microbiome_reads.r \
        --sample_name {params.sam}_1 \
        --fq {input[2]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}


        Rscript {scMIST}/functions/extract_microbiome_reads.r \
        --sample_name {params.sam}_2 \
        --fq {input[3]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}
        '''

rule extract_microbiome_output_species:
    input: 
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.mpa.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.output.txt'
    output:
        'result/step3_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.microbiome.output.txt'
    params:
        sam='{sample}',
        outdir='result/step3_species_bwa_0.95_gs_newlib_unmap_score_6/'
    threads: 2
    shell:
        '''
        Rscript {scMIST}/functions/extract_microbiome_output.r --sample_name {params.sam} \
        --output_file {input[2]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}
        '''


rule bwa_black:
    input:
        'result/step3_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.microbiome.output.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt',
        'result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_2.fa'
    output:
         'result/bwa_black_0.95_gs_newlib_unmap_score_6/{sample}_reads_map_lineage',
         'result/bwa_black_0.95_gs_newlib_unmap_score_6/{sample}_bwa_black_species_id'
    params:
        sam='{sample}',
        outdir='result/bwa_black_0.95_gs_newlib_unmap_score_6'
    threads: 2
    shell:
        '''
# revise ≥score 
Rscript {scMIST}/functions/reads_map_lineage_redatabase_1.r {input[0]} /users/project_name/reference/taxdum_matrix/taxonkit_lineage_matrix.txt \
         /users/project_name/reference/kraken_database/redatabase/kraken_overlap_databases.xlsx  {output[0]} 6
Rscript {scMIST}/functions/del_species_id_bwa_noreserve_droplet.r --sample {params.sam} \
     --fa2 {input[2]} \
     --read_linea  {output[0]}  \
    --reference_fa  /users/project_name/micro/reference/kraken_database/database_ncbi_T2T_split_taxid_score_6/asshuman/bwa_library \
    --condaenv /users/project_name/miniforge/envs/homer/bin \
    --work_path  {Output_dir} \
    --out_file {params.outdir} \
    --bwa_filter {scMIST}/functions/bwa_filter_v2.sh
        '''
########################################################
rule sckmer:
    input: 
        'result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_1.fa',
        'result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_2.fa',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.mpa.txt',
        'result/step3_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.microbiome.output.txt'
    output: 
        'result/step4_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.sckmer.txt'
    params:
        sam='{sample}',
        outdir='result/step4_species_bwa_0.95_gs_newlib_unmap_score_6/'
    threads: 2
    shell:
        '''
       # mkdir {params.outdir}
        Rscript {scMIST}/functions/sckmer_v6.r --sample_name={params.sam} \
        --fa1={input[0]} \
        --fa2={input[1]} \
        --microbiome_output_file={input[4]} \
        --kraken_report={input[2]} \
        --mpa_report={input[3]} \
        --out_path={params.outdir} \
        --cb_len={cb} \
        --umi_len={umi} \
        --ranks=c\\('G','S'\\)  \
        --host=9606 \
        '''

rule barcode_denoising:
    input: 
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt',
        'result/step4_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.sckmer.txt',
        'result/step1_re_unmap_score_6/{sample}.kraken.report.txt'
    output: 
        'result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.barcode.kmer.hits.tsv',
        'result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.rpmm.tsv'
    params: 
        outdir='result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/'
    threads: 2
    shell:
        '''
        #mkdir {params.outdir}
        Rscript {scMIST}/functions/barcode_denoising_v5.r {scMIST}/functions/read_kraken_reports.r {input[0]} {input[1]} {output[0]} {output[1]} {input[2]}
        '''

rule sample_denoising:
    input: expand('result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.rpmm.tsv',sample=sample) 
    output: 'result/sample_denoising_v5_species_bwa_0.95_gs_newlib_unmap_score_6/sample_denoising.txt'
    params:
        outdir='result/sample_denoising_v5_species_bwa_0.95_gs_newlib_unmap_score_6',
        sam=config['sample'],
        step5_path='result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6'
    threads: 2
    shell:
        '''
        #mkdir {params.outdir}
Rscript {scMIST}/functions/sample_denoising_v5.r {params.sam} {params.outdir} {params.step5_path}
        '''



rule cell_line_quantile_test:
    input: 
        'result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.barcode.kmer.hits.tsv',
        'result/step5_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.rpmm.tsv',
        'result/sample_denoising_v5_species_bwa_0.95_gs_newlib_unmap_score_6/sample_denoising.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt'
    output: 
        'result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.cell_line_quantile_hits.tsv',
        'result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.cell_line_quantile_hits_taxa.tsv'
    threads: 2
    shell:
        '''
        Rscript {scMIST}/functions/cell_line_quantile_test_v5.r {input[0]} {input[1]} '' {output[0]} {output[1]} {input[2]} {input[3]}
        '''

rule taxa_counts:
    input: 
        'result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_1.fa',
        'result/step2_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_2.fa',
        'result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.cell_line_quantile_hits_taxa.tsv',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.txt',
        'result/step1_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.kraken.report.mpa.txt'
    output: 
        'result/step7_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.counts.txt'
    params: 
        sam='{sample}',
        outdir='result/step7_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/'
    threads: 1
    shell:
        '''
        #mkdir {params.outdir}
        Rscript {scMIST}/functions/taxa_counts.r \
        --sample_name {params.sam} \
        --fa1 {input[0]} \
        --fa2 {input[1]} \
        --taxa {input[2]} \
        --kraken_report {input[3]} \
        --mpa_report {input[4]} \
        --out_path {params.outdir} \
        --cb_len {cb} \
        --umi_len {umi}
        '''

rule decontamn:
    input: expand('result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.cell_line_quantile_hits.tsv',sample=sample)
    output: 'result/decontamn_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/decontamn.txt'
    params:
        outdir='result/decontamn_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6',
        sam=config['sample'],
        step6_path='result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6',
        kr_path='result/step1_re_unmap_score_6/'
    threads: 2
    shell:
        '''
#        mkdir {params.outdir}
Rscript {scMIST}/functions/decontamn.r {params.sam} {output} 0.1 {params.step6_path} {params.kr_path}
        '''

rule scMIST_cellranger:
    input: 'result/decontamn_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/decontamn.txt',
        'result/cellranger/{sample}.done',
        'result/step7_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.counts.txt',
        'result/step6_sample_v5_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}.cell_line_quantile_hits.tsv',
        'result/cellranger/{sample}/outs/filtered_feature_bc_matrix'
    output: 'result/scMIST_cellranger_species_bwa_0.95_gs_newlib_unmap_score_6/{sample}_scMIST_cellranger.txt'
    threads: 2
    shell:
        '''
Rscript {scMIST}/functions/scMIST_cellranger_v1.r  {input[0]} {input[3]} {input[4]} {input[2]} {scMIST}/black_list.txt {scMIST}/mgnify_black.txt  {output}
        '''

