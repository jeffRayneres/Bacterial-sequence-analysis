#!/bin/bash
# Based on scripts: http://www.mothur.org/wiki/MiSeq_SOP and https://github.com/ekopylova/OTU-clustering.git

#SBATCH -J An_m_J_R
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16000M
#SBATCH --time=48:00:00
#SBATCH --constraint=ib
#SBATCH --mail-user=jeff@dcc.ufmg.br
#SBATCH --mail-type=ALL

mkdir /var/tmp/$SLURM_JOB_ID
echo "Criacao de diretorio temp${SLURM_JOB_ID}" > criacao_diretorio_temp.txt

srun hostname | sort

threads=1

algorithm_abr="an"

## Overlap forward and reverse reads
## Output File Names: 
##  + stability.trim.contigs.fasta
##  + stability.contigs.qual
##  + stability.contigs.report
##  + stability.scrap.contigs.fasta
##  + stability.scrap.contigs.qual
##  + stability.contigs.groups
srun mothur "#make.contigs(file=stability.files, processors=${threads})" > step_1_make_contigs.log

## Filter ambiguous sequences and sequences longer than 550
## Output File Names: 
##  + stability.trim.contigs.good.fasta
##  + stability.trim.contigs.bad.accnos
##  + stability.contigs.good.groups
srun mothur "#screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=550, processors=$threads)" > step_2_screen_seqs.log

## Dereplication
## Output File Names: 
##  + stability.trim.contigs.good.names
##  + stability.trim.contigs.good.unique.fasta
srun mothur "#unique.seqs(fasta=stability.trim.contigs.good.fasta)" > step_3_unique_seqs.log

## Calculate the frequencies of each sequence in each sample
## Output File Names: 
##  + stability.trim.contigs.good.count_table
srun mothur "#count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)" > step_4_count_seqs.log

## Align sequences
## Output File Names:
##  + stability.trim.contigs.good.unique.align
##  + stability.trim.contigs.good.unique.align.report
##  + stability.trim.contigs.good.unique.flip.accnos
srun mothur "#align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=$SLURM_SUBMIT_DIR/silva.v3.fasta, processors=$threads)" > step_5_align_seqs.log

## Run summary.seqs and use output file stability.trim.contigs.good.unique.summary to
## determine the start and end positions of alignments
srun mothur "#summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, processors=$threads)" > step_6_summary_seqs.log

## Filter sequences that fail to overlap from start_pos to end_pos
## Output File Names: 
##  + stability.trim.contigs.good.unique.good.summary
##  + stability.trim.contigs.good.unique.good.align
##  + stability.trim.contigs.good.unique.bad.accnos
##  + stability.trim.contigs.good.good.count_table
srun mothur "#screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=6388, end=25316, maxhomop=8, processors=$threads)" > step_7_screen_seqs.log

## Output File Names: 
##  + stability.filter
##  + stability.trim.contigs.good.unique.good.filter.fasta
srun mothur "#filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=$threads)" > step_8_filter_seqs.log

## Dereplicate aligned sequences
## Output File Names: 
## + stability.trim.contigs.good.unique.good.filter.count_table
## + stability.trim.contigs.good.unique.good.filter.unique.fasta
srun mothur "#unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)" > step_9_unique_seqs.log
srun mothur "#unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, name=stability.trim.contigs.good.names)" > step_9b_unique_seqs_names.log

## Pre-cluster the aligned sequences
## Output File Names: 
## + stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta
## + stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table
## + stability.trim.contigs.good.unique.good.filter.unique.precluster.s1.map
srun mothur "#pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=5, processors=$threads)" > step_10_precluster.log
#cat *.map > single.map

## Identify chimeric sequences and remove them
## Output File Names: 
## + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
srun mothur "#chimera.uchime(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t, processors=$threads); remove.seqs(fasta=current, accnos=current)" > step_11_uchime.log

## Assign taxonomy to clean sequences
## Output File Names: 
## + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.taxonomy
## + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.tax.summary
srun mothur "#classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, reference=$SLURM_SUBMIT_DIR/silva.v3.fasta, taxonomy=$SLURM_SUBMIT_DIR/silva.v3.tax, cutoff=80)" > step_12_classify_seqs.log

## Remove unexpected lineages
## Output File Names: 
##  + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.pick.taxonomy
##  + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
##  + stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table
srun mothur "#remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.taxonomy, taxon=Chloroplast-Mitochondria-Archaea-Eukaryota)" > step_13_remove_lineage.log

# Cluster
# Output File Names: 
#  + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list
srun mothur "#cluster.split(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.15, processors=$threads)" > step_14_cluster_split.log

# Compute number of sequences in each OTU
# Output File Names: 
#  + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.shared
#  + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.s1.rabund
srun mothur "#make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, label=0.03)" > step_15_make_shared.log

srun mothur "#classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.pick.taxonomy, label=0.03)" > step_intermediate_biom

# Create BIOM table
# Output File Names: 
#   + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.biom
srun mothur "#make.biom(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, constaxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy)" > step_16_make_biom.log

## Create representative FASTA file for OTUs
## Output File Names:
##   + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.rep.names
##   + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.0.03.rep.fasta
srun mothur "#get.oturep(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, group=stability.contigs.good.groups, name=stability.trim.contigs.good.names, method=abundance, label=0.03)" > step_17_get_oturep.log

## Get OTU map
## Output File Names:
##  + stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.otu
srun mothur "#get.otulist(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, label=0.03)" > step_18_get_otulist.log

srun mothur "#tree.shared(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared)" > step_19_tree_shared

