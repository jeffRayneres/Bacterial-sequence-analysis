#!/bin/bash

threads=1

algorithm_abr="an"

## Overlap forward and reverse reads
mothur "#make.contigs(file=stability.files, processors=${threads})" > step_1_make_contigs.log

## Filter ambiguous sequences and sequences longer than 550
mothur "#screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=550, processors=$threads)" > step_2_screen_seqs.log

## Dereplication
mothur "#unique.seqs(fasta=stability.trim.contigs.good.fasta)" > step_3_unique_seqs.log

## Calculate the frequencies of each sequence in each sample
mothur "#count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)" > step_4_count_seqs.log

## Align sequences
mothur "#align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v3.fasta, processors=$threads)" > step_5_align_seqs.log

## Run summary.seqs and use output file stability.trim.contigs.good.unique.summary to determine the start and end positions of alignments
mothur "#summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, processors=$threads)" > step_6_summary_seqs.log

## Filter sequences that fail to overlap from start_pos to end_pos
mothur "#screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=1829, end=4088, maxhomop=8, processors=$threads)" > step_7_screen_seqs.log

mothur "#filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=$threads)" > step_8_filter_seqs.log

## Dereplicate aligned sequences
mothur "#unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)" > step_9_unique_seqs.log
mothur "#unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, name=stability.trim.contigs.good.names)" > step_9b_unique_seqs_names.log

## Pre-cluster the aligned sequences
mothur "#pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=5, processors=$threads)" > step_10_precluster.log
cat *.map > single.map

## Identify chimeric sequences and remove them
mothur "#chimera.uchime(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t, processors=$threads); remove.seqs(fasta=current, accnos=current)" > step_11_uchime.log

## Assign taxonomy to clean sequences
mothur "#classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, reference=silva.v3.fasta, taxonomy=silva.v3.tax, cutoff=80)" > step_12_classify_seqs.log

## Remove unexpected lineages
mothur "#remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.taxonomy, taxon=Chloroplast-Mitochondria-Archaea-Eukaryota)" > step_13_remove_lineage.log

# Cluster
mothur "#cluster.split(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.15, processors=$threads)" > step_14_cluster_split.log

## Compute number of sequences in each OTU
mothur "#make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, label=0.03)" > step_15_make_shared.log

mothur "#classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.v3.wang.pick.taxonomy, label=0.03)" > step_pre_biom

## Create .biom file
mothur "#make.biom(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, constaxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy)" > step_16_make_biom.log

## Create representative FASTA file for OTUs
mothur "#get.oturep(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, group=stability.contigs.good.groups, name=stability.trim.contigs.good.names, method=abundance, label=0.03)" > step_17_get_oturep.log

## Get OTU map
mothur "#get.otulist(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.${algorithm_abr}.unique_list.list, label=0.03)" > step_18_get_otulist.log

## Create the tree
mothur "#tree.shared(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared)" > step_19_tree_shared
