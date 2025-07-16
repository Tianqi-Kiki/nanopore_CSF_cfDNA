# CSF_cfDNA_5mC_5hmC_analysis.R
# this file is meant to be used inside Rstudio

# this script takes in promoter methylation (5mC) & intragenic hydroxymethylation (5hmC) % per CpG site with fdr-corrected q values, groups by gene level 
# and plots the promoter 5mC% or intragenic 5hmC% per gene (volcano plot) or per sample (point plot) for both cancer and control

# this file was generated using bedtools intersect with GENCODE v38:
# bedtools intersect -a gencode.v38.genes.sorted.bed -b P1.5mC.sorted.bedgraph -wa -wb -sorted | awk '{print($0 '\tP1')}' > P1.5mC.sorted.genes.intersect.bed
# P1.5mC.sorted.bedgraph is a bedgraph of the modified_bases.5mC.bed from dorado, filtered on column 5 > 0, and extracting only the chrom, start, end, and percent methylation columns
# multiple intersected bed files can be joined using cat, and fed into this script.

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyverse)

# plot differential prmoter 5mCG%  and intragenic 5hmCG% between cancer and control CSF of all genes in volcano plot
# read in merged bed files
# promoter 5mCG 

hyper_pro_5mc <- fread("/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/03_statistical_test/04_volcano_plot/02_samples/gene_summary/hyper.2kbpromoter_summary_5mCG.txt")
hypo_pro_5mc <- fread("/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/03_statistical_test/04_volcano_plot/02_samples/gene_summary/hypo.2kbpromoter_summary_5mCG.txt")

# intragenic 5hmCG

hyper_5hmc <- fread("/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/03_statistical_test/04_volcano_plot/02_samples/gene_summary/hyper.genebody_summary_5hmCG.txt")
hypo_5hmc <- fread("/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/03_statistical_test/04_volcano_plot/02_samples/gene_summary/hypo.genebody_summary_5hmCG.txt")

# gene analysis
# diff_methyl = 5mCG% (cancer ) - 5mCG% (control) per gene

# promoter 5mCG
# use gencodev38 gene names to filter protein_coding genes

gencode <- fread("/mnt/ix1/Projects/M086_221117_TGP_cfDNA_methylation/annotation/gencode.v38.genes_summary.txt")
gencode_proteincoding <- gencode %>% filter(V5=='protein_coding') 
protein_coding_genes <- unique(gencode_proteincoding$V6)

hyper_pro_5mc_mean <- hyper_pro_5mc %>% 
  select(c('V5', 'V6', 'V10', 'V11')) %>%
  rename(gene = V5, gene_id = V6, 
         diff_methyl = V10, qval = V11) %>%
  filter(gene %in% protein_coding_genes) 

hypo_pro_5mc_mean <- hypo_pro_5mc %>% 
  select(c('V5', 'V6', 'V10', 'V11')) %>%
  rename(gene = V5, gene_id = V6, 
         diff_methyl = V10, qval = V11) %>%
  filter(gene %in% protein_coding_genes) 

df_pro_5mc <- rbind(hyper_pro_5mc_mean, hypo_pro_5mc_mean)

# add legend to the plot

df_pro_5mc <- df_pro_5mc %>%
  mutate(group = case_when(
    qval < 0.01 & diff_methyl > 0 ~ "hypermethylated in cancer",
    qval < 0.01 & diff_methyl < 0 ~ "hypomethylated in cancer",
    TRUE ~ "Background"))

# add reported marker labels on the plot

df_marker_pro_5mc <- df_pro_5mc %>% 
  filter(gene %in% c('ATAD3A','CCNH','CD74','DYRK2','IL1R2',
                     'LIPH','MID1','MIF','MMP13','RAD54B',
                     'SETDB1','ZKSCAN1','COL6A6',
                     'EPB41L1','IL6R','PER3')) %>%
  filter(qval<0.01) %>% group_by(gene, diff_methyl,qval) %>% 
  summarize(n = n())

# plot the volcano plot

p1_2 <- ggplot(data=df_pro_5mc, aes(x=diff_methyl, y=-log10(1e-5+qval))) + 
  geom_point(aes(color=group),size=6, alpha=0.6) +
  scale_color_manual(name ='type',
                     values = c("hypermethylated in cancer" = "#CC79A7",       
                                "hypomethylated in cancer" = "#009E73")) +
  geom_text_repel(data = df_marker_pro_5mc, 
                  aes(x= diff_methyl, y=-log10(1e-5+qval), label = gene),
                  min.segment.length = 0,
                  max.overlaps=Inf,
                  size = 8 ) +
  labs(x='CpG promoter 5mCG difference per gene (%)') +
  theme_classic(30) +
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p1_2

# intragenic 5hmCG

hyper_5hmc_mean <- hyper_5hmc %>% 
  select(c('V5', 'V6', 'V10', 'V11')) %>%
  rename(gene_type = V5, gene = V6, 
         diff_methyl = V10, qval = V11) %>%
  filter(gene_type == 'protein_coding')

hypo_5hmc_mean <- hypo_5hmc %>% 
  select(c('V5', 'V6', 'V10', 'V11')) %>%
  rename(gene_type = V5, gene = V6, 
         diff_methyl = V10, qval = V11) %>%
  filter(gene_type == 'protein_coding')

df_5hmc <- rbind(hyper_5hmc_mean, hypo_5hmc_mean)

# add legend to the plot

df_5hmc <- df_5hmc %>%
  mutate(group = case_when(
    qval < 0.01 & diff_methyl > 0 ~ "hyper-hydroxymethylated in cancer",
    qval < 0.01 & diff_methyl < 0 ~ "hypo-hydroxymethylated in cancer",
    TRUE ~ "Background"))

# add reported marker labels on the plot

df_marker_5hmc <- df_5hmc %>% 
  filter(gene %in% c('ADCY9','HDAC9','IL21R','KCNQ1',
                     'MTDH','PBX1','SASH1','TPST1','TWIST2')) %>%
  filter(qval<0.01) %>% group_by(gene, diff_methyl,qval) %>% 
  summarize(n = n())

# plot the volcano plot

p1_2 <- ggplot(data=df_5hmc, aes(x=diff_methyl, y=-log10(1e-6 + qval))) + 
  geom_point(aes(color=group),size=6, alpha=0.6) +
  scale_color_manual(name ='type', values = c("hyper-hydroxymethylated in cancer" = "#CC79A7",
                                              "hypo-hydroxymethylated in cancer" = "#009E73")) +
  geom_text_repel(data = df_marker_5hmc,
                  aes(x = diff_methyl, y = -log10(1e-6 + qval), label = gene),
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size = 8) +
  labs(x='CpG intragenic 5hmCG difference per gene (%)') +
  theme_classic(30) +
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold')) 

p1_2

ggsave("5mCG_promoter_volcano.png", plot = p1_1, dpi = 300, units = "in", height = 12, width = 20)
ggsave("5hmCG_intragenic_volcano.png", plot = p1_2, dpi = 300, units = "in", height = 12, width = 20)


# plot differential prmoter 5mCG%  and intragenic 5hmCG% of marker genes per sample
# these marker genes were defined as sig genes (q<0.01) between cancer vs. control CSF samples, and reported in the literature as biomarkers for NSCLC
# e.g. a sig gene that was hypomethylated in our data and reported as an upregulated marker/overexpressed in NSCLC patients

# promoter 5mCG

p2_1 <- ggplot(marker_5mc_pro, aes(x= gene, y= mean_methyl, color = type)) +
  geom_jitter(size=6, alpha=0.6, width = 0.3, height = 2) +
  #scale_y_continuous(limits = c(0,100), breaks = c(0,50,100)) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,
                            10.5,11.5,12.5,13.5,14.5,15.5), 
             linetype="dotted",size=0.5) +
  labs(x = "Gene" ,y ="CpG promoter 5mCG 
       per sample (%)") +
  theme_classic(25) + 
  theme(axis.text = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = -45, hjust = 0.05)) 
p2_1

# intragenic 5hmCG

p2_2 <- ggplot(marker_5hmc, aes(x= gene, y= mean_methyl, color = type)) +
  geom_jitter(size=6, alpha=0.6, width = 0.3, height = 2) +
  #scale_y_continuous(limits = c(0,100), breaks = c(0,50,100)) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5), 
             linetype="dotted",size=0.5) +
  labs(x = "Gene" ,y ="CpG intragenic 5hmCG 
       per sample (%)") +
  theme_classic(30) + 
  theme(axis.text = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = -45, hjust = 0.05)) 
p2_2

ggsave("5mCG_promoter_marker_plot.png", plot = p2_1, dpi = 300, units = "in", height = 8, width = 20)
ggsave("5hmCG_intragenic_marker_plot.png", plot = p2_2, dpi = 300, units = "in", height = 10, width = 16)


# plot average genome-wide 5hmCG% between cancer and control CSF samples
# read in the merged bedmethyl files for 5hmCG
# these are bedgraphs containing 5hmC% before intersecting with GENCODE

cpg_5hmc <- fread("/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0007_sample_summary_sup_methyl/sorted_bedgraph/5hmCG/merged.genomewide.5hmCG.txt")
cpg_5hmc <- cpg_5hmc %>% mutate(type = ifelse(grepl('control',V5),'Healthy','Cancer'))

# calculate average genome-wide 5hmCG%

avg_5hmc <- cpg_5hmc %>% group_by(V5,type) %>% 
  summarise(mean_5hmCG = mean(V4)) %>% 
  mutate(m = signif(mean_5hmCG,2))

# plot 

p3_1 <- ggplot(avg_5hmc, aes(x= type, y= mean_5hmCG, fill = type)) +
  geom_boxplot(width=0.7) +    
  labs(x = NULL, y ="Average genome-wide 
       5hmCG (%)") +
  #scale_y_continuous(limits = c(0.8,9), breaks = c(1,3,5,7,9)) +
  theme_classic(25) + 
  theme(axis.text = element_text(face = "bold"), 
        axis.title = element_text(face = "bold")) 
p3_1

p3_2 <- ggplot(avg_5hmc, aes(x= m, fill = type)) +
  geom_bar(width=0.1, color ='black') +    
  labs(x = "Average genome-wide 5hmCG (%)", y ="Number of samples", fill = "Type") +
  #scale_x_continuous(limits = c(54,72), breaks = c(54,60,66,72)) +
  #scale_y_continuous(limits = c(0,2), breaks = c(0,1,2)) +
  facet_wrap(~ type, ncol = 1) +
  theme_classic(40) + 
  theme(strip.text = element_text(size = 40),
        axis.text = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"),
        panel.margin = unit(0, "lines")) 
p3_2

# caculate p value for the box plot

p_5hmc <- avg_5hmc %>% summarise(p=t.test(mean_5hmCG[type=='Healthy'], mean_5hmCG[type=='Cancer'])$p.value)

ggsave("genomewide_avg_5hmCG_boxplot.png", plot = p3_1, dpi = 300, units = "in", height = 6, width = 8)
ggsave('genomewide_avg_5hmCG_barplot.png', plot = p3_2, dpi = 300, units = 'in', height = 10, width = 20)



























