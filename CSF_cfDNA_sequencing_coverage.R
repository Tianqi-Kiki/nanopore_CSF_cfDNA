# CSF_cfDNA_sequencing_coverage.R
# this file is meant to be used inside Rstudio
# this file takes a flat table of insert sizes for each read from a bam file, calculates and plots the sequencing coverage of all samples

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)

# read in length files

files_CSF <- Sys.glob('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0007_sample_summary_sup_methyl/align_length/*.align_length.txt')
length1 <- list ()

for (i in 1 : length(files_CSF)) {
  print(i)
  CSF <- files_CSF[i]
  sample <- basename(CSF)
  sample <- gsub(pattern = ".align_length.txt", replacement = "", x=sample)
  length1[[i]] <- fread(files_CSF[i]) %>% 
  filter(align_length != 'NA') %>% 
  mutate(sample = sample)  
}

df_CSF <- do.call(rbind.data.frame, length1)

# calculate sequencing coverage

seq_cov <- df_CSF %>% group_by(sample) %>% 
   summarise(sum_length=sum(align_length, na.rm = TRUE)) %>% 
   mutate(coverage = sum_length/3000000000) 

# normalize seq cov by input CSF volume

seq_cov_normalize <- seq_cov %>% select(c(1,3))%>% 
  mutate(seq_cov_per_ml=c(0.001, 0.033, 0.033, 7.051, 1.78, 0.446, 0.01, 
                          0.66, 0.045, 0.083, 0.006, 0.003, 0.003, 0.015,
                          0.005, 3.371, 4.945, 0.01, 0.043, 0.03, 0.014, 
                          0.065, 0.121, 0.143, 0.023, 0.023)) 

# prepare for plotting

  seq_cov_normalize <- seq_cov_normalize %>% 
mutate(type=ifelse(grepl('control',sample),'Healthy','Cancer')) %>% 
  mutate(depth = case_when(seq_cov_per_ml < 1 ~ 'low', 
                           seq_cov_per_ml > 1 ~ 'high'))

# plot seq cov per mL CSF

p1 <- ggplot(seq_cov_normalize, aes(x=reorder(sample, -seq_cov_per_ml), 
                                y=seq_cov_per_ml, fill=type)) +
  geom_bar(stat = 'identity',width=0.8) +
  facet_grid(depth ~ ., scale='free_y') +
  #scale_y_continuous(trans = 'log10') +
  labs(x = 'Sample',y= 'Sequencing Coverage 
       per mL CSF') +
  #geom_text(aes(label = signif(coverage,2)), nudge_y = 0.3, size=3) +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'),
        axis.text.x = element_text(angle = -45, hjust = 0.05),
        strip.text.y = element_blank())
p1

# compare seq cov cancer vs. control

p2 <- ggplot(seq_cov_normalize, aes(x = type, y = seq_cov_per_ml, 
                                      fill = type)) +
  geom_boxplot(width = 0.5) +
  #scale_y_continuous(limits = c(2,25), breaks = c(5,10,15,20,25)) +
  scale_y_continuous(trans = 'log10') +
  labs(x = NULL,y= 'Sequencing Coverage
       per mL CSF') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p2

pval_seq <- seq_cov_normalize %>% 
  summarise(p=t.test(seq_cov_per_ml[type=='Healthy'],
                     seq_cov_per_ml[type=='Cancer'])$p.value)

ggsave("seq_cov_perml_CSF_barplot.png", plot = p1, dpi = 300, units = "in", height = 6, width = 15)
ggsave("seq_cov_perml_CSF_cancervsctrl.png", plot = p2, dpi = 300, units = "in", height = 5, width = 8)




 
  


