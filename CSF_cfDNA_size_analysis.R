# CSF_cfDNA_size_analysis.R
# this file is meant to be used inside Rstudio
# this file takes a flat table of insert sizes for each read from a bam file, and plots nucleosome ratios and size distribution

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(writexl)

# read in length files

files_CSF <- Sys.glob('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0007_sample_summary_sup_methyl/align_length/*.align_length.txt')
length1_subset <- list () 

# sample 100k aligned reads from each sample

for (i in 1 : length(files_CSF)) {
  print(i)
  CSF <- files_CSF[i]
  sample <- basename(CSF)
  sample <- gsub(pattern = ".align_length.txt", replacement = "", x=sample)
  set.seed(0)
  length1_subset[[i]] <- fread(files_CSF[i]) %>% 
  filter(align_length != 'NA') %>%
  sample_n(100000, replace = TRUE) %>%
  mutate(sample = sample)
}

df_CSF_sub <- do.call(rbind.data.frame, length1_subset)


# classify mono-, di- and tri-nucleosomes
# peak size: mono = 167 bp, di = 330 bp, tri = 500 bp
# mono-di cut-off = 1/2(167+330)= 248.5 bp
# di-tri cut-off = 1/2(330+500)= 415 bp

df_CSF_sub_type <- df_CSF_sub %>% 
  mutate(type = case_when(align_length > 0 & align_length < 248.5 ~ 'mono',
                          align_length > 248.5 & align_length < 415 ~ 'di',
                          align_length > 415 & align_length < 600 ~ 'tri'))
df_CSF_sub_type <- na.omit(df_CSF_sub_type)


# calculate mono/di, mono/tri, di/trinucleosome ratios

df_CSF_ratio <- df_CSF_sub_type %>% group_by(sample) %>%
  summarise(mono = sum(type == 'mono'), di = sum(type == 'di'), 
            tri = sum(type == 'tri'), md_ratio = mono/di, 
            mt_ratio = mono/tri,dt_ratio = di/tri) %>% 
  mutate(status=ifelse(grepl('control',sample),'Healthy','Cancer'))


# check mono- and tri-nucleosome enrichment in cancer vs. control
# calculate the ratio of mono-/all nucleosomes and tri-/all nucleosomes 

df_enrich_ratio <- df_CSF_ratio %>%
  group_by(sample) %>%
  summarise(all = mono + di + tri, 
            mono_to_all = mono/all, tri_to_all = tri/all) %>%
  mutate(type=ifelse(grepl('control',sample),'Healthy','Cancer'))
  
# calculate p value & enrich factor

mono_p <- df_enrich_ratio %>% 
  summarise(p=t.test(mono_to_all[type=='Healthy'],
                     mono_to_all[type=='Cancer'])$p.value)

mono_c <- mean(df_enrich_ratio$mono_to_all[df_enrich_ratio$type == 'Cancer'])
mono_h <- mean(df_enrich_ratio$mono_to_all[df_enrich_ratio$type == 'Healthy'])

mono_enrich <-mono_c/mono_h  

tri_p <- df_enrich_ratio %>% 
  summarise(p=t.test(tri_to_all[type=='Healthy'],
                     tri_to_all[type=='Cancer'])$p.value)

tri_c <- mean(df_enrich_ratio$tri_to_all[df_enrich_ratio$type == 'Cancer'])
tri_h <- mean(df_enrich_ratio$tri_to_all[df_enrich_ratio$type == 'Healthy'])

tri_enrich <-tri_h/tri_c   


# plot

p1_1 <- ggplot(df_enrich_ratio, 
               aes(x = type, y = mono_to_all, fill = type)) +
  geom_boxplot(width = 0.5) +
  #scale_y_continuous(limits = c(0.5,0.85), breaks = c(0.5,0.65,0.8)) +
  labs(x = NULL,y= 'Mononucleosome 
       over all reads') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p1_1

p1_2 <- ggplot(df_enrich_ratio, 
               aes(x = type, y = tri_to_all, fill = type)) +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(0.03,0.2), breaks = c(0.05,0.10,0.15,0.20)) +
  labs(x = NULL,y= 'Trinucleosome 
       over all reads') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p1_2

# pick 1 control and 1 cancer samples to show the nucleosome size distributions

length_641 <- fread("/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0007_sample_summary_sup_methyl/align_length/TP641_control.align_length.txt")   
length_641_sub <- length_641 %>% sample_n(100000)        
length_641_sub$sample <- rep('Healthy_1',dim(length_641_sub)[1])
head(length_641_sub) 

length_432 <- fread("/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0007_sample_summary_sup_methyl/align_length/TP432_cancer.align_length.txt")   
length_432_sub <- length_432 %>% sample_n(100000)        
length_432_sub$sample <- rep('Cancer_1',dim(length_432_sub)[1])
head(length_432_sub) 

df_size <- rbind(length_641_sub, length_432_sub)

# draw plot

p2 <-ggplot(df_size,aes(x= align_length, fill = sample)) +
  geom_density(show.legend = FALSE) +
  scale_x_continuous(limits = c(0,800), breaks = c(0,100,300,500)) +
  scale_y_continuous(limits = c(0,0.015), breaks = c(0.0000, 0.0075, 0.0150)) +
  labs(x= "Insert size (bp)", y="Density") +
  geom_vline(xintercept = c(167,330,500), linetype="dotted",size=1) +       
  facet_wrap(~ sample, ncol = 1) +                                         
  theme_classic(45) +                                                     
  theme(strip.text = element_text(size = 45), 
        strip.background = element_rect(size = 2)) +  
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'),           
        panel.spacing = unit(0, "lines"))                    
p2

ggsave('mono_in_all_enrich_boxplot.png', plot = p1_1, dpi = 300, units = 'in', height = 4, width = 7)
ggsave('tri_in_all_enrich_boxplot.png', plot = p1_2, dpi = 300, units = 'in', height = 4, width = 7)
ggsave('cfDNA_size_exmaple.png', plot = p2, dpi = 300, units = 'in', height = 14, width = 12)









