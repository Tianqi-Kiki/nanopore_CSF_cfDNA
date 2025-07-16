# CSF_plasma_cfDNA_size_comparision.R
# this file is meant to be used inside Rstudio
# this file takes a flat table of insert sizes for each read from a bam file, and plots nucleosome ratios and size distribution

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(writexl)

# read in CSF length files

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

# calculate nucleosome ratios

df_CSF_sub_type <- df_CSF_sub %>% 
  mutate(type = case_when(align_length > 0 & align_length < 248.5 ~ 'mono',
                          align_length > 248.5 & align_length < 415 ~ 'di',
                          align_length > 415 & align_length < 600 ~ 'tri'))
df_CSF_sub_type <- na.omit(df_CSF_sub_type)

df_CSF_ratio <- df_CSF_sub_type %>% group_by(sample) %>%
  summarise(mono = sum(type == 'mono'), di = sum(type == 'di'), 
            tri = sum(type == 'tri'), md_ratio = mono/di, 
            mt_ratio = mono/tri,dt_ratio = di/tri) %>% 
  mutate(status=ifelse(grepl('control',sample),'Healthy','Cancer'))


# read in plasma length files
# cancer plasma

files_cancerplas <- Sys.glob('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0008_met_nsclc_plasma_sup_methyl/align_length/*.txt')
length2_subset <- list () 

# sample 100k aligned reads from each sample

for (i in 1 : length(files_cancerplas)) {
  print(i)
  cancer_plas <- files_cancerplas[i]
  sample <- basename(cancer_plas)
  sample <- gsub(pattern = ".align_length.txt", replacement = "", x=sample)
  set.seed(0)
  length2_subset[[i]] <- fread(files_cancerplas[i]) %>% 
    filter(align_length != 'NA') %>%
    sample_n(100000, replace = TRUE) %>%
    mutate(sample = sample) 
}

df_cancerplas <- do.call(rbind.data.frame, length2_subset) 

# control plasma group 1

files_pcontrolplas <- Sys.glob('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0009_control_plasma_precision_Billy/01_align_length/*length.txt')
length3_subset <- list () 

# sample 100k aligned reads from each sample

for (i in 1 : length(files_pcontrolplas)) {
  print(i)
  control_plas <- files_pcontrolplas[i]
  sample <- basename(control_plas)
  sample <- gsub(pattern = ".align_length.txt", replacement = "", x=sample)
  set.seed(0)
  length3_subset[[i]] <- fread(files_pcontrolplas[i]) %>% 
    filter(align_length != 'NA') %>%
    sample_n(100000, replace = TRUE) %>% 
    mutate(sample = sample) 
}

df_pre_controlplas <- do.call(rbind.data.frame, length3_subset) 

# control plasma group 2
files_dcontrolplas <- Sys.glob('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/00_sample/mp0010_control_plasma_discovery/align_length/*length.txt')
length4_subset <- list () 

# sample 100k aligned reads from each sample

for (i in 1 : length(files_dcontrolplas)) {
  print(i)
  control_plas <- files_dcontrolplas[i]
  sample <- basename(control_plas)
  sample <- gsub(pattern = ".align_length.txt", replacement = "", x=sample)
  set.seed(0)
  length4_subset[[i]] <- fread(files_dcontrolplas[i]) %>% 
    filter(align_length != 'NA') %>%
    sample_n(100000, replace = TRUE) %>% 
    mutate(sample = sample) 
}

df_dis_controlplas <- do.call(rbind.data.frame, length4_subset) 

# calculate nucleosome ratios
# cancer plasma

df_cancerplas_type <- df_cancerplas %>% 
  mutate(type = case_when(align_length > 0 & align_length < 248.5 ~ 'mono',
                          align_length > 248.5 & align_length < 415 ~ 'di',
                          align_length > 415 & align_length < 600 ~ 'tri'))
df_cancerplas_type <- na.omit(df_cancerplas_type)

df_cancerplas_ratio <- df_cancerplas_type %>% group_by(sample) %>%
  summarise(mono = sum(type == 'mono'), di = sum(type == 'di'), 
            tri = sum(type == 'tri'), md_ratio = mono/di, mt_ratio = mono/tri,
            dt_ratio = di/tri) 

# contorl plasma group 1

df_pre_controlplas_type <- df_pre_controlplas %>% 
  mutate(type = case_when(align_length > 0 & align_length < 248.5 ~ 'mono',
                          align_length > 248.5 & align_length < 415 ~ 'di',
                          align_length > 415 & align_length < 600 ~ 'tri'))
df_pre_controlplas_type <- na.omit(df_pre_controlplas_type) 

df_pre_controlplas_ratio <- df_pre_controlplas_type %>% group_by(sample) %>%
  summarise(mono = sum(type == 'mono'), di = sum(type == 'di'), 
            tri = sum(type == 'tri'), md_ratio = mono/di, mt_ratio = mono/tri,
            dt_ratio = di/tri)

# control plasma group 2

df_dis_controlplas_type <- df_dis_controlplas %>% 
  mutate(type = case_when(align_length > 0 & align_length < 248.5 ~ 'mono',
                          align_length > 248.5 & align_length < 415 ~ 'di',
                          align_length > 415 & align_length < 600 ~ 'tri'))
df_dis_controlplas_type <- na.omit(df_dis_controlplas_type) 

df_dis_controlplas_ratio <- df_dis_controlplas_type %>% group_by(sample) %>%
  summarise(mono = sum(type == 'mono'), di = sum(type == 'di'), 
            tri = sum(type == 'tri'), md_ratio = mono/di, mt_ratio = mono/tri,
            dt_ratio = di/tri) 

# combine control and cancer plasma dataset

df_plas_ratio <- rbind(df_cancerplas_ratio, df_pre_controlplas_ratio, df_dis_controlplas_type) %>%
mutate(status=ifelse(grepl('NSCLC',sample),'Cancer','Healthy'))


# plot mono-/di-, mono-/tri-, di-/tri-nucleosome ratios in CSF

p1_1 <- ggplot(df_CSF_ratio, aes(x = status, y = md_ratio, fill = status)) +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(1.5,8), breaks = c(2,4,6,8)) +
  labs(x = NULL,y= 'Mononucleosome to
       dinucleosome ratio') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p1_1

p1_2 <- ggplot(df_CSF_ratio, aes(x = status, y = mt_ratio, fill = status)) +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(2,25), breaks = c(5,10,15,20,25)) +
  labs(x = NULL,y= 'Mononucleosome to
       trinucleosome ratio') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p1_2

p1_3 <- ggplot(df_CSF_ratio, aes(x = status, y = dt_ratio, fill = status)) +
  geom_boxplot(width = 0.5) +
  scale_y_continuous(limits = c(1,5), breaks = c(1,3,5)) +
  #geom_hline(yintercept = 4, linetype="dotted",size=0.5) +
  labs(x = NULL,y= 'Dinucleosome to
       trinucleosome ratio') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p1_3

# calculate p values

p_md_CSF <- df_CSF_ratio %>% summarise(p=t.test(md_ratio[status=='Healthy'],md_ratio[status=='Cancer'])$p.value)
p_mt_CSF <- df_CSF_ratio %>% summarise(p=t.test(mt_ratio[status=='Healthy'],mt_ratio[status=='Cancer'])$p.value)
p_dt_CSF <- df_CSF_ratio %>% summarise(p=t.test(dt_ratio[status=='Healthy'],dt_ratio[status=='Cancer'])$p.value)


# plot mono-/di-, mono-/tri-, di-/tri-nucleosome ratios in plasma

p2_1 <- ggplot(df_plas_ratio, aes(x = status, y = md_ratio, fill = status)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#C77CFF','#7CAE00')) +
  #scale_y_continuous(limits = c(2,25), breaks = c(5,10,15,20,25)) +
  labs(x = NULL,y= 'Mononucleosome to
       dinucleosome ratio') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p2_1

p2_2 <- ggplot(df_plas_ratio, aes(x = status, y = mt_ratio, fill = status)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#C77CFF','#7CAE00')) +
  scale_y_continuous(limits = c(0,200), breaks = c(50,100,150,200)) +
  labs(x = NULL,y= 'Mononucleosome to
       trinucleosome ratio') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p2_2

p2_3 <- ggplot(df_plas_ratio, aes(x = status, y = dt_ratio, fill = status)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c('#C77CFF','#7CAE00')) +
  #scale_y_continuous(limits = c(2,25), breaks = c(5,10,15,20,25)) +
  labs(x = NULL,y= 'Dinucleosome to
       trinucleosome ratio') +
  theme_classic(20) + 
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'))
p2_3

# calculate p values

p_md_plas <- df_plas_ratio %>% summarise(p=t.test(md_ratio[status=='Healthy'], md_ratio[status=='Cancer'])$p.value)
p_mt_plas <- df_plas_ratio %>% summarise(p=t.test(mt_ratio[status=='Healthy'], mt_ratio[status=='Cancer'])$p.value)
p_dt_plas <- df_plas_ratio %>% summarise(p=t.test(dt_ratio[status=='Healthy'], dt_ratio[status=='Cancer'])$p.value)


ggsave('monodi_ratio_CSF_boxplot.png', plot = p1_1, dpi = 300, units = 'in', height = 5, width = 7)
ggsave('monotri_ratio_CSF_boxplot.png', plot = p1_2, dpi = 300, units = 'in', height = 5, width = 7)
ggsave('ditri_ratio_CSF_boxplot.png', plot = p1_3, dpi = 300, units = 'in', height = 5, width = 7)
ggsave('monodi_ratio_plas_boxplot.png', plot = p2_1, dpi = 300, units = 'in', height = 5, width = 7)
ggsave('monotri_ratio_plas_boxplot.png', plot = p2_2, dpi = 300, units = 'in', height = 5, width = 7)
ggsave('ditri_ratio_plas_boxplot.png', plot = p2_3, dpi = 300, units = 'in', height = 5, width = 7)

# plot the cfDNA fragment size distribution
# CSF

p3 <-ggplot(df_CSF_sub,aes(x= align_length, fill = sample)) +
   geom_density(show.legend = FALSE) +
   scale_x_continuous(limits = c(0,800), breaks = c(0,100,300,500)) +
   scale_y_continuous(limits = c(0,0.015), breaks = c(0.0000, 0.0075, 0.0150)) +
  labs(x= "Insert size (bp)", y="Density") +
  geom_vline(xintercept = c(167,330,500), linetype="dotted",size=1) +       
 facet_wrap(~ sample, ncol = 4) +                                         
  theme_classic(50) +                                                     
  theme(strip.text = element_text(size = 45), 
        strip.background = element_rect(size = 2)) +  
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'),          
        panel.spacing = unit(0, "lines"))                    
p3 

# cancer plasma

p4 <-ggplot(df_cancerplas,aes(x= align_length, fill = sample)) +
  geom_density(show.legend = FALSE) +
  scale_x_continuous(limits = c(0,800), breaks = c(0,100,300,500)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02)) +
  labs(x= "Insert size (bp)", y="Density") +
  geom_vline(xintercept = c(167,330,500), linetype="dotted",size=1) +       
  facet_wrap(~ sample, ncol = 4) +                                         
  theme_classic(50) +                                                     
  theme(strip.text = element_text(size = 45), 
        strip.background = element_rect(size = 2)) +  
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'),          
        panel.spacing = unit(0, "lines"))                     
p4

# control plasma group 1 (select 10 samples to plot)

df_10control_pre_plas <- df_pre_controlplas %>% 
  filter(sample %in% c('P10083_30400.D01_control','P10084_30401.D01_control','P10085_30402.D01_control',
           'P10619_27476.D01_control','P10620_27477.D01_control','P10621_27478.D01_control',
           'P10088_30405.D01_control','P10099_30416.D01_control','P10121_30438.D01_control',
           'P10123_30440.D01_control'))

p5 <-ggplot(df_10control_pre_plas,aes(x= align_length, fill = sample)) +
  geom_density(show.legend = FALSE) +
  scale_x_continuous(limits = c(0,800), breaks = c(0,100,300,500)) +
  labs(x= "Insert size (bp)", y="Density") +
  geom_vline(xintercept = c(167,330,500), linetype="dotted",size=1) +       
  facet_wrap(~ sample, ncol = 4) +                                         
  theme_classic(50) +                                                     
  theme(strip.text = element_text(size = 45), 
        strip.background = element_rect(size = 2)) +  
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'),          
        panel.spacing = unit(0, "lines"))                   
p5

# control plasma group 2

p6 <-ggplot(df_dis_controlplas,aes(x= align_length, fill = sample)) +
  geom_density(show.legend = FALSE) +
  scale_x_continuous(limits = c(0,800), breaks = c(0,100,300,500)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02)) +
  labs(x= "Insert size (bp)", y="Density") +
  geom_vline(xintercept = c(167,330,500), linetype="dotted",size=1) +       
  facet_wrap(~ sample, ncol = 4) +                                         
  theme_classic(50) +                                                     
  theme(strip.text = element_text(size = 45), 
        strip.background = element_rect(size = 2)) +  
  theme(axis.text = element_text(face = 'bold'), 
        axis.title = element_text(face = 'bold'),           
        panel.spacing = unit(0, "lines"))                    
p6

ggsave('CSF_cfDNA_fragment_distribution.png', plot = p3, dpi = 300, units = 'in', height = 36, width = 36)
ggsave('cancerplasma_cfDNA_fragment_distribution.png', plot = p4, dpi = 300, units = 'in', height = 25, width = 36)
ggsave('precision_controlplasma_cfDNA_fragment_distribution.png', plot = p5, dpi = 300, units = 'in', height = 18, width = 40)
ggsave('discovery_controlplasma_cfDNA_fragment_distribution.png', plot = p6, dpi = 300, units = 'in', height = 22, width = 38)















