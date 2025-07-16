# cfDNA_enrichr_analysis.R
# this file is meant to be used inside Rstudio
# this script takes in gene lists, connects to the Enrichr databases and output graphs of enriched pathways and proteins

library(enrichR)
library(ggplot2)

# read in significant gene lists (q<0.01)
# intragenic 5hmCG 

gene_5hmc <- fread('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/03_statistical_test/03_graphs/01_all_samples/5hmCG_gb_sig_genes_extend.txt')
sig_genes1 <- gene_5hmc$gene

# promoter 5hmCG 

gene_5hmc_pro <- fread('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/03_statistical_test/03_graphs/01_all_samples/5hmCG_2kbpromoter_sig_genes.txt')
sig_genes2 <- gene_5hmc_pro$gene

# promoter 5mCG

gene_5mc_pro <- fread('/mnt/ix1/Projects/M085_221117_CSF_cfDNA_methylation/03_statistical_test/03_graphs/01_all_samples/5mCG_2kbpromoter_sig_genes_extend.txt')
sig_genes3 <- gene_5mc_pro$gene

listEnrichrSites()
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()

# plot
# intragenic 5hmCG

# bar plot using dbs NCI-Nature 2016

dbs1 <- 'NCI-Nature_2016'
enriched <- enrichr(sig_genes1, dbs1)
enriched[['ncinature']] <- enriched[['NCI-Nature_2016']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  ')) 
View(enriched[['ncinature']])

p1_1 <- plotEnrich(enriched[['ncinature']],showTerms = 10, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=0, size=4) + 
  theme_classic(20)+
  #scale_y_continuous(limits = c(0,2), breaks = c(0,1,2)) +
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched pathways', title = 'NCI-Nature 2016') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),           
        plot.title = element_text(hjust = 0.5, size = 15))
p1_1

# bar plot using dbs PPI Hub Proteins
dbs2 <- 'PPI_Hub_Proteins'
enriched <- enrichr(sig_genes1, dbs2)
enriched[['ppi_hub']] <- enriched[['PPI_Hub_Proteins']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  ')) 
View(enriched[['ppi_hub']])

p1_2 <- plotEnrich(enriched[['ppi_hub']],showTerms = 10, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=-0.1, size=5) + 
  theme_classic(20)+
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched proteins', title = 'PPI_Hub_Proteins') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),           
        plot.title = element_text(hjust = 0.5, size = 15)) 
p1_2

# bar plot using dbs Transcription_Factor_PPIs
dbs3 <- 'Transcription_Factor_PPIs'
enriched <- enrichr(sig_genes1, dbs3)
enriched[['ppi']] <- enriched[['Transcription_Factor_PPIs']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  ')) 
View(enriched[['ppi']])

p1_3 <- plotEnrich(enriched[['ppi']],showTerms = 10, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=-0.1, size=5) + 
  theme_classic(20)+
  #scale_y_continuous(limits = c(0,2), breaks = c(0,1,2)) +
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched proteins', title = 'Transcription_Factor_PPIs') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),           
        plot.title = element_text(hjust = 0.5, size = 15)) 
p1_3

# promoter 5hmCG

# bar plot using dbs Bioplanet 2019

dbs4 <- 'BioPlanet_2019'
enriched <- enrichr(sig_genes2, dbs4)
enriched[['bioplanet']] <- enriched[['BioPlanet_2019']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  '))  
View(enriched[['bioplanet']])

p2_1 <- plotEnrich(enriched[['bioplanet']],showTerms = 14, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=0, size=4) + 
  theme_classic(20)+
  #scale_y_continuous(limits = c(0,7), breaks = c(1,3,5,7)) +
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched pathways', title = 'BioPlanet 2019') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),           
        plot.title = element_text(hjust = 0.5, size = 15)) 
p2_1


# bar plot using dbs PPI Hub Proteins
dbs2 <- 'PPI_Hub_Proteins'
enriched <- enrichr(sig_genes2, dbs2)
enriched[['ppi_hub']] <- enriched[['PPI_Hub_Proteins']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  ')) #combine Term and e_p into 1 column 
View(enriched[['ppi_hub']])

p2_2 <- plotEnrich(enriched[['ppi_hub']],showTerms = 10, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=-0.1, size=5) + 
  theme_classic(20)+
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched proteins', title = 'PPI_Hub_Proteins') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),          
        plot.title = element_text(hjust = 0.5, size = 15)) 
p2_2

# bar plot using dbs Transcription_Factor_PPIs
dbs3 <- 'Transcription_Factor_PPIs'
enriched <- enrichr(sig_genes2, dbs3)
enriched[['ppi']] <- enriched[['Transcription_Factor_PPIs']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  ')) 
View(enriched[['ppi']])

p2_3 <- plotEnrich(enriched[['ppi']],showTerms = 10, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=-0.1, size=5) + 
  theme_classic(20)+
  #scale_y_continuous(limits = c(0,2), breaks = c(0,1,2)) +
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched proteins', title = 'Transcription_Factor_PPIs') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),           
        plot.title = element_text(hjust = 0.5, size = 15)) 
p2_3

# promoter 5mCG

# bar plot using dbs BioPlanet_2019

dbs4 <- 'BioPlanet_2019'
enriched <- enrichr(sig_genes3, dbs4)
enriched[['bioplanet']] <- enriched[['BioPlanet_2019']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  '))  
View(enriched[['bioplanet']])

p3_1 <- plotEnrich(enriched[['bioplanet']],showTerms = 10, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=0, size=4) + 
  theme_classic(20)+
  #scale_y_continuous(limits = c(0,7), breaks = c(1,3,5,7)) +
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched pathways', title = 'BioPlanet 2019') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),           
        plot.title = element_text(hjust = 0.5, size = 15))
p3_1

# bar plot using dbs PPI Hub Proteins
dbs2 <- 'PPI_Hub_Proteins'
enriched <- enrichr(sig_genes3, dbs2)
enriched[['ppi_hub']] <- enriched[['PPI_Hub_Proteins']] %>% 
  mutate(e_p=formatC(P.value, format = "e", digits = 2)) %>%  
  mutate(text=paste(Term,e_p,sep='  ')) 
View(enriched[['ppi_hub']])

p3_2 <- plotEnrich(enriched[['ppi_hub']],showTerms = 8, y = 'Count', orderBy = 'P.value' ) +
  geom_text(aes(label=text), color='black', y=0.01, hjust=-0.1, size=5) + 
  theme_classic(20)+
  scale_fill_gradient(low = 'lightsalmon',
                      high = 'lightskyblue') +
  labs(x = 'Enriched proteins', title = 'PPI_Hub_Proteins') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),          
        plot.title = element_text(hjust = 0.5, size = 15)) 
p3_2

ggsave("ncinature_pathway_5hmc_intragenic.png", plot = p1_1, dpi = 300, units = "in", height = 6, width = 12)
ggsave("bioplanet_pathway_5hmc_promoter.png", plot = p2_1, dpi = 300, units = "in", height = 7, width = 10)
ggsave("bioplanet_pathway_5mCG_promoter.png", plot = p3_1, dpi = 300, units = "in", height = 6, width = 10)

ggsave("ppi_hub_5mCG_promoter.png", plot = p3_2, dpi = 300, units = "in", height = 6, width = 10)
ggsave("ppi_hub_5hmCG_intragenic.png", plot = p1_2, dpi = 300, units = "in", height = 6, width = 10)
ggsave("ppi_transcription_5hmCG_intragenic.png", plot = p1_3, dpi = 300, units = "in", height = 6, width = 10)
ggsave("ppi_hub_5hmCG_promoter.png", plot = p2_2, dpi = 300, units = "in", height = 6, width = 10)
ggsave("ppi_transcription_5hmCG_promoter.png", plot = p2_3, dpi = 300, units = "in", height = 6, width = 10)





