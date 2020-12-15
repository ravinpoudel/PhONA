# installation
library(devtools)

#devtools::install_github("ravinpoudel/PhONA", build_vignettes = TRUE, force = TRUE, auth = "5c5decd148f0378dcb762e7b14c3d1508ef49ba2")

# Browse Vignettes to see the output as below.

#browseVignettes("PhONA")

library(tufte)
library(PhONA)


# GLOBAL VARIABLE
ITERS=2
OTU_OTU_PVALUE = 0.05
OTU_OTU_RVALUE = 0.5
OTU_PHENOTYPE_PVALUE = 1


###### Load the data

otu_data_fungi <- read.mothur.shared("ms_run/data/Fungi1415_trim.contigs.trim.unique.precluster.pick.pick.subsample.nn.unique_list.shared")
rownames(otu_data_fungi) <- paste0("F", rownames(otu_data_fungi))

# upload meta data file
bigmetadata <- read.csv("ms_run/data/block_diversity_selected_tunnel_ww.csv", header = T, row.names = 1, stringsAsFactors = FALSE)
rownames(bigmetadata) <- paste0("F", rownames(bigmetadata))

### select otu_data based on metadata
otu_data <- otu_data_fungi[rownames(bigmetadata), ]
dim(otu_data)
#rownames(otu_data)


otu_data = na.omit(otu_data)

dim(otu_data)


# now select metadata based on od
meta_data <- bigmetadata[rownames(otu_data), , drop = FALSE]

all.equal(row.names(meta_data), row.names(otu_data))


# upload taxanomy file
tax_fungi <- read.mothur.taxonomy("ms_run/data/Fungi1415_trim.contigs.trim.unique.precluster.pick.pick.subsample.nn.unique_list.0.03.cons.taxonomy")



# checking if taxonomy file and count files have same otus
all.equal(colnames(otu_data), row.names(tax_fungi))

#
########## read in phenome data/ Yield data
# read in data
pheno_data <- read.csv("ms_run/data/tomato_yield.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
rownames(pheno_data) <- paste0("F", rownames(pheno_data))

hist(pheno_data$Marketable)
pheno_data_mean <- aggregate(pheno_data$Marketable, by=list(Rootstock=pheno_data$Rootstock), FUN=mean)


pheno_data_sel <- pheno_data %>%
  select(Rootstock, Compartment, Marketable, Study_site, Year, Sample_name)

meta_withpheo = inner_join(meta_data, pheno_data_sel)
rownames(meta_withpheo)<- rownames(meta_data)

# create a phyloseq object~ can combine count, metadata, taxonomy, and phylogentic tree.
tax.mat <- tax_table(as.matrix(tax_fungi))
otu.mat = otu_table(t(otu_data), taxa_are_rows = TRUE)
sample.mat <- sample_data(meta_withpheo)
physeq = phyloseq(otu.mat, tax.mat, sample.mat)


## Assign color to the taxa on the whole phyloseq object so that the same color is assigned for a taxon across treatments
physeq = taxacolor(phyobj = physeq, coloredby = "Phylum")

physeq


#### All the combinations of treatment factor

region <- unique(as.character(data.frame(sample_data(physeq))$Compartment))
region
treatment <- unique(as.character(data.frame(sample_data(physeq))$Rootstock))
treatment



# Endosphere
ng_endo <- subset_samples(physeq, c(Compartment =="Endosphere" & Rootstock=="Nongraft")) %>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
ng_endo_sparcc = otu_table(ng_endo)
write.table(data.frame(OTU_id = rownames(ng_endo_sparcc), ng_endo_sparcc), file = "ms_run/data/sparcc_data/ng_endo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)



sg_endo <- subset_samples(physeq, c(Compartment =="Endosphere" & Rootstock=="Selfgraft")) %>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
sg_endo_sparcc = otu_table(sg_endo)
write.table(data.frame(OTU_id = rownames(sg_endo_sparcc), sg_endo_sparcc), file = "ms_run/data/sparcc_data/sg_endo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)



rst_endo <- subset_samples(physeq, c(Compartment =="Endosphere" & Rootstock=="RST-04-106")) %>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
rst_endo_sparcc = otu_table(rst_endo)
write.table(data.frame(OTU_id = rownames(rst_endo_sparcc), rst_endo_sparcc), file = "ms_run/data/sparcc_data/rst_endo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)



maxi_endo <- subset_samples(physeq, c(Compartment =="Endosphere" & Rootstock=="Maxifort")) %>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
maxi_endo_sparcc = otu_table(maxi_endo)
write.table(data.frame(OTU_id = rownames(maxi_endo_sparcc), maxi_endo_sparcc), file = "ms_run/data/sparcc_data/maxi_endo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)



# Rhizosphere

ng_rhizo <- subset_samples(physeq, c(Compartment =="Rhizosphere" & Rootstock=="Nongraft")) %>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
ng_rhizo_sparcc = otu_table(ng_rhizo)
write.table(data.frame(OTU_id = rownames(ng_rhizo_sparcc), ng_rhizo_sparcc), file = "ms_run/data/sparcc_data/ng_rhizo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)



sg_rhizo <- subset_samples(physeq, c(Compartment =="Rhizosphere" & Rootstock=="Selfgraft")) %>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
sg_rhizo_sparcc = otu_table(sg_rhizo)
write.table(data.frame(OTU_id = rownames(sg_rhizo_sparcc), sg_rhizo_sparcc), file = "ms_run/data/sparcc_data/sg_rhizo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)



rst_rhizo <- subset_samples(physeq, c(Compartment =="Rhizosphere" & Rootstock=="RST-04-106"))%>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
rst_rhizo_sparcc = otu_table(rst_rhizo)
write.table(data.frame(OTU_id = rownames(rst_rhizo_sparcc), rst_rhizo_sparcc), file = "ms_run/data/sparcc_data/rst_rhizo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)



maxi_rhizo <- subset_samples(physeq, c(Compartment =="Rhizosphere" & Rootstock=="Maxifort")) %>%
  prune_taxa(taxa_sums(.) > 2, .)
# Inputfile for running SparCC
maxi_rhizo_sparcc = otu_table(maxi_rhizo)
write.table(data.frame(OTU_id = rownames(maxi_rhizo_sparcc), maxi_rhizo_sparcc), file = "ms_run/data/sparcc_data/maxi_rhizo_sparcc.txt", sep = "\t", row.names = FALSE, quote = FALSE)


########## Running PhONA with SparCC results ##################

ng_endo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Endosphere_Nongraft_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
ng_endo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Endosphere_Nongraft_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)

pdf(file = "ms_run/data/figures/phona_ng_endo.pdf", width = 12, height = 8)
phona_ng_endo <- PhONA(physeqobj = ng_endo, cordata = ng_endo_sparcc.cor,
                       pdata = ng_endo_sparcc.pval, model = "lasso",
                       iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, defineTreatment = "Nongraft",nodesize = 5,
                       PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "Nongraft")
dev.off()






sg_endo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Endosphere_Selfgraft_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
sg_endo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Endosphere_Selfgraft_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)

pdf(file = "ms_run/data/figures/phona_sg_end.pdf", width = 12, height = 8)
phona_sg_endo <-PhONA(physeqobj = sg_endo, cordata = sg_endo_sparcc.cor,
                      pdata = sg_endo_sparcc.pval, model = "lasso",
                      iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, defineTreatment = "Selfgraft",nodesize = 5,
                      PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "Selfgraft")

dev.off()


rst_endo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Endosphere_RST-04-106_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
rst_endo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Endosphere_RST-04-106_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)

pdf(file = "ms_run/data/figures/phona_rst_endo.pdf", width = 12, height = 8)
phona_rst_endo <-PhONA(physeqobj = rst_endo, cordata = rst_endo_sparcc.cor,
                       pdata = rst_endo_sparcc.pval, model = "lasso",
                       iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, defineTreatment = "RST-04-106",nodesize = 5,
                       PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "RST-04-106")
dev.off()


maxi_endo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Endosphere_Maxifort_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
maxi_endo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Endosphere_Maxifort_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)

pdf(file = "ms_run/data/figures/phona_maxi_endo.pdf", width = 12, height = 8)
phona_maxi_endo <-PhONA(physeqobj = maxi_endo, cordata = maxi_endo_sparcc.cor,
                        pdata = maxi_endo_sparcc.pval, model = "lasso",
                        iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, defineTreatment = "Maxifort",nodesize = 5,
                        PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "Maxifort")
dev.off()






ng_rhizo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Rhizosphere_Nongraft_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
ng_rhizo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Rhizosphere_Nongraft_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)
pdf(file = "ms_run/data/figures/phona_ng_rhizo.pdf", width = 12, height = 8)
phona_ng_rhizo <- PhONA(physeqobj = ng_rhizo, cordata = ng_rhizo_sparcc.cor,
                        pdata = ng_rhizo_sparcc.pval, model = "lasso",
                        iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, defineTreatment = "Nongraft",nodesize = 5,
                        PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "Nongraft")


dev.off()

############
sg_rhizo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Rhizosphere_Selfgraft_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
sg_rhizo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Rhizosphere_Selfgraft_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)

pdf(file = "ms_run/data/figures/phona_sg_rhizo.pdf", width = 12, height = 8)
phona_sg_rhizo <- PhONA(physeqobj = sg_rhizo, cordata = sg_rhizo_sparcc.cor,
                        pdata = sg_rhizo_sparcc.pval, model = "lasso",
                        iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, defineTreatment = "Selfgraft",nodesize = 5,
                        PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "Selfgraft")
dev.off()

rst_rhizo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Rhizosphere_RST-04-106_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
rst_rhizo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Rhizosphere_RST-04-106_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)
pdf(file = "ms_run/data/figures/phona_rst_rhizo.pdf", width = 12, height = 8)
phona_rst_rhizo <- PhONA(physeqobj = rst_rhizo, cordata = rst_rhizo_sparcc.cor,
                         pdata = rst_rhizo_sparcc.pval, model = "lasso",
                         iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, defineTreatment = "RST-04-106",nodesize = 5,
                         PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "RST-04-106")

dev.off()

maxi_rhizo_sparcc.cor <- read.delim("ms_run/data/sparcc_output/Rhizosphere_Maxifort_cor_sparcc.out", sep = "\t", header = T, row.names = 1)
maxi_rhizo_sparcc.pval <- read.delim("ms_run/data/sparcc_output/Rhizosphere_Maxifort_pvals.two_sided.txt", sep = "\t", header = T, row.names = 1)
pdf(file = "ms_run/data/figures/phona_maxi_rhizo.pdf", width = 12, height = 8)
phona_maxi_rhizo <- PhONA(physeqobj = maxi_rhizo, cordata = maxi_rhizo_sparcc.cor,
                          pdata = maxi_rhizo_sparcc.pval, model = "lasso",
                          defineTreatment = "Maxifort",nodesize = 5,
                          iters = ITERS,OTU_OTU_pvalue=OTU_OTU_PVALUE, OTU_OTU_rvalue=OTU_OTU_RVALUE, OTU_Phenotype_pvalue=OTU_PHENOTYPE_PVALUE, PhenoNodesize = 12, definePhenotype = "Marketable", PhenoNodelabel = "Maxifort")

dev.off()


# Endosphere

role_df_Endosphere <- rbind(phona_ng_endo$roles,
                            phona_sg_endo$roles,
                            phona_rst_endo$roles,
                            phona_maxi_endo$roles)




fig3b_endo <-ggplot(role_df_Endosphere, aes(participation,connectivity,label=name)) + geom_point(aes(colour = factor(Treatment)))+
  scale_size(guide = 'none')+
  scale_color_manual(breaks= c("Nongraft", "Selfgraft", "RST-04-106", "Maxifort"),values = c("mediumaquamarine","sienna4","#ffae19","blueviolet"))+
  geom_text(aes(label=ifelse(role!="Peripheral",as.character(Order),'')),hjust=0.1,vjust=0.1,size=3,color="black",family="Times New Roman")+
  geom_hline(yintercept=2.5, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0.62, linetype="dashed",color = "blue")+
  annotate("text", x = 0.2, y = 3, label = "Module hubs",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.75, y = 3, label = " Network hubs",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.2, y = -2, label = "Peripherals",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.75, y = -2, label = "Connectors",color="gray20",size=8,family="Times New Roman")+
  labs(x = "Among-module connectivity", y = "Within-module degree", col="Rootstocks") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1.0))+
  theme(axis.text=element_text(size=16,face="bold",family="Times New Roman"),axis.title=element_text(size=16,face="bold",family="Times New Roman"))+
  labs(fill = "Dose (mg)") +
  theme(legend.title=element_text(size=16,face="bold",family="Times New Roman"),
        legend.text=element_text(size=16,family="Times New Roman")) +
  guides(color = guide_legend(override.aes = list(size=5)))

# Save the plot
ggsave("ms_run/data/figures/roles_endo.pdf", width = 8, height = 6, units = "in", device = cairo_pdf)
fig3b_endo
dev.off()


# Endosphere

role_df_Rhizosphere <- rbind(phona_ng_rhizo$roles,
                             phona_sg_rhizo$roles,
                             phona_rst_rhizo$roles,
                             phona_maxi_rhizo$roles)



fig3b_rhizo <-ggplot(role_df_Rhizosphere, aes(participation,connectivity,label=name)) + geom_point(aes(colour = factor(Treatment)))+
  scale_size(guide = 'none')+
  scale_color_manual(breaks= c("Nongraft", "Selfgraft", "RST-04-106", "Maxifort"),values = c("mediumaquamarine","sienna4","#ffae19","blueviolet"))+
  geom_text(aes(label=ifelse(role!="Peripheral",as.character(Order),'')),hjust=0.1,vjust=0.1,size=3,color="black",family="Times New Roman")+
  geom_hline(yintercept=2.5, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0.62, linetype="dashed",color = "blue")+
  annotate("text", x = 0.2, y = 3, label = "Module hubs",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.75, y = 3, label = " Network hubs",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.2, y = -2, label = "Peripherals",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.75, y = -2, label = "Connectors",color="gray20",size=8,family="Times New Roman")+
  labs(x = "Among-module connectivity", y = "Within-module degree", col="Rootstocks") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1.0))+
  theme(axis.text=element_text(size=16,face="bold",family="Times New Roman"),axis.title=element_text(size=16,face="bold",family="Times New Roman"))+
  labs(fill = "Dose (mg)") +
  theme(legend.title=element_text(size=16,face="bold",family="Times New Roman"),
        legend.text=element_text(size=16,family="Times New Roman")) +
  guides(color = guide_legend(override.aes = list(size=5)))

ggsave("ms_run/data/figures/roles_rhizo.pdf", width = 8, height = 6, units = "in", device = cairo_pdf)
fig3b_rhizo
dev.off()


summary_graph_endo_df <- rbind(phona_ng_endo$graph_summary,
                               phona_sg_endo$graph_summary,
                               phona_rst_endo$graph_summary,
                               phona_maxi_endo$graph_summary)
write.csv(summary_graph_endo_df,"ms_run/data/figures/summary_graph_endo.csv")

summary_graph_rhizo_df <- rbind(phona_ng_rhizo$graph_summary,
                                phona_sg_rhizo$graph_summary,
                                phona_rst_rhizo$graph_summary,
                                phona_maxi_rhizo$graph_summary)

write.csv(summary_graph_rhizo_df,"ms_run/data/figures/summary_graph_rhizo.csv")

# save the image
save.image(file = "phona_ms_run.RData")




















