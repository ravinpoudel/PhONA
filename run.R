physeqobj = phyobj
cordata = sparcc.cor
pdata = sparcc.pval
model = c("lm","lasso")
OTU_OTU_pvalue = 0.001
OTU_OTU_rvalue = 0.6
OTU_Phenotype_pvalue = 0.6
definePhenotype="Marketable"
defineTreatment="Maxifort"
coloredby="Phylum"
PhenoNodecolor="yellow"
PhenoNodesize=20
PhenoNodelabel="Yield"
nodesize=10
Pheno2OTUedgecolor = "black"
netlayout=layout.fruchterman.reingold


PhONA(physeqobj, model="lm")
PhONA(physeqobj, model="lasso")
