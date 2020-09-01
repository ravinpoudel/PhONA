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


# Convert to a markdown document
md_document(
  variant = "markdown_strict",
  preserve_yaml = FALSE,
  toc = FALSE,
  toc_depth = 3,
  fig_width = 7,
  fig_height = 5,
  fig_retina = NULL,
  dev = "png",
  df_print = "default",
  includes = NULL,
  md_extensions = NULL,
  pandoc_args = NULL,
  ext = ".md"
)


library(rmarkdown)
render("vignettes/PhONA.Rmd", md_document())

