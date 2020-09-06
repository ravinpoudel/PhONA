physeqobj = phyobj
cordata = sparcc.cor
pdata = sparcc.pval
model = c("lm","lasso")
iters = 2
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
PhONA(physeqobj, model="lasso", iters = 2)

# https://rmarkdown.rstudio.com/docs/reference/md_document.html
# Convert to a markdown document
library(rmarkdown)
md_document(
  variant = "markdown_strict",
  preserve_yaml = FALSE,
  toc = FALSE,
  toc_depth = 3,
  #fig_width = 7,
  #fig_height = 5,
  fig_retina = NULL,
  dev = "png",
  df_print = "default",
  includes = NULL,
  md_extensions = NULL,
  pandoc_args = NULL,
  ext = ".md"
)


render("vignettes/PhONA.Rmd", md_document())



if (model == "lasso"){
  df_list <- list()
  for (i in 1:iters){
    message ("Running Iteration Number::", i)
    tryCatch({df_list[[i]] = model.lasso(x, odata)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  # flat list to df
  df <- do.call(rbind,df_list)
  # filter combined df based on unique otus
  df_uniq = df %>%
    distinct(otus, .keep_all = TRUE)
  df_uniq
  }

install.packages("tictoc")
require(tictoc)
tic()
rnorm(1000,0,1)
toc()

#############
#Introduction to pkgdown
# Run once to configure package to use pkgdown
#usethis::use_pkgdown()
# Run to build the website
document()
build()
pkgdown::build_site()




