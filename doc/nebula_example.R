## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("lhe17/nebula")

## ----echo=TRUE----------------------------------------------------------------
library(nebula)
data(sample_data)

## ----echo=TRUE----------------------------------------------------------------
dim(sample_data$count)

## ----echo=TRUE----------------------------------------------------------------
sample_data$count[1:5,1:5]

## ----echo=TRUE----------------------------------------------------------------
head(sample_data$sid)
table(sample_data$sid)

## ----echo=TRUE----------------------------------------------------------------
head(sample_data$pred)
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
head(df)

## ----echo=TRUE----------------------------------------------------------------
re = nebula(sample_data$count,sample_data$sid,pred=df,ncore=1)
re

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  data_g = group_cell(count=sample_data$count,id=sample_data$sid,pred=df)
#  re = nebula(data_g$count,data_g$id,pred=data_g$pred)

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  library(Matrix)
#  # An example of using the library size of each cell as the scaling factor
#  re = nebula(sample_data$count,sample_data$sid,pred=df,offset=Matrix::colSums(sample_data$count))

## ----echo=TRUE,eval=FALSE-----------------------------------------------------
#  library(nebula)
#  data("sample_seurat")
#  seuratdata <- scToNeb(obj = sample_seurat, assay = "RNA", id = "replicate", pred = c("celltype","tech"), offset="nCount_RNA")
#  ## Make sure that the variables do not contain NA; Otherwise, df would have fewer rows.
#  df = model.matrix(~celltype+tech, data=seuratdata$pred)
#  ## include only the first two cell types in the model to avoid separation due to too many binary variables
#  data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","celltypeactivated_stellate","techcelseq2","techfluidigmc1","techindrop", "techsmartseq2")],offset=seuratdata$offset)
#  re = nebula(data_g$count,data_g$id,pred=data_g$pred,offset=data_g$offset)

## ----eval=TRUE,echo=TRUE------------------------------------------------------
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',ncore=1)
re_hl = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='HL',ncore=1)
## compare the estimated overdispersions
cbind(re_hl$overdispersion,re_ln$overdispersion)

## ----eval=TRUE,echo=TRUE------------------------------------------------------
## compare the p-values for testing the predictors using NEBULA-LN and NEBULA-HL
cbind(re_hl$summary[,10:12],re_ln$summary[,10:12])

## ----eval=TRUE,echo=TRUE------------------------------------------------------
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,model='PMM',ncore=1)

## ----echo=FALSE,results='asis'------------------------------------------------
knitr::kable(re$summary)

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,model='NBLMM',reml=1,ncore=1)

## ----eval=TRUE,echo=TRUE------------------------------------------------------
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',covariance=TRUE,ncore=1)
cov= matrix(NA,4,4)
cov[lower.tri(cov,diag=T)] = as.numeric(re_ln$covariance[1,])
cov[upper.tri(cov)] = t(cov)[upper.tri(cov)]
cov

## ----eval=TRUE,echo=TRUE------------------------------------------------------
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
## the gene to test
gene_i = 1
## output covariance
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',covariance=TRUE,ncore=1)
## recover the covariance matrix
cov= matrix(NA,4,4)
cov[lower.tri(cov,diag=T)] = as.numeric(re_ln$covariance[gene_i,])
cov[upper.tri(cov)] = t(cov)[upper.tri(cov)]
## build the contrast vector
contrast = c(0,1,-1,0)
## testing the hypothesis
eff = sum(contrast*re_ln$summary[gene_i,1:4])
p = pchisq(eff^2/(t(contrast)%*%cov%*%contrast),1,lower.tail=FALSE)
p

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset)
#  pres = nbresidual(re,count=sample_data$count,id=sample_data$sid,pred=df,offset=sample_data$offset)

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,output_re=TRUE)

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  pres = nbresidual(re,count=sample_data$count,id=sample_data$sid,pred=df,offset=sample_data$offset,conditional=TRUE)

