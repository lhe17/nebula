## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("lhe17/nebula")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("nebula", repos="http://R-Forge.R-project.org")

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
re = nebula(sample_data$count,sample_data$sid,pred=df)
re

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  data_g = group_cell(count=sample_data$count,id=sample_data$sid,pred=df)
#  re = nebula(data_g$count,data_g$id,pred=data_g$pred)

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset)

## ----eval=TRUE,echo=TRUE------------------------------------------------------
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN')
re_hl = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='HL')
## compare the estimated overdispersions
cbind(re_hl$overdispersion,re_ln$overdispersion)

## ----eval=TRUE,echo=TRUE------------------------------------------------------
## compare the p-values for testing the predictors using NEBULA-LN and NEBULA-HL
cbind(re_hl$summary[,10:12],re_ln$summary[,10:12])

## ----eval=TRUE,echo=TRUE------------------------------------------------------
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,model='PMM')

## ----echo=FALSE,results='asis'------------------------------------------------
knitr::kable(re$summary)

