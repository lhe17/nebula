-   [NEBULA v1.5.3](#nebula-v1.5.3)
    -   [Overview](#overview)
    -   [Installation](#installation)
        -   [Most recent version](#most-recent-version)
    -   [Functions](#functions)
    -   [Basic usage](#basic-usage)
        -   [Example](#example)
    -   [Specifying scaling factors](#specifying-scaling-factors)
        -   [Example](#example-1)
    -   [Using Seurat/SingleCellExperiment
        Objects](#using-seuratsinglecellexperiment-objects)
        -   [Example](#example-2)
    -   [Difference between NEBULA-LN and
        NEBULA-HL](#difference-between-nebula-ln-and-nebula-hl)
    -   [Filtering low-expression
        genes](#filtering-low-expression-genes)
    -   [Checking convergence for the summary statistics and quality
        control](#checking-convergence-for-the-summary-statistics-and-quality-control)
    -   [Using other mixed models](#using-other-mixed-models)
        -   [Example](#example-3)
    -   [Special attention paid to testing subject-level
        variables](#special-attention-paid-to-testing-subject-level-variables)
        -   [Example](#example-4)
    -   [Testing contrasts](#testing-contrasts)
    -   [Extracting marginal and conditional Pearson
        residuals](#extracting-marginal-and-conditional-pearson-residuals)
    -   [Parallel computing](#parallel-computing)
    -   [References](#references)

# NEBULA v1.5.3

## Overview

The *nebula* package is an R package that provides fast algorithms for
fitting negative binomial and Poisson mixed models for analyzing
large-scale, multi-subject single-cell data. The package *nebula*
accounts for the hierarchical structure of the data by decomposing the
total overdispersion into between-subject and within-subject components
using a negative binomial mixed model (NBMM). Users can utilize the
package for various tasks, such as identifying marker genes, testing
treatment effects, detecting genes with differential expression,
performing cell-level co-expression analysis, and obtaining Pearson
residuals for downstream analyses.

More details can be found in (He et al. 2021)
(<https://www.nature.com/articles/s42003-021-02146-6>).

## Installation

### Most recent version

To install the latest version from github:

``` r
install.packages("devtools")
library(devtools)
install_github("lhe17/nebula")
```

During installation, the *nebula* package may first install the *Rfast*
package, which requires the presence of GSL in the environment. The
installation also requires Rcpp-1.0.7 and has been tested on R-4.1.0.
Starting from version 1.2.0, *nebula* is no longer compatible with R-3.6
or earlier versions of R. Users who have R-3.6 may install version 1.1.8
via R-forge (<https://r-forge.r-project.org/R/?group_id=2407>). However,
it is not recommended to use an older version of *nebula*.

Please contact <hyx520101@gmail.com> for more information.

## Functions

The current version provides the following functions.

-   `nebula`: performs an association analysis using NBMMs given a count
    matrix and subject IDs.
-   `group_cell`: reorders cells to group them by the subject IDs.
-   `nbresidual`: extracts Pearson residuals from the fitted model.
-   `scToNeb`: retrieves data from `Seurat` or `SingleCellExperiment`
    for calling `nebula`.

## Basic usage

We use an example data set to illustrate how to use nebula to perform an
association analysis of multi-subject single-cell data. The example data
set attached to the R package can be loaded as follows.

``` r
library(nebula)
data(sample_data)
```

The example data set includes a count matrix of 6030 cells and 10 genes
from 30 subjects.

``` r
dim(sample_data$count)
#> [1]   10 6176
```

The count matrix can be a matrix object or a sparse dgCMatrix object
(the same format as in `Seurat`). The elements should be integers.

``` r
sample_data$count[1:5,1:5]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>            
#> A . . . . .
#> B . . . . .
#> C . 1 2 . .
#> D . . . . .
#> E . . . . .
```

The subject IDs of each cell are stored in `sample_data$sid`. The
subject IDs can be a character or numeric vector, the length of which
should equal the number of cells.

``` r
head(sample_data$sid)
#> [1] "1" "1" "1" "1" "1" "1"
table(sample_data$sid)
#> 
#>   1  10  11  12  13  14  15  16  17  18  19   2  20  21  22  23  24  25  26  27 
#> 187 230 185 197 163 216 211 195 200 239 196 223 198 202 213 210 199 214 237 200 
#>  28  29   3  30   4   5   6   7   8   9 
#> 205 183 222 191 205 225 211 197 215 207
```

The next step is to build a design matrix for the predictors. The
example data set includes a data frame consisting of three predictors
stored in `sample_data$pred`. To build the design matrix, we can use the
function `model.matrix`. The intercept term must be included in the
design matrix. Each column in the design matrix should have a unique
variable name.

``` r
head(sample_data$pred)
#>           X1        X2      cc
#> 1  0.6155094 0.9759191 control
#> 2  1.4608092 0.9759191    case
#> 3  1.6675054 0.9759191 control
#> 4 -0.1717715 0.9759191    case
#> 5  0.2277492 0.9759191 control
#> 6 -0.2635516 0.9759191 control
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
head(df)
#>   (Intercept)         X1        X2 cccontrol
#> 1           1  0.6155094 0.9759191         1
#> 2           1  1.4608092 0.9759191         0
#> 3           1  1.6675054 0.9759191         1
#> 4           1 -0.1717715 0.9759191         0
#> 5           1  0.2277492 0.9759191         1
#> 6           1 -0.2635516 0.9759191         1
```

The association analysis between the gene expression and the predictors
can then be conducted using the `nebula` function as follows. The count
matrix is an *M* by *N* matrix, where *M* is the number of genes, and
*N* is the number of cells. The function by default fits the negative
binomial gamma mixed model (NBGMM) for each of the genes, and returns a
list of summary statistics including the fold change, p-values, and both
subject-level and cell-level overdispersions (*σ*<sup>2</sup> and
*ϕ*<sup>−1</sup>). The p-values returned by `nebula` are raw p-values
(not adjusted for multiple testing). Users can take advantage of a
multicore CPU by specifying the number of cores to use via the `ncore`
argument.

``` r
re = nebula(sample_data$count,sample_data$sid,pred=df,ncore=1)
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
#> Loading required package: foreach
#> Loading required package: future
#> Loading required package: rngtools
re
#> $summary
#>    logFC_(Intercept)     logFC_X1     logFC_X2 logFC_cccontrol se_(Intercept)
#> 1          -1.902455 -0.016755225 -0.097867225     0.047278197     0.06335820
#> 2          -2.046638 -0.002679074 -0.053812464    -0.022293899     0.06181112
#> 3          -2.033211  0.017954707  0.002398445    -0.048296661     0.08695028
#> 4          -2.008542 -0.005698984 -0.027780387     0.077357703     0.05509711
#> 5          -1.979437  0.011557090 -0.025198987     0.032890493     0.06155853
#> 6          -1.949991  0.013483039 -0.012548791    -0.031590577     0.07440949
#> 7          -1.969248 -0.003531361  0.075230699    -0.009075031     0.06185028
#> 8          -1.964371  0.013639930 -0.061302756    -0.059284665     0.07786361
#> 9          -2.072699 -0.017372176 -0.043828288     0.026624998     0.05737632
#> 10         -2.045646  0.030742876  0.022260805    -0.025516032     0.06842796
#>         se_X1      se_X2 se_cccontrol p_(Intercept)      p_X1      p_X2
#> 1  0.03534659 0.06449424   0.06879634 4.362617e-198 0.6354810 0.1291514
#> 2  0.03787429 0.06255849   0.07385888 2.052788e-240 0.9436079 0.3896819
#> 3  0.03696089 0.09238230   0.07258521 6.275230e-121 0.6271261 0.9792875
#> 4  0.03704556 0.05624824   0.07252600 5.822948e-291 0.8777381 0.6213846
#> 5  0.03750948 0.06101307   0.07331551 7.432319e-227 0.7579977 0.6795995
#> 6  0.03623477 0.07321208   0.07087566 2.257914e-151 0.7098168 0.8639067
#> 7  0.03631619 0.06068697   0.07133730 1.872102e-222 0.9225364 0.2151043
#> 8  0.03551903 0.07955877   0.06969748 1.957495e-140 0.7009654 0.4409831
#> 9  0.03816039 0.05767972   0.07453316 9.307495e-286 0.6489358 0.4473406
#> 10 0.03798694 0.06917485   0.07374591 2.292903e-196 0.4183419 0.7476005
#>    p_cccontrol gene_id gene
#> 1    0.4919443       1    A
#> 2    0.7627706       2    B
#> 3    0.5058082       3    C
#> 4    0.2861434       4    D
#> 5    0.6537089       5    E
#> 6    0.6558008       6    F
#> 7    0.8987718       7    G
#> 8    0.3949916       8    H
#> 9    0.7209245       9    I
#> 10   0.7293432      10    J
#> 
#> $overdispersion
#>       Subject      Cell
#> 1  0.08125256 0.8840821
#> 2  0.07102681 0.9255032
#> 3  0.17159404 0.9266395
#> 4  0.05026165 0.8124118
#> 5  0.07075366 1.2674146
#> 6  0.12086392 1.1096065
#> 7  0.07360445 0.9112956
#> 8  0.13571262 0.7549629
#> 9  0.05541398 0.8139652
#> 10 0.09496649 0.9410035
#> 
#> $convergence
#>  [1] 1 1 1 1 1 1 1 1 1 1
#> 
#> $algorithm
#>  [1] "NBGMM (LN)" "NBGMM (LN)" "NBGMM (LN)" "NBGMM (LN)" "NBGMM (LN)"
#>  [6] "NBGMM (LN)" "NBGMM (LN)" "NBGMM (LN)" "NBGMM (LN)" "NBGMM (LN)"
#> 
#> $covariance
#> NULL
#> 
#> $random_effect
#> NULL
```

The cells in the count matrix need to be grouped by the subjects (that
is, the cells of the same subject should be placed consecutively) before
using as the input to the function `nebula`. If the cells are not
grouped, the function `group_cell` can be used to first reorder the
cells, as shown below. If a scaling factor is specified by the user, it
should also be included in `group_cell`. If the cells are already
grouped, `group_cell` will return *NULL*.

### Example

``` r
data_g = group_cell(count=sample_data$count,id=sample_data$sid,pred=df)
re = nebula(data_g$count,data_g$id,pred=data_g$pred)
```

If `pred` is not specified, `nebula` will fit the model with an
intercept term by default. This can be used when only the
overdispersions are of interest.

## Specifying scaling factors

The scaling factor for each cell is specified in `nebula` using the
argument `offset`. The argument `offset` has to be a vector of length
*N* containing positive values. Note that log(`offset`) will be the
offset term in the NBMM. Common scaling factors can be the library size
of a cell or a normalizing factor adjusted using e.g., TMM. If not
specified, `nebula` will set `offset` as one by default, which means
that each cell is treated equally. If the input count matrix is already
normalized by another tool, e.g., scTransform, then you should not
specify `offset`. However, since `nebula` directly models the raw
counts, it is not recommended to use a normalized count matrix for
`nebula`.

### Example

``` r
library(Matrix)
# An example of using the library size of each cell as the scaling factor
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=Matrix::colSums(sample_data$count))
```

## Using Seurat/SingleCellExperiment Objects

If a single cell data processing package such as `Seurat` or
`SingleCellExperiment` was used, nebula can be easily implemented using
the assistance of the helper function `scToNeb`. Assuming that the
metadata relevant to subject IDs and predictors are available in the
object, `scToNeb` can retrieve and organize these objects and output a
list that is similar to the example data provided in this vignette. For
a `SingleCellExperiment` object, `assay` is not required. For a `Seurat`
object, `assay` can also be specified to fit data from other assays. The
`nebula` package contains a sample Seurat object obtained from (Lab
2019) (<https://github.com/satijalab/seurat-data>) comprised of
pancreatic cells across eight samples.

### Example

``` r
library(nebula)
data("sample_seurat")
seuratdata <- scToNeb(obj = sample_seurat, assay = "RNA", id = "replicate", pred = c("celltype","tech"), offset="nCount_RNA")
## Make sure that the variables do not contain NA; Otherwise, df would have fewer rows.
df = model.matrix(~celltype+tech, data=seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","celltypeactivated_stellate","techcelseq2","techfluidigmc1","techindrop", "techsmartseq2")],offset=seuratdata$offset)
re = nebula(data_g$count,data_g$id,pred=data_g$pred,offset=data_g$offset)
```

The output will be a list with the first element containing `counts`,
the second containing a `data.frame` with all listed predictors, the
third containing a character vector with all subject IDs, and the fourth
containing the normalizing factor. Users can also use other scaling
factors that may be stored within the object’s metadata as a string in
the `offset` argument. If subject ids are un-ordered, `group_cell` can
be used.

## Difference between NEBULA-LN and NEBULA-HL

In *nebula*, a user can choose one of the two algorithms to fit an
NBGMM. NEBULA-LN uses an approximated likelihood based on the law of
large numbers, and NEBULA-HL uses an h-likelihood. A user can select
these methods through `method='LN'` or `method='HL'`. NEBULA-LN is
faster and performs particularly well when the number of cells per
subject (CPS) is large. In addition, NEBULA-LN is much more accurate in
estimating a very large subject-level overdispersion. In contrast,
NEBULA-HL is slower but more accurate in estimating the cell-level
overdispersion.

In the following analysis of the example data set comprising ~200 cells
per subject, the difference of the estimated cell-level overdispersions
between NEBULA-LN and NEBULA-HL is ~5% for most genes.

``` r
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',ncore=1)
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
re_hl = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='HL',ncore=1)
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
## compare the estimated overdispersions
cbind(re_hl$overdispersion,re_ln$overdispersion)
#>       Subject      Cell    Subject      Cell
#> 1  0.08432318 0.9284699 0.08125256 0.8840821
#> 2  0.07455464 0.9726513 0.07102681 0.9255032
#> 3  0.17403263 0.9817569 0.17159404 0.9266395
#> 4  0.05352153 0.8516679 0.05026165 0.8124118
#> 5  0.07480033 1.3254379 0.07075366 1.2674146
#> 6  0.12372424 1.1653129 0.12086392 1.1096065
#> 7  0.07724825 0.9578169 0.07360445 0.9112956
#> 8  0.13797645 0.7991948 0.13571262 0.7549629
#> 9  0.05879495 0.8568850 0.05541398 0.8139652
#> 10 0.09782333 0.9940222 0.09496649 0.9410035
```

Such difference has little impact on testing fixed-effects predictors
under this sample size.

``` r
## compare the p-values for testing the predictors using NEBULA-LN and NEBULA-HL
cbind(re_hl$summary[,10:12],re_ln$summary[,10:12])
#>         p_X1      p_X2 p_cccontrol      p_X1      p_X2 p_cccontrol
#> 1  0.6373036 0.1346298   0.4950795 0.6354810 0.1291514   0.4919443
#> 2  0.9444825 0.3977109   0.7626827 0.9436079 0.3896819   0.7627706
#> 3  0.6282384 0.9787881   0.5087304 0.6271261 0.9792875   0.5058082
#> 4  0.8786074 0.6278827   0.2868256 0.8777381 0.6213846   0.2861434
#> 5  0.7596198 0.6872259   0.6544751 0.7579977 0.6795995   0.6537089
#> 6  0.7134192 0.8656686   0.6576835 0.7098168 0.8639067   0.6558008
#> 7  0.9216994 0.2230964   0.8977251 0.9225364 0.2151043   0.8987718
#> 8  0.7017082 0.4443604   0.3955342 0.7009654 0.4409831   0.3949916
#> 9  0.6505414 0.4561470   0.7238323 0.6489358 0.4473406   0.7209245
#> 10 0.4199828 0.7510837   0.7308108 0.4183419 0.7476005   0.7293432
```

The bias of NEBULA-LN in estimating the cell-level overdispersion gets
larger when the CPS value becomes lower or the gene expression is more
sparse. If the CPS value is \<30, `nebula` will set `method='HL'`
regardless of the user’s input.

When NEBULA-LN is used, the user can opt for better accuracy of
estimating a smaller subject-level overdispersion through the argument
*κ*. NEBULA first fits the data using NEBULA-LN. If the estimated *κ*
for a gene is smaller than the user-defined value, NEBULA-HL will be
used to estimate the subject-level overdispersion for the gene. The
default value of *κ* is 800, which can provide a good estimate of the
subject-level overdispersion as low as ~0.005. Our simulation results
suggest that *κ* = 200 is often sufficient for achieving a well
controlled false positive rate of testing a cell-level predictor. We do
not recommend using a smaller *κ* than 200. Specifying a larger *κ* can
obtain a more accurate estimate of a smaller subject-level
overdispersion when the cell-level overdispersion is large, but will be
computationally slower. On the other hand, testing a subject-level
predictor (i.e., a variable whose values are shared across all cells
from a subject, such as age, sex, treatment, genotype, etc) is more
sensitive to the accuracy of the estimated subject-level overdispersion.
So we recommend using *κ* = 800 (as default) or even larger when testing
a subject-level predictor. Another option to testing a subject-level
predictor is to use a Poisson gamma mixed model, which is extremely fast
(\>50x faster than NEBULA-LN) and will be described below.

## Filtering low-expression genes

NEBULA-HL automatically uses a higher-order Laplace approximation for
low-expressed genes of which the average count per subject is less than
3. The higher-order Laplace approximation substantially increases the
accuracy for estimating the subject-level overdispersion for
low-expressed genes and controls the false positive rate. Nevertheless,
we recommend removing genes with very low expression from the analysis
because there is little statistical power for these genes. Filtering out
low-expressed genes can be specified by `cpc=0.005` (i.e., counts per
cell\<0.5%). The argument `cpc` is defined by the ratio between the
total count of the gene and the number of cells.

## Checking convergence for the summary statistics and quality control

*nebula* reports convergence information about the estimation algorithm
for each gene along with the summary statistics. This is useful and
important information for quality control to filter out genes of which
the estimation procedure potentially does not converge. Generally, a
convergence code ≤ -20 suggests that the algorithm does not converge
well. The results should be interpreted with caution in these cases. The
detailed information about the convergence codes is listed below. The
failure of convergence may occur when the sample size is very small,
there are too few positive counts, or the gene has huge overdispersions.
In these cases, the likelihood can be flat, might reach the maximum at
the infinity, or the optimization is sensitive to the initial values.
For those genes that have a bad convergence code, in many cases, trying
a different negative binomial mixed model (e.g., NBLMM, see below for
more details) may solve the problem.

-   Information about the convergence code:
    -   1: The convergence is reached due to a sufficiently small
        improvement of the function value.
    -   -10: The convergence is reached because the gradients are close
        to zero (i.e., the critical point) and no improvement of the
        function value can be found.
    -   (!) -20: The optimization algorithm stops before the convergence
        because the maximum number of iterations is reached.
    -   (!) -25: The Hessian matrix is either almost singular or not
        positive definite.
    -   (!) -30: The convergence fails because the likelihood function
        returns NaN.  
    -   (!) -40: The convergence fails because the critical point is not
        reached and no improvement of the function value can be found.
    -   (!) -50: Only used for the PMM, indicating a failure of
        convergence.
    -   (!) -60: At least one of the estimated overdispersions reaches
        its upper bound.

Depending on the concrete application, the estimated gene-specific
overdispersions can also be taken into consideration in quality control.
For example, when testing differential expression for a variable, genes
with a very large estimated cell-level overdispersion should be filtered
out because such genes have huge unexplained noises. A large cell-level
overdispersion is generally rare in UMI-based single cell data,
especially among abundantly expressed genes, but more common in e.g.,
SMART-seq2 as PCR duplicates introduce substantial noises. It might be
hard to give a precise cut-off for a large overdispersion because it
also depends on the sample size of the data. Based on the empirical
simulation study in (He et al. 2021), genes with an estimated cell-level
overdispersion \>100 should be removed for a data set with at least 50
cells per subject. On the other hand, if the purpose is to extract
residuals for downstream analysis such as clustering, genes with a large
cell-level overdispersion might be preferable because they have large
variations. If the variable of interest is subject-level, genes with a
very large subject-level overdispersion (\>1) should be removed or
interpreted cautiously as well.

## Using other mixed models

In addition to the NBGMM, the *nebula* package provides efficient
estimation implementation for a Poisson gamma mixed model and a negative
binomial lognormal mixed model (NBLMM). This can be specified through
`model="PMM"` and `model="NBLMM"`, respectively. The NBLMM is the same
model as that adopted in the `glmer.nb` function in the *lme4* R
package, but is computationally much more efficient by setting
`method='LN'`. The only difference between NBGMM and NBLMM is that NBGMM
uses a gamma distribution for the random effects while the NBLMM uses a
lognormal distribution. The PMM is the fastest among these models. Note
that the Poisson mixed model (PMM) should not be used to test a
cell-level predictor because it only estimates the subject-level
overdispersion. Here is an example of using the PMM to fit the example
data set.

### Example

``` r
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,model='PMM',ncore=1)
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
```

| logFC\_(Intercept) |   logFC_X1 |   logFC_X2 | logFC_cccontrol | se\_(Intercept) |     se_X1 |     se_X2 | se_cccontrol | p\_(Intercept) |      p_X1 |      p_X2 | p_cccontrol | gene_id | gene |
|-------:|----:|----:|------:|------:|----:|----:|-----:|-----:|----:|----:|-----:|---:|:--|
|          -1.903559 | -0.0155807 | -0.0976567 |       0.0511051 |       0.0661288 | 0.0329114 | 0.0655551 |    0.0642298 |              0 | 0.6359176 | 0.1363061 |   0.4262299 |       1 | A    |
|          -2.047864 | -0.0032670 | -0.0536887 |      -0.0189269 |       0.0644332 | 0.0355074 | 0.0635450 |    0.0694853 |              0 | 0.9266904 | 0.3981703 |   0.7853239 |       2 | B    |
|          -2.032603 |  0.0179868 |  0.0009264 |      -0.0505295 |       0.0908162 | 0.0345496 | 0.0932444 |    0.0676704 |              0 | 0.6026385 | 0.9920732 |   0.4552440 |       3 | C    |
|          -2.009743 | -0.0055009 | -0.0278638 |       0.0782083 |       0.0573206 | 0.0350744 | 0.0574457 |    0.0686939 |              0 | 0.8753743 | 0.6276435 |   0.2549097 |       4 | D    |
|          -1.980527 |  0.0106340 | -0.0248788 |       0.0312191 |       0.0644293 | 0.0343354 | 0.0621582 |    0.0671645 |              0 | 0.7567817 | 0.6889725 |   0.6420637 |       5 | E    |
|          -1.950454 |  0.0160303 | -0.0134764 |      -0.0345278 |       0.0778201 | 0.0333858 | 0.0738508 |    0.0650410 |              0 | 0.6311185 | 0.8552054 |   0.5955144 |       6 | F    |
|          -1.970271 | -0.0026762 |  0.0750061 |      -0.0063660 |       0.0645989 | 0.0341936 | 0.0615159 |    0.0668723 |              0 | 0.9376159 | 0.2227322 |   0.9241583 |       7 | G    |
|          -1.964322 |  0.0141545 | -0.0610982 |      -0.0578682 |       0.0809950 | 0.0336580 | 0.0800989 |    0.0656801 |              0 | 0.6740920 | 0.4455920 |   0.3782847 |       8 | H    |
|          -2.074035 | -0.0178188 | -0.0436111 |       0.0259748 |       0.0597958 | 0.0362203 | 0.0587687 |    0.0707913 |              0 | 0.6227507 | 0.4580383 |   0.7136780 |       9 | I    |
|          -2.046058 |  0.0307022 |  0.0227234 |      -0.0246096 |       0.0714146 | 0.0354844 | 0.0702255 |    0.0691813 |              0 | 0.3869124 | 0.7462578 |   0.7220452 |      10 | J    |

## Special attention paid to testing subject-level variables

When testing subject-level variables, it should be kept in mind that the
actual sample size is the number of subjects, not the number of cells in
the data set. At least a moderate number of subjects (\>30) are required
for testing a subject-level variable using `nebula` simply because a
small number of subjects are not enough to accurately estimate the
subject-level overdispersion. As shown in the original article (He et
al. 2021), even 30 subjects lead to mild inflated type I errors in most
simulated scenarios. If the number of subjects is very small, methods
designed for small sample size (e.g., DESeq2, edgeR) should be used for
testing subject-level variables.

In addition, when the ratio between the number of subjects and the
number of subject-level variables is small (\<10), it is recommended to
instead use a restricted maximum likelihood (REML) estimate, which is
provided in `nebula` through the argument `reml`. Please see (He et al.
2021) for more details about the formula of REML. This is because the
number of subject-level fixed-effects parameters should be much smaller
than the number of subjects in order to make the maximum likelihood
estimation (MLE) work properly. For example, if the data set has 50
subjects, it is a good practice to keep the number of subject-level
variables below 5 based on our simulation study. Increasing the number
of subject-level parameters will gradually inflate the type I error rate
due to an underestimated overdispersion. When these two numbers are at
the same magnitude, the MLE for the overdispersion will break down and
consequently, the NBMM can degenerate to a negative binomial model. In
contrast, REML takes into account the uncertainty of the estimated fixed
effects and controls the false positive rate even if many subject-level
covariates are included in the model. As shown in the following example,
one could simply specify `reml=1` to use REML, which is supported only
for `model='NBLMM'` in the current version.

### Example

``` r
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,model='NBLMM',reml=1,ncore=1)
```

## Testing contrasts

In some situations, a user may want to test a combination (contrast) of
the log(FC) or perform a global test for multiple variables or levels.
For example, a user may want to test whether the log(FC) of two
variables are the same. Here, we show how `nebula` can be used for this
kind of analysis.

The first step is to tell `nebula` to output the covariance matrix of
the estimated log(FC). This can be done by specifying `covariance=TRUE`
in `nebula`. To save storage, the covariance returned by `nebula` only
contains the elements in the lower triangular part including the
diagonal. Here is an example to recover the covariance matrix from the
output of `nebula`.

``` r
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',covariance=TRUE,ncore=1)
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
cov= matrix(NA,4,4)
cov[lower.tri(cov,diag=T)] = as.numeric(re_ln$covariance[1,])
cov[upper.tri(cov)] = t(cov)[upper.tri(cov)]
cov
#>               [,1]          [,2]         [,3]          [,4]
#> [1,]  4.014261e-03  2.499051e-05 1.384999e-04 -5.197643e-05
#> [2,]  2.499051e-05  1.249382e-03 9.212341e-06 -1.167080e-05
#> [3,]  1.384999e-04  9.212341e-06 4.159507e-03  5.142249e-05
#> [4,] -5.197643e-05 -1.167080e-05 5.142249e-05  4.732936e-03
```

Note that if there are *K* variables, the covariance table in the output
will have *(K+1)K/2* columns. So, for a large *K*, substantial increase
of computational intensity should be expected.

The second step is to build the contrast vector for your hypothesis. In
this example, we want to test whether the log(FC) of *X1* and *X2* are
equal for the first gene. This hypothesis leads to the contrast vector
`(0 1 -1 0)`. Thus, the test can be performed using the following code.

``` r
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
## the gene to test
gene_i = 1
## output covariance
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',covariance=TRUE,ncore=1)
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
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
#>           [,1]
#> [1,] 0.2692591
```

## Extracting marginal and conditional Pearson residuals

Pearson residuals are the distances between the raw count and its
expected value standardized by its standard deviation. Pearson residuals
obtained from fitting the NBMM can be used for normalization and
downstream analyses. The marginal Pearson residuals are obtained by
removing from the raw count the contribution from all fixed-effect
variables included in the model. The conditional Pearson residuals
further remove the subject-level random effects, which capture the
contribution of all other potential subject-level variables that are not
explicitly included in the model. Therefore, the conditional Pearson
residuals are very useful in a situation where one needs to remove the
subject-level batch effects from the normalized residuals for downstream
analyses.

Both Pearson residuals can be easily extracted by using the `nbresidual`
function after successfully running the `nebula` function. To extract
the marginal Pearson residuals, one provides in `nbresidual` the object
returned by `nebula` together with the same arguments including the
count matrix, `id`, `pred` and `offset` used in running the `nebula`
function. Here is an example.

``` r
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset)
pres = nbresidual(re,count=sample_data$count,id=sample_data$sid,pred=df,offset=sample_data$offset)
```

The parameters `count`, `id`, `pred` and `offset` should be the same in
these two functions. Then, the marginal Pearson residuals are available
in the matrix `pres$residuals`. The rows in `pres$residuals` correspond
to the genes in the output of `nebula`, and the columns are the cells in
`count`.

To extract the conditional Pearson residuals, we need to first let
`nebula` output subject-level random effects by setting `output_re=TRUE`
when running `nebula` as shown below.

``` r
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,output_re=TRUE)
```

The returned object will include an *M* by *L* matrix of the random
effects, where *L* is the number of subjects. In the current version,
this option does NOT support `model="PMM"`. Then, the conditional
Pearson residuals can be extracted by running `nbresidual` with
`conditional=TRUE`.

``` r
pres = nbresidual(re,count=sample_data$count,id=sample_data$sid,pred=df,offset=sample_data$offset,conditional=TRUE)
```

## Parallel computing

Starting with version 1.4.0, *nebula* supports parallel computing to
accelerate tasks. To specify the number of logical cores (threads),
users can use the `ncore` argument when running `nebula`. By default,
`nebula` uses two cores if more than two cores are available. However,
if the specified value for `ncore` exceeds the number of available
cores, *nebula* will use available cores detected and issue a warning.
It’s important to note that the maximum number of available processors
on a computing cluster may be restricted by the job scheduler’s
configuration.

While parallel computing significantly improves computational time,
users should monitor memory usage when exploiting multiple cores.
Insufficient memory allocation to a new thread may cause the R session
or job to be terminated by the operating system. In addition, if the
number of cells is small, using too many cores might not improve or even
reduce the efficiency because the overhead for creating new threads
exceeds the speed gain.

## References

He, Liang, Jose Davila-Velderrain, Tomokazu S. Sumida, David A. Hafler,
Manolis Kellis, and Alexander M. Kulminski. 2021. “NEBULA Is a Fast
Negative Binomial Mixed Model for Differential or Co-Expression Analysis
of Large-Scale Multi-Subject Single-Cell Data.” *Communications
Biology*, no. 629 (May). <https://doi.org/10.1038/s42003-021-02146-6>.

Lab, Satija. 2019. *Panc8.SeuratData: Eight Pancreas Datasets Across
Five Technologies*.
