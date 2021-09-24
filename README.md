-   [NEBULA v1.1.8](#nebula-v1.1.8)
    -   [Overview](#overview)
    -   [Installation](#installation)
        -   [Most recent version](#most-recent-version)
    -   [Functions](#functions)
    -   [Basic usage](#basic-usage)
        -   [Example](#example)
    -   [Specifying scaling factors](#specifying-scaling-factors)
        -   [Example](#example-1)
    -   [Selection between NEBULA-LN and
        NEBULA-HL](#selection-between-nebula-ln-and-nebula-hl)
    -   [Filtering low-expressed genes](#filtering-low-expressed-genes)
    -   [Checking convergence for the summary
        statistics](#checking-convergence-for-the-summary-statistics)
    -   [Using other mixed models](#using-other-mixed-models)
        -   [Example](#example-2)
    -   [Testing contrasts](#testing-contrasts)

NEBULA v1.1.8
=============

Overview
--------

The R package, *nebula*, provides fast algorithms for fitting negative
binomial and Poisson mixed models for analyzing large-scale
multi-subject single-cell data. The package *nebula* accounts for the
hierarchical structure of the data by decomposing the total
overdispersion into between-subject and within-subject components using
a negative binomial mixed model (NBMM). The package nebula can be used
for e.g., identifying marker genes, testing treatment effects, detecting
genes with differentail expression, and performing cell-level
co-expression analysis.

More details can be found in the manuscript “NEBULA: a fast negative
binomial mixed model for differential expression and co-expression
analyses of large-scale multi-subject single-cell data”
(<a href="https://www.nature.com/articles/s42003-021-02146-6" class="uri">https://www.nature.com/articles/s42003-021-02146-6</a>).

Installation
------------

### Most recent version

To install the lastest version from github:

``` r
install.packages("devtools")
library(devtools)
install_github("lhe17/nebula")
```

To install the lastest version from R-forge:

``` r
install.packages("nebula", repos="http://R-Forge.R-project.org")
```

Because the package *nebula* uses the R package *Rfast*, the
installation process may first install *Rfast*, which requires that GSL
is installed or available in the environment.

The installation has been tested on R-3.6 and R-3.5. Please contact
<a href="mailto:liang.he@duke.edu" class="email">liang.he@duke.edu</a>
for more information.

Functions
---------

The current version provides the following functions.

-   `nebula`: performs an association analysis using NBMMs given a count
    matrix and subject IDs.
-   `group_cell`: reorders cells to group them by the subject IDs.

Basic usage
-----------

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

The count matrix can be a matrix object or a sparse dgCMatrix object.
The elements should be integers.

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
design matrix.

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
can then be conducted using the function `nebula`. The count matrix is
an *M* by *N* matrix, where *M* is the number of genes, and *N* is the
number of cells.

``` r
re = nebula(sample_data$count,sample_data$sid,pred=df)
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
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
```

The function by default fitted the negative binomial gamma mixed model
(NBGMM) for each of the genes, and return a list of summary statistics
including the fold change, p-values, and both subject-level and
cell-level overdispersions (*σ*<sup>2</sup> and *ϕ*<sup> − 1</sup>). The
cells need to be grouped by the subjects (that is, the cells of the same
subject should be placed consecutively) before using as the input to the
`nebula` function. If the cells are not grouped, the `group_cell`
function can be used to first reorder the cells, as shwon below. If a
scaling factor is specified by the user, it should also be included in
`group_cell`. If the cells are already grouped, `group_cell` will return
*NULL*.

### Example

``` r
data_g = group_cell(count=sample_data$count,id=sample_data$sid,pred=df)
re = nebula(data_g$count,data_g$id,pred=data_g$pred)
```

If `pred` is not specified, `nebula` will fit the model with an intecept
term by default. This can be used when only the overdispersions are of
interest.

Specifying scaling factors
--------------------------

The scaling factor for each cell is specified in `nebula` using the
argument `offset`. The argument `offset` has to be a positive vector of
length *N*. Note that log(`offset`) will be the offset term in the NBMM.
If not specified, `nebula` will set `offset` as 1 by default, which
means that each cell is treated equally. Common scaling factors include
the library size of a cell or a normalizing factor adjusted using e.g.,
TMM.

### Example

``` r
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset)
```

Selection between NEBULA-LN and NEBULA-HL
-----------------------------------------

In *nebula*, a user can choose one of the two algorithms to fit an NBMM.
NEBULA-LN uses an approximated likelihood based on the law of large
numbers, and NEBULA-HL uses an h-likelihood. A user can select these
methods through `method='LN'` or `method='HL'`. NEBULA-LN is faster and
performs particularly well when the number of cells per subject (CPS) is
large. In the following analysis of the example data set comprising
\~200 cells per subject, the difference of the estimated cell-level
overdispersions between NEBULA-LN and NEBULA-HL is \~5% for most genes.

``` r
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN')
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
re_hl = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='HL')
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
## compare the estimated overdispersions
cbind(re_hl$overdispersion,re_ln$overdispersion)
#>       Subject      Cell    Subject      Cell
#> 1  0.08432321 0.9284703 0.08125256 0.8840821
#> 2  0.07455464 0.9726512 0.07102681 0.9255032
#> 3  0.17403276 0.9817570 0.17159404 0.9266395
#> 4  0.05352148 0.8516682 0.05026165 0.8124118
#> 5  0.07480033 1.3254379 0.07075366 1.2674146
#> 6  0.12372426 1.1653128 0.12086392 1.1096065
#> 7  0.07724824 0.9578169 0.07360445 0.9112956
#> 8  0.13797646 0.7991954 0.13571262 0.7549629
#> 9  0.05879492 0.8568854 0.05541398 0.8139652
#> 10 0.09782335 0.9940223 0.09496649 0.9410035
```

Such difference has little impact on testing fixed-effects predictors
under this sample size.

``` r
## compare the p-values for testing the predictors using NEBULA-LN and NEBULA-HL
cbind(re_hl$summary[,10:12],re_ln$summary[,10:12])
#>         p_X1      p_X2 p_cccontrol      p_X1      p_X2 p_cccontrol
#> 1  0.6373037 0.1346298   0.4950795 0.6354810 0.1291514   0.4919443
#> 2  0.9444825 0.3977109   0.7626827 0.9436079 0.3896819   0.7627706
#> 3  0.6282384 0.9787882   0.5087304 0.6271261 0.9792875   0.5058082
#> 4  0.8786074 0.6278826   0.2868256 0.8777381 0.6213846   0.2861434
#> 5  0.7596198 0.6872259   0.6544751 0.7579977 0.6795995   0.6537089
#> 6  0.7134192 0.8656686   0.6576835 0.7098168 0.8639067   0.6558008
#> 7  0.9216994 0.2230964   0.8977251 0.9225364 0.2151043   0.8987718
#> 8  0.7017083 0.4443604   0.3955343 0.7009654 0.4409831   0.3949916
#> 9  0.6505414 0.4561469   0.7238323 0.6489358 0.4473406   0.7209245
#> 10 0.4199828 0.7510837   0.7308108 0.4183419 0.7476005   0.7293432
```

The bias of NEBULA-LN in estimating the cell-level overdispersion gets
larger when the CPS value becomes lower or the gene expression is more
sparse. If the CPS value is \<30, `nebula` will set `method='HL'`
regardless of the user’s input. In contrast, NEBULA-HL is slower, but
its accuracy of estimating the overdispersions depends less on these
factors.

When NEBULA-LN is used, the user can opt for better accuracy of
estimating a smaller subject-level overdispersion through the argument
*κ*. NEBULA first fits the data using NEBULA-LN. If the estimated *κ*
for a gene is smaller than the user-defined value, NEBULA-HL will be
used to estimate the subject-level overdispersion for the gene. The
default value of *κ* is 800, which can provide a good estimate of the
subject-level overdispersion as low as \~0.005. Our simulation results
suggest that *κ* = 200 is often sufficent for achieving a well
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

Filtering low-expressed genes
-----------------------------

NEBULA-HL automatically uses a higher-order Laplace approximation for
low-expressed genes of which the average count per subject is less than
3. The higher-order Laplace approximation substantailly increases the
accuracy for estimating the subject-level overdispersion for
low-expressed genes and controls the false positive rate. Nevertheless,
we recommend removing genes with very low expression from the analysis
because there is little statistical power for these genes. Filtering out
low-expressed genes can be specified by `cpc=0.005` (i.e., counts per
cell\<0.5%). The argument `cpc` is defined by the ratio between the
total count of the gene and the number of cells.

Checking convergence for the summary statistics
-----------------------------------------------

*nebula* reports convergence information about the estimation algorithm
for each gene along with the summary statistics. This is useful and
important information for quality control to filter out genes of which
the estimation procedure potentially does not converge. Generally, a
convergence code \<= -20 suggests that the algorithm does not converge
well. If the convergence code is -30, which indicates a failure of
convergence, their summary statistics should NOT be used. If the
convergence code is -20 or -40, it indicates that the optimization
algorithm stops at the maximum step limit before the complete
convergence. The results should be interpreted with caution in this
case. The failure of convergence may occur when the sample size is very
small, there are too few positive counts, or the gene has huge
overdispersions, in which case the likelihood is flat or the
optimization is sensitive to the initial values. For those genes that
have a bad convergence code, in many cases, trying a different negative
binomial mixed model (e.g., NBLMM, see below for more details) may solve
the problem.

Using other mixed models
------------------------

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
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,model='PMM')
#> Remove  0  genes having low expression.
#> Analyzing  10  genes with  30  subjects and  6176  cells.
```

|  logFC\_(Intercept)|   logFC\_X1|   logFC\_X2|  logFC\_cccontrol|  se\_(Intercept)|     se\_X1|     se\_X2|  se\_cccontrol|  p\_(Intercept)|      p\_X1|      p\_X2|  p\_cccontrol|  gene\_id| gene |
|-------------------:|-----------:|-----------:|-----------------:|----------------:|----------:|----------:|--------------:|---------------:|----------:|----------:|-------------:|---------:|:-----|
|           -1.903571|  -0.0155809|  -0.0976660|         0.0511060|        0.0661297|  0.0329115|  0.0655553|      0.0642299|               0|  0.6359142|  0.1362700|     0.4262222|         1| A    |
|           -2.047864|  -0.0032670|  -0.0536887|        -0.0189269|        0.0644332|  0.0355074|  0.0635450|      0.0694853|               0|  0.9266904|  0.3981703|     0.7853239|         2| B    |
|           -2.032645|   0.0179777|   0.0009387|        -0.0505390|        0.0908196|  0.0345496|  0.0932449|      0.0676706|               0|  0.6028248|  0.9919678|     0.4551611|         3| C    |
|           -2.009746|  -0.0054963|  -0.0278602|         0.0782074|        0.0573209|  0.0350745|  0.0574459|      0.0686939|               0|  0.8754792|  0.6276888|     0.2549156|         4| D    |
|           -1.980528|   0.0106338|  -0.0248791|         0.0312190|        0.0644287|  0.0343355|  0.0621576|      0.0671645|               0|  0.7567865|  0.6889656|     0.6420644|         5| E    |
|           -1.950451|   0.0160341|  -0.0134775|        -0.0345244|        0.0778198|  0.0333858|  0.0738508|      0.0650410|               0|  0.6310363|  0.8551928|     0.5955505|         6| F    |
|           -1.970271|  -0.0026753|   0.0750060|        -0.0063677|        0.0645989|  0.0341936|  0.0615160|      0.0668723|               0|  0.9376369|  0.2227329|     0.9241391|         7| G    |
|           -1.964311|   0.0141532|  -0.0610984|        -0.0578672|        0.0809943|  0.0336579|  0.0800990|      0.0656800|               0|  0.6741201|  0.4455910|     0.3782927|         8| H    |
|           -2.074031|  -0.0178190|  -0.0436094|         0.0259745|        0.0597947|  0.0362203|  0.0587679|      0.0707912|               0|  0.6227459|  0.4580494|     0.7136813|         9| I    |
|           -2.046055|   0.0307026|   0.0227238|        -0.0246112|        0.0714158|  0.0354844|  0.0702268|      0.0691813|               0|  0.3869068|  0.7462578|     0.7220276|        10| J    |

Testing contrasts
-----------------

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
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',covariance=TRUE)
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
this example, we want to test whether the log(FC) of *X1* and
*cccontrol* are equal for the first gene. This hypothesis leads to the
contrast vector `(0 1 -1 0)`. Thus, the test can be performed using the
following code.

``` r
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
## the gene to test
gene_i = 1
## output covariance
re_ln = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset,method='LN',covariance=TRUE)
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
