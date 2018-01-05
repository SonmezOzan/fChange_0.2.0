R package fChange: Change in the Covariance Structure
================
Ozan Sonmez



-   [Introduction](#introduction)
-   [Partial Sum Covariance Function Estimate](#partial-sum-covariance-function-estimate)
-   [Testing changes in the eigenvalues](#testing-changes-in-the-eigenvalues)
    -   [Component-wise](#component-wise)
    -   [Joint](#joint)
    -   [Two step process](#two-step-process)
-   [Testing Changes in the trace](#testing-changes-in-the-trace)

Introduction
============

In this documentation, the updated R package `fChange` will be discussed with the focus on the changes in the covariance structure of the functional data, specially functional time series data. Previously the `fChange` library is explored in the direction of the mean function changes in <a href="https://github.com/SonmezOzan/fChangeExample/blob/master/fChange_vignette.md"> this documentation </a> which was basically implemented from our work <a href="http://onlinelibrary.wiley.com/wol1/doi/10.1111/rssb.12257/full"> dating and detecting structural breaks in functional data without dimension reduction </a>.

This document will focus on testing and detecting changes in the covariance structure of the functional data. Many methods employed in functional data analysis are based on functional principal component analysis (fPCA). These entail projecting the functional observations into a lower dimensional space defined by a few functional principal components that are computed from an empirical covariance operator in the hope that these projections can account for a large percentage of the variability in the sample, and then these lower dimensional projections are easier to work with for subsequent analysis. When data are obtained via randomized experiments, it is reasonable to assume that the covariance structure is homogeneous throughout the sample, and, in this case, the principal components computed from the empirical covariance operator are estimating some “population” principal components that are well known to be optimal as a dimension reduction basis.Often times, it may be unclear whether or not the covariance structure of the data changes during the sampling period. If the covariance does happen to change within the functional time series, fPCA based analysis may be misleading in several ways. Firstly, the basis computed from the sample covariance is not estimating the optimal basis to reduce the dimension of data. Secondly, practitioners frequently choose the number of principal components to analyze based on the “total variance explained” (TVE) approach, which relies on estimating the eigenvalues of the covariance operator. If these quantities are not stable, then too many, or more likely too few, principal components may be considered in the subsequent analysis.

Changes in the covariance might occur in various ways. Here we explore changes in the eigenvalues (both marginally and jointly), eigenfunctions (partially) and the trace of the covariance operator of the functional data. All of the forementioned analyses are based on the partial sum estimate of the covariance operator.

Partial Sum Covariance Function Estimate
========================================

``` r
partial_cov(fdobj, x = NULL)
```

``` r
fdata1 = fun_IID(n=100, nbasis = 21) # IID functional data
eval_90 = partial_cov(fdata1, x=0.9)$eigen_val # Partial Cov estimate when x=0.9
eval_95 = partial_cov(fdata1, x=0.95)$eigen_val # Partial Cov estimate when x=0.95
eval_full = partial_cov(fdata1)$eigen_val # Partial Cov estimate when x=1 (whole sample)
evals = data.frame(eval = c(eval_90, eval_95, eval_full), 
                   x = c(rep(0.9, 21), rep(0.95, 21), rep(1,21)),
                   d = rep(c(1:21),3))
evals$x = as.character(evals$x)
ggplot(data=evals, aes(x=d, y=eval, fill=x)) +
  geom_bar(stat="identity", position=position_dodge())
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-2-1.png) Note that when *x* = 1, the partial sum estimate is identical to the sample covariance estimate that is computed from `pca.fd()` from the `fda` library:

``` r
fPCA = pca.fd(fdata1, nharm = 21, centerfns = T)
eval_fPCA = fPCA$values
# eigenvalues from pca.fd
eval_fPCA
```

    ##  [1] 1.160787428 0.287057779 0.115608648 0.064173834 0.040153045
    ##  [6] 0.026072143 0.021147500 0.017804111 0.010694980 0.008734947
    ## [11] 0.007010943 0.006309423 0.005378429 0.004268792 0.003776888
    ## [16] 0.003154162 0.002892508 0.002470378 0.002216827 0.001612658
    ## [21] 0.001145838

``` r
# eigenvalues from partial_cov
eval_full
```

    ##  [1] 1.160787428 0.287057779 0.115608649 0.064173835 0.040153047
    ##  [6] 0.026072143 0.021147500 0.017804110 0.010694980 0.008734949
    ## [11] 0.007010942 0.006309425 0.005378422 0.004268794 0.003776891
    ## [16] 0.003154163 0.002892500 0.002470379 0.002216832 0.001612657
    ## [21] 0.001145838

Testing changes in the eigenvalues
==================================

One of the most importatnt changes that might occur in the covariance structure of the functional data is the changes in the eigenvalues of the covariance function. testing and detecting these changes might reveal very imortant information. One might be interested in a specific eigenvalue change, which we will refer that as the component wise eigenvalue change, or the focus of interest might also be on the changes in the eigenvalues jointly. Here we analyze these two cases seperatly.

Component-wise
--------------

``` r
eval_component(fdobj, component, h = 2, mean_change = FALSE, delta = 0.1, M = 1000)
```

The test rejects the null if the observed test statistics *T*<sub>*n*1</sub> is larger than the defined quantile of the distribution sup<sub>*δ* ≤ *x* ≤ 1</sub>*W*<sup>2</sup>(*x*), which is obtained by using `M = 1000` Monte-Carlo samples. Above, we have generated an independent functional data with no changes in it. Therefore, applying the above test yields non-significant p value, indicating that there is no change in the defined eigenvalue:

``` r
# testing changes in the first eigenvalue
eval_component(fdata1, 1)
```

    ## $change
    ## [1] 0.66
    ## 
    ## $pvalue
    ## [1] 0.592

``` r
# testing changes in the second eigenvalue
eval_component(fdata1, 2)
```

    ## $change
    ## [1] 0.43
    ## 
    ## $pvalue
    ## [1] 0.831

``` r
# testing changes in the third eigenvalue
eval_component(fdata1, 3)
```

    ## $change
    ## [1] 0.7
    ## 
    ## $pvalue
    ## [1] 0.443

To simulate a case where there is a change in the eigenvalue, we use the following data gerenation process: functional data of size *n* were generated using *D* = 21 Fourier basis functions *v*<sub>1</sub>, …, *v*<sub>*D*</sub> on the unit interval \[0, 1\]. Without loss of generality, the initial mean curve *μ* is assumed to be the zero function.

Therefore, picking *j* = 1 puts the change in the first component of the eigenvalue decay *σ*, which mimics the case of having the change sitting the leading eigenvalue, similarly having *j* = 2 mimics the case of inserting the change in the second leading eigenvalue, and so on. The sensitivity parameter *c* determines the magnitude of the signal or the size of the change that is being inserted. For example, having *c* = 1 is identical to the case of having no change (the null hypothesis).

``` r
D = 21
IID = function(n, Sigma = (1:D)^-1, theta = 0.5, c, j){
  basis = create.fourier.basis(nbasis = 21)
  Sigma1 = Sigma2 = Sigma
  Sigma2[j] = c*Sigma[j]
  x = n*theta
  dat = matrix(0, D, n)
  for (j in 1:x){
    dat[,j] = rnorm(D, 0, Sigma1)
  }
  for (j in (x+1):n){
    dat[,j] = rnorm(D, 0, Sigma2)
  }
  fd(dat, basis)
}
```

Here we generate 3 functional data sets, with no change (*c* = 1), change in the first eigenvalue with a smaller magnitude of change (*c* = 1.5) and a larger magnitude of change (*c* = 3).

``` r
set.seed(1234)
N = 100 # Sample Size
fdata_nochange = IID(n = N, c = 1, j = 1)
fdata_c15 = IID(n = N, c = 1.5, j = 1)
fdata_c3 = IID(n = N, c = 3, j = 1)
eval_component(fdata_nochange, 1)
```

    ## $change
    ## [1] 0.55
    ## 
    ## $pvalue
    ## [1] 0.107

``` r
eval_component(fdata_c15, 1)
```

    ## $change
    ## [1] 0.5
    ## 
    ## $pvalue
    ## [1] 0.011

``` r
eval_component(fdata_c3, 1)
```

    ## $change
    ## [1] 0.5
    ## 
    ## $pvalue
    ## [1] 0.001

Note that as the magnitude of the change *c* gets larger the estimated pvalue gets smaller, indicating a stronger significance for the change in the sepcified eigenvalue. While we put the change in the leading (first) eigenvalue, the remaining eigenvalues should have no changes in them. One can see that by applying the component wise change to the different eigenvalues:

``` r
eval_component(fdata_c3, component = 2)
```

    ## $change
    ## [1] 0.21
    ## 
    ## $pvalue
    ## [1] 0.157

``` r
eval_component(fdata_c15, component = 3)
```

    ## $change
    ## [1] 0.13
    ## 
    ## $pvalue
    ## [1] 0.297

``` r
eval_component(fdata_nochange, component = 4)
```

    ## $change
    ## [1] 0.55
    ## 
    ## $pvalue
    ## [1] 0.423

below plot is the covariance surfaces of the three functional data that is generated above:

``` r
LongRun(fdata_c3, h=0)$contour_plot
LongRun(fdata_c15, h=0)$contour_plot
LongRun(fdata_nochange, h=0)$contour_plot
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-9-1.png)![](Untitled_files/figure-markdown_github/unnamed-chunk-9-2.png)![](Untitled_files/figure-markdown_github/unnamed-chunk-9-3.png)

Joint
-----

The below function `eval_joint()` implemets the joint test for changes in the eigenvalues, where the critical values are computed via Monte-Carlo methods.

``` r
eval_joint(fdobj, d, h = 2, mean_change = FALSE, delta = 0.1, M = 1000)
```

The dimension `d` that the joint test is performed can be selected via total variation explained (TVE) using the function `pick_dim()`. Using the functional data generated previously where the change is only iserted in the first eigenvalue, the joint test is also able to detect the change:

``` r
Joint_Test = eval_joint(fdata_c3, d=1)
Joint_Test
```

    ## $change
    ## [1] 0.5
    ## 
    ## $pvalue
    ## [1] 0.001
    ## 
    ## $eval_before
    ##  [1] 0.9733844630 0.2888920503 0.0696324476 0.0617743931 0.0390480537
    ##  [6] 0.0257334529 0.0208714299 0.0134742350 0.0100685267 0.0085802667
    ## [11] 0.0071743380 0.0055826508 0.0047217868 0.0036357408 0.0031701382
    ## [16] 0.0024962806 0.0021585029 0.0015419946 0.0012349205 0.0011020781
    ## [21] 0.0008714167
    ## 
    ## $eval_after
    ##  [1] 7.3560591927 0.2254591741 0.0910930902 0.0490697015 0.0378175414
    ##  [6] 0.0202579637 0.0169353034 0.0145352668 0.0131276710 0.0111425061
    ## [11] 0.0070605978 0.0051740875 0.0046917159 0.0039029240 0.0025684111
    ## [16] 0.0022999164 0.0018594716 0.0011206153 0.0010483926 0.0009757906
    ## [21] 0.0006607394

Since we now that the change sits in the leading eigenvalue performing the joint test with *d* = 1 will capture the change, where the test comes out to be significant and the location of the change is estimated to be at the true break location that the artificial change is inserted. However, in practice, one will not know where which eigenvalues are changing beforehand, hence defining an appropriate *d* will be crucial for the accurate analysis, since takin small *d* might result in missing the later eigenvalue changes while taking large *d* when only the first few eigenvalues are changing, might result in unreliable analysis, because intriducing large *d* might weaken the signal by having to estimate large number of eigenvaues that do not contribute to the joint change.

``` r
sapply(1:5, function(d) eval_joint(fdata_c3, d=d)$pvalue)
```

    ## [1] 0.001 0.003 0.008 0.021 0.060

``` r
# TVE before the change
TVE_before = cumsum(Joint_Test$eval_before)/sum(Joint_Test$eval_before)
# TVE after the change
TVE_after = cumsum(Joint_Test$eval_after)/sum(Joint_Test$eval_after)
plot(TVE_before, type="o", col="black", main="TVE", ylab="TVE (%)", xlab="d")
points(TVE_after, type="o", col="red")
legend("bottomright", c("before", "after"), col=c("black", "red"), lty=c(1,1))
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-13-1.png)

Note that, in the above example, the change was stronly sitting in the leading eigenvalue, hence the component wise test did also perform pretty well when it was applied to the first component. The power of the joint test, however, comes more appearent when the change is not dominant in one component but rather distributed among mant eigenvalues. In that case the component wise tests will have difficulty detecting these suttle changes while the joint test should be able to detect these changes.

``` r
IID_j = function(n, Sigma = (1:D)^-1, theta = 0.5, c){
  basis = create.fourier.basis(nbasis = 21)
  Sigma1 = Sigma2 = Sigma
  Sigma2[1:3] = c*Sigma[1:3]
  x = n*theta
  dat = matrix(0, D, n)
  for (j in 1:x){
    dat[,j] = rnorm(D, 0, Sigma1)
  }
  for (j in (x+1):n){
    dat[,j] = rnorm(D, 0, Sigma2)
  }
  fd(dat, basis)
}
```

The above function will put the change in the first three eigenvalues whose magnitude is controlled by the sensitivity parameter *c*.

``` r
set.seed(9)
fdat = IID_j(n=100, c=1.25)
eval_joint(fdat, d=3)
```

    ## $change
    ## [1] 0.52
    ## 
    ## $pvalue
    ## [1] 0.037
    ## 
    ## $eval_before
    ##  [1] 0.8764217438 0.1904926859 0.0981095894 0.0674337704 0.0417321499
    ##  [6] 0.0263680925 0.0212646822 0.0142492475 0.0128514680 0.0077086790
    ## [11] 0.0069020211 0.0050693885 0.0037085689 0.0035044318 0.0029932693
    ## [16] 0.0024834659 0.0019617793 0.0014492663 0.0010838865 0.0009375873
    ## [21] 0.0006126148
    ## 
    ## $eval_after
    ##  [1] 1.5855390729 0.4712608533 0.1286102556 0.0620403342 0.0312867149
    ##  [6] 0.0213299867 0.0170935655 0.0128184746 0.0109664321 0.0103616865
    ## [11] 0.0079241232 0.0058745699 0.0043690025 0.0035432553 0.0030953315
    ## [16] 0.0028615543 0.0018049638 0.0015684305 0.0013477261 0.0009788566
    ## [21] 0.0005215008

``` r
sapply(1:3, function(d) eval_component(fdat,d)$pvalue)
```

    ## [1] 0.447 0.017 0.593

Note that while the joint test yields a significant result, the component-wise test results varies.

Two step process
----------------

In practice, one will not have the knowledge of which eigenvalue changes. Therefore, the component wise test seems unpractical. However, these test might be useful to examine which eigenvalues contributes to the joint change, once the joint change rejects the null hyphothesis of no joint eigenvalue change.

We propse a two step process, where the joint eigenvalue test is initially applied to see whether the eigevalues are changing jointly. If the test comes out to be significant, the natural follow up question is how each eigenvalue contributes to the joint test. Hence, the component-wise test can be employed to examine the contribution of individual eigenvalues. Note that one needs to take into account of multiple testing here.

Testing Changes in the trace
============================

Anther importatnt changes to detect is the total variation of the functional data, which is the trace of the covariance operator. Having a nonhomogenous functional variation might effect how certain confidence and prediction intervals are constructed.

The below function tests for the changes in the trace of the covariance operator of the functional data:

``` r
trace_change(fdobj, mean_change = FALSE, delta = 0.1, M = 1000)
```

Here is an example:

``` r
set.seed(12)
# Functional data with No change
funData_1 = IID(n = 100, j=1, c=1) 
# Functional data with change in leading eigenvalue
funData_2 = IID(n = 100, j=1, c=2) 
par(mfrow=c(1,2))
plot(funData_1[1:50], col="red", main="No change, c=1")
```

    ## [1] "done"

``` r
lines(funData_1[51:100], col="blue")
plot(funData_2, col="grey", main="Change with c=3")
```

    ## [1] "done"

``` r
lines(funData_2[51:100], col="blue")
lines(funData_2[1:50], col="red")
```

![](Untitled_files/figure-markdown_github/unnamed-chunk-17-1.png)

Clearly in the second plot, the variation is increased by the magnitude of 2 times the first eigenvalue. This additional variation is clear in the second plot, where the variation before the change (in red) is quite smaller than the variation after the change (in blue).

``` r
# testing changes in the total variation when there is no change
trace_change(funData_1)
```

    ## $change
    ## [1] 0.44
    ## 
    ## $pvalue
    ## [1] 0.516
    ## 
    ## $trace_before
    ## [1] 1.313784
    ## 
    ## $trace_after
    ## [1] 1.811749

``` r
# testing changes in the total variation when the change is inserted
# in the leading eigenvalue with multiplying it by c=2
trace_change(funData_2)
```

    ## $change
    ## [1] 0.5
    ## 
    ## $pvalue
    ## [1] 0
    ## 
    ## $trace_before
    ## [1] 1.524365
    ## 
    ## $trace_after
    ## [1] 5.001732
