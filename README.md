
---
title: "gatars version 0.2.25"
author: ''
output:
  html_document:
    css: ~/Documents/headings.css
    fig_caption: yes
    highlight: default
    number_sections: yes
    theme: readable
    toc: yes
  pdf_document:
    highlight: tango
    number_sections: yes
    pandoc_args: --variable=geometry:margin=0.75in
    toc: yes
  word_document:
    toc: yes
fontsize: 12pt
---
<A NAME="top"> </A>


"2018-02-17 23:28:28 GMT"

# Introduction {#intro}

For comments or questions, please contact Gail Gong at gailgongster@gmail.com

To download the RStudio source file for this webpage, click here: [gatars-follow-me.Rmd](gatars-follow-me.Rmd)

`gatars` ("Genetic Association Tests for Arbitrarily Related Subjects") tests the association between a specified set of genetic markers, called target markers, and a binary or quantitative trait using subjects with any genealogical relationship by calculating the P-values of the following seven statistics: the three basic statistics given by the squared burden $Q_B$, the SKAT $Q_S$ and the trait-based $Q_T$, and their four optimized linear combinations $Q_{BS}$, $Q_{BT}$, $Q_{ST}$, $Q_{BST}$.

Recall the notation from the manuscript: $G$ is the $N \times M$ genotype matrix of target markers, $y$ is the $N \times 1$ column vector of subjects’ coded trait phenotypes, and $\mu$ is the $N \times 1$ column vector of user-specified phenotype predictions. Any linear combination of the three basic statistics can be written as a quadratic form $Q(\alpha) = Z^T A(\alpha) Z$. Here $Z$ is a linear function of the two vectors $W G y$ and $W G \mu$, $W$ is a diagonal matrix whose user-specified diagonal components weight the target markers, $A(\alpha)$ is a square matrix that does not depend on $G$, and $\alpha = \big( \alpha_B, \alpha_S, \alpha_T \big)$  is a triple of coefficients that lie in the unit triangle shown in Figure 1 of the manuscript.  The limiting null distribution of $Z$ when $N$ is large is multivariate normal with mean and covariance matrix that involve the first two moments of $G$.

For any given fixed $\alpha$, the `davies` function from the `CompQuadForm` package can provide the significance level for $Q(\alpha)$, and therefore for any fixed linear combination of the three basic statistics. However the optimized statistics depend on data-driven $\alpha$'s, so the nominal P-values from `davies` do not reflect the their true significance levels.  To estimate the limiting null distribution of an optimized $Q$-statistic, `gatars` uses *genome resampling*:  Simulate  many null versions of the $Q$-statistic corresponding to a large number $K$ of simulated null genotype matrices $\tilde{G}^{(1)}, \cdots, \tilde{G}^{(K)}$.  Specifically, for each $k$, construct the $N \times M$ matrix $\tilde{G}^{(k)}$ by sampling one column from each *sampling set*, described below.  Then, for the $k$-th simulated matrix $\tilde{G}^{(k)}$, find the linear combination $\alpha^{k}$  having minimal nominal significance level $P_{\alpha^{k}}$ and store the transformed value $X_{\alpha^{k}} = - \text{log}_{10} \big( P_{\alpha^{k}} \big)$ . Repeat this procedure $K$ times and estimate the significance level of the observed $Q$-statistic to be the proportion of the $K$ transformed values $\Big\{X_{\alpha^k} \Big\}_{k = 1, \cdots, K}$ that exceed the observed transformed P-value corresponding to the observed $Q$-statistic.

To create sampling sets, `gatars` requires the $N$ subjects' genotypes at roughly 100,000 markers located throughout the autosomal genome. `gatars` chooses from these markers to create $M$ sampling sets, one for each target marker, so that the following two requirements are satisfied: (1) the sampling set for the $m$-th target marker contains markers all of whose minor allele frequencies (MAFs) match the MAF of the $m$-th target marker, and (2)  the markers in the sampling sets are independent of all target markers and all markers known to be trait-associated. (A MAF "matches" $\pi_m$  if it is within $\pi_m \times \big[1 -$ `epsilon` $, 1 +$ `epsilon` $\big]$.)

To perform its calculations `gatars` requires the following input from the user:
1) The binary or quantitative phenotype data $y_n$, $n = 1, \cdots, N$, of the $N$ subjects.
2) A user-specified null trait prediction $\mu_n$ (for example $\mu_n$ can be the mean of the $N$ $y_n$ values or the value predicted from a logistic or linear regression of $y_n$ on nongenetic covariates and/or on principal components of ancestry).
3) Coded genotypes $g_n = \big(g_{n1}, \cdots, g_{nM} \big)$ at $M$ target markers.
4) The pairwise inter-personal correlation coefficients summarized by the $N \times N$ matrix $\Psi$ in the
manuscript.
5) The $N$ subjects' genotypes at roughly 100,000 markers located throughout the autosomal genome.
6) The positions of any autosomal markers known to be trait-associated.

[TOP](#top)

# Installing `gatars` {#install-gatars}

If you have not yet installed R, download the latest version (at least 3.2.4) for your operating system at:

http://www.r-project.org

Run the R application. To install the `gatars` package, enter the following lines of code to the R prompt:

```
install.packages("devtools")
library(devtools)
install_github("gailg/gatars")
if("gatars" %in% rownames(installed.packages())){
  print("gatars installed successfully--you are good to go!")
} else {
  print("something went wrong--ask for help")
}
```
If your installation was successful, you should see the message

```
[1] "gatars installed successfully--you are good to go!"
```

[TOP](#top)

# Data required {#data-required}

`gatars` requires quite large amounts of data, and typically, which can be organized in Plink (http://pngu.mgh.harvard.edu/~purcell/plink/).  We assume that you can bring your data into R and create the following six data sets.

## `bim`

A `data.frame` containing the three columns `chromosome`, `snp`, `bp`, andcontaining $L$ rows corresponding to the $M$ target markers and those used to build the sampling sets.  The $l$-th row of `bim` summarizes the $l$-th marker and corresponds to the $l$-th column of `genotype` (described below). The column labeled `chromosome` contains integers between $1$ and $22$ (other integers may be included, but only chromosomes $1$ through $22$ will be used), the `snp` column contains character strings that name the marker, and the `bp` column contains integers giving the markers' bp positions.  (Because the `gatars` data set `hotspot` (described below) is given in Build hg38/GRCh38, `bp` must also be expressed in Build hg38/GRCh38.) 

## `exclusion_region`

A `data.frame` with one or several rows of the three columns `chromosome`, `start`, and `end`.  Each row of `exclusion_region` reflects one contiguous genomic region known to  be trait-associated and therefore a region used by `gatars` when creating the sampling sets. The column `chromosome` is an integer between `1` and `22` identifying the autosomal chromosome containing the region, and `start` and `end` describe its starting and ending  positions (bp).  If the region consists of a single marker, then `start` and `end` are both equal  to the position of this marker.  `start` and `end` must be expressed in Build hg38/GRCh38.

In addition, `gatars` provides the data set `hotspot` for your convenience. `hotspot contains the columns `chromosome`, `start`, and `end`; `chromosome` gives the  chromosome (an integer between $1$ and $22$) containing a given recombination hotspot;  `start` and `end` give the extent of the recombination hotspot (integers reflecting their Build hg38/GRCh38 bp positions).

## `genotype`

An $N \times L$ matrix, whose $(n,l)$-th element records the $n$-th subject's coded genotype at the $l$-th of $L$ markers. The $l$-th column of `genotype` corresponds to the $l$-th row of `bim`.  The object `genotype` could be obtained by reading in the `.bed` file from plink.  (Here we distinguish `genotype` (the $N \times L$ matrix whose columns correspond to the $M$ target markers AND the additional $L - M$ sampling set marker candidates) from the $N \times M$ matrix $G$, described in the manuscript, containing just the $M$ target markers).

## `phenotype`

A `data.frame` containing $N$ rows and the two columns `y` and `mu` where `y` contains the subjects' coded phenotypes (either `0` or `1` if the trait is binary or a real number if the trait is quantitative) and where `mu` contains the user-specified null trait prediction (a real number) for $y$.

## `Psi` = $\Psi$

`Psi` is the $N \times N$ matrix of user-specified null correlation coefficients between pairs of subjects at any marker.  `Psi` can be estimated using known family pedigree structures and/or from the subjects' genotypes at markers independent of those in the target set in software packages such as KING (Manichaikul et al 2010).




## `target_markers`

A character vector of length $M$ that is a subset of the column `bim$snp`.  This vector identifies the target markers.

[TOP](#top)

# An example data set {#example-data-set}

For illustration I will now use the example data set provided by `gatars`. You can access this example data with the following commands.


```r
library(gatars)
# Preparing the data
bim = gatars_example$bim
exclusion_region = gatars_example$exclusion_region
genotype = gatars_example$genotype
phenotype = gatars_example$phenotype
Psi = gatars_example$Psi
target_markers = gatars_example$target_markers[3:5]
```

## `bim` has 24578 rows and includes the columns `chromosome`, `snp`, and  `bp`


```r
str(bim)
```

```
## 'data.frame':	24578 obs. of  3 variables:
##  $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ snp       : chr  "exm94" "exm210" "exm269" "exm328" ...
##  $ bp        : int  874501 879481 881918 887967 949491 977356 978694 979748 980552 982722 ...
```

## `genotype` has 200 rows and 24578 columns

I chose a relatively small number of columns to keep the example small.


```r
str(genotype)
```

```
##  int [1:200, 1:24578] 0 0 0 1 0 0 0 0 0 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:200] "1" "2" "3" "4" ...
##   ..$ : chr [1:24578] "exm94" "exm210" "exm269" "exm328" ...
```

The column names of `genotype` must match `bim$snp`:


```r
all(colnames(genotype) == bim$snp)
```

```
## [1] TRUE
```

## `target_markers` is a character vector containing 3 names


```r
str(target_markers)
```

```
##  chr [1:3] "exm1061853" "exm1061861" "exm1061863"
```

The elements in `target_markers` must be included in `bim$snp`


```r
all(target_markers %in% bim$snp)
```

```
## [1] TRUE
```

## `phenotype` has 200 rows and includes the columns `y` and `mu`


```r
str(phenotype)
```

```
## 'data.frame':	200 obs. of  3 variables:
##  $ id: int  1 2 3 4 5 6 7 8 9 10 ...
##  $ y : int  0 1 1 1 1 1 1 1 1 1 ...
##  $ mu: num  0.879 0.838 0.798 0.838 0.798 ...
```

## `Psi` is a 200 $\times$ 200  square matrix


```r
str(Psi)
```

```
##  num [1:200, 1:200] 1 0.5 0 0 0 0 0 0 0 0 ...
```

The following figure reflects the fact that the first `100` people are made up of `50` sib pairs, and the last `100` people are independent.


```r
# figure to illustrate Psi
NNN = nrow(phenotype)
first_ten = 1:10
last_ten = NNN - (9:0)
matrix_image_fn(Psi[c(first_ten, last_ten), c(first_ten, last_ten)],
                main = "First and last 10 rows and columns of Psi")
```

<img src="index_files/figure-html/Psi_picture-1.png" width="480" />


## `exclusion_region` is a `data.frame` containing  the columns `chromosome`, `start`, and `end`


```r
str(exclusion_region)
```

```
## 'data.frame':	115 obs. of  3 variables:
##  $ chromosome: num  1 1 1 1 1 2 2 2 2 2 ...
##  $ start     : int  10496040 150685811 154861707 204549714 205788696 9977740 10570604 20688505 43326810 62904596 ...
##  $ end       : int  10496040 150685811 154861707 204549714 205788696 9977740 10570604 20688505 43326810 62904596 ...
```

## `hotspot` is available as soon as you enter `gatars`:


```r
str(hotspot)
```

```
## 'data.frame':	25657 obs. of  4 variables:
##  $ chromosome: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ center    : num  949794 1725398 2021339 2102281 2282281 ...
##  $ start     : int  635391 1722898 1950897 2099781 2273781 2387781 2462781 2835085 2881085 2956085 ...
##  $ end       : int  1264196 1727898 2091781 2104781 2290781 2437781 2493781 2841085 2902085 2971086 ...
```


# Check that your genotype matrix has rank $M$

`genotype_target_markers` is the $N \times M$ genotype matrix.


```r
# Checking the rank of the genotype_target_markers matrix
library(Matrix)
genotype_target_markers = genotype[, target_markers]
list(target_markers = target_markers,
     rank = as.numeric(rankMatrix( genotype_target_markers)))
```

```
## $target_markers
## [1] "exm1061853" "exm1061861" "exm1061863"
## 
## $rank
## [1] 3
```

# gatars-sampling-set {#gatars-sampling-set}

Once you have your data in R, and you have checked that the column rank of your genotype matrix $G =$ `genotype_target_markers` is $M$, use the function `gatars_sampling_set`  to create your sampling sets. 

## The arguments of `gatars_sampling_set`

The following is the function header of `gatars_sampling_set` showing which arguments are needed and in which order.

```
gatars_sampling_set(
  bim,
  epsilon,
  exclusion_region,
  genotype,
  hotspot,
  target_markers
)
```
### `bim`,  `exclusion_region`, `genotype`, `hotspot`, and `target_markers`,

These objects have already been discussed in [An example data set](#example-data-set).

### `epsilon` 

A positive small real number used to  used to determine if a marker can be included in a sampling set.  When creating the $m$-th sampling set for target marker with minor allele frequency $\pi_m$, only those markers whose minor allele frequencies fall within the interval
$\pi_m \times \big[1 -$ `epsilon` $, 1 +$ `epsilon` $\big]$ can be included in the sampling set.

## Example calls to `gatars_sampling_set` 


```r
# Example calls to gatars_sampling_set
set.seed(42)
epsilon = 0.01
exclusion_region = NULL
sampling_set = gatars_sampling_set(
  bim,
  epsilon,
  exclusion_region,
  genotype,
  hotspot,
  target_markers
)
print(sampling_set)
```

```
## $sampling_set_report
##      min     pi    max set_size
## 1 0.0200 0.0200 0.0200      632
## 2 0.0525 0.0525 0.0525      687
## 3 0.0100 0.0100 0.0100      627
## 
## $minimum_sampling_set_size
## [1] 627
```


```r
previous_sampling_set = sampling_set
set.seed(42)
exclusion_region = gatars_example$exclusion_region
head(exclusion_region)
```

```
##    chromosome     start       end
## 1           1  10496040  10496040
## 2           1 150685811 150685811
## 3           1 154861707 154861707
## 4           1 204549714 204549714
## 5           1 205788696 205788696
## 11          2   9977740   9977740
```

```r
sampling_set = gatars_sampling_set(
  bim,
  epsilon,
  exclusion_region,
  genotype,
  hotspot,
  target_markers
)
print(previous_sampling_set)
```

```
## $sampling_set_report
##      min     pi    max set_size
## 1 0.0200 0.0200 0.0200      632
## 2 0.0525 0.0525 0.0525      687
## 3 0.0100 0.0100 0.0100      627
## 
## $minimum_sampling_set_size
## [1] 627
```

```r
print(sampling_set)
```

```
## $sampling_set_report
##      min     pi    max set_size
## 1 0.0200 0.0200 0.0200      527
## 2 0.0525 0.0525 0.0525      570
## 3 0.0100 0.0100 0.0100      509
## 
## $minimum_sampling_set_size
## [1] 509
```

[TOP](#top)

# `gatars_test_size`

## The arguments of `gatars_test_size`

The following is the function header of `gatars_test_size` showing which arguments are needed and in which order.

```
gatars_test_size(phenotype, Psi, sampling_set, N_simulated_nulls, weights)
```

### `phenotype` and `Psi` 

These objects have already been discussed in [An example data set](#example-data-set).

### `sampling_set` 

This object is the result of `gatars_sampling_set` and you may follow the example from [gatars-sampling-set](#gatars-sampling-set)

### `N_simulated_nulls =` $K$

A positive integer which specifies the number of simulated null genotype matrices to generate to determine the p-value of the optimized statistics.


### `weights`

A vector of length equal to $M$ the length of `target_markers`.  The entries in the vector are non-negative real numbers.  The size of the `m`-th entry reflects the importance of the `m`-th target marker. If `weights` is not specified, the `gatars_test_size` assumes you would like equal weights among the $M$ target markers.

[TOP](#top)

## Example calls to `gatars_test_size`


```r
# Calling gatars_test_size 
start = Sys.time()
set.seed(42)
gatars_test_size(phenotype, Psi, sampling_set, N_simulated_nulls = 100, weights = c(5, 3, 2))
```

```
## Warning in davies(-absolute_qqq, lambda = lambda, delta = delta, lim =
## 50000, : Consider playing with 'lim' or 'acc'.
```

```
## $N_simulated_nulls_required
## [1] 100
## 
## $p_value
##            B            S            T           BS           BT 
## 0.0068507014 0.0411525864 0.0106099841 0.0000000000 0.0000000000 
##           ST          BST 
## 0.0000000000 0.0000000000 
## 
## $q
##          B          S          T         BS         BT         ST 
##  48.244722  22.416655 160.533764  48.244722  61.752125  89.103421 
##        BST 
##  61.752125 
## 
## $x
##        BS        BT        ST       BST 
## 2.1642650 2.4065813 2.0275387 2.4065813
```

```r
elapsed_time = Sys.time() - start
paste("100 sim reps used", 
      round(elapsed_time, 1), attributes(elapsed_time)$units)
```

```
## [1] "100 sim reps used 12.3 secs"
```

# All the code in one place


```r
library(gatars)
# Preparing the data
bim = gatars_example$bim
exclusion_region = gatars_example$exclusion_region
genotype = gatars_example$genotype
phenotype = gatars_example$phenotype
Psi = gatars_example$Psi
target_markers = gatars_example$target_markers[3:5]

# figure to illustrate Psi
NNN = nrow(phenotype)
first_ten = 1:10
last_ten = NNN - (9:0)
matrix_image_fn(Psi[c(first_ten, last_ten), c(first_ten, last_ten)],
                main = "First and last 10 rows and columns of Psi")
```

<img src="index_files/figure-html/unnamed-chunk-15-1.png" width="768" />

```r
# Checking the rank of the genotype_target_markers matrix
library(Matrix)
genotype_target_markers = genotype[, target_markers]
list(target_markers = target_markers,
     rank = as.numeric(rankMatrix( genotype_target_markers)))
```

```
## $target_markers
## [1] "exm1061853" "exm1061861" "exm1061863"
## 
## $rank
## [1] 3
```

```r
# Example calls to gatars_sampling_set
set.seed(42)
epsilon = 0.01
exclusion_region = NULL
sampling_set = gatars_sampling_set(
  bim,
  epsilon,
  exclusion_region,
  genotype,
  hotspot,
  target_markers
)
print(sampling_set)
```

```
## $sampling_set_report
##      min     pi    max set_size
## 1 0.0200 0.0200 0.0200      632
## 2 0.0525 0.0525 0.0525      687
## 3 0.0100 0.0100 0.0100      627
## 
## $minimum_sampling_set_size
## [1] 627
```

```r
previous_sampling_set = sampling_set
set.seed(42)
exclusion_region = gatars_example$exclusion_region
head(exclusion_region)
```

```
##    chromosome     start       end
## 1           1  10496040  10496040
## 2           1 150685811 150685811
## 3           1 154861707 154861707
## 4           1 204549714 204549714
## 5           1 205788696 205788696
## 11          2   9977740   9977740
```

```r
sampling_set = gatars_sampling_set(
  bim,
  epsilon,
  exclusion_region,
  genotype,
  hotspot,
  target_markers
)
print(previous_sampling_set)
```

```
## $sampling_set_report
##      min     pi    max set_size
## 1 0.0200 0.0200 0.0200      632
## 2 0.0525 0.0525 0.0525      687
## 3 0.0100 0.0100 0.0100      627
## 
## $minimum_sampling_set_size
## [1] 627
```

```r
print(sampling_set)
```

```
## $sampling_set_report
##      min     pi    max set_size
## 1 0.0200 0.0200 0.0200      527
## 2 0.0525 0.0525 0.0525      570
## 3 0.0100 0.0100 0.0100      509
## 
## $minimum_sampling_set_size
## [1] 509
```

```r
# Calling gatars_test_size 
start = Sys.time()
set.seed(42)
gatars_test_size(phenotype, Psi, sampling_set, N_simulated_nulls = 100, weights = c(5, 3, 2))
```

```
## $N_simulated_nulls_required
## [1] 100
## 
## $p_value
##            B            S            T           BS           BT 
## 0.0068507014 0.0411525864 0.0106099841 0.0000000000 0.0000000000 
##           ST          BST 
## 0.0000000000 0.0000000000 
## 
## $q
##          B          S          T         BS         BT         ST 
##  48.244722  22.416655 160.533764  48.244722  61.752125  89.103421 
##        BST 
##  61.752125 
## 
## $x
##        BS        BT        ST       BST 
## 2.1642650 2.4065813 2.0275387 2.4065813
```

```r
elapsed_time = Sys.time() - start
paste("100 sim reps used", 
      round(elapsed_time, 1), attributes(elapsed_time)$units)
```

```
## [1] "100 sim reps used 12.3 secs"
```

```r
start = Sys.time()
set.seed(42)
gatars_test_size(phenotype, Psi, sampling_set, N_simulated_nulls = 1000, weights = c(5, 3, 2))
```

```
## $N_simulated_nulls_required
## [1] 1000
## 
## $p_value
##            B            S            T           BS           BT 
## 0.0068507014 0.0411525864 0.0106099841 0.0070000000 0.0040000000 
##           ST          BST 
## 0.0130000000 0.0050000000 
## 
## $q
##          B          S          T         BS         BT         ST 
##  48.244722  22.416655 160.533764  48.244722  61.752125  89.103421 
##        BST 
##  61.752125 
## 
## $x
##        BS        BT        ST       BST 
## 2.1642650 2.4065813 2.0275387 2.4065813
```

```r
elapsed_time = Sys.time() - start
paste("1000 sim reps used", 
      round(elapsed_time, 1), attributes(elapsed_time)$units)
```

```
## [1] "1000 sim reps used 2.1 mins"
```

```r
start = Sys.time()
set.seed(42)
gatars_test_size(phenotype, Psi, sampling_set, N_simulated_nulls = 2000, weights = c(5, 3, 2))
```

```
## $N_simulated_nulls_required
## [1] 2000
## 
## $p_value
##            B            S            T           BS           BT 
## 0.0068507014 0.0411525864 0.0106099841 0.0095000000 0.0040000000 
##           ST          BST 
## 0.0120000000 0.0075000000 
## 
## $q
##          B          S          T         BS         BT         ST 
##  48.244722  22.416655 160.533764  48.244722  61.752125  89.103421 
##        BST 
##  61.752125 
## 
## $x
##        BS        BT        ST       BST 
## 2.1642650 2.4065813 2.0275387 2.4065813
```

```r
elapsed_time = Sys.time() - start
paste("2000 sim reps used", 
      round(elapsed_time, 1), attributes(elapsed_time)$units)
```

```
## [1] "2000 sim reps used 4.1 mins"
```










































