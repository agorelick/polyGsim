# polyGsim

Simulate random multi-sample coalescence trees and corresponding poly-guanine genotypes.

## Installation

```r
## install the devtools package if you have not already
install.packages('devtools')

## install the polyGsim package
devtools::install_github('https://github.com/agorelick/polyGsim')
```



## Simulation 1: Poly-G genotypes for impure bulk tissue samples, no cancer cell clone mixing.

The following steps will be used to simulate poly-G genotype data for a given phylogeny of tumor clones and monoclonal metastases. In reality, this phylogeny is not knowable and can only be inferred from the genotype data. In this example, we simulate random impurities in each bulk sample, which  introduced some proportion of non-cancer cells. This example is an ideal case where each sample represents a single cancer cell population. That is, there is no mixture of cancer cells from different clones/metastases.

```r
## Simulate a random phylogeny of primary tumor clones and (monoclonal) metastases.

library(polyGsim)

set.seed(100)

## define following parameters
n_primary_clones <- 3
n_is_avg <- F
met_prob_per_clone <- 0.2
avg_mets_per_clone <- 3
n_markers <- 58
mu <- 1e-4
bdratio <- 1.01

## generate random genotype for the average normal (i.e. zygote)
normal <- get_normal_genotype(n_markers=n_markers)

## get a random number of generations from zygote to last tumor sample taken/patient death
max_gens <- random_generations(bdratio=bdratio)$generations

## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
tree <- get_clone_tree(n_primary_clones=n_primary_clones,
                       n_is_avg=n_is_avg, 
                       met_prob_per_clone=met_prob_per_clone,
                       avg_mets_per_clone=avg_mets_per_clone)

## Fig 1.
plot(tree)
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig1.png" width="500" /> </p>

**Fig 1.** A random phylogeny of primary tumor clones and monoclonal metastases. Specified parameters: Exactly 3 primary clones; for each primary clone there is prop=0.2 of a metastatic clade somewhere on the tree (minimum of 1); each met-clade has an average of 3 monoclonal metastases (minimum of 3).


```r
## Fig 2. plot the true (un-knowable) chronology of the clones
plot_chronology(tree,'Simulated clone/metastasis evolution',max_gens)
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig2.png" width="500" /> </p>

**Fig 2.** Simulated temporal evolution of clones and monoclonal metastases from the first cancer cell. The total number of generations (span of x-axis) is derived from a discrete-time birth-death branching process with ratio of birth-rate to death-rate (bdratio) = 1.01. Branch lengths correspond to exactly number of cell divisions (coalescence times labeled in blue). Primary tumor samples labeled in green (P1-3), metastases in gold (M1-4), normal (germline) in black (N).


```r

## get indels for the given tree
indels <- get_indels_for_tree(tree, max_gens, mu, n_markers) 

## get genotypes for each clone
gt <- get_clone_genotypes(indels, normal)

## get random purities and clonalities
pc <- get_purity_and_clonality(tree)

## for now, we force bulk samples to have 100% clonality for their dominant clone
original_random_clonalities <- pc$clonality
pc$clonality <- 1

## get mixing fractions matrix
mix <- get_mixing_proportions(pc,even_mixing=F)

## Fig 3. Bulk samples with random purity and 100% clonality.
plot_mixtures(mix,'Random purity, 100% clonality')
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig3.png" width="700" /> </p>

**Fig 3.** Each bulk sample is generally simulated to be a mixture of a clonal cancer cell population and other non-clonal populations of tumor cells and admixed normal cells. In this example, we make each dominant clone 100% clonal (i.e. there is _no contribution_ of cells from other tumor clones or metastases). Note: **Purity** is the % of cancer-cells in the bulk sample, drawn from 0-100% and centered at 50%. **Clonality** is the proportion of cancer cells derived from the dominant clone rather than other clones/metastases in the patient, which is usually uniformly distributed between 50-100% (here, entirely 100%).

```r

## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)

## Table 1. The first 10 poly-G markers' mean lengths in each sample
head(t(ml), 10) 

          M1      M2       M3       M4       P1       P2       P3 normal
m1  19.50000 19.5000 19.50000 19.50000 19.50000 19.50000 19.50000   19.5
m2  20.00000 20.0000 20.00000 20.00000 20.00000 20.00000 20.27318   20.0
m3  15.00000 15.0000 15.00000 15.00000 15.00000 15.00000 15.00000   15.0
m4  21.50000 21.5000 21.50000 21.50000 21.50000 21.50000 21.50000   21.5
m5  15.00000 15.0000 15.00000 15.00000 15.00000 15.00000 15.00000   15.0
m6  18.25476 18.0401 18.27913 18.43612 18.32693 18.59848 18.54636   18.0
m7  16.74524 16.9599 16.72087 16.56388 17.00000 17.00000 17.00000   17.0
m8  20.00000 20.0000 19.72087 19.56388 20.00000 20.00000 20.27318   20.0
m9  20.50000 20.5000 20.50000 20.06388 20.50000 20.50000 20.50000   20.5
m10 15.74524 15.9599 15.72087 15.56388 15.67307 16.00000 16.00000   16.0
```

**Table 1.** The mean length of the first 10 poly-G markers in each simulated tumor sample. Mean lengths are calculated from the genotypes and number of copies of each marker in each cell population (currently all diploid), and the proportion of their contribution to the admixed simulated sample. 

```r
## anonymize marker lengths so that each marker's minimum value is 0.
ml <- get_anonymized_marker_lengths(ml)

## calculate the angular distance matrix from the anonymized marker lengths
ad <- get_angular_distance_matrix(ml)

## Fig 4. plot the angular distance neighbor joining tree
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,title='Angular distance tree')
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig4.png" width="700" /> </p>

**Fig 4.** Angular distance neighbor-joining tree for this simulation. Angular distance is generated from the mean poly-G marker lengths, assuming bulk samples consisting of a dominant cancer clone with 100% clonality, and overall impurity due to admixed normal cells.


## Simulation 2: Poly-G genotypes for impure bulk tissue samples with cancer cell clone mixing.

In reality, the cancer cell population in bulk tumor samples may be a mixture of a dominant clone and other subclonal populations derived from different regions in the primary tumor and from other metastases. The resulting within-sample heterogeneity will affect poly-G genotypes for bulk tumor samples. In this example, we modify the above simulation to allow each bulk sample to have a random clonality from 50-100%. The proportion of cancer cells not from the clonal population are drawn from the other cancer clones and metastasis (uniformly distributed), as shown in **Fig 5**. The angular distance tree in **Fig 6** now less clearly resembles the underlying phylogeny in **Fig 2**, but the major qualitative points are still evident. 

```r
## get random purities and clonalities
pc$clonality <- original_random_clonalities

## get mixing fractions matrix
mix <- get_mixing_proportions(pc,even_mixing=F)
plot_mixtures(mix,'Purity [0-100%], clonality [50-100%], random contribution from other clones')
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig5.png" width="700" /> </p>

**Fig 5.** Each bulk sample is now a mixture of a clonal cancer cell population and other non-clonal populations of tumor cells, as well as the original proportion of normal cells. 

```r
## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)
ml <- get_anonymized_marker_lengths(ml)
ad <- get_angular_distance_matrix(ml)
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,title='Angular distance tree, mixed clones')
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig6.png" width="700" /> </p>

**Fig 6.** The angular distance neighbor-joining tree representing a more realistic scenario of impure bulk samples containg a mixture of cells from different clones and metastases throughout the patient's cancer.
