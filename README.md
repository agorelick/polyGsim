# polyGsim
Simulate random multi-sample coalescence trees and corresponding poly-guanine genotypes

--- 

```r
## Fig 1. Simulate a random phylogeny of primary tumor clones and (monoclonal) metastases.

library(polyGsim)
set.seed(100) # reproducable

n_markers <- 58
mu_indel <- 1e-4
mu_cna <- 0.001

tree <- get_clone_tree(n_primary_clones=4, 
                       n_is_avg=F, 
                       met_prob_per_clone=0.2, 
                       avg_mets_per_clone=3)
plot(tree)
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig1.png" width="500" /> </p>

**Fig 1.** A random phylogeny of primary tumor clones and monoclonal metastases. Specified parameters: Exactly 4 primary clones; for each primary clone there is prop=0.2 of a metastatic clade somewhere on the tree; each met-clade has an average of 3 monoclonal metastases.

---

```r
## Fig 2. Get a random number of generations from 1st cancer cell to last tumor sample taken/patient death. 
## Then plot the true (un-knowable) temporal evolution of the clones and metastases.

max_gens <- random_generations(bdratio=1.01)$generations
plot_chronology(tree,'Simulated clone/metastasis evolution',max_gens)
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig2.png" width="500" /> </p>

**Fig 2.** Simulated temporal evolution of clones and monoclonal metastases from the first cancer cell. The total number of generations (span of x-axis) is derived from a discrete-time birth-death branching process with ratio of birth-rate to death-rate (bdratio) = 1.01. Branch lengths correspond to exactly number of cell divisions (coalescence times labeled in blue). Primary tumor samples labeled in green (P1-4), metastases in gold (M1-4), normal (germline) in black (N).

---

```r
## Get allele-specific lengths for each poly-G marker in each clone given indels and CNAs.
## Then plot a heatmap of the marker lengths for each copy in each clone/met.
## Here we use an indel-rate of 1e-4 indels/marker/cell division (with equally likely insertions and deletions).
## We similarly use SCNA rate of 1e-4 copy number deletions or duplications/marker/cell division (at most 4 copies allowed)
gt <- get_genotypes_with_cnas(tree, 
                              n_markers=58, 
                              mu.indel=1e-4, 
                              mu.cna=1e-4,
                              gens_until_first_cancer_cell=1e4)
plot_genotypes(gt)
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig3.png" width="700" /> </p>

**Fig 3.** Simulated genotypes (poly-guanine repeat lengths) for each allelic copy in each sample. Note that Normal always has only 2 copies (diploid). Grey regions have missing values, either due to copy number deletions (in copies 1, 2) or due to no duplication (in copies 3, 4). Each marker has a 50% probability of having copy number alterations early (i.e. after first cancer cell but *before subsequent indels* in all clones) or late (after first cancer cell and *after subsequent indels* in all clones).

---

```r
## get random purities and clonalities
pc <- get_purity_and_clonality(tree)

## for now, replace the default random clonalities to make all samples 80% clonal
pc$clonality <- 0.8 

## get mixing fractions. Uniform mixing assumes the 20% non-clonal cancer cells are evenly from the remaining clones/mets.
mix <- get_mixing_proportions(pc,nonclonal_evenly_distributed=T)
plot_mixtures(mix,'Random purity, 80% clonality, even prop from other clones/mets')
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig4.png" width="600" /> </p>

**Fig 4.** Each clone is simulated to be a mixture of a clonal population and other non-clonal populations of tumor cells and admixed normal cells. These  correspond to simulated 'samples' of each tumor, and will be used to generate 'observed' genotypes for the admixed tumors. Note that **clonality** is the proportion of cancer cells in the clone derived from itself rather than other clonal populations. By default, clonality is random (uniformly distributed between 0-1). In this example, it was set to 80% in each tumor sample. The remaining non-clonal fraction of cancer cells can either be random (uniformly distributed between 0-(1-clonality), or evenly distributed (as shown)).

---

```r
## get "observed" marker lengths, which reflect the true genotypes of each clone, 
## their proportions in the admixed tumor sample, 
## and their number of copies of each poly-G marker in each admixed cell population.

mu <- get_mean_marker_lengths(gt, mix, n_markers)
head(t(mu), 10) ## print the first 10 poly-G markers' mean lengths in each sample

          M1       M2       M3       M4       P1       P2       P3       P4 normal
m1  16.60913 16.64407 16.53508 16.87754 16.75657 16.84734 16.69669 16.60709   16.5
m2  23.00000 23.00000 23.00000 23.00000 23.00000 23.00000 23.00000 23.00000   23.0
m3  15.32108 15.28362 15.42861 15.59829 15.73231 15.14952 15.23775 15.77220   15.5
m4  17.00494 17.00620 17.00178 17.01219 17.00765 17.01159 17.22037 17.00328   17.0
m5  14.17281 14.21683 14.01500 14.56880 14.26774 14.40574 14.27547 14.11482   14.0
m6  19.84563 19.93367 19.62493 20.35365 20.03548 20.31148 20.05093 19.72964   19.5
m7  18.20841 18.26182 18.07509 18.51851 18.32376 18.49260 20.33198 18.13824   18.0
m8  18.17775 18.22303 18.06425 18.76828 18.27539 18.41733 18.28334 18.11810   18.0
m9  18.17281 18.21683 18.06247 18.42682 18.26774 18.40574 18.27547 18.11482   18.0
m10 21.57530 21.43931 21.86346 20.53065 21.20625 20.64819 21.23451 21.66870   22.0
```

**Table 1.** The mean length of the first 10 poly-G markers in each simulated tumor sample. Mean lengths are calculated fron the genotypes and number of copies of each marker in each cell population, and the proportion of their contribution to the admixed simulated sample. Note, in this example, marker m2 had no indels in any clone/met during the cancer evolution. This is evident in Fig 3, where the only variant in m2 is a copy-number deletion in P3 (which does not affect the mean marker length without any indels).
