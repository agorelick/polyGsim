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

