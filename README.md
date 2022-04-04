# polyGsim
Simulate random multi-sample coalescence trees and corresponding poly-guanine genotypes

--- 



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

**Fig 1.** A random phylogeny of primary tumor clones and monoclonal metastases. Specified parameters: Exactly 3 primary clones; for each primary clone there is prop=0.2 of a metastatic clade somewhere on the tree; each met-clade has an average of 3 monoclonal metastases.

---

```r
## Fig 2. plot the true (un-knowable) chronology of the clones
plot_chronology(tree,'Simulated clone/metastasis evolution',max_gens)
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig2.png" width="500" /> </p>

**Fig 2.** Simulated temporal evolution of clones and monoclonal metastases from the first cancer cell. The total number of generations (span of x-axis) is derived from a discrete-time birth-death branching process with ratio of birth-rate to death-rate (bdratio) = 1.01. Branch lengths correspond to exactly number of cell divisions (coalescence times labeled in blue). Primary tumor samples labeled in green (P1-3), metastases in gold (M1-4), normal (germline) in black (N).

---

```r

## get indels for the given tree
indels <- get_indels_for_tree(tree, max_gens, mu, n_markers) 

## get genotypes for each clone
gt <- get_clone_genotypes(indels, normal)

## get random purities and clonalities
pc <- get_purity_and_clonality(tree)

## get mixing fractions matrix
mix <- get_mixing_proportions(pc,uniform_mixing=F)

## Fig 3. Bulk samples with random purity and clonality.
plot_mixtures(mix,'Purity [0-100%], clonality [50-100%], even prop from other clones/mets')

```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig3.png" width="700" /> </p>

**Fig 3.** Each bulk sample is simulated to be a mixture of a clonal cancer cell population and other non-clonal populations of tumor cells and admixed normal cells. These correspond to simulated 'samples' of each tumor, and will be used to generate the 'observed' genotypes for the admixed samples. **Purity** is the % of cancer-cells in the bulk sample, drawn from 0-100% and centered at 50%. **Clonality** is the proportion of cancer cells derived from the dominant clone rather than other clones/metastases in the patient, uniformly distributed between 50-100%. The contributions of the non-dominant clones and mets to remaining proportion of cancer cells is uniformly distributed. 

---

```r

## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)

## Table 1. The first 10 poly-G markers' mean lengths in each sample
head(t(ml), 10) 

          M1       M2       M3       M4       P1       P2       P3 normal
m1  19.50000 19.50000 19.50000 19.50000 19.50000 19.50000 19.50000   19.5
m2  20.02946 20.00098 20.01766 20.00375 20.00858 20.00736 20.24360   20.0
m3  15.00000 15.00000 15.00000 15.00000 15.00000 15.00000 15.00000   15.0
m4  21.50000 21.50000 21.50000 21.50000 21.50000 21.50000 21.50000   21.5
m5  15.00000 15.00000 15.00000 15.00000 15.00000 15.00000 15.00000   15.0
m6  18.30117 18.04219 18.32874 18.45045 18.35371 18.55528 18.51991   18.0
m7  16.79663 16.96242 16.78004 16.58190 16.87044 16.97419 16.97867   17.0
m8  20.00245 19.99904 19.82119 19.59900 19.94589 19.99109 20.23049   20.0
m9  20.49123 20.49920 20.46737 20.09539 20.46424 20.48931 20.49398   20.5
m10 15.79165 15.96200 15.77048 15.57822 15.69985 15.95680 15.97354   16.0

```
**Table 1.** The mean length of the first 10 poly-G markers in each simulated tumor sample. Mean lengths are calculated from the genotypes and number of copies of each marker in each cell population (currently all diploid), and the proportion of their contribution to the admixed simulated sample. 

---

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

**Fig 4.** Angular distance neighbor-joining tree for this simulation. Angular distance is generated from the mean poly-G marker lengths, assuming impure bulk samples consisting of a dominant clone (labeled), admixed normal cells and cells from other non-dominant tumor clones and metastases.
