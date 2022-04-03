# polyGsim
Simulate random multi-sample coalescence trees and corresponding poly-guanine genotypes

```r
library(polyGsim)
set.seed(100)

n_primary_clones <- 4
n_is_avg <- F
met_prob_per_clone <- 0.2
avg_mets_per_clone <- 3
n_markers <- 58
mu_indel <- 1e-4
mu_cna <- 0.001
bdratio <- 1.01

## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
tree <- get_clone_tree(n_primary_clones=4, 
                       n_is_avg=F, 
                       met_prob_per_clone=0.2, 
                       avg_mets_per_clone=3)
plot(tree) ## Fig 1
```

<p align="center"> <img src="https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig1.png" width="500" height="500" /> </p>

**Fig 1.** A random phylogeny of primary tumor clones and monoclonal metastases. Specified parameters: Exactly 4 primary clones; for each primary clone there is prop=0.2 of a metastatic clade somewhere on the tree; each met-clade has an average of 3 monoclonal metastases.
