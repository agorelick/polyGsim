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
tree <- get_clone_tree(n_primary_clones=n_primary_clones, 
                       n_is_avg=n_is_avg, 
                       met_prob_per_clone=met_prob_per_clone, 
                       avg_mets_per_clone=avg_mets_per_clone)
plot(tree) ## Fig 1
```

![alt text](https://github.com/agorelick/polyGsim/blob/main/figures/wiki/fig1.png "Fig 1")
Fig 1. A random phylogeny of primary tumor clones and monoclonal metastases. Specified parameters: Exactly 4 primary clones, each with prob=0.2 of seeding observed metastases, where an average of 3 monoclonal metastases are associated with each met-seeding primary tumor clone.
