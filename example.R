
library(chronoG)

set.seed(42)


## main function.
result <- chronoG(5,F)
gt <- result$gt


## (extra) simulate a cancer with 10^12 cells and return number of generations to reach this size.
tmp <- random_generations(starting_cells=1, s=0.004, k=3, max_gens=1e4, max_cells=1e12)
tmp


## (extra) generate a random coalescence tree where X is proportional to the number of generations the cancer has evolved since the first cancer cell.
tree <- sim_coalescence(5, F)
plot(tree)
nodelabels()


## (extra) internal function 
map <- get_lineage_map(tree)
map


## plot a genotype-based tree
dm <- dist(gt, method='euclidean')
tree2 <- upgma(dm)
tree2 <- nj(dm)
tree2 <- phytools::reroot(tree2, node.number=grep('normal',tree2$tip.label))
plot(tree2)













