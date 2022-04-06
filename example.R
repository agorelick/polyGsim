rm(list=ls())

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
png('figures/wiki/fig1.png')
plot(tree)
dev.off()

## plot the true (un-knowable) chronology of the clones
plot_chronology(tree,'Simulated clone/metastasis evolution',max_gens)
ggsave('../chronoloG/figures/wiki/fig2.png',width=7,height=5)

## get indels for the given tree
indels <- get_indels_for_tree(tree, max_gens, mu, n_markers) 

## get genotypes for each clone
gt <- get_clone_genotypes(indels, normal)


## get random purities and clonalities
pc <- get_purity_and_clonality(tree)
original_random_clonalities <- pc$clonality
pc$clonality <- 1

## get mixing fractions matrix
mix <- get_mixing_proportions(pc,even_mixing=F)
plot_mixtures(mix,'Random purity, 100% clonality') #, even prop from other clones/mets')
ggsave('../chronoloG/figures/wiki/fig3.png',width=7,height=5)


## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)
ang::write_distance_matrix( round(head(t(ml), 10),5), 'figures/wiki/table1.txt')

## anonymize marker lengths so that each marker's minimum value is 0.
ml <- get_anonymized_marker_lengths(ml)

## calculate the angular distance matrix from the anonymized marker lengths
ad <- get_angular_distance_matrix(ml)

## plot the angular distance neighbor joining tree
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,title='Angular distance tree')
ggsave('../chronoloG/figures/wiki/fig4.pdf',width=7,height=5)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# what if dominant clone was not at 100% clonality?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get random purities and clonalities
pc$clonality <- original_random_clonalities

## get mixing fractions matrix
mix <- get_mixing_proportions(pc,even_mixing=F)
plot_mixtures(mix,'Purity [0-100%], clonality [50-100%], random contribution from other clones')
ggsave('../chronoloG/figures/wiki/fig5.png',width=7,height=5)

## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)
ml <- get_anonymized_marker_lengths(ml)
ad <- get_angular_distance_matrix(ml)
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,title='Angular distance tree, mixed clones')
ggsave('../chronoloG/figures/wiki/fig6.png',width=7,height=5)





