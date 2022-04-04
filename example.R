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
mix <- get_mixing_proportions(pc,uniform_mixing=F)
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
ggsave('../chronoloG/figures/wiki/fig4.png',width=7,height=5)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# what if dominant clone was not at 100% clonality?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get random purities and clonalities
pc$clonality <- original_random_clonalities

## get mixing fractions matrix
mix <- get_mixing_proportions(pc,uniform_mixing=F)
plot_mixtures(mix,'Purity [0-100%], clonality [50-100%], random contribution from other clones')
ggsave('../chronoloG/figures/wiki/fig5.png',width=7,height=5)

## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)
ml <- get_anonymized_marker_lengths(ml)
ad <- get_angular_distance_matrix(ml)
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,title='Angular distance tree, mixed clones')
ggsave('../chronoloG/figures/wiki/fig6.pdf',width=7,height=5)










## generate random genotype for the average normal (i.e. zygote)
normal <- get_normal_genotype(n_markers)

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
#pc$clonality <- 1 ## for now, make all samples 80% clonal

## get mixing fractions matrix
mix <- get_mixing_proportions(pc,uniform_mixing=T)
plot_mixtures(mix,'Purity [0-100%], clonality [50-100%], random contribution from other clones')
ggsave('../chronoloG/figures/wiki/fig3.png',width=7,height=5)

## get the admixed marker lengths
mu <- get_mean_marker_lengths(gt, mix, n_markers, anonymize=T)
ang::write_distance_matrix( round(head(t(mu), 10),4), 'figures/wiki/table1.txt')

ad <- get_angular_distance_matrix(mu)
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,title='Angular distance tree')
ggsave('../chronoloG/figures/wiki/fig4.png',width=7,height=5)












## QC
gt1 <- gt[,1:n_markers]
gt2 <- gt[,(n_markers+1):(n_markers*2)]
gt <- (gt1 + gt2) / 2

## plot the genotypes
plot_genotypes <- function(gt) {
    x <- cbind(clone=rownames(gt), as.data.table(gt))
    x <- melt(x, id.vars='clone')
    f <- function(s) {
        ss <- strsplit(s,'[.]')[[1]]
        list(marker=ss[1], copy=ss[2])
    }
    x <- cbind(x, rbindlist(lapply(as.character(x$variable), f)))

    p <- ggplot(x, aes(x=marker,y=value)) +
        geom_point(aes(color=clone),size=0.5,pch=19,position=position_jitter(width=0.1,height=0)) 

}


plot_genotypes(gt)
ggsave('../chronoloG/figures/example_genotypes.pdf',width=9,height=7)



## Example calculation of mean marker length
props <- mix['P4',4:ncol(mix)]; props
lengths <- gt[names(props),grep('^m3[.]',colnames(gt))]; lengths
sum(props * rowSums(lengths,na.rm=T)) / sum(props * rowSums(!is.na(lengths)))
mu['P4','m3']







## QC
gt1 <- gt[,1:n_markers]
gt2 <- gt[,(n_markers+1):(n_markers*2)]
qc <- (gt1 + gt2) / 2

get_anonymized_marker_lengths <- function(gt) {
    marker_min_meanlengths <- apply(gt, 2, min)
    for (mi in 1:n_markers) {
        gt[, mi] <- gt[, mi] - marker_min_meanlengths[mi]
    }
    gt
}

mu2 <- get_anonymized_marker_lengths(qc)
ad <- get_angular_distance_matrix(mu2)
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,title='Angular distance tree (10,000 gens from zygote to cancer ancestor)')









## QCing get genotypes
clones <- tree$tip.label[tree$tip.label!='normal']
n_clones <- length(clones)


get_normal_and_fcc_genotypes <- function(n_markers, mu, gens_until_first_cancer_cell) {
    marker_names <- c(paste0('m',1:n_markers,'.1'), paste0('m',1:n_markers,'.2'))

    ## normal is the average germline
    normal <- t(as.matrix(sample(10:25,replace=T,size=n_markers*2)))
    colnames(normal) <- marker_names

    ## ancestor is based on the germline after large number of divisions, and an early WGD event
    first_cancer_cell <- rcpp_mutate_length_matrix(copy(normal), mu, gens_until_first_cancer_cell)
    colnames(first_cancer_cell) <- marker_names
    list(normal=normal, first_cancer_cell=first_cancer_cell)
}




## get copy number alterations and poly-G indels in single simulation
indels_and_cnas <- get_indels_for_tree(tree, max_gens, mu=1e-4, n_markers, max_ploidy=6) 

## copy number deletions and duplications
cnas <- indels_and_cnas[,(n_markers*5):(n_markers*6)]
cnas[cnas < -1] <- -1
cnas[cnas > 1] <- 1

## copy number deletions and duplications
indels <- indels_and_cnas[,1:(n_markers*4)]










.get_genotypes <- function(indels, normal, first_cancer_cell) {

    clones <- rownames(indels)[rownames(indels)!='normal']
    n_clones <- length(clones)
    nc <- ncol(indels)

    ## gt is the genotype matrix, each tumor starts with 1st cancer cell
    gt <- matrix(nrow=n_clones,ncol=nc)
    for(i in 1:nrow(gt)) gt[i,] <- as.integer(first_cancer_cell)
    rownames(gt) <- clones
    gt <- rbind(gt, normal)
    rownames(gt)[nrow(gt)] <- 'normal'

    ## add the indels encountered to each clone's original genotype
    gt <- gt + indels
    gt
}

gt <- .get_genotypes(indels, normal, first_cancer_cell)







.get_indels_for_tree <- function(tree,max_gens,mu,n_markers) { 

    ## extract tree data
    d <- as.data.table(ggtree(tree)$data)
    d <- d[order(x,decreasing=F),]
    d <- d[-which(d$label=='normal'),]
    xpos <- d$x
    xpos <- xpos / max(xpos)
    d$gen_time <- round(xpos * max_gens)
    d$gens_to_run <- c(0,diff(d$gen_time))
    d <- d[gens_to_run > 0,]


    ## for each branch, run simulations for the according number of generations.
    simulate_indels_over_gens_for_branch <- function(gens,mu,n_markers) {
        gt_0 <- t(as.matrix(rep(0,n_markers*2))) ## initial genotype for each branch
        gt_T <- rcpp_mutate_length_matrix(gt_0, mu, gens)
        out <- as.list(gt_T[1,])
    }
    indels_for_each_branch <- lapply(d$gens_to_run, simulate_indels_over_gens_for_branch, mu, n_markers)
    indels_for_each_branch <- rbindlist(indels_for_each_branch)
    indels_for_each_branch$parent <- d$parent
    indels_for_each_branch$node <- d$node
    tips <- d$label[d$isTip==T]


    ## get the path from the normal to the tip. This includes an extra node for the added normal which we should drop.
    get_node_path_to_normal <- function(to, from, tree, d) {
        node_path_to_normal <- ape::nodepath(tree, from=which(tree$tip.label==from), to=which(tree$tip.label==to))
        node_path_to_normal[node_path_to_normal %in% c(d$parent,d$node[d$label==to])]
    }
    path_to_normal <- lapply(tips, get_node_path_to_normal, 'normal', tree, d)
    names(path_to_normal) <- tips


    ## extract the indels encountered along each branch on the path from the normal to the tip.
    get_indels_for_path <- function(v, indels_for_each_branch) {
        len <- length(v)-1
        f <- function(i, v, indels_for_each_branch) {
            from <- v[i]; to <- v[i+1]
            indels <- indels_for_each_branch[parent==from & node==to,] 
            indels
        }
        indels_along_path <- rbindlist(lapply(1:len, f, v, indels_for_each_branch))
        indels_along_path    
    }
    indels_along_path_to_normal <- lapply(path_to_normal, get_indels_for_path, indels_for_each_branch)
    ## for QC:
    #indels_along_path_to_normal$M2[,100:118]
    #indels_along_path_to_normal$M1[,100:118]
    #indels_along_path_to_normal$P1
    #indels_along_path_to_normal$P2    


    ## for each tip, now combine the indels encountered by summing them over each branch.
    combine_indels_for_each_tip <- function(indels) {
        indels <- indels[,grep('^V[0-9]*',names(indels)),with=F]
        indels <- colSums(indels)     
        as.list(indels)
    }
    indels <- lapply(indels_along_path_to_normal, combine_indels_for_each_tip)
    out <- rbindlist(indels)
    out <- as.matrix(out)
    rownames(out) <- names(indels)

    
    ## add the normal without any indels back in
    toadd <- matrix(nrow=1,ncol=ncol(out))
    toadd[is.na(toadd)] <- 0
    out <- rbind(out, toadd)
    rownames(out)[nrow(out)] <- 'normal'
    out
}











## get allele-specific lengths for each poly-G marker in each clone given indels and CNAs
gt <- get_genotypes_with_cnas(tree, max_gens,  
                              n_markers=58, 
                              mu.indel=1e-3, 
                              mu.cna=0,
                              gens_until_first_cancer_cell=1e3)

## QCing get_genotypes_with_cnas()
gt1 <- gt[,1:n_markers]
gt2 <- gt[,(n_markers+1):(n_markers*2)]
qc <- (gt1 + gt2) / 2

get_anonymized_marker_lengths <- function(gt) {
    marker_min_meanlengths <- apply(gt, 2, min)
    for (mi in 1:n_markers) {
        gt[, mi] <- gt[, mi] - marker_min_meanlengths[mi]
    }
    gt
}

mu2 <- get_anonymized_marker_lengths(qc)
ad <- get_angular_distance_matrix(mu2)
ad_tree <- nj(ad)
plot_simulated_tree(ad_tree,'Angular distance tree')














pindel <- function(k,t,mu) {
    if(length(k) > 1) {
        sapply(k, pindel, t, mu)
    } else {
        a <- exp(-0.5*(k^2)/(mu*t)) / (sqrt(mu*t)*sqrt(2*pi))
        min(c(1,a))
    }
}

x <- seq(-10,10,by=1e-3)
ps <- pindel(x, 10000, 1e-4)
ps <- ps / sum(ps)
tmp <- data.table(x=x,p=ps)
p <- ggplot(tmp, aes(x=x,y=p)) +
    geom_line()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# brownian motion one-step model
# pg 245, eq 15.28 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## probability of k steps displacement in t generations given mu mutation rate
p <- function(t,k,mu) {
    a <- exp(-0.5*(k^2)/(mu*t)) / (sqrt(mu*t)*sqrt(2*pi))
    min(c(1,a))
}

## get the MLE for the number of gens between the samples
get_mle <- function(value, mu) { 
    tmax <- 100*(1/mu)
    gens <- seq(0,tmax,by=1e-5*tmax)
    likelihoods <- sapply(gens, p, value, mu, USE.NAMES=F)
    tmp <- data.table(gen=gens, likelihood=likelihoods)
    tmp <- tmp[!is.na(likelihood),]
    tmp
    #tmp <- tmp[order(likelihood,decreasing=T),]
    #as.list(tmp[1,])
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# deconstruct the angular distance, can I tell how it relates to P?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gt_true <- get_observed_data(gt)

Y4 <- gt_true['P4',] - gt_true['normal',]
Z4 <- Y4 / sqrt(sum(Y4^2))
Y3 <- gt_true['P3',] - gt_true['normal',]
Z3 <- Y3 / sqrt(sum(Y3^2))
acos(pmin(pmax(sum(Z3*Z4), -1), 1)) # 0.2394369

ad_true <- get_angular_distance_matrix(gt_true)
ad_true['P4','P3'] # 0.2394369 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dm <- copy(gt_true)
for(i in 1:nrow(dm)) dm[i,] <- dm[i,] - dm['normal',]
x4 <- dm['P4',]
x3 <- dm['P3',]
l_P4 <- lapply(x4, get_mle, 1e-4)
for(i in 1:58) l_P4[[i]]$marker <- i
l_P3 <- lapply(x3, get_mle, 1e-4)
for(i in 1:58) l_P3[[i]]$marker <- i
P4 <- rbindlist(l_P4)
P4$clone <- 'P4'
P3 <- rbindlist(l_P3)
P3$clone <- 'P3'

## likelihood of the observed value at each marker, given some generation time
P3_ <- dcast(clone + gen ~ marker, value.var='likelihood', data=P3)
P4_ <- dcast(clone + gen ~ marker, value.var='likelihood', data=P4)


gt_obs <- get_observed_data(gt,purities)
dm <- copy(gt_obs)

## subtract normal and normalize each sample
for(i in 1:nrow(dm)) {
    ## this gives the unit vector for each sample's markers
    dm[i,] <- dm[i,] - dm['normal',]
    dm[i,] <- dm[i,] / sqrt(sum(dm[i,]^2))
}


x4 <- dm['P4',]
x3 <- dm['P3',]




P34 <- rbind(P4,P3)

dcast(gen ~ marker, )



l <- rbindlist(l)
dm <- as.data.table(cbind(dm, l))
dm <- dcast(Var1 ~ Var2, value.var='gen', data=dm)
dm <- d2m(dm)
dm <- dm - 10
t <- nj(dm)
plot(t)


tmp <- copy(gt)
for(i in 1:nrow(tmp)) tmp[i,] <- tmp[i,] - gt['normal',]
cor(tmp[1,],tmp[2,],method='pearson')

dm <- as.matrix(dist(gt, method='euclidean'))
x <- dm['normal',]
l2 <- lapply(x, get_mle, 1e-3)
l2 <- rbindlist(l2)
l2$sample <- names(x)
l2 <- l2[order(gen),]






l$sample <- names(x)
l <- l[order(gen),]

dm <- as.matrix(dist(gt, method='euclidean'))
x <- dm['normal',]
l2 <- lapply(x, get_mle, 1e-3)
l2 <- rbindlist(l2)
l2$sample <- names(x)
l2 <- l2[order(gen),]

dm <- as.matrix(dist(gt, method='euclidean'))
x <- dm['normal',]
l3 <- lapply(x, get_mle, 1e-5)
l3 <- rbindlist(l3)
l3$sample <- names(x)
l3 <- l3[order(gen),]









<- dm['normal',]
gens <- seq(0,3000,by=1)
likelihoods <- sapply(gens, get_likelihood, x[1], 1e-4, USE.NAMES=F)
tmp <- data.table(gen=gens, likelihood=likelihoods)
tmp <- tmp[!is.na(likelihood),]
tmp <- tmp[order(likelihood,decreasing=T),]
tmp[1,]





x <- gt[1,] - gt['normal',]
gens <- seq(0,3e3,by=10)
likelihoods <- sapply(gens, get_likelihood, x, 1e-4, USE.NAMES=F)
tmp <- data.table(gen=gens, likelihood=likelihoods)
tmp <- tmp[!is.na(likelihood),]
tmp <- tmp[order(likelihood,decreasing=T),]


## get the most likely distance in gens between two samples
x <- gt[1,] - gt[2,]
gens <- seq(0,3e3,by=1)
likelihoods <- sapply(gens, get_likelihood, x, 1e-4, USE.NAMES=F)
tmp <- data.table(gen=gens, likelihood=likelihoods)
tmp <- tmp[!is.na(likelihood),]
tmp <- tmp[order(likelihood,decreasing=T),]

x <- gt[1,] - gt[3,]
gens <- seq(0,3e3,by=1)
likelihoods <- sapply(gens, get_likelihood, x, 1e-4, USE.NAMES=F)
tmp <- data.table(gen=gens, likelihood=likelihoods)
tmp <- tmp[!is.na(likelihood),]
tmp <- tmp[order(likelihood,decreasing=T),]




pt <- ggplot(tmp, aes(x=gen,y=likelihood)) +
    geom_line()






## euclidean distance NJ tree (pure samples)
gt_true <- get_observed_data(gt)
dm <- dist(gt_true, method='euclidean')
t <- nj(dm)





estimate.mu(t, node.dates, p.tol = 0.05)
estimate.dates(t, node.dates, mu = estimate.mu(t, node.dates),
               min.date = -.Machine$double.xmax, show.steps = 0,
               opt.tol = 1e-8, nsteps = 1000,
               lik.tol = 0, is.binary = is.binary.phylo(t))





d <- sqrt()


library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")
dm  <- dist.ml(primates)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

fit = pml(treeNJ, data=primates)







## angular distance tree with the true mean-marker lengths
gt_true <- get_observed_data(gt)
ad_true <- get_angular_distance_matrix(gt_true)
tree_true <- nj(ad_true)
p_real <- plot_simulated_tree(tree_true,'Angular distance NJ tree from simulated genotypes\n(pure tumors, no noise)', purities)
ggsave('../chronoloG/figures/example_angular_distance_tree_NJ.pdf',width=7,height=7)


## angular distance tree with the admixed tumor+normal mean-marker lengths
gt_admixed_no_noise <- get_observed_data(gt, purities, sd.technical.error=0)
ad_admixed_no_noise <- get_angular_distance_matrix(gt_admixed_no_noise)
true_admixed_no_noise <- nj(ad_admixed_no_noise)
p_admixed_no_noise <- plot_simulated_tree(true_admixed_no_noise,'Admixed no noise',purities)
ggsave('../chronoloG/figures/example_angular_distance_tree_NJ_admixed.pdf',width=7,height=7)


## angular distance tree with the admixed tumor+normal mean-marker lengths
gt_admixed_with_noise <- get_observed_data(gt, purities, sd.technical.error=0.02)
ad_admixed_with_noise <- get_angular_distance_matrix(gt_admixed_with_noise)
true_admixed_with_noise <- nj(ad_admixed_with_noise)
p_admixed_with_noise <- plot_simulated_tree(true_admixed_with_noise,'Admixed with technical error (SD=0.02)',purities)
ggsave('../chronoloG/figures/example_angular_distance_tree_NJ_admixed_with_noise_sd=0.02.pdf',width=7,height=7)


gt_admixed_with_noise <- get_observed_data(gt, purities, sd.technical.error=0.1)
ad_admixed_with_noise <- get_angular_distance_matrix(gt_admixed_with_noise)
true_admixed_with_noise <- nj(ad_admixed_with_noise)
p_admixed_with_noise <- plot_simulated_tree(true_admixed_with_noise,'Admixed with technical error (SD=0.1)',purities)
ggsave('../chronoloG/figures/example_angular_distance_tree_NJ_admixed_with_noise_sd=0.1.pdf',width=7,height=7)








#require(adephylo)
#require(phangorn)
#require(ppcor)
#require(castor)
#require(phytools)
#require(TreeTools)
#source('/Users/alexgorelick/lab_repos/chronoloG/R/func.R')
#sourceCpp('src/simulate.cpp')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3/27/2022
# - To add: 
#   2. create 'samples' which are admixtures of multiple clones and normal




p1 <- ggtree()


plot(tree, type='u')






## generate an amount of time and size for each clone to have grown from the first cell.
dm <- as.matrix(distTips(tree,method='patristic'))
dm <- dm['normal',]
dm <- round(max_gens * dm / max(dm))
run_for_gens <- function(dm) {
    res <- random_generations(max_gens=dm[1], max_cells=1e12)
    list(cells=res$final_cells,gentime=res$generations)    
}
l <- lapply(dm, run_for_gens)
d <- rbindlist(l)
d$clone <- names(dm)
d <- d[clone!='normal',]
d <- d[order(gentime),]
d$add_gentime <- c(0,diff(d$gentime))
#l <- l[2:nrow(l),]







## we will then generate true genotypes for each clone/met. We will simulate the observed samples
## admixed with normals, and also admixed with samples based on how close they are on the tree.



## get lineage-map
map <- get_lineage_map(tree) # rows are the tips. columns are parent nodes







l$detectable <- T
l[cells < 1e6 & clone!='normal',detectable:=F]







## randomly pick one sample on the met tree to be replace be a PT sample

#replace_with_PT <- sample(mtree$tip.label,1)
#mtree$tip.label[mtree$tip.label==replace_with_PT] <- names(met_seeding)
z <- bind.tree(ptree, mtree, where=met_seeding)
plot(z)



## mets
mtree <- sim_chronology(5,n_samples_is_average=F)
mtree$tip.label <- gsub('s','M',mtree$tip.label)
mtree$tip.label <- gsub('normal',paste0('P',2),mtree$tip.label)

#mtree <- rtree(5)
mtree <- ape::root(mtree, outgroup=which(mtree$tip.label=='t1'))
newtree <- ape::bind.tree(ptree, mtree, position=2)


pt_clones <- tree$tip.label[grep('^P',tree$tip.label)]
met_founder_prob <- 0.2
n_met_founders <- rbinom(size=length(pt_clones), n=1, met_founder_prob)



samples <- tree$tip.label[tree$tip.label!='normal']
n_samples <- length(samples) - 1
n_mets <- rbinom(size=n_samples,prob=0.75,n=1)
if(n_mets < 2) n_mets <- 2
if(n_mets >= n_samples - 2) n_mets <- n_samples - 2
mets <- sample(samples,n_mets)
tree$tip.label[tree$tip.label %in% mets] <- gsub('s','M',tree$tip.label[tree$tip.label %in% mets])
tree$tip.label[!tree$tip.label %in% mets] <- gsub('s','P',tree$tip.label[!tree$tip.label %in% mets])




seeds <- seq(100,1e6)
seeds <- sample(seeds, size=1e3, replace=F)
tree <- sim_chronology(5,n_samples_is_average=F)
samples <- tree$tip.label
samples <- samples[samples!='normal']
sample(samples)







s <- 0.004
k <- 5
death_prob <- 0.5*(1-s)^k
birth_prob <- 1 - death_prob
ratio <- birth_prob / death_prob


plot_pearson_against_gens <- function(result) {
    ## get generations separating pairs of samples
    tree <- result$tree
    dm <- as.matrix(adephylo::distTips(tree,method='patristic'))
    n_gens <- result$param$n_gens
    dm <- dm / unique(dm['normal',colnames(dm)!='normal']) ## normalize by the distance to the normal
    dm <- dm * n_gens

    d <- result$gt_obs_centered
    Z <- d['normal',]
    for(r in 1:nrow(d)) d[r,] <- d[r,] - Z

    mat <- cor(t(d),method='pearson')
    mat <- mat[-nrow(mat),-ncol(mat)]
    mat[lower.tri(mat)] <- NA
    for(r in 1:nrow(mat)) mat[r,r] <- NA
    mat2 <- reshape2::melt(mat)

    dm[lower.tri(dm)] <- NA
    for(r in 1:nrow(dm)) dm[r,r] <- NA
    dm2 <- reshape2::melt(dm)

    mat2 <- as.data.table(mat2)
    dm2 <- as.data.table(dm2)
    names(dm2)[3] <- 'gens'
    names(mat2)[3] <- 'r'

    dat <- merge(dm2, mat2, by=c('Var1','Var2'), all=F)
    dat <- dat[!is.na(r) & !is.na(gens),]
    tst <- cor.test(dat$r, dat$gens, method='pearson')
    label <- paste0('Pearson r = ',prettyNum(tst$estimate,digits=2),', P = ',prettyNum(tst$p.value,digits=1))

    m <- -1/(2*n_gens)
    miny <- ifelse(min(dat$r) < 0, min(dat$r), 0)
    p <- ggplot(dat, aes(x=gens,y=r)) +
        geom_abline(slope=m, intercept=1,color='red',linetype='dashed',size=0.25) +
        geom_text(data=dat[1,],x=0.2*n_gens*2,y=0.1,label=label,color='blue',size=4) +
        scale_y_continuous(limits=c(miny,1)) + 
        scale_x_continuous(limits=c(0,n_gens*2)) +
        polyG::theme_ang(base_size=10) +
        geom_point(pch=21,color='black',fill='lightblue',stroke=0.25,size=1.75) +
        labs(x='Number of generations apart',y='Pearson correlation (mean-length minus normal)')
    if(miny < 0) p <- p + geom_hline(yintercept=0,color='black',linetype='dashed',size=0.25)

    p
}


infer_coalescent_tree <- function(gt_obs_centered,normal_name='normal',max_gens=100) {
    #x <- readRDS('~/lab_repos/lung_mets/simulation/simulated_data/rds/k=04_n=35_seed=61512_data.rds'); gt_obs_centered <- x$gt_obs_centered
    #max_gens <- 100
    #normal_name <- 'normal'

    ## original method, predict the generations apart based on extrapolative negative linear correlation.
    d <- copy(gt_obs_centered)
    normal <- d[normal_name,]
    for(i in 1:nrow(d)) d[i,] <- d[i,] - normal
    d <- d[-nrow(d),]
    mat <- cor(t(d),method='pearson')
    mat[lower.tri(mat)] <- NA
    for(i in 1:nrow(mat)) mat[i,i] <- NA
    mat <- reshape2::melt(mat)
    mat <- as.data.table(mat)
    mat <- mat[!is.na(value),]
    mat <- mat[order(value,decreasing=T),]

    intercept <- 1
    pseudo_ngens <- max_gens * 2 ## twice the total number of gens at Y=0
    miny <- -1 
    slope <- (miny - 1) / pseudo_ngens
    mat$predicted <- (mat$value - intercept) / slope
    mat1 <- copy(mat)
    mat2 <- copy(mat)
    names(mat2) <- c('Var2','Var1','value','predicted')
    mat <- rbind(mat1, mat2)
    dm <- dcast(Var1 ~ Var2, value.var='predicted', data=mat)
    dm <- polyG::d2m(dm)
    dm[is.na(dm)] <- 0

    ## make ultrametric tree from the inferred number of generations between samples
    tree <- phangorn::upgma(dm)
    ## add the normal branch to the tree    
    edge.length <- max_gens - 0.5*max(mat$predicted) 
    where <- length(tree$tip.label) + 1 ## where place normal branch
    if(is.null(where)) where <- length(tree$tip)+1
    tip <- list(edge=matrix(c(2,1),1,2), tip.label=normal_name,  edge.length=edge.length, Nnode=1)
    class(tip) <- "phylo"
    tree <- bind.tree(tree,tip,where=where)
    tree_corr_minus_normal_rooted <- phytools::reroot(tree, node.number=grep(normal_name,tree$tip.label))



    ## UPGMA with partial correlation after removing the normal, then added normal
    d <- copy(gt_obs_centered)
    normal <- d['normal',]
    for(r in 1:nrow(d)) d[r,] <- d[r,] - normal
    mat <- t(d)
    mat <- mat[,!colnames(mat) %in% 'normal']
    pc <- pcor(mat, method = "pearson")
    pcmat <- pc$estimate
    rownames(pcmat) <- colnames(mat)
    colnames(pcmat) <- colnames(mat)
    pcmat[lower.tri(pcmat)] <- NA
    for(i in 1:nrow(pcmat)) pcmat[i,i] <- NA
    pcmat <- reshape2::melt(pcmat)
    pcmat <- as.data.table(pcmat)
    pcmat <- pcmat[!is.na(value),]
    pcmat <- pcmat[order(value,decreasing=T),]


    ## NB: because pearson partial correlation scales between -1,1, we use this as the range for the number of gens apart
    intercept <- 1
    pseudo_ngens <- max_gens * 2 ## twice the total number of gens at Y=0
    miny <- -1
    #miny <- ifelse(min(pcmat$value < 0), min(pc$mat), 0)
    slope <- (miny - 1) / pseudo_ngens
    pcmat$predicted <- (pcmat$value - intercept) / slope
    mat1 <- copy(pcmat)
    mat2 <- copy(pcmat)
    names(mat2) <- c('Var2','Var1','value','predicted')
    mat <- rbind(mat1, mat2)
    dm <- dcast(Var1 ~ Var2, value.var='predicted', data=mat)
    dm <- polyG::d2m(dm)
    dm[is.na(dm)] <- 0

    ## make ultrametric tree from the inferred number of generations between samples
    tree <- phangorn::upgma(dm)
    ## add the normal branch to the tree    
    edge.length <- max_gens - 0.5*max(mat$predicted) 
    where <- length(tree$tip.label) + 1 ## where place normal branch
    if(is.null(where)) where <- length(tree$tip)+1
    tip <- list(edge=matrix(c(2,1),1,2), tip.label=normal_name,  edge.length=edge.length, Nnode=1)
    class(tip) <- "phylo"
    tree <- bind.tree(tree,tip,where=where)
    tree_pcorr_minus_normal_rooted <- phytools::reroot(tree, node.number=grep(normal_name,tree$tip.label))



    plot_tree <- function(tree,title) { 
        ## beautify plot
        p <- ggtree(tree, layout='rect') + theme_tree2() + geom_tiplab()
        p$data$x <- 100 * p$data$x / max(p$data$x)
        p$data$node_lab <- as.character(NA)
        p$data$node_lab[p$data$isTip==F & p$data$x < 98] <- round(p$data$x[p$data$isTip==F & p$data$x < 98])
        p <- p + geom_text(aes(label=node_lab),angle=0,size=2,color='blue',hjust=-0.1)
        p <- p + ggplot2::labs(x='% of generations from 1st cancer cell')
        p <- p + ggplot2::ggtitle(title)
        p
    }


   
    corr_minus_normal_rooted = list(plot=plot_tree(tree_corr_minus_normal_rooted,
                                                   'Mean-length minus normal correlation, rooted'),
                                    tree=tree_corr_minus_normal_rooted,
                                    node_distances_from_normal=get_sample_distance_from_normal(tree_corr_minus_normal_rooted),
                                    unique_pairwise_node_distances=get_unique_pairwise_distances(tree_corr_minus_normal_rooted))

    pcorr_minus_normal_rooted = list(plot=plot_tree(tree_pcorr_minus_normal_rooted,
                                                    'Mean-length minus normal partial-correlation, rooted'),
                                    tree=tree_pcorr_minus_normal_rooted,
                                    node_distances_from_normal=get_sample_distance_from_normal(tree_pcorr_minus_normal_rooted),
                                    unique_pairwise_node_distances=get_unique_pairwise_distances(tree_pcorr_minus_normal_rooted))


    list(corr_minus_normal_rooted=corr_minus_normal_rooted, pcorr_minus_normal_rooted=pcorr_minus_normal_rooted)

}


get_sample_distance_from_normal <- function(tree,normal_name='normal') {
    samples <- tree$tip.label; samples <- samples[samples!=normal_name]
    #mrca_node <- getMRCA(tree, samples)
    edge_distance_from_normal <- castor::get_pairwise_distances(tree, samples, rep(normal_name, length(samples)), as_edge_counts=T)
    out <- data.table(sample=samples, edge_distance_from_normal=edge_distance_from_normal)
    out <- out[order(edge_distance_from_normal,decreasing=F),]
    out
}


get_unique_pairwise_distances <- function(tree,normal_name='normal') {
    ## 3-columns output: sample1, sample2, node_distance. sample1 and sample2 are sorted in numeric order.
    dm <- as.matrix(distTips(tree,method='nNodes'))
    fields <- c(paste0('s',1:(nrow(dm)-1)),normal_name)
    fields <- fields[fields %in% rownames(dm)]
    dm <- dm[fields, fields]
    dm[lower.tri(dm)] <- NA
    for(i in 1:nrow(dm)) dm[i,i] <- NA
    dm <- as.data.table(reshape2::melt(dm))
    dm <- dm[!is.na(value),]
    names(dm) <- c('sample1','sample2','node_distance')
    dm
}


compare_node_distances_from_normal <- function(result,partial) {
    obs <- result$true_node_distances_from_normal

    if(partial==F) {
        pred <- result$inferred$corr_minus_normal_rooted$node_distances_from_normal
        title <- 'Distance from normal based on correlation matrix'
    } else {
        pred <- result$inferred$pcorr_minus_normal_rooted$node_distances_from_normal
        title <- 'Distance from N from partial-correlation matrix'
    }

    dat <- merge(obs, pred, by='sample', all=T)
    names(dat) <- c('sample','truth','predicted')
    dat <- dat[!is.na(truth) & !is.na(predicted),]
    r.p=cor(dat$truth, dat$predicted, method='pearson')
    r.s=cor(dat$truth, dat$predicted, method='spearman')
    max_distance <- max(c(dat$truth,dat$predicted))
    p <- ggplot(dat,aes(x=truth,y=predicted)) +
        scale_x_continuous(limits=c(0,max_distance*1.1),breaks=0:max_distance) +
        scale_y_continuous(limits=c(0,max_distance*1.1),breaks=0:max_distance) +
        geom_abline(intercept=0,slope=1,color='#bfbfbf',linetype='dashed',size=0.5) +
        geom_point(pch=21,stroke=0.25,fill='steelblue',size=2,position=position_jitter(width=0.1,height=0.1)) +
        polyG::theme_ang(base_size=12) +
        theme(panel.grid.major=element_line(color='grey92',size=0.25)) +
        labs(x='True node distance from normal',y='Predicted node distance from normal',title=title,
             subtitle=paste0('Pearson r=',prettyNum(r.p,digits=2),', Spearman r=',prettyNum(r.s,digits=2))) 
    p
}


for(seed in seeds) {
    ## use the seed to make the label
    message(seed)

    ## randomly pick the input parameters based on the seed
    set.seed(seed) 
    n_samples <- sample(4:60,size=1)     
    k <- sample(seq(1:5),1)     
    k_label <- ifelse(k < 10, paste0('0',k), as.character(k))
    n_label <- ifelse(n_samples < 10, paste0('0',n_samples), as.character(n_samples))
    label <- paste0('k=',k_label,'_n=',n_label,'_seed=',seed)

    
    ## define output files
    outdir <- '~/lab_repos/lung_mets/simulation/simulated_data'
    out_rds <- file.path(outdir,paste0('rds/',label,'_data.rds'))

    out_chrontree <- file.path(outdir,paste0('pdf/',label,'_chron_tree.pdf'))
    out_chrontree_inferred <- file.path(outdir,paste0('pdf/',label,'_chrontree_inferred.pdf'))
    out_angular_distance_tree <- file.path(outdir,paste0('pdf/',label,'_angular_distance_tree.pdf'))
    out_pearsongens <- file.path(outdir,paste0('pdf/',label,'_pearsongens.pdf'))
    out_distance_from_normal <- file.path(outdir,paste0('pdf/',label,'_distance_from_normal.pdf'))

    ## run the simulation    
    result <- chronoloG(n_samples,n_samples_is_average=F,k=k)
    result$seed <- seed
    chron_tree <- result$tree
    total_gens <- result$param$n_gens
    result$true_node_distances_from_normal <- get_sample_distance_from_normal(chron_tree)
    result$true_unique_pairwise_node_distances <- get_unique_pairwise_distances(chron_tree)
    purity <- result$purities
    purity <- data.frame(label=rownames(purity), purity)
    rownames(purity) <- NULL


    ## infer the chronology tree with various methods, save the tree plots and add the data to the result list
    inferred <- tryCatch({
        infer_coalescent_tree(result$gt_obs_centered)
    },error=function(e) {
        NULL
    })
    if(is.null(inferred)) next


    ## beautify plot
    p <- ggtree(chron_tree, layout='rect') + theme_tree2() + geom_tiplab()
    p$data$x <- total_gens * p$data$x / max(p$data$x)
    p$data$node_lab <- as.character(NA)
    p$data$node_lab[p$data$isTip==F] <- round(p$data$x[p$data$isTip==F])
    p <- p %<+% purity
    p <- p + geom_tippoint(aes(fill=purity),pch=21,stroke=0.25,size=2) 
    p <- p + geom_text(aes(label=node_lab),angle=0,size=2,color='blue',hjust=-0.1)
    p <- p + ggplot2::labs(x='Number of generations from 1st cancer cell')
    p <- p + ggplot2::ggtitle(paste0('Simulated chrono-evolution (gens=',total_gens,', k=',k,', n_samples=',n_samples,', seed=',seed,')'))
    p <- p + ggplot2::scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0.5) 
    p <- p + theme(legend.position='bottom')
    ggsave(out_chrontree,width=8,height=7)
    result$tree_plot <- p    
    ggsave(plot=inferred$corr_minus_normal_rooted$plot, out_chrontree_inferred,width=8,height=7)
    result$inferred <- inferred     


    ## also get angular distance matrix because why not
    tmp <- result$gt_obs_centered
    ad <- get_angular_distance_matrix(tmp)
    excluded_samples <- c()
    bad_samples <- names(which(colMeans(is.na(ad)) > 0.5))
    while(length(bad_samples) > 0) {
        excluded_samples <- c(excluded_samples, bad_samples)
        message(bad_samples,'\n')
        tmp <- tmp[!rownames(tmp) %in% bad_samples,]
        ad <- get_angular_distance_matrix(tmp)
        bad_samples <- names(which(rowSums(is.na(ad)) > 1))
    }
    adtree <- nj(ad)
    adtree <- phytools::reroot(adtree, node.number=grep('normal',adtree$tip.label))
    p_ad <- ggtree(adtree,layout='ape') + geom_tiplab(angle=0)
    if(length(excluded_samples) > 0) {
        title <- paste('Angular distance NJ tree, excluding:',paste(excluded_samples, collapse=', '))
    } else {
        title <- 'Angular distance NJ tree'
    }
    p_ad <- p_ad + ggplot2::ggtitle(title)
    ggsave(plot=p_ad, out_angular_distance_tree,width=9,height=7)
    result$angular_distance_tree <- p_ad
    result$angular_distance_matrix <- ad


    ## get UPGMA angular distance tree
    tree <- phangorn::upgma(ad)
    tree <- phytools::reroot(tree, node.number=grep('normal',tree$tip.label))
    p <- ggtree(tree, layout='rect') + theme_tree2() + geom_tiplab()
    p$data$x <- 100 * p$data$x / max(p$data$x)
    p$data$node_lab <- as.character(NA)
    p$data$node_lab[p$data$isTip==F & p$data$x < 98] <- round(p$data$x[p$data$isTip==F & p$data$x < 98])
    p <- p + geom_text(aes(label=node_lab),angle=0,size=2,color='blue',hjust=-0.1)
    p <- p + ggplot2::labs(x='% of ultrametric mutational-time')
    p <- p + ggplot2::ggtitle(title)
    out_upgma_ad_tree <- file.path(outdir,paste0('pdf/',label,'_angular_distance_upgma.pdf'))
    ggsave(plot=p, out_upgma_ad_tree, width=10, height=8)


    ## extract/save the node distance from each sample to the normal
    dm <- as.matrix(distTips(tree,method='nNodes'))
    out <- dm['normal',]
    out <- data.table(sample=names(out), distance_to_normal=as.integer(out))
    result$angular_distance_nodes_from_normal <- out


    ## plot of node distance from the normal for both 
    p_corr_distance_from_normal <- compare_node_distances_from_normal(result,partial=F)
    ggsave(plot=p_corr_distance_from_normal, out_distance_from_normal,width=7,height=6)
    result$p_corr_distance_from_normal <- p_corr_distance_from_normal


    #p_pcorr_distance_from_normal <- compare_node_distances_from_normal(result,partial=T)
    #distance_from_normal_plots <- list(p_corr_distance_from_normal, p_pcorr_distance_from_normal)
    #pmerged_distance_from_normal <- plot_grid(plotlist=distance_from_normal_plots, nrow=1) 

    ## plot pearson R against generations separating two samples
    #p2 <- plot_pearson_against_gens(result)
    #p2 <- p2 + ggplot2::ggtitle(paste0('Pearson r vs gens apart (gens=',total_gens,', k=',k,', n_samples=',n_samples,', seed=',seed,')'))
    #ggsave(out_pearsongens,width=7,height=5)
    #result$pearsongen_plot <- p2

    saveRDS(result,file=out_rds)
}



## based on the above simulations, how good on average are the inferred node distance from normal?
## how does this relate to the simulated k or n parameters?


tmp <- readRDS('~/lab_repos/lung_mets/simulation/simulated_data/rds/k=04_n=35_seed=61512_data.rds')

rtree_nodes_from_normal <- function(n_samples) {
    tree <- rtree(n_samples, tip.label=paste0('s',1:n_samples))
    branch_lengths <- tree$edge.length
    log10_branch_lengths <- log10(branch_lengths)
    mu <- mean(log10_branch_lengths)
    s <- sd(log10_branch_lengths)
    edge.length <- 10^rnorm(mean=mu,sd=s,n=1)
    where <- length(tree$tip.label) + 1 
    if(is.null(where)) where<-length(tree$tip)+1
    tip<-list(edge=matrix(c(2,1),1,2),
              tip.label='normal',
              edge.length=edge.length,
              Nnode=1)
    class(tip)<-"phylo"
    tree <- bind.tree(tree,tip,where=where)
    tree <- phytools::reroot(tree, node.number=which(tree$tip.label=='normal'))
    dm <- as.matrix(distTips(tree,method='nNodes'))
    out <- dm['normal',]
    out <- data.table(sample=names(out), node_distances_from_normal=as.integer(out))
    out
}


check_correlation <- function(file,method) {
    message(file)
    x <- readRDS(file)
    par <- x$params
    k <- par$k; n <- par$n_samples; gens <- par$n_gens
    truth <- x$true_node_distances_from_normal
    if(method=='corr') {
        pred <- x$inferred$corr_minus_normal_rooted$node_distances_from_normal
    } else if(method=='ad') {
        pred <- x$angular_distance_nodes_from_normal
    } else if(method=='random') {
        pred <- rtree_nodes_from_normal(n)
    }
    d <- merge(truth, pred, by='sample', all=F)
    names(d) <- c('sample','truth','pred')
    d <- d[!is.na(truth) & !is.na(pred),]
    r <- cor(d$truth, d$pred, method='pearson')
    list(k=k, n=n, gens=gens, r=r, seed=x$seed, method=method)
}

library(parallel)
rds_files <- dir('~/lab_repos/lung_mets/simulation/simulated_data/rds',full.names=T)

l_corr <- mclapply(rds_files, check_correlation,'corr', mc.cores=4)
info_corr <- rbindlist(l_corr)
info_corr$method <- 'corr'

l_ad <- mclapply(rds_files, check_correlation,'ad', mc.cores=4)
info_ad <- rbindlist(l_ad)
info_ad$method <- 'ad'

set.seed(42)
l_random <- mclapply(rds_files, check_correlation,'random', mc.cores=4)
info_random <- rbindlist(l_random)
info_random$method <- 'random'


info <- rbind(info_corr, info_ad, info_random)
info[n < 10, cat:='[4,10)']
info[n >= 10 & n < 20, cat:='[10,20)']
info[n >= 20 & n < 30, cat:='[20,30)']
info[n >= 30 & n < 40, cat:='[30,40)']
info[n >= 40 & n < 50, cat:='[40,50)']
info[n >= 50, cat:='[50,60]']
info$cat <- factor(info$cat, levels=c('[4,10)','[10,20)','[20,30)','[30,40)','[40,50)','[50,60]'))

info[method=='corr',method:='Pseudotime']
info[method=='ad',method:='Angular distance']
info[method=='random',method:='Random']
info$method <- factor(info$method, levels=c('Pseudotime','Angular distance','Random'))

p <- ggplot(info, aes(x=method,y=r)) +
    geom_point(position=position_jitter(width=0.1,height=0),pch=21,stroke=0.25,color='black',aes(fill=method)) +
    polyG::theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position='none') +
    geom_boxplot(fill=NA,outlier.shape=NA) +
    facet_wrap(facets=~cat,nrow=1) +
    labs(x='Method used',y='Node distance to Normal\n(Pearson r between True and Predicted)') 
ggsave('simulation/figures/correlation_distance_from_normal.pdf',width=8,height=5)


#info <- dcast(k + n + seed ~ method,data=info, value.var='r')
#p <- ggplot(info, aes(x=ad,y=corr)) +
#    geom_point(pch=21,stroke=0.25,fill='#bfbfbf',size=1.25) +
#    labs(x='Angular distance UPGMA',y='Pseudotime distance UPGMA',subtitle='')





ggsave('simulation/figures/correlation_distance_from_normal.pdf',width=8,height=5.5)



x <- readRDS('~/lab_repos/lung_mets/simulation/simulated_data/rds/k=04_n=35_seed=61512_data.rds')
truth <- x$true_node_distances_from_normal
qs <- quantile(truth$edge_distance_from_normal,c(0.333,0.666),na.rm=T)
truth[edge_distance_from_normal < qs[1],group:='earlier']
truth[edge_distance_from_normal >= qs[1] & edge_distance_from_normal < qs[2],group:='intermediate']
truth[edge_distance_from_normal >= qs[2],group:='later']

pred <- x$inferred$corr_minus_normal_rooted$node_distances_from_normal
qs <- quantile(pred$edge_distance_from_normal,c(0.333,0.666),na.rm=T)
pred[edge_distance_from_normal < qs[1],group:='earlier']
pred[edge_distance_from_normal >= qs[1] & edge_distance_from_normal < qs[2],group:='intermediate']
pred[edge_distance_from_normal >= qs[2],group:='later']

d <- merge(truth[,c('sample','group'),with=F], pred[,c('sample','group'),with=F], by='sample')
names(d) <- c('sample','truth','predicted')
xtabs(~ truth + predicted, data=d)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# try same approach with the cell line data to show 
# that we recover the experimental chrono tree
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- load_marker_files('~/lab_repos/lung_mets/original_data/cell_line_experiments/marker_files/MSI1')
normal_name <- 'MSI1M2'
d <- d[nchar(sample)==10 | sample==normal_name,]
gt_obs_centered <- polyG::d2m(d)
p <- infer_coalescent_tree(gt_obs_centered,normal_name=normal_name,max_gens=100) 
p <- p + ggplot2::xlim(c(0,120))
ggsave('~/lab_repos/lung_mets/processed_data/cell_line_experiments/MSI1_infered_chronology.pdf',width=8,height=6)

d <- load_marker_files('~/lab_repos/lung_mets/original_data/cell_line_experiments/marker_files/MSS1')
normal_name <- 'MSS1M3'
d <- d[nchar(sample)==10 | sample==normal_name,]
gt_obs_centered <- polyG::d2m(d)
p <- infer_coalescent_tree(gt_obs_centered,normal_name=normal_name,max_gens=100) 
p <- p + ggplot2::xlim(c(0,120))
ggsave('~/lab_repos/lung_mets/processed_data/cell_line_experiments/MSS1_infered_chronology.pdf',width=8,height=6)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# real data: C70, other examples with independent primarie?0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subjects <- dir('~/lab_repos/lung_mets/original_data/marker_files')
ann <- fread('~/lab_repos/lung_mets/original_data/annotation_files.txt')

for(this.subject in subjects) { 
    message(this.subject)
    marker_dir <- file.path('~/lab_repos/lung_mets/original_data/marker_files',this.subject)
    ad_file <- paste0('~/lab_repos/lung_mets/original_data/angular_distance_matrices/',this.subject,'.txt')
    d <- load_marker_files(marker_dir)
    info <- ann[Subject==this.subject,]
    cohort <- unique(info$cohort)
    info <- info[,c('Sample_ID','Real_Sample_ID'),with=F]
    d <- merge(info, d, by.x='Sample_ID', by.y='sample', all=F)
    d[,Sample_ID:=NULL]
    d <- polyG::d2m(d)

    ## use the samples in the angular distance matrices
    angular_distance_matrix <- read_distance_matrix(ad_file)
    valid_samples <- rownames(angular_distance_matrix)
    d <- d[valid_samples,]
    normal_name <- c(grep('^N[0-9]',rownames(d),value=T), grep('^Normal[0-9]',rownames(d),value=T))

    ## subset for markers without any NAs
    good_markers <- names(which(colSums(is.na(d)) == 0))
    d <- d[,c(good_markers)]

    #d_tumor <- d[!rownames(d) %in% normals,]
    #d_normals <- d[rownames(d) %in% normals,]
    #if(length(normals) > 1) d_normals <- apply(d_normals, 2, median)
    #d <- rbind(d_tumor, d_normals)
    #rownames(d)[nrow(d)] <- paste(normals,collapse=',')
    #normal_name <- rownames(d)[nrow(d)]

    ## get colors for tissue types
    groups <- group_samples(d,color=T,lun=T,liv=F,per=F)
    setnames(groups,'barcode','label')
    tmp <- groups[!duplicated(groups$group),]
    cols <- tmp$color; names(cols) <- tmp$group

    ## infer chronology and plot with colored tip labels
    p <- infer_coalescent_tree(d,normal_name=normal_name,max_gens=100) 
    p <- p + ggplot2::scale_x_continuous(limits=c(0,110),breaks=seq(0,100,by=25))
    p <- p %<+% groups + geom_tiplab(aes(color=group)) + scale_color_manual(values=cols,name='Tissue') + theme(legend.position='none') + ggtitle(paste0(this.subject,' (',cohort,') inferred chronology'))

    pdf_file <- paste0('~/lab_repos/lung_mets/processed_data/patient_chronologies/',cohort,'/',this.subject,'.pdf')
    ggsave(pdf_file,width=8,height=7)
}



