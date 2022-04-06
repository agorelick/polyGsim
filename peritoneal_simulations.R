rm(list=ls())
library(polyGsim)
set.seed(100)
#source('R/func.R')

## define following parameters
n_primary_clones <- 10
n_is_avg <- T
met_prob_per_clone <- 0.2
avg_mets_per_clone <- 5
n_markers <- 58
mu <- 1e-4
bdratio <- 1.01

set.seed(42)
seeds <- round(runif(0,1e6,n=1e2)) 
seeds <- seeds[!duplicated(seeds)]


compare_trees <- function(i, control, test) {
    require(TreeDist)
    require(Quartet)

    message(i)
    params <- tryCatch({
        params <- test[[i]]
        params <- params[names(params)!='tree']

        ## Generalized RF similarity
        params$grf_similarity <- SharedPhylogeneticInfo( control[[i]]$tree, test[[i]]$tree, normalize=T)

        ## Quartet similarity
        qinfo <- QuartetStatus(trees=test[[i]]$tree, cf=control[[i]]$tree)
        qinfo <- as.data.frame(qinfo)
        params$quartet_similarity <- qinfo$s / qinfo$Q

        params$i <- i
        params
    },error=function(e) {
        NULL
    })
    params
}


run_for_purity_and_clonality <- function(this.seed, purity=NA, clonality=NA) {
    set.seed(this.seed)
    message(this.seed)

    ad_tree <- tryCatch({
        ## generate random genotype for the average normal (i.e. zygote)
        normal <- get_normal_genotype(n_markers=n_markers)

        ## get a random number of generations from zygote to last tumor sample taken/patient death
        max_gens <- random_generations(bdratio=bdratio)$generations

        ## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
        tree <- get_clone_tree(n_primary_clones=n_primary_clones, 
                               n_is_avg=n_is_avg, 
                               met_prob_per_clone=met_prob_per_clone,
                               avg_mets_per_clone=avg_mets_per_clone)

        ## get indels for the given tree
        indels <- get_indels_for_tree(tree, max_gens, mu, n_markers) 

        ## get genotypes for each clone
        gt <- get_clone_genotypes(indels, normal)

        ## get random or specified purity and clonality
        pc <- get_purity_and_clonality(tree)
        if(!is.na(purity)) pc[rownames(pc)!='normal','purity'] <- purity
        if(!is.na(clonality)) pc[rownames(pc)!='normal','clonality'] <- clonality
        
        pvals <- pc[rownames(pc)!='normal','purity']
        p.mu <- mean(pvals,na.rm=T)
        p.sd <- sd(pvals,na.rm=T)
        cvals <- pc[rownames(pc)!='normal','clonality']
        c.mu <- mean(cvals,na.rm=T)
        c.sd <- sd(cvals,na.rm=T)

        mix <- get_mixing_proportions(pc,even_mixing=F)
        mix[is.na(normal),normal:=0]

        ## get the admixed marker lengths
        ml <- get_mean_marker_lengths(gt, mix, n_markers)

        ## anonymize marker lengths so that each marker's minimum value is 0.
        ml <- get_anonymized_marker_lengths(ml)

        ## calculate the angular distance matrix from the anonymized marker lengths
        ad <- get_angular_distance_matrix(ml)

        ## return the angular distance neighbor joining tree
        tree <- nj(ad)
        tree <- phytools::reroot(tree, node.number=which(tree$tip.label=='normal'))
        list(seed=this.seed, tree=tree, size=length(tree$tip.label), p.mu=p.mu, p.sd=p.sd, c.mu=c.mu, c.sd=c.sd)

    }, error=function(e) {
        NULL
    })
    ad_tree
}

## control: perfect purity, perfect clonality
l_control <- lapply(seeds, run_for_purity_and_clonality, purity=1, clonality=1)

## random purity, perfect clonality.
l_random_purity_c1 <- lapply(seeds, run_for_purity_and_clonality, purity=NA, clonality=1)

## random purity, perfect clonality.
l_random_clonality_p1 <- lapply(seeds, run_for_purity_and_clonality, purity=1, clonality=NA)

## random purity, perfect clonality.
l_random_clonality_random_purity <- lapply(seeds, run_for_purity_and_clonality, purity=NA, clonality=NA)



control_vs_random_purity_c1 <- rbindlist(lapply(1:length(l_control), compare_trees, l_control, l_random_purity_c1))
control_vs_random_purity_c1$test <- 'Purity 0-100%, clonality=100%'
control_vs_random_clonality_p1 <- rbindlist(lapply(1:length(l_control), compare_trees, l_control, l_random_clonality_p1))
control_vs_random_clonality_p1$test <- 'Clonality 50-100% (random clone mixing), purity=100%'
control_vs_random_clonality_random_purity <- rbindlist(lapply(1:length(l_control), compare_trees, l_control, l_random_clonality_random_purity))
control_vs_random_clonality_random_purity$test <- 'Clonality 50-100% (random clone mixing), purity=0-100%'
tests <- rbind(control_vs_random_purity_c1, control_vs_random_clonality_p1, control_vs_random_clonality_random_purity)
tests <- melt(tests, id.vars=c('seed','size','p.mu','p.sd','c.mu','c.sd','test','i'))
tests[variable=='grf_similarity', variable:='GRF']
tests[variable=='quartet_similarity', variable:='Quartet']
tests$variable <- factor(tests$variable, levels=c('GRF','Quartet'))
saveRDS(tests,file='figures/test1_n100.rds')

tests <- readRDS('figures/test1_n100.rds')
p <- ggplot(tests, aes(x=test,y=value,group=seed)) +
    scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25)) + 
    geom_boxplot(outlier.shape=NA,aes(group=test,fill=test)) +
    scale_fill_brewer(palette='Accent',name='Simulated tumor samples') +
    facet_wrap(facets=~variable, ncol=2) +
    geom_line(alpha=0.3,size=0.25) + 
    geom_point() +
    ang::theme_ang(base_size=10) +
    labs(x='Comparison to ideal data',y='GRF similarity',subtitle='Similarity to ideal AD tree (100% pure, 100% clonal samples)\nn=100 random phylogenies') +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='right')
p <- ang::extract_gglegend(p)
p2 <- plot_grid(p$plot, p$legend, ncol=1, rel_heights=c(4,1))
ggsave('figures/test1.pdf',width=7,height=7)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we use a brute-force solution to find values of alpha 
# and beta for Beta distribution to achieve desired means all 
# with a standard deviation of 0.1.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

a <- seq(0,20,by=1e-1)
b <- seq(0,20,by=1e-1)
x <- as.data.table(expand.grid(shape1=a, shape2=b))
summarize <- function(i, x) {
    s1 <- x$shape1[i]
    s2 <- x$shape2[i]
    val <- rbeta(n=1e3, shape1=s1, shape2=s2)
    val <- val / 2 + 0.5
    s <- sd(val,na.rm=T)
    mu <- mean(val,na.rm=T)
    list(s1=s1, s2=s2, mu=mu, s=s)
}
info <- rbindlist(lapply(1:nrow(x), summarize, x))
info <- info[s > 0,]

x <- info[ s >= 0.099 & s <= 0.101, ]
x[,diff_from_0.9 := abs(mu - 0.9)]
x[,diff_from_0.8 := abs(mu - 0.8)]
x[,diff_from_0.7 := abs(mu - 0.7)]
x[,diff_from_0.6 := abs(mu - 0.6)]
x[,diff_from_0.5 := abs(mu - 0.5)]

head(x[order(diff_from_0.9,decreasing=F),],1)
head(x[order(diff_from_0.8,decreasing=F),],1)
head(x[order(diff_from_0.7,decreasing=F),],1)
head(x[order(diff_from_0.6,decreasing=F),],1)
head(x[order(diff_from_0.5,decreasing=F),],1)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# next test: run simulations for a variety of clonality distributions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## demonstrate the clonality distributions we will use
get_distribution <- function(n, shape1, shape2) {
    set.seed(42)
    d1 <- data.table(clonality=rbeta(n=n, shape1=shape1, shape2=shape2))
    d1[,clonality:=clonality/2 + 0.5]
    mu <- mean(d1$clonality,na.rm=T)
    s <- sd(d1$clonality, na.rm=T)
    d1$label <- paste0('mu=',round(mu,2),',sd=',round(s,2),', clonality ~ 0.5*Beta(',shape1,',',shape2,')+0.5')
    d1
}
d1 <- get_distribution(10000, 2.3, 0.5) ## 90%
d2 <- get_distribution(10000, 3.4, 2.2) ## 80%
d3 <- get_distribution(10000, 2, 2.9)   ## 70%
d4 <- get_distribution(10000, 0.5, 2.5) ## 60%
dat <- rbind(d1,d2,d3,d4)
dat$label <- factor(dat$label, levels=unique(dat$label))

p <- ggplot(dat, aes(x=clonality)) +
    scale_x_continuous(limits=c(0,1), breaks=seq(0,1,by=0.25)) + 
    geom_histogram(bins=50,color='white',fill='#bfbfbf') + 
    facet_wrap(facets=~label, ncol=1,scale='free_y') +
    ang::theme_ang(base_size=12) +
    labs(x='Clonality',y='Frequency')
ggsave('figures/test2_clonality_distributions.pdf',width=7,height=7)


## run simulations for a variety of clonality distributions
set.seed(42)
seeds <- round(runif(0,1e6,n=1e3)) 
while(any(duplicated(seeds))) {
    message('retrying')
    seeds <- round(runif(0,1e6,n=1e3)) 
}

run_for_clonality_shapes <- function(this.seed, c.shape1, c.shape2) {
    set.seed(this.seed)
    message(this.seed)

    ad_tree <- tryCatch({
        ## generate random genotype for the average normal (i.e. zygote)
        normal <- get_normal_genotype(n_markers=n_markers)

        ## get a random number of generations from zygote to last tumor sample taken/patient death
        max_gens <- random_generations(bdratio=bdratio)$generations

        ## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
        tree <- get_clone_tree(n_primary_clones=n_primary_clones, 
                               n_is_avg=n_is_avg, 
                               met_prob_per_clone=met_prob_per_clone,
                               avg_mets_per_clone=avg_mets_per_clone)

        ## get indels for the given tree
        indels <- get_indels_for_tree(tree, max_gens, mu, n_markers) 

        ## get genotypes for each clone
        gt <- get_clone_genotypes(indels, normal)

        ## get random or specified purity and clonality
        pc <- get_purity_and_clonality(tree)
        if(!is.na(c.shape1) & !is.na(c.shape2)) {
            pc$clonality <- 0.5*rbeta(n=nrow(pc), shape1=c.shape1, shape2=c.shape2) + 0.5
            vals <- 0.5*rbeta(n=1e5, shape1=c.shape1, shape2=c.shape2) + 0.5
            mu <- mean(vals,na.rm=T)
            s <- sd(vals,na.rm=T)
        } else {
            pc$clonality <- 1
            mu <- 1
            s <- 0
        }
        pc['normal','clonality'] <- 0
        
        cvals <- pc[rownames(pc)!='normal','clonality']
        c.mu <- mean(cvals,na.rm=T)
        c.sd <- sd(cvals,na.rm=T)

        mix <- get_mixing_proportions(pc,even_mixing=F)
        mix[is.na(normal),normal:=0]

        ## get the admixed marker lengths
        ml <- get_mean_marker_lengths(gt, mix, n_markers)

        ## anonymize marker lengths so that each marker's minimum value is 0.
        ml <- get_anonymized_marker_lengths(ml)

        ## calculate the angular distance matrix from the anonymized marker lengths
        ad <- get_angular_distance_matrix(ml)

        ## return the angular distance neighbor joining tree
        tree <- nj(ad)
        tree <- phytools::reroot(tree, node.number=which(tree$tip.label=='normal'))
        list(seed=this.seed, tree=tree, size=length(tree$tip.label), c.shape1=c.shape1, c.shape2=c.shape2, c.mu=c.mu, c.sd=c.sd, mu=mu, s=s)

    }, error=function(e) {
        NULL
    })
    ad_tree
}

## control: perfect purity, perfect clonality
l_control <- lapply(seeds[1:100], run_for_clonality_shapes, c.shape1=NA, c.shape2=NA)
l_test_0.9 <- lapply(seeds[1:100], run_for_clonality_shapes, c.shape1=2.3, c.shape2=0.5)
l_test_0.8 <- lapply(seeds[1:100], run_for_clonality_shapes, c.shape1=3.4, c.shape2=2.2)
l_test_0.7 <- lapply(seeds[1:100], run_for_clonality_shapes, c.shape1=2, c.shape2=2.9)
l_test_0.6 <- lapply(seeds[1:100], run_for_clonality_shapes, c.shape1=0.5, c.shape2=2.5)

res_0.9 <- rbindlist(lapply(1:length(l_control), compare_trees, l_control, l_test_0.9))
res_0.9$clonality <- 'mu=0.9,s=0.1'
res_0.8 <- rbindlist(lapply(1:length(l_control), compare_trees, l_control, l_test_0.8))
res_0.8$clonality <- 'mu=0.8,s=0.1'
res_0.7 <- rbindlist(lapply(1:length(l_control), compare_trees, l_control, l_test_0.7))
res_0.7$clonality <- 'mu=0.7,s=0.1'
res_0.6 <- rbindlist(lapply(1:length(l_control), compare_trees, l_control, l_test_0.6))
res_0.6$clonality <- 'mu=0.6,s=0.1'
res <- rbind(res_0.9, res_0.8, res_0.7, res_0.6)
res$clonality <- factor(res$clonality, levels=unique(res$clonality))

tests <- melt(res, id.vars=c('seed','size','c.shape1','c.shape2','c.mu','c.sd','mu','s','i','clonality'))
tests[variable=='grf_similarity', variable:='GRF']
tests[variable=='quartet_similarity', variable:='Quartet']
tests$variable <- factor(tests$variable, levels=c('GRF','Quartet'))
tests$test <- tests$clonality #paste0(tests$clonality,' vs mu=1,s=0')
tests$test <- gsub('s=','SD=',tests$test)
tests$test <- gsub('mu=','Mean=',tests$test)
tests$test <- factor(tests$test, levels=unique(tests$test))
saveRDS(tests,file='figures/test2_data_n100.rds')

#tests <- readRDS('figures/test2_data_n100.rds')
p <- ggplot(tests, aes(x=test,y=value,group=seed)) +
    scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25)) + 
    geom_boxplot(outlier.shape=NA,aes(group=test,fill=test)) + 
    scale_fill_brewer(palette='Accent',name='Simulated tumor samples') +
    facet_wrap(facets=~variable, ncol=2) +
    #geom_line(alpha=0.3,size=0.25) + ## looks weird in this plot
    geom_point(position=position_jitter(width=0.15,height=0),pch=21,fill='grey',color='black',stroke=0.25) +
    ang::theme_ang(base_size=10) +
    labs(x='Beta distribution used for % clonality',y='Similarity',subtitle='Similarity to ideal AD tree with 100% pure, 100% clonal samples\nn=100 random phylogenies with purity 0-100% and clonality from indicated distribution.') +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave('figures/test2_similarity_boxplot.pdf',width=7,height=5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# examples of the worst similarities?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extract_data_for_seed <- function(seed, c.shape1, c.shape2) {
    set.seed(seed)
    normal <- get_normal_genotype(n_markers=n_markers)

    ## get a random number of generations from zygote to last tumor sample taken/patient death
    max_gens <- random_generations(bdratio=bdratio)$generations

    ## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
    tree <- get_clone_tree(n_primary_clones=n_primary_clones, 
                           n_is_avg=n_is_avg, 
                           met_prob_per_clone=met_prob_per_clone,
                           avg_mets_per_clone=avg_mets_per_clone)

    ## get indels for the given tree
    indels <- get_indels_for_tree(tree, max_gens, mu, n_markers) 

    ## get genotypes for each clone
    gt <- get_clone_genotypes(indels, normal)

    ## get random or specified purity and clonality
    pc <- get_purity_and_clonality(tree)
    if(!is.na(c.shape1) & !is.na(c.shape2)) {
        pc$clonality <- 0.5*rbeta(n=nrow(pc), shape1=c.shape1, shape2=c.shape2) + 0.5
        vals <- 0.5*rbeta(n=1e5, shape1=c.shape1, shape2=c.shape2) + 0.5
        mu <- mean(vals,na.rm=T)
        s <- sd(vals,na.rm=T)
    } else {
        pc$clonality <- 1
        mu <- 1
        s <- 0
    }
    pc['normal','clonality'] <- 0

    cvals <- pc[rownames(pc)!='normal','clonality']
    c.mu <- mean(cvals,na.rm=T)
    c.sd <- sd(cvals,na.rm=T)

    mix <- get_mixing_proportions(pc,even_mixing=F)
    mix[is.na(normal),normal:=0]

    ## get the admixed marker lengths
    ml <- get_mean_marker_lengths(gt, mix, n_markers)

    ## anonymize marker lengths so that each marker's minimum value is 0.
    ml <- get_anonymized_marker_lengths(ml)

    ## calculate the angular distance matrix from the anonymized marker lengths
    ad <- get_angular_distance_matrix(ml)

    list(tree=tree, max_gens=max_gens, gt=gt, pc=pc, mix=mix, ml=ml, ad=ad)
}


x <- tests[order(value,decreasing=F),]


index <- x$i[1]
control <- l_control[[index]]
test <- l_test_0.6[[index]]

p1 <- plot_simulated_tree(control$tree,title='Control',legend.position='none')
p2 <- plot_simulated_tree(test$tree,title='Test',legend.position='none')
p <- plot_grid(p1,p2,nrow=1)
ggsave('figures/test3_tree_comparison.pdf',width=7,height=5)

info <- extract_data_for_seed(seed=333072, c.shape1=0.5, c.shape2=2.5)
p <- plot_mixtures(info$mix,'seed=333072')
ggsave('figures/test3_mixtures.pdf',width=7,height=5)



index <- x$i[2]
control <- l_control[[index]]
test <- l_test_0.6[[index]]

p1 <- plot_simulated_tree(control$tree,title='Control',legend.position='none')
p2 <- plot_simulated_tree(test$tree,title='Test',legend.position='none')
p <- plot_grid(p1,p2,nrow=1)
ggsave('figures/test4_tree_comparison.pdf',width=7,height=5)

info <- extract_data_for_seed(seed=515063, c.shape1=0.5, c.shape2=2.5)
p <- plot_mixtures(info$mix,'seed=515063')
ggsave('figures/test4_mixtures.pdf',width=7,height=5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how does allowing insertion/deletion bias affect results?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





