rm(list=ls())
library(polyGsim)

#sourceCpp('src/simulate.cpp')
#source('R/func.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulation 1: poly-G angular distance with impure bulk samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## this must be here only once and at the beginning. 
## I'm not sure why, but otherwise AD trees don't match
set.seed(123) 

## define following parameters
n_clones <- 10
n_markers <- 58
mu <- 1e-4
bdratio <- 1.01

## generate random genotype for the average normal (i.e. zygote)
normal <- get_normal_genotype(n_markers=n_markers)

## get a random number of generations from zygote to last tumor sample taken/patient death
max_gens <- random_generations(bdratio=bdratio)$generations

## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
tree <- get_clone_tree(n_clones)
png('figures/wiki/fig1.png'); plot(tree); dev.off()

## plot the true (un-knowable) chronology of the clones
p <- plot_chronology(tree,'Simulated clone/metastasis evolution',max_gens)
ggsave(plot=p,'figures/wiki/fig2.png',width=7,height=5)

## get random purities
pc <- get_purity(tree)

## generate marker biases (here, set to 0)
bias <- get_marker_bias(n_markers, bias.sd=0)

## get indels for the given tree
indels <- get_indels_for_tree(tree, max_gens, mu, n_markers, bias) 

## get genotypes for each clone
gt <- get_clone_genotypes(indels, normal)

## get mixing fractions matrix
mix <- get_mixing_proportions(pc)
p <- plot_mixtures(mix,'Random purity') 
ggsave(plot=p, 'figures/wiki/fig3.png',width=7,height=5)

## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)
ang::write_distance_matrix( round(head(t(ml), 10),5), 'figures/wiki/table1.txt')

## anonymize marker lengths so that each marker's minimum value is 0.
ml <- get_anonymized_marker_lengths(ml)

z <- get_angular_distance_matrix(ml,return_z=T)



## calculate the angular distance matrix from the anonymized marker lengths
ad <- get_angular_distance_matrix(ml)

## plot the angular distance neighbor joining tree
ad_tree <- nj(ad)
p <- plot_simulated_tree(ad_tree,title='Angular distance tree')
ggsave(plot=p, 'figures/wiki/fig4.pdf',width=7,height=5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2nd simulation:
# Some markers have strong bias to have either deletions or insertions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## this must be here only once and at the beginning. 
## I'm not sure why, but otherwise AD trees don't match
set.seed(123) 

## define following parameters
n_clones <- 10
n_markers <- 58
mu <- 1e-4
bdratio <- 1.01
bias.sd <- 0.03

## generate random genotype for the average normal (i.e. zygote)
normal <- get_normal_genotype(n_markers=n_markers)

## get a random number of generations from zygote to last tumor sample taken/patient death
max_gens <- random_generations(bdratio=bdratio)$generations

## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
tree <- get_clone_tree(n_clones)

## get random purities
pc <- get_purity(tree)

## generate marker biases
bias <- get_marker_bias(n_markers, bias.sd=bias.sd)

## plot the bias [TO DO]
bias$marker <- factor(bias$marker, levels=bias$marker)
p <- ggplot(bias, aes(x=marker, y=bias)) + 
    scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,by=0.25)) +
    geom_bar(stat='identity', fill='#bfbfbf', color='white', size=1.5) + 
    geom_hline(color='red',linetype='solid',size=0.25,yintercept=0.5) + 
    ang::theme_ang(base_size=12) +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    labs(x='Markers',y='Prob of insertions during indel',subtitle=paste0('Simulated marker bias:\nProb(insertion) ~ Beta(mean=0.50, SD=',bias.sd,')'))


## get indels for the given tree
indels <- get_indels_for_tree(tree, max_gens, mu, n_markers, bias) 

## get genotypes for each clone
gt <- get_clone_genotypes(indels, normal)

## get mixing fractions matrix
mix <- get_mixing_proportions(pc)

## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)

## anonymize marker lengths so that each marker's minimum value is 0.
ml <- get_anonymized_marker_lengths(ml)

## calculate the angular distance matrix from the anonymized marker lengths
ad <- get_angular_distance_matrix(ml)

## plot the angular distance neighbor joining tree
ad_tree <- nj(ad)
p <- plot_simulated_tree(ad_tree,title=paste0('Angular distance tree\n(biased markers: Prob insertion ~ Beta(mean=0.50, SD=',bias.sd,')'))
ggsave(plot=p, 'figures/wiki/fig5.pdf',width=7,height=5)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3rd simulation:
# Dominant clone in each sample is not 100% clonal
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(123) 

## define following parameters
n_clones <- 10
n_markers <- 58
mu <- 1e-4
bdratio <- 1.01

## generate random genotype for the average normal (i.e. zygote)
normal <- get_normal_genotype(n_markers=n_markers)

## get a random number of generations from zygote to last tumor sample taken/patient death
max_gens <- random_generations(bdratio=bdratio)$generations

## simulate a random phylogeny of primary tumor clones and (monoclonal) metastases
tree <- get_clone_tree(n_clones)

## plot the true (un-knowable) chronology of the clones
p <- plot_chronology(tree,'Simulated clone/metastasis evolution',max_gens)

## get random purities
pc <- get_purity(tree)

## generate marker biases (here, set to 0)
bias <- get_marker_bias(n_markers, bias.sd=0)

## get indels for the given tree
indels <- get_indels_for_tree(tree, max_gens, mu, n_markers, bias) 

## get genotypes for each clone
gt <- get_clone_genotypes(indels, normal)

## get mixing fractions matrix
mix <- get_mixing_proportions(pc)
p <- plot_mixtures(mix,'Random purity') 

## get the admixed marker lengths
ml <- get_mean_marker_lengths(gt, mix, n_markers)

## anonymize marker lengths so that each marker's minimum value is 0.
ml <- get_anonymized_marker_lengths(ml, n_markers=n_markers)

z <- get_angular_distance_matrix(ml,return_z=T)

mlm <- reshape2::melt(ml)
zm <- reshape2::melt(z)
dat <- merge(mlm, zm, by=c('Var1','Var2'))
dat <- as.data.table(dat)

p <- ggplot(dat, aes(x=value.x,y=value.y)) +
    geom_point() +
    geom_smooth(method='lm',se=T) +
    facet_wrap(facets=~Var1)

get_correlation <- function(dat) {
    x <- dat$value.x
    y <- dat$value.y
    r <- cor(x,y,method='pearson')
    fit <- lm(y ~ x)
    coefs <- summary(fit)$coefficients
    out <- as.list(coefs[2,])
    out <- append(out, list(r=r))
}
cors <- dat[,get_correlation(.SD),by=Var1]
cors <- cors[order(r),]
toadd <- cbind(Var1=rownames(pc), as.data.table(pc))
cors <- merge(cors, toadd, by='Var1')
cors <- cors[order(purity),]
cors[,lwr:=Estimate-1.96*`Std. Error`]
cors[,upr:=Estimate+1.96*`Std. Error`]




fit <- lm(value.y ~ value.x, data=dat[Var1=='P6'])
summary(fit)


library(philentropy)
jsd <- JSD(ml)
rownames(jsd) <- rownames(ml); colnames(jsd) <- rownames(ml)
tree1 <- nj(jsd)

tree2 <- nj(z)








l1 <- l1[!is.na(p),]
    mu <- mean(l1$p)
    qs <- as.numeric(quantile(l1$p,c(0.025,0.975)))
    list(sample=this.sample, mu=mu, lwr=qs[1], upr=qs[2], R=R)






set.seed(42)
samples <- rownames(ml)
samples <- samples[samples!='normal']
l4 <- lapply(samples, run_for_sample, ml, R=1e3, bootstrap=T, version=1)
res4 <- rbindlist(l4)
toadd <- cbind(sample=rownames(pc), as.data.table(pc))
res4 <- merge(res4, toadd, by='sample', all.x=T)
p4 <- ggplot(res4, aes(x=purity, y=mu)) +
    scale_x_continuous(expand=c(0,0),limits=c(0,1)) + scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
    geom_abline(intercept=0, slope=1, color='#bfbfbf', linetype='dashed', size=0.25) + 
    geom_errorbar(aes(min=lwr,max=upr),width=0.05,color='black') +
    geom_point(pch=19,color='black') +
    ang::theme_ang(base_size=12) +
    labs(x='True purity (simulated)', y='Predicted purity (inferred from data)',subtitle='Resampled markers, R=1e4')

p <- plot_grid(p3, p4)






set.seed(42)
samples <- rownames(ml)
samples <- samples[samples!='normal']
l <- lapply(samples, run_for_sample, ml, R=1e4)
res <- rbindlist(l)

toadd <- cbind(sample=rownames(pc), as.data.table(pc))
res <- merge(res, toadd, by='sample', all.x=T)
p <- ggplot(res, aes(x=purity, y=mu)) +
    scale_x_continuous(expand=c(0,0),limits=c(0,1)) + scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
    geom_abline(intercept=0, slope=1, color='#bfbfbf', linetype='dashed', size=0.25) + 
    geom_errorbar(aes(min=lwr,max=upr),width=0.05,color='black') +
    geom_point(pch=19,color='black') +
    ang::theme_ang(base_size=12) +
    labs(x='True purity (simulated)', y='Predicted purity (inferred from data)',subtitle='Resampled markers, R=1e4')




set.seed(42)
samples <- rownames(ml)
samples <- samples[samples!='normal']
l2 <- lapply(samples, run_for_sample, ml, R=1e4, bootstrap=F)
res2 <- rbindlist(l2)

toadd <- cbind(sample=rownames(pc), as.data.table(pc))
res2 <- merge(res2, toadd, by='sample', all.x=T)
p2 <- ggplot(res2, aes(x=purity, y=mu)) +
    scale_x_continuous(expand=c(0,0),limits=c(0,1)) + scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
    geom_abline(intercept=0, slope=1, color='#bfbfbf', linetype='dashed', size=0.25) + 
    geom_errorbar(aes(min=lwr,max=upr),width=0.05,color='black') +
    geom_point(pch=19,color='black') +
    ang::theme_ang(base_size=12) +
    labs(x='True purity (simulated)', y='Predicted purity (inferred from data)',subtitle='R=1e4 runs (no resampling)')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use constrained linear model to estimate the tumor purity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit_model_to_sample <- function(i, ml, penalty, return_corrected=F) {
    t1 <- ml[i,]
    n <- ml['normal',]
    y <- n - t1
    x <- n

    ## constrained gaussian regression: p (slope) must be between 0-1
    require(ConsReg)
    fit = ConsReg(formula = y ~ x, family = 'gaussian', optimizer = 'solnp', LOWER = 0.05, UPPER = 0.95, penalty=penalty)

    coefs <- coef(summary(fit))
    p <- coefs[2,1]
    intercept <- coefs[1,1]
    residuals <- resid(fit)
    x1 <- -1 * (residuals + intercept) / p
    
    if(return_corrected==F) {
        list(sample=rownames(ml)[i], p_infer=p)
    } else {
        list(sample=rownames(ml)[i], p_infer=p, t_pure=x1)
    }
}


run_with_method <- function(i, ml, pc, penalty) {
    if(i %% 10 == 0) message(i)
    true_purity <- pc[rownames(ml),]
    true_purity <- cbind(sample=rownames(true_purity), as.data.table(true_purity))
    setnames(true_purity,'purity','p_true')
    l <- lapply(1:(nrow(ml)-1), fit_model_to_sample, ml, penalty=penalty)
    result <- rbindlist(l)
    result <- merge(result, true_purity, by='sample', all.x=T)
    r = cor(result$p_true, result$p_infer, method='pearson')
    if(i > 0) {
        list(i=i, r=r, penalty=penalty)
    } else {
        result$penalty <- penalty
        result
    }
}

set.seed(42)
l1 <- rbindlist(lapply(1:1000, run_with_method, ml, pc, penalty=1000))
set.seed(42)
res1 <- run_with_method(0, ml, pc, penalty=1000)

p <- ggplot(l1, aes(x=r)) +
    geom_histogram(bins=50)
p











# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(42)
l2 <- rbindlist(lapply(1:10, run_with_method, ml, pc, penalty=10000))

set.seed(42)
l2 <- rbindlist(lapply(1:10, run_with_method, ml, pc, penalty=100000))










fit_model_to_sample <- function(i, ml, method, niter=NA, n.sim=NA) {
    t1 <- ml[i,]
    n <- ml['normal',]
    y <- n - t1
    x <- n

    ## constrained gaussian regression: p (slope) must be between 0-1
    require(ConsReg)
    if(method=='mcmc') {
        fit = ConsReg(formula = y ~ x, family = 'gaussian', optimizer = method, LOWER = 0.05, UPPER = 0.95,penalty=1000, niter=niter)
    } else {
        fit = ConsReg(formula = y ~ x, family = 'gaussian', optimizer = method, LOWER = 0.05, UPPER = 0.95,penalty=1000)
    }

    coefs <- coef(summary(fit))
    p <- coefs[2,1]
    intercept <- coefs[1,1]
    residuals <- resid(fit)
    x1 <- -1 * (residuals + intercept) / p
    list(sample=rownames(ml)[i], p_infer=p) #, t_pure=x1)
}


run_with_method <- function(i, method, niter, ml, pc) {
    if(i %% 10 == 0) message(i)
    true_purity <- pc[rownames(ml),]
    true_purity <- cbind(sample=rownames(true_purity), as.data.table(true_purity))
    setnames(true_purity,'purity','p_true')
    l <- lapply(1:(nrow(ml)-1), fit_model_to_sample, ml, method=method, niter=niter)
    result <- rbindlist(l)
    result <- merge(result, true_purity, by='sample', all.x=T)
    r = cor(result$p_true, result$p_infer, method='pearson')
    if(i > 0) {
        list(i=i, r=r, method=method, niter=niter)
    } else {
        result$method <- method
        result$niter <- niter
        result
    }
}

set.seed(42)
l6 <- rbindlist(lapply(1:10, run_with_method, 'mcmc', niter=1000, ml, pc))
set.seed(42)
res6 <- run_with_method(0, 'mcmc', niter=1000, ml, pc)


set.seed(42)
l8 <- rbindlist(lapply(1:100, run_with_method, 'mcmc', niter=10000, ml, pc))

set.seed(42)
l9 <- rbindlist(lapply(1:10, run_with_method, 'MCMCmetrop', niter=NA, ml, pc))

set.seed(42)
l10 <- rbindlist(lapply(1:10, run_with_method, 'GA', niter=NA, ml, pc))

set.seed(42)
l11 <- rbindlist(lapply(1:10, run_with_method, 'GenSA', niter=NA, ml, pc))




set.seed(42)
l10 <- rbindlist(lapply(1:10, run_with_method, 'GA', niter=NA, ml, pc))
set.seed(42)
res <- run_with_method(0, 'GA', niter=NA, ml, pc)


set.seed(42)
l8 <- rbindlist(lapply(1:10, run_with_method, 'mcmc', niter=20000, ml, pc))
set.seed(42)
res8 <- run_with_method(0, 'mcmc', niter=20000, ml, pc)

p8 <- ggplot(res8,aes(x=p_true,y=p_infer)) +
    geom_abline(intercept=0, slope=1) +
    geom_point()


set.seed(42)
l1 <- rbindlist(lapply(1:10, run_with_method, 'solnp', niter=NA, ml, pc))
set.seed(42)
res1 <- run_with_method(0, 'solnp', niter=NA, ml, pc)

p1 <- ggplot(res1,aes(x=p_true,y=p_infer)) +
    geom_abline(intercept=0, slope=1) +
    geom_point()








residuals - p * t1

residuals - p*n
residuals - p*n




library(nnls)



normal <- ml['normal',]
for(i in 1:nrow(ml)) ml[i,] <- ml[i,] - normal




## calculate the angular distance matrix from the anonymized marker lengths
ad <- get_angular_distance_matrix(ml)

## plot the angular distance neighbor joining tree
ad_tree <- nj(ad)
p <- plot_simulated_tree(ad_tree,title='Angular distance tree')
ggsave(plot=p, 'figures/wiki/fig4.pdf',width=7,height=5)


