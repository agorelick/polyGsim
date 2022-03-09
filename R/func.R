#library(Rcpp)
#library(phytools)
#library(TreeTools)
#library(ape)
#library(data.table)
#library(ggtree)


##' random_generations
##' @export
random_generations <- function(starting_cells=1, s=0.004, k=3, max_gens=1e4, max_cells=1e12) { 
    ## use simple branching process to simulate a number of cell divisions for cancer to reach 10^12 cells
    n_cells <- starting_cells
    gen = 0
    n_extinctions <- 0
    while(gen < max_gens & n_cells < max_cells) {
        death_prob <- 0.5*(1-s)^k
        dead_cells <- rbinom(n=1,size=n_cells,prob=death_prob)
        remaining_cells <- n_cells - dead_cells
        if(remaining_cells > 0) {
            n_cells <- remaining_cells * 2
        } else {
            n_cells <- 1
            gen <- 0
            n_extinctions <- n_extinctions + 1
        }
        gen <- gen+1
    }
    list(generations=gen, final_cells=n_cells, n_extinctions=n_extinctions)
}


##' sim_coalescence
##' @export
sim_coalescence <- function(n_samples,n_samples_is_average=T) {  

    ## pick a random number of samples with the given average lambda
    if(n_samples_is_average) n_samples <- rpois(1,lambda=n_samples)

    ## random coalescence tree using ape::rcoal
    samples <- paste0('s',1:n_samples)
    tree <- rcoal(n_samples, tip.label=samples)

    ## randomly decide a length for the initial tumor branch (from a normal distribution fit to the internal branch lengths)
    branch_lengths <- tree$edge.length
    log10_branch_lengths <- log10(branch_lengths)
    mu <- mean(log10_branch_lengths)
    s <- sd(log10_branch_lengths)
    edge.length <- 10^rnorm(mean=mu,sd=s,n=1)
    where <- length(tree$tip.label) + 1 ## where place normal branch (rcoal always puts the normal at length(tree)+1)

    ## add the normal branch to the tree    
    if(is.null(where)) where<-length(tree$tip)+1
    tip<-list(edge=matrix(c(2,1),1,2),
              tip.label='normal',
              edge.length=edge.length,
              Nnode=1)
    class(tip)<-"phylo"
    tree2 <- bind.tree(tree,tip,where=where)

    ## root the tree at the normal
    tree2 <- phytools::reroot(tree2, node.number=grep('normal',tree2$tip.label))
    return(tree2)
}


##' get_lineage_map
##' @export
get_lineage_map <- function(tree) { 
    samples <- tree$tip.label
    samples <- samples[samples!='normal']
    n_samples <- length(samples)
    ## get a vector of the nodes in the lineage leading to each sample
    get_lineage_for_tip <- function(node,tree,iter.max=1e3) {
        lineage <- node
        i <- 0
        while(i < iter.max & !is.null(node)) {
            trash <- capture.output(node <- phytools::getParent(tree, node))
            lineage <- c(node, lineage)
            i <- i+1
        }
        lineage
    }
    l <- lapply(1:n_samples, get_lineage_for_tip, tree)
    tobinary <- function(x) {
        out <- as.list(rep(1,length(x)))
        names(out) <- x
        out
    }
    tmp <- rbindlist(lapply(l, tobinary),fill=T)
    tmp <- tmp[,order(as.integer(colnames(tmp))),with=F]
    tmp[is.na(tmp)] <- 0
    tmp <- as.data.frame(tmp)
    rownames(tmp) <- samples
    tmp
}


##' chronoG
##' @export
chronoG <- function(n_samples, n_samples_is_average=T, starting_cells=1, s=0.004, k=3, max_gens=1e4, max_cells=1e12, n_markers=58, mu=1e-4) {

    ## simulate the cancer coalescence tree
    #sourceCpp('simulate.cpp')
    tree <- sim_coalescence(n_samples=n_samples,n_samples_is_average=n_samples_is_average); rm(n_samples)
    samples <- tree$tip.label
    samples <- samples[samples!='normal']
    n_samples <- length(samples)
    tree <- Preorder(tree) ## necessary

    ## generate a random number of cell divisions for the above phylogeny
    gens <- random_generations(starting_cells=starting_cells, s=s, k=k, max_gens=max_gens, max_cells=max_cells)
    n_gens <- gens$generations

    ## extract tree data
    d <- as.data.table(ggtree(tree)$data)
    d <- d[order(x,decreasing=F),]
    d <- d[-which(d$label=='normal'),]
    xpos <- d$x
    xpos <- xpos / max(xpos)
    d$gentime <- round(xpos * n_gens)
    d$add_gentime <- c(0,diff(d$gentime))
    d <- d[2:nrow(d),]

    ## get lineage-map
    map <- get_lineage_map(tree) # rows are the tips. columns are parent nodes
    colnames(map) <- paste0('n',colnames(map))

    ## loop across the cell divisions. keep track of the different ongoing genotypes
    ## gm is the generation-matrix tracking changes in lineages
    gm <- matrix(nrow=n_gens,ncol=n_samples)
    gm[is.na(gm)] <- as.integer(0)
    colnames(gm) <- samples
    for(g in 1:(n_gens-1)) {
        current_node <- d$node[d$gentime==g]

        if(length(current_node) > 0) {
            ## get the samples on either side of the bifurcation
            split <- d$node[d$parent==current_node]
            get_samples <- function(node, map) {
                which(map[[paste0('n',node)]]==1)
            }
            split <- lapply(split, get_samples, map)

            ## get the new levels to add to the matrix
            current_level <- max(gm)
            new_levels <- current_level + c(1,2)

            ## update the levels from this generation onward 
            gm[g:nrow(gm),split[[1]]] <- new_levels[1]
            gm[g:nrow(gm),split[[2]]] <- new_levels[2]
        }
    }

    ## gt is the genotype matrix
    gt <- matrix(nrow=n_samples,ncol=n_markers*2)
    gt[is.na(gt)] <- 0
    rownames(gt) <- colnames(gm)

    simulate_mutations_in_levels_for_gen <- function(n_levels, n_markers, mu) {
        tmp <- matrix(nrow=n_levels,ncol=n_markers*2)
        tmp[is.na(tmp)] <- 0
        tmp <- rcpp_mutate_length_matrix(tmp, mu, 1)
        tmp
    }

    limit <- n_gens - 1
    for(g in 1:limit) {
        #message(g)
        levels <- as.integer(gm[g,])
        unique_levels <- unique(levels)
        nl <- length(unique_levels)
        tmp <- simulate_mutations_in_levels_for_gen(nl, n_markers, mu)
        i <- 1
        for(i in 1:nl) { 
            lev <- unique_levels[i]
            affected_samples <- which(levels==lev)
            for(s in affected_samples) {
                gt[s,] <- gt[s,] + tmp[i,]    
            }
        } 
    }


    ## collapse the genotypes to a mean marker length (in tumor) for each sample
    row_to_mean_length <- function(x) {
        n_markers <- length(x) / 2
        x1 <- x[1:n_markers]
        x2 <- x[(n_markers+1):length(x)]
        mu <- (x1 + x2)/2 
        mu
    }
    gt <- t(apply(gt, 1, row_to_mean_length))
    ## add normal genotype (all 0)
    normal <- rep(0,ncol(gt))
    gt <- rbind(gt, normal)

    #dm <- dist(gt, method='euclidean')
    #tree2 <- upgma(dm)
    #tree2 <- nj(dm)
    #tree2 <- phytools::reroot(tree2, node.number=grep('normal',tree2$tip.label))

    params <- list(n_samples=n_samples, n_gens=n_gens, n_samples_is_average=n_samples_is_average,  starting_cells=starting_cells, s=s, k=k, max_gens=max_gens, max_cells=max_cells, n_markers=n_markers, mu=mu) 

    list(tree=tree, gt=gt, gm=gm, map=map, max_gens=max_gens, params=params)

}



