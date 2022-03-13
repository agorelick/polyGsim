
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
    map <- rbindlist(lapply(l, tobinary),fill=T)
    map <- map[,order(as.integer(colnames(map))),with=F]
    map[is.na(map)] <- 0
    map <- as.data.frame(map)
    rownames(map) <- samples
    colnames(map) <- paste0('n',colnames(map))
    map
}


##' chronoG
##' @export
chronoG <- function(n_samples, n_samples_is_average=T, starting_cells=1, s=0.004, k=3, max_gens=1e4, max_cells=1e12, n_markers=58, mu=1e-4, beta_params=c(2,2)) {

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

    ## loop across the cell divisions. keep track of the different ongoing genotypes
    ## gm is the generation-matrix tracking changes in lineages
    gm <- matrix(nrow=n_gens,ncol=n_samples)
    gm[is.na(gm)] <- as.integer(0)
    colnames(gm) <- samples
    for(g in 1:(n_gens-1)) {
        current_node <- d$node[d$gentime==g]

        if(length(current_node) > 0) {
            ## get the samples on either side of the bifurcation
            split <- d$node[d$parent %in% current_node]
            get_samples <- function(node, map) {
                which(map[[paste0('n',node)]]==1)
            }
            split <- lapply(split, get_samples, map)

            ## get the new levels to add to the matrix
            current_level <- max(gm)
            new_levels <- current_level + 1:length(split)

            ## update the levels from this generation onward 
            for(this.level in 1:length(new_levels)) {
                gm[g:nrow(gm),split[[this.level]]] <- new_levels[this.level]
            }
        }
    }

    #browser()
    ## generate a random normal genotype
    normal <- sample(10:25,replace=T,size=n_markers*2)

    ## gt is the genotype matrix
    gt <- matrix(nrow=n_samples,ncol=n_markers*2)
    for(i in 1:nrow(gt)) gt[i,] <- normal #gt[is.na(gt)] <- 0
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
            for(sa in affected_samples) {
                gt[sa,] <- gt[sa,] + tmp[i,]    
            }
        } 
    }
    rm(sa)

    ## add the normal to the last row
    gt <- rbind(gt, normal)
    #browser()

    ## collapse the genotypes to a mean marker length (in tumor) for each sample
    row_to_mean_length <- function(x) {
        n_markers <- length(x) / 2
        x1 <- x[1:n_markers]
        x2 <- x[(n_markers+1):length(x)]
        mu <- (x1 + x2)/2 
        mu
    }
    gt <- t(apply(gt, 1, row_to_mean_length))

    ## randomly generate tumor purities
    purities <- matrix(ncol=1,nrow=n_samples)
    purities[,1] <- rbeta(n=n_samples, shape1=beta_params[1], shape2=beta_params[2])
    rownames(purities) <- samples
    colnames(purities) <- 'purity'

    ## create 'observed' genotype data, which is the underlying tumor marker lengths scaled by their purities. This can be treated as the observed difference in mean-marker-lengths from the normal.
    gt_obs <- copy(gt)
    normal <- as.numeric(gt_obs['normal',])
    for(sa in samples) gt_obs[sa,] <- ( gt[sa,] * purities[sa,1] ) + (normal * (1 - purities[sa,1]))
    rm(sa)

    gt_obs_centered <- copy(gt_obs)
    marker_min_meanlengths <- apply(gt_obs_centered, 2, min)
    for(mi in 1:n_markers) {
        gt_obs_centered[,mi] <- gt_obs_centered[,mi] - marker_min_meanlengths[mi]
    }
    rm(mi)

    params <- list(n_samples=n_samples, n_gens=n_gens, n_cells=gens$final_cells, extinctions=gens$n_extinctions, n_samples_is_average=n_samples_is_average,  starting_cells=starting_cells, s=s, k=k, max_gens=max_gens, max_cells=max_cells, n_markers=n_markers, mu=mu) 

    list(tree=tree, gt=gt, gt_obs=gt_obs, gt_obs_centered=gt_obs_centered, purities=purities, gm=gm, map=map, params=params)
}


##' get_angular_distance_matrix
##' @export
get_angular_distance_matrix <- function(d) {
    ## get meanlength minus normal for sim data
    ## d should be a matrix of mean-lengths where rows are samples and columns are markers

    normal <- d['normal',]
    for(i in 1:nrow(d)) d[i,] <- d[i,] - normal
    colnames(d) <- paste0('m',1:ncol(d))
    markerset <- colnames(d); 
    usedsamps <- rownames(d)
    usedsamps <- usedsamps[usedsamps!='normal']
    meanlen_diff_normal <- as.data.frame(d)

    Zij <- function(usedsamps,markerset,meanlen_diff_normal) {
        #* for each used tumor sample calculate denom = sqrt(sum over markers(Y^2)) and Z = Y/denom
        z <- list()
        denom <- numeric()
        for (t in usedsamps) {
            denom[t] <- 0
            for (m in markerset) {
                denom[t] <- denom[t] + (meanlen_diff_normal[t,m]) ^ 2
            }
            denom[t] <- sqrt(denom[t])
            for (m in markerset) {
                z[[m]][t] <- meanlen_diff_normal[t,m]/denom[t]
            }
        }
        z
    }
    Z <- Zij(usedsamps, markerset, meanlen_diff_normal)

    ang_dist_matrix <- function(Z,usedsamps,ns,markerset,sel_normal_sample,power=1) {
        # for each pair of samples in Z, calculate angular distance (exclude excluded samples)
        # angular distance matrix
        angular_dist <- matrix(0,nrow=ns,ncol=ns)
        for (s1 in 1:(ns-1)) {
            for (s2 in (s1+1):ns) {
                for (m in markerset) {
                    angular_dist[s1,s2] <- angular_dist[s1,s2] + Z[[m]][usedsamps[s1]]*Z[[m]][usedsamps[s2]]
                }
                angular_dist[s1,s2] <- (acos(pmin(pmax(angular_dist[s1,s2],-1.0),1.0)))^power
                angular_dist[s2,s1] <- angular_dist[s1,s2]
            }
        }
        colnames(angular_dist) <- usedsamps
        rownames(angular_dist) <- usedsamps

        #* add a column and row for the normal
        angular_dist_w_root <- rbind(angular_dist,rep((pi/3)^power,ns))
        angular_dist_w_root <- cbind(angular_dist_w_root,c(rep((pi/3)^power,ns),0))
        colnames(angular_dist_w_root) <- c(usedsamps,sel_normal_sample)
        rownames(angular_dist_w_root) <- c(usedsamps,sel_normal_sample)
        angular_dist_w_root
    }

    dm <- ang_dist_matrix(Z, usedsamps, ns=length(usedsamps), markerset, sel_normal_sample='normal', power=1)
    dm
}


##' load_marker_files
##' @export
load_marker_files <- function(marker_dir,return_mean_length=T) {
    polyg_dir <- dir(marker_dir,full.names=T)
    get_lengths <- function(i, polyg_dir) {
        files <- polyg_dir[i]
        load_file <- function(file) {
            marker <- tail(strsplit(file,'[/]')[[1]],1)
            marker <- strsplit(marker,'_')[[1]][1]        
            x <- fread(file)
            x <- cbind(length=1:nrow(x),x)
            x <- melt(x, id.vars=c('length')) 
            x$marker <- marker
            x
        }
        l <- lapply(files, load_file)
        out <- rbindlist(l,fill=T)
        out
    }
    l <- lapply(1:length(polyg_dir), get_lengths, polyg_dir)
    l <- rbindlist(l)

    if(return_mean_length) {
        f=function(s) paste(head(strsplit(s,'[_]')[[1]],-1),collapse='_')
        l$sample <- sapply(as.character(l$variable), f)
        l <- l[!is.na(length),]
        get_mean <- function(l) {
            mu=sum(l$length * l$value,na.rm=T) / sum(l$value,na.rm=T)
            list(mu=mu)
        }
        l <- l[,get_mean(.SD),by=c('sample','marker')]
        d <- dcast(sample ~ marker, value.var='mu', data=l)
    } else {
        d <- copy(l)
    }
    d
}





