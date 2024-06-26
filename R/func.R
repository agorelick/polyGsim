##' get_normal_genotype
##' @export
get_normal_genotype <- function(n_markers, diploid=T) {
    normal1 <- sort(sample(10:25,replace=T,size=n_markers))
    normal2 <- sort(sample(10:25,replace=T,size=n_markers))

    if(diploid) {
        marker_names <- c(paste0('m',1:n_markers,'.1'), paste0('m',1:n_markers,'.2'))
        #normal2 <- normal1 + sample(-3:3, replace=T, size=n_markers) 
        normal <- c(normal1, normal2)
        normal <- t(as.matrix(normal))
    } else {
        marker_names <- paste0('m',1:n_markers,'.1')
        normal <- t(as.matrix(normal1))
    }

    ## normal is the average germline
    colnames(normal) <- marker_names
    normal
}


##' random_generations
##' @export
random_generations <- function(starting_cells=1, bdratio=1.01, max_gens=1e4, max_cells=1e12) { 
    ## use birth-date discrete time branch process to simulate a number of cell divisions for cancer to reach 10^12 cells
    n_cells <- starting_cells
    gen = 0
    n_extinctions <- 0
    death_prob <- 1/(bdratio+1)
    
    while(gen < max_gens & n_cells < max_cells) {
        dead_cells <- rbinom(n=1,size=n_cells,prob=death_prob)
        remaining_cells <- n_cells - dead_cells        
        if(remaining_cells > 0) {
            n_cells <- remaining_cells * 2
            gen <- gen+1
        } else {
            gen <- 0
            n_cells <- 1
            n_extinctions <- n_extinctions + 1
            dead_cells <- 0
            remaining_cells <- 0
        }
    }
    
    list(generations=gen, final_cells=n_cells, n_extinctions=n_extinctions)
}


##' get_clone_tree
##' @export
get_clone_tree <- function(n_clones) { 

    if(n_clones < 4) stop('n_clones must be at least 4!')

    random_chronology <- function(n_clones,ancestor='normal') {  
        ## random tree toplogy
        clones <- paste0('s',1:n_clones)
        tree <- rtree(n_clones, tip.label=clones)

        ## randomly decide a length for the initial tumor branch (from a log-normal distribution fit to the internal branch lengths)
        branch_lengths <- tree$edge.length
        log10_branch_lengths <- log10(branch_lengths)
        mu <- mean(log10_branch_lengths)
        s <- sd(log10_branch_lengths)
        edge.length <- 10^rnorm(mean=mu,sd=s,n=1)
        where <- length(tree$tip.label) + 1 ## where to place normal branch

        ## add the normal branch to the tree    
        if(is.null(where)) where<-length(tree$tip)+1
        tip<-list(edge=matrix(c(2,1),1,2),
                  tip.label=ancestor,
                  edge.length=edge.length,
                  Nnode=1)
        class(tip)<-"phylo"
        tree2 <- bind.tree(tree,tip,where=where)

        ## root the tree at the normal
        tree2 <- phytools::reroot(tree2, node.number=which(tree2$tip.label==ancestor))
        return(tree2)
    }

    ## randomly pick number of PT and Met clones such that always at least 2 PT and 2 mets in every tree
    n_pt_clones <- sample(2:(n_clones-2),1) 
    n_met_clones <- n_clones - n_pt_clones

    ## iteratively select the number of met clones across met-clades to give the total number of met_clones
    remaining_met_clones <- n_met_clones
    met_clade_sizes <- c()
    i <- 0
    while(i < 100 & remaining_met_clones > 0) { 
        n_mets_in_clade <- sample(1:remaining_met_clones,1)
        met_clade_sizes <- c(met_clade_sizes, n_mets_in_clade)
        remaining_met_clones <- remaining_met_clones - n_mets_in_clade
    }
    met_clade_sizes
    n_met_clades <- length(met_clade_sizes)

    ## get the PT clone tumor, with 1 additional tip for each met clade
    ptree <- random_chronology(n_pt_clones+n_met_clades)
    ptree$tip.label <- gsub('s','P',ptree$tip.label)

    ## randomly pick which tips are where the met clade is joined
    tumor_tips <- grep('^P',ptree$tip.label,value=T)
    met_clade_roots <- sample(tumor_tips,n_met_clades)

    ## for each met-seeding clone, get a new random tree of mets, then bind it on 
    growing_tree <- copy(ptree)

    for(i in 1:n_met_clades) { 
        current_met_clade_size <- met_clade_sizes[i]
        current_met_root <- met_clade_roots[i]

        if(current_met_clade_size >= 3) {
            ## generate a new random tree for this clade
            mtree <- ape::rtree(current_met_clade_size,rooted=T)
            mtree$tip.label <- gsub('^t','M',mtree$tip.label)

            ## update its labels so that they don't conflict with the tips already in the growing tree
            tmp <- data.table(sample=mtree$tip.label)
            tmp$pos <- 1:nrow(tmp)
            tmp <- tmp[order(pos),]
            current_met_clones <- grep('^M',growing_tree$tip.label,value=T)
            increment <- ifelse(length(current_met_clones) > 0, max(as.integer(gsub('M','',current_met_clones))), 0)
            tmp[grep('^M',sample),new_name:=paste0('M',pos+increment)]
            tmp[is.na(new_name),new_name:=sample]
            mtree$tip.label <- tmp$new_name       

            ## attach the new met-clade tree onto the growing tree at the right position
            growing_tree <- bind.tree(growing_tree, mtree, where=which(growing_tree$tip.label==current_met_root))

        } else if(current_met_clade_size==2) {
            ## generate a new random tree for this clade
            mtree <- ape::rtree(3,rooted=T)
            n_to_drop <- 3 - current_met_clade_size
            dropped_tip <- sample(mtree$tip.label,n_to_drop)
            mtree <- drop.tip(mtree, dropped_tip)
            mtree$tip.label <- gsub('^t','M',mtree$tip.label)

            ## update its labels so that they don't conflict with the tips already in the growing tree
            tmp <- data.table(sample=mtree$tip.label)
            tmp$pos <- 1:nrow(tmp)
            tmp <- tmp[order(pos),]
            current_met_clones <- grep('^M',growing_tree$tip.label,value=T)
            increment <- ifelse(length(current_met_clones) > 0, max(as.integer(gsub('M','',current_met_clones))), 0)
            tmp[grep('^M',sample),new_name:=paste0('M',pos+increment)]
            tmp[is.na(new_name),new_name:=sample]
            mtree$tip.label <- tmp$new_name       

            ## attach the new met-clade tree onto the growing tree at the right position
            growing_tree <- bind.tree(growing_tree, mtree, where=which(growing_tree$tip.label==current_met_root))

        } else if(current_met_clade_size==1) {
            convert_to_met <- which(growing_tree$tip.label==current_met_root)
            
            tmp <- data.table(sample='M1')
            tmp$pos <- 1:nrow(tmp)
            current_met_clones <- grep('^M',growing_tree$tip.label,value=T)
            increment <- ifelse(length(current_met_clones) > 0, max(as.integer(gsub('M','',current_met_clones))), 0)
            tmp[grep('^M',sample),new_name:=paste0('M',pos+increment)]
            tmp[is.na(new_name),new_name:=sample]
            new_met_name <- tmp$new_name 
            growing_tree$tip.label[convert_to_met] <- new_met_name      
        }
    }

    ## remake PT clone labels (so that no gaps in numbers)
    tmp <- data.table(pt_clone=grep('^P',growing_tree$tip.label,value=T))
    tmp$pos <- 1:nrow(tmp)
    tmp$newpos <- nrow(tmp):1
    tmp$newname <- paste0('p',tmp$newpos)
    for(s in 1:nrow(tmp)) growing_tree$tip.label[which(growing_tree$tip.label==tmp$pt_clone[s])] <- tmp$newname[s]

    updated_tips <- grep('^p',growing_tree$tip.label)    
    growing_tree$tip.label[updated_tips] <- toupper(growing_tree$tip.label[updated_tips])
    growing_tree
}



##' plot_chronology
##' @export
plot_chronology <- function(tree,title,max_gens) {
    groups <- data.table(label=tree$tip.label)
    groups[grep('^P',label),group:='Primary']
    groups[grep('^M',label),group:='Metastasis']
    groups[grep('normal',label),group:='Normal']
    groups$group <- factor(groups$group,levels=c('Normal','Primary','Metastasis'))
    groups[group=='Normal',label:='N1']
    tree$tip.label[tree$tip.label=='normal'] <- 'N1'
    cols <- c('#008c45','#fab31d','black')
    names(cols) <- c('Primary','Metastasis','Normal')

    p <- ggtree(tree,layout='rect') %<+% groups
    p$data$label <- paste0(' ',p$data$label)
    p$data$x <- max_gens * p$data$x / max(p$data$x)
    p <- p + polyGsim_theme(base_size=12) + geom_tiplab(aes(color=group),angle=0)
    p$data$node_lab <- as.character(NA)
    p$data$node_lab[p$data$isTip==F & p$data$x < 0.98*max_gens] <- round(p$data$x[p$data$isTip==F & p$data$x < 0.98*max_gens])
    p <- p + geom_text(aes(label=node_lab),angle=0,size=3,color='blue',hjust=-0.1)
    p <- p + theme(legend.position='none', axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) 
    p <- p + scale_color_manual(values=cols,name='Tissue')
    p <- p + ggplot2::labs(x='Generations from zygote')
    p <- p + ggplot2::ggtitle(title)
    p 
}


##' get_marker_bias
##' @export
get_marker_bias <- function(n_markers, bias.sd) { 
    ## introduce bias for insertions or deletions in each marker
    if(bias.sd > 0.5 | bias.sd < 0 | is.na(bias.sd)) stop('bias.sd must be in [0-0.5]')
    tmp <- data.table(marker=paste0('m',1:n_markers))
    get_bias <- function(n, s) {
        ## calculate alpha/beta shape parameters to ensure mean=0.5 and SD=bias.sd
        v <- s^2
        a <- ((1/v)-4)/8
        b <- a
        sort(rbeta(n=n, shape1=a, shape2=a),decreasing=F)
    }
    tmp$bias <- get_bias(n=n_markers,bias.sd)
    tmp
}


##' get_indels_for_tree
##' @export
get_indels_for_tree <- function(tree,max_gens,n_markers,mu=1e-4,bias=rep(0, n_markers),max_ploidy=2) { 

    ## extract tree data
    d <- as.data.table(ggtree(tree)$data)
    d <- d[order(x,decreasing=F),]
    d <- d[-which(d$label=='normal'),]
    xpos <- d$x
    xpos <- xpos / max(xpos)
    d$gen_time <- round(xpos * max_gens)
    d$gens_to_run <- c(0,diff(d$gen_time))
    #d[gens_to_run==0 & is.Tip==T, gens_to_run:=1]
    d <- d[gens_to_run >= 0,]

    ## for each branch, run simulations for the according number of generations.
    simulate_indels_over_gens_for_branch <- function(gens,mu,n_markers,max_ploidy,biases) {
        gt_0 <- t(as.matrix(rep(0,n_markers*max_ploidy))) ## initial genotype for each branch
        gt_T <- rcpp_mutate_length_matrix(gt_0, mu, gens, biases=c(biases,biases))
        out <- as.list(gt_T[1,])
    }
    indels_for_each_branch <- lapply(d$gens_to_run, simulate_indels_over_gens_for_branch, mu, n_markers, max_ploidy, bias)
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
    toadd[1,1:(n_markers*max_ploidy)] <- 0
    out <- rbind(out, toadd)
    rownames(out)[nrow(out)] <- 'normal'

    ## add marker names
    n <- expand.grid(marker=seq(1,n_markers),copy=seq(1,max_ploidy))
    n$label <- paste0('m',n$marker,'.',n$copy)
    colnames(out) <- n$label
    out
}


 
##' get_clone_genotypes
##' @export
get_clone_genotypes <- function(indels, normal) { 
    clones <- rownames(indels)[rownames(indels)!='normal']
    n_clones <- length(clones)
    nc <- ncol(indels)

    ## gt is the genotype matrix, each tumor starts with 1st cancer cell
    gt <- matrix(nrow=n_clones,ncol=nc)
    for(i in 1:nrow(gt)) gt[i,] <- as.integer(normal)
    rownames(gt) <- clones
    gt <- rbind(gt, normal)
    rownames(gt)[nrow(gt)] <- 'normal'

    ## add the indels encountered to each clone's original genotype
    gt <- gt + indels
    gt
}


##' get_purity
##' @export
get_purity <- function(tree, purity_min=NA, purity_max=NA, clonality_min=NA, clonality_max=NA) {
    samples <- tree$tip.label
    n_samples <- length(samples)

    if(!is.na(purity_min) & !is.na(purity_max)) {
        ## randomly generate tumor purities, uniformly distributed between our allowed range
        purities <- runif(n=n_samples, min=purity_min, max=purity_max)
    } else {
        purities <- rep(1, n_samples)    
    }

    if(!is.na(clonality_min) & !is.na(clonality_max)) {
        ## now we allow mixing of multiple clones within each sample, where the dominant clone has uniformly distributed clonality
        clonalities <- runif(n=n_samples, min=clonality_min, max=clonality_max)
    } else {
        clonalities <- rep(1, n_samples)    
    }

    out <- data.frame(purity=purities, clonality=clonalities)
    rownames(out) <- samples
    out['normal','purity'] <- 0
    out['normal','clonality'] <- 0
    out
}



##' get_mixing_proportions
##' @export
get_mixing_proportions <- function(pc, even_mixing=F) {

    get_sample_mixture <- function(s, pc, even_mixing) {
        ## for each sample, take the non-normal fraction and split it over the tumor samples

        samples <- rownames(pc)
        nonself_samples <- samples[!samples %in% c('normal',s)]
        prop_normal <- 1 - pc[s,'purity']
        prop_cancer <- 1 - prop_normal
        clonality <- pc[s,'clonality']
        prop_self <- clonality * prop_cancer
        prop_others <- prop_cancer - prop_self
        if(even_mixing==F) {
            tmp <- runif(0,1,n=length(nonself_samples))
            prop_others <- prop_others * tmp / sum(tmp)
        } else {
            prop_others <- prop_others / length(nonself_samples)
            prop_others <- rep(prop_others, length(nonself_samples))
        }

        out_samples <- c('normal',s,nonself_samples)   
        props <- c(prop_normal, prop_self, prop_others)
        out <- data.table(sample=out_samples, prop=props)
        out$self <- s
        out$clonality <- clonality
        out$prop_cancer <- prop_cancer
        out <- out[!(sample=='normal' & prop==0),]
        out
    }

    pc['normal','clonality'] <- 0
    samples <- rownames(pc)
    l <- lapply(samples, get_sample_mixture, pc, even_mixing)
    d <- rbindlist(l) 
    out <- dcast( self + prop_cancer + clonality ~ sample, value.var='prop', data=d)
    out$normal[is.na(out$normal)] <- 0
    setnames(out,'prop_cancer','purity')
    out
}


##' plot_mixtures
##' @export
plot_mixtures <- function(x,title='') {
    x <- x[order(purity,decreasing=F),]
    pd <- x[self!='normal']
    self_levels <- pd$self
    sample_levels <- (c('normal',pd$self))
    pd <- melt(pd, id.vars=c('self','purity','clonality'))
    pd$self <- factor(pd$self, levels=self_levels)
    pd$variable <- factor(pd$variable,levels=sample_levels)
    cols <- c('white', rep(brewer.pal(8,'Accent'), length.out=nrow(x)-1))
    names(cols) <- sample_levels
    label_data <- data.table(self=x$self, value=x$purity)
    label_data$self <- factor(label_data$self, levels=self_levels)
    label_data <- label_data[!is.na(self),]
    label_data$label <- paste0(prettyNum(100*label_data$value, digits=2),'%')
    p <- ggplot(pd, aes(x=self,y=value)) +
        scale_y_continuous(expand=c(0,0)) + 
        geom_bar(stat='identity',aes(fill=variable),color='black',linewidth=0.25) +
        geom_text(data=label_data,aes(label=label),vjust=-0.15,size=3) +
        scale_fill_manual(values=cols,name='Origin') +
        polyGsim_theme(base_size=12) +
        labs(x='Sample',y='Mixture proportion', subtitle=title)
    p
}

 

##' get_mean_marker_lengths
##' @export
get_mean_marker_lengths <- function(gt, mix, n_markers, noise_cv=NA) { 

    get_mean_lengths_per_sample <- function(s, samples, gt, x) { 
        ## get admixed observed genotypes based on purity, clonality, and copy number alterations
        props <- as.matrix(x[,(samples),with=F])
        rownames(props) <- x$self
        props <- props[s,]
        props <- data.table(gt=names(props), prop=props) 
        tmp <- cbind(gt=rownames(gt), as.data.table(gt))
        tmp <- merge(props, tmp, by='gt', all=T)
        tmp <- melt(tmp, id.vars=c('gt','prop'))
        tmp <- tmp[!is.na(value),]
        f <- function(s) strsplit(s,'[.]')[[1]][1]
        tmp$marker <- sapply(as.character(tmp$variable), f)
        collapse <- function(tmp) {
            mu <- sum(tmp$value * tmp$prop) / sum(tmp$prop)
            list(mu=mu)
        }
        res <- tmp[,collapse(.SD),by=c('marker')]
        res$sample <- s
        res
    }
    samples <- rownames(gt)
    l <- lapply(samples, get_mean_lengths_per_sample, samples, gt, mix)
    l <- rbindlist(l)

     ## markers vary in their noise, so get a random CV per marker and use it to simulate noise for all samples
    if(!is.na(noise_cv)) {
        #message('Adding noise to mean lengths after mixing.')
        noise_s <- l$mu * noise_cv
        l$mu <- rnorm(mean=l$mu, sd=noise_s, n=nrow(l))
    }

    l <- dcast(sample ~ marker, value.var='mu', data=l)
    markers <- paste0('m',1:n_markers)
    l <- l[,c('sample',markers),with=F]

    ## convert to a matrix and return it
    rows <- l$sample
    l <- l[, c(2:ncol(l)), with = F]
    m <- as.matrix(l)
    rownames(m) <- rows
    m
}


##' get_anonymized_marker_lengths
##' @export
get_anonymized_marker_lengths <- function(d) { 
    ## anonymize by subtracting the minimum length for each marker 

    n_markers <- ncol(d)
    marker_min_meanlengths <- apply(d, 2, min)
    for (mi in 1:n_markers) {
        d[, mi] <- d[, mi] - marker_min_meanlengths[mi]
    }
    d
}


##' get_angular_distance_matrix
##' @export
get_angular_distance_matrix <- function(d,return_z=F) {
    ## get meanlength minus normal for sim data
    ## d should be a matrix of mean-lengths where rows are samples and columns are markers

    normal <- d['normal',]
    for(i in 1:nrow(d)) d[i,] <- d[i,] - normal
    markerset <- colnames(d)
    #markerset <- names(which(colMeans(d==0) != 1))
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

    Z <- Zij(usedsamps, markerset, meanlen_diff_normal)
    if(return_z) {
        markers <- names(Z)
        out <- as.matrix(rbindlist(lapply(Z, as.list)))
        rownames(out) <- markers
        out <- t(out)
        normal <- out[1,]
        normal[normal!=0] <- 0
        out <- rbind(out, normal)
        out
    } else {
        ang_dist_matrix(Z, usedsamps, ns=length(usedsamps), markerset, sel_normal_sample='normal', power=1)
    }
}


##' plot_simulated_tree
##' @export
plot_simulated_tree <- function(tree,layout='ape',title=NA,purities=NULL,legend.position='right', base_size=12) {
    if(!is.rooted(tree)) tree <- phytools::reroot(tree, node.number=which(tree$tip.label=='normal'))
    groups <- data.table(label=tree$tip.label)
    groups[grep('^P',label),group:='Primary']
    groups[grep('^M',label),group:='Metastasis']
    groups[grepl('N1',label) | grepl('normal',label),group:='Normal']
    groups$group <- factor(groups$group,levels=c('Normal','Primary','Metastasis'))

    if(!is.null(purities)) {
        toadd <- data.table(label=rownames(purities),purity=purities[,1])
        groups <- merge(groups, toadd, by='label', all.x=T)
        groups[group=='Normal',purity:=NA]
    }
    tree$tip.label[tree$tip.label=='normal'] <- 'N1'
    groups[label=='normal',label:='N1']
    cols <- c('#008c45','#fab31d','black')
    names(cols) <- c('Primary','Metastasis','Normal')
    p <- ggtree(tree,layout=layout) %<+% groups
    p <- p + geom_tiplab(aes(color=group),angle=0) +
        polyGsim_theme(base_size=base_size) +
        theme(legend.position=legend.position,
              axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank()) + 
        scale_color_manual(values=cols,name='Tissue')
    if(!is.null(purities)) {
        p <- p + geom_tippoint(aes(fill=purity),pch=21,size=2,stroke=0.25)
        p <- p + scale_fill_gradient2(low='blue',mid='white',high='red',na.value='black',midpoint=0.5,name='Purity',limits=c(0,1)) 
    }
    if(!is.na(title)) p <- p + ggtitle(title)
    p 
}


##' get_mean_marker_length_matrix
##' @export
get_mean_marker_length_matrix <- function(gt,n_markers) {
    f <- function(x,n_markers) {
        max_ploidy <- length(x) / n_markers        
        mat <- matrix(x, nrow=max_ploidy, ncol=n_markers, byrow=T)    
        mean_lengths <- colMeans(mat,na.rm=T)
        mean_lengths
    }
    tmp <- t(apply(gt, 1, f, n_markers))
    tmp
}


##' polyGsim_theme
##' @export
polyGsim_theme <- function (base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_bw(base_size = base_size, base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(line = element_line(colour = "black", linewidth = base_line_size, linetype = 1, lineend = "round"), text = element_text(colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = F), axis.text = element_text(colour = "black", size = rel(0.8)), axis.ticks = element_line(colour = "black", linewidth = rel(1)), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", linewidth = rel(1)), legend.key = element_blank(), strip.background = element_blank())
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions to potentially move out/delete
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##' get_ovoid_dimensions
get_ovoid_dimensions <- function(n_cells) { 
    longest_dimension_cm <- get_longest_dimension(n_cells)
    volume_cm3 <- get_volume(longest_dimension_cm)
    surface_area_cm2 <- get_surface_area(longest_dimension_cm)
    cells_per_cm2 <- 215443.5
    n_cells_on_surface <- surface_area_cm2 * cells_per_cm2
    list(longest_dimension_cm=longest_dimension_cm, volume_cm3=volume_cm3, 
         surface_area_cm2=surface_area_cm2, n_cells_on_surface=n_cells_on_surface)
}


##' load_marker_files
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


##' get_genotypes_with_cnas
get_genotypes_with_cnas <- function(tree, max_gens, n_markers, mu.indel, mu.cna, gens_until_first_cancer_cell, prop_markers_with_early_cnas=0.5) {
    ## DOES NOT WORK YET!!

    clones <- tree$tip.label[tree$tip.label!='normal']
    n_clones <- length(clones)
    marker_names <- c(paste0('m',1:n_markers,'.1'), paste0('m',1:n_markers,'.2'))

    ## normal is the average germline
    normal <- t(as.matrix(sample(10:25,replace=T,size=n_markers*2)))
    colnames(normal) <- marker_names

    ## ancestor is based on the germline after large number of divisions, and an early WGD event
    first_cancer_cell <- rcpp_mutate_length_matrix(copy(normal), mu.indel, gens=gens_until_first_cancer_cell)
    colnames(first_cancer_cell) <- marker_names


    indels <- get_indels_for_tree(tree, max_gens, mu.indel, n_markers, max_ploidy=4) 

    ## get scnas for diploid (-1 = deletion of that copy, +1 = duplication of that copy)
    scnas <- get_indels_for_tree(tree, max_gens, mu.cna, n_markers, max_ploidy=2)
    scnas[scnas < -1] <- -1     
    scnas[scnas > 1] <- 1     

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # genotypes for early SCNAs. SCNAs first, then indels
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## get length for each marker (diploid) in each sample with indels
    nc <- ncol(normal)
    gt_early_scnas <- matrix(nrow=n_clones,ncol=nc)
    for(i in 1:nrow(gt_early_scnas)) gt_early_scnas[i,] <- as.integer(first_cancer_cell)
    rownames(gt_early_scnas) <- clones
    gt_early_scnas <- rbind(gt_early_scnas, normal)
    rownames(gt_early_scnas)[nrow(gt_early_scnas)] <- 'normal'

    ## add copy number deletions to gt, then add in the duplications, then finally introduce indels
    gt_early_scnas[scnas==-1] <- NA
    dups <- copy(gt_early_scnas)
    dups[scnas!=1] <- NA
    gt_early_scnas <- cbind(gt_early_scnas, dups)
    colnames(gt_early_scnas) <- colnames(indels)
    gt_early_scnas <- gt_early_scnas + indels

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # genotypes for lates SCNAs: indels first, then SCNAs
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nc <- ncol(normal)
    gt_late_scnas <- matrix(nrow=n_clones,ncol=nc)
    for(i in 1:nrow(gt_late_scnas)) gt_late_scnas[i,] <- as.integer(first_cancer_cell)
    rownames(gt_late_scnas) <- clones
    gt_late_scnas <- rbind(gt_late_scnas, normal)
    rownames(gt_late_scnas)[nrow(gt_late_scnas)] <- 'normal'
    gt_late_scnas <- gt_late_scnas + indels[,1:(n_markers*2)]
    gt_late_scnas[scnas==-1] <- NA

    ## add copy number deletions to gt, then add in the duplications, then finally introduce indels
    dups <- copy(gt_late_scnas)
    dups[scnas!=1] <- NA
    gt_late_scnas <- cbind(gt_late_scnas, dups)
    colnames(gt_late_scnas) <- colnames(gt_early_scnas)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # randomly pick some markers from the early SCNAs and others from late SCNAs
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## split the markers into those that had early SCNAs and others that had late SCNAs
    all_markers <- colnames(gt_early_scnas)
    early_markers <- seq(1:n_markers)[rbinom(size=1,n=n_markers,p=prop_markers_with_early_cnas)==1]
    early_markers <- c(paste0('m',early_markers,'.1'), paste0('m',early_markers,'.2'), paste0('m',early_markers,'.3'), paste0('m',early_markers,'.4'))
    late_markers <- all_markers[!all_markers %in% early_markers]

    gt_with_cnvs <- cbind(gt_early_scnas[,early_markers], gt_late_scnas[,late_markers])
    gt_with_cnvs <- gt_with_cnvs[,all_markers]
    gt_with_cnvs
}


##' plot_genotypes_with_cnas
plot_genotypes_with_cnas <- function(gt) {
    ## plot the true (unknowable) genotypes

    x <- as.data.table(gt)
    x$clone <- rownames(gt)
    x <- melt(x, id.vars='clone')
    f=function(s) {
        s <- strsplit(s,'[.]')[[1]]
        list(marker=s[1],copy=as.integer(s[2]))
    }
    markers <- rbindlist(lapply(as.character(x$variable), f))
    x <- cbind(x, markers)
    x$copy <- paste('copy',x$copy)
    x$copy <- factor(x$copy, levels=paste0('copy ',1:4))
    x$clone <- factor(x$clone, levels=rev(rownames(gt)))
    x$marker <- factor(x$marker, levels=unique(x$marker))
    p <- ggplot(x, aes(x=marker, y=clone)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        geom_tile(data=x[is.na(value),], aes(fill=value)) + 
        geom_tile(data=x[!is.na(value),],aes(fill=value),color='black',size=0.25) + 
        scale_fill_gradient2(low='blue',mid='white',high='red',
                             midpoint=round(median(x$value,na.rm=T)),name='Poly-G length (bp)') +
        facet_wrap(facets=~copy, ncol=1) +
        polyGsim_theme(base_size=10) +
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    p
}


##' get_observed_data
get_observed_data <- function(gt, purities=NULL, sd.technical.error=0) {
    ## create 'admixed' genotype data, which is the underlying tumor marker lengths scaled by their purities. 
    ## This can be treated as the admixed difference in mean-marker-lengths from the normal.

    if(!is.null(purities)) {
        tumors <- rownames(gt)
        n_markers <- ncol(gt)
        normal <- as.numeric(gt['normal',])
        for(sa in tumors) gt[sa,] <- ( gt[sa,] * purities[sa,1] ) + (normal * (1 - purities[sa,1]))
    }

    marker_min_meanlengths <- apply(gt, 2, min)
    for(mi in 1:n_markers) {
        gt[,mi] <- gt[,mi] - marker_min_meanlengths[mi]
    }
   
    ## add noise after admixing the mean lengths according to the sample purities
    if(sd.technical.error > 0) for(i in 1:nrow(gt)) gt[i,] <- gt[i,] + rnorm(mean=0,sd=sd.technical.error,n=ncol(gt))
    gt
}


