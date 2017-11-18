#' Get all alteration binary matrics
#' 
#' @param gt gene centric binary table
#' @param alterations A vector of alterations
#' @param alterations.names Name vector for alterations

get_all_alterations <- function(gt, alterations, alteration.names ) {
  if(length(alterations) != length(alteration.names)) {
    ERROR = sprintf("The types of alterations [ %i ] is not equal to the alterations names [ %i ]", 
                    length(alterations), length(alteration.names) )
    stop(ERROR)
  }
  alteration.list = list()
  for( idx in seq.int(length(alterations))) {
    alter =  alterations[idx]
    obj.name = alteration.names[idx]
    INFO = sprintf("[INFO] Getting %s to %s ...\n", alter, obj.name )
    cat(INFO)
    alteration.list[[obj.name]] = get_alteration(gt, alter)
    
  }
  
  return(alteration.list)
  
}

#' Get a binary matrix of a specific type of alteration from gene-centric table
#'
#' @param gt gene-centric table
#' @param alteration one of alteration types
#' @param fill.NA whether to fill NA as 0
#' 
#' @examples
#' aseMat <- get_alteration(gene_centric_table, "ase_all")
get_alteration <- function(gt, 
                           alteration = c("ase_all","cn", "fusion", "isSplice", "variants", "rna_edit", "expr_outlier", "alt_prom", "alt_pa"),
                           fill.NA = TRUE) {
  require(dplyr)
  require(tidyr)  
  alteration =  match.arg(alteration)
  #if( !any(names(gt) %in% alteration)) {
  #  stop("There is no such alteration", alteration, " in the gene centric table")
  #}
  #alteration = "ase_all" 
  mt = gt %>% 
    select_("hgnc_symbol", "ICGC_DONOR_ID", alteration)  %>% 
    spread_("ICGC_DONOR_ID",  alteration ) %>% 
    as.data.frame()
  rownames(mt) = mt$hgnc_symbol
  mt = as.matrix(mt[, 2:ncol(mt)])
  if (fill.NA) {
    mt[is.na(mt)] = 0
  }
  return(mt)
}

#' Calculate P values for raw hits per gene using resampling matrix
#' 
#' @param raw.hits A named vector indicating hits for each gene
#' @param resample.mat  A p*n matrix , p is gene , n is replicate 

calculate_sampling_pvalue <- function(raw.hits, resample.mat) {
  #-- Reorder make sure
  resample.mat = resample.mat[names(raw.hits),]
  
  .ecdf_fun <- function(x,perc) ecdf(x)(perc)   
  
  perc = sapply( seq.int(length(raw.hits)), function(idx) {
    cat(idx, ",")
    .ecdf_fun(resampleMat[idx,], raw.hits[idx])
  })
  cat("\n")
  p.val = 1- perc
  names(p.val) = names(raw.hits)
  return(p.val)
} 


#' Calculate group size for donors

#' @param N The number of total sampling  
#' @param group  A factor vector indicating different group
#' @param is.balance banlance group size by group or not
#' 
calc_group_size <- function(N, group, is.balance = T) {
  
  if (is.balance) {
    #-----------------------------------------------
    # N  |  CancerType 1 (60)  | CancerType 2 (40) |
    # ----------------------------------------------
    # 10 |      5              |          5        |
    #-----------------------------------------------
    # 20 |      10             |          10       |
    #------------------------------------------------
      total.group = length(unique(group))
      group.size = rep(N %/% total.group, total.group)
      remainder = N %% total.group
      if( remainder != 0) {
        idx = sample(seq.int(total.group), size = remainder)
        group.size[idx] = group.size[idx] + 1
      }
      names(group.size) = unique(group)
  } else {
      #-----------------------------------------------
      # N  |  CancerType 1 (60)  | CancerType 2 (40) |
      # ----------------------------------------------
      # 10 |       6             |          4        |
      #-----------------------------------------------
      # 20 |      12             |          8       |
      #------------------------------------------------
      group.size = unlist(as.list(table(group)))
  }
  
  return(group.size)
    
}


#' @param group.table data.frame column 1 is donor, column 2 is group
#' @param N sampling size 
sampling_by_group <-  function(group.table, N = 100, is.balance = T) {
  #--- calculate group size
  group.size = calc_group_size(N = N, group = group.table[[2]] , is.balance = is.balance)
  
  donors.list.by.group  = split(group.table[[1]], group.table[[2]])
  
  rdonors = Map(function(x, n) {sample(x, n, replace = T)},  donors.list.by.group, group.size)
  return(rdonors)
}




#' Sampling  hits matrix from the orignial hits matrix 

#' @param hits.mat orignial hits matrix
#' @param n the number of donors to sample
#' @param n.repeat the number of replicates to sample
#' @param group.table a data.frame indicating which group a donor belong to 
#' @param balance whether to blance group size difference  
#profvis( {
sampling_hits_matrix <- function(hits.mat, n, n.repeat = 1000, group.table = NULL, seed.use = 319, n.cores = 4, 
                                 is.balance = T) {
  #--- Initialing Matrix
  set.seed(seed.use)
  #res.mat = matrix(0, nrow = nrow(hits.mat), ncol = n.repeat)
  require(Matrix)
  donors = colnames(hits.mat)
  
  #--- Initial Parallel
  require(doSNOW)
  cl = makeCluster(n.cores)
  registerDoSNOW(cl)
  res.mat <- foreach (seq.int(n.repeat), .combine = cBind, .export = ls(envir = globalenv())) %dopar%  {
      
    if(!is.null(group.table)) {
        rdonors = sampling_by_group(group.table = group.table, N = n, is.balance = is.balance  )  
        rdonors = unlist(rdonors)
    } else {
        rdonors = sample(donors, size = n, replace = F)
    }
      Matrix::rowSums(hits.mat[,rdonors])
  }
  stopCluster(cl)
  rownames(res.mat) <- rownames(hits.mat)
  colnames(res.mat) <- paste0("Replicate", seq.int(ncol(res.mat)))
  return(res.mat)
}
#system.time(aa <- sampling_hits_matrix(hits.mat = ase, n = 100, n.repeat = 100, group.table = donor2cancer.use, n.cores = 2, do.parallel = T))
#})  
#cBind()
  


#' convert FPKM to TPM. Borrowed from https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

#' @param fpkm FPKM values for a sample , including total genes
fpkm_to_tpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
