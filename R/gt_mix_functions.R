#' @include utilities.R

#' This is a function to read in the experimental design
#' @title read_experimental_design
#' @param experimental_design_path path to the experimental design .csv file
#' @return a matrix of the experimental design
#' @examples
#' \donttest{
#' read_experimental_design()
#' }
#' @export
read_experimental_design <- function(experimental_design_path){
  message("reading experimental design")
  edesign <- read.csv(experimental_design_path, row.names = 1, header = TRUE)
  return(edesign)
}

#' This is a function to plot the experimental design
#' @param experimental_design an experimental design matrix rownames should be microfluidics channels, colnames should be genotypes
#' @import reshape2 ggplot2
#' @return a ggplot experimental design heatmap
#' @examples
#' \donttest{
#' #' plot_experimental_design()
#' }
#' @export
plot_experimental_design <- function(experimental_design){
  requireNamespace("reshape2")
  channels <- factor(rownames(experimental_design))
  design_df <- reshape2::melt(data.frame(experimental_design, "channel" = rownames(experimental_design)))
  design_df$channel <- factor(design_df$channel, levels = channels)

  pl <-  ggplot(design_df,
                aes(y= channel, x = variable, fill = factor(value))) +
    geom_tile(color = 'grey80')+
    scale_fill_manual(values = c('1' = 'black',
                                 '0' = 'white')) +
    theme_classic() + coord_fixed() +
    theme(legend.position = "none", axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size =7.5, color = 'black', face = 'bold'),
          axis.text = element_text(color = 'black', size = 5, face = 'plain'),
          axis.text.x = element_text(angle = 45, hjust=1)) +
    xlab("Donor") + ylab("10X channel")
  return(pl)
}

#' This is a function to get the number of shared genotypes between channels
#' @param experimental_design design matrix rownames should be microfluidics channels, colnames should be genotypes
#' @return a matrix of channels and genotypes shared
#' @examples
#' \donttest{
#' # get_shared()
#' }
#' @export
get_shared <- function(experimental_design){
  requireNamespace('reshape2')
  mat <- matrix(NA, ncol = nrow(experimental_design), nrow = nrow(experimental_design))
  colnames(mat) <- rownames(mat) <- rownames(experimental_design)
  design <- lapply(rownames(experimental_design), function(x){colnames(experimental_design)[experimental_design[x,] >0]})
  names(design) <- colnames(mat)
  for(i in colnames(mat)){
    for(j in rownames(mat)){
      mat[i, j] <- length(intersect(design[[i]], design[[j]]))
    }
  }
  mat[lower.tri(mat)] = 0 # we can omit the lower tri because distances will be symmetric
  # make the lower tri of this longform
  mat <- reshape2::melt(mat)
  mat <- mat[!mat$value == 0, ]
  mat <- mat[!mat$Var1 == mat$Var2, ]
  colnames(mat) <- c("channel_1", "channel_2", "shared")
  rownames(mat) <- 1:nrow(mat)
return(mat)
}

#' This is a function to read in the locations of Soup Or Cell outputs
#' @param SOC_locations path to the .csv file giving channel names as the first column and file paths to the souporcell output in the second column
#' @return a matrix of SOC file paths and the channel names
#' @examples
#' \donttest{
#' #' read_SOC_locations()
#' }
#' @export
read_SOC_locations <- function(SOC_locations){
  message("reading souporcell file locations")
  locations <- read.csv(SOC_locations, header = TRUE)
  return(locations)
}

#' this is a function that compares genotype calls across two experiments
#' @param experiment_1_path the file path to the first souporcell directory
#' @param experiment_2_path the file path to the second souporcell directory
#' @param experiment_1_name the name of experiment 1
#' @param experiment_2_name the name of experiment 2
#' @param shared numeric- the number of shared genotypes between the experiments
#' @return a data frame of the shared genotype clusters between the two experiments
#' @import vcfR
#' @examples
#' \donttest{
#' #' shared_genotypes()
#' }
#' @export
shared_genotypes <- function(experiment_1_path, experiment_2_path, shared, experiment_1_name, experiment_2_name ){
  library(vcfR, quietly = TRUE)
  vcf1_path <- file.path(experiment_1_path, "cluster_genotypes.vcf")
  vcf2_path <- file.path(experiment_2_path, "cluster_genotypes.vcf")
  vcf1 =  read.vcfR(file.path(experiment_1_path, "cluster_genotypes.vcf"),
                          verbose = FALSE)
  vcf2 =  read.vcfR(file.path(experiment_2_path, "cluster_genotypes.vcf"),
                          verbose = FALSE)
  #get intersection of fix - here we just use
  vcf1_idx <- paste0(vcf1@fix[, 1], "_", vcf1@fix[,2 ], "_", vcf1@fix[,4], "_", vcf1@fix[, 5])
  vcf2_idx <- paste0(vcf2@fix[, 1], "_", vcf2@fix[,2 ], "_", vcf2@fix[,4], "_", vcf2@fix[, 5])
  intersect_idx <-  intersect(vcf1_idx, vcf2_idx)
  #subset vcfs to intersect of fix
  vcf1 <- vcf1[match(intersect_idx, vcf1_idx), ]
  vcf2 <- vcf2[match(intersect_idx, vcf2_idx), ]
  #now get the genotype calls
  gt_1  <- vcf1@gt
  AO_idx_gt_1 <- which(strsplit(gt_1[,1], ":")[[1]] == "AO")
  RO_idx_gt_1 <-  which(strsplit(gt_1[,1], ":")[[1]] == "RO")
  gt_1 <- lapply(2:ncol(gt_1), function(x){
     gt1_mat <-  do.call(rbind, strsplit(gt_1[,x], ":"))[, c(AO_idx_gt_1, RO_idx_gt_1)]
     gt1_mat <- apply(gt1_mat, 2, as.numeric)
     return(gt1_mat)
  })
  names(gt_1) <- 0:(length(gt_1)-1)
  gt_2 <- vcf2@gt
  AO_idx_gt_2 <- which(strsplit(gt_2[,1], ":")[[1]] == "AO")
  RO_idx_gt_2 <-  which(strsplit(gt_2[,1], ":")[[1]] == "RO")
  gt_2 <- lapply(2:ncol(gt_2), function(x){
    gt2_mat <-  do.call(rbind, strsplit(gt_2[,x], ":"))[, c(AO_idx_gt_2, RO_idx_gt_2)]
    gt2_mat <- apply(gt2_mat, 2, as.numeric)
    return(gt2_mat)
  })
  names(gt_2) <- 0:(length(gt_2)-1)

  #subset to variants that have sufficient counts
  rs_1 <- lapply(gt_1, rowSums)
  rs_2 <- lapply(gt_2, rowSums)
  sufficient_counts <- apply(do.call(cbind, c(rs_1, rs_2)), 2, function(x){x > 10}) #set sufficient counts as 10
  sufficient_counts <- apply(sufficient_counts, 1, all)
  #do this subsetting and calculate the vaf
  gt_1 <- lapply(gt_1, function(x){x <- x[sufficient_counts, ]
  x <- x[,2]/(x[,1] + x[,2])
  return(x)})
  gt_2 <- lapply(gt_2, function(x){x <- x[sufficient_counts, ]
  x <- x[,2]/(x[,1] + x[,2])
  return(x)})

  #and do a double for loop to calculate the mean squared error
  df_list <- list()
  p = 1
  for(i in names(gt_1)){
    for(j in names(gt_2)){
      i_geno <- gt_1[[i]]
      j_geno <- gt_2[[j]]
      mean_squared_error <- mean((i_geno-j_geno)^2)
      df <- data.frame("experiment_1" = paste0(experiment_1_name, "_", i), "experiment_2" = paste0(experiment_2_name, "_", j), "mse" = mean_squared_error)
      df_list[[p]] <- df
      p = p+1
    }
  }
  df <- do.call(rbind, df_list)
  df <- df[order(df$mse, decreasing= F), ]
  df <- df[1:shared, ]
  return(df)
}

#' this is a big function that constructs a genotype cluster graph using a set of SOC directories and an experimental design
#' @param experimental_design an experimental design matrix rownames should be microfluidics channels, colnames should be genotypes
#' @param file_locations the file locations
#' @return a list $graph_membership gives you cluster memberships $graph_plot gives a force directed embedding of the graph $membership_plot gives a heatmap of memberships $membership_matrix gives a matrix of channel memberships
#' @import igraph ggraph utils ggplot2
#' @examples
#' \donttest{
#' construct_genotype_cluster_graph()
#' }
#' @export
construct_genotype_cluster_graph <- function(experimental_design, file_locations){
  if(all(file_locations$channel == rownames(experimental_design))){
    message("checking files")
  }else{stop("! the channel column of the locations file does not match the rownames of the experimental design matrix")}
  contrasts <- get_shared(experimental_design)
    #now make a graph matrix!
  n_genotypes <- data.frame("n_genotypes" = rowSums(experimental_design)) #this is the number of clusters in each channel
  ch_use <- rownames(experimental_design)
  mat_names <- unlist(lapply(ch_use, function(x){
    paste0(x, "_", c(0:(n_genotypes[x, 1]-1)))
  })) #these are the names of the clusters in the channels
 #make a graph matrix -  we will fill this in
   graph_matrix <- matrix(0, ncol = length(mat_names), nrow = length(mat_names))
  colnames(graph_matrix) <- rownames(graph_matrix) <- mat_names
  pb <- txtProgressBar(min = 0, max = nrow(contrasts), style = 3)
for(x in rownames(contrasts)){
  #get which experiments we want to compare and get their paths
  selected_contrast <- as.list(contrasts[x, ])
  exp1 = as.character(selected_contrast[[1]])
  exp1_path <- file_locations[file_locations$channel %in% exp1, 2]
  exp2 = as.character(selected_contrast[[2]])
  exp2_path <- file_locations[file_locations$channel %in% exp2, 2]
  #and the number of shared genotypes
  shared_gt = as.integer(selected_contrast[[3]])
  out <- shared_genotypes(experiment_1_path = exp1_path, experiment_2_path = exp2_path, shared = shared_gt,
                          experiment_1_name = exp1, experiment_2_name = exp2)
      #then fill in the values into our graph matrix
    for(p in seq_len(nrow(out))){
      idx_1 <- out[p, 1]
      idx_2 <- out[p, 2]
      graph_matrix[idx_1, idx_2] = 1
      graph_matrix[idx_2, idx_1] = 1
    }
  #set progress bar
  setTxtProgressBar(pb, as.numeric(x))
}
  close(pb)
  #now make a graph using igraph
  requireNamespace('igraph')
  gr <- igraph::graph_from_adjacency_matrix(graph_matrix, mode = "undirected", weighted = NULL)
  #cluster the graph
  graph_membership <- igraph::cluster_walktrap(gr)$membership
  names(graph_membership) <- igraph::V(gr)$name
  cluster_colors <- cluster_color_ramp()(length(unique(graph_membership)))
  names(cluster_colors) <- unique(graph_membership)
  requireNamespace('ggraph')
  set.seed(0)
  membership_graph_plot <- ggraph::ggraph(gr, layout = 'fr') + ggraph::geom_edge_link() +
    #geom_node_point(pch = 21, aes(fill = factor(clusters)), size = 5) +
    theme_void() +
    scale_fill_manual(values = cluster_colors) +
    theme(legend.position = 'none') + ggraph::geom_node_label(aes(label = graph_membership, fill = factor(graph_membership)))
#now get the membership of each of these clusters per channel
  ident_mat <- matrix(0, ncol = length(ch_use), nrow = length(unique(graph_membership)))
  colnames(ident_mat) <- ch_use
  rownames(ident_mat) <- unique(graph_membership)
  for(i in rownames(ident_mat)){
    for(j in colnames(ident_mat) ){
      cl_samples <- colnames(graph_matrix)[graph_membership %in% i]
      ident_mat[i, j] = any(grep(j, cl_samples))*1
    }
  }

  #plot membership by cluster
  requireNamespace("reshape2")
  ident_mat_df <- reshape2::melt(ident_mat)
  membership_plot <- ggplot(ident_mat_df, aes(x= factor(Var1), y = factor(Var2, levels = ch_use), fill = factor(value))) +
    geom_tile(color = 'grey80') + scale_fill_manual(values = c('1' = 'black',
                                                                             '0' = 'white')) +
    theme_classic() + coord_fixed() +
    theme(legend.position = "none", axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size =7.5, color = 'black', face = 'bold'),
          axis.text = element_text(color = 'black', size = 5, face = 'plain')) +
    xlab("Computed genotypes") + ylab("10X channel")

  #now we test whether the graph is a collection  of complete subgraphs

subgraph_checks <- unlist(lapply(unique(graph_membership), function(m){
  requireNamespace('igraph')
    subgraph <- igraph::induced_subgraph(gr, vids = igraph::V(gr)$name %in% names(graph_membership)[graph_membership %in% m] )
    check_complete_graph(subgraph)
  }))

if(all(subgraph_checks)){
  message("the membership graph is a collection of complete subgraphs as expected")
}else{
  warning("the membership graph is NOT a collection of complete subgraphs - donor genotype matching is probably not possible - please double check the experimental design ")
}
  #then return a bunch of stuff
  return(list("genotype_cluster_graph" = gr,
              "graph_membership" = graph_membership,
              "graph_plot" = membership_graph_plot,
              "membership_plot" = membership_plot,
              "membership_matrix" = ident_mat))
}
#' this is a big function that checks if the genotype graph is a set of complete subgraphs
#' @param gr an igraph object
#' @return a message telling you if your graph is a collection of complete subgraphs
#' @import igraph
#' @examples
#' \donttest{
#' check_complete_graph()
#' }
#' @export
check_complete_graph <- function(gr){
  requireNamespace('igraph')
  expected_edges <- igraph::vcount(gr)*(igraph::vcount(gr)-1)/2
  measured_edges <- igraph::ecount(gr)
  return(expected_edges == measured_edges)
}
#' this is a function that maps channel membership to genotype ID
#' @param experimental_design design matrix rownames should be microfluidics channels, colnames should be genotypes
#' @param graph_output the output of construct_genotype_cluster_graph given as x$genotype_cluster_graph
#' @return a matrix of channel memberships for each cluster
#' @import reshape2
#' @examples
#' \donttest{
#' membership_map()
#' }
#' @export
membership_map <- function(experimental_design, graph_output){
  requireNamespace("reshape2")
  graph_membership <- graph_output$graph_membership
  memb_mat <- graph_output$membership_matrix
  ident_mat_df <- reshape2::melt(memb_mat)
  cluster_mapping <- list()
  n=1
  for(i in rownames(memb_mat)){
    for(j in colnames(experimental_design)){
      if(all(experimental_design[,j] == memb_mat[i,])){
        chs <- colnames(memb_mat)[memb_mat[i, ] ==1]
        cluster_mapping[[n]] <- data.frame("channel" = chs, "SOC_cluster" = names(graph_membership[graph_membership == i]), "genotype_cluster" = i, "genotype" = j)
        n=n+1
      }
    }
  }
  cluster_mapping <- do.call(rbind, cluster_mapping)
  #throw a warning message if we can't find all the genotypes
  if(!all(colnames(experimental_design) %in% unique(cluster_mapping$genotype))){
    warning("not all genotypes can be recovered - go back and check the experimental design and the number of clusters called")
  }
return(cluster_mapping)
}
#' this is a function that maps the genotype clusters to genotype names
#' @export
#' @param SOC_locations the locations file
#' @param membership_mat output of membership_map()
#' @return matrix of cells by genotypes
#' @examples
#' \donttest{
#' cells_to_genotypes()
#' }
#' @export
cells_to_genotypes <- function(SOC_locations, membership_mat){
aggregated_clusters  <-  do.call(rbind, lapply(1:nrow(SOC_locations), function(x){
channel <- SOC_locations[x, 1]
clusters <- read.csv(file.path(SOC_locations[x, 2], "clusters.tsv"), sep = "\t")
clusters <- clusters[, 1:3]
clusters$channel <- channel
if(all(clusters$assignment == 0)){
  clusters$genotype <- membership_mat[membership_mat$channel %in% clusters$channel, "genotype"]
  clusters$SOC_cluster <- membership_mat[membership_mat$channel %in% clusters$channel, "SOC_cluster"]
}
else{
SOC_cluster <- strsplit(clusters$assignment, "/")
genotype_mapping <- do.call(rbind, lapply(SOC_cluster, function(soc){
  if(length(soc) == 1){
    p <- paste0(channel, "_", soc)
    gt <- membership_mat[membership_mat$SOC_cluster %in% p, "genotype"]
    df <- data.frame("genotype" = gt, "SOC_cluster" = p)
  }else{
    p1 <- paste0(channel, "_", soc[1])
    p2 <- paste0(channel, "_", soc[2])
    gt <- paste0(membership_mat[membership_mat$SOC_cluster %in% p1, "genotype"], "/", membership_mat[membership_mat$SOC_cluster %in% p2, "genotype"])
    p <- paste0(p1, "/",p2)
    df <- data.frame("genotype" = gt, "SOC_cluster" = p)
  }
  return(df)
}))
clusters <- cbind(clusters, genotype_mapping)
}
return(clusters)
}))
return(aggregated_clusters)}

#' this is a function that maps the genotype clusters to genotype names
#' @param experiment_1_path the SOC directory for the first experiment
#' @param experiment_2_path the SOC directory for the second experiment
#' @param experiment_1_name the name of the first experiment
#' @param experiment_2_name the name of the second experiment
#' @import ggplot2 vcfR
#' @return matrix of cells by genotypes
#' @examples
#' \donttest{
#' cells_to_genotypes()
#' }
#' @export
plot_cross_vaf <- function(experiment_1_path, experiment_2_path, experiment_1_name, experiment_2_name){
  requireNamespace("vcfR")
  vcf1_path <- file.path(experiment_1_path, "cluster_genotypes.vcf")
  vcf2_path <- file.path(experiment_2_path, "cluster_genotypes.vcf")
  vcf1 =  vcfR::read.vcfR(file.path(experiment_1_path, "cluster_genotypes.vcf"),
                          verbose = FALSE)
  vcf2 =  vcfR::read.vcfR(file.path(experiment_2_path, "cluster_genotypes.vcf"),
                          verbose = FALSE)
  #get intersection of fix - here we just use
  vcf1_idx <- paste0(vcf1@fix[, 1], "_", vcf1@fix[,2 ], "_", vcf1@fix[,4], "_", vcf1@fix[, 5])
  vcf2_idx <- paste0(vcf2@fix[, 1], "_", vcf2@fix[,2 ], "_", vcf2@fix[,4], "_", vcf2@fix[, 5])
  intersect_idx <-  intersect(vcf1_idx, vcf2_idx)
  #subset vcfs to intersect of fix
  vcf1 <- vcf1[match(intersect_idx, vcf1_idx), ]
  vcf2 <- vcf2[match(intersect_idx, vcf2_idx), ]
  #now get the genotype calls
  gt_1  <- vcf1@gt
  AO_idx_gt_1 <- which(strsplit(gt_1[,1], ":")[[1]] == "AO")
  RO_idx_gt_1 <-  which(strsplit(gt_1[,1], ":")[[1]] == "RO")
  gt_1 <- lapply(2:ncol(gt_1), function(x){
    gt1_mat <-  do.call(rbind, strsplit(gt_1[,x], ":"))[, c(AO_idx_gt_1, RO_idx_gt_1)]
    gt1_mat <- apply(gt1_mat, 2, as.numeric)
    return(gt1_mat)
  })
  names(gt_1) <- 0:(length(gt_1)-1)
  gt_2 <- vcf2@gt
  AO_idx_gt_2 <- which(strsplit(gt_2[,1], ":")[[1]] == "AO")
  RO_idx_gt_2 <-  which(strsplit(gt_2[,1], ":")[[1]] == "RO")
  gt_2 <- lapply(2:ncol(gt_2), function(x){
    gt2_mat <-  do.call(rbind, strsplit(gt_2[,x], ":"))[, c(AO_idx_gt_2, RO_idx_gt_2)]
    gt2_mat <- apply(gt2_mat, 2, as.numeric)
    return(gt2_mat)
  })
  names(gt_2) <- 0:(length(gt_2)-1)

  #subset to variants that have sufficient counts
  rs_1 <- lapply(gt_1, rowSums)
  rs_2 <- lapply(gt_2, rowSums)
  sufficient_counts <- apply(do.call(cbind, c(rs_1, rs_2)), 2, function(x){x > 10}) #set sufficient counts as 10
  sufficient_counts <- apply(sufficient_counts, 1, all)
  #do this subsetting and calculate the vaf
  gt_1 <- do.call(rbind, lapply(names(gt_1), function(i){
    x <- gt_1[[i]]
    x <- x[sufficient_counts, ]
    x <- x[,2]/(x[,1] + x[,2])
    df <- data.frame("VAF" = x, "genotype" = paste0(experiment_1_name, "_", i), "channel" = experiment_1_name)
    return(df)}))
  gt_2 <- do.call(rbind, lapply(names(gt_2), function(i){
    x <- gt_2[[i]]
    x <- x[sufficient_counts, ]
    x <- x[,2]/(x[,1] + x[,2])
    df <- data.frame("VAF" = x, "genotype" = paste0(experiment_2_name, "_", i), "channel" = experiment_2_name)
    return(df)}))

  #mung this into a big df
  df_list <- list()
  p=1
  for(i in unique(gt_1$genotype)){
    for(j in unique(gt_2$genotype)){
      df_list[[p]] <- data.frame(
        "experiment_1_genotype" = gt_1[gt_1$genotype %in% i, "VAF"],
        "experiment_2_genotype" = gt_2[gt_2$genotype %in% j, "VAF"],
        "channel_1" = i,
        "channel_2" = j)
      p=p+1
    }
  }
  df <- do.call(rbind, df_list)
  pl <- ggplot(df, aes(x = experiment_1_genotype, y = experiment_2_genotype)) +
    geom_point(pch = 19, cex = 0.05) +
    facet_grid(rows = vars(channel_2), cols = vars(channel_1)) + coord_fixed() + ylim(c(0, 1)) + xlim(c(0, 1)) + theme_minimal() + theme(axis.text = element_blank())

  return(pl)
}

#' this is a function produces an experimental design with varying degrees of sparsity.
#' @param n_mixtures the number of mixtures to include
#' @param n_genotypes the number of genotypes to include
#' @param density the density of the design - numeric value from 0->1
#' @return matrix of mixtures by genotypes
#' @examples
#' \donttest{
#' make_overlapping_mixture()
#' }
#' @export
make_overlapping_mixture <- function (n_mixtures, n_genotypes, density = 1) {
  max_genotypes = 2^n_mixtures - 1
  if (n_genotypes > max_genotypes) {
    message("too many genotypes")
  }
  else {
    exp_design <- do.call(cbind, lapply(1:n_mixtures, function(m) {
      mixture_names <- paste0("mixtures_", 1:n_mixtures)
      comb <- combn(mixture_names, m = m, simplify = TRUE)
      p = 1
      combo_list <- list()
      for (i in 1:ncol(comb)) {
        vec <- mixture_names %in% comb[, i]
        names(vec) <- mixture_names
        combo_list[[p]] <- vec
        p = p + 1
      }
      return(do.call(cbind, combo_list) * 1)
    }))
    colnames(exp_design) <- paste0("genotype_", 1:max_genotypes)
    windows <- 1:(ncol(exp_design) - n_genotypes)
    windows <- lapply(windows, function(x) {
      colnames(exp_design)[c(x:(x + n_genotypes - 1))]
    })
    #at this point check that all the windows incorporate all the channels - this can be a problem if the n mixtures is high and n genotypes is low and the density is low
    window_rs <- lapply(windows, function(x){
      all(rowSums(exp_design[, x]) > 0)
    } )
    windows <- windows[unlist(window_rs)]
    window_density <- 1:length(windows)/length(windows)
    if (n_genotypes == max_genotypes) {
      window_select <- windows[[1]]
    }
    else {
      window_select <- windows[[which.min(abs(window_density -
                                                density))]]
    }
    exp_design <- exp_design[, window_select]
    colnames(exp_design) <- paste0("genotype_", 1:n_genotypes)
    return(exp_design)
  }
}

