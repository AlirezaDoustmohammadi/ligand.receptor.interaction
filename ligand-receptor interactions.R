library(CellChat)
library(patchwork)
library(ggplot2)
library(igraph)


read_files <- function(input_dir, GEO.number) {
  
  #' Read normalized count matrix and metadata files for a GEO number
  #' This function loads the normalized count matrix and metadata files for a given GEO number. 
  #' It also reorders the rows of the metadata file to match the order of the cells in the count matrix.
  #' @param input_dir the input directory where the files are located
  #' @param GEO.number the GEO number for which to load the files
  
  # Load the normalized count matrix and metadata
  # genes should be in rows with rownames and cells in columns with colnames
  normalized.csv.file <- paste0(input_dir, GEO.number,'/', GEO.number, 
                                "_Normalized_Gene_count_table.csv")
  print(normalized.csv.file)
  count.matrix <- as.matrix(read.csv(normalized.csv.file, row.names=1, header = TRUE))
  
  # a matrix (rows are cells with rownames) consisting of cell information,
  # which will be used for defining cell groups.
  metadata.file <- paste0(input_dir, GEO.number,'/', GEO.number, 
                          "_metadata.csv")
  print(metadata.file)
  metadata <- as.matrix(read.csv(metadata.file, row.names=1, header = TRUE))
  
  
  # Reorder the rows of the count.matrix & metadata 
  # to match the order of the cells in the count.matrix
  colnames_input <- colnames(count.matrix)
  rownames_metadata <- rownames(metadata)
  
  # Find the indices of the metadata rows that match the column names of the count matrix
  match_indices <- match(colnames_input, rownames_metadata)
  
  # Reorder the rows of the metadata file using the match indices
  metadata_reordered <- metadata[match_indices, ]
  
  return(list(count.matrix, metadata_reordered)) 
}



create_cellchat_object <- function(count.matrix, metadata_reordered, selected.group) {
  
  #' This function creates a CellChat object for performing 
  #' cell-cell communication analysis on the given count matrix and metadata, 
  #' and sets the ligand-receptor interaction database for the analysis.
  
  # Create a CellChat object
  cellchat <- createCellChat(object = count.matrix, 
                             meta = metadata_reordered, group.by = selected.group)
  
  # Set the ligand-receptor interaction database
  # Load the mouse CellChat database
  CellChatDB <- CellChatDB.mouse
  # use a subset of CellChatDB for cell-cell communication analysis 
  # This step is necessary even if using the whole database
  # "use Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact" (all interaction annotation types)
  # "Secreted Signaling": secrete autocrine/paracrine signaling interactions (59.9 %)
  # "ECM-Receptor": extracellular matrix (ECM)-receptor interactions (21.4 %)
  CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact" )) 
  
  # use all CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  # subset the expression data of signaling genes for saving computation cost
  # This step is necessary even if using the whole count matrix
  cellchat <- subsetData(cellchat) 
  
  return(cellchat)
}


cell_cell_communication <- function(cellchat, workers) {
  
  #' This function performs cell-cell communication analysis. 
  #' The function performs various steps including identification of over-expressed signaling genes, 
  #' over-expressed ligand-receptor interactions, computation of communication probabilities, 
  #' and filtering of communication events.
  
  # set configuration for parallel computing
  future::plan("multiprocess", workers = workers) 
  
  # To infer the cell state-specific communications, over-expressed ligands or receptors in one cell group 
  # has been identified and then  over-expressed ligand-receptor interactions has been identified 
  # if either ligand or receptor is over-expressed.
  # Identify over-expressed signaling genes associated with each cell group
  # Threshold of p-values = 0.05
  cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.05)
  # Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # CellChat infers the biologically significant cell-cell communication by assigning 
  # each interaction with a probability value and peforming a permutation test.
  # CellChat models the probability of cell-cell communication by integrating gene expression 
  # as well as prior known knowledge of the interactions between 
  # signaling ligands, receptors and their cofactors using the law of mass action.
  # triMean is used for calculating the average gene expression per cell group.
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.1)
  
  # # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # CellChat computes the communication probability on signaling pathway level by 
  # summarizing the communication probabilities of all ligands-receptors interactions 
  # associated with each signaling pathway.
  cellchat <- computeCommunProbPathway(cellchat)
  # calculation of the aggregated cell-cell communication network 
  # by counting the number of links or summarizing the communication probability.
  cellchat <- aggregateNet(cellchat)  
  
  return(cellchat)
}



aggregated_cell_cell_communication_network <- function(cellchat, output_dir, GEO.number, selected_group) {
  
  # Set up PDF file
  pdf(paste0(output_dir, GEO.number, "/", GEO.number, "_", selected_group, 
             ".aggregated cell-cell communication network.pdf"), 
      width = 11.69, height = 8.27, paper = "a4r")
  
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                   weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
                   weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  # Close PDF file
  dev.off()  

  }


pathways_in_detail <- function(cellchat, output_dir, GEO.number, selected_group) {
  
  # Visualization of cell-cell communication network
  # detected pathway according to count matrix
  pathways.show = cellchat@netP[["pathways"]]
  
  # Create the new directory
  dir.create(paste0(output_dir, GEO.number, "/", selected_group))
  dir.create(paste0(output_dir, GEO.number, "/", selected_group,  "/signalling pathway network"))
  dir.create(paste0(output_dir, GEO.number, "/", selected_group,  "/signalling pathway network/circle_plot"))
  dir.create(paste0(output_dir, GEO.number, "/", selected_group,  "/signalling pathway network/chord_plot"))
  
  # Loop through pathways.show and generate plot for each element
  for (i in seq_along(pathways.show)) {
    pathway <- pathways.show[[i]]
    
    # "circle_plot"
    # Set up PDF file
    pdf(paste0(output_dir, GEO.number, "/", selected_group,  "/signalling pathway network/circle_plot/", 
               pathway, " signalling pathway network.pdf"), 
        width = 11.69, height = 8.27, paper = "a4r")
    
    # Generate plot
    netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
    
    # Close PDF file
    dev.off()
    
    
    "chord_plot"
    # Set up PDF file
    pdf(paste0(output_dir, GEO.number, "/", selected_group,  "/signalling pathway network/chord_plot/", 
               pathway, " signalling pathway network.pdf"), 
        width = 11.69, height = 8.27, paper = "a4r")
    
    # Generate plot
    netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
    
    # Close PDF file
    dev.off()
    
  }
  
}



contribution_l_R <- function(cellchat, output_dir, GEO.number, selected_group) {
  
  # Visualization of cell-cell communication network
  # detected pathway according to count matrix
  pathways.show = cellchat@netP[["pathways"]]
  
  # Create the new directory
  dir.create(paste0(output_dir, GEO.number, "/", selected_group,  "/contribution of each ligand-receptor"))
  
  # Loop through pathways.show and generate plot for each element
  for (i in seq_along(pathways.show)) {
    pathway <- pathways.show[[i]]
    
    # Generate plot
    # Compute the contribution of each ligand-receptor pair to the overall signaling pathway
    plt <- netAnalysis_contribution(cellchat, signaling = pathway)
    
    ggsave(paste0(output_dir, GEO.number, "/", selected_group,  "/contribution of each ligand-receptor/", 
                  pathway, " contribution of each ligand-receptor.pdf"), plot = plt)
  }
  
  return(cellchat)
}


functional_similarity <- function(cellchat, output_dir, GEO.number, selected_group) {
  # Identify signaling groups based on their functional similarity
  
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netEmbedding(cellchat, type = "functional")
  #> Manifold learning of the signaling networks for a single dataset
  cellchat <- netClustering(cellchat, type = "functional")
  # Classification learning of the signaling networks for a single dataset
  
  # Visualization in 2D-space
  plt <- netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  
  ggsave(paste0(output_dir, GEO.number, "/", selected_group,  "/functional similarity.pdf"), plot = plt)
  
  return(cellchat)
}



structural_similarity <- function(cellchat, output_dir, GEO.number, selected_group) {
  # Identify signaling groups based on their structural similarity
  
  # Identify signaling groups based on their structural similarity
  
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural")
  #> Manifold learning of the signaling networks for a single dataset
  cellchat <- netClustering(cellchat, type = "structural")
  # Classification learning of the signaling networks for a single dataset
  
  # Visualization in 2D-space
  plt <- netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  
  ggsave(paste0(output_dir, GEO.number, "/", selected_group,  "/structural similarity.pdf"), plot = plt)
  
  return(cellchat)
}

Centrality_analysis <- function(cellchat, output_dir, GEO.number, selected_group) {
  
  # Visualization of cell-cell communication network
  # detected pathway according to count matrix
  pathways.show = cellchat@netP[["pathways"]]
  
  # Compute the network centrality scores
  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  # Create the new directory
  dir.create(paste0(output_dir, GEO.number, "/", selected_group,  "/Centrality analysis"))
  
  # Loop through pathways.show and generate plot for each element
  for (i in seq_along(pathways.show)) {
    pathway <- pathways.show[[i]]
    
    # Generate plot
    # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
    pdf(paste0(output_dir, GEO.number, "/", selected_group,  "/centrality analysis/", 
               pathway, " centrality analysis.pdf"), 
        width = 11.69, height = 8.27, paper = "a4r")
    netAnalysis_signalingRole_network(cellchat, signaling = pathway, width = 8, height = 2.5, font.size = 10)
    
    dev.off()
  }
  
  return(cellchat)
  
}


save_interactions <- function(cellchat, output_dir, GEO.number, selected_group) {
  # founded intercations
  interactions = cellchat@LR[["LRsig"]]
  write.csv(interactions, paste0(output_dir, GEO.number, "/", selected_group,
                                 "/", GEO.number,  "_interactions.csv"))

}

# GEO numbers vector
GEO.vector.number = c('GSE109125')

# input directory
input_dir <- "input/"
# output directory to save results
output_dir <- "output/"

group.vector <- c("Cell_type", "TISSUE")


for (i in 1:length(GEO.vector.number)) {
  GEO.number <- GEO.vector.number[i]
  for (j in 1: length(group.vector)) {
    selected.group <- group.vector[j]
    
    # read files
    readed.file <- read_files(input_dir, GEO.number)
    
    count.matrix <- readed.file[[1]]
    metadata_reordered <- readed.file[[2]]
    
    # create cellchat object
    cellchat <- create_cellchat_object(count.matrix, metadata_reordered, selected.group)
    # finding strong communication
    cellchat <- cell_cell_communication(cellchat, workers = 6)
    
    # aggregated cell cell communication network
    aggregated_cell_cell_communication_network(cellchat, output_dir, GEO.number, selected.group)
    
    
    # pathways detail
    pathways_in_detail(cellchat, output_dir, GEO.number, selected.group)
    # ligand receptor contribution
    contribution_l_R(cellchat, output_dir, GEO.number, selected.group)
    
    # Identify signaling groups based on their functional similarity
    cellchat <- functional_similarity(cellchat, output_dir, GEO.number, selected.group)
    # Identify signaling groups based on their structural similarity
    cellchat <- structural_similarity(cellchat, output_dir, GEO.number, selected.group)
    
    # centrality analysis
    cellchat <- Centrality_analysis(cellchat, output_dir, GEO.number, selected.group)
    
    # save intercations
    save_interactions(cellchat, output_dir, GEO.number, selected.group)
    
    # Save the CellChat object
    saveRDS(cellchat, file = paste0(output_dir, GEO.number, "/", selected.group,
                                    "/", GEO.number,  "_cellchat_LS.rds"))
    
  }
}
