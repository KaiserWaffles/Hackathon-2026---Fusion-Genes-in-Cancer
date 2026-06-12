#Script for finding the human orthologs for mouse genes and creating TSV files with the ortholog and fusion gene columns

#Load libraries
library(gprofiler2)
library(stringr)

#Create function to separate genes in the same cell
gene_separator <- function(genes) {
  gene_list <- str_extract_all(genes, "[A-Za-z0-9]+(?=\\()")
  result <- unlist(gene_list)
  if(length(result) == 0) {
    result <- trimws(unlist(strsplit(genes, ",")))
  }
  return(result)
}

#Create function to find human orthologs of mouse genes 
#Write orthologs and the fused genes to a new TSV file
ortholog_finder <- function(tsv_file, output_file){
  final_list <- c()
  
  #add empty columns into tsv
  tsv$gene1_ortholog <- character(nrow(tsv))
  tsv$gene2_ortholog <- character(nrow(tsv))
  tsv$fusion_genes <- character(nrow(tsv))
  
  #iterate through rows
  for(row in 1:nrow(tsv)){
    genes1 = gene_separator(tsv$gene1[row])
    genes2 = gene_separator(tsv$gene2[row])
    
    fusion_gene_list = c()
    gene1_ortho_list = c()
    gene2_ortho_list = c()
    
    for(gene1 in genes1) {
      for(gene2 in genes2) {
        
        hs_gene1 = mouse_to_human[gene1]
        hs_gene2 = mouse_to_human[gene2]
        
        if(is.na(hs_gene1)){
          hs_gene1 ='no_ortholog' #label as no_ortholog if none are found
        }
        if(is.na(hs_gene2)){
          hs_gene2 = 'no_ortholog'
        }
        
        gene1_ortho_list = append(gene1_ortho_list, hs_gene1)
        gene2_ortho_list = append(gene2_ortho_list, hs_gene2)
        
        fusion_gene = paste(hs_gene1, '::', hs_gene2)
        fusion_gene_list = append(fusion_gene_list, fusion_gene)
      }
    }
    #Paste genes into the new columns
    tsv$gene1_ortholog[row] <- paste(unique(gene1_ortho_list), collapse = ",")
    tsv$gene2_ortholog[row] <- paste(unique(gene2_ortho_list), collapse = ",")
    tsv$fusion_genes[row] <- paste(unique(fusion_gene_list), collapse = ",")
    
    final_list[[row]] = fusion_gene_list
  }
  
  #Custom order for columns
  original_cols <- names(tsv)[!names(tsv) %in% c("gene1_ortholog", "gene2_ortholog", "fusion_genes")]
  new_order <- c(original_cols[1:2], "gene1_ortholog", "gene2_ortholog", "fusion_genes", original_cols[3:length(original_cols)])
  tsv <- tsv[, new_order]
  
  #Write to new file
  write.table(tsv, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
  
  return(tsv)
}

#Create list of TSV files within the mouse folder
tsv_files <- list.files(path = "PRJNA1062304_arriba_mouse_ifn", pattern = "\\.tsv$", full.names = TRUE)

#Loop through TSV files and apply function
for(file in tsv_files){

  tsv_file=file
  
  base_name = basename(tsv_file)
  dir_name = dirname(tsv_file)
  
  #Ensure that ortholog_results folder is present beforehand
  output_file = file.path(dir_name,'ortholog_results',paste0("ortholog_", base_name))
  
  tsv=read.delim(tsv_file,sep='\t')
  names(tsv)[1]<-'gene1'
  names(tsv)[2]<-'gene2'
  
  # Create list of all mouse genes in TSV file
  all_mouse_genes <- unique(c(unlist(sapply(tsv$gene1, gene_separator)), unlist(sapply(tsv$gene2, gene_separator))))
  print(paste("Found", length(all_mouse_genes), "unique mouse genes in file", base_name))
  
  # Get orthologs for all mouse genes
  orthologs <- gorth(query = all_mouse_genes, source_organism = "mmusculus", target_organism = "hsapiens",filter_na = TRUE)
  
  # Create table for mouse-human orthologs
  mouse_to_human <- setNames(orthologs$ortholog_name, orthologs$input)
  
  ortholog_finder(tsv, output_file)

}
