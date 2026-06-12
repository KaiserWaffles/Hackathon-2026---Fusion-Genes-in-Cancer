#Creating fusion genes column within the human data TSV files

#Load library
library(stringr)

#Create function to separate genes 
gene_separator <- function(genes) {
  gene_list <- str_extract_all(genes, "[A-Za-z0-9]+(?=\\()")
  result <- unlist(gene_list)
  if(length(result) == 0) {
    result <- trimws(unlist(strsplit(genes, ",")))
  }
  return(result)
}

# !!! Substitute this with the other antiEGFR filepath to generate resulting TSV files !!!
tsv_files <- list.files(path = "PRJNA1269539_arriba_human_bile_salts", pattern = "\\.tsv$", full.names = TRUE)

#Loop through TSV files within each folder
for(file in tsv_files){
  
  tsv_file=file

  base_name = basename(tsv_file)
  dir_name = dirname(tsv_file)
  
  #Ensure subdirectory named 'fusion_results' is already present before running the code
  output_file = file.path(dir_name,'fusion_results',paste0("fused_genes_", base_name))
  
  tsv=read.delim(tsv_file,sep='\t')
  names(tsv)[1]<-'gene1'
  names(tsv)[2]<-'gene2'
  
  #Iterate through rows of TSV file
  for (row in 1:nrow(tsv)){
    
    genes1 = gene_separator(tsv$gene1[row])
    genes2 = gene_separator(tsv$gene2[row])
    
    fusion_gene_list = c()
    
    #Create fusion gene name and write to new column
    for(gene1 in genes1) {
      for(gene2 in genes2) {
        
        fusion_gene = paste(genes1, '::', genes2)
        fusion_gene_list = append(fusion_gene_list, fusion_gene)
        tsv$fusion_genes[row] <- paste(unique(fusion_gene_list), collapse = ",")
      }
    }
    
    original_cols <- names(tsv)[!names(tsv) %in% c("fusion_genes")]
    
    #create custom order for columns
    new_order <- c(original_cols[1:2], "fusion_genes", original_cols[3:length(original_cols)])
    tsv <- tsv[, new_order]
    
    #Write new TSV
    write.table(tsv, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
  }
}
