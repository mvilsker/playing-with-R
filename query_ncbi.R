#!/usr/bin/env Rscript
library(modules)
library(readr)

query_ncbi <- module({
  import(rentrez)
  import(jsonlite)

  rna_sequencing_handler <- function(acc) {
    #A handler for rna sequencing datasets that includes srp identifier, sra link and srr.
    
    gse_search <- entrez_search(db="gds", term=paste(acc," AND gse[entry type]"), use_history = TRUE)
    gse_sum <- entrez_summary(db="gds", web_history = gse_search$web_history)
    srp <- gse_sum$extrelations$targetobject
    
    srr_search_res <- entrez_search(db="sra", term=srp, use_history = TRUE)
    accs <- entrez_fetch(db="sra", web_history = srr_search_res$web_history, rettype="acclist")
    
    ex <- extract_from_esummary(gse_sum,c("extrelations"))
    accs <- gsub("\n\n","\n",accs)
    ex$srr <- toJSON(stringr::str_split(accs, '\n'), auto_unbox = TRUE)

    return(ex)
  }
  
  experiment_summary <- function(acc) {
      #experiment_summary: accession id, gpls, suppfile and ftplink. Use GSE59847
      search_res <- entrez_search(db="gds", term=paste(acc," AND gse[entry type]"))
      summ_result <- entrez_summary(db="gds", id=search_res$ids)
      
      df <- extract_from_esummary(summ_result,c("accession", "gpl", "suppfile", "ftplink"))
      return(df)
  }
})
