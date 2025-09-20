library(plumber)
library(dplyr)
library(readr)

data <- read_tsv("../data/interface_data.tsv")

#* @apiTitle HapscoreDB API
#* @apiDescription API interface for the data presented in HapScoreDB

#* @param genes:[character] Ensembl gene ID of genes of interest
#* @param variants:[character] rsid of variants of interest 
#* @param format:[character] Either 'csv' or 'json'
#* @get /api/data

function(res, genes = NULL, variants = NULL, format = "json"){
        
        df <- data
        
        # Filtering
        if (!is.null(genes)) {
                # Handle comma-separated values
                gene_list <- unlist(strsplit(genes, ","))
                df <- df %>% dplyr::filter(gene_id %in% gene_list)
        }
        if (!is.null(variants)) {
                # Handle comma-separated values
                variants_list <- unlist(strsplit(variants, ","))
                variants_pattern <- paste(variants_list, collapse = "|")
                df <- df %>% dplyr::filter(grepl(variants_pattern, rsid))
        }
        
        # Return JSON (default)
        if (tolower(format) == "json") {
                res$setHeader("Content-Type", "application/json")
                return(df)
        }
        
        # Return CSV
        if (tolower(format) == "csv") {
                res$setHeader("Content-Type", "text/csv")
                res$setHeader("Content-Disposition", 'attachment; filename="data.csv"')
                return(paste(capture.output(write.csv(df, row.names = FALSE)), collapse = "\n"))
        }
        
        # Unsupported format
        res$status <- 400
        return(list(error = "Invalid format. Use 'json' or 'csv'."))
}