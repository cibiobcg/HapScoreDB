library(tidyverse)

# Load variant rsids
df <- read_tsv(file = "../data/interface_data.tsv") %>% select(rsid)
variants <- do.call(c, list(df$rsid))
variants <- strsplit(variants, ",")
variants <- unlist(variants) 
variants <- variants %>% unique()
variants <- variants[variants!="wt"]
variants <- data.frame("rsid" = variants)

# Load variants present in clinvar
clinvar <- read_delim("../data/clinvar_variants.txt", delim = "\n", col_names = FALSE)
colnames(clinvar) <- "rsid"
clinvar <- clinvar %>% mutate(rsid = paste0("rs", rsid))
# Only keep variants in hapscoredb
clinvar <- clinvar %>% dplyr::filter(rsid %in% variants$rsid)

# Load variants present in gwas catalog
gwas_cat <- read_tsv("../data/gwas_catalog_v1.0-associations_e115_r2025-09-15.tsv") %>% select(SNPS)
# Only keep variants in hapscoredb
gwas_cat <- gwas_cat %>% filter(SNPS %in% variants$rsid)

# Create variants links to clinvar
variants$link <- lapply(variants$rsid, function(rsid){
        if (rsid %in% clinvar$rsid) {
                return(sprintf('<a href="https://www.ncbi.nlm.nih.gov/clinvar/?term=%s" target="_blank">%s</a>(ClinVar)', rsid, rsid))
        }else if (rsid %in% gwas_cat$SNPS) {
                return(sprintf('<a href="https://www.ebi.ac.uk/gwas/variants/%s" target="_blank">%s</a>(GWAS Catalog)', rsid, rsid))
        }else{
                return(sprintf('<a href="https://www.ncbi.nlm.nih.gov/snp/%s" target="_blank">%s</a>(dbSNP)', rsid, rsid))
        }
})

# save link pan
saveRDS(variants, file = "../data/link_map.rds")
