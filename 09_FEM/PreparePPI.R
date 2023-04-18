# Create PPT network.

library(STRINGdb)
library(igraph)
library(biomaRt)

message("getSTRINGdb for Zebrafish")
string_db <- STRINGdb$new(version="11", species=7955)
zebrafish_graph <- string_db$get_graph()

message("get edges with high confidence score")
edge.scores <- E(zebrafish_graph)$combined_score
ninetyth.percentile <- quantile(edge.scores, 0.9)
thresh <- data.frame(name='90th percentile',
                     val=ninetyth.percentile)
zebrafish_graph <- subgraph.edges(zebrafish_graph,
                              E(zebrafish_graph)[combined_score > ninetyth.percentile])

message("create adjacency matrix")
adj_matrix <- as_adjacency_matrix(zebrafish_graph)


# 4. map gene ids to protein ids

message("get gene/protein ids via Biomart")
mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='drerio_gene_ensembl')

message("extract protein ids from the zebrafish network")
protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'), function(x) x[2])

message("get protein to gene id mappings")
mart_results <- getBM(attributes = c("ensembl_gene_id",
                                     "ensembl_peptide_id"),
                      filters = "ensembl_peptide_id", values = protein_ids,
                      mart = mart)

message("replace protein ids with gene ids")
ix <- match(protein_ids, mart_results$ensembl_peptide_id)
ix <- ix[!is.na(ix)]

newnames <- protein_ids
newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
    mart_results[ix, 'ensembl_gene_id']
rownames(adj_matrix) <- newnames
colnames(adj_matrix) <- newnames

ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
nullrows <- Matrix::rowSums(ppi)==0
ppi <- ppi[!nullrows,!nullrows] ## ppi is the network with gene ids

message("Save Generated PPI information into rda file in Data Folder")
if (!file.exists("./Data")) dir.create("./Data")
save(ppi, file="./Data/ppi.rda")

