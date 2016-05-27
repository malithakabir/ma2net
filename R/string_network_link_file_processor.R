


dir("./input")
a1=readLines("./input/9606.protein.links.v10.txt")
# > length(a1)
# [1] 8548003

a2=strsplit(a1, " ")
# > length(a2)
# [1] 8548003

rm("a1")
a3=lapply(a2, function(x){as.numeric(x[3])>700})

a4=a2[unlist(a3)]
# > length(a4)
# [1] 639779

rm("a2")
rm("a3")

a5=do.call(rbind.data.frame, a4)
# > dim(a5)
# [1] 639778      3

rm("a4")

colnames(a5)=c("protein1", "protein2", "combined_score")

# write.csv(a5, "./input/unannotated_fitered_string_dataset.csv")

library(igraph)
g <- graph.data.frame(d = a5, directed = FALSE)
# > length(V(g)$name)
# [1] 15478
# #
# V(g)$name[1]
# V(g)$name[duplicated(V(g)$name)]
# Graph created and no duplicate id found

V(g)$name=gsub("9606.","",V(g)$name)

library(org.Hs.eg.db)
columns(org.Hs.eg.db)

annotation=select(org.Hs.eg.db, V(g)$name, "SYMBOL", keytype="ENSEMBLPROT")

if(nrow(annotation)==length(V(g)$name)){
  g1=g
  V(g1)$name=annotation$SYMBOL
  g2=delete_vertices(g1, which(is.na(V(g1)$name)==TRUE))
} else {
  multiple_gene_symbol_for_same_protein=annotation$ENSEMBLPROT[duplicated(annotation$ENSEMBLPROT)]
  g1=delete_vertices(g, multiple_gene_symbol_for_same_protein)
  # > length(V(g)$name)>length(V(g1)$name)
  # [1] TRUE
  annotation1=select(org.Hs.eg.db, V(g1)$name, "SYMBOL", keytype="ENSEMBLPROT")
  # > nrow(annotation1)==length(V(g1)$name)
  # [1] TRUE
  V(g1)$name=annotation1$SYMBOL
  # > length(V(g1)$name)
  # [1] 15328
  g2=delete_vertices(g1, which(is.na(V(g1)$name)==TRUE))
  # > length(V(g2)$name)
  # [1] 14505
}


# write.csv(get.data.frame(g1), "./data/annotated_filtered_string_dataset.csv")

saveRDS(global_network_data, "./data/global_network_data.Rds")
# rm(list=ls())
# global_network_data <- readRDS("./data/global_network_data.Rds")

