# 查询和校正植物学名
taxon_info <- Taxonstand::TPL(c(
  "Salix futura Seemen",
  "Salix hukaoana Kimura"
))
write.csv(taxon_info, "taxon_info.csv")