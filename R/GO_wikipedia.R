library(dplyr)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(fgsea)
options(stringsAsFactors = F)

###Load necessary files####
wikipedia_GO_gmt <- fread("data/wikidata_GO_terms.tsv")

wikipedia_GO_df <- wikipedia_GO_gmt %>%
  dplyr::filter(!duplicated(.[,-7]))

GO_gmt <- clusterProfiler::read.gmt("data/c5.go.bp.v7.5.1.symbols.gmt")
GO_gmt$term <- tolower(gsub(gsub(as.character(GO_gmt$term),pattern = "GOBP_",replacement = ""),
                            pattern = "_",replacement = " "))

wikipedia_GO_fgsea <- split(wikipedia_GO_df$gene_symbol,
                            wikipedia_GO_df$itemLabel)

GO_fgsea <- split(GO_gmt$gene,
                  GO_gmt$term)

nterms_GO <- length(unique(GO_gmt$term))
nterms_wikipedia <- length(unique(wikipedia_GO_df$itemLabel))

#Terms GO = 7658
terms_GO <- unique(GO_gmt$term) %>%
  gsub(pattern = "GOBP_",replacement = "") %>%
  gsub(pattern = "_",replacement = " ") %>%
  tolower()

#Terms Wikidata = 269
terms_wikidata <- unique(wikipedia_GO_df$itemLabel)

#ACE2 overexpressing cells infected with MOI=2.0 SARS-CoV-2 for 24h
COVID_DEGs <- fread("data/DEGs_series16.csv")

COVID_DEGs$log2FC_FDR <- COVID_DEGs$log2FC * (-log(COVID_DEGs$FDR))

ranks <- COVID_DEGs$log2FC
names(ranks) <- COVID_DEGs$external_gene_name

#run fgsea
set.seed(1231)
fgsea_wikipedia <- fgsea::fgseaSimple(pathways = wikipedia_GO_fgsea,stats = ranks,
                                      nperm = 1000)

fgsea_GO <- fgsea::fgseaSimple(pathways = GO_fgsea,stats = ranks,
                               nperm = 1000)

#Filter enrichment by padj < 0.1
fgsea_wikipedia_sig <- fgsea_wikipedia %>%
  filter(padj < 0.1) %>%
  mutate(db="wikipedia")

p <- fgsea_wikipedia_sig %>%
  mutate(dir=ifelse(NES>0,"up","down")) %>%
  group_by(dir) %>%
  top_n(n = 6,wt = abs(NES)) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(pathway=factor(pathway,levels = unique(pathway))) %>%
  ggplot(aes(x=NES,y = pathway,fill=NES))+
  geom_col()+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red3")+
  theme_minimal()+
  ggtitle("Wikipedia GO terms")
p

fgsea_GO_sig <- fgsea_GO %>%
  filter(padj < 0.1) %>%
  mutate(db="GO")

p <- fgsea_GO_sig %>%
  mutate(dir=ifelse(NES>0,"up","down")) %>%
  group_by(dir) %>%
  top_n(n = 6,wt = abs(NES)) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(pathway=factor(pathway,levels = unique(pathway))) %>%
  ggplot(aes(x=NES,y = pathway,fill=NES))+
  geom_col()+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red3")+
  theme_minimal()+
  ggtitle("All GO terms")
p

fgsea_all <- rbind(fgsea_wikipedia_sig,fgsea_GO_sig)

p <- fgsea_all %>%
  mutate(dir=ifelse(NES>0,"up","down")) %>%
  group_by(dir,db) %>%
  top_n(n = 6,wt = abs(NES)) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(pathway=factor(pathway,levels = unique(pathway))) %>%
  ggplot(aes(x=NES,y = pathway,fill=NES))+
  geom_col()+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red3")+
  theme_minimal()+
  facet_wrap(facets = ~db,scales = "free")
pdf(file = "figures/go_wikipedia_fgsea_COVID.pdf",width = 9,height = 3)
print(p)
dev.off()

#Schizophrenia DEGs Gandal 2018b
SCZ_DEGs <- fread("data/gandal_DEGs_clean.csv")
SCZ_DEGs <- SCZ_DEGs[SCZ_DEGs$gene_name!="NOG",]

ranks <- SCZ_DEGs$SCZ.log2FC
names(ranks) <- SCZ_DEGs$gene_name

#run fgsea
set.seed(1231)
fgsea_wikipedia <- fgsea::fgseaSimple(pathways = wikipedia_GO_fgsea,stats = ranks,
                                      nperm = 1000)

fgsea_GO <- fgsea::fgseaSimple(pathways = GO_fgsea,stats = ranks,
                               nperm = 1000)

#Filter enrichment by padj < 0.1
fgsea_wikipedia_sig <- fgsea_wikipedia %>%
  filter(padj < 0.1) %>%
  mutate(db="wikipedia")

fgsea_GO_sig <- fgsea_GO %>%
  filter(padj < 0.1) %>%
  mutate(db="GO")

fgsea_all <- rbind(fgsea_wikipedia_sig,fgsea_GO_sig)

p <- fgsea_all %>%
  mutate(dir=ifelse(NES>0,"up","down")) %>%
  group_by(dir,db) %>%
  top_n(n = 6,wt = abs(NES)) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(pathway=factor(pathway,levels = unique(pathway))) %>%
  ggplot(aes(x=NES,y = pathway,fill=NES))+
  geom_col()+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red3")+
  theme_minimal()+
  facet_wrap(facets = ~db,scales = "free")
pdf(file = "figures/go_wikipedia_fgsea_SCZ.pdf",width = 9,height = 3)
print(p)
dev.off()
