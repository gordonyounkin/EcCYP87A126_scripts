####### Script to plot differential expression between leaves and roots vs collinum z-score ########
####### Gordon Younkin, April 17, 2023 #######
####### Script written for CYP87A126 paper, figure 2A #######

library("edgeR")

##### Load Data #####
# Load tissue gene expression data #
expr <- read.csv("tissue_gene_counts_normalized.matrix", sep="\t")
# Load species gene expression data #
spec <- read.csv("species_gene_counts_normalized_vst_transformed.matrix", sep="\t")
# load gene id to Arabidopsis thaliana function for stringtie IDs #
func <- data.frame(readxl::read_xlsx("ECE_STRG_to_Artha_function.xlsx"))
# load gene id to Artha funciton for original IDs
func_orig <- read.csv("Eche_v2.0_genes_with_annotated_function.txt", sep="\t")
func_orig$gene_id <- gsub("\\..*", "", func_orig$gene_id)

#### Organize and prepare for analysis #####
# make vector with names of young leaves and roots samples
samples <- names(expr)[grepl("YLC", names(expr)) | grepl("RtC", names(expr))]
# genes: get full list of P450s
# first from new stringtie genes based on BLAST to Arabidopsis
genes <- func[grepl("P450", func$At_annotation) & func$geneid %in% row.names(expr) & !grepl("antisense", func$geneid), "geneid"]
# then based on original genome annotation
genes_orig <- func_orig[grepl("P450", func_orig$annotation) & func_orig$gene_id %in% row.names(expr), "gene_id"]
# and combine
genes_all <- unique(c(genes, genes_orig))


########## edgeR for differential expression analysis between young leaves and roots ########################
# get P450 expression in leaves and roots
expr_subset <- expr[genes_all,samples]
# remove samples with 0 expression
expr_subset <- expr_subset[apply(expr_subset, 1, max) > 0, ]
groups <- numeric(length(samples))
groups[grepl("YL", samples)] <- 2
groups[grepl("Rt", samples)] <- 1

y <- DGEList(expr_subset, group=groups)
y <- estimateDisp(y)
et <- exactTest(y)
et_fdr <- topTags(et, nrow(et$table))
#############################################################################################################

########## Calculate z-score of expression in E. collinum relative to expression in all other species ########
# get P450s
spec_subset <- spec[genes_all, ]
# calculate z-score
spec_subset$mean <- rowMeans(spec_subset)
spec_subset$sd <- apply(spec_subset, 1, sd)
spec_subset$col_z <- (spec_subset$PAS_L2 - spec_subset$mean) / spec_subset$sd
# spec_subset$col_z contains colinum z score to be used for plotting 

########## Plot logFC vs colinum z-score ####################################################################
# combine FC and z-score into dataframe
toplot <- merge(et_fdr, spec_subset[, c("mean","sd","col_z")], by="row.names")
# convert ln to log10
toplot$log10fc <- log10(exp(toplot$logFC))
# label points of interest
tolabel <- toplot[toplot$log10fc > 0 & toplot$col_z < -3, ]

# plot
pdf("P450s_LfRt_colinumZ_NOCOLOR.pdf", h=6, w=5)
plot(toplot$col_z ~ toplot$log10fc, pch=16, las=1,
     xlab=expression(paste("log"[10],"fold-change expression in ", italic("E. cheiranthoides"), " leaves relative to roots")),
     ylab=expression(paste("z-score CPM in ", italic("E. colinum"), "relative to other ", italic("Erysimum"), " species")))
# add lines for x and y axes
abline(v=0)
abline(h=0)
# add labels
text(x=c(tolabel$log10fc)-c(0,0.3,0), y=tolabel$col_z, 
     c("CYP716A418", "CYP87A126", "CYP71B132"), pos=3, cex=0.8)
dev.off()
