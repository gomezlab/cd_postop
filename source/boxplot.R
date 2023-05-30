# %%
## Ben's pipeline for RNA-seq instructions 
library("AnnotationDbi")
library("org.Hs.eg.db")
library(edgeR)
library(RUVSeq)

# This is useful for when you want to run in a different directory
# Should be updated based on where you are running this
source("downstreamAnalysis_RNAseqFunctions.R")
# %%

# %%
#CHANGE THIS GENE LIST TO PLOT
gene_list <- c('HCAR3')
# %%

# %%
# Read in data
cts = as.matrix(read.csv("../data/cts_i_symbol.csv", sep=',', row.names='symbol'))
ruvdata = read.csv('../results/recur/i/data_ruvg.csv', sep=',', row.names='X')

cts = cts[,-which(colnames(cts) == 'PB17_236_257_A5')]
cts = cts[,-which(colnames(cts) == 'PB600298_A4')]
cts = cts[,-which(colnames(cts) == 'PB43300_A1')]

dim(ruvdata)
dim(cts)
# %%


# %%
factor_cols = c('sex', 'recur')
for (col in factor_cols) {
  ruvdata[,col] = as.factor(ruvdata[,col])
}
sex = ruvdata$sex
source = ruvdata$source
# tissue = ruvdata$tissue
recur = ruvdata$recur
inflammation_status = ruvdata$inflammation_status
tin = ruvdata$TIN

sample = rownames(ruvdata)

#round counts to integers
cts <- round(cts);

#check if row names of ruvdata match colnames of cts
all(rownames(ruvdata) == colnames(cts))

# %%

# %%
formula_list = list(
~ sex + W_1 + recur, 
~ sex + W_1 + W_2 + recur, 
~ sex + W_1 + W_2 + W_3 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + W_18 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + W_18 + W_19 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + W_18 + W_19 + W_20 + recur)

#use 18 covariates
design_formula = formula(formula_list[[18]])
design_formula
dds = DESeqDataSetFromMatrix(countData=cts, colData=ruvdata, design = design_formula)
#filter genes for low expression with edgeR
y <- DGEList(counts=counts(dds), group=recur)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds <- dds[keep,]
#run DESeq
dds <- DESeq(dds)
res <- results(dds)

#use RemoveBatchEffect to adjust for the identified covariates
#first create variables for each of the RUV factors
ruv1 = ruvdata$W_1
ruv2 = ruvdata$W_2
ruv3 = ruvdata$W_3
ruv4 = ruvdata$W_4
ruv5 = ruvdata$W_5
ruv6 = ruvdata$W_6
ruv7 = ruvdata$W_7
ruv8 = ruvdata$W_8
ruv9 = ruvdata$W_9
ruv10 = ruvdata$W_10
ruv11 = ruvdata$W_11
ruv12 = ruvdata$W_12
ruv13 = ruvdata$W_13
ruv14 = ruvdata$W_14
ruv15 = ruvdata$W_15
ruv16 = ruvdata$W_16
ruv17 = ruvdata$W_17
ruv18 = ruvdata$W_18

#run RemoveBatchEffect using the RUV factors and sex
vsd = vst(dds)
mat= assay(vsd)
assay(vsd) <- removeBatchEffect(assay(vsd), batch=sex, covariates=cbind(ruv1, ruv2, ruv3, ruv4, ruv5, ruv6, ruv7, ruv8, ruv9, ruv10, ruv11, ruv12, ruv13, ruv14, ruv15, ruv16, ruv17, ruv18))
# %%

# %%
#CHANGE THIS GENE LIST TO THE GENES OF INTEREST
#define a list of genes to plot

#plot the boxplot for each gene and save the results to the results/boxplots folder
for (gene in gene_list){
    d <- plotCounts(dds, gene=gene, intgroup='recur', returnData=TRUE)
    ggplot(d, aes(x=recur, y=count, fill=recur)) +
        geom_boxplot(alpha = 0) +
        geom_point(position=position_jitter(w=0.1,h=0)) +
        scale_y_log10(breaks=c(100, 200, 500, 1000)) +
        labs(x="Subtype", y="Normalized/Transformed Expression", title=gene, fill='Subtype') +
        scale_color_manual(values=c("#00BFC4", "#F8766D")) +
        theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),    
                legend.title=element_blank(), legend.key=element_blank())
    ggsave(paste0("../results/boxplots/recur/", gene, "_i.png"), width=6, height=6, dpi=300)
}

# %%

# %%
# Read in data
cts = as.matrix(read.csv("../data/cts_c_symbol.csv", sep=',', row.names='symbol'))
ruvdata = read.csv('../results/recur/c/data_ruvg.csv', sep=',', row.names='X')

cts = cts[,-which(colnames(cts) == 'PB17_236_257_A1')]
cts = cts[,-which(colnames(cts) == 'PB600052_A7')]
cts = cts[,-which(colnames(cts) == 'PB600298_A8')]

dim(ruvdata)
dim(cts)
# %%


# %%
factor_cols = c('sex', 'recur')
for (col in factor_cols) {
  ruvdata[,col] = as.factor(ruvdata[,col])
}
sex = ruvdata$sex
source = ruvdata$source
# tissue = ruvdata$tissue
recur = ruvdata$recur
inflammation_status = ruvdata$inflammation_status
tin = ruvdata$TIN

sample = rownames(ruvdata)

#round counts to integers
cts <- round(cts);

#check if row names of ruvdata match colnames of cts
all(rownames(ruvdata) == colnames(cts))

# %%

# %%
formula_list = list(
~ sex + W_1 + recur, 
~ sex + W_1 + W_2 + recur, 
~ sex + W_1 + W_2 + W_3 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + recur, 
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + W_18 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + W_18 + W_19 + recur,
~ sex + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10 + W_11 + W_12 + W_13 + W_14 + + W_15 + W_16 + W_17 + W_18 + W_19 + W_20 + recur)

#use 18 covariates
design_formula = formula(formula_list[[18]])
design_formula
dds = DESeqDataSetFromMatrix(countData=cts, colData=ruvdata, design = design_formula)
#filter genes for low expression with edgeR
y <- DGEList(counts=counts(dds), group=recur)
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds <- dds[keep,]
#run DESeq
dds <- DESeq(dds)
res <- results(dds)

#use RemoveBatchEffect to adjust for the identified covariates
#first create variables for each of the RUV factors
ruv1 = ruvdata$W_1
ruv2 = ruvdata$W_2
ruv3 = ruvdata$W_3
ruv4 = ruvdata$W_4
ruv5 = ruvdata$W_5
ruv6 = ruvdata$W_6
ruv7 = ruvdata$W_7
ruv8 = ruvdata$W_8
ruv9 = ruvdata$W_9
ruv10 = ruvdata$W_10
ruv11 = ruvdata$W_11
ruv12 = ruvdata$W_12
ruv13 = ruvdata$W_13
ruv14 = ruvdata$W_14
ruv15 = ruvdata$W_15
ruv16 = ruvdata$W_16
ruv17 = ruvdata$W_17
ruv18 = ruvdata$W_18

#run RemoveBatchEffect using the RUV factors and sex
vsd = vst(dds)
mat= assay(vsd)
assay(vsd) <- removeBatchEffect(assay(vsd), batch=sex, covariates=cbind(ruv1, ruv2, ruv3, ruv4, ruv5, ruv6, ruv7, ruv8, ruv9, ruv10, ruv11, ruv12, ruv13, ruv14, ruv15, ruv16, ruv17, ruv18))
# %%

# %%
#plot the boxplot for each gene and save the results to the results/boxplots folder
for (gene in gene_list){
    d <- plotCounts(dds, gene=gene, intgroup='recur', returnData=TRUE)
    ggplot(d, aes(x=recur, y=count, fill=recur)) +
        geom_boxplot(alpha = 0) +
        geom_point(position=position_jitter(w=0.1,h=0)) +
        scale_y_log10(breaks=c(100, 200, 500, 1000)) +
        labs(x="Subtype", y="Normalized/Transformed Expression", title=gene, fill='Subtype') +
        scale_color_manual(values=c("#00BFC4", "#F8766D")) +
        theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), axis.title=element_text(face="bold", size=20),
                axis.text=element_text(size=14, color="black"),
                panel.background = element_rect(fill="white",color="black"),
                legend.title=element_blank(), legend.key=element_blank())
    ggsave(paste0("../results/boxplots/recur/", gene, "_c.png"), width=6, height=6, dpi=300)
}

# %%