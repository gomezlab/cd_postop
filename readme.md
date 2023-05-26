JAK/STAT Signaling and Oxidative Phosphorylation Underlying Post-Operative Infectious Complications and Crohn’s Disease Recurrence: A Transcriptomic Analysis

This code performs differential expression analysis for Crohn's Disease recurrence and post-operative colorectal surgical site infection. There are two sections in /source for each of the two databases and outcomes. All code is in R except for the lr.ipynb notebook.

deseq_newruv_..._R.ipynb - performs the DESeq analysis using RUVg factors, also generate heatmaps using normalized/transformed counts
gsea_R.ipynb - performs gene set enrichment analysis using the Hallmark gene sets
ssgsea_R.ipynb - performs single sample gene set enrichment using the Hallmark gene sets
volcano_R.ipynb - creates a volcano plot from the DESeq output

