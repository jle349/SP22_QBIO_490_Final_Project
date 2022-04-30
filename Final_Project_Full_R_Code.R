library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(survival)
library(survminer)
library(maftools)

# Download Clinical Data
query_clin <- GDCquery(project = "TCGA-COAD", 
                       data.category = "Clinical",
                       file.type = "xml")
GDCdownload(query_clin)
clinic <- GDCprepare_clinic(query_clin, clinical.info = "patient")
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"
colnames(clinic)[colnames(clinic) == "gender"] <- "sex"
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

#replace NA days to death with days to last follow up
clinic$days_to_death = ifelse(is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death)

clinic$death_event = ifelse(clinic$vital_status == "Dead", 1, 0)

#produce survival plot comparing female vs male patients
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

sex_fit <- surv_fit( surv_object ~ clinic$sex, data = clinic)

survplot = ggsurvplot(sex_fit, 
                      pval=TRUE, xlab = "Time (Days)",
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

#rna gene counts analysis
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "STAR - Counts"
)

GDCdownload(query, method="api")
sum_exp = GDCprepare(query)


#maf - mutation analysis
query_maf <- GDCquery(project = "TCGA-COAD",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      legacy = F)
GDCdownload(query_maf)

maf_prep <- GDCprepare(query_maf)
maf_object <- read.maf(maf = maf_prep,
                       clinicalData = clinic,
                       isTCGA = TRUE)


# cooncoplot
clinic <- maf_object@clinical.data

name <- ifelse(clinic$gender=="MALE", TRUE, FALSE)
male_patient_ids <- clinic[name, Tumor_Sample_Barcode]
male_maf <- subsetMaf(maf = maf_object,
                      tsb = male_patient_ids)

name2 <- ifelse(clinic$gender=="FEMALE", TRUE, FALSE)
female_patient_ids <- clinic[name2, Tumor_Sample_Barcode]
female_maf = subsetMaf(maf = maf_object,
                       tsb = female_patient_ids)


coOncoplot(m1 = male_maf, 
           m2 = female_maf, 
           genes = c("APC", "TP53", "TTN", "KRAS", "PIK3CA"),
           m1Name = "Male Patients", 
           m2Name = "Female Patients")

#boxplots
#APC
geneA_id_mask <- rowData(sum_exp)$gene_name == "APC"
apc_counts = assays(sum_exp)$"unstranded"[geneA_id_mask, ]
boxplot(apc_counts ~ colData(sum_exp)$gender,
        main = "APC Gene Counts Comparison of Both Sexes",
        col = "green",
        xlab = "Sex",
        ylab = "Gene Counts"
)

#TP53
geneB_id_mask <- rowData(sum_exp)$gene_name == "TP53"
tp53_counts = assays(sum_exp)$"unstranded"[geneB_id_mask, ]
boxplot(tp53_counts ~ colData(sum_exp)$gender,
        main = "TP53 Gene Counts Comparison of Both Sexes",
        col = "green",
        xlab = "Sex",
        ylab = "Gene Counts"
)

#TTN
geneC_id_mask <- rowData(sum_exp)$gene_name == "TTN"
ttn_counts = assays(sum_exp)$"unstranded"[geneC_id_mask, ]
boxplot(ttn_counts ~ colData(sum_exp)$gender,
        main = "TTN Gene Counts Comparison of Both Sexes",
        col = "green",
        xlab = "Sex",
        ylab = "Gene Counts"
)

#KRAS
geneD_id_mask <- rowData(sum_exp)$gene_name == "KRAS"
kras_counts = assays(sum_exp)$"unstranded"[geneD_id_mask, ]
boxplot(kras_counts ~ colData(sum_exp)$gender,
        main = "KRAS Gene Counts Comparison of Both Sexes",
        col = "green",
        xlab = "Sex",
        ylab = "Gene Counts"
)

#PIK3CA
geneE_id_mask <- rowData(sum_exp)$gene_name == "PIK3CA"
pik3ca_counts = assays(sum_exp)$"unstranded"[geneE_id_mask, ]
boxplot(pik3ca_counts ~ colData(sum_exp)$gender,
        main = "PIK3CA Gene Counts Comparison of Both Sexes",
        col = "green",
        xlab = "Sex",
        ylab = "Gene Counts"
)


#get counts data
counts = assays(sum_exp)$unstranded

gender_na_mask = is.na(colData(sum_exp)$gender)
counts_no_na = counts[,!gender_na_mask] #remove NA from counts
colData_no_na = sum_exp@colData[!gender_na_mask,] #remove NAs from gender column in colData

female_mask = colData_no_na$gender == "female"
male_mask = colData_no_na$gender == "male"

female_counts = counts_no_na[,female_mask] #counts for female patients
male_counts = counts_no_na[,male_mask] #counts for male patients

#temp1 = colnames(counts[,!gender_na_mask])
#temp2 = rownames(sum_exp@colData[!gender_na_mask,])

#prune lowly expressed genes
counts_row_sums = rowSums(counts_no_na)
low_counts_mask = ifelse(counts_row_sums >= 10, TRUE, FALSE)
sum(low_counts_mask)

counts_no_na = counts_no_na[low_counts_mask,]

#DESeq analysis
dds = DESeqDataSetFromMatrix(countData = counts_no_na,
                             colData = colData_no_na,
                             design = ~gender)
dds_obj = DESeq(dds)
resultsNames(dds_obj)
results = results(dds_obj, format = "DataFrame", contrast = c("gender", "male", "female"))

row_order = order(results$padj)
results = results[row_order,]
results = results[!is.na(results$padj),]

#change gene ensembl ID to gene names
get_gene_names <- function(gene_id, sum_exp_row_data) {
  
  name_mask = (sum_exp_row_data$gene_id == gene_id)
  name = sum_exp_row_data$gene_name[name_mask]
  
  return(name)
}
gene_names = lapply(rownames(results), get_gene_names, rowData(sum_exp))
rownames(results) = gene_names

female_order = order(results$log2FoldChange)
female_results = results[female_order,] #top 5 most expressed genes in female patients

male_order = order(results$log2FoldChange, decreasing = TRUE)
male_results = results[male_order,] #top 5 most expressed genes in male patients

#show DESeq results for top 5 most mutated genes
key_genes = c('APC', 'TP53', 'TTN', 'KRAS', 'PIK3CA')
key_genes_results = results[key_genes,]

#KDM61, KDM5C, GENE has log fold change 0.49 higher in females
#IL2 has 0.25 log fold change higher in males


log2FoldChange_threshold = 1
padj_threshold = 0.05

genes_pass_log_threshold = ((results$log2FoldChange > log2FoldChange_threshold) | (results$log2FoldChange < (-1 * log2FoldChange_threshold)))

genes_pass_padj_threshold = (results$padj < padj_threshold) & !is.na(results$padj)

results_subset = results[genes_pass_log_threshold & genes_pass_padj_threshold,]

fc_threshold = 1  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in males",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in males", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Male/Female)",
       y = "-log10 Adjusted p-value")

volcano_plot

log_p = -log10(p_threshold)

#volcano plot of all male vs female genes
EnhancedVolcano(results, lab = rownames(results), x = 'log2FoldChange', 
                y = 'padj', selectLab = c('ZFY', 'AC010086.3','EIF1AY', 'TMSB4Y','TXLNGY',
                                          'XIST', 'APOA4', 'MAGEB2', 'Z84478.1', 'MS4A10',
                                          'APC', 'TP53', 'TTN', 'KRAS', 'PIK3CA'), 
                title = "Male vs. Female CRC Patients", xlab = bquote(~Log[2]~ 'fold change'), 
                pCutoff = log_p, legendLabels = c("NS", "Significant Log2FC", "Significant p-value", "Significant p-value and log2FC"),
                FCcutoff = 1.0, pointSize = 1.5, labSize = 3, labCol = 'black', 
                labFace = 'bold', boxedLabels = TRUE, colAlpha = 4/5, legendPosition = 'right', 
                legendLabSize = 9, legendIconSize = 3.0, drawConnectors = TRUE, hline = log_p, 
                widthConnectors = 0.50, colConnectors = 'black', maxoverlapsConnectors = 30)


tumor_mask = colData(sum_exp)$sample_type == "Primary Tumor"
normal_mask = colData(sum_exp)$sample_type == "Solid Tissue Normal"

female_mask = colData_no_na$gender == "female"
male_mask = colData_no_na$gender == "male"

female_counts = counts_no_na[,female_mask] #counts for female patients
male_counts = counts_no_na[,male_mask] #counts for male patients

female_colData = colData_no_na[female_mask,]
male_colData = colData_no_na[male_mask,]

#DESeq analysis comparing normal tissue to cancer tissue for male and female patients

#Analyis for Female Patients
female_dds = DESeqDataSetFromMatrix(countData = female_counts,
                                    colData = female_colData,
                                    design = ~sample_type)

female_dds_obj = DESeq(female_dds)
resultsNames(female_dds_obj)
#comparing female normal to cancer tissue
female_comp_results = results(female_dds_obj, format = "DataFrame", contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))

female_row_order = order(female_comp_results$log2FoldChange)

#overexpressed in normal tissue
female_comp_results = female_comp_results[female_row_order,]
female_comp_results = female_comp_results[!is.na(female_comp_results$log2FoldChange),]

#overexpressed in cancer tissue
female_row_order2 = order(female_comp_results$log2FoldChange, decreasing = TRUE)
female_comp_results2 = female_comp_results[female_row_order2,]

#rename genes
female_gene_names = lapply(rownames(female_comp_results), get_gene_names, rowData(sum_exp))
rownames(female_comp_results) = female_gene_names

#comparing male normal to cancer tissue
male_dds = DESeqDataSetFromMatrix(countData = male_counts,
                                  colData = male_colData,
                                  design = ~sample_type)

male_dds_obj = DESeq(male_dds)
resultsNames(male_dds_obj)
male_comp_results = results(male_dds_obj, format = "DataFrame", contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))

#rename genes
male_gene_names = lapply(rownames(male_comp_results), get_gene_names, rowData(sum_exp))
rownames(male_comp_results) = male_gene_names

#overexpressed in normal tissue
male_row_order = order(male_comp_results$log2FoldChange)
male_comp_results = male_comp_results[male_row_order,]
male_comp_results = male_comp_results[!is.na(male_comp_results$log2FoldChange),]

#overexpressed in cancer tissue
male_row_order2 = order(male_comp_results$log2FoldChange, decreasing = TRUE)
male_comp_results2 = male_comp_results[male_row_order2,]

#top 5 most mutated genes
key_genes_female_results = female_comp_results2[key_genes,]
key_genes_male_results = male_comp_results2[key_genes,]

EnhancedVolcano(female_comp_results, lab = rownames(female_comp_results), x = 'log2FoldChange', 
                y = 'padj', selectLab = c(), 
                title = "Male vs. Female CRC Patients", xlab = bquote(~Log[2]~ 'fold change'), 
                pCutoff = -log10(0.05), legendLabels = c("NS", "Significant Log2FC", "Significant p-value", "Significant p-value and log2FC"),
                FCcutoff = 1.0, pointSize = 1.5, labSize = 1.5, labCol = 'black', 
                labFace = 'bold', boxedLabels = TRUE, colAlpha = 4/5, legendPosition = 'right', 
                legendLabSize = 9, legendIconSize = 3.0, drawConnectors = TRUE, 
                widthConnectors = 0.50, colConnectors = 'black', maxoverlapsConnectors = 30)

#genes that are overexpresed in female cancer tissue
female_overexpressed_mask = (female_comp_results$log2FoldChange >= 1)
female_overexpressed = female_comp_results[female_overexpressed_mask,]

#genes that are overexpresed in male cancer tissue
male_overexpressed_mask = (male_comp_results$log2FoldChange >= 1)
male_overexpressed = male_comp_results[male_overexpressed_mask,]

male_overexpressed_genes = rownames(male_overexpressed)
female_overexpressed_genes = rownames(female_overexpressed)

#get genes that are only overexpressed in males and those that are only overexpressed in females
overexpressed_genes_in_male_only_mask = !(male_overexpressed_genes %in% female_overexpressed_genes)
overexpressed_genes_in_male_only = male_overexpressed[overexpressed_genes_in_male_only_mask,]
overexpressed_genes_in_female_only_mask = !(female_overexpressed_genes %in% male_overexpressed_genes)
overexpressed_genes_in_female_only = female_overexpressed[overexpressed_genes_in_female_only_mask,]

#list of genes that are only overexpressed in males vs females
#first 5 are only those in male, second row are only those in female
key_genes2 = c("SLC30A2", "RNA5SP225", "BPIFB1", "CASC8", "RNA5-8SP5", 
               "NDUFB4P1", "AL121759.1", "MTCYBP3", "AC092666.1", "VSIG1")

#show comparison of these genes
male_comp_results[key_genes2,]
female_comp_results[key_genes2,]

#SLC30A2, CASC8, VSIG1, BPIFB1

#get genes that are underexpressed in cancer tissue vs normal tissue
female_underexpressed_mask = (female_comp_results$log2FoldChange <= -1)
female_underexpressed = female_comp_results[female_underexpressed_mask,]

male_underexpressed_mask = (male_comp_results$log2FoldChange <= -1)
male_underexpressed = male_comp_results[male_underexpressed_mask,]

male_underexpressed_genes = rownames(male_underexpressed)
female_underexpressed_genes = rownames(female_underexpressed)

underexpressed_genes_in_male_only_mask = !(male_underexpressed_genes %in% female_underexpressed_genes)
underexpressed_genes_in_male_only = male_underexpressed[underexpressed_genes_in_male_only_mask,]
underexpressed_genes_in_female_only_mask = !(female_underexpressed_genes %in% male_underexpressed_genes)
underexpressed_genes_in_female_only = female_underexpressed[underexpressed_genes_in_female_only_mask,]

#key genes: BTNL2, CYP2A8, PCDH11Y,  THBS4**, OR51E2**

#list of genes that are only underexpressed in males and females
#first 5 are only those in male, second row are only those in female
key_genes3 = c("PCDH11Y", "EVX2", "CPB1", "THBS4", "OR51E2",
               "TREH", "XPNPEP2", "FADS6", "BTNL2", "CYP2C8")

male_comp_results[key_genes3,]
female_comp_results[key_genes3,]

#volcano plot of normal vs cancer tissue in male patients
EnhancedVolcano(male_comp_results, lab = rownames(male_comp_results), x = 'log2FoldChange', 
                y = 'padj', selectLab = c("SLC30A2", "RNA5SP225", "BPIFB1", "CASC8", "RNA5-8SP5", 
                                          "PCDH11Y", "EVX2", "CPB1", "THBS4", "OR51E2"), 
                title = "Cancer Tissue vs. Normal Tissue in Male CRC Patients", xlab = bquote(~Log[2]~ 'fold change'), 
                pCutoff = log_p, legendLabels = c("NS", "Significant Log2FC", "Significant p-value", "Significant p-value and log2FC"),
                FCcutoff = 1.0, pointSize = 1.5, labSize = 4, labCol = 'black', 
                labFace = 'bold', boxedLabels = TRUE, colAlpha = 4/5, legendPosition = 'right', 
                legendLabSize = 9, legendIconSize = 3.0, drawConnectors = TRUE, hline = log_p, 
                widthConnectors = 0.50, colConnectors = 'black', maxoverlapsConnectors = 30)

#volcano plot of normal vs cancer tissue in female patients
EnhancedVolcano(female_comp_results, lab = rownames(female_comp_results), x = 'log2FoldChange', 
                y = 'padj', selectLab = c("NDUFB4P1", "AL121759.1", "MTCYBP3", "AC092666.1", "VSIG1", 
                                          "TREH", "XPNPEP2", "FADS6", "BTNL2", "CYP2C8"), 
                title = "Cancer Tissue vs. Normal Tissue in Female CRC Patients", xlab = bquote(~Log[2]~ 'fold change'), 
                pCutoff = log_p, legendLabels = c("NS", "Significant Log2FC", "Significant p-value", "Significant p-value and log2FC"),
                FCcutoff = 1.0, pointSize = 1.5, labSize = 5, labCol = 'black', 
                labFace = 'bold', boxedLabels = TRUE, colAlpha = 4/5, legendPosition = 'right', 
                legendLabSize = 9, legendIconSize = 3.0, drawConnectors = TRUE, hline = log_p, 
                widthConnectors = 0.50, colConnectors = 'black', maxoverlapsConnectors = 30)











