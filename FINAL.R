#####Install Libraries#####
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
BiocManager::install("chromVAR")
BiocManager::install("MotifDb")

#####Load Libraries#####
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(dorothea))
suppressPackageStartupMessages(library(progeny))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(JASPAR2020))
suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(chromVAR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(calibrate))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(cicero))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(BuenColors))
library(readr)
library(RColorBrewer)
library(hdf5r)
library(DESeq2)
library(BiocParallel)
library(zellkonverter)
library(SingleCellExperiment)
library(ggrepel)
library(SingleR)
library(celldex)
library(scran)
library(Azimuth)
library(monocle3)
library(ComplexHeatmap)
library(CellChat)
library(Nebulosa)
library(tidyplots)

#####Initiate#####
set.seed(1234)
knitr::opts_chunk$set(echo = TRUE)

workingDir <- getwd()
dir.create(paste0(workingDir, "session_info"))
sink(paste0(workingDir, "session_info/sessionInfo.txt"))
sessionInfo()
sink()

#####Build SeuratObject#####
# Load raw H5 count matrix 
counts <- Read10X_h5("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/filtered_feature_bc_matrix.h5")
fragpath <- paste0(workingDir, "/atac_fragments.tsv.gz")

# Load annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

#Combine all per_barcode_metrics file to one csv file
#We need to change the barcodes so that the last digit reflects the sample number
replace_last_digit <- function(string, new_digit) {
  gsub("\\d$", new_digit, string)
}

per_barcode_metrics_s005 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S005.csv")
per_barcode_metrics_s006 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S006.csv")
per_barcode_metrics_s007 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S007.csv")
per_barcode_metrics_s008 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S008.csv")
per_barcode_metrics_s009 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S009.csv")
per_barcode_metrics_s010 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S010.csv")
per_barcode_metrics_s013 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S013.csv")
per_barcode_metrics_s014 <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/per_barcode_metrics_S014.csv")

# Use mutate to apply the function to the entire column
#the number after the barcode is the position not the sample number so see aggr file for positions
per_barcode_metrics_s005 <- per_barcode_metrics_s005 |> 
  mutate(barcode = replace_last_digit(barcode, 1)) 

per_barcode_metrics_s006 <- per_barcode_metrics_s006 |> 
  mutate(barcode = replace_last_digit(barcode, 4)) 

per_barcode_metrics_s007 <- per_barcode_metrics_s007 |> 
  mutate(barcode = replace_last_digit(barcode, 3)) 

per_barcode_metrics_s008 <- per_barcode_metrics_s008 |> 
  mutate(barcode = replace_last_digit(barcode, 6)) 

per_barcode_metrics_s009 <- per_barcode_metrics_s009 |> 
  mutate(barcode = replace_last_digit(barcode, 5)) 

per_barcode_metrics_s010 <- per_barcode_metrics_s010 |> 
  mutate(barcode = replace_last_digit(barcode, 7)) 

per_barcode_metrics_s013 <- per_barcode_metrics_s013 |> 
  mutate(barcode = replace_last_digit(barcode, 2)) 

per_barcode_metrics_s014 <- per_barcode_metrics_s014 |> 
  mutate(barcode = replace_last_digit(barcode, 8)) 


per_barcode_metrics <- rbind(per_barcode_metrics_s005, per_barcode_metrics_s006,
                             per_barcode_metrics_s007, per_barcode_metrics_s008,
                             per_barcode_metrics_s009, per_barcode_metrics_s010,
                             per_barcode_metrics_s013, per_barcode_metrics_s014)

#Add other metadata based on position
per_barcode_metrics$position <- as.numeric(sub(".*-(\\d+)$", "\\1", per_barcode_metrics$barcode))
aggr <- read_csv("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/Multiome/HPV-X-Neutra/Data/aggr.csv")
per_barcode_metrics_merged <- merge(per_barcode_metrics, aggr, by = 'position', all = TRUE)
per_barcode_metrics_merged <- per_barcode_metrics_merged |> column_to_rownames("barcode")

per_barcode_metrics <- per_barcode_metrics_merged

rm(per_barcode_metrics_s005, per_barcode_metrics_s006,
   per_barcode_metrics_s007, per_barcode_metrics_s008,
   per_barcode_metrics_s009, per_barcode_metrics_s010,
   per_barcode_metrics_s013, per_barcode_metrics_s014,
   per_barcode_metrics_merged, replace_last_digit)

# Create a twins object containing the RNA data
twins <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  meta.data=per_barcode_metrics)

# Create ATAC assay and add it to the object
twins[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# removing peaks that are in scaffolds, is not so useful for analysis and cause error
DefaultAssay(twins) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(twins))) %in% main.chroms)
twins[["ATAC"]] <- subset(twins[["ATAC"]], features = rownames(twins[["ATAC"]])[keep.peaks])

#####QC DNA Accessibility#####
twins$pct_reads_in_peaks <- twins$atac_peak_region_fragments / twins$atac_fragments.x * 100
twins$blacklist_fraction <- FractionCountsInRegion(object = twins, assay = "ATAC",regions = blacklist_hg38)
twins <- NucleosomeSignal(twins)
twins <- TSSEnrichment(twins)

# View quality control - DNA accessibility metrics
VlnPlot(twins, features = c("nCount_ATAC", "nucleosome_signal","TSS.enrichment", "pct_reads_in_peaks", 
                            "atac_peak_region_fragments", "blacklist_fraction"), pt.size = 0) & 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank())

VlnPlot(twins, features = "nCount_ATAC", pt.size = 0, log = T)
upper_limit <- 5000  # Change this to your desired cutoff value
VlnPlot(twins, features = "nCount_ATAC", pt.size = 0) +
  scale_y_continuous(breaks = seq(0, upper_limit, by = 200), limits = c(0, upper_limit))& 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_line(color = "grey80"),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted"))  # Minor gridlines

VlnPlot(twins, features = "atac_peak_region_fragments", pt.size = 0, log = T)
upper_limit <- 5000  # Change this to your desired cutoff value
VlnPlot(twins, features = "atac_peak_region_fragments", pt.size = 0) +
  scale_y_continuous(breaks = seq(0, upper_limit, by = 200), limits = c(0, upper_limit))& 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_line(color = "grey80"),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted"))  # Minor gridlines

VlnPlot(twins, features = "nucleosome_signal", pt.size = 0, log = T)
upper_limit <- 20  # Change this to your desired cutoff value
VlnPlot(twins, features = "TSS.enrichment", pt.size = 0) +
  scale_y_continuous(breaks = seq(0, upper_limit, by = 1), limits = c(0, upper_limit))& 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_line(color = "grey80"),  # Major gridlines
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted"))  # Minor gridlines

VlnPlot(twins, features = "percent.mt", pt.size = 0, log = F)  +
  theme_bw() +
  scale_fill_manual(values = rep("grey", length(unique(Idents(twins))))) +  # Set all violins to grey
  ylab("Mitochondrial genes (%)") + 
  theme(
    text = element_text(family = "Helvetica", size = 12),  # Base font settings
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),  # X-axis text color
    axis.text.y = element_text(size = 12, colour = "black"),  # Y-axis text color
    axis.title.x = element_blank(),  # X-axis title color
    axis.title.y = element_text(size = 12, colour = "black"),  # Y-axis title color
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

#####QC Gene expression#####
DefaultAssay(twins) <- "RNA"
Idents(twins) <- "orig.ident"

twins[["percent.mt"]] <- PercentageFeatureSet(twins, pattern = "^MT-")
twins[["percent.ribo"]] <- PercentageFeatureSet(twins, pattern = "^RP[SL]")

# View quality control - Gene expression metrics
V1 <- VlnPlot(object = twins, features =  c("nCount_RNA"), group.by = "subject",pt.size = 0, log = T) & 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) 
V2 <- VlnPlot(object = twins, features =  c("nFeature_RNA"), pt.size = 0) & 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) 
V3 <- VlnPlot(object = twins, features =  c("percent.mt", "percent.ribo"), pt.size = 0) & 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
(V1 | V2) / V3

#####Visualize QC#####
VlnPlot(object = twins, features = c("nCount_RNA"), group.by = "subject", 
        pt.size = 0, log = T) & 
  scale_y_continuous(name = "log10(nCount_RNA)", transform = "log10") & 
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) & 
  scale_fill_manual(values = rep(NA, 8)) & 
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0, face = "bold"))

VlnPlot(object = twins, features = c("nFeature_RNA"), group.by = "subject", 
        pt.size = 0, log = T) & 
  scale_y_continuous(name = "log10(nFeature_RNA)", transform = "log10") & 
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) & 
  scale_fill_manual(values = rep(NA, 8)) & 
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0, face = "bold"))

VlnPlot(object = twins, features = c("percent.mt"), group.by = "subject", 
        pt.size = 0, log = F) & 
  scale_y_continuous(name = "Mitochondrial genes (%)") & 
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) & 
  scale_fill_manual(values = rep(NA, 8)) & 
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0, face = "bold"))

VlnPlot(object = twins, features = c("nCount_ATAC"), group.by = "subject", 
        pt.size = 0, log = F) & 
  scale_y_continuous(name = "log10(nCount_ATAC)", transform = "log10") & 
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) & 
  scale_fill_manual(values = rep(NA, 8)) & 
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0, face = "bold"))

VlnPlot(object = twins, features = c("nucleosome_signal"), group.by = "subject", 
        pt.size = 0, log = F) & 
  scale_y_continuous(name = "nucleosome signal") & 
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) & 
  scale_fill_manual(values = rep(NA, 8)) & 
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0, face = "bold"))

VlnPlot(object = twins, features = c("TSS.enrichment"), group.by = "subject", 
        pt.size = 0, log = F) & 
  scale_y_continuous(name = "TSS enrichment") & 
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) & 
  scale_fill_manual(values = rep(NA, 8)) & 
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0, face = "bold"))

# Extract metadata
meta <- twins@meta.data

# Select and combine relevant columns
combined_qc_table <- meta[, c("subject", 
                              "nCount_RNA", 
                              "nFeature_RNA", 
                              "percent.mt", 
                              "nCount_ATAC", 
                              "nucleosome_signal", 
                              "TSS.enrichment")]
# Optionally, save to CSV
write.csv(combined_qc_table, "combined_QC_metrics.csv", row.names = TRUE)

#####Doublets#####
## Doublet identification
library(DoubletFinder)
DefaultAssay(twins) <- "RNA"
twins.diet <- DietSeurat(twins, assays = 'RNA')
twins.list <- SplitObject(twins.diet, split.by = 'subject') #DoubletFinder always per sample
twins.list <- lapply(twins.list, 
                     function(x) {
                       x <- NormalizeData(x) |> 
                         FindVariableFeatures() |> 
                         ScaleData() |> 
                         RunPCA()
                       x <- FindNeighbors(x, dims = 1:15)
                       x <- FindClusters(x)
                       x <- RunUMAP(x, dims = 1:15)
                       
                       sweep.res.list <- paramSweep(x, PCs = 1:15)
                       sweep.res.nsclc <- summarizeSweep(sweep.res.list)
                       bcmvn_nsclc <- find.pK(sweep.res.nsclc)
                       
                       pK <- bcmvn_nsclc %>% 
                         filter(BCmetric == max(BCmetric)) %>% 
                         pull(pK)
                       pK <- as.numeric(as.character(pK[[1]]))
                       
                       annotations <- x@meta.data$seurat_clusters
                       homotic.prop <- modelHomotypic(annotations)
                       nExp_poi <- round(0.05*nrow(x@meta.data))
                       nExp_poi.adj <- round(nExp_poi*(1-homotic.prop))
                       
                       doubletFinder(x, pK = pK, PCs = 1:15, nExp = nExp_poi.adj)})

meta.data.list <- lapply(twins.list, function(x) {
  colnames(x@meta.data) <- sub('_0.25.*', '', colnames(x@meta.data))
  return(x@meta.data)
})

metadata.df <- bind_rows(meta.data.list)

twins <- AddMetaData(twins, metadata.df)
rm(twins.list)
rm(twins.diet)
gc()

metadata.df <- as.data.frame(twins@meta.data)

# Calculate the frequency of 'single' and 'doublet' per 'sample_id'
frequency_table <- metadata.df %>%
  group_by(subject, DF.classifications) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total counts per 'sample_id' and 'predicted.id'
total_counts <- frequency_table %>%
  group_by(subject) %>%
  summarise(total = sum(count), .groups = 'drop')

# Calculate the percentage of each classification
percentage_table <- frequency_table %>%
  left_join(total_counts, by = c("subject")) %>%
  mutate(percentage = (count / total) * 100)

write.csv(percentage_table, "Doublets.csv")

# Filter to include only doublet percentages
doublet_percentage <- percentage_table %>%
  filter(DF.classifications == "Doublet")

# Create the plot for doublet percentages per cluster
ggplot(doublet_percentage, aes(x = subject, y = percentage)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Percentage of Doublet Cells per Subject",
       x = "Subject",
       y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

frequency_table <- frequency_table[frequency_table$DF.classifications == "Doublet", ]
ggplot(frequency_table, aes(x = subject, y = count)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Number of Doublet Cells per Subject",
       x = "",
       y = "Count") +
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, face = "bold", color = "black"))

ggplot(frequency_table, aes(x = subject, y = count)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  theme_bw() +
  labs(title = "Number of Doublet Cells per Subject",
       x = "",
       y = "Count") +
  scale_x_discrete(labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", "Subject 009", "Subject 010", "Subject 013", "Subject 014")) +
  theme(
    text = element_text(family = "Helvetica", size = 12),  # Base font settings
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),  # X-axis text color
    axis.text.y = element_text(size = 10, colour = "black"),  # Y-axis text color
    axis.title.x = element_blank(),  # X-axis title color
    axis.title.y = element_text(size = 12, colour = "black"),  # Y-axis title color
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold")
  )

rm(doublet_frequency, doublet_percentage, frequency_pivot, frequency_table, metadata, percentage_table, total_counts)

#Save before subsetting
dir.create("SeuratObjects")
saveRDS(twins, "SeuratObjects/multiome_pre_processed.RDS")
twins <- readRDS("SeuratObjects/multiome_pre_processed.RDS")

######Filter low quality cells#####
#Filter Low Quality Cells
CellTotal <- print(paste0("Before QC : [",nrow(twins),",",ncol(twins)," ]"))

twins <- subset(
  x = twins,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 15000 & #See figure, I put it on 15000 which should be enough, 25000 was max nCount
    nCount_ATAC > 200 & #I put this on 200 because 1000 was too high, see figure with upperlimit 5000
    nCount_RNA > 100 & #I put this on 100 because 300 was too high, see figure with log scale
    nucleosome_signal < 2 &
    TSS.enrichment > 2 &
    percent.mt < 25 &
    atac_peak_region_fragments > 200 & #I put this on 200 because 2000 was too high, see figure with upperlimit 5000
    atac_peak_region_fragments < 50000 &
    pct_reads_in_peaks > 20 &
    blacklist_fraction < 0.05
)
twins

CellFiltered <- print(paste0("After QC : [",nrow(twins),",",ncol(twins)," ]"))
CellFilter <- rbind(CellTotal, CellFiltered)
write.csv(CellFilter, "Cell_Low_Quality_Filtered.csv")

#Filter Doublets
twins <- subset(twins, subset = DF.classifications == "Singlet")

#####Add HPV nAB titers#####
#Add Ab titers
#Read the Excel file with antibody titers
titer_data <- read_excel("titer_data.xlsx")

#Extract the metadata from the Seurat object
metadata <- twins@meta.data

#Merge the titer data into the Seurat metadata based on the "subject" column
#BE SURE TO transfer ROW NAMES OF METADATA TWINS, these are your Cells!
merged_metadata <- left_join(metadata, titer_data, by = "subject")
rownames(merged_metadata) <- rownames(metadata)

#Update the Seurat object with the new metadata
twins@meta.data <- merged_metadata

# Save the updated Seurat object if needed
saveRDS(twins, "SeuratObjects/HPV-X-Neutra2.RDS")
twins <- readRDS("SeuratObjects/HPV-X-Neutra.RDS")

#####Gene Expression Data processing#####
DefaultAssay(twins) <- "RNA"

# split per twin for later integration datasets
twins[["RNA"]] <- split(twins[["RNA"]], f = twins$twin)

# Normalization
twins <- SCTransform(twins, conserve.memory=TRUE)
# Perform dimensionality reduction by Principal Component Analysis (PCA) 
twins <- RunPCA(twins)
# look at how much variation there is along each of the 50 PCs we computed
ElbowPlot(twins, ndims=50)
# UMAP (uniform manifold approximation and projection)
#We will perform dimensional reduction on both assays independently, using standard approaches for RNA and ATAC-seq data
twins <- RunUMAP(twins, dims = 1:50, reduction.name = 'umap.rnaNI', reduction.key = 'rnaNIUMAP_')
twins <- FindNeighbors(twins, dims = 1:30)
twins <- FindClusters(twins, resolution = 0.5)
p1 <- DimPlot(twins, reduction = "umap.rnaNI", group.by = "subject")
p2 <- DimPlot(twins, reduction = "umap.rnaNI", group.by = "vaccine")
p3 <- DimPlot(twins, reduction = "umap.rnaNI", group.by = "twin")
p1+p2+p3

plot <- DimPlot(twins, reduction = "umap.rnaNI", group.by = "subject", 
                #split.by = "vaccine", 
                label = F, repel = T, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12),# Adjust subtitle size if necessary
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
        ) 
plot <- plot + annotate("segment", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]), 
                        xend = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]), 
                        yend = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]), 
           xend = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]), 
           y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]), 
           yend = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot

#Integration
twins <- IntegrateLayers(object = twins, method = HarmonyIntegration,
                         orig.reduction = "pca", new.reduction = 'harmony',
                         assay = "SCT", verbose = FALSE)
twins <- FindNeighbors(twins, reduction = "harmony", dims = 1:30)
twins <- FindClusters(twins, resolution = 0.6)
twins <- RunUMAP(twins, dims = 1:30, reduction = "harmony", reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

pa <- DimPlot(twins, reduction = "umap.rna", group.by = "subject")
pb <- DimPlot(twins, reduction = "umap.rna", group.by = "vaccine")
pc <- DimPlot(twins, reduction = "umap.rna", group.by = "twin")
pa+pb+pc
gc()

######DNA accessibility Data processing#####
DefaultAssay(twins) <- "ATAC"
twins <- FindTopFeatures(twins, min.cutoff = 5)
twins <- RunTFIDF(twins)
twins <- RunSVD(twins)
twins <- RunUMAP(twins, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atacNI", reduction.key = "atacNIUMAP_")
p1 <- DimPlot(twins, reduction = "umap.atacNI", group.by = "subject")
p2 <- DimPlot(twins, reduction = "umap.atacNI", group.by = "vaccine")
p3 <- DimPlot(twins, reduction = "umap.atacNI", group.by = "twin")
p1+p2+p3

plot <- DimPlot(twins, reduction = "umap.atacNI", group.by = "subject", 
                #split.by = "vaccine", 
                label = F, repel = T, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12),# Adjust subtitle size if necessary
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  ) 
plot <- plot + annotate("segment", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]), 
                        xend = min(twins@reductions$umap.atacNI@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]), 
                        yend = min(twins@reductions$umap.atacNI@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]), 
           xend = min(twins@reductions$umap.atacNI@cell.embeddings[,1]), 
           y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]), 
           yend = min(twins@reductions$umap.atacNI@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot

# Save reductions before splitting
reduction_list <- list(
  pca = twins[["pca"]],
  umap_rnaNI = twins[["umap.rnaNI"]],
  harmony = twins[["harmony"]],
  umap_rna = twins[["umap.rna"]],
  lsi = twins[["lsi"]],
  umap_atacNI = twins[["umap.atacNI"]]
)

saveRDS(twins, "SeuratObjects/intermediatebeforeatacintegration.RDS")

DefaultAssay(twins) <- "ATAC"
twins.list <- SplitObject(twins, split.by = "twin")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = twins.list,
  anchor.features = rownames(twins),
  reduction = "rlsi",
  dims = 2:30
)

saveRDS(integration.anchors, "SeuratObjects/IntegrationAnchorsATAC.RDS")
rm(twins.list)
gc()

# integrate LSI embeddings
twins <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = twins[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# Restore the saved reductions
twins[["pca"]] <- reduction_list$pca
twins[["umap.rnaNI"]] <- reduction_list$umap_rnaNI
twins[["harmony"]] <- reduction_list$harmony
twins[["umap.rna"]] <- reduction_list$umap_rna
twins[["lsi"]] <- reduction_list$lsi
twins[["umap.atacNI"]] <- reduction_list$umap_atacNI

rm(integration.anchors, twins.list, reduction_list)

# create a new UMAP using the integrated embeddings
twins <- RunUMAP(twins, reduction = "integrated_lsi", dims = 2:30, reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')
p1 <- DimPlot(twins, reduction = "umap.atac", group.by = "subject")
p2 <- DimPlot(twins, reduction = "umap.atac", group.by = "vaccine")
p3 <- DimPlot(twins, reduction = "umap.atac", group.by = "twin")
p1+p2+p3

saveRDS(twins, "SeuratObjects/HPV-X-Neutra.RDS")

######Build a joint neighbor graph using both assays#####
twins <- FindMultiModalNeighbors(
  object = twins,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
twins <- RunUMAP(twins, nn.name = "weighted.nn",reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = TRUE)

#Clustering
twins <- FindNeighbors(twins, dims = 1:30)
# You should make sure your assay is set correctly (the assay that you originally run PCA).
DefaultAssay(twins) <- "SCT"
# Now we cluster the resultant graph.
twins  <- FindClusters(twins, resolution = 0.5)

#Cluster visualization
# visualize the UMAP
p1 <- DimPlot(twins, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(twins, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(twins, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p4 <- DimPlot(twins, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Clusters")
p1 <- DimPlot(twins, reduction = "wnn.umap", group.by = "subject", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Subject")
p2 <- DimPlot(twins, reduction = "wnn.umap", group.by = "vaccine", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Vaccine")
p3 <- DimPlot(twins, reduction = "wnn.umap", group.by = "twin", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Twin")
(p1+p3)/(p2+p4)

DimPlot(twins, reduction = "umap.rnaNI", group.by = "subject", raster = F) 

plot <- DimPlot(twins, reduction = "umap.rnaNI", group.by = "subject",
                label = F, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12)) # Adjust subtitle size if necessary
plot <- plot + annotate("segment", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]), 
                        xend = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]), 
                        yend = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]), 
           xend = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]), 
           y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]), 
           yend = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins@reductions$umap.rnaNI@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(twins@reductions$umap.rnaNI@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot

plot <- DimPlot(twins, reduction = "umap.atacNI", group.by = "subject",
                label = F, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12)) # Adjust subtitle size if necessary
plot <- plot + annotate("segment", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]), 
                        xend = min(twins@reductions$umap.atacNI@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]), 
                        yend = min(twins@reductions$umap.atacNI@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]), 
           xend = min(twins@reductions$umap.atacNI@cell.embeddings[,1]), 
           y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]), 
           yend = min(twins@reductions$umap.atacNI@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins@reductions$umap.atacNI@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(twins@reductions$umap.atacNI@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot

plot <- DimPlot(twins, reduction = "wnn.umap", group.by = "seurat_clusters",
                label = F, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12)) # Adjust subtitle size if necessary
plot <- plot + annotate("segment", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
                        xend = min(twins@reductions$wnn.umap@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
                        yend = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
           xend = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
           yend = min(twins@reductions$wnn.umap@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot

cluster_frq <- twins@meta.data %>%
  group_by(seurat_clusters, vaccine) %>%
  summarise(count=n()) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(data_set = "HPV-X-Neutra")

cluster_frq

freq_plot <- ggplot(cluster_frq, aes(x = seurat_clusters, y = relative_freq)) +
  geom_col(color="black", aes(fill=vaccine), position = position_stack(reverse = TRUE)) +
  ylab("Count") +
  xlab("Seurat Clusters") +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73"), 
                    labels = c("Cervarix", "Gardasil")) + 
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = 12),  # Base font settings
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y=element_text(size = 10, colour = "black"),
    axis.text.x=element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    axis.ticks.x=element_blank(),
    axis.line = element_blank()) +
  guides(fill = guide_legend(title = "Vaccine")) 

freq_plot

cluster_frq <- twins@meta.data %>%
  group_by(seurat_clusters, subject) %>%
  summarise(count=n()) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(data_set = "HPV-X-Neutra")

cluster_frq

freq_plot <- ggplot(cluster_frq, aes(x = seurat_clusters, y = relative_freq)) +
  geom_col(color="black", aes(fill=subject), position = position_stack(reverse = TRUE)) +
  ylab("Count") +
  xlab("Seurat Clusters") +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                               "#0072B2", "#D55E00", "#CC79A7", "#999999"),
                    labels = c("Subject 005", "Subject 006", "Subject 007", "Subject 008", 
                               "Subject 009", "Subject 010", "Subject 013", "Subject 014")) + 
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = 12),  # Base font settings
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y=element_text(size = 10, colour = "black"),
    axis.text.x=element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    axis.ticks.x=element_blank(),
    axis.line = element_blank()) +
  guides(fill = guide_legend(title = "Subject")) 

freq_plot

saveRDS(twins, "SeuratObjects/HPV-X-Neutra.RDS")

#####Annotation Azimuth#####
# load reference
reference <- LoadH5Seurat("/Users/valentino/Library/CloudStorage/OneDrive-UGent/Systems Vaccinology/Bioinformatics/pbmc_multimodal.h5seurat")

DefaultAssay(twins) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = twins,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictionsl1 <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l1,
  weight.reduction = twins[['pca']],
  dims = 1:50
)

# Get the column names of predictions
pred_cols <- colnames(predictionsl1)

# Add prefix to each column name
pred_cols_new <- paste("PBMC_l1_", pred_cols, sep = "")

# Rename columns in predictions object
colnames(predictionsl1) <- pred_cols_new

predictionsl2 <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = twins[['pca']],
  dims = 1:50
)

# Get the column names of predictions
pred_cols <- colnames(predictionsl2)

# Add prefix to each column name
pred_cols_new <- paste("PBMC_l2_", pred_cols, sep = "")

# Rename columns in predictions object
colnames(predictionsl2) <- pred_cols_new

predictionsl3 <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l3,
  weight.reduction = twins[['pca']],
  dims = 1:50
)

# Get the column names of predictions
pred_cols <- colnames(predictionsl3)

# Add prefix to each column name
pred_cols_new <- paste("PBMC_l3_", pred_cols, sep = "")

# Rename columns in predictions object
colnames(predictionsl3) <- pred_cols_new

twins <- AddMetaData(
  object = twins,
  metadata = c(predictionsl1, predictionsl2, predictionsl3)
)

# set the cell identities to the cell type predictions
Idents(twins) <- "PBMC_l2_predicted.id"

# remove cell types with lower than 15 cells
remCells <- twins@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(PBMC_l2_predicted.id) %>%
  mutate(total_cells = n_distinct(cell)) %>%
  mutate(removeCells = ifelse(total_cells > 16, total_cells, 0)) %>%
  dplyr::select(cell, removeCells) %>% 
  column_to_rownames("cell")

twins <- AddMetaData(twins, metadata = remCells)
twins <- subset(twins, subset = removeCells > 1 )
twins

Idents(twins) <- "PBMC_l2_predicted.id"
levels(twins@active.ident)

# View annotations
p1 <- DimPlot(twins, reduction = "wnn.umap", group.by = "PBMC_l1_predicted.id", label = TRUE, label.size = 3, repel = TRUE) + ggtitle("Cells level 1")
p2 <- DimPlot(twins, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE) + ggtitle("Twins")
p3 <- DimPlot(twins, reduction = "wnn.umap", group.by = "PBMC_l3_predicted.id", label = TRUE, label.size = 3, repel = TRUE) + ggtitle("Cells level 3")
p1/p2/p3

saveRDS(twins, "HPV-X-Neutra.RDS")

#####Check Annotations#####
# Check if level 2 and level 1 corresponds
B <- subset(twins, subset = PBMC_l1_predicted.id == "B")
DimPlot(B, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE)

CD4 <- subset(twins, subset = PBMC_l1_predicted.id == "CD4 T")
DimPlot(CD4, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE)

Mono <- subset(twins, subset = PBMC_l1_predicted.id == "Mono")
DimPlot(Mono, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE)

CD8 <- subset(twins, subset = PBMC_l1_predicted.id == "CD8 T")
DimPlot(CD8, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE)

NK <- subset(twins, subset = PBMC_l1_predicted.id == "NK")
DimPlot(NK, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE)

DC <- subset(twins, subset = PBMC_l1_predicted.id == "DC")
DimPlot(DC, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE)

otherT <- subset(twins, subset = PBMC_l1_predicted.id == "other T")
DimPlot(otherT, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id", label = TRUE, label.size = 3, repel = TRUE)

#Finding markers for each cluster
#Check if markers correspond with annotated cells
#I also checked that level 1 and level 2 correspond
DefaultAssay(twins) <- "SCT"
Idents(twins) <- "PBMC_l2_predicted.id"
twins <- PrepSCTFindMarkers(twins, assay = "SCT")
gc()
cluster_markers <- FindAllMarkers(twins,
                                  logfc.threshold = 0, #To test all possible genes (very slow), set this to 0
                                  only.pos = F, #and set this to F
                                  recorrect_umi = F)
# And lets look at the top markers (by logFC) for each cluster. 
cluster_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

# Let"s plot the top gene for every cluster
# Plot options: Feature plot, Violin plots, Ridge plots, Heatmap, Dot Plot
top_genes <- cluster_markers %>% 
  group_by(cluster) %>% 
  top_n(1, avg_log2FC) %>% 
  pull(gene)

plotLS <- list()

plotLS[[1]] <- FeaturePlot(twins, features=top_genes, cols = c("lightgrey", "red"), pt.size=0.25, reduction = "wnn.umap",
                           combine = T)

plotLS[[2]] <- VlnPlot(twins, features="CD8B", pt.size=0.25) + NoLegend()

plotLS[[3]] <- RidgePlot(twins, features="IGHA2")+ NoLegend()

top_genes <- cluster_markers %>% 
  group_by(cluster) %>% 
  top_n(5, avg_log2FC) %>% 
  pull(gene) %>%
  unique() #in case any markers are duplicated

plotLS[[4]] <- DoHeatmap(twins, features=top_genes, label = T, size = 1)

plotLS[[5]] <- DotPlot(twins, features=top_genes, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size=2))

headings <- c("1: Feature Plot", "2: Violin plot", 
              "3: Ridge plots", "4: Heatmap",
              "5: Dot Plot")

for (i in 1:length(plotLS)) {
  cat("##### ",headings[i],"\n")
  print(plotLS[[i]])
  cat('\n\n')
}

#Per cluster
#Extract Celtypes
celtypes <- levels(twins@active.ident)
celtypes <- data.frame(celtypes)
celtypes

CD4 <- c("CD4 Naive", 'CD4 TCM', 'Treg', 'CD4 TEM', 'CD4 Proliferating', 'CD4 CTL')
CD4 <- celtypes[celtypes$celtypes %in% CD4, , drop = FALSE]
CD4 <- unique(CD4$celtypes)
CD8 <- c('CD8 TEM', 'CD8 Naive', 'CD8 TCM', 'CD8 Proliferating')
CD8 <- celtypes[celtypes$celtypes %in% CD8, , drop = FALSE]
CD8 <- unique(CD8$celtypes)
otherT <- c('gdT', 'MAIT', 'dnT')
otherT <- celtypes[celtypes$celtypes %in% otherT, , drop = FALSE]
otherT <- unique(otherT$celtypes)
NK <- c('NK', 'NK_CD56bright', 'NK Proliferating')
NK <- celtypes[celtypes$celtypes %in% NK, , drop = FALSE]
NK <- unique(NK$celtypes)
B <- c('B naive', 'B memory', 'B intermediate', 'Plasmablast')
B <- celtypes[celtypes$celtypes %in% B, , drop = FALSE]
B <- unique(B$celtypes)
Mono <- c('CD14 Mono', 'CD16 Mono')
Mono <- celtypes[celtypes$celtypes %in% Mono, , drop = FALSE]
Mono <- unique(Mono$celtypes)
DC <- c('pDC', 'cDC1', 'cDC2', 'ASDC')    
DC <- celtypes[celtypes$celtypes %in% DC, , drop = FALSE]
DC <- unique(DC$celtypes)
other <- c('ILC', 'HSPC', 'Platelet', 'Eryth')
other <- celtypes[celtypes$celtypes %in% other, , drop = FALSE]
other <- unique(other$celtypes)

#Topgenes per category
top_genes_CD4 <- cluster_markers %>% 
  filter(cluster %in% CD4) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

top_genes_CD8 <- cluster_markers %>% 
  filter(cluster %in% CD8) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

top_genes_otherT <- cluster_markers %>% 
  filter(cluster %in% otherT) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

top_genes_NK <- cluster_markers %>% 
  filter(cluster %in% NK) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

top_genes_B <- cluster_markers %>% 
  filter(cluster %in% B) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

top_genes_Mono <- cluster_markers %>% 
  filter(cluster %in% Mono) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

top_genes_DC <- cluster_markers %>% 
  filter(cluster %in% DC) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

top_genes_other <- cluster_markers %>% 
  filter(cluster %in% other) %>%
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>% 
  pull(gene)%>%
  unique()

#Subset per celltype
CD4 <- subset(twins, idents = CD4)
CD8 <- subset(twins, idents = CD8)
otherT <- subset(twins, idents = otherT)
NK <- subset(twins, idents = NK)
B <- subset(twins, idents = B)
Mono <- subset(twins, idents = Mono)
DC <- subset(twins, idents = DC)
other <- subset(twins, idents = other)

#Dotplot per celltype
p1 <- DotPlot(CD4, features=top_genes_CD4, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("CD4 cells")

p2 <- DotPlot(CD8, features=top_genes_CD8, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("CD8 cell")

p3 <- DotPlot(NK, features=top_genes_NK, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("Natural Killer cells")

p4 <- DotPlot(B, features=top_genes_B, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("B cells")

p5 <- DotPlot(Mono, features=top_genes_Mono, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("Monocytes")

p6 <- DotPlot(otherT, features=top_genes_otherT, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("Other T cells")

p7 <- DotPlot(DC, features=top_genes_DC, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("Dendritic Cells")

p8 <- DotPlot(other, features=top_genes_other, cols=c("lightgrey", "red")) + 
  theme(axis.text = element_text(angle=45, hjust=1, size = 5)) + ggtitle("Other cells")

(((p1+p2)/(p3+p6))/(p5+p7))/(p4+p8)

#####Final UMAP#####
plot <- DimPlot(twins, reduction = "wnn.umap", group.by = "subject", 
                #split.by = "vaccine", 
                label = F, repel = T, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
        )
plot <- plot + annotate("segment", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
                        xend = min(twins@reductions$wnn.umap@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
                        yend = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
           xend = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
           yend = min(twins@reductions$wnn.umap@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot

plot <- DimPlot(twins, reduction = "wnn.umap", group.by = "PBMC_l1_predicted.id", 
                #split.by = "vaccine", 
                label = T, repel = T, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
        )
plot <- plot + annotate("segment", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
                        xend = min(twins@reductions$wnn.umap@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
                        yend = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
           xend = min(twins@reductions$wnn.umap@cell.embeddings[,1]), 
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]), 
           yend = min(twins@reductions$wnn.umap@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins@reductions$wnn.umap@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(twins@reductions$wnn.umap@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot
ggsave(filename = "Figure2D_final_UAMP.pdf", plot = plot, width = 40, height = 40, dpi = 320, units = "cm", device = "pdf")

#####Frequencies######
cluster_frq <- twins@meta.data %>%
  group_by(vaccine, PBMC_l2_predicted.id) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(relative_freq_percent = relative_freq * 100)

#Calculate per subject
cluster_frq <- twins@meta.data %>%
  group_by(twin, subject, vaccine, PBMC_l2_predicted.id) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(relative_freq_percent = relative_freq * 100)

# Filter for specific cell types
selected_celltypes <- c("CD14 Mono", "CD16 Mono", "NK", "NK_CD56bright","ASDC", "cDC1", "cDC2", "pDC",
                        "CD4 Naive", "CD4 Proliferating", "CD4 TEM", "CD4 TCM", "Treg", "CD4 CTL",
                        "CD8 Naive", "CD8 Proliferating", "CD8 TEM", "CD8 TCM", 
                        "B naive", "B intermediate", "B memory", "Plasmablast")
#"dnT", "gdT", "MAIT", "ILC",
#"Platelet", "Eryth","HSPC")  
cluster_frq_subset <- cluster_frq %>%
  filter(PBMC_l2_predicted.id %in% selected_celltypes)

# Set the order of the 'PBMC_l2_predicted.id' factor
cluster_frq_subset$PBMC_l2_predicted.id <- factor(cluster_frq_subset$PBMC_l2_predicted.id, levels = selected_celltypes)

# Plot bar charts with significance markers
vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")

cluster_frq_subset %>% 
  tidyplot(x = vaccine, y = relative_freq, color = vaccine) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_data_labels_repel(label = twin) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "fisher_test", p.adjust.method = "bonferroni", label = "p.signif", hide_info = T) %>%
  add_test_pvalue(method = "t_test", p.adjust.method = "bonferroni",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Proportion") %>%
  adjust_legend_title("Vaccine") %>%
  split_plot(by = PBMC_l2_predicted.id)

ggsave(filename = "Figure_2B_CellFreq_Differences.pdf", plot = plot, width = 40, height = 20, dpi = 320, units = "cm", device = "pdf")

# Save the results to a CSV file
write.csv(cluster_frq, "RelativeCellFrequencies_PerSubject.csv", row.names = FALSE)

#####Subsets of cells#####
Idents(twins) <- "PBMC_l2_predicted.id"
levels(twins@active.ident)

otherT <- c('gdT', 'MAIT', 'dnT')
other <- c('HSPC', 'Platelet', 'Eryth')
CD8 <- c('CD8 TEM', 'CD8 Naive', 'CD8 TCM', 'CD8 Proliferating')

CD4 <- c("CD4 Naive", 'CD4 TCM', 'Treg', 'CD4 TEM', 'CD4 Proliferating', 'CD4 CTL')
twins.cd4 <- subset(twins, idents = CD4)
DimPlot(twins.cd4, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id")
twins.cd4 <- FindMultiModalNeighbors(
  object = twins.cd4,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
twins.cd4 <- RunUMAP(twins.cd4, nn.name = "weighted.nn", reduction.name = 'wnn.umapcd4', reduction.key = 'rnaUMAPcd4_')
#Clustering
ElbowPlot(twins.cd4, ndims=50)
twins.cd4 <- FindNeighbors(twins.cd4, dims = 1:30)
# You should make sure your assay is set correctly (the assay that you originally run PCA).
DefaultAssay(twins.cd4) <- "SCT"
# Now we cluster the resultant graph.
twins.cd4 <- FindClusters(twins.cd4, resolution = 0.5)
DimPlot(twins.cd4, reduction = "wnn.umapcd4", group.by = "subject")
saveRDS(twins.cd4, "SeuratObjects/CD4.RDS")

B <- c('B naive', 'B memory', 'B intermediate', 'Plasmablast')
twins.b <- subset(twins, idents = B)
DimPlot(twins.b, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id")
twins.b <- FindMultiModalNeighbors(
  object = twins.b,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
twins.b <- RunUMAP(twins.b, nn.name = "weighted.nn", reduction.name = 'wnn.umapb', reduction.key = 'rnaUMAPb_')
DimPlot(twins.b, reduction = "wnn.umapb", group.by = "PBMC_l2_predicted.id")
ElbowPlot(twins.b, ndims=50)
twins.b <- FindNeighbors(twins.b, dims = 1:30)
# You should make sure your assay is set correctly (the assay that you originally run PCA).
DefaultAssay(twins.b) <- "SCT"
# Now we cluster the resultant graph.
twins.b <- FindClusters(twins.b, resolution = 0.5)
DimPlot(twins.b, reduction = "wnn.umapb", group.by = "seurat_clusters")
saveRDS(twins.b, "SeuratObjects/B.RDS")

Innate <- c('CD14 Mono', 'CD16 Mono','pDC', 'cDC1', 'cDC2', 'ASDC','NK', 'NK_CD56bright','ILC')
twins.innate <- subset(twins, idents = Innate)
DimPlot(twins.innate, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id")
twins.innate <- FindMultiModalNeighbors(
  object = twins.innate,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
twins.innate <- RunUMAP(twins.innate, nn.name = "weighted.nn", reduction.name = 'wnn.umapinnate', reduction.key = 'rnaUMAPinnate_')
DimPlot(twins.innate, reduction = "wnn.umapinnate", group.by = "PBMC_l2_predicted.id")
saveRDS(twins.innate, "SeuratObjects/Innate.RDS")

DCB <- c('pDC', 'cDC1', 'cDC2', "B naive", "B memory", "B intermediate", "Plasmablast")
twins.DC <- subset(twins, idents = DCB)
DimPlot(twins.DC, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id")
twins.DC <- FindMultiModalNeighbors(
  object = twins.DC,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
twins.DC <- RunUMAP(twins.DC, nn.name = "weighted.nn", reduction.name = 'wnn.umapDC', reduction.key = 'rnaUMAPDC_')
DimPlot(twins.DC, reduction = "wnn.umapDC", group.by = "PBMC_l2_predicted.id")
saveRDS(twins.DC, "SeuratObjects/DC_B.RDS")
twins.DC <- readRDS("SeuratObjects/DC_B.RDS")
twins <- readRDS("SeuratObjects/HPV-X-Neutra.RDS")

#####Monocyte and Dendritic Cell Annotation Check####
MonoDCB <- c('CD14 Mono', 'CD16 Mono','pDC', 'cDC1', 'cDC2', "B naive", "B memory", "B intermediate", "Plasmablast")
twins.DC <- subset(twins, idents = MonoDCB)
DimPlot(twins.DC, reduction = "wnn.umap", group.by = "PBMC_l2_predicted.id")
twins.DC <- FindMultiModalNeighbors(
  object = twins.DC,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
twins.DC <- RunUMAP(twins.DC, nn.name = "weighted.nn", reduction.name = 'wnn.umapDC', reduction.key = 'rnaUMAPDC_')
DimPlot(twins.DC, reduction = "wnn.umapDC", group.by = "PBMC_l2_predicted.id")
saveRDS(twins.DC, "SeuratObjects/Mono_DC_B.RDS")

dc_markers <- c("ITGAM", #Myeloid
                "CD86", "CD80", #General
                "HLA-DRA", "HLA-DPA1", "HLA-DQB1",  #Antigen Presentation
                "ITGAX", #Classical
                "FLT3", "CD24","IRF8", #cDC1
                "CD1C", "FCER1A","CLEC10A", "TLR8",#cDC2
                "IL3RA","NRP1","IRF7","IRF4", "TLR7","TLR9", "LAMP3","TCF4", #pDC
                "CD14", "FCGR3A", "CCR2", #Monocyte-derived DC (also CD1C, ITGAX, and HLA-DRA)
                "CD163", #DC3 (also CD1C, FCGR3A and CD14
                "CCR7") #Migratory DC (also CD86, dc0 and HLA-DRA)

classical <- c("CD14", "LYZ", "S100A9", "FCGR3A", "CCR2", "TREM1") #CD14++ CD16-, Inflammatory and phagocytotic FCGR3A = CD16
intermediate <- c("CD14", "FCGR3A", "HLA-DRA", "ITGAX") #CD14High CD16 Pos, Ag Presentation, transitional
nonclassical <- c("CD14", "FCGR3A", "CX3CR1", "ITGAL","S100A8", "S100A9") #CD14+ CD16++, antiviral, tissue repair, patrolling
macrophage <- c("CD14", "FCGR3A","CD163", "F13A1", "APOE", "C1QA", "C1QB")
M1 <- c("IL1B", "TNF", "NOS2", "CXCL9", "CXCL10") #inflammatory
M2 <- c("MRC1", "CD163", "IL10", "TREM2", "ARG1") #anti-inflammatory
macro <- c(macrophage, M1, M2) %>% unique()
all <- c(classical, intermediate, nonclassical, macro) %>% unique()

library(viridis)
DotPlot(twins.DC, features = dc_markers, group.by = "PBMC_l2_predicted.id") + 
  scale_color_viridis(option = "H", discrete = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "Dendritic Cell Markers", y = "Seurat Clusters")

#####Final UMAP DC and B#####

plot <- DimPlot(twins.DC, reduction = "wnn.umapDC", group.by = "PBMC_l2_predicted.id", 
                #split.by = "vaccine", 
                label = T, repel = T, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  )
plot <- plot + annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
                        xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) + 0.8, #This is length of x arrow
                        y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
                        yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) -1, 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) -2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 2 ,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 1, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2.5, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 0.5, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot
ggsave(filename = "final_UAMP_DCandB.pdf", plot = plot, width = 18, height = 18, dpi = 320, units = "cm", device = "pdf")

#####Frequencies######
#All
cluster_frq <- twins@meta.data %>%
  group_by(vaccine, PBMC_l2_predicted.id) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(relative_freq_percent = relative_freq * 100)

#Calculate per subject
cluster_frq <- twins@meta.data %>%
  group_by(twin, subject, vaccine, PBMC_l2_predicted.id) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(relative_freq_percent = relative_freq * 100)

# Filter for specific cell types
selected_celltypes <- c("CD14 Mono", "CD16 Mono", "NK", "NK_CD56bright","ASDC", "cDC1", "cDC2", "pDC",
                        "CD4 Naive", "CD4 Proliferating", "CD4 TEM", "CD4 TCM", "Treg", "CD4 CTL",
                        "CD8 Naive", "CD8 Proliferating", "CD8 TEM", "CD8 TCM", 
                        "B naive", "B intermediate", "B memory", "Plasmablast")
#"dnT", "gdT", "MAIT", "ILC",
#"Platelet", "Eryth","HSPC")  
cluster_frq <- cluster_frq %>%
  filter(PBMC_l2_predicted.id %in% selected_celltypes)
# Set the order of the 'PBMC_l2_predicted.id' factor
cluster_frq$PBMC_l2_predicted.id <- factor(cluster_frq$PBMC_l2_predicted.id, levels = selected_celltypes)

vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")

plot <- cluster_frq %>% 
  tidyplot(x = vaccine, y = relative_freq, color = vaccine) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "fisher_test", p.adjust.method = "bonferroni", label = "p.signif", hide_info = T) %>%
  add_test_pvalue(method = "t_test", p.adjust.method = "bonferroni",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Proportion", family = "helvetica", fontsize = 8, color = "black") %>%
  adjust_y_axis(limits = c(0,NA), padding = c(0,0.15)) %>%
  adjust_legend_title("Vaccine", family = "helvetica", fontsize = 10, color = "black") %>%
  adjust_font(family = "Helvetica", fontsize = 8, color = "black", face = "bold") %>%
  adjust_size(width = 100, height = 100, unit = "mm") %>%
  split_plot(by = PBMC_l2_predicted.id)
plot

#DC and B
cluster_frq <- twins.DC@meta.data %>%
  group_by(vaccine, PBMC_l2_predicted.id) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(relative_freq_percent = relative_freq * 100)

#Calculate per subject
cluster_frq <- twins.DC@meta.data %>%
  group_by(twin, subject, vaccine, PBMC_l2_predicted.id) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(relative_freq_percent = relative_freq * 100)

# Filter for specific cell types
selected_celltypes <- c("cDC1", "cDC2", "pDC",
                        "B naive", "B intermediate", "B memory", "Plasmablast")

# Set the order of the 'PBMC_l2_predicted.id' factor
cluster_frq$PBMC_l2_predicted.id <- factor(cluster_frq$PBMC_l2_predicted.id, levels = selected_celltypes)

vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")

plot <- cluster_frq %>% 
  tidyplot(x = vaccine, y = relative_freq, color = vaccine) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "fisher_test", p.adjust.method = "bonferroni", label = "p.signif", hide_info = T) %>%
  add_test_pvalue(method = "t_test", p.adjust.method = "bonferroni",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Proportion", family = "helvetica", fontsize = 8, color = "black") %>%
  adjust_y_axis(limits = c(0,NA), padding = c(0,0.15)) %>%
  adjust_legend_title("Vaccine", family = "helvetica", fontsize = 10, color = "black") %>%
  adjust_font(family = "Helvetica", fontsize = 8, color = "black", face = "bold") %>%
  adjust_size(width = 100, height = 100, unit = "mm") %>%
  split_plot(by = PBMC_l2_predicted.id)

plot
ggsave(filename = "DCandB_CellFreq_Differences.pdf", plot = plot, width = 20, height = 20, dpi = 320, units = "cm", device = "pdf")

# Save the results to a CSV file
write.csv(cluster_frq, "RelativeCellFrequencies_PerSubject.csv", row.names = FALSE)

#####DEG and GSEA per cell group#####
library(tidyverse)
library(clusterProfiler)
library(RCy3)
library(enrichplot)
library(pheatmap)

#DEG
DefaultAssay(twins.DC) <- "SCT"
Idents(twins.DC) <- "PBMC_l2_predicted.id"
subsets <- levels(twins.DC@active.ident)

twins.DC$celltype.vaccine <- paste(twins.DC$PBMC_l2_predicted.id, twins.DC$vaccine, sep = "_")
Idents(twins.DC) <- "celltype.vaccine"
levels(twins.DC@active.ident)
DefaultAssay(twins.DC) <- "SCT"

degs.list <- list()
for(i in seq_along(subsets)) {
  degs.list[[i]] <- FindMarkers(twins.DC,
                                ident.1 = paste(subsets[i], "CERVARIX", sep = "_" ),
                                ident.2 = paste(subsets[i], "GARDASIL", sep = "_" ), recorrect_umi = F)
}
names(degs.list) <- subsets

# Combine all data frames with an identifier and save as CSV
# Add gene column and remove rownames in each subset
for (i in seq_along(degs.list)) {
  degs.list[[i]]$gene <- rownames(degs.list[[i]])
  rownames(degs.list[[i]]) <- NULL
}
combined_degs <- bind_rows(degs.list, .id = "subset")
write.csv(combined_degs, file = "DEGs_All_Cells.csv", row.names = FALSE)

# Initialize an empty list to store the plots
plots_list <- list()

# Loop through each subset in degs.list
for(i in seq_along(degs.list)) {
  # Extract the current data frame and reset row names
  degs.list[[i]]$gene <- rownames(degs.list[[i]])
  rownames(degs.list[[i]]) <- NULL
  
  # Categorize the points
  degs.list[[i]]$category <- ifelse(degs.list[[i]]$p_val < 0.05 & degs.list[[i]]$avg_log2FC < -1, "Significant < -1",
                                    ifelse(degs.list[[i]]$p_val < 0.05 & degs.list[[i]]$avg_log2FC > 1, "Significant > 1",
                                           "Not Significant"))
  
  # Create the ggplot for the current subset
  p <- ggplot(degs.list[[i]], aes(avg_log2FC, -log10(p_val))) + 
    geom_point(aes(color = category), size = 0.5, alpha = 0.5) + 
    scale_color_manual(values = c("Significant < -1" = "blue", "Significant > 1" = "red", "Not Significant" = "grey")) +
    theme_bw() +
    ylab("-log10(p-value)") + 
    xlab("log2(Fold Change)") +
    geom_text_repel(aes(label = ifelse(p_val < 0.01, gene, "")), 
                    colour = "black", size = 3) +
    geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme(legend.position = "none") +
    ggtitle(names(degs.list)[i])
  
  # Store the plot in the list
  plots_list[[names(degs.list)[i]]] <- p
}

volcano_innate <- (plots_list[[3]] + plots_list[[9]]) / 
  (plots_list[[10]] + plots_list[[19]]) / 
  (plots_list[[24]] + plots_list[[8]]) / 
  (plots_list[[21]] + plots_list[[28]])

volcano_cd4 <- (plots_list[[2]] + plots_list[[20]]) / 
  (plots_list[[5]] + plots_list[[6]]) / 
  (plots_list[[11]] + plots_list[[22]])

volcano_b <- (plots_list[[1]] + plots_list[[18]]) / 
  (plots_list[[14]] + plots_list[[27]])

#GSEA
# Load gene set - REACTOME
gene_set_r <- msigdbr::msigdbr(category = 'C2', subcategory = 'CP:REACTOME') %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  unique()

# Ranked list of genes
ranks.list <- list()
for (i in seq_along(degs.list)) {
  ranks.list[[i]] <- degs.list[[i]] %>%
    mutate(p_val = ifelse(p_val == 0, .Machine$double.xmin*row_number(), p_val)) %>% 
    mutate(rank = -log10(p_val)*sign(avg_log2FC)) %>% 
    dplyr::select(rank) %>% 
    rownames_to_column() %>% 
    arrange(dplyr::desc(rank)) %>% 
    deframe() 
}
names(ranks.list) <- names(degs.list)

# Gene set enrichment analysis
gsea.list <- compareCluster(ranks.list, fun = 'GSEA', 
                            TERM2GENE = gene_set_r, eps = 0)

#Save
GSEA_df <- gsea.list@compareClusterResult
GSEA_df$Cluster <- as.factor(GSEA_df$Cluster)
levels(GSEA_df$Cluster)
# Split the data frame by the 'Cluster' column
clusters <- split(GSEA_df, GSEA_df$Cluster)
# Write the list of data frames to an Excel file
write.xlsx(clusters, file = "GSEA_DC_and_B.xlsx", rowNames = T)
levels(twins@active.ident)
# Convert the x-axis variable to a factor if it is not already
gsea.list@compareClusterResult$Cluster <- factor(gsea.list@compareClusterResult$Cluster, levels = subsets)

# Dotplot with enriched pahtways per cell type
dotplot(gsea.list, color = 'NES', label_format = 35, showCategory = Inf) +
  scale_fill_viridis_c(option = 'H') +
  RotatedAxis() +
  theme(axis.text.y = element_text(size = 10)) # Adjust y-axis text size

ggsave(filename = "DotPlot_Pathways_Reactome.pdf", width = 20, height = 120, dpi = 320, units = "cm", device = "pdf")

# Define the pathways of interest (antigen presentation and immune memory)
antigen_immune <- c(
  #Antigen presentation
  "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
  "REACTOME_TCR_SIGNALING",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  #Memory formation
  "REACTOME_TCR_SIGNALING",
  "REACTOME_B_CELL_RECEPTOR_SIGNALING",
  "REACTOME_ADAPTIVE_IMMUNE_SYSTEM",
  #Metabolism and cell response to stress
  "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
  "REACTOME_SELENOEAMINO_ACID_METABOLISM",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  #Phagocytosis
  "REACTOME_DAP12_SIGNALING",
  "REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION",
  "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
  "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS",
  #Cytokine
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
  "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
  "REACTOME_TNFR1_INDUCED_NFKAPPAB_SIGNALING_PATHWAY",
  "REACTOME_IKK_COMPLEX_RECRUITMENT_MEDIATED_BY_RIP1",
  "REACTOME_TAK1_ACTIVATES_NFKB_BY_PHOSPHORYLATION_AND_ACTIVATION_OF_IKKS_COMPLEX",
  "REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION"
)

# Define the pathways of interest (antigen presentation and immune memory)
cell_cycle_translation <- c(
  #CellCycle
  "REACTOME_CYCLIN_D_ASSOCIATED_EVENTS_IN_G1",
  #Cyclin D is crucial for progression through the G1 phase of the cell cycle, which is important for cell proliferation.
  "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_ADDITIONAL_CELL_CYCLE_GENES",
  #TP53 plays a critical role in cell cycle regulation, particularly in response to stress, DNA damage, or immune activation.
  "REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_METHYLATION",
  #This pathway controls p53 activity via methylation, which is important for regulating immune cell proliferation.
  "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B",
  #Involved in mitotic exit and cell cycle progression, this pathway indicates how immune cells divide.
  "REACTOME_CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY",
  #This pathway marks the initiation of DNA replication, another important event in immune cell proliferation.
  "REACTOME_DEVELOPMENTAL_BIOLOGY",
  "REACTOME_MAP2K_AND_MAPK_ACTIVATION",
  "REACTOME_PROLONGED_ERK_ACTIVATION_EVENTS",
  "REACTOME_SIGNALING_BY_NOTCH4",
  "REACTOME_SIGNALING_BY_BRAF_AND_RAF1_FUSIONS",
  #RNA translation
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  #Involves the initiation of mRNA translation, a key step for producing proteins that are important for immune cell function.
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  #Represents the elongation phase of translation, where polypeptide chains are synthesized, crucial for immune activation.
  "REACTOME_RRNA_PROCESSING",
  #The production and modification of ribosomal RNA, necessary for ribosome function and protein synthesis.
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  #This pathway is involved in targeting proteins to the membrane during translation, critical for producing receptors and other membrane-bound proteins essential for immune function.
  "REACTOME_NON_SENSE_MEDIATED_DECAY_NMD",
  #This pathway degrades faulty mRNAs, maintaining the fidelity of protein synthesis during immune responses.
  "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S",
  "REACTOME_DNA_DAMAGE_RECOGNITION_IN_GG_NER",
  "REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL",
  "REACTOME_TRANSLATION",
  "REACTOME_METABOLISM_OF_RNA",
  "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_THE_AP_2_TFAP2_FAMILY_OF_TRANSCRIPTION_FACTORS"
)

celladhesion <- c(
  "REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL",
  "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
  "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", #This pathway describes how integrins mediate cell adhesion to the extracellular matrix and are crucial for migration.
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", #Covers how cells interact with and organize the extracellular matrix, which is key for cell adhesion and migration.
  "REACTOME_CELL_JUNCTION_ORGANIZATION", #Involves the formation and regulation of cell junctions, important for maintaining adhesion between cells.
  "REACTOME_FOCAL_ADHESION", #Details how focal adhesions form and facilitate cellular attachment to the substrate, essential for cell migration.
  "REACTOME_LEUKOCYTE_TRANSMIGRATION", #Explains the migration of leukocytes across endothelial barriers during inflammation.
  "REACTOME_CELL_ADHESION_MOLECULES_CAMs" #Focuses on molecules that mediate the adhesion of cells to each other and to the extracellular matrix.
)

# Subset the gsea.list for b cells to include only antigen presentation pathways
GSEA_filtered <- GSEA_df[GSEA_df$Description %in% antigen_immune, ]
GSEA_filtered <- GSEA_df[GSEA_df$Description %in% cell_cycle_translation, ]
GSEA_filtered <- GSEA_df[GSEA_df$Description %in% celladhesion, ]

gsea.list@compareClusterResult <- GSEA_filtered

# Create the DotPlot for b cells
dotplot(gsea.list, color = 'NES', label_format = 35, showCategory = Inf) +
  scale_fill_viridis_c(option = 'H') +
  RotatedAxis() +
  theme(axis.text.y = element_text(size = 10)) # Adjust y-axis text size

ggsave(filename = "DotPlot_Pathways_B_Reactome_CellAdhesion.pdf", width = 20, height = 15, dpi = 320, units = "cm", device = "pdf")
#####Visualize Genes of interest####
library(tidyverse)
library(pheatmap)
library(ggplot2)

bmem <- as.data.frame(degs.list[["B memory"]])
bmem$gene <- rownames(bmem)

#Select top genes based on rank (same ranking as GSEA)
top <- bmem %>%
  mutate(rank = -log10(p_val) * sign(avg_log2FC)) %>% 
  dplyr::select(rank, gene, p_val, avg_log2FC) %>% 
  rownames_to_column() %>%
  #filter(p_val <0.05 & abs(avg_log2FC) > 1) %>% 
  arrange(dplyr::desc(rank))

#Selectpathways of interest
gene_set_r <- msigdbr::msigdbr(category = 'C2', subcategory = 'CP:REACTOME') %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  unique()

paths <- GSEA_df %>% pull(ID) %>% unique()

gene_set_r_filtered <- gene_set_r %>%
  filter(gs_name %in% paths)

#Select Leading edge Genes from pathways that are also in top
#Use them to label volcano plot
selected_genes <- top %>%
  filter(gene %in% combined_genes)

genes_to_label <- selected_genes$gene

# Categorize the points
bmem$category <- ifelse(bmem$p_val < 0.1 & bmem$avg_log2FC < -1, "Significant < -1",
                        ifelse(bmem$p_val < 0.1 & bmem$avg_log2FC > 1, "Significant > 1",
                               "Not Significant"))

# Create the ggplot for the current subset
p <- ggplot(bmem, aes(avg_log2FC, -log10(p_val))) + 
  geom_point(aes(color = category), size = 0.5, alpha = 0.5) + 
  scale_color_manual(values = c("Significant < -1" = "blue", "Significant > 1" = "red", "Not Significant" = "grey")) +
  theme_bw() +
  ylab("-log10(p-value)") + 
  xlab("log2(Fold Change)") +
  geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), 
                  colour = "black", size = 3,  max.overlaps = 500) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey") +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),  # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 10),  # X-axis tick labels font and size
        panel.grid = element_blank()
  )
p

ggsave(filename = "BMemory_VolcanoPlot.pdf", plot = p, width = 10, height = 8, dpi = 300)

genesB <- selected_genes

cdc1 <- as.data.frame(degs.list[["cDC1"]])
cdc1$gene <- rownames(cdc1)

#Select top genes based on rank (same ranking as GSEA)
top <- cdc1 %>%
  mutate(rank = -log10(p_val) * sign(avg_log2FC)) %>% 
  dplyr::select(rank, gene, p_val, avg_log2FC) %>% 
  rownames_to_column() %>%
  #filter(p_val <0.05 & abs(avg_log2FC) > 1) %>% 
  arrange(dplyr::desc(rank))

#Select Leading edge Genes from pathways that are also in top
#Use them to label volcano plot
selected_genes <- top %>%
  filter(gene %in% combined_genes)

genes_to_label <- selected_genes$gene

# Categorize the points
cdc1$category <- ifelse(cdc1$p_val < 0.1 & cdc1$avg_log2FC < -1, "Significant < -1",
                       ifelse(cdc1$p_val < 0.1 & cdc1$avg_log2FC > 1, "Significant > 1",
                              "Not Significant"))

# Create the ggplot for the current subset
p <- ggplot(cdc1, aes(avg_log2FC, -log10(p_val))) + 
  geom_point(aes(color = category), size = 0.5, alpha = 0.5) + 
  scale_color_manual(values = c("Significant < -1" = "blue", "Significant > 1" = "red", "Not Significant" = "black")) +
  theme_bw() +
  ylab("-log10(p-value)") + 
  xlab("log2(Fold Change)") +
  geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), 
                  colour = "black", size = 3,  max.overlaps = 500) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey") +
  theme_bw() + 
  theme(legend.position = "none",
              axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
              axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),  # Y-axis title font and size
              axis.text.x = element_text(family = "Helvetica", size = 10),  # X-axis tick labels font and size
              panel.grid = element_blank()
              )
p

ggsave(filename = "cdc1_VolcanoPlot.pdf", plot = p, width = 10, height = 8, dpi = 300)

genescdc1 <-selected_genes

pdc <- as.data.frame(degs.list[["pDC"]])
pdc$gene <- rownames(pdc)

#Select top genes based on rank (same ranking as GSEA)
top <- pdc %>%
  mutate(rank = -log10(p_val) * sign(avg_log2FC)) %>% 
  dplyr::select(rank, gene, p_val, avg_log2FC) %>% 
  rownames_to_column() %>%
  #filter(p_val <0.05 & abs(avg_log2FC) > 1) %>% 
  arrange(dplyr::desc(rank))

#Select Leading edge Genes from pathways that are also in top
#Use them to label volcano plot
selected_genes <- top %>%
  filter(gene %in% combined_genes)

genes_to_label <- selected_genes$gene

# Categorize the points
pdc$category <- ifelse(pdc$p_val < 0.1 & pdc$avg_log2FC < -1, "Significant < -1",
                        ifelse(pdc$p_val < 0.1 & pdc$avg_log2FC > 1, "Significant > 1",
                               "Not Significant"))

# Create the ggplot for the current subset
p <- ggplot(pdc, aes(avg_log2FC, -log10(p_val))) + 
  geom_point(aes(color = category), size = 0.5, alpha = 0.5) + 
  scale_color_manual(values = c("Significant < -1" = "blue", "Significant > 1" = "red", "Not Significant" = "grey")) +
  theme_bw() +
  ylab("-log10(p-value)") + 
  xlab("log2(Fold Change)") +
  geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), 
                  colour = "black", size = 3,  max.overlaps = 500) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey") +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),  # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 10),  # X-axis tick labels font and size
        panel.grid = element_blank()
  )
p

ggsave(filename = "pDC_VolcanoPlot.pdf", plot = p, width = 10, height = 8, dpi = 300)

genespdc <-selected_genes

cells <- c("cDC1", "pDC", "B memory")
twins2 <- subset(twins.DC, subset = PBMC_l2_predicted.id %in% cells)

# Find all genes between the selected genes across the 4 cell types
genescdc1top <- c("RPS6", "RPS23", "RPS17", "RPSA", "EIF3H", "EIF3A", "RPS5", "RPS27", "RPS8", "RPS27A", "EIF3L", "EIF4A1", "RPS4X", "EIF4H", "RPS19", "RPS11", "RPS21", "PABPC1", "RPS15A", "RPS24", "FAU", "RPS2", 
                  "RBPJ", "UBC", "MAML2", "CREB1", "RPS27A", "NOTCH2", "APH1A", "EP300", 
                  "EEF1A1", "EEF2", "CAMKMT", "HSPA8")
genespdctop <- c("RPL41", "RPS12", "RPS10", "RPLP1", "RPL37A", "RPLP2", "RPS11", "RPL23", "RPL13", "RPS28", 
  "RPL34", "RPS19", "RPS25", "RPL19", "RPL39", "FAU", "RPS15A", "RPS21", "RPL28", "RPL10A", 
  "RPSA", "RPS6", "RPL35A", "RPL27A", "RPS3A", "RPL37", "RPL26", "RPL36A", "RPS15", "RPS4X", 
  "RPS14", "RPL13A", "RPS23", "TRAM1", "RPS27", "RPL32", "RPS18", "RPL23A", "RPL10", "UBA52", 
  "RPL17", "RPS27A", "RPLP0", "RPL9", "RPL22", "RPS24", "EEF1G", "RPL18A", "RPS3", "RPS16", 
  "RPS29", "RPS5", "RPL12", "RPL31", "EEF1A1", "RPS8", "RPL7A", "RPL24", "RPS20", "RPS9", 
  "RPL38", "RPL14", "RPS2", "PABPC1", "EIF3I", "EIF4B", "PSMB9", "PSMB4", "PSMC6", "CXCR4", 
  "LAMTOR3", "FLCN", "IMPACT", "ATF3", "GSR", "POLR2L", "MRPS21", "SEC61A1", "SSR3", "MRPL2", 
  "OXA1L", "NOL6", "MRM3", "RPL6", "WDR12", "RPS13", "RPL11", "RPL36", "RIOK3", "GLS", "OAZ1", 
  "GOT2", "AMT", "SAT1", "NDUFAB1"
)
genescd4top <- c(
  "RPS14", "RPL32", "UBA52", "RPL17", "RPS18", "RPS15A", "RPS12", "RPLP1", "RPS20", 
  "RPL21", "RPL12", "RPS7", "RPL9", "RPS24", "RPLP2", "RPL30", "RPL19", "RPL26", "RPL28", 
  "RPS8", "RPS21", "RPL13A", "RPL10", "RPS23", "RPL13", "RPS15", "RPL37A", "RPL27A", 
  "RPL18", "RPS10", "RPL35A", "RPL31", "RPL7", "RPL23", "RPS25", "RPS11", "RPS16", "RPL29", 
  "RPL18A", "RPL6", "RPS17", "RPS13", "RPS4X", "RPS26", "RPL4", "RPS29", "RPL39", "RPS6", 
  "RNPS1", "RPSA", "RPL37", "RPL34", "RPL5", "MAGOH", "RPS19", "RPL23A", "RPS3A", 
  "EIF3H", "PSMB4", "ELOB", "PSMD7", "MRPL10", "PTPN11", "RDX", "LRPPRC", "FKBP5", 
  "GSTP1", "COX6B1", "MAPKAPK5", "DIS3", "HSBP1", "DNAJB9", "RPS19BP1", "P4HB", "GOSR2", 
  "MT-CO3", "RPL10A", "CEBPG", "SESN1", "POLR2J", "MTAP", "EXOSC5", "NR3C1", "RPS15", 
  "RPS13", "RPS4X", "RPL34", "MRPS26", "MRPL34", "EIF2B2", "RPL5", "MRPL22", "MRPL28", 
  "RPL7A", "RPL27A", "RPL6", "RPL9", "RPS25", "RPS3A", "RPSA"
)
genesBtop <- c(
  "RPLP1", "RPL41", "RPL23A", "RPL35A", "RPL10", "RPS12", "RPL7A", "RPS18", "RPL17", "RPS27", 
  "RPS20", "RPS10", "RPL12", "RPL13", "RPL29", "RPL5", "RPS27A", "RPSA", "RPL32", "RPLP2", 
  "RPS15", "RPS25", "RPL9", "RPL10A", "MAGOH", "RPS6", "RPL26", "RPL28", "RPL27A", "RPS9", 
  "RPS15A", "RPL36A", "RPL34", "RPL30", "RPS14", "RPL15", "PSMA1", "RPS7", "RPL11", "RPS23", 
  "RPS27L", "RPL37", "RPL36AL", "RPS11", "RPS2", "RPL18", "RPL35", "RPL27", "RPS21", "RPL19", 
  "RPS19", "PSMA4", "RPS4X", "RPL14", "RPL22", "RPS26", "RPL23", "ROBO1", "RPS16", "FAU", "RPL7", 
  "PSMA7", "RPL8", "ETF1", "PSMB10", "RPS29", "RPLP0", "RPL21", "PSMD9", "RPL18A", "PSMB1", 
  "RPS8", "RPS24", "RPS13", "POLR2B", "GRSF1", "NDC1", "NUP153", "RPS3A", "RPL24", "RPS3", "RPL38", 
  "KPNA3", "RPS5", "FNIP1", "MTOR", "ATP6V1C1", "BMT2", "NOP10", "MTERF4", "WDR43", "DIS3", "PRORP", 
  "RRP9", "EXOSC9", "TSR3", "PES1", "MPHOSPH10", "MRPL32", "MRPS18A", "MRPL42", "MRPL54", "RARS2", 
  "MRPL45", "TARS2", "MRPL48", "MTFMT", "MRPL50", "MRPS2", "MRPL22", "MRPL44"
)

genes <- c(genescdc1top, genespdctop, genescd4top,genesBtop)
genes <- genes %>% unique()

# Extract the expression data for the selected genes from the Seurat object
expression_data <- twins2@assays$SCT@data[genes, ]

# Add metadata columns to group by cell type and vaccine condition
metadata <- twins2@meta.data %>%
  dplyr::select(PBMC_l2_predicted.id, vaccine)  # adjust metadata names as needed

# Add cell type and vaccine condition to the expression data
expression_data <- cbind(metadata, t(expression_data))

# Aggregate the expression by cell type and vaccine condition
aggregated_expression <- expression_data %>%
  group_by(PBMC_l2_predicted.id, vaccine) %>%
  summarise(across(all_of(genes), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

write_csv(aggregated_expression, "aggregated_expression.csv")

# Reshape the data for heatmap
heatmap_data <- aggregated_expression %>%
  pivot_longer(cols = all_of(genes),
               names_to = "gene",
               values_to = "expression") %>%
  unite("cell_vaccine", PBMC_l2_predicted.id, vaccine, sep = "_") %>%
  spread(key = gene, value = expression)

# Convert the heatmap data to a numeric matrix, excluding the first column (which contains row labels)
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$cell_vaccine

# Create row annotations for cell type and vaccine
row_annotations <- data.frame(
  celltype = gsub("_.*", "", rownames(heatmap_matrix)),  # Extract cell type from rowname
  vaccine = gsub(".*_", "", rownames(heatmap_matrix))    # Extract vaccine condition from rowname
)

# Convert to an annotation object for the heatmap
row_annotation_obj <- rowAnnotation(
  Celltype = row_annotations$celltype,
  Vaccine = row_annotations$vaccine,
  col = list(Celltype = c("cDC1" = "pink","pDC" = "#4DAF4A","B memory" = "#FFFF33"),
             Vaccine = c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")),
  gp = gpar(col = "white"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Helvetica"),  # Set legend font size and font
    labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica", face = "bold")
  )
)


# Plot the heatmap
scaled_matrix <- scale(heatmap_matrix)
plot <- Heatmap(scaled_matrix, 
                name = "z-score", 
                row_names_gp = gpar(fontsize = 12, family = "Helvetica",  face = "bold", color = "black"),
                column_names_gp = gpar(fontsize = 8, family = "Helvetica", color = "black"),
                cluster_rows = F, 
                cluster_columns = T,
                show_column_names = TRUE, 
                show_row_names = F, 
                col = colorRampPalette(c("blue", "white", "red"))(100),
                na_col = "white", 
                column_title = "",
                border = T,
                border_gp = gpar(col = "white"),
                rect_gp = gpar(col = "white", lwd = 0.5),
                right_annotation = row_annotation_obj,
                width = unit(45, "cm"),
                height = unit(2.5, "cm"),
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 12, fontfamily = "Helvetica", fontface = "bold"),  # Customize legend title font
                  labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica")  # Customize legend label font
                )
                )

plot
#####Average Expression NOTCH Genes#####
library(tidyplots)

notch_genes <- c(
  "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",  # NOTCH family receptors
  "DLL1", "DLL3", "DLL4",  # Delta-like ligands
  "JAG1", "JAG2",  # Jagged ligands
  "MAML1",  # Mastermind-like 1
  "RBPJ",  # Transcriptional regulator of NOTCH signaling
  "HES1", "HES5",  # HES transcription factors
  "HEY1", "HEY2", "HEYL",  # Hey family transcription factors
  "CSL"  # CBF1, also known as CSL (receptor binding)
)

ets_genes <- c(
  "ETS1", "ETS2", "ELK1", "ELK3", "ELK4", "ERG", "FLI1", "PU.1", "SPI1", "EGR1", "EGR2", "EGR3", "EGR4"  # ETS transcription factors
)

# Combine both lists
combined_genes <- c(notch_genes, ets_genes)

combined_genes <- c("RPS6", "RPS23", "RPS17", "RPSA", "EIF3H", "EIF3A", "RPS5", "RPS27", "RPS8", "RPS27A", "EIF3L", "EIF4A1", "RPS4X", "EIF4H", "RPS19", "RPS11", "RPS21", "PABPC1", "RPS15A", "RPS24", "FAU", "RPS2", 
                    "RBPJ", "UBC", "MAML2", "CREB1", "RPS27A", "NOTCH2", "APH1A", "EP300", 
                    "EEF1A1", "EEF2", "CAMKMT", "HSPA8", combined_genes) %>% unique()

cdc <- subset(twins.DC, subset = PBMC_l2_predicted.id == "cDC1")
DefaultAssay(cdc) <- "SCT"

genes <- as.data.frame(rownames(cdc@assays$SCT))
colnames(genes) <- "gene"
genes <- genes[genes$gene %in% combined_genes, ]

# Extract the expression matrix for the top DEGs
expr_data <- GetAssayData(cdc, slot = "data")[genes, , drop = FALSE]

# Add subject and cell type information to the metadata
metadata <- cdc@meta.data
metadata$subject <- gsub("_", "-", metadata$subject)
metadata$subject_celltype_vaccine <- paste(metadata$subject, metadata$PBMC_l2_predicted.id, metadata$vaccine,sep = "_")

# Calculate average expression per subject per cell type
aggregated_expr_data <- expr_data %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(subject_celltype_vaccine = metadata$subject_celltype_vaccine) %>%
  group_by(subject_celltype_vaccine) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("subject_celltype_vaccine") %>%
  t()

write.csv(aggregated_expr_data, "aggregated_expression.csv")

# Scale the data for heatmap visualization
scaled_expr_matrix <- t(scale(t(aggregated_expr_data)))

annotation_col <- data.frame(
  Subject = gsub("_.*", "", colnames(scaled_expr_matrix)),  # Extract subject (everything before the first '_')
  CellType = gsub("^[^_]+_([^_]+)_.*", "\\1", colnames(scaled_expr_matrix)),  # Extract cell type (everything after the first '_' but before the last)
  Vaccine = gsub(".*_", "", colnames(scaled_expr_matrix))  # Extract vaccine (everything after the last '_')
)

rownames(annotation_col) <- colnames(scaled_expr_matrix)

#Remove Subject
annotation_col <- annotation_col[, !colnames(annotation_col) %in% "Subject"]

# Order columns by cell type, then by Vaccine
annotation_col <- annotation_col %>%
  arrange(CellType, Vaccine)

# Reorder the columns in `scaled_expr_matrix` to match this sorted order
scaled_expr_matrix <- scaled_expr_matrix[, rownames(annotation_col)]

# Reshape the expression matrix into a long format for ggplot
aggregated_expr_data_long <- as.data.frame(aggregated_expr_data) %>%
  tibble::rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "subject_celltype_vaccine", values_to = "expression")

# Split the `subject_celltype_vaccine` column into separate columns for `subject`, `celltype`, and `vaccine`
aggregated_expr_data_long <- aggregated_expr_data_long %>%
  separate(subject_celltype_vaccine, into = c("subject", "celltype", "vaccine"), sep = "_")

aggregated_expr_data_long$subject <- as.factor(aggregated_expr_data_long$subject)
aggregated_expr_data_long$celltype <- as.factor(aggregated_expr_data_long$celltype)
aggregated_expr_data_long$vaccine <- as.factor(aggregated_expr_data_long$vaccine)

# Create the violin plot using ggplot2
library(ggplot2)

notch2_signaling_reactome <- c("RBPJ", "UBC", "MAML2", "CREB1", "RPS27A", "NOTCH2", "APH1A", "EP300")

# Filter the data to include only the selected genes
filtered_data <- aggregated_expr_data_long %>%
  filter(gene %in% notch2_signaling_reactome)

vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")
plot <- ggplot(filtered_data, aes(x = celltype, y = expression, fill = vaccine)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +  # Facet by gene for individual plots
  scale_fill_manual(values = vaccine_colors, labels = c("Cervarix", "Gardasil")) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = " ",
       x = "", y = "Expression", fill = "Vaccine")
plot

aggregated_expr_data_long %>% 
  tidyplot(x = vaccine, y = expression, color = vaccine) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "wilcoxon", p.adjust.method = "bonferroni", label = "p.adj.signif", hide_info = T) %>%
  add_test_pvalue(method = "wilcoxon", p.adjust.method = "bonferroni",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Average Expression") %>%
  adjust_legend_title("Vaccine") %>%
  split_plot(by = gene) 

filtered_data %>% 
  tidyplot(x = vaccine, y = expression, color = vaccine) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "wilcoxon", p.adjust.method = "bonferroni", label = "p.adj.signif", hide_info = T) %>%
  add_test_pvalue(method = "wilcoxon", p.adjust.method = "bonferroni",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Average Expression", family = "helvetica", fontsize = 8, color = "black") %>%
  adjust_y_axis(limits = c(0,NA), padding = c(0.1,0.15)) %>%
  adjust_legend_title("Vaccine",, family = "helvetica", fontsize = 8, color = "black", face = "bold") %>%
  adjust_font(family = "Helvetica", fontsize = 10, color = "black", face = "bold") %>%
  adjust_size(width = 100, height = 100, unit = "mm") %>%
  split_plot(by = gene) 

ggsave(filename = "Figure_BoxPlot_Notch2_in_cDC1.pdf", plot = plot, width = 15, height = 12, dpi = 320, units = "cm", device = "pdf")
#####Visualize Pathways of interest#####
library(tidyverse)
library(pheatmap)
library(ggplot2)

filtered_GSEA <- GSEA_df %>%
  dplyr::filter(Cluster %in% c("cDC1", "cDC2", "pDC")) %>%
  pull(ID)

filtered_GSEA_df <- GSEA_df %>%
  dplyr::filter(ID %in% filtered_GSEA)

# Reshape the data to get a matrix where rows are pathways and columns are cell types
heatmap_data <- filtered_GSEA_df %>%
  dplyr::select(Cluster, ID, NES) %>%
  spread(key = Cluster, value = NES)  # Create a wide format with cell types as columns

# Ensure pathways are in the rows and cell types are in the columns
heatmap_data2 <- heatmap_data[, -1]  # Remove the Description column (now row names)
rownames(heatmap_data2) <- heatmap_data$ID # Set pathway names as row names

# Optional: Scale the data (z-score) for better visualization
#heatmap_data_scaled <- t(scale(t(heatmap_data2)))

col_order <- c("cDC1", "cDC2", "pDC",
               "B naive", "B intermediate", "B memory", "Plasmablast")
colnames(heatmap_data2)

heatmap_data2 <- heatmap_data2[, col_order]
heatmap_data_matrix <- as.matrix(heatmap_data2)

#heatmap_data_scaled <- t(scale(t(heatmap_data2)))
#heatmap_data_scaled[is.nan(heatmap_data_scaled)] <- NA

write.xlsx(GSEA_df, "All_GSEA.xlsx")

# Define custom font
font_settings <- par(family = "Helvetica", cex = 1)

library(ComplexHeatmap)

row_annotation_colors <- ifelse(rownames(heatmap_data_matrix) == "SIGNALING BY NOTCH2", "Orange", "Black")

# Define colors for each cell type
celltype_colors <- c("cDC1" = "pink",  # Red
                     "cDC2" = "#377EB8",  # Blue
                     "pDC" = "#4DAF4A",   # Green
                     "B naive" = "#984EA3",  # Purple
                     "B intermediate" = "#FF7F00",  # Orange
                     "B memory" = "#FFFF33",  # Yellow
                     "Plasmablast" = "#A65628")  # Brown

column_annot <- HeatmapAnnotation(
  Cell_Type = factor(colnames(heatmap_data_matrix), levels = col_order),  # Ensure correct order
  col = list(Cell_Type = celltype_colors),
  show_legend = TRUE,
  annotation_legend_param = list(title = "Cell Type")
)

library(circlize)
# Define color scale centered at 0
min(heatmap_data_matrix, na.rm = TRUE)
max(heatmap_data_matrix, na.rm = TRUE)

col_fun <- colorRamp2(c(min(heatmap_data_matrix, na.rm = TRUE), 0, max(heatmap_data_matrix, na.rm = TRUE)), 
                      c("blue", "white", "red")) # Adjust colors as needed

# Generate heatmap
plot <- Heatmap(heatmap_data_matrix, 
                name = "NES Scores", 
                row_names_gp = gpar(fontsize = 8, family = "Helvetica", color = "black"),
                column_names_gp = gpar(fontsize = 12, family = "Helvetica", color = "black", rotation = 45),
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                show_column_names = FALSE, 
                show_row_names = TRUE, 
                col = col_fun,  # Apply color function here
                na_col = "white", 
                column_title = "",
                row_names_side = "left", 
                rect_gp = gpar(col = "white", lwd = 2),
                column_names_rot = 45,
                width = unit(2.5, "cm"), 
                height = unit(20, "cm"),
                top_annotation = column_annot,
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 12, fontfamily = "Helvetica", fontface = "bold"),  # Customize legend title font
                  labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica")  # Customize legend label font
                ))   

print(plot)

# Save the heatmap if desired
ggsave(filename = "Heatmap_Pathways_NES.pdf", plot = plot, width = 16, height = 5, dpi = 300)

#####Correlate GSEA pathways with HPV titers#####
unique_paths <- length(unique(GSEA_df$ID))
print(unique_paths)
unique_paths <- unique(GSEA_df$ID)

#Select pathways of interest
gene_set_r <- msigdbr::msigdbr(category = 'C2', subcategory = 'CP:REACTOME') %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  unique()

gene_set_r_filtered <- gene_set_r %>%
  filter(gs_name %in% unique_paths)

# For each unique pathway, extract the corresponding gene list
gene_lists <- lapply(unique_paths, function(pathway) {
  # Filter for genes corresponding to the current pathway
  gene_list <- gene_set_r_filtered %>%
    filter(gs_name == pathway) %>%
    pull(gene_symbol) # Extract the gene symbols for this pathway
  return(gene_list)
})

# To check the first few gene lists for the first pathway
names(gene_lists) <- unique_paths

# Now, you can use AddModuleScore with these gene lists
module_scores <- sapply(1:length(gene_lists), function(i) {
  gene_list <- gene_lists[[i]]
  # Filter out genes that are not in the Seurat object
  valid_genes <- gene_list[gene_list %in% rownames(twins.DC)]
  if(length(valid_genes) > 0) {
    # Use AddModuleScore with the valid genes only
    AddModuleScore(twins.DC, features = list(valid_genes), name = paste0("ModuleScore_", unique_paths[i]), ctrl = 5)
  } else {
    message(paste("No valid genes found for pathway:", unique_paths[i]))
  }
})

# Assuming your list of Seurat objects is called 'seurat_list' and the existing Seurat object is 'existing_seurat'
# Create the column names dynamically
column_names <- paste0("ModuleScore_", unique_paths, "1")

# Initialize an empty list to store the extracted dataframes
extracted_dataframes <- list()

# Loop through each name in column_names
for (i in 1:length(column_names)) {
  # Extract the corresponding dataframe
  column_name <- column_names[i]
  
  # Assuming module_scores is a list of dataframes, you can extract the desired column
  df_name <- paste0("module_scores[[", i, "]]$", column_name)
  extracted_dataframes[[column_name]] <- as.data.frame(eval(parse(text = df_name)))
}

# Combine all dataframes into one dataframe (assuming each is a column-wise dataframe)
combined_df <- do.call(cbind, extracted_dataframes)

# If you want to give specific names to each column (e.g., "p1", "p2", ..., "p42"), you can set the column names
colnames(combined_df) <- paste0("pa", 1:92) #in the order of "colnames"

# Assuming your Seurat object is called 'seurat_obj', add the combined dataframe to its meta.data
twins.DC@meta.data <- cbind(twins.DC@meta.data, combined_df)

rm(module_scores)

#Correlation of Module Scores with HPV Titers
module_scores <- as.data.frame(twins.DC@meta.data) %>%
  dplyr::select(PBMC_l2_predicted.id, 
                HPV6.y, HPV16.y, HPV18.y, HPV31.y, HPV33.y, HPV45.y, HPV52.y, HPV58.y,
                pa1, pa2, pa3, pa4, pa5, pa6, pa7, pa8, pa9, pa10, pa11, pa12, pa13, pa14, pa15, 
                pa16, pa17, pa18, pa19, pa20, pa21, pa22, pa23, pa24, pa25, pa26, pa27, pa28, 
                pa29, pa30, pa31, pa32, pa33, pa34, pa35, pa36, pa37, pa38, pa39, pa40, pa41, pa42,
                pa43, pa44, pa45, pa46, pa47, pa48, pa49, pa50, pa51, pa52, pa53, pa54, pa55, pa56, 
                pa57, pa58, pa59, pa60, pa61, pa62, pa63, pa64, pa65, pa66, pa67, pa68, pa69, pa70, 
                pa71, pa72, pa73, pa74, pa75, pa76, pa77, pa78, pa79, pa80, pa81, pa82, pa83, pa84, 
                pa85, pa86, pa87, pa88, pa89, pa90, pa91, pa92)

# Log10-transform the HPV titers
module_scores <- module_scores %>%
  mutate(across(starts_with("HPV"), ~ log10(. + 1)))

# List of HPV strains and Module Scores
hpv_strains <- c("HPV6.y", "HPV16.y", "HPV18.y", "HPV31.y", "HPV33.y", "HPV45.y", "HPV52.y", "HPV58.y")
module_names <- c("pa1", "pa2", "pa3", "pa4", "pa5", "pa6", "pa7", "pa8", "pa9", "pa10", "pa11", "pa12", "pa13", "pa14", "pa15", 
                  "pa16", "pa17", "pa18", "pa19", "pa20", "pa21", "pa22", "pa23", "pa24", "pa25", "pa26", "pa27", "pa28", 
                  "pa29", "pa30", "pa31", "pa32", "pa33", "pa34", "pa35", "pa36", "pa37", "pa38", "pa39", "pa40", "pa41", "pa42", 
                  "pa43", "pa44", "pa45", "pa46", "pa47", "pa48", "pa49", "pa50", "pa51", "pa52", "pa53", "pa54", "pa55", "pa56", 
                  "pa57", "pa58", "pa59", "pa60", "pa61", "pa62", "pa63", "pa64", "pa65", "pa66", "pa67", "pa68", "pa69", "pa70", 
                  "pa71", "pa72", "pa73", "pa74", "pa75", "pa76", "pa77", "pa78", "pa79", "pa80", "pa81", "pa82", "pa83", "pa84", 
                  "pa85", "pa86", "pa87", "pa88", "pa89", "pa90", "pa91", "pa92")
conditions <- c("CERVARIX", "GARDASIL")

# Create an empty list to store the results
correlation_results_list <- list()

library(broom)
# Loop through each condition, Module Score, and HPV strain
for (module in module_names) {
  for (strain in hpv_strains) {
    
    # Filter data by the current condition
    correlation_data <- module_scores %>%
      mutate(module_score = get(module), antibody_titer_log10 = get(strain))  # Use the pre-transformed titers
    
    # Perform correlation analysis
    correlation_results <- correlation_data %>%
      do(tidy(cor.test(.$module_score, .$antibody_titer_log10, method = "kendall"))) %>%
      mutate(Module = module, HPV_strain = strain)
    
    # Store the results in the list
    correlation_results_list[[paste(module, strain, sep = "_")]] <- correlation_results
  }
}

# Combine the results from all module score and strain combinations into a single data frame
final_correlation_results <- bind_rows(correlation_results_list)

# Remove rows with NA in the 'estimate' column
final_correlation_results <- final_correlation_results %>%
  filter(!is.na(estimate))

# View the combined correlation results
head(final_correlation_results)

column_names

# Heatmap of correlations for each Module Score and HPV strain
custom_labels <- c(
  "pa1" = "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "pa2" = "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  "pa3" = "REACTOME_SELENOAMINO_ACID_METABOLISM",
  "pa4" = "REACTOME_TRANSLATION",
  "pa5" = "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  "pa6" = "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
  "pa7" = "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
  "pa8" = "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
  "pa9" = "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "pa10" = "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
  "pa11" = "REACTOME_RRNA_PROCESSING",
  "pa12" = "REACTOME_INFLUENZA_INFECTION",
  "pa13" = "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
  "pa14" = "REACTOME_DEVELOPMENTAL_BIOLOGY",
  "pa15" = "REACTOME_NERVOUS_SYSTEM_DEVELOPMENT",
  "pa16" = "REACTOME_CELLULAR_RESPONSES_TO_STIMULI",
  "pa17" = "REACTOME_METABOLISM_OF_RNA",
  "pa18" = "REACTOME_INFECTIOUS_DISEASE",
  "pa19" = "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S",
  "pa20" = "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
  "pa21" = "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS",
  "pa22" = "REACTOME_FCERI_MEDIATED_CA_2_MOBILIZATION",
  "pa23" = "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
  "pa24" = "REACTOME_GPVI_MEDIATED_ACTIVATION_CASCADE",
  "pa25" = "REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES",
  "pa26" = "REACTOME_INOSITOL_PHOSPHATE_METABOLISM",
  "pa27" = "REACTOME_RAB_GERANYLGERANYLATION",
  "pa28" = "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "pa29" = "REACTOME_CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY",
  "pa30" = "REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "pa31" = "REACTOME_FORMATION_OF_TC_NER_PRE_INCISION_COMPLEX",
  "pa32" = "REACTOME_REGULATION_OF_RUNX3_EXPRESSION_AND_ACTIVITY",
  "pa33" = "REACTOME_COMPLEX_I_BIOGENESIS",
  "pa34" = "REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL",
  "pa35" = "REACTOME_ACTIVATION_OF_IRF3_IRF7_MEDIATED_BY_TBK1_IKK_EPSILON",
  "pa36" = "REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION",
  "pa37" = "REACTOME_INFECTION_WITH_MYCOBACTERIUM_TUBERCULOSIS",
  "pa38" = "REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY",
  "pa39" = "REACTOME_PINK1_PRKN_MEDIATED_MITOPHAGY",
  "pa40" = "REACTOME_RESPONSE_OF_MTB_TO_PHAGOCYTOSIS",
  "pa41" = "REACTOME_E3_UBIQUITIN_LIGASES_UBIQUITINATE_TARGET_PROTEINS",
  "pa42" = "REACTOME_RAS_PROCESSING",
  "pa43" = "REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS",
  "pa44" = "REACTOME_LATE_ENDOSOMAL_MICROAUTOPHAGY",
  "pa45" = "REACTOME_PROTEIN_UBIQUITINATION",
  "pa46" = "REACTOME_GLYCOGEN_METABOLISM",
  "pa47" = "REACTOME_IRON_UPTAKE_AND_TRANSPORT",
  "pa48" = "REACTOME_ASSEMBLY_OF_THE_HIV_VIRION",
  "pa49" = "REACTOME_BUDDING_AND_MATURATION_OF_HIV_VIRION",
  "pa50" = "REACTOME_SYNTHESIS_OF_ACTIVE_UBIQUITIN_ROLES_OF_E1_AND_E2_ENZYMES",
  "pa51" = "REACTOME_ONCOGENE_INDUCED_SENESCENCE",
  "pa52" = "REACTOME_TRANSLATION_OF_SARS_COV_2_STRUCTURAL_PROTEINS",
  "pa53" = "REACTOME_LISTERIA_MONOCYTOGENES_ENTRY_INTO_HOST_CELLS",
  "pa54" = "REACTOME_APC_CDC20_MEDIATED_DEGRADATION_OF_NEK2A",
  "pa55" = "REACTOME_NEGATIVE_REGULATORS_OF_DDX58_IFIH1_SIGNALING",
  "pa56" = "REACTOME_SPRY_REGULATION_OF_FGF_SIGNALING",
  "pa57" = "REACTOME_NUCLEOTIDE_SALVAGE",
  "pa58" = "REACTOME_MITOPHAGY",
  "pa59" = "REACTOME_ENDOSOMAL_SORTING_COMPLEX_REQUIRED_FOR_TRANSPORT_ESCRT",
  "pa60" = "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B",
  "pa61" = "REACTOME_DAP12_INTERACTIONS",
  "pa62" = "REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_METHYLATION",
  "pa63" = "REACTOME_TRAF6_MEDIATED_INDUCTION_OF_TAK1_COMPLEX_WITHIN_TLR4_COMPLEX",
  "pa64" = "REACTOME_SARS_COV_2_INFECTION",
  "pa65" = "REACTOME_NEGATIVE_REGULATION_OF_FGFR1_SIGNALING",
  "pa66" = "REACTOME_NEGATIVE_REGULATION_OF_FGFR2_SIGNALING",
  "pa67" = "REACTOME_NEGATIVE_REGULATION_OF_FGFR3_SIGNALING",
  "pa68" = "REACTOME_NEGATIVE_REGULATION_OF_FGFR4_SIGNALING",
  "pa69" = "REACTOME_SARS_COV_INFECTIONS",
  "pa70" = "REACTOME_TNFR1_INDUCED_NFKAPPAB_SIGNALING_PATHWAY",
  "pa71" = "REACTOME_SIGNALING_BY_EGFR_IN_CANCER",
  "pa72" = "REACTOME_DOWNREGULATION_OF_SMAD2_3_SMAD4_TRANSCRIPTIONAL_ACTIVITY",
  "pa73" = "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "pa74" = "REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS",
  "pa75" = "REACTOME_SIGNALING_BY_FGFR1",
  "pa76" = "REACTOME_SIGNALING_BY_FGFR4",
  "pa77" = "REACTOME_P75NTR_SIGNALS_VIA_NF_KB",
  "pa78" = "REACTOME_CIRCADIAN_CLOCK",
  "pa79" = "REACTOME_NEGATIVE_REGULATION_OF_MAPK_PATHWAY",
  "pa80" = "REACTOME_NEGATIVE_REGULATION_OF_MET_ACTIVITY",
  "pa81" = "REACTOME_CELL_DEATH_SIGNALLING_VIA_NRAGE_NRIF_AND_NADE",
  "pa82" = "REACTOME_CYCLIN_D_ASSOCIATED_EVENTS_IN_G1",
  "pa83" = "REACTOME_MAP3K8_TPL2_DEPENDENT_MAPK1_3_ACTIVATION",
  "pa84" = "REACTOME_SIGNALING_BY_FGFR2",
  "pa85" = "REACTOME_PEROXISOMAL_PROTEIN_IMPORT",
  "pa86" = "REACTOME_TERMINATION_OF_TRANSLESION_DNA_SYNTHESIS",
  "pa87" = "REACTOME_IRAK1_RECRUITS_IKK_COMPLEX",
  "pa88" = "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "pa89" = "REACTOME_TICAM1_RIP1_MEDIATED_IKK_COMPLEX_RECRUITMENT",
  "pa90" = "REACTOME_MYOGENESIS",
  "pa91" = "REACTOME_SIGNALING_BY_NOTCH2",
  "pa92" = "REACTOME_PROTEIN_METHYLATION"
  )

custom_labels2 <- c("HPV6.y" = "HPV6",
                    "HPV16.y" = "HPV16",
                    "HPV18.y" = "HPV18",
                    "HPV31.y" = "HPV31",
                    "HPV33.y" = "HPV33",
                    "HPV45.y" = "HPV45",
                    "HPV52.y" = "HPV52",
                    "HPV58.y" = "HPV58")
                    
final_correlation_results <- final_correlation_results %>%
  mutate(Module = custom_labels[Module], 
         HPV_strain = custom_labels2[HPV_strain])

write.xlsx(final_correlation_results, "Kendall_Correlation_Pathways_Titers.xlsx")

p <- ggplot(final_correlation_results, aes(x = HPV_strain, y = Module, fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", name="Spearman\nCorrelation") +
  scale_y_discrete(labels = custom_labels) +  # Apply custom labels to x-axis
  #facet_wrap(~ Condition, ncol = 2) +  # Separate heatmaps for each condition
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, size = 8),
        axis.text.y = element_text(hjust = 1, size = 8)) +
  labs(title = "",
       x = "HPV type",
       y = "Module Score")
p
ggsave(filename = "Correlation_ModuleScores_Log10(Titres)_ByVaccine.png", plot = p, width = 30, height = 10, dpi = 320, units = "cm", device = "png", bg = "white")

# Create a correlation matrix for clustering
correlation_matrix <- final_correlation_results %>%
  dplyr::select(Module, HPV_strain, estimate) %>%
  pivot_wider(names_from = HPV_strain, values_from = estimate) %>%
  column_to_rownames(var = "Module") %>%
  as.matrix()
rownames(correlation_matrix) <- custom_labels

# Create a data frame for row annotations based on the Vaccine column
# Subset annotation_row to match the correlation_matrix
annotation_row <- final_correlation_results %>%
  dplyr::select(Module, Vaccine) %>%
  filter(ID %in% rownames(correlation_matrix)) %>%
  distinct() %>%
  column_to_rownames(var = "ID")

rownames(correlation_matrix) <- gsub("^REACTOME_", "", rownames(correlation_matrix))
rownames(correlation_matrix)[rownames(correlation_matrix) == "RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS"] <- "RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS"

# Define color scale centered at 0
min(correlation_matrix, na.rm = TRUE)
max(correlation_matrix, na.rm = TRUE)

library(colorRamp2)
col_fun <- colorRamp2(c(min(correlation_matrix, na.rm = TRUE), 0, max(correlation_matrix, na.rm = TRUE)), 
                      c("blue", "white", "red")) # Adjust colors as needed

# Perform hierarchical clustering
hc_col <- hclust(dist(correlation_matrix), method = "complete")
hc_row <- hclust(dist(t(correlation_matrix)), method = "complete")

# Create a heatmap with clustering
pheatmap(t(correlation_matrix),
         #scale = "column",
         cluster_rows = hc_row, 
         cluster_cols = hc_col, 
         cutree_rows = 4,
         #cutree_cols = 4,
         display_numbers = T, 
         color = col_fun,
         main = " ",
         name = "Kendall's Tau rank",
         fontsize_row = 12,       # Adjust font size for row labels
         fontface_row = "bold",
         fontsize_col = 10,       # Adjust font size for column labels
         number_color = "black",  # Color of the numbers displayed
         #annotation_row = annotation_row,  # Add annotation for rows
         #annotation_colors = list(Vaccine = c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")),  # Specify colors for vaccines
         border_color = "white" ,
         na_col = "white",
         cellwidth = 20,
         cellheight = 20,
         angle_col = "45"
)

####NOTCH2 Module Score Correlation####
# Extract the metadata and create a dataframe
metadata <- twins.DC@meta.data
metadata <- metadata[,1:406]

# Aggregate by subject
aggregated_scores <- metadata %>%
  group_by(subject) %>%
  summarise(mean_module_score = mean(pa91, na.rm = TRUE))

merged_data <- merge(aggregated_scores, titer_data, by = "subject")
merged_data$HPV6_log <- log10(merged_data$HPV6)
merged_data$HPV16_log <- log10(merged_data$HPV16)
merged_data$HPV18_log <- log10(merged_data$HPV18)
merged_data$HPV31_log <- log10(merged_data$HPV31)
merged_data$HPV33_log <- log10(merged_data$HPV33)
merged_data$HPV45_log <- log10(merged_data$HPV45)
merged_data$HPV52_log <- log10(merged_data$HPV52)
merged_data$HPV58_log <- log10(merged_data$HPV58)
merged_data$Vaccine <- c("CERVARIX", "GARDASIL",
                             "CERVARIX", "GARDASIL",
                             "GARDASIL", "CERVARIX", 
                             "CERVARIX", "GARDASIL")

# Fit linear model
lm_model <- lm(HPV58_log ~ mean_module_score, data = merged_data)

# Extract R-squared and p-value
r_squared <- summary(lm_model)$r.squared
p_value <- summary(lm_model)$coefficients[2, 4]  # Extract p-value for slope

# Scatter plot of mean Module Score vs Titer
plothpv58 <- ggplot(merged_data, aes(x = mean_module_score, y = HPV58_log)) +
  geom_point(aes(color = Vaccine)) +  # Color only the points
  geom_smooth(method = "lm", se = TRUE, fill = "grey", alpha = 0.2, color = "black") +  # Single regression line
  labs(title = "",
       x = "NOTCH2 signaling Module Score",
       y = "log10(HPV58 Titer)") +
  theme_bw() +  # Keep the minimal theme
  scale_color_manual(values = vaccine_colors) +  # Custom colors for the Vaccine categories
  theme(
    #panel.border = element_rect(color = "black", size = 1),  # Add black border around the plot
    #plot.background = element_rect(fill = "white"),  # Add background color for the entire plot
    panel.background = element_rect(fill = "white"),  # Ensure the plot panel has a white background
    text = element_text(family = "Helvetica", size = 12),  # Base font settings
    axis.text.x = element_text(size = 12, colour = "black"),  # X-axis text color
    axis.text.y = element_text(size = 10, colour = "black"),  # Y-axis text color
    axis.title.x = element_text(size = 12, colour = "black"),  # X-axis title color
    axis.title.y = element_text(size = 12, colour = "black"),  # Y-axis title color
    legend.text = element_text(size = 12, colour = "black"),  # Legend text color
    panel.grid = element_blank()) +
  annotate("text", x = min(merged_data$mean_module_score), 
           y = max(merged_data$HPV58_log)+2, 
           label = paste0("R² = ", round(r_squared, 3), 
                          "\np = ", signif(p_value, 3)),
           hjust = 0, vjust = 1, size = 5, color = "black")

plothpv6 <- plothpv6
plothpv16 <- plothpv16 + scale_y_continuous(limits = c(0, NA))
plothpv18 <- plothpv18 + scale_y_continuous(limits = c(0, NA),n.breaks = 6)
plothpv31 <- plothpv31 + scale_y_continuous(limits = c(0, NA), n.breaks = 4)
plothpv33 <- plothpv33 + scale_y_continuous(limits = c(0, 6), n.breaks = 4)
plothpv45 <- plothpv45 + scale_y_continuous(limits = c(0, 6), n.breaks = 4)
plothpv52 <- plothpv52 + scale_y_continuous(limits = c(0, 6), n.breaks = 4)
plothpv58 <- plothpv58 + scale_y_continuous(limits = c(0, 6), n.breaks = 4)

plot <- (plothpv6 + plothpv16 + plothpv18) / 
  (plothpv31 + plothpv33 + plothpv45) / 
  (plothpv52 + plothpv58 + plot_spacer()) 

plot + plot_layout(guides = "collect")

plot

#####Breadth Score#####
titers <- as.data.frame(twins@meta.data) %>%
  dplyr::select(vaccine, PBMC_l2_predicted.id, HPV6.y, HPV16.y, HPV18.y, HPV31.y, HPV33.y, HPV45.y, HPV52.y, HPV58.y)

titers <- titers %>%
  mutate(
    breadth2 = rowSums(across(starts_with("HPV") & !ends_with("_flag")))
  )

summary(titers$breadth2)

titers <- titers %>%
  mutate(
    breadth2a = ifelse(breadth2 > median(breadth2), 1,0)
  )

titers$breadth2a <- as.factor(titers$breadth2a)
levels(titers$breadth2a)
summary(titers$breadth2a)

twins$breadth2a <- titers$breadth2a

Idents(twins) <- "PBMC_l2_predicted.id"
subsets <- levels(twins@active.ident)

twins$celltype.breadth <- paste(twins$PBMC_l2_predicted.id, twins$breadth2a, sep = "_")
Idents(twins) <- "celltype.breadth"
levels(twins@active.ident)
DefaultAssay(twins) <- "SCT"

degs.list <- list()
for(i in seq_along(subsets)) {
  degs.list[[i]] <- FindMarkers(twins,
                                ident.1 = paste(subsets[i], "1", sep = "_" ),
                                ident.2 = paste(subsets[i], "0", sep = "_" ), recorrect_umi = F)
}
names(degs.list) <- subsets

# Combine all data frames with an identifier and save as CSV
# Add gene column and remove rownames in each subset
for (i in seq_along(degs.list)) {
  degs.list[[i]]$gene <- rownames(degs.list[[i]])
  rownames(degs.list[[i]]) <- NULL
}
combined_degs <- bind_rows(degs.list, .id = "subset")
write.csv(combined_degs, file = "DEGs_Breadth.csv", row.names = FALSE)

# Initialize an empty list to store the plots
plots_list <- list()

# Loop through each subset in degs.list
for(i in seq_along(degs.list)) {
  # Extract the current data frame and reset row names
  degs.list[[i]]$gene <- rownames(degs.list[[i]])
  rownames(degs.list[[i]]) <- NULL
  
  # Categorize the points
  degs.list[[i]]$category <- ifelse(degs.list[[i]]$p_val < 0.05 & degs.list[[i]]$avg_log2FC < -1, "Significant < -1",
                                    ifelse(degs.list[[i]]$p_val < 0.05 & degs.list[[i]]$avg_log2FC > 1, "Significant > 1",
                                           "Not Significant"))
  
  # Create the ggplot for the current subset
  p <- ggplot(degs.list[[i]], aes(avg_log2FC, -log10(p_val))) + 
    geom_point(aes(color = category), size = 0.5, alpha = 0.5) + 
    scale_color_manual(values = c("Significant < -1" = "blue", "Significant > 1" = "red", "Not Significant" = "grey")) +
    theme_bw() +
    ylab("-log10(p-value)") + 
    xlab("log2(Fold Change)") +
    geom_text_repel(aes(label = ifelse(p_val < 0.01, gene, "")),
                    colour = "black", size = 3, max.overlaps = 20) +
    geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme(legend.position = "none",
          axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
          axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),  # Y-axis title font and size
          axis.text.x = element_text(family = "Helvetica", size = 10),  # X-axis tick labels font and size
          panel.grid = element_blank(),
          plot.title = element_text(family = "Helvetica", size = 14, face = "bold", hjust = 0)
          ) +
    ggtitle(names(degs.list)[i])

  # Store the plot in the list
  plots_list[[names(degs.list)[i]]] <- p
}

volcano_dc <- (plots_list[[24]] /plots_list[[21]] / plots_list[[8]])

volcano_b <- (plots_list[[1]] / plots_list[[14]] / plots_list[[18]] / plots_list[[27]])

# Ranked list of genes
ranks.list <- list()
for (i in seq_along(degs.list)) {
  ranks.list[[i]] <- degs.list[[i]] %>%
    mutate(p_val = ifelse(p_val == 0, .Machine$double.xmin*row_number(), p_val)) %>% 
    mutate(rank = -log10(p_val)*sign(avg_log2FC)) %>% 
    dplyr::select(rank) %>% 
    rownames_to_column() %>% 
    arrange(dplyr::desc(rank)) %>% 
    deframe() 
}
names(ranks.list) <- names(degs.list)

gene_set_r <- msigdbr::msigdbr(category = 'C2', subcategory = 'CP:REACTOME') %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  unique()

# Gene set enrichment analysis
gsea.list <- compareCluster(ranks.list, fun = 'GSEA', 
                            TERM2GENE = gene_set_r, eps = 0)

#Save
GSEA_df <- gsea.list@compareClusterResult
GSEA_df$Cluster <- as.factor(GSEA_df$Cluster)
levels(GSEA_df$Cluster)
# Split the data frame by the 'Cluster' column
clusters <- split(GSEA_df, GSEA_df$Cluster)
# Write the list of data frames to an Excel file
write.xlsx(clusters, file = "GSEA_High.vs.Low.Breadth_Reactome.xlsx", rowNames = T)

# Convert the x-axis variable to a factor if it is not already
gsea.list@compareClusterResult$Cluster <- factor(gsea.list@compareClusterResult$Cluster, levels = subsets)

# Dotplot with enriched pahtways per cell type
dotplot(gsea.list, color = 'NES', label_format = 35, showCategory = Inf) +
  scale_fill_viridis_c(option = 'H') +
  RotatedAxis() +
  theme(axis.text.y = element_text(size = 10)) # Adjust y-axis text size

ggsave(filename = "DotPlot_Pathways_High.vs.Low.Breadth.pdf", width = 20, height = 120, dpi = 320, units = "cm", device = "pdf")

Significant_Pathways <- c(
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  "REACTOME_SELENOAMINO_ACID_METABOLISM",
  "REACTOME_TRANSLATION",
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
  "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
  "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
  "REACTOME_RRNA_PROCESSING",
  "REACTOME_INFLUENZA_INFECTION",
  "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
  "REACTOME_DEVELOPMENTAL_BIOLOGY",
  "REACTOME_NERVOUS_SYSTEM_DEVELOPMENT",
  "REACTOME_CELLULAR_RESPONSES_TO_STIMULI",
  "REACTOME_METABOLISM_OF_RNA",
  "REACTOME_INFECTIOUS_DISEASE",
  "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S",
  "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS",
  "REACTOME_FCERI_MEDIATED_CA_2_MOBILIZATION",
  "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
  "REACTOME_GPVI_MEDIATED_ACTIVATION_CASCADE",
  "REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES",
  "REACTOME_INOSITOL_PHOSPHATE_METABOLISM",
  "REACTOME_RAB_GERANYLGERANYLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY",
  "REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_FORMATION_OF_TC_NER_PRE_INCISION_COMPLEX",
  "REACTOME_REGULATION_OF_RUNX3_EXPRESSION_AND_ACTIVITY",
  "REACTOME_COMPLEX_I_BIOGENESIS",
  "REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL",
  "REACTOME_ACTIVATION_OF_IRF3_IRF7_MEDIATED_BY_TBK1_IKK_EPSILON",
  "REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION",
  "REACTOME_INFECTION_WITH_MYCOBACTERIUM_TUBERCULOSIS",
  "REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY",
  "REACTOME_PINK1_PRKN_MEDIATED_MITOPHAGY",
  "REACTOME_RESPONSE_OF_MTB_TO_PHAGOCYTOSIS",
  "REACTOME_E3_UBIQUITIN_LIGASES_UBIQUITINATE_TARGET_PROTEINS",
  "REACTOME_RAS_PROCESSING",
  "REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS",
  "REACTOME_LATE_ENDOSOMAL_MICROAUTOPHAGY",
  "REACTOME_PROTEIN_UBIQUITINATION",
  "REACTOME_GLYCOGEN_METABOLISM",
  "REACTOME_IRON_UPTAKE_AND_TRANSPORT",
  "REACTOME_ASSEMBLY_OF_THE_HIV_VIRION",
  "REACTOME_BUDDING_AND_MATURATION_OF_HIV_VIRION",
  "REACTOME_SYNTHESIS_OF_ACTIVE_UBIQUITIN_ROLES_OF_E1_AND_E2_ENZYMES",
  "REACTOME_ONCOGENE_INDUCED_SENESCENCE",
  "REACTOME_TRANSLATION_OF_SARS_COV_2_STRUCTURAL_PROTEINS",
  "REACTOME_LISTERIA_MONOCYTOGENES_ENTRY_INTO_HOST_CELLS",
  "REACTOME_APC_CDC20_MEDIATED_DEGRADATION_OF_NEK2A",
  "REACTOME_NEGATIVE_REGULATORS_OF_DDX58_IFIH1_SIGNALING",
  "REACTOME_SPRY_REGULATION_OF_FGF_SIGNALING",
  "REACTOME_NUCLEOTIDE_SALVAGE",
  "REACTOME_MITOPHAGY",
  "REACTOME_ENDOSOMAL_SORTING_COMPLEX_REQUIRED_FOR_TRANSPORT_ESCRT",
  "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B",
  "REACTOME_DAP12_INTERACTIONS",
  "REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_METHYLATION",
  "REACTOME_TRAF6_MEDIATED_INDUCTION_OF_TAK1_COMPLEX_WITHIN_TLR4_COMPLEX",
  "REACTOME_SARS_COV_2_INFECTION",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR1_SIGNALING",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR2_SIGNALING",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR3_SIGNALING",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR4_SIGNALING",
  "REACTOME_SARS_COV_INFECTIONS",
  "REACTOME_TNFR1_INDUCED_NFKAPPAB_SIGNALING_PATHWAY",
  "REACTOME_SIGNALING_BY_EGFR_IN_CANCER",
  "REACTOME_DOWNREGULATION_OF_SMAD2_3_SMAD4_TRANSCRIPTIONAL_ACTIVITY",
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS",
  "REACTOME_SIGNALING_BY_FGFR1",
  "REACTOME_SIGNALING_BY_FGFR4",
  "REACTOME_P75NTR_SIGNALS_VIA_NF_KB",
  "REACTOME_CIRCADIAN_CLOCK",
  "REACTOME_NEGATIVE_REGULATION_OF_MAPK_PATHWAY",
  "REACTOME_NEGATIVE_REGULATION_OF_MET_ACTIVITY",
  "REACTOME_CELL_DEATH_SIGNALLING_VIA_NRAGE_NRIF_AND_NADE",
  "REACTOME_CYCLIN_D_ASSOCIATED_EVENTS_IN_G1",
  "REACTOME_MAP3K8_TPL2_DEPENDENT_MAPK1_3_ACTIVATION",
  "REACTOME_SIGNALING_BY_FGFR2",
  "REACTOME_PEROXISOMAL_PROTEIN_IMPORT",
  "REACTOME_TERMINATION_OF_TRANSLESION_DNA_SYNTHESIS",
  "REACTOME_IRAK1_RECRUITS_IKK_COMPLEX",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_TICAM1_RIP1_MEDIATED_IKK_COMPLEX_RECRUITMENT",
  "REACTOME_MYOGENESIS",
  "REACTOME_SIGNALING_BY_NOTCH2",
  "REACTOME_PROTEIN_METHYLATION"
)

# Subset the gsea.list to include only pathways significant between vaccines
#then filter for Cell specific pathways
GSEA_filtered <- GSEA_df[GSEA_df$Description %in% Significant_Pathways, ]

DC_Pathways <- GSEA_filtered %>%
  dplyr::filter(Cluster %in% c("cDC1", "cDC2", "pDC", "B naive", "B intermediate", "B memory", "Plasmablast")) %>%
  pull(ID)

GSEA_filtered <- GSEA_filtered %>%
  dplyr::filter(ID %in% DC_Pathways)

gsea.list@compareClusterResult <- GSEA_filtered

selected_clusters <- c("cDC1", "cDC2", "pDC", "B naive", "B intermediate", "B memory", "Plasmablast") # Replace with actual names
filtered_gsea <- gsea.list[gsea.list@compareClusterResult$Cluster %in% selected_clusters, ]
cluster_frq_subset$Thelp <- factor(cluster_frq_subset$Thelp, levels = selected_celltypes)
filtered_gsea$Cluster <- factor(filtered_gsea$Cluster, levels = selected_clusters)
gsea.list@compareClusterResult <- filtered_gsea

# Create the DotPlot
dotplot(gsea.list, color = 'NES', label_format = 35, showCategory = Inf) +
  scale_fill_viridis_c(option = 'H') +
  RotatedAxis() +
  theme(axis.title.x = element_blank(),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 10),  # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis tick labels font and size
        #panel.grid = element_blank()
  )

ggsave(filename = "DotPlot_Breadth_Reactome_CellPathways.pdf", width = 25, height = 75, dpi = 320, units = "cm", device = "pdf")

#####CellChat#####
Idents(twins.DC) <- "PBMC_l2_predicted.id"
levels(twins.DC@active.ident)

#Select CERVARIX or GARDASIL
query_CER <- subset(x=twins.DC, subset = vaccine == "CERVARIX")
dim(query_CER)
query_GAR <- subset(x=twins.DC, subset = vaccine == "GARDASIL")
dim(query_GAR)

#CERVARIX
data.input <- query_CER[["SCT"]]$data 
labels <- Idents(query_CER)
meta <- data.frame(labels = labels, row.names = names(labels))
cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellChat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 3 * 1024^3)  # 3 GiB
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 576

cellChat <- computeCommunProb(cellChat, type = "triMean")

cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)

groupSize <- as.numeric(table(cellChat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

celltype_colors <- c("cDC1" = "pink",  # Red
                     "cDC2" = "#377EB8",  # Blue
                     "pDC" = "#4DAF4A",   # Green
                     "B naive" = "#984EA3",  # Purple
                     "B intermediate" = "#FF7F00",  # Orange
                     "B memory" = "#FFFF33",  # Yellow
                     "Plasmablast" = "#A65628")  # Brown

selected_clusters <- c("cDC1", "cDC2", "pDC", "B naive", "B intermediate", "B memory", "Plasmablast") # Replace with actual names
cellChat@idents <- factor(cellChat@idents, levels = selected_clusters, ordered = TRUE)

debug(netAnalysis_signalingRole_heatmap)
undebug(netAnalysis_signalingRole_heatmap)
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing", color.use = celltype_colors,
                                         width = 5, height = 10, font.size = 12, font.size.title = 16,
                                         title = "Cervarix")
ht1
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming", height = 15)
ht2
                                                                                   
saveRDS(cellChat, file = "SeuratObjects/cellChat_query_Cervarix.rds")
cellchat.CER <- cellChat

#Gardasil
data.input <- query_GAR[["SCT"]]$data 
labels <- Idents(query_GAR)
meta <- data.frame(labels = labels, row.names = names(labels))

cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human 
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellChat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 3 * 1024^3)  # 3 GiB
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 593

cellChat <- computeCommunProb(cellChat, type = "triMean")

cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)

groupSize <- as.numeric(table(cellChat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellChat@idents <- factor(cellChat@idents, levels = selected_clusters, ordered = TRUE)
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing", color.use = celltype_colors,
                                         width = 5, height = 10, font.size = 12, font.size.title = 16,
                                         title = "Gardasil-4")
ht1
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming", height = 15)
ht2

saveRDS(cellChat, file = "SeuratObjects/cellChat_query_Gardasil.rds")
cellchat.GAR <- cellChat

#Compare 
cellchat.CER <- readRDS("cellChat_query_Cervarix.rds")
cellchat.GAR <- readRDS("cellChat_query_Gardasil.rds")

object.list <- list(GAR = cellchat.GAR, CER = cellchat.CER)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

save(object.list, file = "SeuratObjects/cellchat_object.list.RData")
save(cellchat, file = "SeuratObjects/cellchat_merged.RData")

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",
                           color.use = c("1" = "#007A73", "2" = "#E4AD2E"), color.alpha = 0.5,
                           digits = 2, size.text = 14, width = 0.5,
                           xlabel = "Vaccine", remove.xtick = F) +
  theme_bw() + 
  theme(legend.position = "none",
    axis.text.x.top = element_blank(),
    axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
    axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"), 
    axis.text.y = element_text(family = "Helvetica", size = 10, color = "black"), # Y-axis title font and size
    axis.text.x = element_text(family = "Helvetica", size = 10, color = "black"),  # X-axis tick labels font and size
    panel.grid = element_blank(),
    legend.text = element_text(family = "Helvetica", size = 10, color = "black"),
    legend.title = element_text(family = "Helvetica", size = 12, face = "bold")
  )
  
gg2

dc <- c("cDC1", "cDC2", "pDC")
b <- c("B naive", "B intermediate", "B memory","Plasmablast")

gg1 <- netVisual_heatmap(cellchat, row.show = dc, col.show = b, cluster.rows = T, cluster.cols = T)
gg1
gg2 <- netVisual_heatmap(cellchat, measure = "weight",row.show = dc, col.show = b, 
                         cluster.rows = T, cluster.cols = T, color.use = celltype_colors,
                         #width = 2, height = 2, font.size = 12, font.size.title = 16
                         ) 
gg2

#####NicheNet#####
library(nichenetr)

twins.DC <- SetIdent(twins.DC, value = "PBMC_l2_predicted.id")
levels(twins.DC@active.ident)

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network <- lr_network %>% distinct(from, to)

#Sender-agnostic
receiver = "B memory"
expressed_genes_receiver <- get_expressed_genes(receiver, twins.DC, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#Sender-focused
sender_celltypes <- c("cDC1", "pDC") 

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, twins.DC, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 8102
length(potential_ligands)
## [1] 472
length(potential_ligands_focused)
## [1] 89

#Gene set of interest
condition_oi <-  "CERVARIX"
condition_reference <- "GARDASIL"

twins_receiver <- subset(twins.DC, idents = receiver)

DE_table_receiver <-  FindMarkers(object = twins_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "vaccine",
                                  min.pct = 0.05, recorrect_umi = FALSE) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define the background
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 5144
length(geneset_oi)
## [1] 465

#Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(dplyr::desc(aupr_corrected)))
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

#Sender-focused approach
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme_bw() + 
  theme(#legend.position = "none",
        axis.text.x.top = element_blank(),
        axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"), 
        axis.text.y = element_text(family = "Helvetica", size = 10, color = "black"), # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 10, color = "black"),  # X-axis tick labels font and size
        panel.grid = element_blank(),
        legend.text = element_text(family = "Helvetica", size = 10, color = "black"),
        legend.title = element_text(family = "Helvetica", size = 12, face = "bold")
  )+
  coord_fixed(ratio = 0.5)

p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 40) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple") +
  theme_bw() + 
  theme(#legend.position = "none",
        #axis.text.x.top = element_blank(),
        axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"), 
        axis.text.y = element_text(family = "Helvetica", size = 10, color = "black"), # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 10, color = "black"),  # X-axis tick labels font and size
        panel.grid = element_blank(),
        legend.text = element_text(family = "Helvetica", size = 10, color = "black"),
        legend.title = element_text(family = "Helvetica", size = 12, face = "bold")
  )

p_ligand_target


#####NicheNet NOTCH focused#####
library(nichenetr)
library(Seurat)
library(tidyverse)

twins <- SetIdent(twins, value = "PBMC_l2_predicted.id")
levels(twins@active.ident)

genes <- c("RPS27A", "ITCH", "PSEN2", "NOTCH1", "CNTN1", "MIB1", "JAG1", "NUMB", "APH1A", "DTX2", "PSENEN", "PSEN1", 
           "RPS27A", "PSMB5", "PSME3", "PSMB2", "PSMA5", "PSMD3", "PSMD14", "PSMA2", "PSMD10", "PSMA3", "PSMB6", 
           "PSMD8", "PSMC2", "PSMB4", "PSMD4", "PSMD9", "PSMB7", "PSMD6", "PSMB3", "PSMD12", "PSMF1", "PSMC1", 
           "PSMB1", "FBXW7", "PSMD1", "PSMD13", "PSMC6", "PSMA7", "PSMC5", "PSMD7", "NOTCH4", "PSMA1", "PSMC4", 
           "PSMA4", "RBX1", "PSMA6", "HDAC9", "HDAC8", "NOTCH1", "NCOR2", "HDAC7", "TBL1XR1", "CREBBP", "NOTCH4", 
           "HDAC2", "TBL1X", "NBEA", "TLE3", "RPS27A", "TLE1", "CCNC", "HDAC9", "HDAC8", "MYC", "HEY1", "NOTCH1", 
           "NCOR2", "RPS27A", "PSEN2", "CNTN1", "MIB1", "JAG1", "APH1A", "PSENEN", "PSEN1", "RPS27A", "PSEN2", "MIB1", 
           "JAG1", "EGF", "APH1A", "WWP2", "PSENEN", "PSEN1", "DLGAP5", "PBX1", "HEY1", "PTCRA", "NOTCH1", "WWC1", 
           "FABP7", "EP300", "CREBBP", "HEY1", "NOTCH1", "ACTA2", "EP300", "CREBBP", "NOTCH4", "SEL1L", "TFDP1", 
           "JUN", "ATP2A2", "TMED2", "NOTCH1", "B4GALT1", "PRKCI", "RFNG", "RAB6A", "ATP2A3", "MFNG", "AGO3", "TP53", 
           "ST3GAL3", "E2F1", "EP300", "CREBBP", "ST3GAL6", "LFNG", "NOTCH4", "POFUT1", "RUNX1", "SEL1L", "ATP2A2", 
           "TMED2", "NOTCH1", "B4GALT1", "RFNG", "RAB6A", "ATP2A3", "MFNG", "ST3GAL3", "ST3GAL6", "LFNG", "NOTCH4", 
           "SEL1L", "DLGAP5", "NBEA", "TLE3", "RPS27A", "TLE1", "CCNC", "PSMB5", "HDAC9", "PSME3", "PBX1", "ITCH", 
           "TFDP1", "HDAC8", "JUN", "ATP2A2", "PSMB2", "MYC", "HEY1", "TMED2", "PSMA5", "PTCRA", "PSEN2", "NOTCH1", 
           "CNTN1", "MIB1", "PSMD3", "PSMD14", "B4GALT1", "PRKCI", "PSMA2", "FCER2", "JAG1", "PSMD10", "NCOR2", "PSMA3", 
           "WWC1", "RFNG", "RAB6A", "PSMB6", "NUMB", "ATP2A3", "MFNG", "PSMD8", "ACTA2", "EGF", "APH1A", "PSMC2", 
           "PSMB4", "PSMD4", "PSMD9", "WWP2", "PSMB7", "PSMD6", "DTX2", "PSENEN", "HDAC7", "PSMB3", "PSMD12", "PSMF1", 
           "PSMC1", "FABP7", "PSMB1", "AGO3", "FBXW7", "TBL1XR1", "PSEN1", "TP53", "PSMD1", "ST3GAL3", "PSMD13", 
           "PSMC6", "E2F1", "EP300", "CREBBP", "PSMA7", "PSMC5", "DTX1", "PSMD7", "ST3GAL6", "LFNG", "NOTCH4", "PSMA1", 
           "POFUT1", "PSMC4", "RUNX1", "HDAC2", "PSMA4", "RBX1", "NBEA", "TLE3", "RPS27A", "TLE1", "CCNC", "HDAC9", 
           "ITCH", "HDAC8", "MYC", "HEY1", "PSEN2", "NOTCH1", "CNTN1", "MIB1", "JAG1", "NCOR2", "NUMB", "APH1A", "DTX2", 
           "PSENEN", "HDAC7", "FBXW7", "TBL1XR1", "PSEN1", "EP300", "CREBBP", "RPS27A", "NOTCH1", "MIB1", "JAG1", 
           "RPS27A", "CCNC", "HDAC9", "HDAC8", "MYC", "HEY1", "PSEN2", "NOTCH1", "MIB1", "JAG1", "NCOR2", "APH1A", 
           "PSENEN", "HDAC7", "FBXW7", "TBL1XR1", "PSEN1", "EP300", "CREBBP", "RPS27A", "PSEN2", "CNTN1", "MIB1", 
           "FCER2", "JAG1", "APH1A", "PSENEN", "PSEN1", "EP300", "DLGAP5", "RPS27A", "PBX1", "HEY1", "PTCRA", "PSEN2", 
           "NOTCH1", "MIB1", "JAG1", "WWC1", "EGF", "APH1A", "WWP2", "PSENEN", "FABP7", "PSEN1", "EP300", "CREBBP", 
           "RPS27A", "PSMB5", "PSME3", "PSMB2", "HEY1", "PSMA5", "PSEN2", "NOTCH1", "PSMD3", "PSMD14", "PSMA2", "JAG1", 
           "PSMD10", "PSMA3", "PSMB6", "PSMD8", "ACTA2", "APH1A", "PSMC2", "PSMB4", "PSMD4", "PSMD9", "PSMB7", "PSMD6", 
           "PSENEN", "PSMB3", "PSMD12", "PSMF1", "PSMC1", "PSMB1", "FBXW7", "PSEN1", "PSMD1", "PSMD13", "PSMC6", "EP300", 
           "CREBBP", "PSMA7", "PSMC5", "PSMD7", "NOTCH4", "PSMA1", "PSMC4", "PSMA4", "RBX1", "PSMA6", "NOTCH2" ,"NOTCH3", "CD40") %>% unique()

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network <- subset(lr_network, from %in% genes | to %in% genes)
#lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction" & database == "kegg")
#lr_network <- subset(lr_network, from %in% c("IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNB1","IFNE","IFNG","IFNK","IFNL1","IFNL2","IFNL3","IFNW1","IL10","IL11","IL12A","IL12B","IL13","IL15","IL16","IL17A","IL17B","IL17C","IL17F","IL18","IL19","IL1A","IL1B","IL1F10","IL1RN","IL2","IL20","IL21","IL22","IL23A","IL24","IL25","IL27","IL3","IL31","IL33","IL34","IL36A","IL36B","IL36G","IL4","IL5","IL6","IL7","IL9","TGFA","TGFB1","TGFB2","TGFB3","TNF"))

ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
weighted_networks <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to"))

#Sender-agnostic
receiver = "B memory" 
expressed_genes_receiver <- get_expressed_genes(receiver, twins.DC, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#Sender-focused
sender_celltypes <- c("cDC1", "pDC")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, twins.DC, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 8102
length(potential_ligands)
## [1] 17
length(potential_ligands_focused)
## [1]6

#Gene set of interest
condition_oi <-  "CERVARIX"
condition_reference <- "GARDASIL"

twins_receiver <- subset(twins.DC, idents = receiver)

DE_table_receiver <-  FindMarkers(object = twins_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "vaccine",
                                  min.pct = 0.05, recorrect_umi = FALSE) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define the background
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 5144 (best between 5000 - 10000 and sufficiently larger than geneset_oi)
length(geneset_oi)
## [1] 465

#Perform NicheNet ligand activity analysis
#Sender-agnostic
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(dplyr::desc(aupr_corrected)))

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange", bins = 100)  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected) %>% base::as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

#Infer target genes and receptors of top-ranked ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 10) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)
## [1] 555
head(active_ligand_target_links_df)
## # A tibble: 6 × 3
##   ligand target  weight
##   <chr>  <chr>    <dbl>
## 1 Ifna1  Ddx58    0.247
## 2 Ifna1  Eif2ak2  0.246
## 3 Ifna1  Gbp2     0.192
## 4 Ifna1  Gbp7     0.195
## 5 Ifna1  H2-D1    0.206
## 6 Ifna1  H2-K1    0.206

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

nrow(active_ligand_target_links)
## [1] 543
head(active_ligand_target_links)
##        Ifna13 Ifna2 Ifna6 Ifna15 Ifna7 Ifna5 Ifnab Ifna9 Ifna11 Ifna12 Ifna16 Ifna4 Ifna14 Ptprc        Tnf      Il36g       Il10      Il21        Osm
## Irf1        0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.27692301 0.07400782 0.07722567 0.1342983 0.16962803
## Ddx60       0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.11281871 0.00000000 0.05478472 0.0000000 0.08116101
## Parp14      0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.07003101 0.06895448 0.00000000 0.0000000 0.08011593
## Ddx58       0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.24433255 0.06891134 0.00000000 0.0000000 0.08862524
## Parp12      0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.19298997 0.06687691 0.05621734 0.0000000 0.07252823
## Tap1        0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.25076038 0.07099514 0.00000000 0.0280451 0.15842935
##             Il27     Ifna1      Ifnb1      Ifng       Ifnk       Ifne      Lrtm2      Ifnl3       Ebi3      Ifnl2
## Irf1   0.3393635 0.2493108 0.25825704 0.2864755 0.04430047 0.04063913 0.03483987 0.10021371 0.11317944 0.07408235
## Ddx60  0.1596771 0.1218225 0.13569911 0.1171453 0.02629924 0.02771157 0.00000000 0.08217529 0.04882999 0.03746171
## Parp14 0.1563348 0.1269487 0.07891102 0.1142710 0.00000000 0.00000000 0.00000000 0.07074981 0.04568410 0.03135479
## Ddx58  0.2265024 0.2467722 0.21469112 0.2480807 0.00000000 0.00000000 0.00000000 0.07125327 0.04024212 0.02705719
## Parp12 0.1580405 0.1844883 0.14626411 0.1851128 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
## Tap1   0.1949126 0.1937607 0.16257249 0.2563000 0.00000000 0.02717588 0.00000000 0.08010402 0.04710621 0.03769675

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

#Receptors of top-ranked ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


#Sender-focused approach
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 40) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

#Visualizing expression and log-fold change in sender cells
# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(twins, PBMC_l2_predicted.id %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

#Check upregulation of ligands in sender cells
celltype_order <- levels(Idents(twins)) 

library(purrr)
DE_table_top_ligands <- DE_table_top_ligands %>%  purrr::reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ]) 

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc

(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))

#Summary Visualizations
figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot
#####Thelp#####
cd4 <- readRDS("SeuratObjects/CD4.RDS")
DefaultAssay(cd4) <- "SCT"
Idents(cd4) <- "seurat_clusters"

DimPlot(cd4, reduction = "wnn.umapcd4", group.by = "seurat_clusters", label = T)

plot <- DimPlot(cd4, reduction = "wnn.umapcd4", group.by = "seurat_clusters", 
                #split.by = "vaccine", 
                label = F, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  )
plot <- plot + annotate("segment", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]), 
                        xend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]), 
                        yend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]), 
           xend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]), 
           y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]), 
           yend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis
plot
ggsave(filename = "Figure_5A_final_UMAP_CD4.pdf", plot = plot, width = 20, height = 20, dpi = 320, units = "cm", device = "pdf")

t_cell_markers <- c("CD3E", "CD3D", "CD3G", "sct_CD4", "CD8A", "CD8B", "NCR1", "KLRD1", 
                    "NCAM1", "TBX21", "CXCR3", "IFNG", "CCR4", "GATA3", "PTGDR2", 
                    "IL4", "RORC", "IL17A", "CCR6", "FOXP3", "IL7R", "IL2RA", "TGFB1", 
                    "TGFB3", "AHR", "IL22", "SPI1", "BCL6", "CXCR5", "ICOS", 
                    "PDCD1", "IL21", "B3GAT1", "CD69", "GZMB", "PRF1", "SELL", "CCR7", 
                    "CD27", "CD28", "KLRG1", "EOMES", "TCF7", "IKZF1", "MKI67", 
                    "HAVCR2", "TOX", "TOX2", "TIGIT")

DefaultAssay(cd4) <- "SCT"
genes <- as.data.frame(rownames(twins@assays$SCT))

library(viridis)
DotPlot(cd4, features = t_cell_markers, group.by = "seurat_clusters") + 
  scale_color_viridis(option = "H", discrete = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "T Cell Markers", y = "Seurat Clusters")

DotPlot(cd4, features = c("GATA3", "IL4", "IL5", "IL13", "CCR4", "PTGDR2"), group.by = "seurat_clusters") + 
  scale_color_viridis(option = "H", discrete = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "T Cell Markers", y = "Seurat Clusters")

DotPlot(cd4, features = c("BCL6", "CXCR5", "ICOS", "PDCD1", "IL21", "B3GAT1"), group.by = "seurat_clusters") + 
  scale_color_viridis(option = "H", discrete = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "T Cell Markers", y = "Seurat Clusters")

library(Nebulosa)
p <- plot_density(cd4, reduction = "wnn.umapcd4", features = c("GATA3", "IL5", "IL13", "CCR4", "PTGDR2"), joint = T, raster = F)
p

p <- plot_density(cd4, reduction = "wnn.umapcd4", features = c("BCL6", "CXCR5", "ICOS", "PDCD1", "IL21", "B3GAT1"), joint = T, raster = F)
p

cd4 <- AddModuleScore(
  object = cd4,
  features = list(c("BCL6", "CXCR5", "ICOS", "PDCD1", "IL21", "B3GAT1")),
  ctrl = 100,
  name = 'Tfh'
)
FeaturePlot(cd4, reduction = "wnn.umapcd4", features = "Tfh1", order = T, min.cutoff = 0.05)

Idents(cd4) <- "seurat_clusters"
levels(cd4@active.ident)
PrepSCTFindMarkers(cd4)
th2_cluster <- FindAllMarkers(cd4, only.pos = F, logfc.threshold = 0.0,
                              return.thresh = 1,
                              recorrect_umi = F)
print(th2_cluster)
top_genes <- th2_cluster %>% 
  group_by(cluster) %>% 
  top_n(25, avg_log2FC) %>% 
  pull(gene) %>%
  unique()
DoHeatmap(cd4, features=top_genes, label = T, size = 10)

#Cluster 8 is Th2
#Cluster 9 is Th1
#Cluster 7 is Th17
#Cluster 6 = Treg
#Part of Cluster 2 is Tfh
tfh_threshold <- 0.1  # Adjust this cutoff as necessary

# Assign "Tfh" to cells with a high Tfh score
cd4$Thelp <- ifelse(cd4$seurat_clusters == 2 & cd4$Tfh1 > tfh_threshold, "Tfh", NA)

# Assign "Th2" to cells in cluster 8 (assuming 'seurat_clusters' column is available)
cd4$Thelp[cd4$seurat_clusters == 8] <- "Th2"
cd4$Thelp[cd4$seurat_clusters == 9] <- "Th1"
cd4$Thelp[cd4$seurat_clusters == 7] <- "Th17"
cd4$Thelp[cd4$seurat_clusters == 6] <- "Treg"

# Assign "CD4" to all remaining cells that are not Tfh or Th2
cd4$Thelp[is.na(cd4$Thelp)] <- "other CD4"

# Check the distribution of cell types
table(cd4$Thelp)

DimPlot(cd4, 
        reduction = "wnn.umapcd4",  # Specify reduction
        group.by = "Thelp",  # The new metadata column
        order = TRUE)  # Order cells by their label

# Compare the two metadata columns using table
table(cd4$PBMC_l2_predicted.id, cd4$Thelp)
# Check for exact matches between the two columns
overlap <- cd4$PBMC_l2_predicted.id == cd4$Thelp
table(overlap)  # This will give you TRUE/FALSE counts
# Check the number of cells per condition in the 'condition' metadata column
table(cd4$PBMC_l2_predicted.id)
table(cd4$Thelp)

plot <- DimPlot(cd4, reduction = "wnn.umapcd4", group.by = "Thelp", 
                #split.by = "vaccine", 
                label = F, raster = F) + 
  ggtitle("") +
  theme_void() +  # Remove default axes
  theme(strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  )
plot <- plot + annotate("segment", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]), 
                        xend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]) + 3, #This is length of x arrow
                        y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]), 
                        yend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]), 
                        arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]), 
           xend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]), 
           y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]), 
           yend = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]) + 3,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]) + 2, #Horizontal shift
           y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]) - 0.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,1]) - 0.8, #Horizontal shift
           y = min(cd4@reductions$wnn.umapcd4@cell.embeddings[,2]) + 2, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis

plot
ggsave(filename = "Figure_5C_final_UMAP_CD4_Thelp_subsets.pdf", plot = plot, width = 20, height = 20, dpi = 320, units = "cm", device = "pdf")


Idents(cd4) <- "Thelp" 
levels(cd4@active.ident)

library(dplyr)

cluster_frq <- cd4@meta.data %>%
  group_by(subject, vaccine, Thelp) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(relative_freq = count/sum(count)) %>%
  mutate(relative_freq_percent = relative_freq * 100)

# Filter for specific cell types
selected_celltypes <- c("Th1", "Th2", "Th17", "Treg","Tfh")

cluster_frq_subset <- cluster_frq %>%
  filter(Thelp %in% selected_celltypes)

# Set the order of the 'PBMC_l2_predicted.id' factor
cluster_frq_subset$Thelp <- factor(cluster_frq_subset$Thelp, levels = selected_celltypes)

# Plot bar charts with significance markers
vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")

plot <- cluster_frq_subset %>% 
  tidyplot(x = vaccine, y = relative_freq, color = vaccine) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "fisher_test", p.adjust.method = "bonferroni", label = "p.adj.signif", hide_info = T) %>%
  add_test_pvalue(method = "t_test", p.adjust.method = "bonferroni",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Proportion",  fontsize = 8, color = "black") %>%
  adjust_y_axis(limits = c(0,NA), padding = c(0,0.15)) %>%
  adjust_legend_title("Vaccine",  fontsize = 10, color = "black") %>%
  adjust_font(family = "Helvetica", fontsize = 8, color = "black", face = "bold") %>%
  adjust_size(width = 100, height = 100, unit = "cm") %>%
  split_plot(by = Thelp) 

plot
ggsave(filename = "Figure_5C_Th1_Th2_Th17_Treg_Tfh_BoxPlot.pdf", plot = plot, width = 40, height = 20, dpi = 320, units = "cm", device = "pdf")

saveRDS(cd4,"SeuratObjects/CD4.RDS")

#DEG
Idents(cd4) <- "Thelp"
subsets <- levels(cd4@active.ident)

cd4$celltype.vaccine <- paste(cd4$Thelp, cd4$vaccine, sep = "_")
Idents(cd4) <- "celltype.vaccine"
levels(cd4@active.ident)
DefaultAssay(cd4) <- "SCT"

degs.list <- list()
for(i in seq_along(subsets)) {
  degs.list[[i]] <- FindMarkers(cd4,
                                ident.1 = paste(subsets[i], "CERVARIX", sep = "_" ),
                                ident.2 = paste(subsets[i], "GARDASIL", sep = "_" ), recorrect_umi = F)
}
names(degs.list) <- subsets

gene_set_r <- msigdbr::msigdbr(category = 'C2', subcategory = 'CP:REACTOME') %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  unique()

# Ranked list of genes
ranks.list <- list()
for (i in seq_along(degs.list)) {
  ranks.list[[i]] <- degs.list[[i]] %>%
    mutate(p_val = ifelse(p_val == 0, .Machine$double.xmin*row_number(), p_val)) %>% 
    mutate(rank = -log10(p_val)*sign(avg_log2FC)) %>% 
    dplyr::select(rank) %>% 
    rownames_to_column() %>% 
    arrange(dplyr::desc(rank)) %>% 
    deframe() 
}
names(ranks.list) <- names(degs.list)

# Gene set enrichment analysis
gsea.list <- compareCluster(ranks.list, fun = 'GSEA', 
                            TERM2GENE = gene_set_r, eps = 0)

#Save
GSEA_df <- gsea.list@compareClusterResult
GSEA_df$Cluster <- as.factor(GSEA_df$Cluster)
levels(GSEA_df$Cluster)
# Split the data frame by the 'Cluster' column
clusters <- split(GSEA_df, GSEA_df$Cluster)
# Write the list of data frames to an Excel file
write.xlsx(clusters, file = "GSEA_Thelp_all_Reactome.xlsx", rowNames = T)
levels(cd4@active.ident)
# Convert the x-axis variable to a factor if it is not already
gsea.list@compareClusterResult$Cluster <- factor(gsea.list@compareClusterResult$Cluster, levels = subsets)

th_gsea_list <- gsea.list
th_gsea_list@compareClusterResult <- subset(gsea.list@compareClusterResult, Cluster %in% selected_celltypes)
th_gsea_list@compareClusterResult$Cluster <- factor(th_gsea_list@compareClusterResult$Cluster, levels = selected_celltypes, ordered = TRUE)

# Dotplot with enriched pahtways per cell type
dotplot(th_gsea_list, color = 'NES', label_format = 35, showCategory = Inf) +
  scale_fill_viridis_c(option = 'H') +
  RotatedAxis() +
  theme(axis.title.x = element_blank(),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 10),  # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis tick labels font and size
        #panel.grid = element_blank()
  )

notch <- c("REACTOME_NOTCH2_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
           "REACTOME_NOTCH3_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
           "REACTOME_SIGNALING_BY_NOTCH1_HD_DOMAIN_MUTANTS_IN_CANCER",
           "REACTOME_SIGNALING_BY_NOTCH2",
           "REACTOME_ACTIVATED_NOTCH1_TRANSMITS_SIGNAL_TO_THE_NUCLEUS",
           "REACTOME_SIGNALING_BY_NOTCH",
           "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
           "REACTOME_KEAP1_NFE2L2_PATHWAY")

GSEA_filtered <- GSEA_df[GSEA_df$Description %in% notch, ]

gsea.list@compareClusterResult <- GSEA_filtered

filtered_gsea <- gsea.list[gsea.list@compareClusterResult$Cluster %in% selected_celltypes, ]
filtered_gsea$Cluster <- factor(filtered_gsea$Cluster, levels = selected_celltypes, order = T)
gsea.list@compareClusterResult <- filtered_gsea

dotplot(gsea.list, color = 'NES', label_format = 35, showCategory = Inf) +
  scale_fill_viridis_c(option = 'H') +
  RotatedAxis() +
  theme(axis.title.x = element_blank(),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 10),  # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis tick labels font and size
        #panel.grid = element_blank()
  )

ggsave(filename = "DotPlot_Pathways_Th_all_Reactome.pdf", width = 20, height = 120, dpi = 320, units = "cm", device = "pdf")

#####IgG4 Production#####
twins.b <- readRDS("SeuratObjects/B.RDS")
DefaultAssay(twins.b) <- "SCT"
Idents(twins.b) <- "PBMC_l2_predicted.id"

DimPlot(twins.b, reduction = "wnn.umapb", group.by = "PBMC_l2_predicted.id")
DimPlot(twins.b, reduction = "wnn.umapb", group.by = "seurat_clusters", label = T)
DimPlot(twins.b, reduction = "wnn.umapb", group.by = "vaccine")

igg <- c("IGHG1","IGHG2","IGHG3","IGHG4")

igg4 <- c("IGHG4", "AICDA", "UNG", "CD40", "CD27", #Genes Involved in Class Switch Recombination (CSR) to IgG4
          "BCL6", "IRF4", "PRDM1", "XBP1", #Genes Involved in B Cell Activation and Differentiation
          "IL4R", "STAT6", "IL21R", "TNFRSF13B", "TNFRSF13C", #Th2 Cytokine and Signaling Pathways
          "SLAMF7", "HSPA5", "SDC1", "CXCR4", #Plasmablast Survival and Antibody Secretion
          "TNFAIP3", "SOCS1", "SOCS3", "FOXP1",#Other Genes Supporting IgG4 Production
          "IGHG1","IGHG2","IGHG3") 

b <- c("CD19", "MS4A1", "SDC1", "TNFRSF17")

DotPlot(twins.b, features = igg, group.by = "seurat_clusters") + 
  scale_color_viridis(option = "H", discrete = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "B Cell Markers", y = "Seurat Clusters")

DotPlot(twins.b, features = igg4, group.by = "seurat_clusters") + 
  scale_color_viridis(option = "H", discrete = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "B Cell Markers", y = "Seurat Clusters")

DotPlot(twins.b, features = b, group.by = "seurat_clusters") + 
  scale_color_viridis(option = "H", discrete = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "B Cell Markers", y = "Seurat Clusters")

library(pheatmap)

DefaultAssay(twins.b) <- "SCT"

# Extract the expression matrix for the top DEGs
expr_data <- GetAssayData(twins.b, slot = "data")[igg4, , drop = FALSE]

# Add subject and cell type information to the metadata
metadata <- twins.b@meta.data
metadata$subject <- gsub("_", "-", metadata$subject)
metadata$subject_celltype_vaccine <- paste(metadata$subject, metadata$PBMC_l2_predicted.id, metadata$vaccine,sep = "_")

# Calculate average expression per subject per cell type
aggregated_expr_data <- expr_data %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(subject_celltype_vaccine = metadata$subject_celltype_vaccine) %>%
  group_by(subject_celltype_vaccine) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("subject_celltype_vaccine") %>%
  t()

# Scale the data for heatmap visualization
scaled_expr_matrix <- t(scale(t(aggregated_expr_data)))

annotation_col <- data.frame(
  Subject = gsub("_.*", "", colnames(scaled_expr_matrix)),  # Extract subject (everything before the first '_')
  CellType = gsub("^[^_]+_([^_]+)_.*", "\\1", colnames(scaled_expr_matrix)),  # Extract cell type (everything after the first '_' but before the last)
  Vaccine = gsub(".*_", "", colnames(scaled_expr_matrix))  # Extract vaccine (everything after the last '_')
)

rownames(annotation_col) <- colnames(scaled_expr_matrix)

#Remove Subject
annotation_col <- annotation_col[, !colnames(annotation_col) %in% "Subject"]

# Order columns by cell type, then by Vaccine
annotation_col <- annotation_col %>%
  arrange(CellType, Vaccine)

# Reorder the columns in `scaled_expr_matrix` to match this sorted order
scaled_expr_matrix <- scaled_expr_matrix[, rownames(annotation_col)]

# Create the heatmap
hm.parameters <- list(scaled_expr_matrix,
                      scale = "row",
                      cellwidth = 5, cellheight= 5,
                      color = colorRampPalette(c("blue4","blue1","white","red1","brown1"))(100),
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = F,
                      fontsize_row = 5,
                      main = "",
                      clustering_method = "ward.D2",
                      cluster_rows = T, cluster_cols = F,
                      clustering_distance_rows = "euclidean",
                      clustering_distance_cols = "euclidean",
                      annotation_col = annotation_col)

kmean.hm <- do.call("pheatmap", hm.parameters)
print(kmean.hm)

#bmemory <- subset(twins.b, subset = PBMC_l2_predicted.id == "B memory")
#VlnPlot(bmemory, features = "IGHG4")

# Reshape the expression matrix into a long format for ggplot
aggregated_expr_data_long <- as.data.frame(aggregated_expr_data) %>%
  tibble::rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "subject_celltype_vaccine", values_to = "expression")

# Split the `subject_celltype_vaccine` column into separate columns for `subject`, `celltype`, and `vaccine`
aggregated_expr_data_long <- aggregated_expr_data_long %>%
  separate(subject_celltype_vaccine, into = c("subject", "celltype", "vaccine"), sep = "_")

aggregated_expr_data_long$subject <- as.factor(aggregated_expr_data_long$subject)
aggregated_expr_data_long$celltype <- as.factor(aggregated_expr_data_long$celltype)
aggregated_expr_data_long$vaccine <- as.factor(aggregated_expr_data_long$vaccine)

# Filter the data to only include genes of interest from the `igg4` list
aggregated_expr_data_long <- aggregated_expr_data_long %>%
  filter(gene %in% igg4)

# Filter for the specific cell type (e.g., "Treg")
aggregated_expr_data_long <- aggregated_expr_data_long %>%
  filter(celltype == "B memory")

# Create the violin plot using ggplot2
library(ggplot2)

IGG <- c("IGHG1","IGHG2","IGHG3","IGHG4")
CSR <- c("CD40", "CD27", "AICDA", "UNG")
TFs <- c("BCL6", "IRF4", "PRDM1", "XBP1")
Th2signal <- c("IL4R", "STAT6", "IL21R", "TNFRSF13B", "TNFRSF13C","TNFAIP3")

# Filter the data to include only the selected genes
filtered_data <- aggregated_expr_data_long %>%
  filter(gene %in% IGG)

vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")
plot <- ggplot(filtered_data, aes(x = celltype, y = expression, fill = vaccine)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") +  # Facet by gene for individual plots
  scale_fill_manual(values = vaccine_colors, labels = c("Cervarix", "Gardasil")) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = " ",
       x = "", y = "Expression", fill = "Vaccine")
plot

filtered_data %>% 
  tidyplot(x = vaccine, y = expression, color = vaccine) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  add_test_asterisks(method = "wilcoxon", p.adjust.method = "bonferroni", label = "p.adj.signif", hide_info = T) %>%
  add_test_pvalue(method = "wilcoxon", p.adjust.method = "bonferroni",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Average Expression",  fontsize = 8, color = "black") %>%
  adjust_y_axis(limits = c(0,NA), padding = c(0.1,0.15)) %>%
  adjust_legend_title("Vaccine",  fontsize = 10, color = "black") %>%
  adjust_font(family = "Helvetica", fontsize = 8, color = "black", face = "bold") %>%
  adjust_size(width = 100, height = 100, unit = "cm") %>%
  split_plot(by = gene) 

ggsave(filename = "Figure_5F_Increased_IgG4_in_BMem.pdf", plot = plot, width = 40, height = 20, dpi = 320, units = "cm", device = "pdf")

#DEG
Idents(twins.b) <- "PBMC_l2_predicted.id"
subsets <- levels(twins.b@active.ident)

twins.b$celltype.vaccine <- paste(twins.b$PBMC_l2_predicted.id, twins.b$vaccine, sep = "_")
Idents(twins.b) <- "celltype.vaccine"
levels(twins.b@active.ident)
DefaultAssay(twins.b) <- "SCT"

degs.list <- list()
for(i in seq_along(subsets)) {
  degs.list[[i]] <- FindMarkers(twins.b,
                                ident.1 = paste(subsets[i], "CERVARIX", sep = "_" ),
                                ident.2 = paste(subsets[i], "GARDASIL", sep = "_" ), recorrect_umi = F)
}
names(degs.list) <- subsets

degsBMEM <- degs.list[["B memory"]]
degsBMEM$gene <- rownames(degsBMEM)
BMEM_genes <- degsBMEM[degsBMEM$gene %in% igg4, ]

write.csv(BMEM_genes, file = "Supp_Table14_DEGs_IgG4_Production_BMemory.csv", row.names = FALSE)

#####CellChat Th to B Mem#####
# Add Tfh and Th annotations from the 'cd4' object to the 'twins' object
# Ensure 'twins' contains the metadata column 'PBMC_l2_predicted.id'

# Extract cell IDs for Tfh and Th2 from the 'cd4' object
selected_celltypes <- c("Th1", "Th2", "Th17", "Treg","Tfh")
th_cells <- colnames(cd4)[cd4$Thelp %in% selected_celltypes]

# Add annotations to the 'twins' object
twins$Thelp <- twins$PBMC_l2_predicted.id  # Start with the original annotation

# Update cells that are in 'cd4' and labeled 'Tfh' or 'Th2'
twins$Thelp[colnames(twins) %in% th_cells] <- cd4$Thelp[cd4$Thelp %in% selected_celltypes]

table(twins$Thelp)
table(cd4$Thelp)

saveRDS(twins, "SeuratObjects/HPV-X-Neutra.RDS")
twins <- readRDS("SeuratObjects/HPV-X-Neutra.RDS")

Idents(twins) <- "Thelp"
levels(twins@active.ident)

cells <- c("Th1", "Th2", "Th17",
           "B naive", "B intermediate", "B memory","Plasmablast")

BT <- subset(twins, subset = Thelp %in% cells)
DimPlot(BT, reduction = "wnn.umap", group.by = "Thelp")
BT <- FindMultiModalNeighbors(
  object = BT,
  reduction.list = list("harmony", "integrated_lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
BT <- RunUMAP(BT, nn.name = "weighted.nn", reduction.name = 'wnn.umapcd4b', reduction.key = 'rnaUMAPcd4b_')
#Clustering
ElbowPlot(BT, ndims=50)
BT <- FindNeighbors(BT, dims = 1:30)
# You should make sure your assay is set correctly (the assay that you originally run PCA).
DefaultAssay(BT) <- "SCT"
# Now we cluster the resultant graph.
#BT <- FindClusters(BT, resolution = 0.5)
DimPlot(BT, reduction = "wnn.umapcd4b", group.by = "Thelp")
saveRDS(BT, "SeuratObjects/CD4_B.RDS")
BT <- readRDS("SeuratObjects/CD4_B.RDS")

Idents(BT) <- "Thelp"
levels(BT@active.ident)

#Select CERVARIX or GARDASIL
query_CER <- subset(x=BT, subset = vaccine == "CERVARIX")
dim(query_CER)
query_GAR <- subset(x=BT, subset = vaccine == "GARDASIL")
dim(query_GAR)

#Cervarix
data.input <- query_CER[["SCT"]]$data # normalized data matrix
labels <- Idents(query_CER)
meta <- data.frame(labels = labels, row.names = names(labels))

cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellChat@DB <- CellChatDB.use

#subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 3 * 1024^3)  # 3 GiB
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 699

cellChat <- computeCommunProb(cellChat, type = "triMean")

cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)

groupSize <- as.numeric(table(cellChat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

debug(netAnalysis_signalingRole_heatmap)
undebug(netAnalysis_signalingRole_heatmap)
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing", height = 15)
ht1
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming", height = 15)
ht2

saveRDS(cellChat, file = "SeuratObjects/cellChat_query_Cervarix.rds")
cellchat.CER <- cellChat

#Gardasil
data.input <- query_GAR[["SCT"]]$data # normalized data matrix
labels <- Idents(query_GAR)
meta <- data.frame(labels = labels, row.names = names(labels))

cellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
cellChat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 3 * 1024^3)  # 3 GiB
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 716

cellChat <- computeCommunProb(cellChat, type = "triMean")

cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)

groupSize <- as.numeric(table(cellChat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing", height = 15)
ht1
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming", height = 15)
ht2

saveRDS(cellChat, file = "SeuratObjects/cellChat_query_Gardasil.rds")
cellchat.GAR <- cellChat

#Compare Cer and Gar
cellchat.CER <- readRDS("cellChat_query_Cervarix.rds")
cellchat.GAR <- readRDS("cellChat_query_Gardasil.rds")

object.list <- list(GAR = cellchat.GAR, CER = cellchat.CER)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

save(object.list, file = "SeuratObjects/cellchat_object.list.RData")
save(cellchat, file = "SeuratObjects/cellchat_merged.RData")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

cd4 <- c("Th1", "Th2", "Th17", "Treg","Tfh")
b <- c("B naive", "B intermediate", "B memory","Plasmablast")
cells <- c("Th1", "Th2", "Th17", "Treg","Tfh", "B naive", "B intermediate", "B memory","Plasmablast")

debug(netVisual_heatmap)       
undebug(netVisual_heatmap)

celltype_colors <- c("Th1" = "blue",  # Red
                     "Th2" = "green",  # Blue
                     "Th17" = "purple",   # Green,
                     "Treg" = "turquoise",
                     "Tfh" = "lightgrey",
                     "B naive" = "#984EA3",  # Purple
                     "B intermediate" = "#FF7F00",  # Orange
                     "B memory" = "#FFFF33",  # Yellow
                     "Plasmablast" = "#A65628")  # Brown

cellChat@idents <- factor(cellChat@idents, levels = cells, ordered = TRUE)

gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",
                           color.use = c("1" = "#007A73", "2" = "#E4AD2E"), color.alpha = 0.5,
                           digits = 2, size.text = 14, width = 0.5,
                           xlabel = "Vaccine", remove.xtick = F) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x.top = element_blank(),
        axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"), 
        axis.text.y = element_text(family = "Helvetica", size = 10, color = "black"), # Y-axis title font and size
        axis.text.x = element_text(family = "Helvetica", size = 10, color = "black"),  # X-axis tick labels font and size
        panel.grid = element_blank(),
        legend.text = element_text(family = "Helvetica", size = 10, color = "black"),
        legend.title = element_text(family = "Helvetica", size = 12, face = "bold")
  )

gg2

gg1 <- netVisual_heatmap(cellchat, row.show = cd4, col.show = b, cluster.rows = T, cluster.cols = T)
gg1
gg2 <- netVisual_heatmap(cellchat, measure = "weight",row.show = cd4, col.show = b, 
                         cluster.rows = T, cluster.cols = T, color.use = celltype_colors,
                         width = 2, height = 2, font.size = 12, font.size.title = 16
) 
gg2


####NicheNet Cytokines to B Mem####
library(nichenetr)
DefaultAssay(twins) <- "SCT"
twins <- SetIdent(twins, value = "Thelp")
levels(twins@active.ident)

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

genes_list <- c(
  "IL1B", "IL1R1", "IL1RN", "IL1R2", "IL1RL1", "IL1RL2", "IL2", "IL2RA", "IL2RB", "IL2RG", 
  "IL3", "IL3RA", "IL4", "IL4R", "IL5", "IL5RA", "IL6", "IL6R", "IL6ST", "IL7", "IL7R", 
  "IL9", "IL9R", "IL10", "IL10RA", "IL10RB", "IL11", "IL11RA", "IL12A", "IL12B", "IL12RB1", 
  "IL12RB2", "IL13", "IL13RA1", "IL13RA2", "IL15", "IL15RA", "IL17A", "IL17F", "IL17RA", 
  "IL17RB", "IL17RC", "IL18", "IL18BP", "IL18R1", "IL18RAP", "IL22", "IL22RA1", "IL22RA2", 
  "IL23A", "IL23R", "IL27", "IL27RA", "IL31", "IL31RA", "IL33", "IL34", "TSLP", 
  "IFNA", "IFNB1", "IFNG", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1", "IFITM1", 
  "IFITM2", "IFITM3", "IFIT1", "IFIT2", "IFIT3", "ISG15", "IRF3", "IRF9", 
  "TNF", "TNFRSF1A", "TNFRSF1B", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF4", "TNFSF8", 
  "TNFSF9", "TNFSF11", "TNFRSF11A", "TNFRSF11B", "TNFRSF12A", "TNFSF12", "TNFRSF18", 
  "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", 
  "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCL11", "CCL19", "CCL20", 
  "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CXCL9", "CXCL10", "CXCL11", 
  "FASLG", "CD4", "CD27", "CD36", "CD40LG", "CD70", "CSF1", "CSF1R", "CSF2", 
  "CSF2RA", "CSF2RB", "CSF3", "CSF3R", "CNTF", "CNTFR", "LIF", "LIFR", "OSM", 
  "OSMR", "HGF", "PTGS2", "MAIT"
)

filtered_values <- lr_network[lr_network$from %in% genes_list, ]

lr_network <- filtered_values %>% distinct(from, to)
head(lr_network)

#Sender-agnostic
receiver = "B memory"
expressed_genes_receiver <- get_expressed_genes(receiver, twins, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#Sender-focused
sender_celltypes <- c("Th1", "Th2", "Th17", "Tfh")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, twins, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 6729
length(potential_ligands)
## [1] 31
length(potential_ligands_focused)
## [1] 3

#Gene set of interest
condition_oi <-  "CERVARIX"
condition_reference <- "GARDASIL"

twins_receiver <- subset(twins, idents = receiver)

DE_table_receiver <-  FindMarkers(object = twins_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "vaccine",
                                  min.pct = 0.05, recorrect_umi = FALSE) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define the background
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 5144
length(geneset_oi)
## [1] 465

#Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(dplyr::desc(aupr_corrected)))
best_upstream_ligands <- ligand_activities %>% top_n(40, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange", bins = 50)  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(40, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% dplyr::arrange(aupr_corrected) %>% base::as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange", size = 5) + 
    theme(axis.text.x.top = element_blank()))  

#Infer target genes and receptors of top-ranked ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix, n = 100) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)
## [1] 90
head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.1) 

nrow(active_ligand_target_links)
## [1] 29
head(active_ligand_target_links)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes", #[c(1,2,4,5,7,8,12,20,21,24),c(2,11,14,15)]
                    color = "purple", legend_title = "Regulatory potential", size = 4) +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple") +
  theme_bw() + 
  theme(#legend.position = "none",
    #axis.text.x.top = element_blank(),
    axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
    axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"), 
    axis.text.y = element_text(family = "Helvetica", size = 10, color = "black"), # Y-axis title font and size
    axis.text.x = element_text(family = "Helvetica", size = 10, color = "black", , angle = 45, hjust = -0.01),  # X-axis tick labels font and size
    panel.grid = element_blank(),
    legend.text = element_text(family = "Helvetica", size = 10, color = "black"),
    legend.title = element_text(family = "Helvetica", size = 12, face = "bold")
  )

df <- as.data.frame(vis_ligand_target[c(1,2,4,5,7,8,12,20,21,24),c(2,11,14,15)])

#####Differentially Accessable Regions#####
DefaultAssay(twins) <- "ATAC"
Idents(twins) <- "celltype.vaccine"
dar_pdc <- FindMarkers(
  twins,
  ident.1 = "cDC1_CERVARIX",     # Cells treated with CERVARIX in the pDC population
  ident.2 = "cDC1_GARDASIL",     # Cells treated with GARDASIL in the pDC population
  test.use = "LR",              # Logistic regression for DAR analysis
  latent.vars = "atac_peak_region_fragments", # Correct for sequencing depth
  min.pct = 0.1                 # Only consider peaks accessible in ≥10% of cells
)

#Annotate and Write to an xlsx file with default parameters
annotated_dar <- ClosestFeature(twins, regions = rownames(dar_pdc), sep = c(":", "-"))
dar_pdc <- cbind(dar_pdc, gene = annotated_dar$gene_name, distance = annotated_dar$distance)
write.xlsx(dar_pdc, file = "cDC1_DifferentiallyAccessbileRegions.xlsx", rowNames = T)

dar_pdc_sig <- dar_pdc %>%
  rownames_to_column(var = "coord") %>%
  dplyr::filter(avg_log2FC > 0) %>%  # Filter for regions with positive log2FC
  dplyr::distinct(gene, .keep_all = TRUE) %>% 
  dplyr::mutate(celltype = "pDC")  # Annotate as pDC peaks

dar_pdc_top200 <- dar_pdc_sig %>%
  top_n(200, wt = avg_log2FC)

write.xlsx(dar_pdc_top200, file = "cDC1_Top200_DifferentiallyAccessbileRegions.xlsx", rowNames = T)

twins_pdc <- subset(twins, idents = c("cDC1_CERVARIX", "cDC1_GARDASIL"))

# Extract the accessibility values for these top 200 regions
top200_regions <- dar_pdc_top200$coord
accessibility_data <- FetchData(twins_pdc, vars = c(top200_regions, "subject")) #Only 6 subjects: S10 and S14 don't have cDC1

region_gene_mapping <- dar_pdc_top200 %>%
  dplyr::select(coord, gene) 

# Aggregate by subject: calculate mean accessibility per region per subject
aggregated_data <- accessibility_data %>%
  group_by(subject) %>%
  summarise(across(all_of(top200_regions), mean, na.rm = TRUE))

aggregated_matrix <- as.matrix(aggregated_data[, -1])  # Exclude the 'subject' column
rownames(aggregated_matrix) <- aggregated_data$subject
colnames(aggregated_matrix) <- region_gene_mapping$gene
aggregated_matrix <- t(aggregated_matrix)

# Classify regions based on distance
classify_region_distance <- function(distance) {
  if (distance >= -2000 & distance <= 500) {
    return("promoter")
  } else if (distance < -10000 | distance > 10000) {
    return("trans")
  } else {
    return("distal")
  }
}

# Apply this function to classify each region based on its distance
dar_pdc_top200$region_type <- sapply(dar_pdc_top200$distance, classify_region_distance)

region_annot <- data.frame(
  Region = dar_pdc_top200$region_type
)

# Ensure the annotation rows match the DAR coordinates in the heatmap
row.names(region_annot) <- dar_pdc_top200$gene

# Create a color scheme for the region types
region_colors <- c("promoter" = "#D6A6CD", "distal" = "#93BFC0", "trans" = "#A3DAB1")

annot <- data.frame(
  Subject = c("Subject_005", "Subject_006", "Subject_007", "Subject_008", "Subject_009", "Subject_013"), #"Subject 010",, "Subject 014"
  Vaccine = c("CERVARIX", "GARDASIL", "CERVARIX", "GARDASIL", "GARDASIL",  "CERVARIX") #"CERVARIX",, "GARDASIL"
)
rownames(annot) <- annot$Subject
annot$Subject <- NULL
annot$Vaccine <- as.factor(annot$Vaccine)

annotation_colors <- list(
  Vaccine = c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73"),
  Region = region_colors
)
# Create a annotation_row# Create a heatmap using pheatmap
plot <- pheatmap(aggregated_matrix, 
                 cluster_rows = F, 
                 cluster_cols = T, 
                 show_rownames = FALSE, 
                 show_colnames = FALSE, 
                 scale = "row", 
                 main = "",
                 annotation_col = annot,
                 annotation_row = region_annot,
                 annotation_colors = annotation_colors, 
                 cellwidth = 15,
                 cellheight = 1,
                 name = "z-score"
)

plot

library(ComplexHeatmap)
library(circlize)

# Function to scale rows (z-score normalization)
scale_rows <- function(mat) {
  t(apply(mat, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
}

# Scale the matrix
scaled_matrix <- scale_rows(aggregated_matrix)

# Define colors
region_colors <- c("promoter" = "#D6A6CD", "distal" = "#93BFC0", "trans" = "#A3DAB1")
vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")

# Create row annotation for region type
row_ha <- rowAnnotation(
  Region = region_annot$Region,
  col = list(Region = region_colors),
  annotation_legend_param = list(title = "Region Type")
)

# Create column annotation for vaccine
col_ha <- HeatmapAnnotation(
  Vaccine = annot$Vaccine,
  col = list(Vaccine = vaccine_colors),
  annotation_legend_param = list(title = "Vaccine")
)

# Plot heatmap
plot <- Heatmap(
  scaled_matrix,
  name = "z-score",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  #row_names_gp = gpar(fontsize = 8, family = "Helvetica", color = "black"),
  #column_names_gp = gpar(fontsize = 12, family = "Helvetica", color = "black", rotation = 45),
  top_annotation = col_ha,
  left_annotation = row_ha,
  na_col = "white", 
  column_title = "",
  #row_names_side = "left", 
  #rect_gp = gpar(col = "white", lwd = 2),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  column_names_rot = 45,
  width = unit(2.5, "cm"), 
  height = unit(20, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Helvetica", fontface = "bold"),  # Customize legend title font
    labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica")  # Customize legend label font
  ))   
plot

ggsave(filename = "Heatmap_Top200_Significant_DAR_cDC1.pdf", plot = plot, width = 20, height = 20, dpi = 320, units = "cm", device = "pdf")

notch_genes <- c(
  "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",  # NOTCH family receptors
  "DLL1", "DLL3", "DLL4",  # Delta-like ligands
  "JAG1", "JAG2",  # Jagged ligands
  "MAML1",  # Mastermind-like 1
  "RBPJ",  # Transcriptional regulator of NOTCH signaling
  "HES1", "HES5",  # HES transcription factors
  "HEY1", "HEY2", "HEYL",  # Hey family transcription factors
  "CSL"  # CBF1, also known as CSL (receptor binding)
)

ets_genes <- c(
  "ETS1", "ETS2", "ELK1", "ELK3", "ELK4", "ERG", "FLI1", "PU.1", "SPI1", "EGR1", "EGR2", "EGR3", "EGR4"  # ETS transcription factors
)

selected_tfs <- c("NRF1", "KLF15", "KLF2", "EGR1", "ZBTB14", "KLF14", "HINFP", "SP9",
                  "EGR3", "KLF3", "Klf12", "SP3", "SP4", "SP1", "SP2", "EGR2", "KLF6",
                  "TFAP2B", "TFAP2E", "KLF4", "TFDP1", "E2F6", "KLF10", "KLF11", "ZNF148",
                  "CTCFL", "TCFL5", "ZBTB33", "TFAP2C", "SP8", "Zfx", "KLF5", "Wt1",
                  "ATF7", "EGR4")

# Combine both lists
combined_genes <- c(notch_genes, ets_genes, selected_tfs)

combined_genes <- c("RPS6", "RPS23", "RPS17", "RPSA", "EIF3H", "EIF3A", "RPS5", "RPS27", "RPS8", "RPS27A", "EIF3L", "EIF4A1", "RPS4X", "EIF4H", "RPS19", "RPS11", "RPS21", "PABPC1", "RPS15A", "RPS24", "FAU", "RPS2", 
                    "RBPJ", "UBC", "MAML2", "CREB1", "RPS27A", "NOTCH2", "APH1A", "EP300", 
                    "EEF1A1", "EEF2", "CAMKMT", "HSPA8", combined_genes) %>% unique()

notch <- c("TCFL5" ,"APH1A", "SP3", "KLF15", "RPS4X", "RPS15A", "EEF2",
           "FLI1", "EGR4", "RPS21", "DLL4", "EIF4A1", "KLF14", "KLF10",
           "HEY1", "RPS11", "HES5", "EIF3A", "NOTCH1", "KLF11", "ELK1",
           "RPS17", "ELK4", "JAG1")

dar_pdc_notch <- dar_pdc %>%
  rownames_to_column(var = "coord") %>%
  dplyr::filter(gene %in% combined_genes) %>%  # Filter for regions with positive log2FC
  dplyr::distinct(gene, .keep_all = TRUE) %>%  
  dplyr::mutate(celltype = "pDC") 

write.xlsx(dar_pdc_notch, file = "pDC1_Notch_DifferentiallyAccessbileRegions.xlsx", rowNames = T)

notch_regions <- dar_pdc_notch$coord
accessibility_data <- FetchData(twins_pdc, vars = c(notch_regions, "subject"))

region_gene_mapping <- dar_pdc_notch %>%
  dplyr::select(coord, gene) 

# Aggregate by subject: calculate mean accessibility per region per subject
aggregated_data <- accessibility_data %>%
  group_by(subject) %>%
  summarise(across(all_of(notch_regions), mean, na.rm = TRUE))

aggregated_matrix <- as.matrix(aggregated_data[, -1])  # Exclude the 'subject' column
rownames(aggregated_matrix) <- aggregated_data$subject
colnames(aggregated_matrix) <- region_gene_mapping$gene
aggregated_matrix <- t(aggregated_matrix)

# Classify regions based on distance
classify_region_distance <- function(distance) {
  if (distance >= -2000 & distance <= 500) {
    return("promoter")
  } else if (distance < -10000 | distance > 10000) {
    return("trans")
  } else {
    return("distal")
  }
}

# Apply this function to classify each region based on its distance
dar_pdc_notch$region_type <- sapply(dar_pdc_notch$distance, classify_region_distance)

region_annot <- data.frame(
  Region = dar_pdc_notch$region_type
)

# Ensure the annotation rows match the DAR coordinates in the heatmap
row.names(region_annot) <- dar_pdc_notch$gene

# Create a color scheme for the region types 
region_colors <- c("promoter" = "#D6A6CD", "distal" = "#93BFC0", "trans" = "#A3DAB1")

annot <- data.frame(
  Subject = c("Subject_005", "Subject_007", "Subject_013", "Subject_006",  "Subject_008", "Subject_009"), #,, 
  Vaccine = c("CERVARIX","CERVARIX" , "CERVARIX", "GARDASIL", "GARDASIL", "GARDASIL") #,
)
rownames(annot) <- annot$Subject
annot$Subject <- NULL
annot$Vaccine <- as.factor(annot$Vaccine)

annotation_colors <- list(
  Vaccine = c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73"),
  Region = region_colors
)

# Assuming ordered_subjects is a vector containing the desired column order
ordered_subjects <- c("Subject_005", "Subject_007", "Subject_013", "Subject_006",  "Subject_008", "Subject_009") #"Subject_010",", "Subject_014" )
aggregated_matrix <- aggregated_matrix[, ordered_subjects]

# Create a annotation_row# Create a heatmap using pheatmap
plot <- pheatmap(aggregated_matrix, 
                 cluster_rows = T, 
                 cluster_cols = F, 
                 show_rownames = T, 
                 show_colnames = FALSE, 
                 scale = "row", 
                 main = "",
                 annotation_col = annot,
                 annotation_row = region_annot,
                 annotation_colors = annotation_colors, 
                 cellwidth = 10,
                 cellheight = 10,
                 #cutree_rows = 2,
                 #cutree_cols = 1,
                 name = "z-score"
)

plot

# Scale the matrix
scaled_matrix <- scale_rows(aggregated_matrix)

# Define colors
region_colors <- c("promoter" = "#D6A6CD", "distal" = "#93BFC0", "trans" = "#A3DAB1")
vaccine_colors <- c("CERVARIX" = "#E4AD2E", "GARDASIL" = "#007A73")

# Create row annotation for region type
row_ha <- rowAnnotation(
  Region = region_annot$Region,
  col = list(Region = region_colors),
  annotation_legend_param = list(title = "Region Type")
)

# Create column annotation for vaccine
col_ha <- HeatmapAnnotation(
  Vaccine = annot$Vaccine,
  col = list(Vaccine = vaccine_colors),
  annotation_legend_param = list(title = "Vaccine")
)

# Plot heatmap
plot <- Heatmap(
  scaled_matrix,
  name = "z-score",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = T,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8, family = "Helvetica", color = "black"),
  #column_names_gp = gpar(fontsize = 12, family = "Helvetica", color = "black", rotation = 45),
  top_annotation = col_ha,
  left_annotation = row_ha,
  na_col = "white", 
  column_title = "",
  row_names_side = "right", 
  #rect_gp = gpar(col = "white", lwd = 2),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  column_names_rot = 45,
  width = unit(2.5, "cm"), 
  height = unit(10, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12, fontfamily = "Helvetica", fontface = "bold"),  # Customize legend title font
    labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica")  # Customize legend label font
  ))   
plot

ggsave(filename = "Heatmap_NOTCH_Significant_DAR_pDC.pdf", plot = plot, width = 15, height = 15, dpi = 320, units = "cm", device = "pdf")

DefaultAssay(twins.DC) <- "ATAC"
peaks_of_interest <- c("chr8-79767610-79768539", "chr15-40941526-40942423", #HEY1, DLL4
                       "chr1-2529080-2529966", "chr9-136539451-136540222", #HES5, NOTCH1
                       "chrX-47650288-47651211", "chr20-10673834-10674648") #ELK1, JAG1

p <- plot_density(twins.DC, reduction = "wnn.umapDC", 
             features = peaks_of_interest, 
             joint = FALSE, raster = FALSE)
p1 <- p[[1]] + theme_bw() + ggtitle("HEY1") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),   # Remove axis text (tick labels)
        axis.ticks = element_blank(),  # Remove axis ticks
        panel.grid = element_blank(),   # Remove gridlines
        strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  ) +
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
             xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) + 0.8, #This is length of x arrow
             y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
             yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) -1, 
             arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) -2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 2 ,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 1, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2.5, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 0.5, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis

p2 <- p[[2]] + theme_bw() + ggtitle("DLL4") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),   # Remove axis text (tick labels)
        axis.ticks = element_blank(),  # Remove axis ticks
        panel.grid = element_blank(),   # Remove gridlines
        strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  ) +
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) + 0.8, #This is length of x arrow
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) -1, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) -2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 2 ,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 1, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2.5, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 0.5, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis

p3 <- p[[3]] + theme_bw() + ggtitle("HES5") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),   # Remove axis text (tick labels)
        axis.ticks = element_blank(),  # Remove axis ticks
        panel.grid = element_blank(),   # Remove gridlines
        strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  ) +
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) + 0.8, #This is length of x arrow
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) -1, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) -2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 2 ,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 1, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2.5, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 0.5, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis

p4 <- p[[4]] + theme_bw() + ggtitle("NOTCH1") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),   # Remove axis text (tick labels)
        axis.ticks = element_blank(),  # Remove axis ticks
        panel.grid = element_blank(),   # Remove gridlines
        strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  ) +
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) + 0.8, #This is length of x arrow
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) -1, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) -2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 2 ,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 1, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2.5, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 0.5, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis

p5 <- p[[5]] + theme_bw() + ggtitle("ELK1") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),   # Remove axis text (tick labels)
        axis.ticks = element_blank(),  # Remove axis ticks
        panel.grid = element_blank(),   # Remove gridlines
        strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  ) +
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) + 0.8, #This is length of x arrow
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) -1, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) -2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 2 ,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 1, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2.5, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 0.5, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis

p6 <- p[[6]] + theme_bw() + ggtitle("JAG1") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),   # Remove axis text (tick labels)
        axis.ticks = element_blank(),  # Remove axis ticks
        panel.grid = element_blank(),   # Remove gridlines
    strip.text = element_text(size = 12),
        text = element_text(family = "Helvetica", size = 12, face = "bold") 
  ) +
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) + 0.8, #This is length of x arrow
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) -1, 
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is x arrow >
  annotate("segment", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) -2, 
           xend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2, 
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1, 
           yend = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 2 ,#This is length of y arrow
           arrow = arrow(length = unit(0.2, "cm")), color = "black") + #This is y arrow >
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 1, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) - 1.8, #Vertical shift
           label = "UMAP1", size = 4) + #This is UMAP1 text on x axis
  annotate("text", x = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,1]) - 2.5, #Horizontal shift
           y = min(twins.DC@reductions$wnn.umapDC@cell.embeddings[,2]) + 0.5, #Vertical shift
           label = "UMAP2", size = 4, angle = 90) #This is UMAP2 text on y axis

plot <- (p2 + p4 + p6) / (p1 + p3 + p5)
plot

#####Enriched Motifs######
BiocManager::install("motifmatchr")
library(motifmatchr)
DefaultAssay(twins) <- "ATAC"
Idents(twins) = "PBMC_l2_predicted.id"

# Get a list of motif position weight matrices from the JASPAR database.
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

# Add the Motif object to the assay.
twins <- AddMotifs(
  object = twins,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Scan the DNA sequence of each peak for the presence of each motif.
twins <- RegionStats(
  object = twins,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

# Create a new Mofif object to store the results
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(twins), sep = c(":", "-")),
  pwm = pfm,
  genome = "BSgenome.Hsapiens.UCSC.hg38",
  sep = c(":", "-")
)

motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Compute motif activities using chromvar
twins <- RunChromVAR(
  object = twins,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix = motif.matrix
)

#Finding overrepresented motifs
Idents(twins) <- "PBMC_l2_predicted.id"
GetChromvarActivities <- function(cluster, seurat_aggregate, motif) {
  print(paste0("Finding chromVAR activities for: ", cluster))
  twins <- seurat_aggregate
  DefaultAssay(twins) <- "chromvar"
  dam <- FindMarkers(twins, 
                     ident.1 = cluster,
                     test.use = "LR",
                     latent.vars = "nCount_ATAC",
                     logfc.threshold = 0)
  motifLookup <- rownames(dam)
  motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
  return(cbind(dam, gene = motifNames))
}

idents <- levels(Idents(twins))
list.cluster.dam <- lapply(idents, function(x) GetChromvarActivities(x, seurat_aggregate = twins, motif = motif))
names(list.cluster.dam) <- idents
write.xlsx(list.cluster.dam, file = paste0(workingDir, "/ATAC/chromVAR.celltypes.xlsx"), sheetName = idents, rowNames = T)

library(openxlsx)
library(readxl)
workingDir <- getwd()
file_path <- paste0(workingDir, "/ATAC/chromVAR.celltypes.xlsx")

sheet_names <- excel_sheets(file_path)
print(sheet_names)

# Read each sheet as a list
chromvar_data <- lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})

# Assign names based on sheet names
names(chromvar_data) <- sheet_names

# Check the structure of the data
str(chromvar_data)

Idents(twins) <- "celltype.vaccine"
levels(twins@active.ident)
da_peaks <- FindMarkers(
  object = twins,
  ident.1 = 'pDC_CERVARIX',
  ident.2 = 'pDC_GARDASIL',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.05 & da_peaks$pct.1 > 0.2 & da_peaks$avg_log2FC > 0.5, ])

#OR use prevously identified DARs
top200 <- dar_pdc_top200$coord
notch <- dar_pdc_notch$coord

# test enrichment
enriched.motifs <- FindMotifs(
  object = twins,
  features = notch
)

write.xlsx(enriched.motifs, file = paste0(workingDir, "/ATAC/EnrichedMotifs_cDC2.xlsx"), rowNames = F)

install.packages('ggseqlogo')
library(ggseqlogo)

MotifPlot(
  object = twins,
  motifs = head(rownames(enriched.motifs))
) +
  theme_bw() + 
  theme(#legend.position = "none",
    axis.text.x.top = element_blank(),
    axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),  # X-axis title font and size
    axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"), 
    axis.text.y = element_text(family = "Helvetica", size = 10, color = "black"), # Y-axis title font and size
    axis.text.x = element_text(family = "Helvetica", size = 10, color = "black"),  # X-axis tick labels font and size
    panel.grid = element_blank(),
    legend.text = element_text(family = "Helvetica", size = 10, color = "black"),
    legend.title = element_text(family = "Helvetica", size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold", colour = "black"), # Facet label color
    strip.background = element_rect(fill = "lightgrey", color = "black")  # Change background color of facet titles
  )

#####CoveragePlots#####
Idents(twins.DC) <- "celltype.vaccine"
levels(twins.DC@active.ident)
DefaultAssay(twins.DC) <- "ATAC"

CoveragePlot(
  object = twins.DC,
  region = "NOTCH2",
  #features = "NOTCH2",
  #expression.assay = "SCT",
  extend.upstream = 50,
  extend.downstream = 50,
  links = F,
  idents = c("cDC1_CERVARIX", "cDC1_GARDASIL"))
  #or use group.by

#####GEO Search#####
BiocManager::install("GEOquery")
library(GEOquery)
library(limma)
library(oligo)
library(affy)
library(edgeR)
library(tidyverse)
library(hgu133plus2.db)
library(annotate)
library(org.Hs.eg.db) 
library(pheatmap)
library(fgsea)
library(msigdbr)
library(ggplot2)

#GSE52245
# load series and platform data from GEO
gset <- getGEO("GSE52245", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL13158", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- paste0("XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX",
               "0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0XX0")
sml <- strsplit(gsms, split="")[[1]]

# samples selection within a series
sel <- which(sml != "X")  # eliminate samples marked as "X
sml <- sml[sel]
gset <- gset[ ,sel]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

groups <- c("Day7")
gs <- factor(sml)
levels(gs) <- groups

metadata <- pData(gset)
geo <- as.data.frame(metadata[, c("geo_accession", "source_name_ch1")])
#Import antibody_data manually
metadata <- merge(metadata, antibody_data, by = "geo_accession")
median_IgG <- median(metadata$IgG, na.rm = TRUE)
metadata$Ab <- ifelse(metadata$IgG > median_IgG, "higher", "lower")

pData(gset) <- metadata

group <- metadata$Ab  # Replace "Group" with the relevant column name in your metadata
group <- factor(group) 

design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(group)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(group[1], group[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

write.csv(tT, "DEGs_GSE52245_MeningococcalVaccine.csv")

#Memory B cell genes HPV twins
gene_signature <- c("RPLP1", "RPL41", "RPL23A", "RPL35A", "RPL10", "MTRNR2L12", "RPS12", "RPL7A", 
                    "EEF1B2", "TPST1", "JUNB", "RPS18", "ALG13", "RPL17", "RPS27", "MT-ATP6", "MT-ND5", 
                    "RPS20", "PLCG2", "RPS10", "RPL12", "VAPA", "HVCN1", "RPL13", "OPA1", "HLA-A", 
                    "RPL29", "RPL5", "GPCPD1", "RPS27A", "RPSA", "RPL32", "NACA", "SEPHS1", "RPLP2", 
                    "BTF3", "PEX6", "AFF3", "SHISA8", "RPS15", "CBFA2T3", "CDK7", "NEXN", "RPS25", 
                    "PSPC1", "MT-ND1", "RPL9", "MGLL", "DIPK2A", "MMUT", "RHOQ", "RPL10A", "PMF1", 
                    "MAGOH", "AC087500.1", "RPS6", "GTF2H1", "TUG1", "CARD16", "AC027018.1", "AC026979.2", 
                    "NUFIP2", "LY9", "INTS9", "MBNL3", "GOLGA8B", "YBX1", "CSTB", "RPL26", "RPL28", 
                    "DHTKD1", "ARRB1", "C5orf67", "DARS", "PIGP", "DCXR", "EHMT1", "RPL27A", "INPP4A", 
                    "KBTBD8", "SUPV3L1", "NDUFB1", "PPWD1", "PRR11", "COX16", "SRRT", "MRTFB", "SNRPB", 
                    "HIPK2", "HSPD1", "MT-ND4", "MZT2B", "TSPAN14", "RPS9", "SIPA1L3", "RWDD1", "RDX", 
                    "GRIPAP1", "WAPL", "MT-CYB")
gene_signature_full <- c("RPLP1", "RPL41", "RPL23A", "RPL35A", "RPL10", "MTRNR2L12", "RPS12", "RPL7A", 
                         "EEF1B2", "TPST1", "JUNB", "RPS18", "ALG13", "RPL17", "RPS27", "MT-ATP6", "MT-ND5", 
                         "RPS20", "PLCG2", "RPS10", "RPL12", "VAPA", "HVCN1", "RPL13", "OPA1", "HLA-A", 
                         "RPL29", "RPL5", "GPCPD1", "RPS27A", "RPSA", "RPL32", "NACA", "SEPHS1", "RPLP2", 
                         "BTF3", "PEX6", "AFF3", "SHISA8", "RPS15", "CBFA2T3", "CDK7", "NEXN", "RPS25", 
                         "PSPC1", "MT-ND1", "RPL9", "MGLL", "DIPK2A", "MMUT", "RHOQ", "RPL10A", "PMF1", 
                         "MAGOH", "AC087500.1", "RPS6", "GTF2H1", "TUG1", "CARD16", "AC027018.1", "AC026979.2", 
                         "NUFIP2", "LY9", "INTS9", "MBNL3", "GOLGA8B", "YBX1", "CSTB", "RPL26", "RPL28", 
                         "DHTKD1", "ARRB1", "C5orf67", "DARS", "PIGP", "DCXR", "EHMT1", "RPL27A", "INPP4A", 
                         "KBTBD8", "SUPV3L1", "NDUFB1", "PPWD1", "PRR11", "COX16", "SRRT", "MRTFB", "SNRPB", 
                         "HIPK2", "HSPD1", "MT-ND4", "MZT2B", "TSPAN14", "RPS9", "SIPA1L3", "RWDD1", "RDX", 
                         "GRIPAP1", "WAPL", "MT-CYB", "RPS15A", "CMC1", "ATOX1", "CARHSP1", "RPL36A", "STK17A", 
                         "CNST", "TOMM6", "PHF11", "RPL34", "KLHL28", "ATXN7", "IFNGR2", "CD53", "RPL30", 
                         "ACBD3", "RPS14", "SIPA1L1", "AC093157.1", "KLHL7", "LRRC28", "CSK", "SH3YL1", "BIRC2", 
                         "RPL15", "AL137009.1", "AC073343.2", "TCAF1", "MED22", "C1QTNF6", "SYNE3", "PSMA1", 
                         "GCNT7", "DUSP12", "F5", "AC010978.1", "TMEM44", "FOXO3B", "AC023421.2", "IQGAP2", 
                         "FNIP1", "SUGP2", "MRPL32", "SEC14L1", "AC092747.4", "COQ5", "DPH1", "PTP4A2", "ECE1", 
                         "TANK", "IL4R", "MED27", "YIPF4", "LRRC37A3", "RPS7", "GSAP", "RPL11", "TIAL1", 
                         "ZNF737", "MIR4435-2HG", "CMTR1", "FRYL", "RPS23", "PDHX", "TUBGCP6", "TFCP2", "SEL1L3", 
                         "WNT5B", "CDIP1", "RAD51C", "LY86-AS1", "BTK", "TRIP12", "LMNB1", "SCARB1", "TNIP1", 
                         "AC132153.1", "ALKBH8", "LINC00624", "PCNA", "RIC8A", "ACTR8", "AC025171.2", "MRPS18A", 
                         "LINC00265", "HSD17B10", "TUT7", "LYSMD3", "RPS27L", "GCA", "GALNT2", "STX5", "RPL37", 
                         "EFCAB13", "ST6GALNAC3", "MAGED1", "KLF8", "RAD17", "ANO10", "RPL36AL", "FARP2", "MAP3K7", 
                         "BCL2", "AFF1", "LAMP2", "PTPRC", "JAZF1", "COQ8A", "TRIM52", "MDH2", "GABPB1-IT1", 
                         "NLRP1", "CD40", "RB1", "MAP3K4", "CDKN1B", "COX17", "CHD1", "RPS11", "AC083837.1", 
                         "BLK", "MRPL42", "DNAJC7", "COX7C", "SLC1A4", "LINC01138", "NXPE3", "MICAL2", "P4HA3", 
                         "RPH3AL", "MVP", "ZNF281", "AC105052.4", "POLK", "MRPL54", "SLC35B2", "TUBE1", "ABCC10", 
                         "ZNF32", "AL049840.2", "VAT1", "DGKE", "CHD2", "LINC01876", "C2CD2L", "HELLPAR", "THRA1/BTR", 
                         "AC009961.1", "SLC43A3", "KREMEN2", "MYPOP", "DBP", "ZNF444", "TCFL5", "DFFA", "IKBKE", 
                         "KLHL12", "POLR3E", "RPS2", "PNRC1", "WDR25", "RAB11B", "MTOR", "ATF7IP", "PAFAH1B1", 
                         "SPCS2", "ZBTB20-AS5", "THEMIS2", "OSTC", "NT5DC3", "TMED2", "JSRP1", "LZTR1", "PIGF", 
                         "RPL18", "AP002075.1", "FAM160B1", "STK38L", "RPL35", "ODF2L", "CTSH", "OPTN", "IFITM2", 
                         "RAB27A", "MYO18A", "BUB3", "MAT2A", "CAP1", "FKBP15", "ARPC5L", "MAPKAPK5-AS1", "FBXO33", 
                         "RPL27", "SGCE", "BBS2", "RPS21", "TUBB", "CHORDC1", "RPL19", "PDLIM1", "POLR2B", "HEXD", 
                         "LRCH3", "SMDT1", "HELB", "SLC35E2A", "TMCO1", "FAM185A", "DOP1A", "KLHL36", "AC025442.2", 
                         "PUM1", "TOMM7", "ZNF800", "IGHM", "UTRN", "RARS2", "COMMD6", "AHCTF1", "SYPL1", "TMTC3", 
                         "CBFB", "EEA1", "CLYBL", "FOXK2", "RPS19", "VPS9D1", "TMEM50A", "METTL6", "SSBP2", "MNAT1", 
                         "QRICH1", "RTTN", "FOCAD", "DNM3", "USP36", "CD84", "MAP4K3", "AC073111.4", "AC044849.1", 
                         "AL390774.2", "AC010320.2", "SHKBP1", "TIMM10", "NOP10", "CARM1", "ZNF701", "RDH14", 
                         "ZNF713", "DEXI", "AL354740.1", "CCNB1IP1", "TXNL4B", "TIMM50", "LINC01266", "TPRA1", 
                         "OSMR", "RUNDC3B", "FAM160B2", "PDCD1LG2", "LRRC27", "FANCF", "AVIL", "DTD2", "GLRX5", 
                         "TMEM266", "AC022960.1", "CYTOR", "ADPGK", "MIB", "TMEM131L", "TLE1", "CRY1", "SNHG32", "NIPSNAP3B", "HMGA1", "OSBPL2", "TMEM208", 
                         "ERBIN", "ANKRD28", "PSMA4", "NCOA2", "ARL5A", "ZFAT", "FAF2", "BTN3A2", "MOSMO", 
                         "VPREB3", "UBE2E3", "CPEB4", "UBE2E2", "ILF3", "AL136962.1", "TRMT10C", "NT5E", 
                         "AC024084.1", "SMG9", "SYT11", "C5orf22", "POLH", "DDX51", "AC011939.2", "SNHG17", 
                         "LDLRAD4", "DENND6A", "TTC19", "SLK", "MPP7", "STAG1", "ATP5IF1", "TEX10", "LAMC1", 
                         "CDKL3", "MBTPS1", "ATF4", "TMEM123", "MTERF4", "FAM3C", "ZNF76", "AL592183.1", 
                         "WWC3", "SAMD9", "CDC37L1", "MTG2", "EIF5", "TIMM44", "LINC01934", "DSTYK", "BCL10", 
                         "SYNGR2", "JMY", "GOLPH3L", "RC3H2", "SKI", "SERPINB1", "AC246817.1", "LY75", "CSAD", 
                         "KIAA1958", "ARFIP1", "DERL1", "MPRIP", "LMTK2", "UCKL1", "NUBPL", "PLEKHO1", "EAF2", 
                         "TTC39B", "CCDC141", "POU2F1", "AC108863.2", "NOTCH1", "AC087239.1", "FAM234A", 
                         "AC090772.2", "USPL1", "EVI2B", "NADK2", "SUPT16H", "AARS", "POGK", "SLC2A11", "MT-CO2", 
                         "MICU3", "NHEJ1", "MT-CO1", "FDX1", "KIAA1143", "PDE4DIP", "FXYD5", "MED23", "SOCS4", 
                         "THYN1", "KCTD12", "AL109628.2", "LAGE3", "DNAJC2", "RPS4X", "MFN2", "AC124016.1", 
                         "CABP4", "MAP3K12", "AC087190.1", "FAM209A", "ARFRP1", "RGL4", "YAF2", "ARRDC2", 
                         "CD99L2", "EEF2", "CD48", "TRIR", "CAMK4", "NFATC2IP", "PTER", "B4GALT4", "AC078845.1", 
                         "QKI", "VDAC2", "CARMIL1", "ETAA1", "MMADHC", "NDUFS3", "TMEM120B", "LRRC59", "AAR2", 
                         "UBE2G1", "SUB1", "PNRC2", "SMARCD1", "POLR3F", "EFCAB2", "CLCN7", "TRDC", "SPIN1", 
                         "AC128687.3", "POLE4", "DHFR2", "CPM", "ISCA2", "SPATA21", "LINC01355", "MTX1", 
                         "AC110769.2", "FLACC1", "CAPN10", "PDIA5", "SEPSECS-AS1", "RELL2", "CYP3A5", "FZD3", 
                         "NDUFA8", "ZBTB34", "TMEM141", "FUOM", "TNFRSF1A", "DNAAF4", "MCRIP2", "AC015853.3", 
                         "AC090617.4", "SYCE2", "AC022098.1", "ZNF628", "VPS16", "SLC12A2", "PKNOX1", "RNF169", 
                         "MUC20", "AC246817.2", "SFXN2", "LYRM9", "AL022328.4", "AL445685.1", "CACNA1E", 
                         "AL512306.3", "CMPK2", "AC130814.1", "MFSD2B", "SLC8A1", "MARS2", "RNF25", "ZNF619", 
                         "AC108519.1", "FNIP2", "HIST1H4H", "AL354892.1", "GRID2IP", "AC000065.1", "TRPA1", 
                         "CASC9", "AC011773.4", "DEPTOR", "AL157938.3", "ADO", "AC127035.1", "CARNS1", "HEPHL1", 
                         "LOH12CR2", "AL356966.1", "C13orf46", "LTB4R", "LINC01146", "AL163932.1", "RPUSD2", 
                         "CELF6", "ZSCAN2", "LYSMD4", "AC020765.2", "LLGL1", "AC018628.2", "AC011933.2", "UHRF1", 
                         "ZNF576", "ZNF616", "AC010320.3", "SLC13A3", "SRRD", "AL031595.2", "SYP", "FAM226B", 
                         "AC244090.1", "ICE1", "AL031428.1", "YARS", "COLGALT1", "WDTC1", "SELENOF", "CPPED1", 
                         "WDR43", "STAU1", "UHRF1BP1", "PHF14", "PDS5B", "SDAD1", "LPAR6", "CERT1", "AC090945.1", 
                         "SERPINI1", "MED19", "MRPL45", "CD2AP", "ZNF429", "TRIM26", "RBM15-AS1", "AC015971.1", 
                         "AC013264.1", "NME6", "PIGZ", "AC083862.3", "TTC26", "TMEM203", "FIBP", "DVL2", "NACC1", 
                         "CDK16", "MAN1B1", "ZNF350", "PPP1R21", "SLC16A10", "SP140L", "CCDC126", "PHB", "GRSF1", 
                         "SRA1", "SHMT2", "TRMU", "EIF3I", "ARF3", "ZNF587", "KLHL5", "RBFOX2", "AKIRIN2", "RPL14", 
                         "NDUFV2", "PILRB", "RBM5", "TBCA", "SCAP", "TSC1", "JAM3", "FUBP3", "ATXN3", "UXS1", 
                         "TXNDC12", "CEP126", "GAS5", "RFC3", "AQR", "STX7", "CXXC1", "RPL22", "PRKN", "COMMD2", 
                         "POLR2J3", "TBC1D1", "RIPOR2", "ZDHHC2", "CNPPD1", "CNOT6L", "C22orf34", "AC008875.3", 
                         "DNASE1L1", "CPQ", "SSH2", "LAP3", "HPS3", "ZBTB7A", "COMMD9", "ATP11C", "PHTF1", 
                         "RPS26", "CTSO", "TRIP11", "CSNK1G1", "VPS45", "PLEKHM1", "RPL23", "BTD", "ZEB2", 
                         "RNF146", "PPHLN1", "CHMP1B", "ANKRD6", "AFDN", "BNIP3L", "DDX60", "TSTD2", "ZNF765", 
                         "STX18-AS1", "CNOT9", "NAF1", "GTF2H2", "FIS1", "MT-CO3", "ABHD10", "AC061958.1", 
                         "PDCD7", "MARCH7", "TRIM59", "YWHAB")

#Visualize DEGs
deg_results2 <- tT

# Filter DEGs (adjusted p-value < 0.05 and |logFC| > 1)
significant_degs <- deg_results2[deg_results2$Gene.symbol %in% gene_signature, ]
#OR
significant_degs2 <- deg_results2[deg_results2$P.Value < 0.05 & abs(deg_results2$logFC) > 1, ]

duplicated_genes <- deg_results2$Gene.symbol[duplicated(deg_results2$Gene.symbol)]
unique_duplicated_genes <- unique(duplicated_genes)

# Extract unique genes in significant_degs
unique_significant_genes <- unique(significant_degs$Gene.symbol)

# Find genes that are in significant_degs but not in gene_signature
extra_genes <- setdiff(unique_significant_genes, gene_signature)

# Remove duplicate rows based on Gene.symbol
deg_results2_unique <- deg_results2[!duplicated(deg_results2$Gene.symbol), ]

# Filter again with unique data
significant_degs <- deg_results2_unique[deg_results2_unique$Gene.symbol %in% gene_signature, ]

#significant_degs2 <- significant_degs[significant_degs$P.Value < 0.05, ]

# Extract expression levels for significant DEGs
#significant_expr <- exprs(gset)[rownames(exprs(gset)) %in% rownames(significant_degs), ]
significant_expr <- exprs(gset)[fData(gset)$Gene.symbol %in% significant_degs2$Gene.symbol, ]

# Ensure unique Gene.symbol in both datasets
unique_genes_in_degs <- unique(significant_degs2$Gene.symbol)  # Unique significant genes

# Filter the expression matrix using unique matching Gene.symbol
matching_rows <- fData(gset)$Gene.symbol %in% unique_genes_in_degs  # Logical vector for matches
significant_expr <- exprs(gset)[matching_rows, ]

# Update rownames with unique Gene.symbols for clarity
rownames(significant_expr) <- fData(gset)$Gene.symbol[matching_rows]

# Remove duplicate gene symbols, keeping the first occurrence
significant_expr <- significant_expr[!duplicated(rownames(significant_expr)), ]

annotation_col <- data.frame(Ab = pData(gset)$Ab, treatment = pData(gset)$`treatment:ch1`, IgG = pData(gset)$IgG)

# Ensure row names in the annotation match the column names in the heatmap
rownames(annotation_col) <- colnames(significant_expr)

pheatmap(
  significant_expr,
  cluster_rows = TRUE,         # Cluster rows (genes)
  cluster_cols = TRUE,         # Cluster columns (samples)
  scale = "row",               # Normalize rows to z-scores
  show_rownames = TRUE,        # Show gene names
  show_colnames = TRUE,        # Show sample names
  annotation_col = annotation_col,
  cutree_rows = 2,
  cutree_cols = 2,
  main = "", # Title of the heatmap
  color = colorRampPalette(c("blue", "white", "red"))(50) # Color scheme
)

#GSEA
# Create a ranked gene list (assuming tT is your limma result table)
#Do not select significant genes, use all
geneList <- tT$logFC
names(geneList) <- tT$Gene.symbol  # Replace with the actual gene symbols or Entrez IDs
geneList <- sort(geneList, decreasing = TRUE)

# Fetch MSigDB gene sets (C2 category for example)
geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Prepare the gene sets
geneSets <- geneSets %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# Run GSEA with fgsea
gsea_result <- fgsea(pathways = geneSets, 
                     stats = geneList, 
                     minSize = 15, 
                     maxSize = 500, 
                     nperm = 1000)

# View the results
head(gsea_result)
gsea_result_clean <- gsea_result %>%
  mutate(across(where(is.list), ~ sapply(., toString)))
write.csv(gsea_result_clean, "GSEA_GSE52245_MeningococcalVaccine.csv")

# Filter for significant results based on adjusted p-value (e.g., < 0.05)
significant_results <- gsea_result %>%
  filter((padj < 0.05 & NES > 2.5) | (pval < 0.45 & NES < -1))

significant_results$pathway <- factor(significant_results$pathway, 
                                      levels = significant_results$pathway[order(significant_results$NES, decreasing = TRUE)])

# Create a bar chart
plot <- ggplot(significant_results, aes(x = NES, y = pathway, fill = padj < 0.05)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey", "red")) +  # Gray for non-significant, Red for significant
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway", 
       fill = "Significance",
       title = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5))

plot
ggsave(filename = "GSEA_GSE52245_MeningococcalVaccine.pdf", plot = plot, width = 50, height = 20, dpi = 320, units = "cm", device = "pdf")

#GSE169159
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE169159", "file=GSE169159_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

gsms <- paste0("XX0XXXXX0XXXXX0XXXX0XXXXXXXXXX0XXXXX0XXXXXXXXXXX0X",
               "XXXX0XXXX0XXXXX0XXXXX0XXXXXX0XXXX0XXXXX0XXXXX0XXX0",
               "XXXX0XXXXX0XXXXX0XXXXX0XXXXX0XXXXX0XXXX0XXXXXX0XXX",
               "XX0XXXXXXXXXX0XXXXX0XXXXX0XXXXX0XXX")  # sample selection
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)

#Import antibody data manually
antibody_data_mRNA <- antibody_data_mRNA %>%
  mutate(
    breadth = rowSums(dplyr::select(., WA:RSA), na.rm = TRUE)
  )

summary(antibody_data_mRNA$breadth)

antibody_data_mRNA <- antibody_data_mRNA %>%
  mutate(
    breadth2 = ifelse(breadth > median(breadth), "higher", "lower")
  )

antibody_data_mRNA$breadth2 <- as.factor(antibody_data_mRNA$breadth2)
levels(antibody_data_mRNA$breadth2)
summary(antibody_data_mRNA$breadth2)

# Check if metadata SampleID matches column names of tbl
all(antibody_data_mRNA$geo_accession %in% colnames(dat))

# Reorder metadata to match count matrix columns
antibody_data_mRNA <- antibody_data_mRNA[match(colnames(tbl), antibody_data_mRNA$geo_accession), ]
all(colnames(dat) == antibody_data_mRNA$geo_accession)

# Check the row names of `antibody_data_mRNA` and column names of `dat`
rownames(antibody_data_mRNA) <- antibody_data_mRNA$geo_accession

# Verify alignment
all(rownames(antibody_data_mRNA) == colnames(dat)) 

group <- antibody_data_mRNA$breadth2  # Replace "Group" with the relevant column name in your metadata
group <- factor(group) 

library(Biobase)
gset <- ExpressionSet(assayData = as.matrix(dat), phenoData = AnnotatedDataFrame(antibody_data_mRNA))

# Pre-filter the genes (e.g., remove rows with all zeros)
gset <- gset[rowSums(exprs(gset)) > 0, ]  # Keep genes with at least one count

validObject(gset)

design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(group)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(levels(group)[1], levels(group)[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# Merge the annotation data (gene symbols) with tT by the GeneID
tT$Gene.symbol <- annot[rownames(tT), "Symbol"]

# If you want to also include the Gene title, you can add it similarly:
tT$Gene.title <- annot[rownames(tT), "Description"]

write.csv(tT, "DEGs_GSE169159_mRNACOVIDVaccine.csv")

#Memory B cell genes HPV twins
gene_signature <- c("RPLP1", "RPL41", "RPL23A", "RPL35A", "RPL10", "MTRNR2L12", "RPS12", "RPL7A", 
                    "EEF1B2", "TPST1", "JUNB", "RPS18", "ALG13", "RPL17", "RPS27", "MT-ATP6", "MT-ND5", 
                    "RPS20", "PLCG2", "RPS10", "RPL12", "VAPA", "HVCN1", "RPL13", "OPA1", "HLA-A", 
                    "RPL29", "RPL5", "GPCPD1", "RPS27A", "RPSA", "RPL32", "NACA", "SEPHS1", "RPLP2", 
                    "BTF3", "PEX6", "AFF3", "SHISA8", "RPS15", "CBFA2T3", "CDK7", "NEXN", "RPS25", 
                    "PSPC1", "MT-ND1", "RPL9", "MGLL", "DIPK2A", "MMUT", "RHOQ", "RPL10A", "PMF1", 
                    "MAGOH", "AC087500.1", "RPS6", "GTF2H1", "TUG1", "CARD16", "AC027018.1", "AC026979.2", 
                    "NUFIP2", "LY9", "INTS9", "MBNL3", "GOLGA8B", "YBX1", "CSTB", "RPL26", "RPL28", 
                    "DHTKD1", "ARRB1", "C5orf67", "DARS", "PIGP", "DCXR", "EHMT1", "RPL27A", "INPP4A", 
                    "KBTBD8", "SUPV3L1", "NDUFB1", "PPWD1", "PRR11", "COX16", "SRRT", "MRTFB", "SNRPB", 
                    "HIPK2", "HSPD1", "MT-ND4", "MZT2B", "TSPAN14", "RPS9", "SIPA1L3", "RWDD1", "RDX", 
                    "GRIPAP1", "WAPL", "MT-CYB")
gene_signature_full <- c("RPLP1", "RPL41", "RPL23A", "RPL35A", "RPL10", "MTRNR2L12", "RPS12", "RPL7A", 
                         "EEF1B2", "TPST1", "JUNB", "RPS18", "ALG13", "RPL17", "RPS27", "MT-ATP6", "MT-ND5", 
                         "RPS20", "PLCG2", "RPS10", "RPL12", "VAPA", "HVCN1", "RPL13", "OPA1", "HLA-A", 
                         "RPL29", "RPL5", "GPCPD1", "RPS27A", "RPSA", "RPL32", "NACA", "SEPHS1", "RPLP2", 
                         "BTF3", "PEX6", "AFF3", "SHISA8", "RPS15", "CBFA2T3", "CDK7", "NEXN", "RPS25", 
                         "PSPC1", "MT-ND1", "RPL9", "MGLL", "DIPK2A", "MMUT", "RHOQ", "RPL10A", "PMF1", 
                         "MAGOH", "AC087500.1", "RPS6", "GTF2H1", "TUG1", "CARD16", "AC027018.1", "AC026979.2", 
                         "NUFIP2", "LY9", "INTS9", "MBNL3", "GOLGA8B", "YBX1", "CSTB", "RPL26", "RPL28", 
                         "DHTKD1", "ARRB1", "C5orf67", "DARS", "PIGP", "DCXR", "EHMT1", "RPL27A", "INPP4A", 
                         "KBTBD8", "SUPV3L1", "NDUFB1", "PPWD1", "PRR11", "COX16", "SRRT", "MRTFB", "SNRPB", 
                         "HIPK2", "HSPD1", "MT-ND4", "MZT2B", "TSPAN14", "RPS9", "SIPA1L3", "RWDD1", "RDX", 
                         "GRIPAP1", "WAPL", "MT-CYB", "RPS15A", "CMC1", "ATOX1", "CARHSP1", "RPL36A", "STK17A", 
                         "CNST", "TOMM6", "PHF11", "RPL34", "KLHL28", "ATXN7", "IFNGR2", "CD53", "RPL30", 
                         "ACBD3", "RPS14", "SIPA1L1", "AC093157.1", "KLHL7", "LRRC28", "CSK", "SH3YL1", "BIRC2", 
                         "RPL15", "AL137009.1", "AC073343.2", "TCAF1", "MED22", "C1QTNF6", "SYNE3", "PSMA1", 
                         "GCNT7", "DUSP12", "F5", "AC010978.1", "TMEM44", "FOXO3B", "AC023421.2", "IQGAP2", 
                         "FNIP1", "SUGP2", "MRPL32", "SEC14L1", "AC092747.4", "COQ5", "DPH1", "PTP4A2", "ECE1", 
                         "TANK", "IL4R", "MED27", "YIPF4", "LRRC37A3", "RPS7", "GSAP", "RPL11", "TIAL1", 
                         "ZNF737", "MIR4435-2HG", "CMTR1", "FRYL", "RPS23", "PDHX", "TUBGCP6", "TFCP2", "SEL1L3", 
                         "WNT5B", "CDIP1", "RAD51C", "LY86-AS1", "BTK", "TRIP12", "LMNB1", "SCARB1", "TNIP1", 
                         "AC132153.1", "ALKBH8", "LINC00624", "PCNA", "RIC8A", "ACTR8", "AC025171.2", "MRPS18A", 
                         "LINC00265", "HSD17B10", "TUT7", "LYSMD3", "RPS27L", "GCA", "GALNT2", "STX5", "RPL37", 
                         "EFCAB13", "ST6GALNAC3", "MAGED1", "KLF8", "RAD17", "ANO10", "RPL36AL", "FARP2", "MAP3K7", 
                         "BCL2", "AFF1", "LAMP2", "PTPRC", "JAZF1", "COQ8A", "TRIM52", "MDH2", "GABPB1-IT1", 
                         "NLRP1", "CD40", "RB1", "MAP3K4", "CDKN1B", "COX17", "CHD1", "RPS11", "AC083837.1", 
                         "BLK", "MRPL42", "DNAJC7", "COX7C", "SLC1A4", "LINC01138", "NXPE3", "MICAL2", "P4HA3", 
                         "RPH3AL", "MVP", "ZNF281", "AC105052.4", "POLK", "MRPL54", "SLC35B2", "TUBE1", "ABCC10", 
                         "ZNF32", "AL049840.2", "VAT1", "DGKE", "CHD2", "LINC01876", "C2CD2L", "HELLPAR", "THRA1/BTR", 
                         "AC009961.1", "SLC43A3", "KREMEN2", "MYPOP", "DBP", "ZNF444", "TCFL5", "DFFA", "IKBKE", 
                         "KLHL12", "POLR3E", "RPS2", "PNRC1", "WDR25", "RAB11B", "MTOR", "ATF7IP", "PAFAH1B1", 
                         "SPCS2", "ZBTB20-AS5", "THEMIS2", "OSTC", "NT5DC3", "TMED2", "JSRP1", "LZTR1", "PIGF", 
                         "RPL18", "AP002075.1", "FAM160B1", "STK38L", "RPL35", "ODF2L", "CTSH", "OPTN", "IFITM2", 
                         "RAB27A", "MYO18A", "BUB3", "MAT2A", "CAP1", "FKBP15", "ARPC5L", "MAPKAPK5-AS1", "FBXO33", 
                         "RPL27", "SGCE", "BBS2", "RPS21", "TUBB", "CHORDC1", "RPL19", "PDLIM1", "POLR2B", "HEXD", 
                         "LRCH3", "SMDT1", "HELB", "SLC35E2A", "TMCO1", "FAM185A", "DOP1A", "KLHL36", "AC025442.2", 
                         "PUM1", "TOMM7", "ZNF800", "IGHM", "UTRN", "RARS2", "COMMD6", "AHCTF1", "SYPL1", "TMTC3", 
                         "CBFB", "EEA1", "CLYBL", "FOXK2", "RPS19", "VPS9D1", "TMEM50A", "METTL6", "SSBP2", "MNAT1", 
                         "QRICH1", "RTTN", "FOCAD", "DNM3", "USP36", "CD84", "MAP4K3", "AC073111.4", "AC044849.1", 
                         "AL390774.2", "AC010320.2", "SHKBP1", "TIMM10", "NOP10", "CARM1", "ZNF701", "RDH14", 
                         "ZNF713", "DEXI", "AL354740.1", "CCNB1IP1", "TXNL4B", "TIMM50", "LINC01266", "TPRA1", 
                         "OSMR", "RUNDC3B", "FAM160B2", "PDCD1LG2", "LRRC27", "FANCF", "AVIL", "DTD2", "GLRX5", 
                         "TMEM266", "AC022960.1", "CYTOR", "ADPGK", "MIB", "TMEM131L", "TLE1", "CRY1", "SNHG32", "NIPSNAP3B", "HMGA1", "OSBPL2", "TMEM208", 
                         "ERBIN", "ANKRD28", "PSMA4", "NCOA2", "ARL5A", "ZFAT", "FAF2", "BTN3A2", "MOSMO", 
                         "VPREB3", "UBE2E3", "CPEB4", "UBE2E2", "ILF3", "AL136962.1", "TRMT10C", "NT5E", 
                         "AC024084.1", "SMG9", "SYT11", "C5orf22", "POLH", "DDX51", "AC011939.2", "SNHG17", 
                         "LDLRAD4", "DENND6A", "TTC19", "SLK", "MPP7", "STAG1", "ATP5IF1", "TEX10", "LAMC1", 
                         "CDKL3", "MBTPS1", "ATF4", "TMEM123", "MTERF4", "FAM3C", "ZNF76", "AL592183.1", 
                         "WWC3", "SAMD9", "CDC37L1", "MTG2", "EIF5", "TIMM44", "LINC01934", "DSTYK", "BCL10", 
                         "SYNGR2", "JMY", "GOLPH3L", "RC3H2", "SKI", "SERPINB1", "AC246817.1", "LY75", "CSAD", 
                         "KIAA1958", "ARFIP1", "DERL1", "MPRIP", "LMTK2", "UCKL1", "NUBPL", "PLEKHO1", "EAF2", 
                         "TTC39B", "CCDC141", "POU2F1", "AC108863.2", "NOTCH1", "AC087239.1", "FAM234A", 
                         "AC090772.2", "USPL1", "EVI2B", "NADK2", "SUPT16H", "AARS", "POGK", "SLC2A11", "MT-CO2", 
                         "MICU3", "NHEJ1", "MT-CO1", "FDX1", "KIAA1143", "PDE4DIP", "FXYD5", "MED23", "SOCS4", 
                         "THYN1", "KCTD12", "AL109628.2", "LAGE3", "DNAJC2", "RPS4X", "MFN2", "AC124016.1", 
                         "CABP4", "MAP3K12", "AC087190.1", "FAM209A", "ARFRP1", "RGL4", "YAF2", "ARRDC2", 
                         "CD99L2", "EEF2", "CD48", "TRIR", "CAMK4", "NFATC2IP", "PTER", "B4GALT4", "AC078845.1", 
                         "QKI", "VDAC2", "CARMIL1", "ETAA1", "MMADHC", "NDUFS3", "TMEM120B", "LRRC59", "AAR2", 
                         "UBE2G1", "SUB1", "PNRC2", "SMARCD1", "POLR3F", "EFCAB2", "CLCN7", "TRDC", "SPIN1", 
                         "AC128687.3", "POLE4", "DHFR2", "CPM", "ISCA2", "SPATA21", "LINC01355", "MTX1", 
                         "AC110769.2", "FLACC1", "CAPN10", "PDIA5", "SEPSECS-AS1", "RELL2", "CYP3A5", "FZD3", 
                         "NDUFA8", "ZBTB34", "TMEM141", "FUOM", "TNFRSF1A", "DNAAF4", "MCRIP2", "AC015853.3", 
                         "AC090617.4", "SYCE2", "AC022098.1", "ZNF628", "VPS16", "SLC12A2", "PKNOX1", "RNF169", 
                         "MUC20", "AC246817.2", "SFXN2", "LYRM9", "AL022328.4", "AL445685.1", "CACNA1E", 
                         "AL512306.3", "CMPK2", "AC130814.1", "MFSD2B", "SLC8A1", "MARS2", "RNF25", "ZNF619", 
                         "AC108519.1", "FNIP2", "HIST1H4H", "AL354892.1", "GRID2IP", "AC000065.1", "TRPA1", 
                         "CASC9", "AC011773.4", "DEPTOR", "AL157938.3", "ADO", "AC127035.1", "CARNS1", "HEPHL1", 
                         "LOH12CR2", "AL356966.1", "C13orf46", "LTB4R", "LINC01146", "AL163932.1", "RPUSD2", 
                         "CELF6", "ZSCAN2", "LYSMD4", "AC020765.2", "LLGL1", "AC018628.2", "AC011933.2", "UHRF1", 
                         "ZNF576", "ZNF616", "AC010320.3", "SLC13A3", "SRRD", "AL031595.2", "SYP", "FAM226B", 
                         "AC244090.1", "ICE1", "AL031428.1", "YARS", "COLGALT1", "WDTC1", "SELENOF", "CPPED1", 
                         "WDR43", "STAU1", "UHRF1BP1", "PHF14", "PDS5B", "SDAD1", "LPAR6", "CERT1", "AC090945.1", 
                         "SERPINI1", "MED19", "MRPL45", "CD2AP", "ZNF429", "TRIM26", "RBM15-AS1", "AC015971.1", 
                         "AC013264.1", "NME6", "PIGZ", "AC083862.3", "TTC26", "TMEM203", "FIBP", "DVL2", "NACC1", 
                         "CDK16", "MAN1B1", "ZNF350", "PPP1R21", "SLC16A10", "SP140L", "CCDC126", "PHB", "GRSF1", 
                         "SRA1", "SHMT2", "TRMU", "EIF3I", "ARF3", "ZNF587", "KLHL5", "RBFOX2", "AKIRIN2", "RPL14", 
                         "NDUFV2", "PILRB", "RBM5", "TBCA", "SCAP", "TSC1", "JAM3", "FUBP3", "ATXN3", "UXS1", 
                         "TXNDC12", "CEP126", "GAS5", "RFC3", "AQR", "STX7", "CXXC1", "RPL22", "PRKN", "COMMD2", 
                         "POLR2J3", "TBC1D1", "RIPOR2", "ZDHHC2", "CNPPD1", "CNOT6L", "C22orf34", "AC008875.3", 
                         "DNASE1L1", "CPQ", "SSH2", "LAP3", "HPS3", "ZBTB7A", "COMMD9", "ATP11C", "PHTF1", 
                         "RPS26", "CTSO", "TRIP11", "CSNK1G1", "VPS45", "PLEKHM1", "RPL23", "BTD", "ZEB2", 
                         "RNF146", "PPHLN1", "CHMP1B", "ANKRD6", "AFDN", "BNIP3L", "DDX60", "TSTD2", "ZNF765", 
                         "STX18-AS1", "CNOT9", "NAF1", "GTF2H2", "FIS1", "MT-CO3", "ABHD10", "AC061958.1", 
                         "PDCD7", "MARCH7", "TRIM59", "YWHAB")

#Visualize DEGs
deg_results2 <- tT

# Filter DEGs (adjusted p-value < 0.05 and |logFC| > 1)
significant_degs <- deg_results2[deg_results2$Gene.symbol %in% gene_signature, ]
#OR
significant_degs2 <- deg_results2[deg_results2$P.Value < 0.05 & abs(deg_results2$logFC) > 1, ]

# Extract expression levels for significant DEGs
significant_expr <- exprs(gset)[rownames(exprs(gset)) %in% rownames(significant_degs), ]
rownames(significant_expr) <- significant_degs$Gene.symbol

annotation_col <- data.frame(Breadth = pData(gset)$breadth2)

# Ensure row names in the annotation match the column names in the heatmap
rownames(annotation_col) <- colnames(significant_expr)

pheatmap(
  significant_expr,
  cluster_rows = TRUE,         # Cluster rows (genes)
  cluster_cols = TRUE,         # Cluster columns (samples)
  scale = "row",               # Normalize rows to z-scores
  show_rownames = TRUE,        # Show gene names
  show_colnames = TRUE,        # Show sample names
  annotation_col = annotation_col,
  cutree_rows = 2,
  cutree_cols = 4,
  main = "", # Title of the heatmap
  color = colorRampPalette(c("blue", "white", "red"))(50) # Color scheme
)

#GSEA
# Create a ranked gene list (assuming tT is your limma result table)
#Do not select significant genes, use all
geneList <- tT$logFC
names(geneList) <- tT$Gene.symbol  # Replace with the actual gene symbols or Entrez IDs
geneList <- sort(geneList, decreasing = TRUE)

# Fetch MSigDB gene sets (C2 category for example)
geneSets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Prepare the gene sets
geneSets <- geneSets %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# Run GSEA with fgsea
gsea_result <- fgsea(pathways = geneSets, 
                     stats = geneList, 
                     minSize = 15, 
                     maxSize = 500, 
                     nperm = 1000)

# View the results
head(gsea_result)
gsea_result_clean <- gsea_result %>%
  mutate(across(where(is.list), ~ sapply(., toString)))
write.csv(gsea_result_clean, "GSEA_GSE169159_mRNACOVIDVaccine.csv")

# Filter for significant results based on adjusted p-value (e.g., < 0.05)
significant_results <- gsea_result %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>%  # Sort by NES in descending order
  slice_head(n = 20)     # Select the top 20 pathways

significant_results$pathway <- factor(significant_results$pathway, 
                                      levels = significant_results$pathway[order(significant_results$NES, decreasing = TRUE)])

# Create a bar chart
plot <- ggplot(significant_results, aes(x = NES, y = pathway)) +
  geom_bar(stat = "identity", fill = "#D6A6CD") +
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway", 
       fill = "Significance",
       title = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5))

plot
ggsave(filename = "GSEA_GSE169159_mRNACOVIDVaccine.pdf", plot = plot, width = 25, height = 20, dpi = 320, units = "cm", device = "pdf")

#Combine
# Load necessary libraries
library(fgsea)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(circlize)

# Step 1: Read both csv files
file1 <- "/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/MSD_Cytokines/MSD Cytokines/GSEA_GSE52245_MeningococcalVaccine.csv"
file2 <- "/Users/valentino/Library/CloudStorage/OneDrive-UGent/Research Projects/HPV Cross Neutra study/PublicData/GSEA_GSE169159_mRNACOVIDVaccine.csv"

fgsea_result1 <- read.csv(file1)
fgsea_result2 <- read.csv(file2)

fgsea_result1$X <- NULL
fgsea_result2$X <- NULL

cell_cycle_translation <- c(
  #CellCycle
  "REACTOME_CYCLIN_D_ASSOCIATED_EVENTS_IN_G1",
  #Cyclin D is crucial for progression through the G1 phase of the cell cycle, which is important for cell proliferation.
  "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_ADDITIONAL_CELL_CYCLE_GENES",
  #TP53 plays a critical role in cell cycle regulation, particularly in response to stress, DNA damage, or immune activation.
  "REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_METHYLATION",
  #This pathway controls p53 activity via methylation, which is important for regulating immune cell proliferation.
  "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B",
  #Involved in mitotic exit and cell cycle progression, this pathway indicates how immune cells divide.
  "REACTOME_CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY",
  #This pathway marks the initiation of DNA replication, another important event in immune cell proliferation.
  "REACTOME_DEVELOPMENTAL_BIOLOGY",
  "REACTOME_MAP2K_AND_MAPK_ACTIVATION",
  "REACTOME_PROLONGED_ERK_ACTIVATION_EVENTS",
  "REACTOME_SIGNALING_BY_NOTCH4",
  "REACTOME_SIGNALING_BY_BRAF_AND_RAF1_FUSIONS",
  #RNA translation
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  #Involves the initiation of mRNA translation, a key step for producing proteins that are important for immune cell function.
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  #Represents the elongation phase of translation, where polypeptide chains are synthesized, crucial for immune activation.
  "REACTOME_RRNA_PROCESSING",
  #The production and modification of ribosomal RNA, necessary for ribosome function and protein synthesis.
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  #This pathway is involved in targeting proteins to the membrane during translation, critical for producing receptors and other membrane-bound proteins essential for immune function.
  "REACTOME_NON_SENSE_MEDIATED_DECAY_NMD",
  #This pathway degrades faulty mRNAs, maintaining the fidelity of protein synthesis during immune responses.
  "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S",
  "REACTOME_DNA_DAMAGE_RECOGNITION_IN_GG_NER",
  "REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL",
  "REACTOME_TRANSLATION",
  "REACTOME_METABOLISM_OF_RNA",
  "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_THE_AP_2_TFAP2_FAMILY_OF_TRANSCRIPTION_FACTORS"
)


# Step 2: Extract significantly enriched pathways (adjust p-value threshold as needed)
sig_pathways1 <- fgsea_result1 %>% filter(fgsea_result1$pathway %in% cell_cycle_translation)
sig_pathways2 <- fgsea_result2 %>% filter(fgsea_result2$pathway %in% cell_cycle_translation)

#OR Filter for significant pathways based on adjusted p-value
#sig_pathways1 <- fgsea_result1 %>% filter(padj < 0.05 & abs(NES) > 2)
#sig_pathways2 <- fgsea_result2 %>% filter(padj < 0.05 & abs(NES) > 2)

# Step 3: Select overlapping pathways between the two sets
#overlapping_pathways <- intersect(sig_pathways1$pathway, sig_pathways2$pathway)

# Extract the NES values for overlapping pathways
#sig_pathways1_overlap <- sig_pathways1 %>% filter(pathway %in% overlapping_pathways)
#sig_pathways2_overlap <- sig_pathways2 %>% filter(pathway %in% overlapping_pathways)

#Heatmap
# Combine the NES values from both results for the overlapping pathways
combined_data <- merge(sig_pathways1[, c("pathway", "NES")],
                       sig_pathways2[, c("pathway", "NES")],
                       by = "pathway", suffixes = c("_file1", "_file2"))

# Step 4: Plot the heatmap with y-axis the selected pathways and colour by NES
# Prepare data for heatmap
heatmap_data <- combined_data[, c("NES_file1", "NES_file2")]
rownames(heatmap_data) <- combined_data$pathway
colnames(heatmap_data)[1] <- "Meningococcal conjugated polysaccharide vaccine"
colnames(heatmap_data)[2] <- "COVID-19 mRNA vaccine"

col_annotation <- data.frame(
  Vaccine = c("Meningococcal conjugated polysaccharide vaccine", 
              "COVID-19 mRNA vaccine")
)
rownames(col_annotation) <- colnames(heatmap_data)

annotation_colors <- list(
  Vaccine = c(
    "Meningococcal conjugated polysaccharide vaccine" = "#B4EC86",
    "COVID-19 mRNA vaccine" = "#FFC5F5"
  )
)

# Plot heatmap
plot <- pheatmap(heatmap_data, 
                 cluster_rows = T, 
                 cluster_cols = F, 
                 scale = "none", 
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 main = "",
                 cellwidth = 25,
                 cellheight = 25,
                 annotation_col = col_annotation,
                 annotation_colors = annotation_colors,
                 labels_col = "",
                 fontsize = 5) 

plot
ggsave(filename = "GSEA_OtherVaccines.pdf", plot = plot, width = 25, height = 20, dpi = 320, units = "cm", device = "pdf")