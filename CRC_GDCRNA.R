project<-'TCGA-CRC'
rnadir<-paste(project,"RNAseq",sep="/")
mirdir<-paste(project,"miRNAs",sep="/")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-CRC',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir)
gdcRNADownload(project.id = "TCGA-CRC",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir)
#Download clinical data
clinicaldir<-paste(project,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-CRC",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_CRC<-gdcParseMetadata(project.id = "TCGA-CRC",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_CRC<-gdcFilterDuplicate(metaMatrix.RNA_CRC)
metaMatrix.RNA_CRC<-gdcFilterSampleType(metaMatrix.RNA_CRC)

####Parse miRNA metadata
metaMatrix.MIR_CRC<-gdcParseMetadata(project.id = "TCGA-CRC",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_CRC<-gdcFilterDuplicate(metaMatrix.MIR_CRC)
metaMatrix.MIR_CRC<-gdcFilterSampleType(metaMatrix.MIR_CRC)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_CRC<-gdcRNAMerge(metadata = metaMatrix.RNA_CRC,path = rnadir,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_CRC<-gdcRNAMerge(metadata = metaMatrix.MIR_CRC,path = mirdir,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_CRC<-gdcClinicalMerge(path = clinicaldir,key.info = TRUE)
ClinicalDa_CRC[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_CRC<-gdcVoomNormalization(counts = rnaCounts_CRC,filter = FALSE)
####Normalization of miRNAs data
mirExpr_CRC<-gdcVoomNormalization(counts = mirCounts_CRC,filter = FALSE)
####Differential gene expression analysis
DEGAll_CRC<-gdcDEAnalysis(counts = rnaCounts_CRC,group = metaMatrix.RNA_CRC$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
data("DEGAll")
####All DEGs
deALL_CRC<-gdcDEReport(deg = DEGAll_CRC,gene.type = "all")
####DE long-noncoding
deLNC_CRC<-gdcDEReport(deg = DEGAll_CRC,gene.type = "long_non_coding")
####DE protein coding genes
dePC_CRC<-gdcDEReport(deg=DEGAll_CRC,gene.type = "protein_coding")

####DEG visualisation
gdcVolcanoPlot(DEGAll_CRC)
gdcBarPlot(deg = deALL_CRC,angle = 45,data.type = "RNAseq")
degName=rownames(deALL_CRC)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_CRC,rna.expr = rnaExpr_CRC)

##ceRNAs network analysis

ceOutput_CRC<- gdcCEAnalysis(lnc= rownames(deLNC_CRC),pc= rownames(dePC_CRC),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_CRC,mir.expr    = mirExpr_CRC)

####Network Visualization in Cytoscape

ceOutput2_CRC<-ceOutput_CRC[ceOutput_CRC$hyperPValue<0.01 & ceOutput_CRC$corPValue<0.01 & ceOutput_CRC$regSim !=0,]
edges_CRC<-gdcExportNetwork(ceNetwork = ceOutput2_CRC,net = "edges")
nodes_CRC<-gdcExportNetwork(ceNetwork = ceOutput2_CRC,net = "nodes")

write.table(edges_CRC,file = "edges_CRC.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_CRC,file="nodes_CRC.txt",sep="\t",quote = FALSE,row.names = FALSE)

head(edges_CRC)

shinyCorPlot(gene1    = rownames(deLNC_CRC), 
             gene2    = rownames(dePC_CRC), 
             rna.expr = rnaExpr_CRC, 
             metadata = metaMatrix.RNA_CRC)
####Survival Analysis
#####Cox-PH Model
survOutput_CRC_COXPH<- gdcSurvivalAnalysis(gene     = rownames(deALL_CRC), 
                                            method   = 'coxph', 
                                            rna.expr = rnaExpr_CRC, 
                                            metadata = metaMatrix.RNA_CRC)

#####Kaplan Meir Plot
survOutput_CRC_KM<- gdcSurvivalAnalysis(gene     = rownames(deALL_CRC), 
                                         method   = 'KM', 
                                         rna.expr = rnaExpr_CRC, 
                                         metadata = metaMatrix.RNA_CRC, 
                                         sep      = 'median')

gdcKMPlot(gene     = 'ENSG00000137331',
          rna.expr = rnaExpr_CRC,
          metadata = metaMatrix.RNA_CRC,
          sep      = 'median')











