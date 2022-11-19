if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("GDCRNATools")
project<-'TCGA-PRAD'
rnadir<-paste(project,'RNAseq',sep="/")
mirdir<-paste(project,"miRNAs",sep="/")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-PRAD',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir)
gdcRNADownload(project.id = "TCGA-PRAD",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir)
#Download clinical data
clinicaldir<-paste(project,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-PRAD",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_PRAD<-gdcParseMetadata(project.id = "TCGA-PRAD",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_PRAD<-gdcFilterDuplicate(metaMatrix.RNA_PRAD)
metaMatrix.RNA_PRAD<-gdcFilterSampleType(metaMatrix.RNA_PRAD)

####Parse miRNA metadata
metaMatrix.MIR_PRAD<-gdcParseMetadata(project.id = "TCGA-PRAD",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_PRAD<-gdcFilterDuplicate(metaMatrix.MIR_PRAD)
metaMatrix.MIR_PRAD<-gdcFilterSampleType(metaMatrix.MIR_PRAD)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_PRAD<-gdcRNAMerge(metadata = metaMatrix.RNA_PRAD,path = rnadir,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_PRAD<-gdcRNAMerge(metadata = metaMatrix.MIR_PRAD,path = mirdir,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_PRAD<-gdcClinicalMerge(path = clinicaldir,key.info = TRUE)
ClinicalDa_PRAD[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_PRAD<-gdcVoomNormalization(counts = rnaCounts_PRAD,filter = FALSE)
####Normalization of miRNAs data
mirExpr_PRAD<-gdcVoomNormalization(counts = mirCounts_PRAD,filter = FALSE)
####Differential gene expression analysis
DEGAll_PRAD<-gdcDEAnalysis(counts = rnaCounts_PRAD,group = metaMatrix.RNA$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
DEGMIR_PRAD2<-gdcDEAnalysis(counts = mirCounts_PRAD,group = metaMatrix.MIR_PRAD$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
data("DEGAll")
####All DEGs
deALL_PRAD<-gdcDEReport(deg = DEGAll_PRAD,gene.type = "all")
####DE long-noncoding
deLNC_PRAD<-gdcDEReport(deg = DEGAll_PRAD,gene.type = "long_non_coding")
####DE protein coding genes
dePC_PRAD<-gdcDEReport(deg=DEGAll_PRAD,gene.type = "protein_coding")
dePseudo_PRAD<-gdcDEReport(deg = DEGAll_PRAD,gene.type="pseudogene")
####DEG visualisation
gdcVolcanoPlot(DEGAll_PRAD)
gdcBarPlot(deg = deALL_PRAD,angle = 45,data.type = "RNAseq")
degName=rownames(deALL_PRAD)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_PRAD,rna.expr = rnaExpr_PRAD)

##ceRNAs network analysis
ceOutput_PRAD<- gdcCEAnalysis(lnc= rownames(deLNC_PRAD),pc= rownames(dePC_PRAD),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_PRAD,mir.expr    = mirExpr_PRAD)
ceOutput_PRAD2<- gdcCEAnalysis(lnc= rownames(deLNC_PRAD),pc= rownames(dePC_PRAD),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_PRAD,mir.expr    = mirExpr_PRAD)
ceOutput_PRAD3<- gdcCEAnalysis(lnc= rownames(deLNC_PRAD),pc= rownames(dePC_PRAD),lnc.targets = 'miRcode',pc.targets= 'starBase',rna.expr= rnaExpr_PRAD,mir.expr    = mirExpr_PRAD)
ceOutput2_PRAD3<-ceOutput_PRAD3[ceOutput_PRAD3$hyperPValue<0.01 & ceOutput_PRAD3$cor>0.4,]




####Network Visualization in Cytoscape
ceOutput2_PRAD<-ceOutput_PRAD[ceOutput_PRAD$hyperPValue<0.01 & ceOutput_PRAD$corPValue<0.01 & ceOutput_PRAD
                              $regSim !=0,]
ceOutput2_PRAD2<-ceOutput_PRAD2[ceOutput_PRAD2$hyperPValue<0.01 & ceOutput_PRAD2$corPValue<0.01 & ceOutput_PRAD2$regSim !=0,]
edges_PRAD<-gdcExportNetwork(ceNetwork = ceOutput2_PRAD,net = "edges")
nodes_PRAD<-gdcExportNetwork(ceNetwork = ceOutput2_PRAD,net = "nodes")

write.table(edges_PRAD,file = "edges_PRAD.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_PRAD,file="nodes_PRAD.txt",sep="\t",quote = FALSE,row.names = FALSE)

###Enrichment analysis
library(GDCRNATools)
enrichOutput_PRAD<- gdcEnrichAnalysis(gene = rownames(deALL_PRAD), simplify = TRUE)
gdcEnrichPlot(enrichOutput_PRAD, type = 'bar', category = 'GO', num.terms = 10)
gdcEnrichPlot(enrichOutput_PRAD, type='bubble', category='GO
              ', num.terms = 10)
shinyCorPlot(gene1    = rownames(deLNC_PRAD), 
             gene2    = rownames(dePC_PRAD), 
             rna.expr = rnaExpr_PRAD, 
             metadata = metaMatrix.RNA_PRAD)
####Survival Analysis
#####Cox-PH Model
survOutput_PRAD_COXPH<- gdcSurvivalAnalysis(gene     = rownames(deALL_PRAD), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr_PRAD, 
                                  metadata = metaMatrix.RNA_PRAD)

#####Kaplan Meir Plot
survOutput_PRAD_KM<- gdcSurvivalAnalysis(gene     = rownames(demiR_PRAD), 
                                  method   = 'KM', 
                                  rna.expr = mirExpr_PRAD, 
                                  metadata = metaMatrix.MIR_PRAD, 
                                  sep      = 'median')

gdcKMPlot(gene     = 'hsa-miR-204-5P',
          rna.expr = mirExpr_PRAD,
          metadata = metaMatrix.MIR_PRAD,
          sep      = 'median')

shinyKMPlot(gene = rownames(demiR_PRAD), rna.expr = mirExpr_PRAD, 
            metadata = metaMatrix.MIR_PRAD)









