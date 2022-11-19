project<-'TCGA-UCEC'
rnadir<-paste(project,"RNAseq",sep="/")
mirdir<-paste(project,"miRNAs",sep="/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("GDCRNATools",lib="C:/Users/n10136142/AppData/Local/Temp/RtmpMhnjtL/downloaded_packages")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-UCEC',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir)
gdcRNADownload(project.id = "TCGA-UCEC",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir)
#Download clinical data
clinicaldir<-paste(project,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-UCEC",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_UCEC<-gdcParseMetadata(project.id = "TCGA-UCEC",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_UCEC<-gdcFilterDuplicate(metaMatrix.RNA_UCEC)
metaMatrix.RNA_UCEC<-gdcFilterSampleType(metaMatrix.RNA_UCEC)

####Parse miRNA metadata
metaMatrix.MIR_UCEC<-gdcParseMetadata(project.id = "TCGA-UCEC",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_UCEC<-gdcFilterDuplicate(metaMatrix.MIR_UCEC)
metaMatrix.MIR_UCEC<-gdcFilterSampleType(metaMatrix.MIR_UCEC)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_UCEC<-gdcRNAMerge(metadata = metaMatrix.RNA_UCEC,path = rnadir,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_UCEC<-gdcRNAMerge(metadata = metaMatrix.MIR_UCEC,path = mirdir,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_UCEC<-gdcClinicalMerge(path = clinicaldir,key.info = TRUE)
ClinicalDa_UCEC[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_UCEC<-gdcVoomNormalization(counts = rnaCounts_UCEC,filter = FALSE)
####Normalization of miRNAs data
mirExpr_UCEC<-gdcVoomNormalization(counts = mirCounts_UCEC,filter = FALSE)
####Differential gene expression analysis
DEGAll_UCEC<-gdcDEAnalysis(counts = rnaCounts_UCEC,group = metaMatrix.RNA_UCEC$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
DEGMIR_UCEC<-gdcDEAnalysis(counts = mirCounts_UCEC,group = metaMatrix.MIR_UCEC$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma")
####As I need log2FC for all set of genes I added a modification for gdcDEAnlaysis
DEGMIR_UCEC2<-gdcDEAnalysis(counts = mirCounts_UCEC,group = metaMatrix.MIR_UCEC$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
setDT(DEGMIR_UCEC2, keep.rownames = TRUE)[]
demiR_UCEC2<-gdcDEReport(DEGMIR_UCEC2,fc=0.05)

data("DEGAll")
####All DEGs
deALL_UCEC<-gdcDEReport(deg = DEGAll_UCEC,gene.type = "all")
####DE long-noncoding
deLNC_UCEC<-gdcDEReport(deg = DEGAll_UCEC,gene.type = "long_non_coding")
####DE protein coding genes
dePC_UCEC<-gdcDEReport(deg=DEGAll_UCEC,gene.type = "protein_coding")
##DE pseudogenes
dePseudo_UCEC<-gdcDEReport(deg=DEGAll_UCEC,gene.type = "pseudogene")
demiR_UCEC<-gdcDEReport(DEGMIR_UCEC)
################################################################################################
DEGALL_UCEC_MIR<-gdcDEAnalysis(counts = mirCounts_UCEC,group = metaMatrix.MIR_UCEC$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
gdcVolcanoPlot(DEGALL_UCEC_MIR)
gdcHeatmap(deg.id = rownames(DEGALL_UCEC_MIR),metadata = metaMatrix.MIR_UCEC,rna.expr = mirExpr_UCEC)
demiR_UCEC<-gdcDEReport(deg=DEGALL_UCEC_MIR,gene.type = "miRNAs")
####DEG visualisation
gdcVolcanoPlot(DEGAll_UCEC)
gdcBarPlot(deg = deALL_UCEC,angle = 45,data.type = "RNAseq")
degName=rownames(deALL_UCEC)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_UCEC,rna.expr = rnaExpr_UCEC)
gdcBarPlot(deg=demiR_UCEC,angle=45,data.type="miRNAs")
t1<-table(which(demiR_UCEC$logFC<0))
dim(t1)
t2<-table(which(demiR_UCEC$logFC>0))
dim(t2)
##ceRNAs network analysis

ceOutput_UCEC<- gdcCEAnalysis(lnc= rownames(deLNC_UCEC),pc= rownames(dePC_UCEC),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_UCEC,mir.expr    = mirExpr_UCEC)
ceOutput_UCEC2<- gdcCEAnalysis(lnc= rownames(deLNC_UCEC),pc= rownames(dePC_UCEC),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_UCEC,mir.expr    = mirExpr_UCEC)
MC_HC_ceOutput_UCEC<-ceOutput_UCEC2[ceOutput_UCEC2$hyperPValue<0.01 & ceOutput_UCEC2$cor>0.4,]
ceOutput_UCEC3<- gdcCEAnalysis(lnc= rownames(deLNC_UCEC),pc= rownames(dePC_UCEC),lnc.targets = 'miRcode',pc.targets= 'starBase',rna.expr= rnaExpr_UCEC,mir.expr    = mirExpr_UCEC)
ceOutput2_UCEC3<-ceOutput_UCEC3[ceOutput_UCEC3$hyperPValue<0.01 & ceOutput_UCEC3$cor>0.4,]
####Network Visualization in Cytoscape

ceOutput2_UCEC<-ceOutput_UCEC[ceOutput_UCEC$hyperPValue<0.01 & ceOutput_UCEC$corPValue<0.01 & ceOutput_UCEC$regSim !=0,]
ceOutput2_UCEC2<-ceOutput_UCEC2[ceOutput_UCEC2$hyperPValue<0.01 & ceOutput_UCEC2$corPValue<0.01 & ceOutput_UCEC2$regSim !=0,]
edges_UCEC<-gdcExportNetwork(ceNetwork = ceOutput2_UCEC,net = "edges")
nodes_UCEC<-gdcExportNetwork(ceNetwork = ceOutput2_UCEC,net = "nodes")

write.table(edges_UCEC,file = "edges_UCEC.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_UCEC,file="nodes_UCEC.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(ceOutput2_UCEC,file="ceOutput_UCEC.txt",sep="\t",quote = FALSE,row.names=FALSE)
head(edges_UCEC)


#########I have sent this predictions to my third aim analysis. The ceOutput_UCEC2 contains miRNA-target predictions for 
##First i will extract gene-miRNA predictions
miRNA_UCEC<-strsplit(ceOutput_UCEC2$miRNAs,split = ",")
mRNA_miRNA_UCEC<-data.frame(V1=rep(ceOutput_UCEC2$Genes,sapply(miRNA_UCEC,length)),V2=unlist(miRNA_UCEC))
lncRNA_miRNA_UCEC<-data.frame(V1=rep(ceOutput_UCEC2$lncRNAs,sapply(miRNA_UCEC,length)),V2=unlist(miRNA_UCEC))

##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
mRNA_UCEC<-mRNA_miRNA_UCEC$V1
symbol_mRNA_UCEC<-ensembldb::select(EnsDb.Hsapiens.v79,keys =mRNA_UCEC,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(mRNA_miRNA_UCEC,symbol_mRNA_UCEC,by=c("V1"="GENEID"))
mRNA_miRNA_UCEC<-data.frame(d1$V2,d1$SYMBOL)
names(mRNA_miRNA_UCEC)<-c("miRNA","target")
View(mRNA_miRNA_UCEC)

###Then, I will extract lncRNA-miRNA predictions
##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
lncRNA_UCEC<-lncRNA_miRNA_UCEC$V1
symbol_lncRNA_UCEC<-ensembldb::select(EnsDb.Hsapiens.v79,keys =lncRNA_UCEC,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(lncRNA_miRNA_UCEC,symbol_lncRNA_UCEC,by=c("V1"="GENEID"))
lncRNA_miRNA_UCEC<-data.frame(d1$V2,d1$SYMBOL)
names(lncRNA_miRNA_UCEC)<-c("miRNA","target")
lncRNA_miRNA_UCEC$target<-paste("lncRNA",lncRNA_miRNA_UCEC$target,sep = "-")
View(lncRNA_miRNA_UCEC)
#Binding mRNA and lncRNA predictions
mRNA_lncRNA_miRNA_UCEC<-rbind(mRNA_miRNA_UCEC,lncRNA_miRNA_UCEC)
miRNA.target.interactions<-mRNA_lncRNA_miRNA_UCEC
write.table(miRNA.target.interactions,file = "Z:/GENETICS/DNA_Methylation/UCEC/miRNA.target.interactions.txt",sep = "\t",quote=FALSE)


