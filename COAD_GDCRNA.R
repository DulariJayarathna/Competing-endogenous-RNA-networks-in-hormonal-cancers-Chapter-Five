project<-'TCGA-COAD'
rnadir<-paste(project,"RNAseq",sep="/")
mirdir<-paste(project,"miRNAs",sep="/")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-COAD',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir)
gdcRNADownload(project.id = "TCGA-COAD",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir)
#Download clinical data
clinicaldir<-paste(project,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-COAD",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_COAD<-gdcParseMetadata(project.id = "TCGA-COAD",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_COAD<-gdcFilterDuplicate(metaMatrix.RNA_COAD)
metaMatrix.RNA_COAD<-gdcFilterSampleType(metaMatrix.RNA_COAD)

####Parse miRNA metadata
metaMatrix.MIR_COAD<-gdcParseMetadata(project.id = "TCGA-COAD",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_COAD<-gdcFilterDuplicate(metaMatrix.MIR_COAD)
metaMatrix.MIR_COAD<-gdcFilterSampleType(metaMatrix.MIR_COAD)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_COAD<-gdcRNAMerge(metadata = metaMatrix.RNA_COAD,path = rnadir,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_COAD<-gdcRNAMerge(metadata = metaMatrix.MIR_COAD,path = mirdir,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_COAD<-gdcClinicalMerge(path = clinicaldir,key.info = TRUE)
ClinicalDa_COAD[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_COAD<-gdcVoomNormalization(counts = rnaCounts_COAD,filter = FALSE)
####Normalization of miRNAs data
mirExpr_COAD<-gdcVoomNormalization(counts = mirCounts_COAD,filter = FALSE)
####Differential gene expression analysis
DEGAll_COAD<-gdcDEAnalysis(counts = rnaCounts_COAD,group = metaMatrix.RNA_COAD$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")

DEGMIR_COAD<-gdcDEAnalysis(counts = mirCounts_COAD,group = metaMatrix.MIR_COAD$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma")
####As I need log2FC for all set of genes I added a modification for gdcDEAnlaysis
DEGMIR_COAD2<-gdcDEAnalysis(counts = mirCounts_COAD,group = metaMatrix.MIR_COAD$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
setDT(DEGMIR_COAD2, keep.rownames = TRUE)[]
demiR_COAD2<-gdcDEReport(DEGMIR_COAD2,fc=0.05)
data("DEGAll")
####All DEGs
deALL_COAD<-gdcDEReport(deg = DEGAll_COAD,gene.type = "all")
####DE long-noncoding
deLNC_COAD<-gdcDEReport(deg = DEGAll_COAD,gene.type = "long_non_coding")
####DE protein coding genes
dePC_COAD<-gdcDEReport(deg=DEGAll_COAD,gene.type = "protein_coding")
##DE pseudogenes
dePseudo_COAD<-gdcDEReport(deg=DEGAll_COAD,gene.type = "pseudogene")
demiR_COAD<-gdcDEReport(DEGMIR_COAD)
################################################################################################
DEGALL_COAD_MIR<-gdcDEAnalysis(counts = mirCounts_COAD,group = metaMatrix.MIR_COAD$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
gdcVolcanoPlot(DEGALL_COAD_MIR)
gdcHeatmap(deg.id = rownames(DEGALL_COAD_MIR),metadata = metaMatrix.MIR_COAD,rna.expr = mirExpr_COAD)
demiR_COAD<-gdcDEReport(deg=DEGALL_COAD_MIR,gene.type = "miRNAs")
####DEG visualisation
gdcVolcanoPlot(DEGAll_COAD)
gdcBarPlot(deg = deALL_COAD,angle = 45,data.type = "RNAseq")
degName=rownames(deALL_COAD)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_COAD,rna.expr = rnaExpr_COAD)
gdcBarPlot(deg=demiR_COAD,angle=45,data.type="miRNAs")
t1<-table(which(demiR_COAD$logFC<0))
dim(t1)
t2<-table(which(demiR_COAD$logFC>0))
dim(t2)
##ceRNAs network analysis

ceOutput_COAD<- gdcCEAnalysis(lnc= rownames(deLNC_COAD),pc= rownames(dePC_COAD),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_COAD,mir.expr    = mirExpr_COAD)
ceOutput_COAD2<- gdcCEAnalysis(lnc= rownames(deLNC_COAD),pc= rownames(dePC_COAD),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_COAD,mir.expr    = mirExpr_COAD)
MC_HC_ceOutput_COAD<-ceOutput_COAD2[ceOutput_COAD2$hyperPValue<0.01 & ceOutput_COAD2$cor>0.4,]

##Network analysis
############################################
ceOutput_COAD3<- gdcCEAnalysis(lnc= rownames(deLNC_COAD),pc= rownames(dePC_COAD),lnc.targets = 'miRcode',pc.targets= 'starBase',rna.expr= rnaExpr_COAD,mir.expr    = mirExpr_COAD)
ceOutput2_COAD3<-ceOutput_COAD3[ceOutput_COAD3$hyperPValue<0.01 & ceOutput_COAD3$cor>0.4,]
####Network Visualization in Cytoscape
ceOutput2_COAD<-ceOutput_COAD[ceOutput_COAD$hyperPValue<0.01 & ceOutput_COAD$corPValue<0.01 & ceOutput_COAD$regSim !=0,]
ceOutput2_COAD2<-ceOutput_COAD2[ceOutput_COAD2$hyperPValue<0.01 & ceOutput_COAD2$corPValue<0.01 & ceOutput_COAD2$regSim !=0,]

edges_COAD<-gdcExportNetwork(ceNetwork = ceOutput2_COAD,net = "edges")
nodes_COAD<-gdcExportNetwork(ceNetwork = ceOutput2_COAD,net = "nodes")

write.table(edges_COAD,file = "edges_COAD.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_COAD,file="nodes_COAD.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(ceOutput2_COAD,file="ceOutput_COAD.txt",sep="\t",quote = FALSE,row.names=FALSE)
head(edges_COAD)




#########I have sent this predictions to my third aim analysis. The ceOutput_COAD2 contains miRNA-target predictions for 
##First i will extract gene-miRNA predictions
miRNA_COAD<-strsplit(ceOutput_COAD2$miRNAs,split = ",")
mRNA_miRNA_COAD<-data.frame(V1=rep(ceOutput_COAD2$Genes,sapply(miRNA_COAD,length)),V2=unlist(miRNA_COAD))
lncRNA_miRNA_COAD<-data.frame(V1=rep(ceOutput_COAD2$lncRNAs,sapply(miRNA_COAD,length)),V2=unlist(miRNA_COAD))

##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
mRNA_COAD<-mRNA_miRNA_COAD$V1
symbol_mRNA_COAD<-ensembldb::select(EnsDb.Hsapiens.v79,keys =mRNA_COAD,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(mRNA_miRNA_COAD,symbol_mRNA_COAD,by=c("V1"="GENEID"))
mRNA_miRNA_COAD<-data.frame(d1$V2,d1$SYMBOL)
names(mRNA_miRNA_COAD)<-c("miRNA","target")
View(mRNA_miRNA_COAD)

###Then, I will extract lncRNA-miRNA predictions
##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
lncRNA_COAD<-lncRNA_miRNA_COAD$V1
symbol_lncRNA_COAD<-ensembldb::select(EnsDb.Hsapiens.v79,keys =lncRNA_COAD,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(lncRNA_miRNA_COAD,symbol_lncRNA_COAD,by=c("V1"="GENEID"))
lncRNA_miRNA_COAD<-data.frame(d1$V2,d1$SYMBOL)
names(lncRNA_miRNA_COAD)<-c("miRNA","target")
lncRNA_miRNA_COAD$target<-paste("lncRNA",lncRNA_miRNA_COAD$target,sep = "-")
View(lncRNA_miRNA_COAD)
#Binding mRNA and lncRNA predictions
mRNA_lncRNA_miRNA_COAD<-rbind(mRNA_miRNA_COAD,lncRNA_miRNA_COAD)
miRNA.target.interactions_COAD<-mRNA_lncRNA_miRNA_COAD
write.table(miRNA.target.interactions,file = "Z:/GENETICS/DNA_Methylation/COLCA/miRNA.target.interactions_COAD.txt",sep = "\t",quote=FALSE)

