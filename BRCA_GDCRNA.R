project_BRCA<-'TCGA-BRCA'
rnadir_BRCA<-paste(project_BRCA,'RNAseq',sep="/")
mirdir_BRCA<-paste(project_BRCA,"miRNAs",sep="/")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-BRCA',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir_BRCA)
gdcRNADownload(project.id = "TCGA-BRCA",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir_BRCA)
#Download clinical data
clinicaldir_BRCA<-paste(project_BRCA,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-BRCA",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir_BRCA)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_BRCA<-gdcParseMetadata(project.id = "TCGA-BRCA",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_BRCA<-gdcFilterDuplicate(metaMatrix.RNA_BRCA)
metaMatrix.RNA_BRCA<-gdcFilterSampleType(metaMatrix.RNA_BRCA)

####Parse miRNA metadata
metaMatrix.MIR_BRCA<-gdcParseMetadata(project.id = "TCGA-BRCA",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_BRCA<-gdcFilterDuplicate(metaMatrix.MIR_BRCA)
metaMatrix.MIR_BRCA<-gdcFilterSampleType(metaMatrix.MIR_BRCA)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_BRCA<-gdcRNAMerge(metadata = metaMatrix.RNA_BRCA,path = rnadir_BRCA,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_BRCA<-gdcRNAMerge(metadata = metaMatrix.MIR_BRCA,path = mirdir_BRCA,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_BRCA<-gdcClinicalMerge(path = clinicaldir_BRCA,key.info = TRUE)

ClinicalDa_BRCA[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_BRCA<-gdcVoomNormalization(counts = rnaCounts_BRCA,filter = FALSE)
####Normalization of miRNAs data
mirExpr_BRCA<-gdcVoomNormalization(counts = mirCounts_BRCA,filter = FALSE)
####Differential gene expression analysis
DEGAll_BRCA<-gdcDEAnalysis(counts = rnaCounts_BRCA,group = metaMatrix.RNA_BRCA$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
DEGMIR_BRCA<-gdcDEAnalysis(counts = mirCounts_BRCA,group = metaMatrix.MIR_BRCA$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma")
####As I need log2FC for all set of genes I added a modification for gdcDEAnlaysis
DEGMIR_BRCA2<-gdcDEAnalysis(counts = mirCounts_BRCA,group = metaMatrix.MIR_BRCA$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
setDT(DEGMIR_BRCA2, keep.rownames = TRUE)[]
demiR_BRCA2<-gdcDEReport(DEGMIR_BRCA2,fc=0.05)
data("DEGAll_BRCA")
####All DEGs
deALL_BRCA<-gdcDEReport(deg = DEGAll_BRCA,gene.type = "all")
####DE long-noncoding
deLNC_BRCA<-gdcDEReport(deg = DEGAll_BRCA,gene.type = "long_non_coding")
####DE protein coding genes
dePC_BRCA<-gdcDEReport(deg=DEGAll_BRCA,gene.type = "protein_coding")
##DE pseudogenes
dePseudo_BRCA<-gdcDEReport(deg=DEGAll_BRCA,gene.type = "pseudogene")
demiR_BRCA<-gdcDEReport(deg = DEGMIR_BRCA)
################################################################################################
DEGALL_BRCA_MIR<-gdcDEAnalysis(counts = mirCounts_BRCA,group = metaMatrix.MIR_BRCA$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
gdcVolcanoPlot(DEGALL_BRCA_MIR)
gdcHeatmap(deg.id = rownames(DEGALL_BRCA_MIR),metadata = metaMatrix.MIR_BRCA,rna.expr = mirExpr_BRCA)

####DEG visualisation
gdcVolcanoPlot(DEGMIR_BRCA)
demiR_BRCA<-gdcDEReport(deg=DEGMIR_BRCA,gene.type = "all")
gdcBarPlot(deg=demiR_BRCA,angle = 45,data.type = "miRNAs")
gdcBarPlot(deg = deALL_BRCA,angle = 45,data.type = "RNAseq")
gdcHeatmap(deg.id = rownames(demiR_BRCA),metadata = metaMatrix.MIR_BRCA,rna.expr = mirExpr_BRCA)
degName=rownames(deALL_BRCA)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_BRCA,rna.expr = rnaExpr_BRCA)
gdcBarPlot(deg=demiR_BRCA,angle=45,data.type="miRNAs")
t1<-table(which(demiR_BRCA$logFC<0))
dim(t1)
t2<-table(which(demiR_BRCA$logFC>0))
dim(t2)
##ceRNAs network analysis

ceOutput_BRCA<- gdcCEAnalysis(lnc= rownames(deLNC_BRCA),pc= rownames(dePC_BRCA),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_BRCA,mir.expr= mirExpr_BRCA)
ceOutput_BRCA2<- gdcCEAnalysis(lnc= rownames(deLNC_BRCA),pc= rownames(dePC_BRCA),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_BRCA,mir.expr= mirExpr_BRCA)
MC_HC_ceOutput_BRCA<-ceOutput_BRCA2[ceOutput_BRCA2$hyperPValue<0.01 & ceOutput_BRCA2$cor>0.4,]
ceOutput_BRCA3<- gdcCEAnalysis(lnc= rownames(deLNC_BRCA),pc= rownames(dePC_BRCA),lnc.targets = 'miRcode',pc.targets= 'starBase',rna.expr= rnaExpr_BRCA,mir.expr    = mirExpr_BRCA)
ceOutput2_BRCA3<-ceOutput_BRCA3[ceOutput_BRCA3$hyperPValue<0.01 & ceOutput_BRCA3$cor>0.4,]


#########I have sent this predictions to my third aim analysis. The ceOutput_BRCA2 contains miRNA-target predictions for 
##First i will extract gene-miRNA predictions
miRNA_BRCA<-strsplit(ceOutput_BRCA2$miRNAs,split = ",")
mRNA_miRNA_BRCA<-data.frame(V1=rep(ceOutput_BRCA2$Genes,sapply(miRNA_BRCA,length)),V2=unlist(miRNA_BRCA))
lncRNA_miRNA_BRCA<-data.frame(V1=rep(ceOutput_BRCA2$lncRNAs,sapply(miRNA_BRCA,length)),V2=unlist(miRNA_BRCA))

##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
mRNA_BRCA<-mRNA_miRNA_BRCA$V1
symbol_mRNA_BRCA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =mRNA_BRCA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(mRNA_miRNA_BRCA,symbol_mRNA_BRCA,by=c("V1"="GENEID"))
mRNA_miRNA_BRCA<-data.frame(d1$V2,d1$SYMBOL)
names(mRNA_miRNA_BRCA)<-c("miRNA","target")
View(mRNA_miRNA_BRCA)

###Then, I will extract lncRNA-miRNA predictions
##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
lncRNA_BRCA<-lncRNA_miRNA_BRCA$V1
symbol_lncRNA_BRCA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =lncRNA_BRCA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(lncRNA_miRNA_BRCA,symbol_lncRNA_BRCA,by=c("V1"="GENEID"))
lncRNA_miRNA_BRCA<-data.frame(d1$V2,d1$SYMBOL)
names(lncRNA_miRNA_BRCA)<-c("miRNA","target")
lncRNA_miRNA_BRCA$target<-paste("lncRNA",lncRNA_miRNA_BRCA$target,sep = "-")
View(lncRNA_miRNA_BRCA)
#Binding mRNA and lncRNA predictions
mRNA_lncRNA_miRNA_BRCA<-rbind(mRNA_miRNA_BRCA,lncRNA_miRNA_BRCA)
miRNA.target.interactions<-mRNA_lncRNA_miRNA_BRCA
write.table(miRNA.target.interactions,file = "Z:/GENETICS/DNA_Methylation/BRCA/miRNA.target.interactions.txt",sep = "\t",quote=FALSE)
####Network Visualization in Cytoscape

ceOutput2_BRCA<-ceOutput_BRCA[ceOutput_BRCA$hyperPValue<0.01 & ceOutput_BRCA$corPValue<0.01 & ceOutput_BRCA$regSim !=0,]
ceOutput2_BRCA2<-ceOutput_BRCA2[ceOutput_BRCA2$hyperPValue<0.01 & ceOutput_BRCA2$corPValue<0.01 & ceOutput_BRCA2$regSim !=0,]
edges_BRCA<-gdcExportNetwork(ceNetwork = ceOutput2_BRCA,net = "edges")
nodes_BRCA<-gdcExportNetwork(ceNetwork = ceOutput2_BRCA,net = "nodes")

write.table(edges_BRCA,file = "edges_BRCA.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_BRCA,file="nodes_BRCA.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(ceOutput2_BRCA,file="ceOutput_BRCA.txt",sep="\t",quote = FALSE,row.names=FALSE)
head(edges)


#########I have sent this predictions to my third aim analysis. The ceOutput_BRCA2 contains miRNA-target predictions for 
##First i will extract gene-miRNA predictions
miRNA_BRCA<-strsplit(ceOutput_BRCA2$miRNAs,split = ",")
mRNA_miRNA_BRCA<-data.frame(V1=rep(ceOutput_BRCA2$Genes,sapply(miRNA_BRCA,length)),V2=unlist(miRNA_BRCA))
lncRNA_miRNA_BRCA<-data.frame(V1=rep(ceOutput_BRCA2$lncRNAs,sapply(miRNA_BRCA,length)),V2=unlist(miRNA_BRCA))

##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
mRNA_BRCA<-mRNA_miRNA_BRCA$V1
symbol_mRNA_BRCA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =mRNA_BRCA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(mRNA_miRNA_BRCA,symbol_mRNA_BRCA,by=c("V1"="GENEID"))
mRNA_miRNA_BRCA<-data.frame(d1$V2,d1$SYMBOL)
names(mRNA_miRNA_BRCA)<-c("miRNA","target")
View(mRNA_miRNA_BRCA)

###Then, I will extract lncRNA-miRNA predictions
##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
lncRNA_BRCA<-lncRNA_miRNA_BRCA$V1
symbol_lncRNA_BRCA<-ensembldb::select(EnsDb.Hsapiens.v79,keys =lncRNA_BRCA,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(lncRNA_miRNA_BRCA,symbol_lncRNA_BRCA,by=c("V1"="GENEID"))
lncRNA_miRNA_BRCA<-data.frame(d1$V2,d1$SYMBOL)
names(lncRNA_miRNA_BRCA)<-c("miRNA","target")
lncRNA_miRNA_BRCA$target<-paste("lncRNA",lncRNA_miRNA_BRCA$target,sep = "-")
View(lncRNA_miRNA_BRCA)
#Binding mRNA and lncRNA predictions
mRNA_lncRNA_miRNA_BRCA<-rbind(mRNA_miRNA_BRCA,lncRNA_miRNA_BRCA)
miRNA.target.interactions_BRCA<-mRNA_lncRNA_miRNA_BRCA
write.table(miRNA.target.interactions,file = "Z:/GENETICS/DNA_Methylation/BRCA/miRNA.target.interactions_BRCA.txt",sep = "\t",quote=FALSE)









