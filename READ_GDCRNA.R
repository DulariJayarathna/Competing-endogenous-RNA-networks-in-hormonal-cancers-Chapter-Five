project<-'TCGA-READ'
rnadir<-paste(project,"RNAseq",sep="/")
mirdir<-paste(project,"miRNAs",sep="/")
library(GDCRNATools)
##Download RNA and mature miRNA expression data
gdcRNADownload(project.id = 'TCGA-READ',data.type = "RNAseq",write.manifest = FALSE,method = 'gdc-client',directory = rnadir)
gdcRNADownload(project.id = "TCGA-READ",data.type = "miRNAs",write.manifest = FALSE,method = "gdc-client",directory = mirdir)
#Download clinical data
clinicaldir<-paste(project,"clinical",sep="/")
gdcClinicalDownload(project.id = "TCGA-READ",write.manifest = FALSE,method = 'gdc-client',directory = clinicaldir)

##Data Organization and DE analysis
###Parse metadata
####Parse RNAseq metadata
metaMatrix.RNA_READ<-gdcParseMetadata(project.id = "TCGA-READ",data.type = "RNAseq",write.meta = FALSE)
####Filter duplicated samples in RNAseq metadata
metaMatrix.RNA_READ<-gdcFilterDuplicate(metaMatrix.RNA_READ)
metaMatrix.RNA_READ<-gdcFilterSampleType(metaMatrix.RNA_READ)

####Parse miRNA metadata
metaMatrix.MIR_READ<-gdcParseMetadata(project.id = "TCGA-READ",data.type = "miRNAs",write.meta = FALSE)
####Filter duplicated samples in miRNA metadata
metaMatrix.MIR_READ<-gdcFilterDuplicate(metaMatrix.MIR_READ)
metaMatrix.MIR_READ<-gdcFilterSampleType(metaMatrix.MIR_READ)

###Merge raw counts data
####Merge RNAseq data
rnaCounts_READ<-gdcRNAMerge(metadata = metaMatrix.RNA_READ,path = rnadir,organized = FALSE,data.type = "RNAseq")
####Merge miRNA data
mirCounts_READ<-gdcRNAMerge(metadata = metaMatrix.MIR_READ,path = mirdir,organized = FALSE,data.type = "miRNAs")
####Merge clinical data
ClinicalDa_READ<-gdcClinicalMerge(path = clinicaldir,key.info = TRUE)
ClinicalDa_READ[1:6,5:10]

###TMM Normalization and voom transformation

####Normalization of RNAseq data
rnaExpr_READ<-gdcVoomNormalization(counts = rnaCounts_READ,filter = FALSE)
####Normalization of miRNAs data
mirExpr_READ<-gdcVoomNormalization(counts = mirCounts_READ,filter = FALSE)
####Differential gene expression analysis
DEGAll_READ<-gdcDEAnalysis(counts = rnaCounts_READ,group = metaMatrix.RNA_READ$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
DEGMIR_READ<-gdcDEAnalysis(counts = mirCounts_READ,group = metaMatrix.MIR_READ$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma")
####As I need log2FC for all set of genes I added a modification for gdcDEAnlaysis
DEGMIR_READ2<-gdcDEAnalysis(counts = mirCounts_READ,group = metaMatrix.MIR_READ$sample_type,comparison ="PrimaryTumor-SolidTissueNormal",method = "limma",filter = FALSE)
setDT(DEGMIR_READ2, keep.rownames = TRUE)[]

data("DEGAll")
####All DEGs
deALL_READ<-gdcDEReport(deg = DEGAll_READ,gene.type = "all")
####DE long-noncoding
deLNC_READ<-gdcDEReport(deg = DEGAll_READ,gene.type = "long_non_coding")
####DE protein coding genes
dePC_READ<-gdcDEReport(deg=DEGAll_READ,gene.type = "protein_coding")
##DE pseudogenes
dePseudo_READ<-gdcDEReport(deg=DEGAll_READ,gene.type = "pseudogene")
demiR_READ<-gdcDEReport(DEGMIR_READ)
demiR_READ2<-gdcDEReport(DEGMIR_READ2,fc=0.05)
################################################################################################
DEGALL_READ_MIR<-gdcDEAnalysis(counts = mirCounts_READ,group = metaMatrix.MIR_READ$sample_type,comparison = "PrimaryTumor-SolidTissueNormal",method = "limma")
gdcVolcanoPlot(DEGALL_READ_MIR)
gdcHeatmap(deg.id = rownames(DEGALL_READ_MIR),metadata = metaMatrix.MIR_READ,rna.expr = mirExpr_READ)
####DEG visualisation
gdcVolcanoPlot(DEGAll_READ)
gdcBarPlot(deg = deALL_READ,angle = 45,data.type = "RNAseq")
degName=rownames(deALL_READ)
gdcHeatmap(deg.id = degName,metadata = metaMatrix.RNA_READ,rna.expr = rnaExpr_READ)
gdcBarPlot(deg=demiR_READ,angle=45,data.type="miRNAs")
t1<-table(which(demiR_READ$logFC<0))
dim(t1)
t2<-table(which(demiR_READ$logFC>0))
dim(t2)
##ceRNAs network analysis

ceOutput_READ<- gdcCEAnalysis(lnc= rownames(deLNC_READ),pc= rownames(dePC_READ),lnc.targets = 'starBase',pc.targets= 'starBase',rna.expr= rnaExpr_READ,mir.expr    = mirExpr_READ)
ceOutput_READ2<- gdcCEAnalysis(lnc= rownames(deLNC_READ),pc= rownames(dePC_READ),lnc.targets = 'miRcode',pc.targets= 'miRcode',rna.expr= rnaExpr_READ,mir.expr    = mirExpr_READ)
MC_HC_ceOutput_READ<-ceOutput_READ2[ceOutput_READ2$hyperPValue<0.01 & ceOutput_READ2$cor>0.4,]
ceOutput_READ3<- gdcCEAnalysis(lnc= rownames(deLNC_READ),pc= rownames(dePC_READ),lnc.targets = 'miRcode',pc.targets= 'starBase',rna.expr= rnaExpr_READ,mir.expr    = mirExpr_READ)
ceOutput2_READ3<-ceOutput_READ3[ceOutput_READ3$hyperPValue<0.01 & ceOutput_READ3$cor>0.4,]
####Network Visualization in Cytoscape

ceOutput2_READ<-ceOutput_READ[ceOutput_READ$hyperPValue<0.01 & ceOutput_READ$corPValue<0.01 & ceOutput_READ$regSim !=0,]
ceOutput2_READ2<-ceOutput_READ2[ceOutput_READ2$hyperPValue<0.01 & ceOutput_READ2$corPValue<0.01 & ceOutput_READ2$regSim !=0,]
edges_READ<-gdcExportNetwork(ceNetwork = ceOutput2_READ,net = "edges")
nodes_READ<-gdcExportNetwork(ceNetwork = ceOutput2_READ,net = "nodes")

write.table(edges_READ,file = "edges_READ.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(nodes_READ,file="nodes_READ.txt",sep="\t",quote = FALSE,row.names = FALSE)
write.table(ceOutput2_READ,file="ceOutput_READ.txt",sep="\t",quote = FALSE,row.names=FALSE)
head(edges_READ)



#########I have sent this predictions to my third aim analysis. The ceOutput_READ2 contains miRNA-target predictions for 
##First i will extract gene-miRNA predictions
miRNA_READ<-strsplit(ceOutput_READ2$miRNAs,split = ",")
mRNA_miRNA_READ<-data.frame(V1=rep(ceOutput_READ2$Genes,sapply(miRNA_READ,length)),V2=unlist(miRNA_READ))
lncRNA_miRNA_READ<-data.frame(V1=rep(ceOutput_READ2$lncRNAs,sapply(miRNA_READ,length)),V2=unlist(miRNA_READ))

##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
mRNA_READ<-mRNA_miRNA_READ$V1
symbol_mRNA_READ<-ensembldb::select(EnsDb.Hsapiens.v79,keys =mRNA_READ,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(mRNA_miRNA_READ,symbol_mRNA_READ,by=c("V1"="GENEID"))
mRNA_miRNA_READ<-data.frame(d1$V2,d1$SYMBOL)
names(mRNA_miRNA_READ)<-c("miRNA","target")
View(mRNA_miRNA_READ)

###Then, I will extract lncRNA-miRNA predictions
##Convert ensembel ID to gene symbol
library(EnsDb.Hsapiens.v79)
lncRNA_READ<-lncRNA_miRNA_READ$V1
symbol_lncRNA_READ<-ensembldb::select(EnsDb.Hsapiens.v79,keys =lncRNA_READ,keytype = "GENEID",columns = c("SYMBOL","GENEID"))
d1<-dplyr::left_join(lncRNA_miRNA_READ,symbol_lncRNA_READ,by=c("V1"="GENEID"))
lncRNA_miRNA_READ<-data.frame(d1$V2,d1$SYMBOL)
names(lncRNA_miRNA_READ)<-c("miRNA","target")
lncRNA_miRNA_READ$target<-paste("lncRNA",lncRNA_miRNA_READ$target,sep = "-")
View(lncRNA_miRNA_READ)
#Binding mRNA and lncRNA predictions
mRNA_lncRNA_miRNA_READ<-rbind(mRNA_miRNA_READ,lncRNA_miRNA_READ)
miRNA.target.interactions_READ<-mRNA_lncRNA_miRNA_READ

###We joined predictions for COAD and read
temp3 <- rbind(miRNA.target.interactions_COAD, miRNA.target.interactions_READ)
library(data.table)
temp3 <- setkey(temp3, NULL)
temp3 <- unique(temp3)
head(temp3)
mirna.target.interactions_COLCA<-temp3
write.table(mirna.target.interactions_COLCA,file = "Z:/GENETICS/DNA_Methylation/COLCA/miRNA.target.interactions.txt",sep = "\t",quote=FALSE)





