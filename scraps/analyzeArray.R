library(beadarray)
library(affy)

targets = read.table("../mirnaSamples_cebpa.txt",sep="\t",header=TRUE,as.is=TRUE)
targets = read.table("../mirnaSamples_hnf1.txt",sep="\t",header=TRUE,as.is=TRUE)
targets = read.table("../mirnaSamples_hnf4af.txt",sep="\t",header=TRUE,as.is=TRUE)
targets = read.table("../mRNAsamples_hnf4a.txt",sep="\t",header=TRUE,as.is=TRUE)

BLData = readIllumina(arrayNames = targets$ArrayName, textType=".txt", targets = targets , backgroundMethod = "none")
BLData = readIllumina(arrayNames = targets$ArrayName, textType=".txt", targets = targets , useImages=FALSE)
par(mfrow=c(4,3),mai=c(0.2,0.2,0.2,0.2))
boxplotBeads(BLData,las=2,ylim=c(4,12),main="Foreground",outline=FALSE)
boxplotBeads(BLData,las=2,ylim=c(4,12),whatToPlot="Gb",main="Background",outline=FALSE)
BLData.bc = backgroundCorrect(BLData,method="subtract")
boxplotBeads(BLData.bc,las=2,ylim=c(4,12),main="Foreground Corrected")

bashOut = BASH(BLData.bc,array=1:36)
BLData.bash = setWeights(BLData.bc,bashOut$wts,array=1:36)

# outliers
outliers = NULL
for (i in 1:10) { outliers[i] = length(findAllOutliers(BLData,array=i))}


# look at images
par(mfrow=c(2,1))
an = arrayNames(BLData.bc)
for (i in 1:9) { imageplot(BLData.bc,array=i,nrow=20,ncol=200,zlim=range(5,8),what="G",main=an[i])}
for (i in 1:18) { imageplot(BLData.bc,array=i,nrow=20,ncol=200,zlim=range(5,8),what="G",main=an[i])}

BSData = createBeadSummaryData(BLData.bc)
BSData = createBeadSummaryData(BLData.bc,imagesPerArray=2)
par(mfrow=c(1,2))
boxplot(as.data.frame(log2(exprs(BSData))),las=2,ylab="log2(intensity)")
boxplot(as.data.frame(NoBeads(BSData)),las=2,ylab="number of beads")

BSData.q = normaliseIllumina(BSData,method="quantile",transform="log2")

par(mfrow=c(2,1))
plotDensity(log2(exprs(BSData)), main="Non-normalized")
plotDensity(exprs(BSData.q), main="Normalized")

samples = as.factor(pData(BSData.q)$Origin)
d=dist(t(exprs(BSData.q)))
plot(hclust(d),labels=samples)

# or...
d <- dist(t(exprs(BSData.q)),method="euclidian")
plclust(hclust(d,method="complete"),main="Normalized Euclidian Clustering",labels=an)
plclust(hclust(d,method="complete"),main="Normalized Euclidian Clustering",labels=samples)

# or ... cluster on most variable
evar <- apply(exprs(BSData.q),1,var)
evarind <- evar > summary(evar)[5]
d <- dist(t(exprs(BSData.q)[evarind,]),method="euclidian")
plclust(hclust(d,method="complete"),main="Normalized Euclidian Clustering -- Most Variable",labels=samples)

boxplot(as.data.frame(exprs(BSData.q)),las=2)
design = model.matrix(~0+samples)
colnames(design) = levels(samples)
fit = lmFit(exprs(BSData.q),design)
cont.matrix = makeContrasts(cebpa = cEBPa_KO - cEBPa_WT,levels = design)
fit2 = contrasts.fit(fit,cont.matrix)
ebFit = eBayes(fit2)
topTable(ebFit,coef=1,number=20)

logWTmean = as.matrix(rowMeans(exprs(BSData.q)[,1:3]))
logKOmean = as.matrix(rowMeans(exprs(BSData.q)[,4:9]))
WTmean = 2^logWTmean
KOmean = 2^logKOmean


mouseMI = read.table("../mouseMIannot.txt", sep="\t", quote="", header=TRUE,as.is=TRUE,comment.char="")
mouseWG = read.table("../mouseWGannot.txt", sep="\t", quote="", header=TRUE,as.is=TRUE,comment.char="")

rownames(mouseMI) = mouseMI$Array_Address_Id
rownames(mouseWG) = mouseWG$Array_Address_Id

probeIds = rownames(exprs(BSData.q))
annotMap = match(probeIds,mouseMI$Array_Address_Id)
annotMap = match(probeIds,mouseWG$Array_Address_Id)

annot = cbind(arrayAddress = probeIds,
              Gene = as.character(mouseMI[annotMap,"ILMN_Gene"]),
              Mature = as.character(mouseMI[annotMap,"Mature_miRNA_seq"]),
              Chromosome = as.character(mouseMI[annotMap,"Chromosome"]),
              Coordinates = as.character(mouseMI[annotMap,"Probe_Coordinates"]),
              Strand = as.character(mouseMI[annotMap,"Probe_Chr_Orientation"]),
              logKOmean = as.character(logKOmean),
              logWTmean = as.character(logWTmean),
              KOmean = as.character(KOmean),
              WTmean = as.character(WTmean))
annot = cbind(arrayAddress = probeIds,
              Gene = as.character(mouseWG[annotMap,"Symbol"]),
              Accession = as.character(mouseWG[annotMap,"Accession"]),
              Refseq = as.character(mouseWG[annotMap,"RefSeq_Id"]),
              Unigene = as.character(mouseWG[annotMap,"Unigene_Id"]),
              logKOmean = as.character(logKOmean),
              logWTmean = as.character(logWTmean),
              KOmean = as.character(KOmean),
              WTmean = as.character(WTmean),
              Definition = as.character(mouseWG[annotMap,"Definition"]),
              Chromosome = as.character(mouseWG[annotMap,"Chromosome"]),
              Strand = as.character(mouseWG[annotMap,"Probe_Chr_Orientation"]),
              Coordinates = as.character(mouseWG[annotMap,"Probe_Coordinates"]),
              GO_Component = as.character(mouseWG[annotMap,"Ontology_Component"]),
              GO_Process = as.character(mouseWG[annotMap,"Ontology_Process"]),
              GO_Function = as.character(mouseWG[annotMap,"Ontology_Function"]))

ebFit$genes = annot
              
cebpaRaw = topTable(ebFit,coef=1,number=2000)
cebpa = cebpaRaw[!isNA(cebpaRaw$Gene),]
cebpa$FC = 2^cebpa$logFC

#microRNA
reorder = c(1,2,17,15,14,9,10,3,4,5,6,11,7,8,13,16,12)
#mRNA
reorder = c(1,2,21,19,6,7,3,8,18,9,4,5,9,10,11,14,13,12,15,16,17,20)

cebpaR = cebpa[,reorder]
write.table(cebpaR,file="../miRNAcebpa_results.txt")
