source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(affy)
library(limma)
library(statmod)
setwd("C:/Users/Justin/Desktop/MicroArray")
countsCoding = readTargets("140625_Klijn_counts_coding.txt")
dim(countsCoding)

fit = princomp(data,cor=TRUE)
summary(fit)
loadings(fit)
plot(fit, type = "lines")
fit$scores
biplot(fit)

#information about the gene ID's for each of the samples
geneToTranscript = read.table("140625_Klijn_geneToTranscript.txt")

source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment")

#Make a summaried experiment by hand
#rownames is column 1
dim(countsCoding)
#number of rows is number of geneID
#number of columns is number of cancer cell lines
as(countsCoding, "RangedSummarizedExperiment")
ExpRowRanges = 
ExpColumnData = 
se = SummarizedExperiment(assays = SimpleList(counts = countsCoding[,2:676]),
                         rowRanges = ExpRowRanges,
                         colData = ExpColumnData)

columnsData = read.csv("cell line data E-MTAB-2706.csv")
even_indicies = seq(2, nrow(columnsData), by = 2)
columnsData = columnsData[even_indicies,]
rownames(columnsData) = paste("Sample",1:nrow(columnsData))
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
countsCoding2 = countsCoding[,-1]
rownames(countsCoding2) = countsCoding[,1]

dds = DESeqDataSetFromMatrix(countData = countsCoding2,
                             colData = columnsData,
                             design = ~1)
head(dds)
head(countsCoding2)
head(assay(dds))
#Check fragments that uniquely align to the genes
boxplot(colSums(assay(dds)), ylab = "read count", main = "Boxplot of total read counts for 675 Human Cancer Cell Lines")

colData(dds)
#Can enter genomic ranges into the row ranges
rowRanges(dds)
#Need to normalize the variance, because variance increases with mean
# 2 ways to do this: vsn() and rlog()

#try rlog()
#takes a year to run
#rld = rlog(dds)
#head(assay(rld))

#try variance stabilizing transformation
vsd = varianceStabilizingTransformation(dds)
head(assay(vsd))
#this also takes a lifetime to run
VSDCounts = readTargets("140625_Klijn_VSD_coding.txt")
VSDCountsCoding2 = VSDCounts[,-1]
rownames(VSDCountsCoding2) = VSDCounts[,1]

vsd = DESeqDataSetFromMatrix(countData = VSDCountsCoding2,
                             colData = columnsData,
                             design = ~1)

#still really slow
# instead use getVarianceStabilizedData

#Sample Distances
sampleDistances = dist(t(VSDCountsCoding2))

#name of cell lines


library("gplots")
library("RColorBrewer")
DistMatrix = as.matrix(sampleDistances)
#rownames(DistMatrix) = columnsData$Characteristics.cell.line.
rownames(DistMatrix) = columnsData$Characteristics.organism.part.
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc = hclust(sampleDistances)
heatmap.2(DistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE, labRow = FALSE, RowSideColors = as.vector(poop$Colors), ColSideColors = as.vector(poop$Colors))

dat = as.data.frame(rownames(DistMatrix))
colnames(dat) = "cellline"
n <- nlevels(as.factor(rownames(DistMatrix)))
dat.col <- data.frame(cellline =unique(as.factor(rownames(DistMatrix))),
                      Colors =rainbow(n))  ## you can also use rainbow(n)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = dat.col$cellline, # category labels
       col = dat.col$Colors,  # color key
       lty= 1,             # line style
       lwd = 5,            # line width
       cex = 0.5
)

poop = merge(dat,dat.col)

#More Hierarchical Clustering
library(dendextend)
library(circlize)
hc2 = hc
hc2$labels = rownames(DistMatrix)
dend = as.dendrogram(hc2)
dend = dend %>%
  color_branches(k=4) %>%
  color_labels
par(mar = rep(0,4))
circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)

#Interesting to see the patterns, the squares are likely similar cell lines/types
#We can take a look based on the type of cancer and the location of the cancer
#Lazy way out is to simply create a distance matrix and use the heatmap to plot

bodyParts = columnsData$Factor.Value.organism.part.
bodyPartsMatrix = matrix(rep(0,(length(bodyParts)*length(bodyParts))), 
                         nrow = length(bodyParts), 
                         ncol = length(bodyParts))
for (i in 1:length(bodyParts)){
  bodyPartsMatrix[i,] = as.numeric(1 - 1 *(bodyParts[i] == bodyParts))
}

rownames(bodyPartsMatrix) = columnsData$Characteristics.cell.line.
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#hc2 = hclust(bodyPartsMatrix)
##herpderp need a disimilarity object
bodyPartsMatrixDistances = as.dist(bodyPartsMatrix)
hc2 = hclust(bodyPartsMatrixDistances)
bodyPartsMatrixDistances = as.matrix(bodyPartsMatrixDistances)

heatmap.2(bodyPartsMatrixDistances, Rowv=as.dendrogram(hc2),
          symm=TRUE, trace="none", col=colors,
          margins=c(2,10), labCol=FALSE)

cancerType = columnsData$Factor.Value.disease.

###This doesn't actually work like I thought it 
# did oops.





#PCA Plot
plotPCA(VSDCountsCoding2)


fit = princomp(VSDCountsCoding2,cor=TRUE)
summary(fit)
loadings(fit)
plot(fit, type = "lines")
fit$scores
biplot(fit)
#shoot i actually want the transpose
transposedCounts = t(VSDCountsCoding2)
dim(transposedCounts)
rownames(transposedCounts) = bodyParts

fit = prcomp(transposedCounts,cor=TRUE)
summary(fit)
loadings(fit)
plot(fit, type = "lines")
fit$scores
biplot(fit)


bodyPlaces = as.character(columnsData$Characteristics.tissue.supergroup.)
library(devtools)
#install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(fit, obs.scale = 1, var.scale = 1, 
              groups = bodyPlaces, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
#printing other principle components
g <- ggbiplot(fit, choices = c(1,3), obs.scale = 1, var.scale = 1, 
              groups = bodyPlaces, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g <- ggbiplot(fit, choices = c(2,3), obs.scale = 1, var.scale = 1, 
              groups = bodyPlaces, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g <- ggbiplot(fit, choices = c(3,4), obs.scale = 1, var.scale = 1, 
              groups = bodyPlaces, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

g <- ggbiplot(fit, choices = c(5,6), obs.scale = 1, var.scale = 1, 
              groups = bodyPlaces, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)




#Let's try knn
library(FNN)

bodyCategories = columnsData$Characteristics.tissue.supergroup.

knnSolution = knn(train = trainingSet, test = testingSet, cl = y, k = 3)
#head(knnSolution)
#table(bodyCategories[-trainingFlag],knnSolution)
#cbind(bodyCategories [-trainingFlag],knnSolution)
#It kinda seems that these things don't cluster together particularly well
#based on where the cancer originates from. While some of the knn results might make sense
# based on the fact that for example skin and lung cancers are both epithelial, some of the other
# classifications are really confusing. Might be interesting to perform 
# similar analysis on cell lines that are healthy.

CV2 = matrix(rep(0,17*17),nrow = 17, ncol = 17)
rownames(CV2)=levels(bodyCategories)
colnames(CV2)=levels(bodyCategories)
for(i in 1:100){
  trainingFlag = sample(1:nrow(transposedCounts),450,replace=FALSE)
  trainingSet = transposedCounts[trainingFlag,]
  testingSet = transposedCounts[-trainingFlag,]
  y = bodyCategories [trainingFlag]
  knnSolution = knn(train = trainingSet, test = testingSet, cl = y, k = 5)
  
  
  knnSolution3 = c(knnSolution, 1:17)
  test = c(bodyCategories[-trainingFlag] , 1:17)
  CV = table(knnSolution3, test)
  for(i in 1: 17){
    CV[i,i] = CV[i,i]-1
  }
  
  rownames(CV)=levels(bodyCategories)
  colnames(CV)=levels(bodyCategories)
  #CV2 = CV
  CV2 = CV2 + CV
}

CV2
t(CV2)#transpose is what we want here
a = sum(diag(CV2))
b = sum(CV2)
a/b







#Testing stuff
library("airway")
data("airway")
head(airway)
head(data)



#Text datasets analysis
setwd("C:/Users/Justin/Desktop/TextData")
source("https://bioconductor.org/biocLite.R")
readtext = function(filename){
  textFileString = readChar(paste(filename, ".txt", sep = "")
                            , file.info(paste(filename, ".txt", sep = ""))$size)
  textFileString = gsub("\r?\n|\r", " ", textFileString)
  textFileString = gsub("<.*?>", " ", textFileString)
  wordFrequenciesTable = table(unlist(strsplit(tolower(textFileString), "\\W")))
  wordFrequenciesTable.df = as.data.frame(wordFrequenciesTable)
  colnames(wordFrequenciesTable.df) = c('word', 'freq')
  return(as.data.frame(wordFrequenciesTable.df))
}

#Fitzgerald Stuff
setwd("C:/Users/Justin/Desktop/TextData/Fitzgerald")
f1 = readtext("winterdreams")
f2 = readtext("bernicebobsherhair")
f3 = readtext("headandshoulders")
f4 = readtext("theicepalace")
f5 = readtext("theoffshorepirate")
f6 = readtext("thecuriouscaseofbenjaminbutton")
f7 = readtext("thebabyparty")
f8 = readtext("thefreshestboy")
f9 = readtext("thebridalparty")
f10 = readtext("anewleaf")
f11 = readtext("babylonrevisited")
f12 = readtext("crazysunday")
f13 = readtext("MayDay")
f14 = readtext("TheJellyBean")
f15 = readtext("TheCamelsBack")
f16 = readtext("PorcelainAndPink")
f17 = readtext("TheLeesOfHappiness")
f18 = readtext("MrIcky")
f19 = readtext("JeminaTheMountainGirl")
f20 = readtext("TheCutGlassBowl")
setwd("C:/Users/Justin/Desktop/TextData")

#Chekhov Stuff
setwd("C:/Users/Justin/Desktop/TextData/Chekhov")
c1 = readtext("adoctorsvisit")
c2 = readtext("theladywiththedog")
c3 = readtext("alivingchattel")
c4 = readtext("aboutlove")
c5 = readtext("thegrasshopper")
c6 = readtext("whowastoblame")
c7 = readtext("ananonymousstory")
c8 = readtext("volodya")
c9 = readtext("theprincess")
c10 = readtext("thedeathofagovernmentclerk")
c11 = readtext("theninny")
c12 = readtext("thebetrothed")
c13 = readtext("AStoryWithoutATitle")
c14 = readtext("Peasants")
c15 = readtext("AnArtistsStory")
c16 = readtext("AtHome")
c17 = readtext("AtChristmasTime")
c18 = readtext("TheDarling")
c19 = readtext("ANervousBreakdown")
c20 = readtext("Sleepy")
c21 = readtext("Boys")
c22 = readtext("Uprooted")
c23 = readtext("Happiness")
c24 = readtext("BadWeather")
c25 = readtext("Darkness")
c26 = readtext("TheLotteryTicket")
c27 = readtext("AMystery")
c28 = readtext("Frost")
c29 = readtext("Champagne")
c30 = readtext("Drunk")
c31 = readtext("Enemies")
setwd("C:/Users/Justin/Desktop/TextData")

#Pop Music
setwd("C:/Users/Justin/Desktop/TextData/PopMusic")
p1 = readtext("work")
p2 = readtext("7years")
p3 = readtext("whereisthelove")
p4 = readtext("OneDance")
p5 = readtext("Californication")
p6 = readtext("StairwaytoHeaven")
p7 = readtext("BohemianRhapsody")
p8 = readtext("HotelCalifornia")
p9 = readtext("Chandelier")
p10 = readtext("Heroes")
p11 = readtext("ManInTheMirror")
p12 = readtext("MyHeartWillGoOn")
p13 = readtext("CircleOfLife")
p14 = readtext("AWholeNewWorld")
p15 = readtext("BeautyAndTheBeast")
p16 = readtext("OneDayMore")
p17 = readtext("Memory")
p18 = readtext("MoMoneyMoProblems")
p19 = readtext("IDreamedADream")
p20 = readtext("Changes")
p21 = readtext("DearMama")
p22 = readtext("KeepYaHeadUp")
p23 = readtext("HeyJude")
p24 = readtext("WithOrWithoutYou")
p25 = readtext("BeautifulDay")
p26 = readtext("PaintItBlack")
p27 = readtext("Desperado")
p28 = readtext("EyeOfTheTiger")
p29 = readtext("NovemberRain")
p30 = readtext("IWillAlwaysLoveYou")
p31 = readtext("Human")
p32 = readtext("21Guns")

setwd("C:/Users/Justin/Desktop/TextData")

#Seminal Papers
setwd("C:/Users/Justin/Desktop/TextData/SeminalPapers")
sp1 = readtext("ridgeregression")
sp2 = readtext("rulefit")
sp3 = readtext("iPSC")
sp4 = readtext("proteinmeasurement")
sp5 = readtext("TheProblemOfSocialCost")
sp6 = readtext("InactivationEColi")
sp7 = readtext("GenomeCluster")
sp8 = readtext("MicroarraysIonization")
sp9 = readtext("ProteinDyeLabeling")
sp10 = readtext("DNASequencing")
sp11 = readtext("RNAIsolation")
sp12 = readtext("WesternBlot")
sp13 = readtext("ElectronDensity")
sp14 = readtext("ThermoChemistry")
sp15 = readtext("BasicLocalAlignment")
sp16 = readtext("SHELX")
sp17 = readtext("BLAST")
sp18 = readtext("SocialMedia")
sp19 = readtext("CRISPR")
sp20 = readtext("GeneticAlgorithms")

setwd("C:/Users/Justin/Desktop/TextData")

#Hemingways
setwd("C:/Users/Justin/Desktop/TextData/Hemingway")
h1 = readtext("TheShortHappyLifeofFrancisMacomber")
h2 = readtext("TheCapitaloftheWorld")
h3 = readtext("TheSnowsofKilimanjaro")
h4 = readtext("OldManattheBridge")
h5 = readtext("UpinMichigan")
h6 = readtext("OntheQuaiatSmyrna")
h7 = readtext("IndianCamp")
h8 = readtext("TheDoctorandtheDoctorsWife")
h9 = readtext("TheEndofSomething")
h10 = readtext("The Three-Day Blow")
h11 = readtext("TheBattler")
h12 = readtext("AVeryShortStory")
h13 = readtext("SoldiersHome")
h14 = readtext("TheRevolutionist")
h15 = readtext("MrandMrsElliot")
h16 = readtext("CatintheRain")               
h17 = readtext("OutofSeason")
h18 = readtext("CrossCountrySnow")
h19 = readtext("MyOldMan")
h20 = readtext("BigTwoHeartedRiverPart1")
h21 = readtext("BigTwoHeartedRiverPart2")
h22 = readtext("TheUndefeated")
h23 = readtext("InAnotherCountry")
h24 = readtext("HillsLikeWhiteElephants")
h25 = readtext("TheKillers")
h26 = readtext("CheTiDiceLaPatria")
h27 = readtext("FiftyGrand")
h28 = readtext("ASimpleEnquiry")
h29 = readtext("TenIndians")
h30 = readtext("ACanaryforOne")
h31 = readtext("AnAlpineIdyll")
h32 = readtext("APursuitRace")
h33 = readtext("BanalStory")
h34 = readtext("NowILayMe")
h35 = readtext("acleanwelllightedplace")
setwd("C:/Users/Justin/Desktop/TextData")

#O'Connor
setwd("C:/Users/Justin/Desktop/TextData/OConnor")
oc1 = readtext("TheGeranium")
oc2 = readtext("TheBarber")
oc3 = readtext("Wildcat")
oc4 = readtext("TheCrop")
oc5 = readtext("TheTurkey")
oc6 = readtext("TheTrain")
oc7 = readtext("ThePeeler")
oc8 = readtext("TheHeartOfThePark")
oc9 = readtext("AStrokeOfGoodFortune")
oc10 = readtext("EnochAndTheGorilla")
oc11 = readtext("AGoodManIsHardToFind")
oc12 = readtext("ALateEncounterWithTheEnemy")
oc13 = readtext("TheLifeYouSaveMayBeYourOwn")
oc14 = readtext("TheRiver")
oc15 = readtext("ACircleInTheFire")
oc16 = readtext("TheDisplacedPerson")
oc17 = readtext("ATempleOfTheHolyGhost")
oc18 = readtext("TheArtificialNigger")
oc19 = readtext("GoodCountryPeople")
oc20 = readtext("YouCantBeAnyPoorerThanDeath")
oc21 = readtext("GreenLeaf")
oc22 = readtext("AViewOfTheWoods")
oc23 = readtext("TheEnduringChill")
oc24 = readtext("TheComfortsOfHome")
oc25 = readtext("EverythingThatRisesMustConverge")
oc26 = readtext("ThePartridgeFestival")
oc27 = readtext("TheLameShallEnterFirst")
oc28 = readtext("WhyDoTheHeathenRage")
oc29 = readtext("Revelation")
oc30 = readtext("ParkersBack")
oc31 = readtext("JudgementDay")
setwd("C:/Users/Justin/Desktop/TextData")

wordFrequenciesDataframe = f1
for(i in 2:20){
  wordFrequenciesDataframe = merge(wordFrequenciesDataframe, 
                                   eval(parse(text=(paste("f", i, sep="")))),
                                   by = "word", 
                                   suffixes = c(paste("f", 1:i, sep="")),
                                   all = TRUE)
}
for(i in 1:31){
  wordFrequenciesDataframe = merge(wordFrequenciesDataframe, 
                                   eval(parse(text=(paste("c", i, sep="")))),
                                   by = "word", 
                                   suffixes = c(paste("c", 1:i, sep="")),
                                   all = TRUE)
}

for(i in 1:32){
  wordFrequenciesDataframe = merge(wordFrequenciesDataframe, 
                                   eval(parse(text=(paste("p", i, sep="")))),
                                   by = "word", 
                                   suffixes = c(paste("p", 1:i, sep="")),
                                   all = TRUE)
}

for(i in 1:20){
  wordFrequenciesDataframe = merge(wordFrequenciesDataframe, 
                                   eval(parse(text=(paste("sp", i, sep="")))),
                                   by = "word", 
                                   suffixes = c(paste("sp", 1:i, sep="")),
                                   all = TRUE)
}

for(i in 1:35){
  wordFrequenciesDataframe = merge(wordFrequenciesDataframe, 
                                   eval(parse(text=(paste("h", i, sep="")))),
                                   by = "word", 
                                   suffixes = c(paste("h", 1:i, sep="")),
                                   all = TRUE)
}

for(i in 1:31){
  wordFrequenciesDataframe = merge(wordFrequenciesDataframe, 
                                   eval(parse(text=(paste("oc", i, sep="")))),
                                   by = "word", 
                                   suffixes = c(paste("oc", 1:i, sep="")),
                                   all = TRUE)
}

wordFrequenciesDataframe[is.na(wordFrequenciesDataframe)] = 0
colnames(wordFrequenciesDataframe) = c("word", paste("freq.f", 1:20, sep=""), 
                                       paste("freq.c", 1:31, sep=""),
                                       paste("freq.p", 1:32, sep=""), 
                                       paste("freq.sp", 1:20, sep=""),
                                       paste("freq.h", 1:35, sep=""),
                                       paste("freq.oc", 1:31, sep=""))
saveRDS(wordFrequenciesDataframe, file = "wordFreqDF.RDS")
#####
#EDA
textDataClasses = c(rep("Fitzgerald", 20),
                    rep("Chekhov", 31),
                    rep("Pop Music", 32),
                    rep("Science Papers", 20),
                    rep("Hemingway", 35),
                    rep("OConnor", 31))
p =colSums(wordFrequenciesDataframe[-1,-1])
mean(p[1:20])
sqrt(var(as.numeric(p[1:20])))
mean(p[21:51])
sqrt(var(p[21:51]))
mean(p[52:83])
sqrt(var(p[52:83]))
mean(p[84:103])
sqrt(var(p[84:103]))
mean(p[104:138])
sqrt(var(p[104:138]))
mean(p[139:169])
sqrt(var(p[139:169]))
                 
boxplot(p ~ textDataClasses, main = "Boxplot of Number of Words in Texts by Author", ylab = "Number of Words")                 
n = wordFrequenciesDataframe$word
wf.df = as.data.frame(t(wordFrequenciesDataframe[,-1]))
colnames(wf.df) = n
#wf.df$myfactor <- factor(row.names(wf.df))
#as.data.frame(colnames(wordFrequenciesDataframe))
#any(wf.df < 0)
columnDataUsed = as.data.frame(cbind(
  colnames(wordFrequenciesDataframe)[-1],
  c(rep("Fitzgerald", 20),
    rep("Chekhov", 31),
    rep("Pop Music", 32),
    rep("Science Paper", 20),
    rep("Hemingway", 35),
    rep("OConnor", 31))))
colnames(columnDataUsed) = c("Frequency", "Author")
library(DESeq2)
#library(affycoretools)
dds = DESeqDataSetFromMatrix(countData = t(wf.df),
                             colData = columnDataUsed,
                             design = ~1)


###
#Plotting PCA
diagdds = dds
diagdds = estimateSizeFactors(diagdds)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dds2 = DESeq(dds)
plotPCA(assay(diagdds), intgroup = c("Author"))
rld <- rlog(dds2, blind=FALSE)
DESeq2::plotPCA(dds2, intgroup = c("Author"))
dim(diagvst)
plotPCA(rld, intgroup = c("Author"))

numberOfSamples = ncol(wordFrequenciesDataframe)



#Let's look at hierarchical clustering
#Sample Distances
sampleDistances = dist(t(diagvst))
library("gplots")
library("RColorBrewer")
DistMatrix = as.matrix(sampleDistances)
rownames(DistMatrix) = unlist(strsplit(rownames(wf.df), split = '.', fixed=TRUE))[seq(2,(numberOfSamples-1)*2, by = 2)]
colnames(DistMatrix) = rownames(DistMatrix)
#colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc = hclust(sampleDistances)
heatmap.2(DistMatrix, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", #col=colors,
          margins=c(2,10))#, labCol=TRUE )

#More Hierarchical Clustering
library(dendextend)
library(circlize)
hc2 = hc
hc2$labels = colnames(DistMatrix)
dend = as.dendrogram(hc2)
dend = dend %>%
  color_branches(k=6) %>%
  color_labels
par(mar = rep(0,4))
circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
#source("https://bioconductor.org/biocLite.R")
#biocLite("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(diagdds,normalized=TRUE)),decreasing=TRUE)[1:20]
nt <- normTransform(diagdds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(diagdds))

rld <- rlog(dds, blind=FALSE)


pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
#pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)

plotPCA(rld, intgroup = c("Author"))

#library(ggplot2)
#data <- plotPCA(rld, intgroup=c("Author", "Frequency"), returnData=TRUE)
#percentVar <- round(1000 * attr(data, "percentVar"))
#ggplot(data, aes(PC1, PC2, color = Author))
#  geom_point(size=3) 
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) 
#ylab(paste0("PC2: ",percentVar[2],"% variance"))


Authors = columnDataUsed$Author
fit = prcomp(t(diagvst),center = TRUE, scale. = TRUE)
plot(fit, type = "lines")
#print(fit)
#summary(fit)
#loadings(fit)
#plot(fit, type = "lines")
#fit$scores
#biplot(fit)

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(fit, obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
#printing other principle components
g <- ggbiplot(fit, choices = c(1,3), obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g <- ggbiplot(fit, choices = c(2,3), obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g <- ggbiplot(fit, choices = c(3,4), obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

g <- ggbiplot(fit, choices = c(5,6), obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

library(rgl)
plot3d(fit$scores[,1:3])



unnormalizedTransposedTermFreqMatrix = as.matrix(wf.df)
unnormalizedTransposedTermFreqMatrix = unnormalizedTransposedTermFreqMatrix/colSums(unnormalizedTransposedTermFreqMatrix)
#Removing the Science Papers
Authors = columnDataUsed$Author[c(1:83, 104:169)]
diagvstNoSP = unnormalizedTransposedTermFreqMatrix[,c(1:83, 104:169)]
fit = prcomp(t(diagvstNoSP),center = TRUE, scale. = TRUE)

plot(fit, type = "lines")

g <- ggbiplot(fit, obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
#printing other principle components
g <- ggbiplot(fit, choices = c(1,3), obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g <- ggbiplot(fit, choices = c(2,3), obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g <- ggbiplot(fit, choices = c(2,4), obs.scale = 1, var.scale = 1, 
              groups = Authors, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)



#########
#KNN
#Let's try knn
library(FNN)
transposedTermFreqMatrix = t(diagvst)

CV2 = matrix(rep(0,36),nrow = 6, ncol = 6)
rownames(CV2)=c("f","c","p","sp","h","oc")
colnames(CV2)=c("f","c","p","sp","h","oc")
for(i in 1:100){
  trainingFlag = c(sample(1:20,15,replace=FALSE),
                   sample(21:51,15,replace=FALSE),
                   sample(52:83,15,replace=FALSE),
                   sample(84:103,15,replace=FALSE),
                   sample(104:138,15,replace=FALSE),
                   sample(139:169,15,replace=FALSE)
  )
  trainingSet = transposedTermFreqMatrix[trainingFlag,]
  testingSet = transposedTermFreqMatrix[-trainingFlag,]
  textDataClasses = c(rep("1", 20),
                      rep("2", 31),
                      rep("3", 32),
                      rep("4", 20),
                      rep("5", 35),
                      rep("6", 31))
  
  y = textDataClasses[trainingFlag]
  knnSolution = knn(train = trainingSet, test = testingSet, cl = y, k = 2)
  #head(knnSolution)
  #cbind(textDataClasses[-trainingFlag],knnSolution)
  knnSolution3 = c(knnSolution, 1,2,3,4,5,6)
  test = c(textDataClasses[-trainingFlag] , 1,2,3,4,5,6)
  CV = table(knnSolution3, test)
  CV[1,1] = CV[1,1]-1
  CV[2,2] = CV[2,2]-1
  CV[3,3] = CV[3,3]-1
  CV[4,4] = CV[4,4]-1
  CV[5,5] = CV[5,5]-1
  CV[6,6] = CV[6,6]-1
  rownames(CV)=c("f","c","p","sp","h","oc")
  colnames(CV)=c("f","c","p","sp","h","oc")
  #CV2 = CV
  CV2 = CV2 + CV
}

CV2
t(CV2)#transpose is what we want here
a = sum(diag(CV2))
b = sum(CV2)
a/b
#1
#0.5868354
#2
#0.6631646
#3
#0.343038
#4
#0.415443
#5
#
#6
#

#cbind(textDataClasses[-trainingFlag],knnSolution)
#table(textDataClasses[-trainingFlag],knnSolution)
#### without variance stabilization

unnormalizedTransposedTermFreqMatrix = as.matrix(wf.df)
unnormalizedTransposedTermFreqMatrix = unnormalizedTransposedTermFreqMatrix/colSums(unnormalizedTransposedTermFreqMatrix)


CV2 = matrix(rep(0,36),nrow = 6, ncol = 6)
rownames(CV2)=c("f","c","p","sp","h","oc")
colnames(CV2)=c("f","c","p","sp","h","oc")
for(i in 1:100){
trainingFlag = c(sample(1:20,15,replace=FALSE),
                 sample(21:51,15,replace=FALSE),
                 sample(52:83,15,replace=FALSE),
                 sample(84:103,15,replace=FALSE),
                 sample(104:138,15,replace=FALSE),
                 sample(139:169,15,replace=FALSE)
)
trainingSet = unnormalizedTransposedTermFreqMatrix[trainingFlag,]
testingSet = unnormalizedTransposedTermFreqMatrix[-trainingFlag,]
textDataClasses = c(rep("1", 20),
                    rep("2", 31),
                    rep("3", 32),
                    rep("4", 20),
                    rep("5", 35),
                    rep("6", 31))

y = textDataClasses[trainingFlag]
knnSolution = knn(train = trainingSet, test = testingSet, cl = y, k = 1)
knnSolution3 = c(knnSolution, 1,2,3,4,5,6)
test = c(textDataClasses[-trainingFlag] , 1,2,3,4,5,6)
CV = table(knnSolution3, test)
CV[1,1] = CV[1,1]-1
CV[2,2] = CV[2,2]-1
CV[3,3] = CV[3,3]-1
CV[4,4] = CV[4,4]-1
CV[5,5] = CV[5,5]-1
CV[6,6] = CV[6,6]-1
rownames(CV)=c("f","c","p","sp","h","oc")
colnames(CV)=c("f","c","p","sp","h","oc")
#CV2 = CV
CV2 = CV2 + CV
}

CV2
t(CV2)
a = sum(diag(CV2))
b = sum(CV2)
a/b
#1
#0.5301951
#2
#0.4736709
#3
#0.5075949
#4
#0.4870886
#5
#0.4970886
#6
#0.5067089


#########
#SVM
library(e1071)
##transposedTermFreqMatrix = t(diagvst)

unnormalizedTransposedTermFreqMatrix = as.matrix(wf.df)
unnormalizedTransposedTermFreqMatrix = unnormalizedTransposedTermFreqMatrix/colSums(unnormalizedTransposedTermFreqMatrix)

CV2 = matrix(rep(0,36),nrow = 6, ncol = 6)
rownames(CV2)=c("f","c","p","sp","h","oc")
colnames(CV2)=c("f","c","p","sp","h","oc")
for(i in 1:100){
  trainingFlag = c(sample(1:20,15,replace=FALSE),
                   sample(21:51,15,replace=FALSE),
                   sample(52:83,15,replace=FALSE),
                   sample(84:103,15,replace=FALSE),
                   sample(104:138,15,replace=FALSE),
                   sample(139:169,15,replace=FALSE)
  )
  #trainingSet = transposedTermFreqMatrix[trainingFlag,]
  #testingSet = transposedTermFreqMatrix[-trainingFlag,]
  trainingSet = unnormalizedTransposedTermFreqMatrix[trainingFlag,]
  testingSet = unnormalizedTransposedTermFreqMatrix[-trainingFlag,]
  textDataClasses = c(rep("1", 20),
                      rep("2", 31),
                      rep("3", 32),
                      rep("4", 20),
                      rep("5", 35),
                      rep("6", 31))
  
  y = factor(textDataClasses[trainingFlag])
  svmSolution = svm(trainingSet,y)
  #head(knnSolution)
  #cbind(textDataClasses[-trainingFlag],knnSolution)
  pred <- predict(svmSolution, testingSet)
  
  svmSolution3 = c(pred, 1,2,3,4,5,6)
  test = c(textDataClasses[-trainingFlag] , 1,2,3,4,5,6)
  CV = table(svmSolution3, test)
  CV[1,1] = CV[1,1]-1
  CV[2,2] = CV[2,2]-1
  CV[3,3] = CV[3,3]-1
  CV[4,4] = CV[4,4]-1
  CV[5,5] = CV[5,5]-1
  CV[6,6] = CV[6,6]-1
  rownames(CV)=c("f","c","p","sp","h","oc")
  colnames(CV)=c("f","c","p","sp","h","oc")
  #CV2 = CV
  CV2 = CV2 + CV
}

CV2
t(CV2)#transpose is what we want here
a = sum(diag(CV2))
b = sum(CV2)
a/b
#variance stabilized
#0.9008861


#############
# Annotating Texts

annotateText = function(filename){
  anno = annotateFile(paste(filename, ".txt", sep = ""))
  saveRDS(anno, file = paste(filename, ".RDS", sep = ""))
}

#Fitzgerald Stuff
setwd("C:/Users/Justin/Desktop/TextData/Fitzgerald")
f1 = annotateText("winterdreams")
f2 = annotateText("bernicebobsherhair")
f3 = annotateText("headandshoulders")
f4 = annotateText("theicepalace")
f5 = annotateText("theoffshorepirate")
f6 = annotateText("thecuriouscaseofbenjaminbutton")
f7 = annotateText("thebabyparty")
f8 = annotateText("thefreshestboy")
f9 = annotateText("thebridalparty")
f10 = annotateText("anewleaf")
f11 = annotateText("babylonrevisited")
f12 = annotateText("crazysunday")
f13 = annotateText("MayDay")
f14 = annotateText("TheJellyBean")
f15 = annotateText("TheCamelsBack")
f16 = annotateText("PorcelainAndPink")
f17 = annotateText("TheLeesOfHappiness")
f18 = annotateText("MrIcky")
f19 = annotateText("JeminaTheMountainGirl")
f20 = annotateText("TheCutGlassBowl")
setwd("C:/Users/Justin/Desktop/TextData")

#Chekhov Stuff
setwd("C:/Users/Justin/Desktop/TextData/Chekhov")
c1 = annotateText("adoctorsvisit")
c2 = annotateText("theladywiththedog")
c3 = annotateText("alivingchattel")
c4 = annotateText("aboutlove")
c5 = annotateText("thegrasshopper")
c6 = annotateText("whowastoblame")
c7 = annotateText("ananonymousstory")
c8 = annotateText("volodya")
c9 = annotateText("theprincess")
c10 = annotateText("thedeathofagovernmentclerk")
c11 = annotateText("theninny")
c12 = annotateText("thebetrothed")
c13 = annotateText("AStoryWithoutATitle")
c14 = annotateText("Peasants")
c15 = annotateText("AnArtistsStory")
c16 = annotateText("AtHome")
c17 = annotateText("AtChristmasTime")
c18 = annotateText("TheDarling")
c19 = annotateText("ANervousBreakdown")
c20 = annotateText("Sleepy")
c21 = annotateText("Boys")
c22 = annotateText("Uprooted")
c23 = annotateText("Happiness")
c24 = annotateText("BadWeather")
c25 = annotateText("Darkness")
c26 = annotateText("TheLotteryTicket")
c27 = annotateText("AMystery")
c28 = annotateText("Frost")
c29 = annotateText("Champagne")
c30 = annotateText("Drunk")
c31 = annotateText("Enemies")
setwd("C:/Users/Justin/Desktop/TextData")

#Pop Music
setwd("C:/Users/Justin/Desktop/TextData/PopMusic")
p1 = annotateText("work")
p2 = annotateText("7years")
p3 = annotateText("whereisthelove")
p4 = annotateText("OneDance")
p5 = annotateText("Californication")
p6 = annotateText("StairwaytoHeaven")
p7 = annotateText("BohemianRhapsody")
p8 = annotateText("HotelCalifornia")
p9 = annotateText("Chandelier")
p10 = annotateText("Heroes")
p11 = annotateText("ManInTheMirror")
p12 = annotateText("MyHeartWillGoOn")
p13 = annotateText("CircleOfLife")
p14 = annotateText("AWholeNewWorld")
p15 = annotateText("BeautyAndTheBeast")
p16 = annotateText("OneDayMore")
p17 = annotateText("Memory")
p18 = annotateText("MoMoneyMoProblems")
p19 = annotateText("IDreamedADream")
p20 = annotateText("Changes")
p21 = annotateText("DearMama")
p22 = annotateText("KeepYaHeadUp")
p23 = annotateText("HeyJude")
p24 = annotateText("WithOrWithoutYou")
p25 = annotateText("BeautifulDay")
p26 = annotateText("PaintItBlack")
p27 = annotateText("Desperado")
p28 = annotateText("EyeOfTheTiger")
p29 = annotateText("NovemberRain")
p30 = annotateText("IWillAlwaysLoveYou")
p31 = annotateText("Human")
p32 = annotateText("21Guns")

setwd("C:/Users/Justin/Desktop/TextData")

#Seminal Papers
setwd("C:/Users/Justin/Desktop/TextData/SeminalPapers")
sp1 = annotateText("ridgeregression")
sp2 = annotateText("rulefit")
sp3 = annotateText("iPSC")
sp4 = annotateText("proteinmeasurement")
sp5 = annotateText("TheProblemOfSocialCost")
sp6 = annotateText("InactivationEColi")
sp7 = annotateText("GenomeCluster")
sp8 = annotateText("MicroarraysIonization")
sp9 = annotateText("ProteinDyeLabeling")
sp10 = annotateText("DNASequencing")
sp11 = annotateText("RNAIsolation")
sp12 = annotateText("WesternBlot")
sp13 = annotateText("ElectronDensity")
sp14 = annotateText("ThermoChemistry")
sp15 = annotateText("BasicLocalAlignment")
sp16 = annotateText("SHELX")
sp17 = annotateText("BLAST")
sp18 = annotateText("SocialMedia")
sp19 = annotateText("CRISPR")
sp20 = annotateText("GeneticAlgorithms")

setwd("C:/Users/Justin/Desktop/TextData")

#Hemingways
setwd("C:/Users/Justin/Desktop/TextData/Hemingway")
h1 = annotateText("TheShortHappyLifeofFrancisMacomber")
h2 = annotateText("TheCapitaloftheWorld")
h3 = annotateText("TheSnowsofKilimanjaro")
h4 = annotateText("OldManattheBridge")
h5 = annotateText("UpinMichigan")
h6 = annotateText("OntheQuaiatSmyrna")
h7 = annotateText("IndianCamp")
h8 = annotateText("TheDoctorandtheDoctorsWife")
h9 = annotateText("TheEndofSomething")
h10 = annotateText("The Three-Day Blow")
h11 = annotateText("TheBattler")
h12 = annotateText("AVeryShortStory")
h13 = annotateText("SoldiersHome")
h14 = annotateText("TheRevolutionist")
h15 = annotateText("MrandMrsElliot")
h16 = annotateText("CatintheRain")               
h17 = annotateText("OutofSeason")
h18 = annotateText("CrossCountrySnow")
h19 = annotateText("MyOldMan")
h20 = annotateText("BigTwoHeartedRiverPart1")
h21 = annotateText("BigTwoHeartedRiverPart2")
h22 = annotateText("TheUndefeated")
h23 = annotateText("InAnotherCountry")
h24 = annotateText("HillsLikeWhiteElephants")
h25 = annotateText("TheKillers")
h26 = annotateText("CheTiDiceLaPatria")
h27 = annotateText("FiftyGrand")
h28 = annotateText("ASimpleEnquiry")
h29 = annotateText("TenIndians")
h30 = annotateText("ACanaryforOne")
h31 = annotateText("AnAlpineIdyll")
h32 = annotateText("APursuitRace")
h33 = annotateText("BanalStory")
h34 = annotateText("NowILayMe")
h35 = annotateText("acleanwelllightedplace")
setwd("C:/Users/Justin/Desktop/TextData")

#O'Connor
setwd("C:/Users/Justin/Desktop/TextData/OConnor")
oc1 = annotateText("TheGeranium")
oc2 = annotateText("TheBarber")
oc3 = annotateText("Wildcat")
oc4 = annotateText("TheCrop")
oc5 = annotateText("TheTurkey")
oc6 = annotateText("TheTrain")
oc7 = annotateText("ThePeeler")
oc8 = annotateText("TheHeartOfThePark")
oc9 = annotateText("AStrokeOfGoodFortune")
oc10 = annotateText("EnochAndTheGorilla")
oc11 = annotateText("AGoodManIsHardToFind")
oc12 = annotateText("ALateEncounterWithTheEnemy")
oc13 = annotateText("TheLifeYouSaveMayBeYourOwn")
oc14 = annotateText("TheRiver")
oc15 = annotateText("ACircleInTheFire")
oc16 = annotateText("TheDisplacedPerson")
oc17 = annotateText("ATempleOfTheHolyGhost")
oc18 = annotateText("TheArtificialNigger")
oc19 = annotateText("GoodCountryPeople")
oc20 = annotateText("YouCantBeAnyPoorerThanDeath")
oc21 = annotateText("GreenLeaf")
oc22 = annotateText("AViewOfTheWoods")
oc23 = annotateText("TheEnduringChill")
oc24 = annotateText("TheComfortsOfHome")
oc25 = annotateText("EverythingThatRisesMustConverge")
oc26 = annotateText("ThePartridgeFestival")
oc27 = annotateText("TheLameShallEnterFirst")
oc28 = annotateText("WhyDoTheHeathenRage")
oc29 = annotateText("Revelation")
oc30 = annotateText("ParkersBack")
oc31 = annotateText("JudgementDay")
setwd("C:/Users/Justin/Desktop/TextData")



