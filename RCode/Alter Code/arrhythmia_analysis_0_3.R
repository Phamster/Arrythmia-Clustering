require("psych")
require("AER")
require("mlbench")
require("cluster")
require("zoo")
require("dendextend")
require("factoextra")

# install.packages("missForest") für gemischt Datensätze
##  Daten einlesen ##
setwd("C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt/Daten")
rawdata <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/arrhythmia/arrhythmia.data", header = FALSE, na.strings = "?" )
arrhythmia <- read.csv("scaled.arrhythmia.csv",header = TRUE)
unscaled.arrhythmia <- read.csv("unscaled.arrhythmia.csv",header = TRUE)
# Verteilung der Herzrythmusstörungsklassen 1 == "keine Störung"
# extrahiere Ergebnisspalte als vorbereitung zur Hauptkomponentenanalyse

trueclus <- as.numeric(arrhythmia$class)
table(trueclus)
arrhythmia$class <- NULL
raw.trueclus <- as.numeric(rawdata$V280)
rawdata$V280 <- NULL
unscaled.arrhythmia$class <- NULL

arrhythmia.sick <- arrhythmia[(trueclus>1)]
sick.trueclus <- trueclus[(trueclus>1)]

##  Hierachical Clusteranalyse ohne normalisierung und ohne Merkmalsreduktion  ##
# berechne Distanzmatrix dd
par(mfrow = c(1,1))
dd <- daisy(arrhythmia,metric="gower")
hier.ave <- hclust(dd, method="average")
cor_cophenetic(hier.ave, dd)

hier.complete <- hclust(dd, method="complete")
cor_cophenetic(hier.complete, dd)

hier.wardD2 <- hclust(dd, method="ward.D2")
cor_cophenetic(hier.wardD2, dd)

plot(hier.ave)
plot(hier.complete)
plot(hier.wardD2)


cluslist.wardD2 <- lapply(1:15, 
                   function(obj)
                     cutree(hier.wardD2, obj))

funclus <- function(x, k){
  clus <- list(cluster=cluslist.wardD2[[k]])
}

gap.wardD2 <- clusGap(arrhythmia, funclus, K.max=15)
plot(gap.wardD2)
print(gap.wardD2)

clusters=cutree(hier.wardD2, 12)
bannerplot(gap.wardD2)
sil.wardD2 <- silhouette(clusters, dd)
fviz_silhouette(sil.wardD2)

print(gap.wardD2, method="Tibs2001SEmax")

# keine sinnvolle Clusterung
table(cluslist.wardD2[[12]], trueclus)

table(cluslist.wardD2[[6]], trueclus)


cluslist.complete <- lapply(1:13, 
                   function(obj)
                     cutree(hier.complete, obj))

funclus <- function(x, k){
  clus <- list(cluster=cluslist.complete[[k]])
}

gap.complete <- clusGap(arrhythmia, funclus, K.max=13)
plot(gap.complete)
print(gap.complete)

clusters=cutree(gap.complete, 3)
bannerplot(gap.complete)
sil.complete <- silhouette(clusters, dd)
fviz_silhouette(sil.complete)

print(gap.complete, method="Tibs2001SEmax")

# keine sinnvolle Clusterung
table(cluslist.complete[[4]], trueclus)

table(cluslist.complete[[11]], trueclus)

## Teste Clusterung ohne gesunde Patienten ##

arrhythmia.sick <- arrhythmia[(trueclus>1),]
sick.trueclus <- trueclus[(trueclus>1)]

dd <- daisy(arrhythmia.sick,metric="gower")

hier.ave <- hclust(dd, method="average")
cor_cophenetic(hier.ave, dd)

hier.complete <- hclust(dd, method="complete")
cor_cophenetic(hier.complete, dd)

hier.wardD2 <- hclust(dd, method="ward.D2")
cor_cophenetic(hier.wardD2, dd)

plot(hier.ave)
plot(hier.complete)
plot(hier.wardD2)

##

cluslist <- lapply(1:13, 
                   function(obj)
                     cutree(hier.wardD2, obj))

funclus <- function(x, k){
  clus <- list(cluster=cluslist[[k]])
}

gap <- clusGap(arrhythmia.sick, funclus, K.max=13)
plot(gap)

print(gap)

# keine sinnvolle Clusterung
table(cluslist[[7]], sick.trueclus)


# keine sinnvolle Clusterung
table(cluslist[[3]], trueclus)

##  kmean Clusteranalyse mit normalisierung und ohne Merkmalsreduktion  ##
'
k-means assume the variance of the distribution of each attribute (variable) is spherical;

all variables have the same variance;

the prior probability for all k clusters are the same, i.e. each cluster has roughly equal number of observations;
If any one of these 3 assumptions is violated, then k-means will fail.'
kmliste <- lapply(2:16, function(obj) 
  kmeans(arrhythmia, obj, nstart=25))
str(kmliste[[2]])
names(kmliste) <- 2:10
## Scree plot nach Backhaus et al.
plot(2:16, sapply(kmliste, 
                  function(obj) obj$tot.withinss), type="b")
## kmeans gegen tatsaechliche Klassifizierung ##
table(kmliste[[2]]$cluster, trueclus)
table(kmliste[[15]]$cluster, trueclus)


##  Dimensionsreduktion mit PCA ##
PCarrhythmia <- prcomp(arrhythmia, scale.=TRUE)
EVarrhythmia <- eigen(cor(arrhythmia))
VSS.scree(arrhythmia)
plot(PCarrhythmia)
unrot <- principal(arrhythmia,29,rotate="none")
unrot

unrot27 <- principal(arrhythmia,27,rotate="none")
round(100*unrot27$communality, 2)


cumprob <- cumsum(PCarrhythmia$sdev^2 / sum(PCarrhythmia$sdev^2))
plot(cumprob[0:112], xlab = "HK #", ylab = "Erklaerte Varianz",
     main = "kumulierte Varianz")


rot9 <- principal(arrhythmia,9,rotate="varimax")

rot9

pcsummary <- summary(PCarrhythmia)
pcsummary$rotation[,1:2]
biplot(PCarrhythmia, cex=0.6, xpd=NA)

rot2 <- principal(arrhythmia, nfactors=2,rotate="varimax")
biplot(rot2, labels=rownames(arrhythmia), cex=0.6, xpd=NA, pch=1)

HKWerte <- predict(PCarrhythmia) 
cor(arrhythmia, HKWerte[,1:2]) ## Ladungen (Begriff loadings wird uneinheitlich genutzt)
rowSums(cor(arrhythmia, HKWerte[,1:2])^2) ## Kommunalitäten der ersten beiden HKen

pca_50 <- principal(arrhythmia,34,rotate="varimax") ## 6 factors, Kommunalit?ten 1
heatmap((pca_50$loadings))

screeplot(PCarrhythmia, type = "l", npcs = 112)

principal(arrhythmia,13,rotate="none") ## 4 factors, höhere Kommunalitäten
HKWerte <- predict(PCarrhythmia) 
HKWerte <- PCarrhythmia$x





dd <- dist(HKWerte)
hier.ave <- hclust(dd,method="average")
plot(hier.ave)
kmliste <- lapply(2:16, function(obj) 
  kmeans(HKWerte[,1:82], obj, nstart=25))
str(kmliste[[2]])
names(kmliste) <- 2:10
## Scree plot nach Backhaus et al.
plot(2:16, sapply(kmliste, 
                  function(obj) obj$tot.withinss), type="b")
## kmeans gegen tatsaechliche Klassifizierung
table(kmliste[[2]]$cluster, trueclus)
table(kmliste[[15]]$cluster, trueclus)