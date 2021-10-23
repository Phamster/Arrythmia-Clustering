require("psych")
require("AER")
require("mlbench")
require("cluster")
require("zoo")
require("dendextend")
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

arrhythmia.sick <- arrhythmia[(trueclus>1),]
sick.trueclus <- trueclus[(trueclus>1)]

##  Hierachical Clusteranalyse ohne normalisierung und ohne Merkmalsreduktion  ##
# berechne Distanzmatrix dd
par(mfrow = c(1,1))
dd <- daisy(arrhythmia,metric="gower")

hier.wardD2 <- agnes(dd^2,method="ward")
dendextend::cor_cophenetic(hier.wardD2,dd)

hier.wardD2 <- hclust(dd, method="ward.D2")
cor_cophenetic(hier.wardD2, dd)

bannerplot(hier.wardD2)
pltree(hier.wardD2)

hier.ave <- agnes(dd,method="average")
bannerplot(hier.ave)
pltree(hier.ave)
dendextend::cor_cophenetic(hier.ave,dd)

hier.complete <- agnes(dd,method="complete")
bannerplot(hier.complete)
pltree(hier.complete)
dendextend::cor_cophenetic(hier.complete,dd)


cluslist <- lapply(1:13, 
                   function(obj)
                     cutree(hier.ave, obj))

funclus <- function(x, k){
  clus <- list(cluster=cluslist[[k]])
}

#gap <- clusGap(arrhythmia, funclus, K.max=13)
plot(gap)

print(gap)

# keine sinnvolle Clusterung
table(cluslist[[9]], trueclus)


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
plot(PCarrhythmia)
VSS.scree(arrhythmia)
plot(PCarrhythmia)

biplot(PCarrhythmia, cex=0.6, xpd=NA)

HKWerte <- predict(PCarrhythmia) 
cor(arrhythmia, HKWerte[,1:2]) ## Ladungen (Begriff loadings wird uneinheitlich genutzt)
rowSums(cor(arrhythmia, HKWerte[,1:2])^2) ## Kommunalitäten der ersten beiden HKen

pca_50 <- principal(arrhythmia,34,rotate="varimax") ## 6 factors, Kommunalit?ten 1
heatmap((pca_50$loadings))

screeplot(PCarrhythmia, type = "l", npcs = 112)
cumprob <- cumsum(PCarrhythmia$sdev^2 / sum(PCarrhythmia$sdev^2))
plot(cumprob[0:112], xlab = "HK #", ylab = "Erklaerte Varianz",
     main = "kumulierte Varianz")

principal(arrhythmia,112,rotate="varimax") ## 4 factors, höhere Kommunalitäten
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