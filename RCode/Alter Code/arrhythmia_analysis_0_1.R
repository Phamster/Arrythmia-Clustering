require("psych")
require("AER")
require("mlbench")
require("cluster")
require("zoo")
library(ggplot2)
library(gridExtra)

# install.packages("missForest") für gemischt Datensätze
##  Daten einlesen ##
##setwd("C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt")
arrhythmia <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/arrhythmia/arrhythmia.data", header = FALSE, na.strings = "?" )
rawdata <- arrhythmia

# allgemeine Strukturanalyse
names(arrhythmia)
summary(arrhythmia)
describe(arrhythmia)
str(arrhythmia)

##  Datenbereinigung  ##

sumfun <- function(x){sum(abs(x))} 
test <- sapply(arrhythmia,sumfun)
# Missing Value Analyse
table(is.na(arrhythmia))
table(complete.cases(arrhythmia))
count.na <- function(x){sum(is.na(x))} 

# V14 376 Missing Values => wird gedroppt
missing <- sapply(arrhythmia,count.na)
missing[missing>0]

# lösche V14 
arrhythmia$V14 <- NULL

# bearbeite Zeilen mit NA
arrhythmia <- na.omit(arrhythmia)
# missforest fuer gemischte Datentypen
#arrhythmia <- na.approx(arrhythmia)
describe(arrhythmia)

# entferne Spalten mit Standardabweichung = 0
arrhythmia.omit <- arrhythmia[,apply(arrhythmia, 2, sd)==0] 
arrhythmia <- arrhythmia[,apply(arrhythmia, 2, sd)>0]
describe(arrhythmia)

# Verteilung der Herzrythmusstörungsklassen 1 == "keine Störung"
# extrahiere Ergebnisspalte als vorbereitung zur Hauptkomponentenanalyse
trueclus <- as.numeric(arrhythmia$V280)
table(trueclus)
arrhythmia$V280 <- NULL

# Daten Standardisierung
scaled.arrhythmia <- scale(arrhythmia)
colMeans(scaled.arrhythmia)
apply(scaled.arrhythmia, 2, sd)

# arrythmia <- arrhythmia[,1:14]
scaled.arrhythmia

##  Hierachical Clusteranalyse ohne normalisierung und ohne Merkmalsreduktion  ##

# berechne Distanzmatrix dd
dd <- daisy(arrhythmia,metric="gower")
hier.ave <- hclust(dd, method="ward.D2")
arrhythmiaagnes <- agnes(dd)
bannerplot(arrhythmiaagnes)
pltree(arrhythmiaagnes)
dendextend::cor_cophenetic(arrhythmiaagnes,dd)


cluslist <- lapply(1:16, 
                   function(obj)
                     cutree(hier.ave, obj))

funclus <- function(x, k){
  clus <- list(cluster=cluslist[[k]])
}

#gap <- clusGap(arrhythmia, funclus, K.max=5)
plot(gap)

print(gap)

# keine sinnvolle Clusterung

df <- table(cluslist[[16]], trueclus)


##  Hierachical Clusteranalyse mit normalisierung und ohne Merkmalsreduktion  ##

# berechne Distanzmatrix dd
dd <- dist(scaled.arrhythmia)
hier.ave <- hclust(dd, method="single")
arrhythmiaagnes <- agnes(dd)
bannerplot(arrhythmiaagnes)
pltree(arrhythmiaagnes)
dendextend::cor_cophenetic(arrhythmiaagnes, dd)


cluslist <- lapply(1:16, 
                   function(obj)
                     cutree(hier.ave, obj))

funclus <- function(x, k){
  clus <- list(cluster=cluslist[[k]])
}

#gap <- clusGap(scaled.arrhythmia, funclus, K.max=3)
plot(gap)

# keine sinnvolle Clusterung
table(cluslist[[16]], trueclus)

##  kmean Clusteranalyse mit normalisierung und ohne Merkmalsreduktion  ##

kmliste <- lapply(2:16, function(obj) 
  kmeans(scaled.arrhythmia, obj, nstart=25))
str(kmliste[[2]])
names(kmliste) <- 2:10
## Scree plot nach Backhaus et al.
plot(2:16, sapply(kmliste, 
                  function(obj) obj$tot.withinss), type="b")
## kmeans gegen tatsaechliche Klassifizierung
table(kmliste[[2]]$cluster, trueclus)
table(kmliste[[15]]$cluster, trueclus)

##  Dimensionsreduktion mit PCA ##
PCarrhythmia <- prcomp(scaled.arrhythmia, scale.=TRUE)
EVarrhythmia <- eigen(cor(scaled.arrhythmia))
plot(PCarrhythmia)
VSS.scree(scaled.arrhythmia)
plot(PCarrhythmia)

biplot(PCarrhythmia, cex=0.6, xpd=NA)

HKWerte <- predict(PCarrhythmia) 
cor(scaled.arrhythmia, HKWerte[,1:2]) ## Ladungen (Begriff loadings wird uneinheitlich genutzt)
rowSums(cor(scaled.arrhythmia, HKWerte[,1:2])^2) ## Kommunalitäten der ersten beiden HKen

pca_50 <- principal(scaled.arrhythmia,50,rotate="varimax") ## 6 factors, Kommunalit?ten 1
heatmap((pca_50$loadings))

screeplot(PCarrhythmia, type = "l", npcs = 100)
cumprob <- cumsum(PCarrhythmia$sdev^2 / sum(PCarrhythmia$sdev^2))
plot(cumprob[0:100], xlab = "HK #", ylab = "Erklaerte Varianz",
     main = "kumulierte Varianz")

principal(scaled.arrhythmia,79,rotate="varimax") ## 4 factors, höhere Kommunalitäten
HKWerte <- predict(PCarrhythmia) 
HKWerte <- PCarrhythmia$x

##  kmean Clusteranalyse mit normalisierung   ##

kmliste <- lapply(2:16, function(obj) 
  kmeans(HKWerte[,1:50], obj, nstart=25))
str(kmliste[[2]])
names(kmliste) <- 2:10
## Scree plot nach Backhaus et al.
plot(2:16, sapply(kmliste, 
                  function(obj) obj$tot.withinss), type="b")
## kmeans gegen tatsaechliche Klassifizierung
table(kmliste[[2]]$cluster, trueclus)
table(kmliste[[15]]$cluster, trueclus)
