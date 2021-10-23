require("psych")
require("AER")
require("mlbench")
require("cluster")
require("zoo")
require("plyr")
require("df.table")
require("factoextra")
require("matrixStats") ## for rowSums, rowMeans....

require("outliers") ## for grubbs test

##  Daten einlesen ##
##setwd("C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt")
arrhythmia <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/arrhythmia/arrhythmia.data", header = FALSE, na.strings = "?" )
rawdata <- arrhythmia
arrhythmia <- rawdata ## dim(arrhythmia) [1] 452 280
## allgemeine Strukturanalyse
'names(arrhythmia)
summary(arrhythmia)
describe(arrhythmia)
str(arrhythmia)'
farben <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#B0E2FF","#EE799F","#CCEBC5")

## ergaenze relevante Spaltennamen
attnames <- c("Age","Sex","Height","Weight","QRS_duration","P_R_interval",
              "Q_T_interval","T_interval","P_interval","QRS","T","P","QRST","J","Heartrate","class")
names(arrhythmia)[c(1:15,280)] <- attnames 

psych::describe(arrhythmia)

## Untersuchung Klassifizierung

trueclus <- as.numeric(arrhythmia$class)
table(trueclus)

temp <- as.numeric(arrhythmia$class)
isSick <- c()
  for(i in 1:length(temp)){
    if(temp[i] > 1){
      isSick[i] <- 1
    }
    else{
      isSick[i] <- 0
    }
  }
arrhythmia$sick <- isSick

boxplot(arrhythmia$Age~arrhythmia$sick, names=c("Nein","Ja"), col=farben[1:2], 
        main="Krankheitsfaelle nach Alter",
        xlab="Herzrythmusstoerung erkannt", ylab="Alter")
barplot(table(arrhythmia$Sex,isSick),
        col=farben[4:5] ,main="Krankheitsfaelle nach Geschlecht",names=c("Nein","Ja"), 
        xlab="Herzrythmusstoerung erkannt", ylab="Haeufigkeit"); legend("topright",
         title="Geschlecht", c("Maenner","Frauen"),cex=0.8, fill = farben[4:5])
boxplot(rm~data$UNIVERSITY, col=farben[1:2],names=c("UOC","UPD"), xlab="UniversitÃ¤t", ylab="Durchschnitt",main="Durchschnitt je UniversitÃ¤t" )

##  Datenbereinigung  ##

## Missing Value Analyse
table(is.na(arrhythmia))
table(complete.cases(arrhythmia))
count.na <- function(x){sum(is.na(x))} 

## J 376 Missing Values
## lösche J Auspraegung nach der Zusammenfassung um Verrueckung zu vermeiden
## ersetze J durch Nullzeile vor na.omit um die TOTALE DEZIMIERUNG der Daten zu vermeiden lol

missing <- sapply(arrhythmia,count.na)
missing[missing>0]

arrhythmia$J <- rep(0,nrow(arrhythmia))


## loesche Zeilen mit NA
## oder approximieren
#arrhythmia <- na.omit(arrhythmia)
## dim(arrhythmia) [1] 420 280
arrhythmia <- na.omit(arrhythmia)

## Channel Zusammenfassung  ##

varName <- c("Qwidth","Rwidth","Swidth","R_width","S_width","nIntrinsDefl")

varName1 <- c("D1_aVL","D23_aVF","V1_V2","V3_V4","V5_V6","aVR")

varName2 <- c("JJamp","Qamp","Ramp","Samp","R_amp","S_amp","Pamp","Tamp","QRSAamp","QRSTAamp")

varName3 <- c("Rraggged","Rderiv","Praggged","Pderiv","Traggged","Tderiv")

## in der Praxis nimmt eine Aussage erst fuer wahr wenn sie von 2 Channel bestaetigt wird
evalTrue <- function(df){
  temp <-rowSums(as.matrix(df))
  vec <- c()
  for(i in 1:length(temp)){
    if(temp[i] >= 2){
      vec[i] <- 1
    }
    else{
      vec[i] <- 0
    }
  }
  return(vec)
}
## Zusammenfassung der fachlich zusammenhaengende Channels und Umbenennung der aVR bezogenen Zeiten
for(i in 1:10){
  if(i <= 6){
    ## Zusammenfassung der Wellenbreiten
    arrhythmia$x <-rowMeans(as.matrix(arrhythmia[c(15+i,63+i)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[1],".",varName[i],".mean")
    
    arrhythmia$x <-rowMeans(as.matrix(arrhythmia[c(27+i,39+i,75+i)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[2],".",varName[i],".mean")
    
    arrhythmia$x <-rowMeans(as.matrix(arrhythmia[c(87+i,99+i)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[3],".",varName[i],".mean")
    
    arrhythmia$x <-rowMeans(as.matrix(arrhythmia[c(111+i,123+i)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[4],".",varName[i],".mean")
    
    arrhythmia$x <-rowMeans(as.matrix(arrhythmia[c(135+i,147+i)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[5],".",varName[i],".mean")
    
    arrhythmia$x <-arrhythmia[,(51+i)]
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[6],".",varName[i])
  }
    ## Zusammenfassung der Amplituden
    
    ## wir berechnen das betragsmässige Maximum fuer jede Gruppierung da die Amplituden abhaengig vom Kanal und Individuum neg bzw pos sein koennen
    
    arrhythmia$x <- arrhythmia[c(159+i,199+i)][cbind(1:nrow(arrhythmia), max.col(abs(arrhythmia[c(159+i,199+i)])))]
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[1],".",varName2[i],".sigMax")
    
    arrhythmia$x <- arrhythmia[c(169+i,179+i,209+i)][cbind(1:nrow(arrhythmia), max.col(abs(arrhythmia[c(169+i,179+i,209+i)])))]
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[2],".",varName2[i],".sigMax")
    
    arrhythmia$x <- arrhythmia[c(219+i,229+i)][cbind(1:nrow(arrhythmia), max.col(abs(arrhythmia[c(219+i,229+i)])))]
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[3],".",varName2[i],".sigMax")
    
    arrhythmia$x <- arrhythmia[c(239+i,249+i)][cbind(1:nrow(arrhythmia), max.col(abs(arrhythmia[c(239+i,249+i)])))]
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[4],".",varName2[i],".sigMax")
    
    arrhythmia$x <- arrhythmia[c(259+i,269+i)][cbind(1:nrow(arrhythmia), max.col(abs(arrhythmia[c(259+i,269+i)])))]
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[5],".",varName2[i],".sigMax")
    
    arrhythmia$x <-arrhythmia[,(189+i)]
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName1[6],".",varName2[i])
}
for(i in 1:6){
  ## Auswertung der binaeren Aussagen wird getrennt ausgefuehrt um sie leichter abtrennen zu koennen
  arrhythmia$x <- evalTrue(arrhythmia[seq(21+i,159,12)])
  names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varName3[i],".eval")
}

## loesche zusammengefasste Spalten und J Spalte
## dim(arrhythmia) 420 117
arrhythmia[,c(16:279)] <- NULL
arrhythmia$J <- NULL


## entferne Spalten mit Standardabweichung = 0
## dim(arrhythmia) 420  113
arrhythmia.omit <- arrhythmia[,apply(arrhythmia, 2, sd)==0] 
arrhythmia <- arrhythmia[,apply(arrhythmia, 2, sd)>0]

## Symmetrie Transformation
## Eine Symmetrie Transformation der Merkmale die eine besonders schiefe Univariate Verteilung aufweisen ist sinnvoll
## umm eine leichtere Interpretation der Merkmale in der Cluster- und Hauptkomponentenanalyse
## zu gewaehrleisten und die Unterscheidung mithilfe der Distanzmatrix im Bereichn niedriger Werte 
## zu verbessern

par(mfrow = c(2, 3)) 

## binaere Werte uninteressant daher Ausschluss
arrhythmia.num <- arrhythmia[,c(1,3:(length(arrhythmia)-6))]
arrhythmia.bin <- arrhythmia[,-c(1,3:(length(arrhythmia)-6))]

for(i in 1:ncol(arrhythmia.num)){
  col <- paste0(names(arrhythmia.num[i])," distribution")
  hist(arrhythmia.num[,i],breaks=12,col="red",main=col)
}
## beim Hoehen Merkmal sind zwei eindeutige falsche Daten erkennbar Koerpergroesse 608 780  
## fachlicher Ausschluss nach nochmaligem Test mit Rohdaten rawdata[(rawdata$V3>300),]
arrhythmia.bin <- arrhythmia.bin[!(arrhythmia.num$Height>250),]
arrhythmia.num <- arrhythmia.num[!(arrhythmia.num$Height>250),] ## dim(arrhythmia) 418  99

col <- paste0(names(arrhythmia.num[2])," distribution")
hist(arrhythmia.num[,2],col="red",main=col)

## im Plot sind links- und rechtsschiefe Verteilungen zu erkennen
## bevor weitere Transformationen durchgefuehrt werden sollten jedoch Merkmale mit negativen Auspraegungen
## in positive Werte uebersetzt werden mithilfe einer einfachen Linearverschiebung
is.negative <- function(vec){
  for(i in vec){
    if(i < 0){
      return(TRUE)
    }
  }
  return(FALSE)
}
linScalePos <- function(df){
  for(i in 1:length(df)){
    if(is.negative(df[,i])){
      df[,i] <- df[,i]- min(df[,i])
    }
  }
  return(df)
}
arrhythmia.num <- linScalePos(arrhythmia.num)
## in den Univariaten Verteilungsplots kann man einige schiefe Verteilungen erkennen
## abhaengig vom Grade der schiefe werden sie als naechstes transformiert

## rechts Schief: stark 1/(x+1) mittel log10(x+1) schwach sqrt()
## links Schief: stark x^3 mittel x^2 schwach

hist((arrhythmia.num$QRS_duration),breaks=12,col="red")
hist(log10(arrhythmia.num$QRS_duration+1),breaks=12,col="red")
arrhythmia.num$QRS_duration <- log10(arrhythmia.num$QRS_duration+1)

hist((arrhythmia.num$P_R_interval),breaks=12,col="red")
hist(sqrt(arrhythmia.num$P_R_interval),breaks=12,col="red")
arrhythmia.num$P_R_interval <- sqrt(arrhythmia.num$P_R_interval)

hist((arrhythmia.num$T_interval),breaks=12,col="red")
hist(1/(arrhythmia.num$T_interval + 1),breaks=12,col="red")
arrhythmia.num$T_interval <- (1/(arrhythmia.num$T_interval + 1))

hist((arrhythmia.num$QRS),breaks=12,col="red")
hist((arrhythmia.num$QRS)^2,breaks=12,col="red")
arrhythmia.num$QRS <- (arrhythmia.num$QRS)^2

hist((arrhythmia.num$P),breaks=12,col="red")
hist((arrhythmia.num$P)^2,breaks=12,col="red")
arrhythmia.num$P <- (arrhythmia.num$P)^2

hist((arrhythmia.num$D1_aVL.Qwidth.mean),breaks=12,col="red")
hist(log10(arrhythmia.num$D1_aVL.Qwidth.mean +1 ),breaks=12,col="red")
arrhythmia.num$D1_aVL.Qwidth.mean <- log10(arrhythmia.num$D1_aVL.Qwidth.mean +1 )

hist((arrhythmia.num$D23_aVF.Qwidth.mean),breaks=12,col="red")
hist(log10(arrhythmia.num$D23_aVF.Qwidth.mean +1),breaks=12,col="red")
arrhythmia.num$D23_aVF.Qwidth.mean <- log10(arrhythmia.num$D23_aVF.Qwidth.mean +1)

hist((arrhythmia.num$aVR.Qwidth),breaks=12,col="red")
hist(log(arrhythmia.num$aVR.Qwidth+1),breaks=12,col="red")
arrhythmia.num$aVR.Qwidth <- log(arrhythmia.num$aVR.Qwidth+1)

hist((arrhythmia.num$D1_aVL.Swidth.mean),breaks=12,col="red")
hist(sqrt(arrhythmia.num$D1_aVL.Swidth.mean),breaks=12,col="red")
arrhythmia.num$D1_aVL.Swidth.mean <- sqrt(arrhythmia.num$D1_aVL.Swidth.mean)

hist((arrhythmia.num$V5_V6.Swidth.mean),breaks=12,col="red")
hist(sqrt(arrhythmia.num$V5_V6.Swidth.mean),breaks=12,col="red")
arrhythmia.num$V5_V6.Swidth.mean <- sqrt(arrhythmia.num$V5_V6.Swidth.mean)

hist((arrhythmia.num$D1_aVL.Ramp.sigMax),breaks=12,col="red")
hist(sqrt(arrhythmia.num$D1_aVL.Ramp.sigMax),breaks=12,col="red")
arrhythmia.num$D1_aVL.Ramp.sigMax <- sqrt(arrhythmia.num$D1_aVL.Ramp.sigMax)

hist((arrhythmia.num$D23_aVF.Ramp.sigMax),breaks=12,col="red")
hist(sqrt(arrhythmia.num$D23_aVF.Ramp.sigMax),breaks=12,col="red")
arrhythmia.num$D23_aVF.Ramp.sigMax <- sqrt(arrhythmia.num$D23_aVF.Ramp.sigMax)

hist((arrhythmia.num$V1_V2.Ramp.sigMax),breaks=12,col="red")
hist(log10(arrhythmia.num$V1_V2.Ramp.sigMax+1),breaks=12,col="red")
arrhythmia.num$V1_V2.Ramp.sigMax <- log10(arrhythmia.num$V1_V2.Ramp.sigMax+1)

hist((arrhythmia.num$V3_V4.Ramp.sigMax),breaks=12,col="red")
hist(sqrt(arrhythmia.num$V3_V4.Ramp.sigMax),breaks=12,col="red")
arrhythmia.num$V3_V4.Ramp.sigMax <- sqrt(arrhythmia.num$V3_V4.Ramp.sigMax)

hist((arrhythmia.num$D1_aVL.Samp.sigMax),breaks=12,col="red")
hist((arrhythmia.num$D1_aVL.Samp.sigMax)^3,breaks=12,col="red")
arrhythmia.num$D1_aVL.Samp.sigMax <- (arrhythmia.num$D1_aVL.Samp.sigMax)^3

hist((arrhythmia.num$D23_aVF.Samp.sigMax),breaks=12,col="red")
hist((arrhythmia.num$D23_aVF.Samp.sigMax)^3,breaks=12,col="red")
arrhythmia.num$D23_aVF.Samp.sigMax <- (arrhythmia.num$D23_aVF.Samp.sigMax)^3

hist((arrhythmia.num$V1_V2.Samp.sigMax),breaks=12,col="red")
hist((arrhythmia.num$V1_V2.Samp.sigMax)^2,breaks=12,col="red")
arrhythmia.num$V1_V2.Samp.sigMax <- (arrhythmia.num$V1_V2.Samp.sigMax)^2

hist((arrhythmia.num$V3_V4.Samp.sigMax),breaks=12,col="red")
hist((arrhythmia.num$V3_V4.Samp.sigMax)^3,breaks=12,col="red")
arrhythmia.num$V3_V4.Samp.sigMax <- (arrhythmia.num$V3_V4.Samp.sigMax)^3

hist((arrhythmia.num$V5_V6.Samp.sigMax),breaks=12,col="red")
hist((arrhythmia.num$V5_V6.Samp.sigMax)^3,breaks=12,col="red")
arrhythmia.num$V5_V6.Samp.sigMax <- (arrhythmia.num$V5_V6.Samp.sigMax)^3

hist((arrhythmia.num$D1_aVL.nIntrinsDefl.mean),breaks=12,col="red")
hist(sqrt(arrhythmia.num$D1_aVL.nIntrinsDefl.mean),breaks=12,col="red")
arrhythmia.num$D1_aVL.nIntrinsDefl.mean <- sqrt(arrhythmia.num$D1_aVL.nIntrinsDefl.mean)

hist((arrhythmia.num$V1_V2.nIntrinsDefl.mean),breaks=12,col="red")
hist(log(arrhythmia.num$V1_V2.nIntrinsDefl.mean +1),breaks=12,col="red")
arrhythmia.num$V1_V2.nIntrinsDefl.mean <- log(arrhythmia.num$V1_V2.nIntrinsDefl.mean+1)

hist((arrhythmia.num$D23_aVF.QRSAamp.sigMax),breaks=12,col="red")
hist((arrhythmia.num$D23_aVF.QRSAamp.sigMax)^2,breaks=12,col="red")
arrhythmia.num$D23_aVF.QRSAamp.sigMax <- (arrhythmia.num$D23_aVF.QRSAamp.sigMax)^2

hist((arrhythmia.num$aVR.QRSAamp),breaks=12,col="red")
hist((arrhythmia.num$aVR.QRSAamp)^2,breaks=12,col="red")
arrhythmia.num$aVR.QRSAamp <- (arrhythmia.num$aVR.QRSAamp)^2

hist((arrhythmia.num$D1_aVL.Tamp.sigMax),breaks=12,col="red")
hist((arrhythmia.num$D1_aVL.Tamp.sigMax)^2,breaks=12,col="red")
arrhythmia.num$D1_aVL.Tamp.sigMax <- (arrhythmia.num$D1_aVL.Tamp.sigMax)^2

hist((arrhythmia.num$aVR.QRSAamp),breaks=12,col="red")
hist((arrhythmia.num$aVR.QRSAamp)^2,breaks=12,col="red")
arrhythmia.num$aVR.QRSAamp <- (arrhythmia.num$aVR.QRSAamp)^2

hist((arrhythmia.num$V1_V2.QRSTAamp.sigMax),breaks=12,col="red")
hist(sqrt(arrhythmia.num$V1_V2.QRSTAamp.sigMax+1),breaks=12,col="red")
arrhythmia.num$V1_V2.QRSTAamp.sigMax <- sqrt(arrhythmia.num$V1_V2.QRSTAamp.sigMax)

hist((arrhythmia.num$aVR.QRSTAamp),breaks=12,col="red")
hist((arrhythmia.num$aVR.QRSTAamp)^2,breaks=12,col="red")
arrhythmia.num$aVR.QRSTAamp <- (arrhythmia.num$aVR.QRSTAamp)^2

trueclus <- arrhythmia.num$class
arrhythmia.num <-scale(arrhythmia.num)
arrhythmia.num <- as.data.frame(arrhythmia.num)
arrhythmia.num$class  <- trueclus

arrhythmia <- cbind(arrhythmia.num,arrhythmia.bin)

## ein Ausreisser Test waere eventuell angebracht aber aufgrund des fehlenden Sachverstaendnisses fuer die Thematik und der bereits geringen
## Datensatzgroesse belasse ich es hier bei einer haendischen Filterung

write.csv(arrhythmia,paste0("C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt/Daten/arrhythmia_cat_",1,".csv"), row.names = FALSE)
