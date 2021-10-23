require("psych")
require("AER")
require("mlbench")
require("cluster")
require("zoo")
require("plyr")
require(data.table)
rq
install.packages("dgCMatrix") 
##  Daten einlesen ##
##setwd("C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt")
arrhythmia <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/arrhythmia/arrhythmia.data", header = FALSE, na.strings = "?" )
rawdata <- arrhythmia
arrhythmia <- rawdata
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

# lösche V14 nach channel Zusammenfassung um Verrueckung zu vermeiden
arrhythmia$V14 <- rep(0,nrow(arrhythmia))
# loesche Zeilen mit NA
arrhythmia <- na.omit(arrhythmia)

# Fasse channel Zusammen

conc <- function(data,FUN, col){
  old_names <- names(data)
  x <- FUN(data[,col])
  A <- cbind(data,x)
  var <- paste0("V",ncol(A))
  colnames(A)[ncol(A)] <- var
  return(A)
}

rowmax <- function(data,x){
  apply(data,1,function(x) max(x))
}
Qwidth.col <- seq(16,159,12)
Rwidth.col <- seq(17,159,12)
Swidth.col <- seq(18,159,12)
R_width.col <- seq(19,159,12)
S_width.col <- seq(20,159,12)

JJamp.col <- seq(160,279,10)
Qamp.col <- seq(161,279,10)
Ramp.col <- seq(162,279,10)
Samp.col <- seq(163,279,10)
R_amp.col <- seq(164,279,10)
S_amp.col <- seq(165,279,10)
Pamp.col <- seq(166,279,10)
Tamp.col <- seq(167,279,10)
QRSAamp.col <- seq(168,279,10)
QRSTAamp.col <- seq(169,279,10)

arrhythmia <- conc(arrhythmia,rowMeans,Qwidth.col) 
arrhythmia <- conc(arrhythmia,rowMeans,Rwidth.col) 
arrhythmia <- conc(arrhythmia,rowMeans,Swidth.col) 
arrhythmia <- conc(arrhythmia,rowMeans,R_width.col) 
arrhythmia <- conc(arrhythmia,rowMeans,S_width.col) 

arrhythmia <- conc(arrhythmia,rowmax,JJamp.col) 
arrhythmia <- conc(arrhythmia,rowmax,Qamp.col) 
arrhythmia <- conc(arrhythmia,rowmax,Ramp.col) 
arrhythmia <- conc(arrhythmia,rowmax,Samp.col) 

arrhythmia <- conc(arrhythmia,rowmax,R_amp.col) 
arrhythmia <- conc(arrhythmia,rowmax,S_amp.col) 
arrhythmia <- conc(arrhythmia,rowmax,Pamp.col) 
arrhythmia <- conc(arrhythmia,rowmax,Tamp.col) 
arrhythmia <- conc(arrhythmia,rowmax,QRSAamp.col) 
arrhythmia <- conc(arrhythmia,rowmax,QRSTAamp.col) 

'arrhythmia$Qwidth <- rowmax(arrhythmia[seq(16,159,12)])
arrhythmia$Rwidth <- rowmax(arrhythmia[seq(17,159,12)])
arrhythmia$Swidth <- rowmax(arrhythmia[seq(18,159,12)])
arrhythmia$R_width <- rowmax(arrhythmia[seq(19,159,12)])
arrhythmia$S_width <- rowmax(arrhythmia[seq(20,159,12)])
arrhythmia$int_num<- rowmax(arrhythmia[seq(21,159,12)])

arrhythmia$Rraggged <- rowmax(arrhythmia[seq(22,159,12)])
arrhythmia$Rderiv <- rowmax(arrhythmia[seq(23,159,12)])
arrhythmia$Praggged <- rowmax(arrhythmia[seq(24,159,12)])
arrhythmia$Pderiv <- rowmax(arrhythmia[seq(25,159,12)])
arrhythmia$Traggged <- rowmax(arrhythmia[seq(26,159,12)])
arrhythmia$Tderiv <- rowmax(arrhythmia[seq(27,159,12)])

arrhythmia$JJamp <- rowmax(arrhythmia[seq(160,279,10)])
arrhythmia$Qamp <- rowmax(arrhythmia[seq(161,279,10)])
arrhythmia$Ramp <- rowmax(arrhythmia[seq(162,279,10)])
arrhythmia$Samp <- rowmax(arrhythmia[seq(163,279,10)])
arrhythmia$R_amp <- rowmax(arrhythmia[seq(164,279,10)])
arrhythmia$S_amp <- rowmax(arrhythmia[seq(165,279,10)])
arrhythmia$Pamp <- rowmax(arrhythmia[seq(166,279,10)])
arrhythmia$Tamp <- rowmax(arrhythmia[seq(167,279,10)])
arrhythmia$QRSAamp <- rowmax(arrhythmia[seq(168,279,10)])
arrhythmia$QRSTAamp <- rowmax(arrhythmia[seq(169,279,10)])'

arrhythmia[,c(Qwidth.col,Rwidth.col,Swidth.col,R_width.col,S_width.col,
              160:279)] <- NULL

# entferne Spalten mit Standardabweichung = 0
arrhythmia.omit <- arrhythmia[,apply(arrhythmia, 2, sd)==0] 
arrhythmia <- arrhythmia[,apply(arrhythmia, 2, sd)>0]


write.csv(arrhythmia,paste0(  "C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt/arrhythmia_cat_",1,".csv"), row.names = FALSE)
