require("psych")
require("AER")
require("mlbench")
require("cluster")
require("zoo")
require("plyr")
require(data.table)
require(factoextra)
require(matrixStats)

##  Daten einlesen ##
##setwd("C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt")
arrhythmia <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/arrhythmia/arrhythmia.data", header = FALSE, na.strings = "?" )
rawdata <- arrhythmia
arrhythmia <- rawdata
## allgemeine Strukturanalyse
'names(arrhythmia)
summary(arrhythmia)
describe(arrhythmia)
str(arrhythmia)'

## ergaenze relevante Spaltennamen
attnames <- c("Age","Sex","Height","Weight","QRS duration","P-R interval",
              "Q-T interval","T interval","P interval","QRS","T","P","QRST","J","Heart rate","class")
names(arrhythmia)[c(1:15,280)] <- attnames 

##  Datenbereinigung  ##

# Missing Value Analyse
table(is.na(arrhythmia))
table(complete.cases(arrhythmia))
count.na <- function(x){sum(is.na(x))} 

# J 376 Missing Values => wird gedroppt
missing <- sapply(arrhythmia,count.na)
missing[missing>0]

# lösche V14 nach channel Zusammenfassung um Verrueckung zu vermeiden
arrhythmia$J <- rep(0,nrow(arrhythmia))
# loesche Zeilen mit NA
arrhythmia <- na.omit(arrhythmia)

## Channel Zusammenfassung  ##

rowmax <- function(data,x){
  apply(data,1,function(x) max(x))
}
rowmin <- function(data,x){
  apply(data,1,function(x) min(x))
}
rowVars(as.matrix(arrhythmia[seq(16,159,12)]))

rowMeans
varname <- c("Qwidth","Rwidth","Swidth","R_width","S_width","int_num",
               "Rraggged","Rderiv","Praggged","Pderiv","Traggged","Tderiv")

varname2 <- c("JJamp","Qamp","Ramp","Samp","R_amp","S_amp","Pamp","Tamp","QRSAamp","QRSTAamp")

for(i in 1:12){
  if(i < 7){
    arrhythmia$x <-rowMaxs(as.matrix(arrhythmia[seq(15+i,159,12)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname[i],".max")
    
    arrhythmia$x <-rowMins(as.matrix(arrhythmia[seq(15+i,159,12)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname[i],".min")
    
    arrhythmia$x <-rowVars(as.matrix(arrhythmia[seq(15+i,159,12)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname[i],".var")
    
    arrhythmia$x <-rowMeans(as.matrix(arrhythmia[seq(15+i,159,12)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname[i],".mean")
  }
  else{
    arrhythmia$x <-rowMaxs(as.matrix(arrhythmia[seq(15+i,159,12)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname[i],".max")
  }
  if(i < 11){
    arrhythmia$x <-rowMaxs(as.matrix(arrhythmia[seq(159+i,279,10)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname2[i],".max")
    
    arrhythmia$x <-rowMins(as.matrix(arrhythmia[seq(159+i,279,10)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname2[i],".min")
    
    arrhythmia$x <-rowVars(as.matrix(arrhythmia[seq(159+i,279,10)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname2[i],".var")
    
    arrhythmia$x <-rowMeans(as.matrix(arrhythmia[seq(159+i,279,10)]))
    names(arrhythmia)[names(arrhythmia) == "x"] <- paste0(varname[i],".mean")
  }
}

arrhythmia[,c(16:279)] <- NULL
arrhythmia$J <- NULL

# entferne Spalten mit Standardabweichung = 0
arrhythmia.omit <- arrhythmia[,apply(arrhythmia, 2, sd)==0] 
arrhythmia <- arrhythmia[,apply(arrhythmia, 2, sd)>0]

write.csv(arrhythmia,paste0(  "C:/Users/David/Documents/Uni/Explorative Datenanalyse/1. Projekt/arrhythmia_cat_",1,".csv"), row.names = FALSE)
