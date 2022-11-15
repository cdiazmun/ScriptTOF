#install.packages(readxl)   # install.packages is only needed once, if you have not installed these packages yet
#install.packages(tidyverse)
#install.packages(openxlsx)

library(readxl)  # every time you open R you have to use the library function (to activate the packages)
library(tidyverse)
library(openxlsx)

setwd("C:/Users/hadeca/Vrije Universiteit Brussel/IMDO - Hannes Decadt/Documents/overig/DarioTofScript") # adapt to the folder where your peak area xls files are saved

namen<-list.files(pattern = ' .xls')

alles<-c() # create 
t<-1 # to count all the xls files, increases with every new file processed
for (j in namen) { 
  M<-as.data.frame(read_xls(j, sheet = 2, skip=7)) # skip 7: the first 7 rows don't contain data and can be skipped
  
  # remove siloxanes in 2 steps
  os <- c(grep('sil', M[,1]),grep('Sil', M[,1])) # find all siloxanes 
  M<-M[-os,]  #remove all rows containing siloxanes 
  
  # find compounds that are more than one time found in one sample and give them a '2'
  d<-duplicated(M[,1]);    w<-which(d)
  while (length(w)>0){for (k in w){ M[k,1]<-str_c(M[k,1],"2")}; d<-duplicated(M[k,1]);    w<-which(d)}
  
  sampleX<-data.frame(compound = M[,1], RT = as.numeric(M[,4]),area = as.numeric(M[,13]), sample = t); t<-t+1 # in column 1 are the compound names, in column 4 is the RT, in column 14 the area
  alles<-rbind(alles,sampleX)
}

co<-unique(alles$compound) # find all unique compounds
co2<-as.data.frame(co)
# make matrix containing peak areas (m) and RTs (n)
m <- data.frame(matrix(0, ncol = length(co), nrow = length(namen)),row.names = namen); colnames(m)<-co
n<-m

#transform alles in a useful matrix
for (za in 1:length(co)){
  gevonden<-which(alles$compound==co2[za,1])
  m[alles$sample[gevonden],za]<-alles$area[gevonden]
  n[alles$sample[gevonden],za]<-alles$RT[gevonden]  }

write.xlsx(m,'Areas.xlsx',rowNames=TRUE)
write.xlsx(n,'RTs.xlsx',rowNames=TRUE)

