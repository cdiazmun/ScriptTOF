library(readxl)  
library(tidyverse)
library(openxlsx)

# all samples are treated independently, you can process yourself the duplicates/triplicates afterwards

setwd("C:/Users/hadeca/Vrije Universiteit Brussel/IMDO - Hannes Decadt/Documents/Methodes en resultaten/metabolites/SPME-GC-MS/Pekel")
pad<-"C:/Users/hadeca/Vrije Universiteit Brussel/IMDO - Hannes Decadt/Documents/Methodes en resultaten/metabolites/SPME-GC-MS/Pekel"

lijst<-list.files(path=pad,pattern=".xls")

besta<-read_xls(lijst[1], skip=32) # skip 32: 32 first rows are always empty, so can be skipped

#now only the best match is chosen
u<-as.numeric(besta[,12]=="DBC TIC")
index<-which(u %in% 1)
b1<-data.frame(compound = as.data.frame(besta[index+1,1]), RT = as.data.frame(besta[index,1]), area = as.data.frame(besta[index,22]))
colnames(b1)<-c("Compound","RT","Area")
d<-duplicated(b1$Compound);    w<-which(d) # some compounds can be found multiple times: list all of them up
while (length(w)>0){ #while loop: adds a '2' to every compound that is found twice (two 2's if it is found 3 times etc)
  for (k in w){ b1$Compound[k]<-str_c(b1$Compound[k],"2")}; d<-duplicated(b1$Compound);    w<-which(d)
}
colnames(b1)<-c("Compound","RT_1","Area_1")

for (j in c(2:length(lijst))) { # j goes to every first xls sheet of each cheese sample (default: if you have 3 different cheese samples, you make for each sample 3 different vials and start with all first vials of all cheese samples, then all second, then all thirth. i is default 3, number of vials run of same cheese sample
  
  besta<-read_xls(lijst[j], skip=32) # skip 32: 32 first rows are always empty, so can be skipped
  
  #now only the best match is chosen
  u<-as.numeric(besta[,12]=="DBC TIC")
  index<-which(u %in% 1)
  b2<-data.frame(compound = as.data.frame(besta[index+1,1]), RT = as.data.frame(besta[index,1]), area = as.data.frame(besta[index,22]))
  colnames(b2)<-c("Compound","RT","Area")
  d<-duplicated(b2$Compound);    w<-which(d)
  while (length(w)>0){for (k in w){ b2$Compound[k]<-str_c(b2$Compound[k],"2")}; d<-duplicated(b2$Compound);    w<-which(d)}
  
  colnames(b2)<-c("Compound",str_c("RT_",j),str_c("Area_",j))
  
  b1<-merge(b1,b2,by='Compound',all=TRUE)
}

write.xlsx(b1,'All.xlsx',rowNames=TRUE)
