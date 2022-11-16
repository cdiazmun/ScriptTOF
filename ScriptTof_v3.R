# This script is based on ScriptTof_v2.R from Hannes DeCadt and it is updated by Cristian Díaz-Muñoz

# This script takes the output files from the GC_TOF-MS and process them doing:
# 1) Reshape the dataframe and get the best hit for each RT.
# 2) Flag all duplicated compounds and merge all files together.
# 3) Eliminate those compounds present in only 1 out of 3 samples.
# 4) Calculate the average of the triplicates.

# Your samples need to be in a separate folder inside your working directory (Sample reports)
# Your samples need to be named as follows: 
# sample1_a.xls, sample1_b.xls, sample1_c.xls,
# sample2_a.xls, sample2_b.xls, sample2_c.xls, etc.

# Import all needed libraries to run the script
library(readxl) # To open excel files
library(tidyverse) # Type of R syntax, to make code easier
library(reshape2)

# Set the working directory
setwd("C:/Users/cdiazmun/OneDrive - Vrije Universiteit Brussel/IMDO/R scripts IMDO/GC-TOF-MS/")
samples_wd <- "C:/Users/cdiazmun/OneDrive - Vrije Universiteit Brussel/IMDO/R scripts IMDO/GC-TOF-MS/Sample reports/"

# List all files in the folder (all your samples)
files <- list.files(path = samples_wd, pattern=".xls")
names <- as.vector(sapply(strsplit(basename(files), ".xls"), `[`, 1))


##########################################################
#1) Reshape the dataframe and get the best hit for each RT
##########################################################

for (i in 1:length(files)) {
  # Import all files
  df <- read_excel(path = paste("Sample reports/",files[i], sep = ""), skip=32) # 32 first rows are always empty, so can be skipped
  # Select the RT rows
  rt <- as.numeric(df[,12]=="DBC TIC") 
  # In which position they are
  index <- which(rt %in% 1) 
  # Create a data frame with all information needed
  df1 <- data.frame(Compound = as.data.frame(df[index + 1, 1]), 
                    CAS = as.data.frame(df[index + 1, 17]),
                    RT = as.data.frame(df[index, 1]),
                    MF = as.data.frame(df[index + 1, 19]),
                    area = as.data.frame(df[index, 22])
  )
  # Modify the column names
  colnames(df1) <- c("Compound", "CAS", "RT", "MF", "Area")
  # Write the table in tsv format (easy to work with in R and Excel)
  write.table(file = paste("Sample reports/",names[i],".txt", sep = ""), 
              x = df1, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
}


##############################################################
#2) Flag all duplicated compounds and merge all files together
##############################################################

# List all files in the folder (all your samples)
files2 <- list.files(path = samples_wd, pattern=".txt")

df_all <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df_all) <- c("Compound", "CAS", "RT", "MF", "Area", "Sample","Duplicated")

for (i in 1:length(files2)) {
  # Import all the files
  df1 <- read.delim(file = paste("Sample reports/",files2[i], sep = ""))
  # Put the sample name
  df1$Sample <- rep(names[i], nrow(df1))
  # Some compounds can be found multiple times: list all of them up
  df1$Duplicated <- duplicated(df1$Compound)
  # In which position they are
  index2 <- which(duplicated(df1$Compound))
  # Add an '*' to every compound that is found twice (two *'s if it is found 3 times, etc.)
  while (length(index2)>0){ 
    for (k in index2){
      df1$Compound[k] <- str_c(df1$Compound[k],"*")
    }
    index2 <- which(duplicated(df1$Compound))
  }
  df_all <- rbind(df_all,df1)
}

# Write the table in tsv format (easy to work with in R and Excel)
write.table(file = "Output reports/all_samples.txt", 
            x = df_all, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


################################################################
#3) Eliminate those compounds present in only 1 out of 3 samples
################################################################

# Add a column with the triplicate info (change according to the name of the file)
df_tri <- df_all
df_tri$Rep <- sapply(strsplit(basename(as.character(df_tri$Sample)), "_"), `[`, 2)
df_tri$Sample <- sapply(strsplit(basename(as.character(df_tri$Sample)), "_"), `[`, 1)

# Calculate the standard deviation for each compound in the triplicate. 
df_tri <- df_tri %>%
  group_by(Compound, Sample) %>%
  mutate(Area_SD = sd(Area))

# Eliminate all the compounds with non-existing SD (present in only 1 sample)
df_cur <- subset(df_tri, df_tri$Area_SD > 0)

# Count how many duplicates there are still in the samples
length(df_all$Compound[df_all$Duplicated == TRUE])
length(df_cur$Compound[df_cur$Duplicated == TRUE])

# Write the table in tsv format (easy to work with in R and Excel)
write.table(file = "Output reports/all.samples.curated.txt", 
            x = df_cur, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


##############################################
# 4) Calculate the average of the triplicates
##############################################

# Calculate the mean for every numeric column and remove the replicates
df_av <- df_cur %>%
  group_by(Compound, Sample) %>%
  mutate(RT_SD = sd(RT)) %>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

# Write the table in tsv format (easy to work with in R and Excel)
write.table(file = "Output reports/all.samples.curated.averaged.txt", 
            x = df_av, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
