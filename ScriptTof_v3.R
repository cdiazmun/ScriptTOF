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
setwd("C:/Users/mmorenoc/Vrije Universiteit Brussel/Cristian DIAZ MUNOZ - Mery/TOF")
samples_wd <- "C:/Users/mmorenoc/Vrije Universiteit Brussel/Cristian DIAZ MUNOZ - Mery/TOF/SSR/"

# List all files in the folder (all your samples)
files <- list.files(path = samples_wd, pattern=".xls")
#names <- as.vector(sapply(strsplit(basename(files), ".xls"), `[`, 1))
names_2 <- c(rep(c("a", "b","c"),27))
names_3 <- c(rep("Captain_1", 3),rep("Captain_2", 3),rep("Captain_3", 3),
             rep("Karma_1", 3),rep("Karma_2", 3),rep("Karma_3", 3),
             rep("Thylbert_1", 3),rep("Thylbert_2", 3),rep("Thylbert_3", 3),
             rep("Rish_1", 3),rep("Rish_2", 3),rep("Rish_3", 3),
             rep("Gremline_1", 3), rep("Gremline_2", 3), rep("Gremline_3", 3),
             rep("Yaya_1", 3),rep("Yaya_2", 3),rep("Yaya_3", 3),
             rep("Smile_1", 3),rep("Smile_2", 3),rep("Smile_3", 3),
             rep("Yuguen_1", 3),rep("Yuguen_2", 3),rep("Yuguen_3", 3),
             rep("Bekombucha_1", 3),rep("Bekombucha_2", 3),rep("Bekombucha_3", 3))
names_f <- paste(names_3,names_2, sep = "_")

##########################################################
#1) Reshape the dataframe and get the best hit for each RT
##########################################################

for (i in 1:length(files)) {
  # Import all files
  df <- read_excel(path = paste("SSR/",files[i], sep = ""), skip=32) # 32 first rows are always empty, so can be skipped
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
  write.table(file = paste("SSR/",names_f[i],".txt", sep = ""), 
              x = df1, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
}


##############################################################
#2) Flag all duplicated compounds and merge all files together
##############################################################

# List all files in the folder (all your samples)
files2 <- list.files(path = samples_wd, pattern=".txt")
names_f2 <- sort(names_f)

df_all <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df_all) <- c("Compound", "CAS", "RT", "MF", "Area", "Sample","Duplicated")

for (i in 1:length(files2)) {
  # Import all the files
  df2 <- read.delim(file = paste("SSR/",files2[i], sep = ""))
  # Put the sample name
  df2$Sample <- rep(names_f2[i], nrow(df2))
  # Some compounds can be found multiple times: list all of them up
  df2$Duplicated <- duplicated(df2$Compound)
  # In which position they are
  index2 <- which(duplicated(df2$Compound))
  # Add an '*' to every compound that is found twice (two *'s if it is found 3 times, etc.)
  while (length(index2)>0){ 
    for (k in index2){
      df2$Compound[k] <- str_c(df2$Compound[k],"*")
    }
    index2 <- which(duplicated(df2$Compound))
  }
  df_all <- rbind(df_all,df2)
}

# Write the table in tsv format (easy to work with in R and Excel)
write.table(file = "Output reports/all_samples.txt", 
            x = df_all, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

# Check IS (Toluene-D8)
check_list <- df_all$Sample[df_all$Compound == "Toluene-D8"]
print(paste("There are ", length(files2) - length(check_list), " files with no Toluene-D8"))

################################################################
#3) Eliminate those compounds present in only 1 out of 3 samples
################################################################

# Add a column with the triplicate info (change according to the name of the file)
df_tri <- df_all
df_tri$Rep <- sapply(strsplit(basename(as.character(df_tri$Sample)), "_"), `[`, 3)
df_tri$Sample <- paste(sapply(strsplit(basename(as.character(df_tri$Sample)), "_"), `[`, 1),
                       sapply(strsplit(basename(as.character(df_tri$Sample)), "_"), `[`, 2),
                       sep = " ")

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


##############################################
# 5) Plot the results as a PCA
##############################################

# Create the matrix of areas with compounds in the columns and samples in the rows
tof <- df_av[c(1,2,5)] %>%
  spread(key = Compound, value = Area)

tof[is.na(tof)] <- 0

# Write the table
write.table(file = "Output reports/all.samples.curated.averaged.matrix.txt", 
            x = tof, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

# Eliminate Bekombucha and divide all columns by Toluene-D8
tof2 <- tof[-(1:3),] %>%
  mutate_if(is.numeric, funs(./`Toluene-D8`))

tof2$Bottle <- sapply(strsplit(basename(as.character(tof2$Sample)), " "), `[`, 1)

# Remove the Toluene-D8 column

tof2 <- tof2[-327]

library(ggfortify) # Package needed for plotting the PCA

# If the variances of the variables in the data vary by orders of magnitude, 
# scale is appropriate
ar.pca <- prcomp(tof2[2:(ncol(tof2)-1)],
                 center = TRUE,
                 scale. = FALSE)

# Explore the factor loadings of the PCA
print(ar.pca)
# Explore the contribution of each PC to the total variance in the dataset
screeplot(ar.pca)

# Write the factor loadings to a csv file, if you want
write.table(ar.pca,"pca_loading_factors.txt")

# Plot the PCA (score plot)
g <- autoplot(ar.pca, data = tof2, colour = "Bottle") +
  theme_bw() + scale_colour_brewer(type = "qual", palette = "Set1")

print(g)

# Save the PCA in 300 DPI, fitting publication standards
ggsave("PCA_fertemp_all_curated.tiff", width = 6, height = 4, units = "in", dpi = 300)

# Plot the factor loadings (score plot + loading plot)
h <- autoplot(ar.pca, data = na.omit(df), colour = "Brewery", shape = "Brewery", 
              loadings = 5, loadings.colour = "black", loadings.label = TRUE, 
              loadings.label.size = 3, loadings.label.colour = "Black",
              loadings.label.repel = TRUE, lable = T) +
  theme_bw()

print(h)

# Save the PCA in 300 DPI, fitting publication standards
ggsave("PCA_metabolite_loadings_all.tiff", width = 6, height = 4, units = "in", dpi = 300)






