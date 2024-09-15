library(DESeq2)
library(tidyverse)
library(dplyr)

#Code written with assistance from online tutorials by Khushbu Patel (Bioinformagician: https://www.youtube.com/@Bioinformagician)
#Code written in collaboration with Professor Walden Ai, University of South Carolina

#Read dataset
Counts_data<-read.table('C:/Users/Walden.Ai/Documents/R analysis/Gencodecounts.txt')
names(Counts_data) <- Counts_data[1,]
Counts_data <- Counts_data[-1,]

#save geneid names for later
write.csv(Counts_data$Geneid, 'C:/Users/Walden.Ai/Documents/R analysis/Geneid.csv')

#Read conditions data file
colData <- read.csv("C:/Users/Walden.Ai/Documents/R analysis/Co_culture_conditions2.csv", header = TRUE, 
                    fill = TRUE)

#subset Counts_data to have just the samples
subset_cols <- c('S01Aligned.bam', 'S02Aligned.bam', 'S03Aligned.bam', 'S04Aligned.bam', 
                 'S05Aligned.bam', 'S06Aligned.bam')
Counts_data <- Counts_data[subset_cols]


#converting from character to number
Counts_data$S01Aligned.bam <- as.numeric(as.character(unlist(Counts_data$S01Aligned.bam)))

Counts_data$S02Aligned.bam <- as.numeric(as.character(unlist(Counts_data$S02Aligned.bam)))

Counts_data$S03Aligned.bam <- as.numeric(as.character(unlist(Counts_data$S03Aligned.bam)))

Counts_data$S04Aligned.bam <- as.numeric(as.character(unlist(Counts_data$S04Aligned.bam)))

Counts_data$S05Aligned.bam <- as.numeric(as.character(unlist(Counts_data$S05Aligned.bam)))

Counts_data$S06Aligned.bam <- as.numeric(as.character(unlist(Counts_data$S06Aligned.bam)))

#Creating Condition_medium dataframe 
subset_cols <- c('S01Aligned.bam', 'S02Aligned.bam', 'S03Aligned.bam', 'S04Aligned.bam')
Conditioned_medium <- Counts_data[subset_cols]


Conditioned_medium_coldata <- colData %>% filter(
  Samples == "S01Aligned.bam" | Samples == "S02Aligned.bam" | Samples == 'S03Aligned.bam' |
    Samples == 'S04Aligned.bam')

Conditioned_medium_coldata <- Conditioned_medium_coldata %>% remove_rownames %>% column_to_rownames(var="Samples")

#Create dds object
dds <- DESeqDataSetFromMatrix(countData = round(Conditioned_medium),
                              colData = Conditioned_medium_coldata,
                              design = ~ Co_culture)

# set the factor level
dds$Co_culture <- relevel(dds$Co_culture, ref = "control")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res
write.csv(res,"C:/Users/Walden.Ai/Documents/R analysis/Output(July 29).csv")
