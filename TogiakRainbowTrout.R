#Data analysis for Togiak NWR Rainbow Trout study that examines genetic diversity and population structure of RBT in SW Alaska
#Use base R and tidyverse functions

# Open packages--------------------------------------------------------------------------------------------------------------
library(adegenet)
library(allelematch)
library(ape)
library(cowplot)
library(genepop)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(gtools)
library(hierfstat)
library(pegas)
library(PopGenKit)
library(poppr)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(viridis)

# Set print to 200 lines-----------------------------------------------------------------------------------------------------
options(tibble.print_max = 200, tibble.print_min = 200)

# Open functions-------------------------------------------------------------------------------------------------------------

# Tibble to Genepop 
# Create a dataframe (tibble) with Individual and Population in first two columns followed by columns for each locus in two or three digit format ("0103", "001003").  
# The dataframe (tibble) must be organize as follows: 
#     First column contains individual identifiers and is named IndList.
#     Second column contains population identifiers and is named PopList.
#     The remaining columns contain the genotypes for each locus in either two (e.g., 0102) or three (e.g., 001002) digit code with no symbol separating alleles.
#     The locus columns are named for each locus
# The individual and population columns are both used here and one will be selected for use in GenePop (see subset step below) depending upon analysis.
Tib2Genepop <- function(x, y, z) {
  inputGP <- x
  # Creates blank dataframe with same number of columns as inputGP and no. columns - 1 rows.
  NewGP <-  data.frame(matrix(nrow = length(inputGP) - 1, ncol = length(inputGP))) 
  # Assign columns one and two of NewGP the column names of inputGP column 2 through end.  
  NewGP[,1:2] <-  names(inputGP[,2:length(inputGP)]) 
  # Assign the column names of inputGP to NewGP
  names(NewGP) <-  names(inputGP) 
  # Assign the first row, column 1 of NewGP the string "Genepop Input'.
  NewGP[1,1:2] <- 'Genepop Input' 
  # Create LociCount vector with one element that is the number of columns of NewGP minus 2.
  LociCount <-  length(NewGP)-2 
  # Add a comma to each character in the column IndList in the dataframe inputGP.  The comma is needed for genepop. 
  inputGP$IndList <- paste(inputGP$IndList,",", sep = "")
  # Add a comma to each character in the column PopList in the dataframe inputGP. The comma is needed for genepop. 
  inputGP$PopList <- paste(inputGP$PopList,",", sep = "")
  # Create a one row empty dataframe 
  Pop <-  as.data.frame(t(c(rep(NA, length(inputGP))))) 
  # Assign the column names of inputGP to Pop
  names(Pop) <-  names(inputGP) 
  # Assign the first row, column 1 and 2 of Pop the string 'Pop'
  Pop[1:2] <-  'Pop' 
  # Create a PopN vector with the names of each population from the the column inputGP$PopList in inputGP.
  PopN <-  unique(inputGP$PopList) 
  # for loop appends the datafram NewGP with the Pop dataframe and each population, in turn, from the dataframe inputGP.
  for(i in 1:length(PopN)) {
    tempGP <-  subset(inputGP, inputGP$PopList == PopN[i])
    NewGP <-  rbind(NewGP, Pop, tempGP)  
  }
  # Remove either column titled IndList or PopList from NewGP
  #NewGP <-  subset(NewGP, select = -c(PopList)) #Retaining the Indlist column for intra-watersded assignment test.
  NewGP_Ind <-  subset(NewGP, select = -c(PopList)) #Retaining the Poplist column for the TNWRRBT study
  NewGP_Pop <-  subset(NewGP, select = -c(IndList)) #Retaining the Indlist column for the TNWRRBT study.  
  # Write NewGP_Ind and NewGP_Pop to file.
  write.table(NewGP_Ind,y,quote=FALSE,row.names=F,col.names=F,sep="\t",na="")   
  write.table(NewGP_Pop,z,quote=FALSE,row.names=F,col.names=F,sep="\t",na="")    
}

#============================================================================================================================


# 1=PREPARE DATA=============================================================================================================

# 1.1-Read and arrange input file(s)-----------------------------------------------------------------------------------------
# The data files represent SNP genotypes from TaqMan assays.
# The data files are in the folder "DataFiles" in the TogiakRainbowTrout R-project.
# The data file format is from the quantstudio real time PCR instrument.
# Read input as a tibble.
# The object "InputDataDir" is the directory path to input data file folder "InputFiles" in the TogiakRainbowTrout R-project.
InputDataDir <- "./InputFiles/"
# The object "files" is a character vector of the file names in "InputFiles"
files = dir(path = InputDataDir, pattern="*.csv")
# Define input
# Read the input data files from the "InputFiles" folder in the TogiakRainbowTrout R-project and create tibble "input"
input <-  paste(InputDataDir,files,sep="") %>% 
  # Use the map function in purrr to read the input data files from the quantstudia as lists. The first 17 rows are skipped and the four columns are type "character".
  map(read_csv, skip = 17, col_types = 'cccc') %>% 
  # Use the reduce function in purrr to combine the lists (data files) by row using rbind.
  reduce(rbind) %>% 
  # rename the columns
  rename(Ind = 'Sample ID', Locus = 'Assay Name or ID', Allele1 = 'Allele 1 Call', Allele2 = 'Allele 2 Call') %>% 
  # use the mutate function in dplyr to replace periods and dashes with underscore in strings in the columns "Ind" and "Locus"
  mutate(Ind = str_replace_all(Ind, "[.]|-", "_")) %>% 
  mutate(Locus = str_replace_all(Locus, "[.]|-", "_")) %>%   
  # Use filter in dplyr to remove rows with individual in the Ind column beginning with "Bonn" or "NTC" (The Bonniville collection from CRTFC and the no template controls)
  filter(!(str_detect(Ind, "^Bonn|^NTC")))

# 1.2-Rename samples that were mis-ID'd to collection based on evaluation of GPS data in the file "TNWR_MetaData.csv created in step 4.9, below.
# The GPS data was evaluated in Google Earth Pro to indentify samples that were mis-ID'd to collection.
input <- mutate(input, Ind = str_replace_all(Ind, c("RGOO0901_090" = "RKUK0906_090", "REAG0904_195" = "REAG0904_165", "RGEC0003_001" = "RPUN0002_024", 
                                                    "RGEC0003_002" = "RPUN0002_025", "RGEC0003_003" = "RPUN0002_026", "RGEC0003_004" = "RPUN0002_027", 
                                                    "RGEC0003_005" = "RPUN0002_028")))

# 1.3-Define vector CheckScores as the distinct scores in the infile.--------------------------------------------------------  
# Use CheckScores to confirm expected scores and create vector snpScore (below) for numerical recoding of scores.
CheckScores <- pivot_longer(input, c("Allele1","Allele2"), names_to = "Allele", values_to = "Score") %>% 
  distinct(Score)

# 1.4-Compare the collection names (codes) in the data file input to those expected.----------------------------------------- 
# Highlights any possible error in code transcription during input for genotyping. 
# Collections.csv includes the expected names (codes) of the samples in input.
CollectionNames <- read.csv("Collections.csv")
# Compare the collection names (code) in the input data to the expected code (name) in CollectionNames.
CheckNames <- input %>% 
  # Use mutate to remove the underscore and individual number from the last four digits in the name in the Ind column.  The remaining name is the collection name.
  mutate(Ind = str_sub(Ind, end=-5)) %>% 
  distinct(Ind) %>% 
  # Use mutate to add column to CheckNames that shows TRUE/FALSE if a collection name (code) in input is one of the expected names (codes) in CollectionNames.
  mutate(InList = Ind %in% CollectionNames$Code) 
  # Collection names in the Ind column in CheckNames that are FALSE in the InList column should be checked for typing (transcription) errors.

# 1.5-Create vectors for numerical recoding of SNP scores.-------------------------------------------------------------------
# Vector of all possible SNP scores from the input file (see the CheckScores object above).
snpScore <- c('A','C','G','T','-','VIC','FAM','UND','INV','NOAMP','PRA') 
# Vector of numeric codes for all scores in vector "snpScore". UND, INV, NOAMP, PRA are all given the code "00".
snpCode <- c(paste('0',c(1:7),sep=''),rep('00',4))
# Create named vector linking snpCode to snpScore.
snpCon <- setNames(snpCode,snpScore)

# 1.6-Recode genotypes and remove duplicate samples.-------------------------------------------------------------------------
input <- input %>%
  # Use the mutate and recode functions to change snp scores to snp codes in "Allele1Num" and "Allele2Num" using the character vector snpCon. The 3 exclamations are needed for splicing the snpCon vector (see recode help).
  mutate(Allele1Num = recode(input$Allele1,!!!snpCon)) %>% 
  mutate(Allele2Num = recode(input$Allele2,!!!snpCon)) %>% 
  # Use the mutate and paste functions to add a column variable "Genotype" that merges the allele1 and allele 2 scores (e.g., 01 and 02 to 01/02 ).
  mutate(Genotype = paste(Allele1Num,"/",Allele2Num, sep = "")) %>% 
  # Use the mutate and str_sub functions to add a column "Collection" that identifies the collection code from the "Ind" column for each sample by removing the last four digits (underscore plus three numbers).
  mutate(Collection = str_sub(Ind, end=-5)) %>% 
  # The next two lines of code used results from CheckNames above to identify collection names (codes) that were written incorrectly. 
  # Use the mutate and str_replace_all functions to correct the collection code name for three collections (some individual samples are incorrect in dataset).
  mutate(Collection = str_replace_all(Collection, c("RSAN0914" = "RSNA0914", "TDOG0904" = "TDOG9904", "RNUG0302" = "RUNG0302",
                                                    "rbtKUS97/2__525" = "rbtKUS97", "rbtBB97_10_552" = "rbtBB97"))) %>% 
  # Use the mutate and str_replace_all functions to correct the individual code name for three collections (some individual samples are incorrect in dataset). 
  mutate(Ind = str_replace_all(Ind, c("RSAN0914" = "RSNA0914", "TDOG0904" = "TDOG9904", "RNUG0302" = "RUNG0302",
                                      "rbtKUS97/2__525" = "rbtKUS97", "rbtBB97_10_552" = "rbtBB97"))) %>% 
  # Use the mutate and duplicated functions to add a column "dups" that returns TRUE or FALSE if Ind x Locus sample is duplicated in the tibble.
  mutate(dups = duplicated(paste(Ind,Locus))) %>% 
  # Use the filter function to retain only samples with FALSE in the "dups" column (in other words, remove all duplicates).
  filter(dups == FALSE) %>% 
  # Use the select function to retain columns "Ind", "Collection", "Locus", "Genotype".
  select(Ind,Collection,Locus,Genotype)
  
# 1.7-Create lists of the Individuals and Collections and count the number of each in input.--------------------------------
Indiv <- distinct(input,Ind)
NumIndiv <- count(Indiv)
Collec <- distinct(input,Collection)
NumCollec <- count(Collec)

# 1.8-Count the number of observations for each locus and retain loci that were genotyped in all individuals------------------ 
# (this was done because "input" includes two data sets that have some different loci). 
SharedLoci <- input %>% 
  group_by(Locus) %>% 
  tally() %>% 
  filter(n == NumIndiv$n)

# 1.9-Find rows in input that have a matching Locus in SharedLoci-------------------------------------------------------------                              
#(SharedLoci lists loci that were genotyped in both raw data files).
input <- input %>% 
  semi_join(SharedLoci,by="Locus")

# 1.10-Split the Kanektok River collection into upper and lower river samples to test for differentiation----------------------
# 42 Samples are segregated into upper and lower reaches based on radio tagging results
# 8 samples are removed because they had no reach location information.
# Read csv file with reach data for Kanektok samples and replace period in Ind name with underscore
KanektokByReach <- read_csv("KanektokByReach.csv") %>% 
  mutate(Ind = str_replace_all(Ind, "[.]", "_"))
# Create named vector linking subcollection to individual for Kanektok samples
ReachCodeInd <- set_names(KanektokByReach$SubCollection,KanektokByReach$Ind)
# Recode collection code with collection subcode for Kanektok samples. Kanektok samples not assigned a subcollection are filtered out
input <- input %>% 
  mutate(Collection = recode(Ind, !!!ReachCodeInd, .default = Collection)) %>% 
  filter(!(Collection == "RKAN0901"))

# 1.11-Read the spreadsheet of TNWR collections and define TNWRcoll.---------------------------------------------------------- 
# Use cbind to add column "InDataTF" that uses binary TRUE/FALSE for collections that are/are not in inputv3.
# This step was done to confirm that the input file contained the intended collections.  
TNWRcoll <- read_csv("TNWR_Collections.csv")
TNWRcoll <- cbind(TNWRcoll, InDataTF = TNWRcoll$Collection %in% input$Collection)

# 1.12-Create named vectors linking Watershed, river code and river number to collection.-------------------------------------
Watershed <- setNames(TNWRcoll$WatershedCode,TNWRcoll$Collection)
RiverCode <- setNames(TNWRcoll$RiverCode, TNWRcoll$Collection)
RiverNo <- setNames(TNWRcoll$RiverNo, TNWRcoll$Collection)

# 1.13-Modify input to include only rows that have a matching collection in TNWRcoll.-----------------------------------------
# This step removes collections that are not part of TNWR rainbow trout study.
input <- input %>% 
  semi_join(TNWRcoll,by="Collection") %>% 
  # Use the mutate and recode functions to add a column "WatershedCode" that has the watershed code for each collection (uses the character vector Watershed and the 3 exclamations are needed for splicing the Watershed vector (see recode help).
  mutate(WatershedCode = recode(Collection, !!!Watershed)) %>% 
  mutate(RiverCode = recode(Collection, !!!RiverCode)) %>% 
  mutate(RiverNo = recode(Collection, !!!RiverNo)) %>% 
  arrange(RiverNo,Ind,Locus)

# 1.14-Identify individuals in input that are missing genotypes at more than 25% of loci.-------------------------------------
DropIndividuals <- group_by(input,Ind,Genotype,Locus) %>% 
  tally() %>% 
  summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(Genotype == "00/00", freq > 0.25)
  
# 1.15-Identify Loci in input that are missing genotypes at more than 25% of individuals.-------------------------------------
DropLoci_missingGT <- group_by(input,Locus,Genotype,Ind) %>% 
  tally() %>% 
  summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(Genotype == "00/00", freq > 0.25)

# 1.16-Identify loci in any river collection (RiverCode) that are missing genotypes at more than 25% of individuals.----------
# Did not use the code for this study.
# GenoCount00 <- group_by(input,RiverCode,Locus,Genotype) %>% 
  # tally() %>% 
  # mutate(freq = n/sum(n)) %>% 
  # filter(Genotype == "00/00", freq > 0.25)

# 1.17-Identify loci in input with a cumulative minor allele frequency (over all collections) less than 0.01.---------------
AlleleByLocus <- input %>% 
  filter(!(Genotype == "00/00")) %>% 
  separate(Genotype, into=(c("A1","A2"))) %>% 
  select(Locus, A1, A2) %>% 
  pivot_longer(c("A1","A2"), names_to = "Allele", values_to = "Score")

DropLoci_lowMAF <- AlleleByLocus %>% 
  group_by(Locus, Score) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(freq < 0.01)

# 1.18-Remove Loci and Individuals from input that are missing genotypes at 25% or more individuals or loci----------------
input <- input %>% 
  filter(!(Ind %in% DropIndividuals$Ind)) %>% 
  filter(!(Locus %in% DropLoci_missingGT$Locus)) %>%
  filter(!(Locus %in% DropLoci_lowMAF$Locus)) %>% 
  select(Ind,Collection,WatershedCode,RiverCode,RiverNo,everything())

# 1.19-Pivot input so rows are individual genotypes and columns are loci.----------------------------------------------------
input_pivot <- input %>% 
  pivot_wider(names_from = Locus, values_from = Genotype) 

saveRDS(DropLoci_missingGT, "DropLoci_missingGT.rds")
saveRDS(DropIndividuals, "DropIndividuals.rds")
saveRDS(DropLoci_lowMAF, "DropLoci_lowMAF.rds")


#============================================================================================================================


# 2=IDENTIFY AND REMOVE MATCHING GENOTYPES===================================================================================

# Use package allelematch to identify matching individuals (geneotypes) presumably sampled multiple times

# 2.1-Use the unite and then separate functions to separate alleles for each locus in input_pivot--------------------------
# into separate columns (for use in the package allelematch).
LocusNames <- sort(rep(names(input_pivot[-c(1:5)]),2))
NewLocNames <- paste(LocusNames,c("A1","A2"),sep="_")
input_Allele <- input_pivot %>% 
  unite("All",LocusNames,sep="/") %>% 
  separate(All,NewLocNames,sep="/")

# 2.2-Use the amDataset function in allelematch to produce an input dataset for allelematch from input_Allele.-------------
AM_input <- amDataset(input_Allele, missingCode = "00", indexColumn = 1, ignoreColumn = c(2:5))

# 2.3-Use the amUniqueProfile function to determine the optimum setting------------------------------------------------------ 
# for the alleleMismatch parameter (tested values 1-15).
amUniqueProfile(AM_input, alleleMismatch = c(1:15), doPlot = TRUE)

# 2.4-Use the amUnique function to identify unique genotypes (individuals) and matching individuals.------------------------- 
# (used alleleMismatch = 6 rather than 14 from amUniqueProfile.  Thus, allowing for ~ 5% allele mismatch)
PW_unique <- amUnique(AM_input, alleleMismatch = 6)
# Send PW_unique results to web browser as html (see allelematch documentation for amUnique function in R).
summary(PW_unique, html = TRUE)
# Send PW_unique results to local directory as csv file (see allelematch documentation for amUnique function in R).
summary(PW_unique, csv = "myUnique.csv")
# Define UniqueOut
UniqueOut <- read_csv("myUnique.csv")

# 2.5-Further analyze individuals that were unclassified or with multiple matches.-------------------------------------------
# (see tutorial Example 3 in allelematch supplementary documentation).
# Preview individuals that were unclassified.
select(filter(UniqueOut, rowType == "UNCLASSIFIED"), c(uniqueGroup:matchIndex,Psib))
# If there are unclassified samples then use amPairwise to evaluate those samples (see tutorial Example 3 in allelematch supplementary documentation).
# The following code is used for unclassified samples but first see the tutorial Example 3 for details.
    PW_unclassified <- amPairwise(PW_unique$unclassified, PW_unique$unique, alleleMismatch = 7)
    #Send PW_unclassified results to web browser as html to inspect unclassified individual.
    summary(PW_unclassified, html = TRUE)
# Preview individuals with multiple matches.
select(filter(UniqueOut, rowType == "MULTIPLE_MATCH"), c(uniqueGroup:matchIndex,Psib))
# Because individuals with multiple matches were present, use amPairwise to evaluate each match (see tutorial Example 3 in allelematch supplementary documentation).
PW_mismatch <- amPairwise(PW_unique$multipleMatches, PW_unique$unique, alleleMismatch = 6)
#Send PW_mismatch results to web browser as html to inspect individuals with multiple matches.
summary(PW_mismatch, html = TRUE)
# PW_mismatch results indicated that most Snake River and Eagle Creek putative matches included mismatches, not missing genotypes. Chose to retain most Snake and Eagle samples.

# 2.6-Create vector of samples to remove because of matches.-----------------------------------------------------------------
# Use the pull function to pull out a column.
SampRemove <- UniqueOut %>% 
  filter(nUniqueGroup > 1) %>% 
  select(uniqueGroup:score) %>% 
  filter(!(uniqueIndex == matchIndex)) %>% 
  filter(!(rowType == "MULTIPLE_MATCH")) %>% 
  filter(!(Psib == "NA")) %>% 
  pull(matchIndex)

saveRDS(SampRemove, "SampRemove.rds")
SampRemove <- readRDS("SampRemove.rds")

# 2.7-Remove matching samples from input and input_pivot.--------------------------------------------------------------------
#Use filter function to remove samples in sampRemove from inputv5_pivot.
input <- input %>%
  filter(!(Ind %in% SampRemove))
input_pivot <- input_pivot %>% 
  filter(!(Ind %in% SampRemove))

# 2.8-Meta Data 1: Latitude and longitude for samples in input final (if available in metadata).---------------------------------
# Not included with publication
# Use for import into ArcGIS pro and google earth pro.
MetaData1 <- read_csv("TNWR_MetaData.csv", skip_empty_rows = TRUE) %>% 
  select(Watershed, SurveyDate, Length, Sex, River:Ind, Latitude:Longitude) %>% #added length and sex temporarily
  arrange(Watershed, Collection) %>% 
  mutate(Ind = str_replace_all(Ind, "[.]", "_")) %>% 
  mutate_at(vars(Latitude:Longitude), as.character) %>% 
  mutate_at(vars(Latitude:Longitude), ~str_pad(., 8, "right", "0")) %>% 
  mutate(Longitude = str_replace_all(Longitude, "^1", "-1")) %>% 
  semi_join(input_pivot,by="Ind") %>% 
  mutate(Collection = recode(Ind, !!!ReachCodeInd, .default = Collection)) %>% 
  rename(Aggregation = River, Date = SurveyDate, Length_mm = Length) %>% 
  select(Ind, Watershed, Aggregation, Date, Year, Length_mm, Sex, Latitude, Longitude)
# Use if want to remove samples from table that have NA for lat and long
# filter(!is.na(Latitude) | !is.na(Longitude))

# Modify aggregation code for Kanektok to reflect upper and lower drainage sections
KanektokBR <- read_csv("KanektokByReach.csv") %>% 
  mutate(Ind = str_replace_all(Ind, "[.]", "_")) %>% 
  mutate(Aggregation = paste("Kanektok_River", str_sub(SubCollection, start = -1L, end = -1L), sep = ""))
# Create named vector linking subcollection to individual for Kanektok samples
AggregationCodeInd <- set_names(KanektokBR$Aggregation,KanektokBR$Ind)
# Recode collection code with collection subcode for Kanektok samples. Kanektok samples not assigned a subcollection are filtered out
MetaData1 <- MetaData1 %>% 
  mutate(Aggregation = recode(Ind, !!!AggregationCodeInd, .default = Aggregation))

# 2.8.1-Write Meta Data 1 to directory-------------------------------
write.table(MetaData1,"MetaData1_TNWRRBT.txt",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",") 

# 2.8.2-Remove NA's from Lat and Long in Meta Data 1 and write to directory for use in ArcGIS Pro-------------------
TNWR_LatLong <- MetaData1 %>% drop_na(Longitude) %>% drop_na(Latitude)
write.table(TNWR_LatLong,"TNWR_LatLong.txt",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")

# 2.9-Meta Data 2: Genotype data--------------------------------------
# Uploaded to Dryad, not included with paper.
MetaData2 <- input %>% rename(AggregationCode = RiverCode) %>% 
    mutate(ID = group_indices(., factor(Ind, levels = unique(Ind)))) %>% 
    select(ID, everything(), -(RiverNo))

# 2.9.1-Write Meta Data 2 to directory-----------------------------------
write.table(MetaData2,"MetaData2_TNWRRBT.txt",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",") 

# 2.10-Create named vector linking Ind to ID-------------------------------------
ID_Ind <- setNames(TableS2$ID, TableS2$Ind)
saveRDS(ID_Ind, "ID_Ind.rds")
ID_Ind <- read_rds("ID_Ind.rds")
# open the modified version of Table S2 perpared for the manuscript
TableS2 <- read_csv("TableS2_TNWRRBT.txt", skip = 3)

#============================================================================================================================


# 3=PREPARE INPUT FOR FIRST LEVEL ANALYSIS USING HIERFSTAT, GENEPOP==========================================================

# 3.1-Determine sample size by Collection code and River code.--------------------------------------------------------------- 
# There are multiple collections for some rivers.
group_by(input_pivot, Collection) %>% tally()
group_by(input_pivot, RiverCode) %>% tally()

# 3.2-Simplify locus names in input_pivot from original (long) to short code (e.g., Loc1, Loc2...)-------------------------
# Define vector of original (long) locus names
LocNameLong <- names(select(input_pivot,-(Ind:RiverNo)))
# Define vector of short vector names
LocNameCode <- paste("Loc",rep(c(1:length(LocNameLong))),sep="")
# Define vector associating long locus names to codes. This object is only for reference, not used in subsequent code.
LocNameCon <- setNames(LocNameCode,LocNameLong)
saveRDS(LocNameCon, "LocNameCon.rds")
# Use rename_at function to rename multiple columns (all locus names). The tilda identifies function
input_pivot <- input_pivot %>% 
  rename_at(vars(LocNameLong), ~ LocNameCode)

# 3.3-Prepare input data for use in hierfsat packages------------------------------------------------------------------------
# Vectors of populations using Collection code (CollList) and River code (PopList) to allow for analyses by Collection and by River (for adegenet input)
CollList <- input_pivot$Collection
PopList <- input_pivot$RiverCode
# Vector of individuals (for adegenet input)
IndList <- input_pivot$Ind
# Vector of loci
LociList <- LocNameCode
# Use 'select' to remove columns without genotypes.
input4AD <- select(input_pivot, -(Ind:RiverNo)) 
# Convert to adegenet input
AD <- df2genind(input4AD,sep = '/', ncode=4, ind.names=IndList, pop=PopList, NA.char="00/00", type="codom") 
# Convert to hierfsat input
HF <- genind2hierfstat(AD) 

# 3.4-Prepare input data for use in genepop package--------------------------------------------------------------------------  
# Create a dataframe (tibble) with Individual and Population in first two columns followed by columns for each locus in two digit format ('0103').  
# The individual and population columns are both used here and one will be selected for use in GenePop depending upon analysis.
input2GP <- pivot_longer(input_pivot, c("Loc1":"Loc55"), names_to = "Loci", values_to = "Genotype") %>% 
  #remove slash "/" separating alleles
  mutate(Genotype = str_replace_all(Genotype, "/", "")) %>% 
  pivot_wider(names_from = "Loci", values_from = "Genotype") %>% 
  #two columns (Ind and RiverCode) and all loci
  select(Ind, RiverCode, Loc1:Loc55) %>% 
  rename(IndList = Ind, PopList = RiverCode)

# Use function Tib2Genepop to convert input2GP tibble to genepop input files for individual and populations.
Tib2Genepop(input2GP, "TNWRRBT_GP_Ind_Final.txt", "TNWRRBT_GP_Pop_Final.txt")

saveRDS(input_pivot, "input_pivot.rds")
saveRDS(input, "input.rds")

#============================================================================================================================


# 4-FIRST LEVEL ANALYSIS=====================================================================================================
input <- readRDS("input.rds")
input_pivot <- readRDS("input_pivot")

# 4.1-Use hierFstat to compute Ho, Hs, Fis, Fst for all loci and populations-------------------------------------------------
B <-  basic.stats(HF) 
PopHs_mean <- round(colMeans(B$Hs,na.rm=TRUE),3)
WC <- wc(HF)

# 4.2-Use Poppr to estimate the Shannon-Weiner index for each locus----------------------------------------------------------
# Used to select most informative locus of locus pair if linked
SW_Nuc <- locus_table(AD, index = 'shannon')  
# Coerce SW_Nuc to tibble
SW_Nuc <- as_tibble(rownames_to_column(as.data.frame(SW_Nuc), var = "Locus")) 
# Use "filter" in dplyr package to remove last row (says "mean" in Locus column) showing mean SW estimates for all loci
SW_Nuc <- filter(SW_Nuc, Locus != "mean") 

#4.3-Use Genepop to test HWP and GD----------------------------------------------------------------------------------------
locinfile <-  "TNWRRBT_GP.txt"
basic_info(locinfile, outputFile = "TNWRRBT_GP_Out.txt")
# Test HWP (Hardy_Weinberg proportions) for all loci and collections
test_HW(locinfile, outputFile = "TNWRRBT_GP_HW.txt", enumeration = TRUE, dememorization = 10000, batches = 100, iterations = 5000)
# Test GD (gametic disequilibrium) for all pairs of loci across collections
test_LD(locinfile, outputFile = "TNWRRBT_GP_GD.txt", dememorization = 10000, batches = 100, iterations = 5000)

# 4.4-Evaluate results of HWP tests from genepop-----------------------------------------------------------------------------

# First, test assumption of global HWP across all loci and collections
# use cumulative binomial probability distribution to derive 95% CI to the number of table-wide significant tests (p < 0.05). See Waples 2015. 
# Modify the test_HW output using python script (outside of R) to create a locus x collection matrix of p-values. 
# Read data file matrix (created using python script) of HWP p-values 24 pops. The col_types argument explicitly specifies the column data type (c=character, d=double).
inputHWP <- read_csv("TNWRRBT_GP_HW_Pvalues.txt", col_types = paste("c",strrep("d",24), sep=""), na = "-") 
# Create tibble from inputHWP that sums the number p-values < 0.05 for each locus.
LocLTalfa <- tibble('Loc' = inputHWP$Loc, 'NLocLTalfa' = rowSums(inputHWP[-1] < 0.05, na.rm = TRUE)) 
# Create tibble from inputHWP that sums the number p-values < 0.05 for each population.
PopLTalfa <- tibble('Pop' = names(inputHWP[-1]), 'NPopLTalfa' = colSums(inputHWP[-1] < 0.05, na.rm=TRUE)) 
# Table-wide number of HWP tests counts number of tests (p-values) in inputHWP (without counting unmeaningful tests denoted by NA) 
NoTestsHWP <- inputHWP %>% 
  select(-(Loc)) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  mutate(Tot = rowSums(.)) %>% 
  select(Tot)
NoTestsSigHWP <- sum(LocLTalfa$NLocLTalfa)
binomtestHWP <- dbinom(1:200, size = NoTestsHWP$Tot, prob = 0.05)
binomtestHWP <- round(binomtestHWP,4)
# cumulative binomial probability distribution to determine if the number of significant tests N (e.g., N = sum(LocLTalfa$NLocLTalfa)) is within (NS) the 0.025 < N < 0.975 interval.
binomtestcumHWP <- cumsum(binomtestHWP) 

# Second, test assumption of HWP for each collection across loci and each locus across collections

# Population-level test of HWP across loci. Genepop data file summarized using python script.
PopSigTest <-  read_csv("TNWRRBTGenepopHWP_PopSigTest.txt", col_types = paste('c',strrep('d',6), sep=''), na = '-')
PopBinomTest <- dbinom(0:length(LocLTalfa$Loc), size = length(LocLTalfa$Loc), prob = 0.05)
PopBinomTestTable <- tibble('NumSigLoci' = 0:length(LocLTalfa$Loc), 'PopBinomExp' = round(PopBinomTest * length(PopLTalfa$Pop), 2), 'PopBinomObs' = as.vector(table(factor(PopSigTest$ST, levels = 0:length(LocLTalfa$Loc)))))
PopBinomTestTableTidy <- gather(PopBinomTestTable, 'PopBinomExp', 'PopBinomObs', key = 'ObsOrExp', value = 'NumAggregations')
# Plot results for collections across loci
PopBinomTestPlot <- ggplot(PopBinomTestTableTidy) + 
  geom_bar(aes(NumSigLoci, NumAggregations, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray"), labels = c(" Expected", " Observed")) + 
  theme_bw() + 
  ggtitle("Aggregations") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,55,5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.2, 0.8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(vjust = -10, hjust = 0.5))

# Locus-level test of HWP across populations. Genepop data file summarized using python script.
LocSigTest <-  read_csv("TNWRRBTGenepopHWP_LocSigTest.txt", col_types = paste('c',strrep('d',6), sep=''), na = '-')
LocBinomTest <- dbinom(0:length(PopLTalfa$Pop), size = length(PopLTalfa$Pop), prob = 0.05)
LocBinomTestTable <- tibble('NumSigAggregations' = 0:length(PopLTalfa$Pop), 'LocBinomExp' = round(LocBinomTest * length(LocLTalfa$Loc), 2), 'LocBinomObs' = as.vector(table(factor(LocSigTest$ST, levels = 0:length(PopLTalfa$Pop)))))
LocBinomTestTableTidy <- gather(LocBinomTestTable, 'LocBinomExp', 'LocBinomObs', key = 'ObsOrExp', value = 'NumLoci')
# Plot results for loci across collections
LocBinomTestPlot <- ggplot(LocBinomTestTableTidy) + 
  geom_bar(aes(NumSigAggregations, NumLoci, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray"), labels = c(" Expected", " Observed")) + 
  theme_bw() + 
  ggtitle("Loci") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,35,5), limits = c(0,35)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,24,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.2,0.8), 
        axis.title = element_text(size=12), axis.text = element_text(size=10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(), legend.key.size = unit(0.3,"cm"),
        plot.title = element_text(vjust = -10, hjust = 0.5)) 

# 4.5-FIGURE S1: Binomial distribution of HWP results for loci and populations-----------------------------------------------
# Print Figure S1 to directory
tiff('FigS1_TNWRRBTPaper.tiff', height = 6.5, width = 8.5, units = 'in', compression = 'none', res = 300) 
ggarrange(PopBinomTestPlot, LocBinomTestPlot, ncol=1, nrow=2, align = 'v')
dev.off()

# 4.6-Evaluate results of GD tests from genepop------------------------------------------------------------------------------

# read data file matrix of GD p-values for all loci pairs for all 24 populations.  Data from genepop summarized using python script.
inputGD <- read_csv("TNWRRBTGenepopGD_Pvalues.txt", col_types = paste('c', strrep('d',24), strrep('i',2), sep=''), na='NA') 
NoTestsGD <- length(inputGD$LocXLoc)*length(PopLTalfa$Pop)
# Table-wide number of GD tests counts number of tests (p-values) in inputGD (without counting unmeaningful tests denoted by NA) 
NoTestsGD <- inputGD %>% 
  select(-c(LocXLoc,ST,BST)) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  mutate(Tot = rowSums(.)) %>% 
  select(Tot)
NoTestsSigGD <- sum(inputGD$ST)
binomtestGD <- dbinom(1:10000, size = NoTestsGD$Tot, prob = 0.05)
binomtestGD <- round(binomtestGD,4)
# cumulative binomial probability distribution to determine if the number of significant tests N (e.g., N = sum(LocLTalfa$NLocLTalfa)) is within (NS) the 0.025 < N < 0.975 interval.
binomtestcumGD <- cumsum(binomtestGD) 
arrange(select(inputGD, LocXLoc, ST), desc(ST))

GDBinomTest <- dbinom(0:length(PopLTalfa$Pop), size = length(PopLTalfa$Pop), prob = 0.05)
GDBinomTestTable <- tibble('NumSigAggregations' = 0:length(PopLTalfa$Pop), " Expected" = round(GDBinomTest * length(inputGD$LocXLoc), 2), " Observed" = as.vector(table(factor(inputGD$ST, levels = 0:length(PopLTalfa$Pop)))))
GDBinomTestTableTidy <- gather(GDBinomTestTable, " Expected", " Observed", key = 'ObsOrExp', value = 'NumLocPairs')

# plot of count of locus pairs versus number of significant collections
GDBinomTestPlot <- ggplot(GDBinomTestTableTidy) + 
  ggtitle("Locus Pairs") +
  geom_bar(aes(NumSigAggregations, NumLocPairs, fill = ObsOrExp), position = 'dodge',  stat = 'identity') + 
  scale_fill_manual(values = c("black", "gray")) +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,850)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,24,4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.3,0.8),
        axis.title = element_text(size=12), axis.text = element_text(size=10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
        legend.title = element_blank(),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=10),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title = element_text(size = 14, vjust = -10, hjust = 0.5)) +
  annotate("text", x = c(13), y = 75, label = c("OMS00154 x\nOmy_arp_630"), size = 4) +
  annotate("segment", x = 13, xend = 13, y = 40, yend = 5, colour = "black", size=1, alpha=1, arrow=arrow(length = unit(0.3, "cm")))

# 4.7-FIGURE S2: GD results showing number of locus pairs versus number of significant collections---------------------------
# Print Figure S2 to directory.
tiff('FigS2_TNWRRBTPaper.tiff', height = 5.5, width = 8.5, units = 'in', compression = 'none', res = 300) 
print(GDBinomTestPlot)
dev.off()


# 4.8-Compute Mean MAF per locus---------------------------------------------------------------------------------------------
input <- readRDS("input.rds")
input_pivot <- readRDS("input_pivot.rds")

MAF <- pivot_longer(input_pivot, -(Ind:RiverNo), names_to = "Locus", names_ptypes = "d", values_to = "Genotype") %>% 
  separate(Genotype, into=(c("A1","A2"))) %>% 
  select(Locus, A1, A2) %>% 
  pivot_longer(c("A1","A2"), names_to = "Allele", values_to = "Score") %>% 
  filter(!(Score == "00")) %>% 
  group_by(Locus, Score) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  filter(freq == min(freq))

# 4.9-TABLE S1: descriptive statistics of loci tested in first level analysis-----------------------------------------------
# Tibble of locus short and long names.
LocusName <- tibble(Locus = names(select(input_pivot, -(Ind:RiverNo))), LocusLong = pull(distinct(input,Locus)))
# Create dataframe of descriptive stats for all loci from hierFstat and add mean MAF computed above
ByLoc <-  B$perloc %>% 
  rownames_to_column(var = "Locus") %>% 
  # coerce ByLoc to tibble
  as_tibble() %>% 
  # create columns of statistics for each locus
  mutate(Ref = c(rep("a", 55)), FisWC = WC$per.loc$FIS, FstWC = WC$per.loc$FST, SWH = SW_Nuc$H) %>% 
  # move Loc17 (OMS00154) in row 18 to bottom. This locus is linked to Loc30 (Omy_arp_630) and was removed for final analysis.
  slice(c(1:16,18:55,17)) %>% 
  inner_join(MAF, by = "Locus") %>% 
  inner_join(LocusName, by = "Locus") %>% 
  rename("MAF" = "freq") %>% 
  select(LocusLong, Ref, MAF, Ho, Hs, FisWC, FstWC, SWH) %>% 
  rename("Locus" = "LocusLong") %>% 
  # add a row of column means for diversity stats 
  bind_rows(summarise_all(., ~ifelse(is.numeric(.), mean(.), "mean"))) %>% 
  # use "mutate_if" in dplyr package to round number to 3 decimals
  mutate_if(is.double, round, 3) %>% 
  mutate_at(vars(MAF:SWH), as.character) %>% 
  mutate_at(vars(MAF:SWH), ~str_pad(.,5,"right","0")) %>% 
  # add loci that were excluded from the data set because mMAF was less than 0.01
  full_join(DropLoci_lowMAF, by = "Locus") %>%  
  # add locus that was exluded from the set because the assay failed in more than 25% of individuals.
  full_join(DropLoci_missingGT, by = "Locus") %>% 
  select(Locus, MAF, Ho, Hs, FisWC, FstWC, SWH)

# write Table S1 to directory.
write_csv(ByLoc, "TableS1_TNWRRBT.txt", col_names = TRUE) 

# 4.10-TABLE S2: descriptive statistics of aggregations tested in first level analysis-----------------------------------------------
ByPopHo <- round(colMeans(B$Ho,na.rm=TRUE),3) 
ByPopHs <- round(colMeans(B$Hs,na.rm=TRUE),3) 
ByPopFis <- round(colMeans(B$Fis,na.rm=TRUE),3) 
# Create dataframe of descriptive stats for all aggregations from hierFstat
ByPop <- tibble(Ho = ByPopHo, Hs = ByPopHs, FisWC = ByPopFis) %>% 
  rbind(c(mean(ByPopHo), mean(ByPopHs), mean(ByPopFis))) %>% 
  mutate_if(is.double, round, 3) %>% 
  mutate(Acode = c(names(ByPopHo), "mean")) %>% 
  select(Acode, everything()) %>% 
  mutate_at(vars(Ho:FisWC), as.character) %>% 
  mutate_at(vars(Ho:FisWC), ~str_pad(.,5,"right","0"))

# write Table S2 to directory.
write_csv(ByPop, "TableS2_TNWRRBT.txt", col_names = TRUE) 

# 4.11-Remove Locus 17 from data set and sort by River Number.----------------------------------------------------------------
# GD test indicated Loci 17 and 29 (OMS00154 and Omy_arp_630) are linked.
# Removed locus 17 because it has the lower overall mMAF.
inputFinal <- input_pivot %>% 
  select(-(Loc17))
input <- input %>%
  filter(!(Locus == "OMS00154"))

# 4.12-Save inputFinal and input as rds files for use in next steps.---------------------------------------------------
saveRDS(inputFinal, "inputFinal.rds")
saveRDS(input, "input.rds")


#============================================================================================================================


# 5=PREPARE INPUT FOR SECOND LEVEL ANALYSIS==================================================================================

inputFinal <- readRDS("inputFinal.rds")

# 5.1-Define WatershedCodes and RiverCodes as factors ordered geographically in input----------------------------------------
WatershedFactors <- c("KW", "KN", "AR", "GN", "OS", "TG", "UN", "IG", "SN", "WD")
RiverFactors <- c("KW1", "KN1", "KN2", "AR1", "GN1", "GN2", "GN3", "GN4", "OS1", "TG1", "TG2", "UN1", "UN2", 
                  "IG1", "IG2", "SN1", "SN2", "WD1", "WD2", "WD3", "WD4", "WD5", "WD6", "WD7")
inputFinal <- inputFinal %>% 
  mutate(WatershedCode = factor(inputFinal$WatershedCode, levels = WatershedFactors)) %>% 
  mutate(RiverCode = factor(inputFinal$RiverCode, levels = RiverFactors)) 

# 5.2-TABLE 1 Summary of collections for this study--------------------------------------------------------------------------
SampSize <- group_by(inputFinal,RiverCode) %>% 
  tally() %>% 
  rename("N" = "n")
PutPop <- select(inputFinal, WatershedCode, RiverCode) %>% 
  group_by(WatershedCode, RiverCode) %>% 
  tally() %>% 
  group_by(WatershedCode) %>% 
  tally() %>% 
  rename("Groups" = "n")
TNWRCollections <- read_csv("TNWR_Collections.csv") %>% 
  mutate(WatershedCode = factor(WatershedCode, levels = WatershedFactors)) %>% 
  mutate(RiverCode = factor(RiverCode, levels = RiverFactors)) %>% 
  select(Watershed, WatershedCode, River, RiverCode, WatershedSize_km2, Year, Month) %>% 
  group_by(Watershed, WatershedCode, River, RiverCode) %>% 
# group_by(Watershed, WatershedCode, River, RiverCode) %>% 
  summarize(Year = paste(sort(unique(Year)),collapse="|"), Month = paste(Month, collapse="|")) %>% 
  left_join(SampSize, by = "RiverCode") %>% 
  #left_join(SampSize, by = "WatershedCode") %>% 
  #left_join(PutPop, by = "WatershedCode") %>% 
  ungroup() %>%
  arrange(WatershedCode, RiverCode) %>% 
  rename(Wcode = WatershedCode, Aggregation = River, Acode = RiverCode) %>% 
  mutate(Watershed = str_replace_all(Watershed, "_River", "")) %>% 
  mutate(Aggregation = str_replace_all(Aggregation, "_River", "")) %>% 
  mutate(Aggregation = str_replace_all(Aggregation, "_Creek", ""))
# rename(Code = WatershedCode, Size_km2 = WatershedSize_km2) 

# Write Table 1 to directory
write_delim(TNWRCollections, "Table1_TNWRRBT.txt", delim = "\t", quote_escape = F, col_names=T)

# 5.3-Prepare final input data for use in adegenet and hierfsat packages-----------------------------------------------------
PopList <- inputFinal$RiverCode
# Vector of individuals (for adegenet input)
IndList <- inputFinal$Ind
# Vector of loci
LociList <- names(select(inputFinal,-(Ind:RiverNo)))
# Use 'select' to remove columns without genotypes.
input4AD <- select(inputFinal, -(Ind:RiverNo)) 
# Convert to adegenet input
AD <- df2genind(input4AD,sep = '/', ncode=4, ind.names=IndList, pop=PopList, NA.char="00/00", type="codom") 
# Convert to hierfsat input
HF<- genind2hierfstat(AD) 

# 5.4-Prepare final input data for use in genepop package--------------------------------------------------------------------
# Create a dataframe with Individual and Population in first two columns followed by columns for each locus in two digit format ('0103').  
# The individual and population columns are both used here and one will be selected for use in GenePop (see subset step below) depending upon analysis.

input2GP <- pivot_longer(inputFinal, c("Loc1":"Loc55"), names_to = "Loci", values_to = "Genotype") %>% 
  mutate(Genotype = str_replace_all(Genotype, "/", "")) %>% 
  pivot_wider(names_from = "Loci", values_from = "Genotype") %>% 
  select(Ind, RiverCode, Loc1:Loc55) %>% 
  rename(IndList = Ind, PopList = RiverCode)

# Input and output for function Tib2Genepop
Tib2Genepop(input2GP, "TNWRRBT_GP_Ind_Final.txt", "TNWRRBT_GP_Pop_Final.txt")

# 5.5-Save inputFinal, AD and HF input as rds files for use in next steps.---------------------------------------------------
saveRDS(inputFinal, "inputFinal.rds")
saveRDS(AD, "AD.rds")
saveRDS(HF, "HF.rds")
saveRDS(WatershedFactors, "WatershedFactors.rds")
saveRDS(RiverFactors, "RiverFactors.rds")

#============================================================================================================================


# 6=SECOND LEVEL ANALYSIS====================================================================================================

inputFinal <- readRDS("inputFinal.rds")
WatershedFactors <- readRDS("WatershedFactors.rds")
RiverFactors <- readRDS("RiverFactors.rds")
AD <- readRDS("AD.rds")
# HF by population
HF <- readRDS("HF.rds")
# HF by watershed
HF_watershed <- mutate(HF, pop = factor(str_sub(pop, end=2), levels(inputFinal$WatershedCode)))



# 6.1-Create an input file for STRUCTURE from inputFinal---------------------------------------------------------------------
# Used STRUCTURE outside or R as an alternative to DAPC. Used DAPC results only in first draft MS. 
WatershedCon <- setNames(c(1:10), levels(inputFinal$WatershedCode))
inputFinal_str <-   mutate(inputFinal, WatershedNo = recode(WatershedCode, !!!WatershedCon)) %>% 
  select(Ind, WatershedNo, everything()) %>% 
  pivot_longer(-c(1:6), names_to = "Loci", values_to = "Genotype") %>% 
  separate(Genotype, into = c("A1", "A2")) %>% 
  pivot_longer(c("A1", "A2"), names_to = "Allele", values_to = "Genotype") %>% 
  mutate(Genotype = str_replace_all(Genotype, c("00" = "-9"))) %>% 
  pivot_wider(names_from = "Loci", values_from = "Genotype") %>% 
  select(-(Collection))
write.table(inputFinal_str,"TNWRRBT_str.txt",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ") 



# 6.2-TABLE S3: Pairwise Fst values and Edwards genetic distance for all aggregation pairs-------------------------------------------
# Compute Weir and Cockerham pairwise Fst and Edward's genetic distance for all aggregation pairs, round to three digits. 
WCpw <- pairwise.WCfst(HF) %>% 
  round(3)
DistEdwards <- edwards.dist(genind2genpop(AD)) %>% 
  as.matrix() %>% 
  round(3)
# TableS3 is matrix with pairwise Fst in lower triangle and pairwise Edward's genetic distance in upper triangle.
TableS5 <- WCpw
TableS5[upper.tri(TableS5)] <- DistEdwards[upper.tri(DistEdwards)]
# Write TableS5 to file as csv file using base R (write_csv from readr does not work for a matrix)
write.csv(TableS5, "TableS3_TNWRRBT.csv") 


IntraPW <- WCpw
IntraPW[upper.tri(IntraPW)] <- NA
IntraPW_melt <- as_tibble(melt(IntraPW)) %>% 
filter(str_sub(Var1, start = 1L, end = 2L) == str_sub(Var2, start = 1L, end = 2L), !value == "NA") 

InterPW <- WCpw
InterPW[upper.tri(InterPW)] <- NA
InterPW_melt <- as_tibble(melt(InterPW)) %>% 
  filter(str_sub(Var1, start = 1L, end = 2L) != str_sub(Var2, start = 1L, end = 2L), !value == "NA") 

DistEdwardsLong <- DistEdwards
DistEdwardsLong[upper.tri(DistEdwardsLong)] <- NA
DistEdwardLong_melt <- as_tibble(melt(DistEdwardsLong)) %>% 
  filter(!Var1 == Var2, !value == "NA") 

# 6.3-Create input file for Arlequin in GenePop file format--------------------------------------------------------------------------
# Tibble for conversion to genepop file for use in Arlequin (three digit code for alleles)
input4Arlequin <- pivot_longer(inputFinal, Loc1:Loc55, names_to = "Loci", values_to = "Genotype") %>% 
  mutate(Genotype = str_replace_all(Genotype, c("^0" = "00", "/0" = "00"))) %>% 
  pivot_wider(names_from = Loci, values_from = Genotype) %>% 
  select(Ind, RiverCode, Loc1:Loc55) %>% 
  rename(IndList = Ind, PopList = RiverCode) 

# Use function Tib2Genepop to create genepop input for use in Arlequin.
Tib2Genepop(input4Arlequin, "GP4Arlequin_Ind.txt", "GP4Arlequin_Pop.txt")

# Write list of Pops to file for appending to Arlequin input file after conversion from genepop file.
write.table(distinct(input4Arlequin, RiverCode), "TNWRRBT_Arlequin_Pops.txt", quote = F, row.names = F, col.names = F)

# 6.3.1-FIGURE 2: AMOVA results from Arelquin--------------------------------------------------------------------------------------
F_Statistic <- c("Fst", "Fct", "Fsc")
Value <- c(0.363, 0.350, 0.020)
Fig2table <- tibble(F_Statistic, Value) %>% 
  mutate(F_Statistic = factor(F_Statistic, levels = c("Fst", "Fct", "Fsc")))

Fig2Plot <- ggplot(data = Fig2table, aes(x = F_Statistic, y = Value)) +
  #Plot type  
  geom_bar(stat = "identity", color = "black", fill = "gray60") +
  geom_text(data = Fig2table, aes(label = c("**"), angle = 0), 
            show.legend = FALSE, nudge_y = 0.03, size = 4) +
  #layout
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), axis.text = element_text(color = "black", size = 8), 
        axis.line = element_line(colour = "black"), panel.border = element_blank()) +
  #overrides  
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.5,0.1), limits = c(0,0.55))

# Print FIGURE 2 to directory
ggsave(filename = "Fig2_TNWRRBTPaper.tiff", plot = Fig2Plot, device = "tiff",
       height=2.5, width = 3.5, units="in", dpi=600)



# 6.4-FIGURE 3: Neighbor-joining Phylogram of Edwards genetic distance using Adegenet, Poppr and Ape-------------------------------
# Use Adegenet to convert genind object to genpop object
ADpop <- genind2genpop(AD)
# Use Poppr to make NJ tree based on Nei's distance and showing bootstrap values for nodes
NJtreeNei <- aboot(ADpop, strata = NULL, tree = "nj", distance = edwards.dist, 
                   sample = 1000, cutoff = 95, showtree = TRUE, missing = "mean", mcutoff = 0, 
                   quiet = FALSE, root = FALSE)

# 6.4.1-Print FIGURE 3 to directory
# Define colors for branch tip labels based on watersheds
NJWatershedCols <- c(brewer.pal(8, "Dark2")[c(1,2,2,3,4,4,4,4,5,6,6,7,7,8,8)], brewer.pal(8, "Set1")[c(1,1,2,2,2,2,2,2,2)])
# Use Ape to plot unrooted tree from Poppr
tiff('Fig3_TNWRRBTpaper.tiff', height = 3.5, width = 3.5, units = 'in', compression = 'none', res = 600) #creates a high resolution TIFF file for publication
plot(NJtreeNei, type = "unrooted", show.node.label = TRUE, edge.width = 0.5, 
     tip.color = NJWatershedCols, cex = 0.6, font = 1, lab4ut = "axial", no.margin = T)
add.scale.bar(cex = 0.6)
dev.off()

# Optional code if need to break long branches
plotBreakLongEdges(NJtreeNei, n = 1, type = "unrooted", show.node.label = TRUE, edge.width = 0.5, 
                   tip.color = NJWatershedCols, cex = 0.6, font = 1, lab4ut = "axial", no.margin = T)
# Optional code for adding a legend - did not complete
legend("bottom", legend = WatershedFactors, pch=15, pt.cex=1, cex=0.75, bty='n',
       col = WatershedCols, horiz = T, yjust = -2)



# 6.5-FIGURE 4: Use DAPC to describe relationship of individuals (and membership probability) in watersheds separated by saltwater.-------
# Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
# Initial cross-validation done to narrow range of the number of PCs evaluated.
xval1 <- xvalDapc(tab(AD, NA.method = "zero"), inputFinal$WatershedCode, training.set = 0.8)
# Final cross-validation evaluates narrow range of number of PCs.
xval2 <- xvalDapc(tab(AD, NA.method = "zero"), inputFinal$WatershedCode, training.set = 0.8,
                  result = "groupMean", center = TRUE, scale = FALSE,
                  n.pca.max = 50, n.rep = 30, xval.plot = TRUE)
# Perform dapc using 20 PCs based on xval2 results. 7 discriminant functions were chosen.
dapc_Watersheds <- dapc(tab(AD, NA.method = "zero"), inputFinal$WatershedCode, n.pca = 20)

# 6.5.1-DAPC results-------------------------------------------------------------------------------------------------
# Combined plot showing dapc scatter plots for first three LD's and membership probability histograms.
# Define colors for 10 watersheds
WatershedCols <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1")[c(1,2)])
# print dapc plot to plot screen for initial evaluation
scatter(dapc_Watersheds, 1, 2, col = WatershedCols, posi.pca = "bottomleft", scree.pca = TRUE, cstar = 0)
# Tibble of LD coordinates from dapc results.
dapc_Watersheds_Ind <-  as_tibble(dapc_Watersheds$ind.coord[,1:4], rownames = "Ind")
dapc_Watersheds_Ind$Group <-  dapc_Watersheds$grp
# dapc scatter plot of LD1 vs LD2
dapc_Watersheds_Ind_Plot1 <- ggplot(dapc_Watersheds_Ind, aes(x = LD1, y = LD2, color = Group, shape = Group, size = Group)) + 
  geom_point(alpha = 0.5, size = 2) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        legend.key = element_rect(size = 0.5), legend.key.size = unit(0.01, "lines"), legend.position = "none", 
        legend.title = element_text(size = 6), legend.text = element_text(size = 5), panel.border = element_rect(fill=NA, colour = "black", size=1.25)) +
  scale_color_manual(values = WatershedCols) +
  scale_shape_manual(values = c(rep(19, 10))) +
  scale_size_manual(values = c(rep(3, 10))) +
  scale_y_continuous(breaks = seq(-9, 9, by = 3)) +
  scale_x_continuous(breaks = seq(-9, 9, by = 3)) +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  #stat_ellipse(type = "t", show.legend = F, size = 0.75) +
  annotate('label', x = dapc_Watersheds$grp.coord[,1], y = dapc_Watersheds$grp.coord[,2], 
           label = levels(inputFinal$WatershedCode), cex = 4, fill = WatershedCols, color = "white")
# dapc scatter plot LD1 vs LD3
dapc_Watersheds_Ind_Plot2 <- ggplot(dapc_Watersheds_Ind, aes(x = LD1, y = LD3, color = Group, shape = Group, size = Group)) + 
  geom_point(alpha = 0.5, size = 2) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        legend.key = element_rect(size = 0.5), legend.key.size = unit(0.01, "lines"), legend.position = "none", 
        legend.title = element_text(size = 6), legend.text = element_text(size = 5), panel.border = element_rect(fill=NA, colour = "black", size=1.25)) +
  scale_color_manual(values = WatershedCols) +
  scale_shape_manual(values = c(rep(19, 10))) +
  scale_size_manual(values = c(rep(3, 10))) +
  scale_y_continuous(breaks = seq(-9, 9, by = 3)) +
  scale_x_continuous(breaks = seq(-9, 9, by = 3)) +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  #stat_ellipse(type = "t", show.legend = F, size = 0.75) +
  annotate('label', x = dapc_Watersheds$grp.coord[,1], y = dapc_Watersheds$grp.coord[,3], 
           label = levels(inputFinal$WatershedCode), cex = 4, fill = WatershedCols, color = "white")
# dapc scatter plot of LD1 vs LD4
dapc_Watersheds_Ind_Plot3 <- ggplot(dapc_Watersheds_Ind, aes(x = LD1, y = LD4, color = Group, shape = Group, size = Group)) + 
  geom_point(alpha = 0.5, size = 2) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        legend.key = element_rect(size = 0.5), legend.key.size = unit(0.01, "lines"), legend.position = "none", 
        legend.title = element_text(size = 6), legend.text = element_text(size = 5), panel.border = element_rect(fill=NA, colour = "black", size=1.25)) +
  scale_color_manual(values = WatershedCols) +
  scale_shape_manual(values = c(rep(19, 10))) +
  scale_size_manual(values = c(rep(3, 10))) +
  scale_y_continuous(breaks = seq(-9, 9, by = 3)) +
  scale_x_continuous(breaks = seq(-9, 9, by = 3)) +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  #stat_ellipse(type = "t", show.legend = F, size = 0.75) +
  annotate("label", x = dapc_Watersheds$grp.coord[,1], y = dapc_Watersheds$grp.coord[,4], 
           label = levels(inputFinal$WatershedCode), cex = 4, fill = WatershedCols, color = "white")

# Table of posterior probability results for each individual for each watershed.
MemProb <- as_tibble(round(dapc_Watersheds$posterior, 4), rownames = "Ind") %>% 
  mutate(WatershedCode = inputFinal$WatershedCode) %>% 
  select(Ind, WatershedCode, everything()) %>% 
  pivot_longer(cols = -c("Ind", "WatershedCode"), names_to = "Group", values_to = "PostProb") %>%
  mutate(Group = factor(Group, levels = levels(WatershedCode)))
# Identify the individuals that are assigned to a different watershed.
MisMatch <- group_by(MemProb, Ind) %>% filter(PostProb == max(PostProb)) %>% filter(Group != WatershedCode)
# Identify the individuals that are assigned to a different watershed at posterior probability >= 0.90
MisMatch90 <- group_by(MemProb, Ind) %>% filter(PostProb == max(PostProb)) %>% filter(Group != WatershedCode) %>% filter(PostProb >= 0.90)
# Summarize the individuals that are admixed (< 0.90 probability in a single watershed)
AdMix90Cum <- group_by(MemProb, Ind) %>% 
  filter(PostProb == max(PostProb)) %>% 
  filter(PostProb > 0) %>% 
  group_by(WatershedCode) %>% 
  summarize(Count = n(), AdMix = sum(PostProb < 0.90)) %>% 
  mutate(AdmixProp = round(AdMix/Count, 2))
# Mean posterior membership probability for each watershed to itself.  
MeanMemProb <- filter(MemProb, WatershedCode == Group) %>% group_by(WatershedCode) %>% summarise(PostProbMean = mean(PostProb))
#Assignment to watershed based on posterior membership probability from dapc.
watershedAssignment <- select(inputFinal, Ind:WatershedCode) %>% mutate(assignedWatershed = dapc_Watersheds$assign) %>% group_by(WatershedCode, assignedWatershed) %>% tally() %>% pivot_wider(names_from = assignedWatershed, values_from = n)

# Three histograms of membership probability (watersheds 1:4, 5:8, 9:10)
# Plot1 (KW, KN, AR, GN)
dapc_WatershedsPlot1 <- ggplot(filter(MemProb, WatershedCode %in% levels(MemProb$WatershedCode)[1:4]), 
                               aes(x = Ind, y = PostProb, fill = Group, color = Group)) + 
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none",
        legend.title = element_blank(), legend.key.width = unit(1.0, "cm"), legend.key.height = unit(0.5, "cm")) +
  guides(color = guide_legend(label.position = "bottom", nrow = 1), fill = guide_legend(label.position = "bottom", nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0001)) +
  scale_fill_manual(values = WatershedCols, aesthetics = c("color", "fill")) +
  facet_grid(. ~ WatershedCode, scales = "free_x", space = "free") +
  labs(fill = "Watershed", color = "Watershed")
# Plot2 (OS, TG, UN, IG)
dapc_WatershedsPlot2 <- ggplot(filter(MemProb, WatershedCode %in% levels(MemProb$WatershedCode)[5:8]), 
                               aes(x = Ind, y = PostProb, fill = Group, color = Group)) + 
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none") +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0001)) +
  scale_fill_manual(values = WatershedCols, aesthetics = c("color", "fill")) +
  facet_grid(. ~ WatershedCode, scales = "free_x", space = "free") +
  labs(fill = "Watershed", color = "Watershed")
#Plot3 (SN, WD)
dapc_WatershedsPlot3 <- ggplot(filter(MemProb, WatershedCode %in% levels(MemProb$WatershedCode)[9:10]), 
                               aes(x = Ind, y = PostProb, fill = Group, color = Group)) + 
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none") +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0001)) +
  scale_fill_manual(values = WatershedCols, aesthetics = c("color", "fill")) +
  facet_grid(. ~ WatershedCode, scales = "free_x", space = "free") +
  labs(fill = "Watershed", color = "Watershed")

# Combine membership probabilities histograms into single plot using ggarrange, 
dapc_Watershedsplots1to3 <- ggarrange(dapc_WatershedsPlot1, dapc_WatershedsPlot2, dapc_WatershedsPlot3, ncol = 1, nrow = 4, 
                                      heights = c(6, 6, 6, 1), common.legend = TRUE, legend = "bottom")  %>% 
  annotate_figure(left = text_grob("Posterior Membership Probability", rot = 90, hjust = 0.30), bottom = text_grob("Individuals", vjust = -6))

# 6.5.2-Print FIGURE 4 to directory-------------------------------------------------------------------------------------------------
# dapc bar plot only
ggsave(filename = "Fig4_TNWRRBTPaper.tiff", plot = dapc_Watershedsplots1to3, device = "tiff",
       height = 4.5, width = 7.5, units="in", dpi=600)

# Optional FIGURE 4 combines scatter plots and histogram 
# Combined dapc scatter plots and combined membership probability histograms plots using ggarrange. 
Fig4Plot <- ggarrange(ggarrange(dapc_Watersheds_Ind_Plot1, dapc_Watersheds_Ind_Plot2, dapc_Watersheds_Ind_Plot3, ncol = 3, labels = c("a", "b", "c")), 
                      dapc_Watershedsplots1to3, nrow = 2, labels = c("", "d"), heights = c(1,2))

tiff('Fig4_TNWRRBTPaper.tiff', width = 14, height = 11, units = 'in', compression = 'none', res = 600) 
Fig2Plot
dev.off()



# 6.6-TABLE 2: Use GeneClass (outside of R) to identify probable first generation migrants---------------------------------------------------
# Used file "TNWRRBT_GP_GeneClass.txt" from step 5.4 above as input for GeneClass.
# Read GeneClass output.
GeneClassOut <- read_delim("TNWRRBT_GeneClass_Bayes_out2.csv", delim = ";", skip = 14, trim_ws = TRUE) %>% 
  mutate(Ind = inputFinal$Ind) %>% 
  select(Ind, home:`Nb. of loci`) %>% 
  rename("logHomeLogMax" = 3) %>% 
  rename("N_Loci" = 15) %>% 
  pivot_longer(cols = -c("Ind":"probability", "N_Loci"), names_to = "AssignedWatershed", values_to = "LogL") %>% 
  group_by(Ind) %>% filter(LogL == min(LogL))
# Table of first generation migrants with assignment probability <= 0.01.
FirstGenMig <- filter(GeneClassOut, home != AssignedWatershed) %>% filter(probability <= 0.01)
# Read vector relating individual code to ID number from Table S2 above
ID_Ind <- readRDS("ID_Ind.rds")
# Add membership probability results from dapc to first generation migrant results from GeneClass. 
Table3 <- filter(MemProb, Ind %in% FirstGenMig$Ind) %>% 
  filter(PostProb > 0) %>% 
  filter(!WatershedCode == Group) %>% 
  full_join(FirstGenMig, by = "Ind") %>% 
  mutate(Support = if_else(PostProb >= 0.9, "**", "*" )) %>% 
  select(Ind, WatershedCode, AssignedWatershed, probability, PostProb, Support) %>% 
  rename(To = WatershedCode, From = AssignedWatershed, MPr = PostProb, P_nm = probability) %>% 
  arrange(desc(Support))%>% 
  mutate(ID = recode(Ind,!!!ID_Ind)) %>% 
  select(ID, everything(), -(Ind)) %>% 
  mutate(MPr = round(MPr, 2)) %>% 
  filter(!(To == "SN" | To == "WD"))

# 6.6.1-Print TABLE 2 to directory.-------------------------------------------
write_csv(Table3,"Table3_TNWRRBT.txt",quote_escape = F, col_names=T, na = "--")



# 6.7-FIGURE S3: Use hierFstat to compute Hs for all aggregations and drainages--------------------------------------------------------------
# mean Hs by population
B <-  basic.stats(HF) 
PopHs_mean <- round(colMeans(B$Hs,na.rm=TRUE),3)
# Watershed stats
B_watersheds <- basic.stats(HF_watershed) 
# mean Hs by watershed
WatershedHs_mean <- round(colMeans(B_watersheds$Hs,na.rm=TRUE),3) %>% 
  tibble::enframe(name = "Watershed")

Watershed_Hs <- as.data.frame(B_watersheds$Hs) %>% rownames_to_column("Locus") %>% as_tibble() %>% 
  pivot_longer(-c("Locus"), names_to = "Watershed", values_to = "Hs") %>% 
  mutate(Watershed = factor(Watershed, levels = WatershedFactors))

Watershed_km2 <- read_csv("TNWR_Collections.csv") %>% 
  group_by(WatershedCode) %>% 
  summarize(Size_km2 = max(WatershedSize_km2)) %>% 
  mutate(WatershedCode = factor(WatershedCode, levels = WatershedFactors)) %>% 
  arrange(WatershedCode) %>% 
  rename(Watershed = WatershedCode)

Watershed_Hs_km2 <- left_join(Watershed_Hs, Watershed_km2)
saveRDS(Watershed_Hs_km2, "Watershed_Hs_km2.rds")
Watershed_Hs_km2 <- readRDS("Watershed_Hs_km2.rds")
  

Watershed_km2 <- distinct(Watershed_Hs_km2, Watershed, .keep_all = TRUE) %>% select(Watershed, Size_km2)

lm_Watershed_Hs_lnkm2 <- lm(Hs ~ Size_km2, test)
lm_Watershed_Hs_lnkm2_noSN <- lm(Hs ~ Size_km2, filter(Watershed_Hs_km2, !(Watershed == c("SN"))))

#Figure S3 plot
FigS3_Hs_km2_plot <- ggplot(data = Watershed_Hs_km2, aes(x = Size_km2, y = Hs)) +
#Data visualization  
  geom_boxplot(aes(group = Watershed)) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="black", fill="black") +
#  geom_smooth(data = filter(Watershed_Hs_km2, !Watershed == c("SN")), aes(x = Size_km2, y = Hs, color = "noSN", linetype = "noSN"), 
#              method=lm, se = FALSE) +
#  geom_smooth(aes(x = Size_km2, y = Hs, color = "All", linetype = "All"), method=lm, se = FALSE) +
#  geom_text(data = filter(Watershed_km2, !Watershed %in% c("OS", "UN", "WD")), aes(label = Watershed, y = 0.53, angle = 0), 
#            show.legend = FALSE, size = 3) +
  geom_text(data = filter(Watershed_km2, Watershed %in% c("OS", "SN", "IG", "KW")), aes(label = Watershed, y = 0.53, angle = 0), 
            show.legend = FALSE, nudge_x = -200, nudge_y = -0.03, size = 2.5) +
  geom_text(data = filter(Watershed_km2, Watershed %in% c("UN", "AR", "KN", "WD")), aes(label = Watershed, y = 0.53, angle = 0), 
            show.legend = FALSE, nudge_x = 50, size = 2.5) +
  geom_text(data = filter(Watershed_km2, Watershed %in% c("GN", "TG")), aes(label = Watershed, y = 0.53, angle = 0), 
            show.legend = FALSE, size = 2.5) +
#layout  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(vjust = -1), 
        axis.text = element_text(color = "black"), legend.position = c(0.2, 0.94), legend.background = element_blank(), 
        legend.title = element_blank(), legend.text = element_text(size = 6), 
        legend.direction = "vertical", legend.key.size = unit(0.35, "cm"), legend.key.width = unit(1.0, "cm"),
        legend.key = element_blank()) +
#overrides  
  scale_linetype_manual(name = "", values = c("All" = 1, "noSN" = 2)) +
  scale_colour_manual(name = "", values = c("All" = "black", "noSN" = "black")) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.6,0.1), limits = c(0,0.65)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,5000,1000), limits = c(0, 5250)) +
  labs(x = "Drainage Size (km2)")

#  annotate("text", x = 1000, y = 0.6, label = c("R2 = 0.003\nP = 0.27"), size = 4) 

# 6.7.1-Print FIGURE S3 to directory------------------------------------------------------------------------------------------
ggsave(filename = "FigS3_TNWRRBTPaper.tiff", plot = FigS3_Hs_km2_plot, device = "tiff",
       height=2.5, width = 3.5, units="in", dpi=600)

# 6.7.2-Write Meta Data 3 to directory------------------------------------------------------------------------------------------
LocNameCon <- readRDS("LocNameCon.rds")
LocNameConRev <- setNames(names(LocNameCon), LocNameCon)
MetaData3 <- Watershed_Hs_km2 %>%
  mutate(Locus = recode(Locus,!!!LocNameConRev))
write_csv(MetaData3, "MetaData3_TNWRRBT.csv", col_names = TRUE) 



# 6.8-Use GenePop to test for genotypic frequency heterogeniety bewtween pairs of collections.-------------------------------
locinfile <-  "TNWRRBT_GP_Pop_Final.txt"
test_diff(locinfile, genic = FALSE, pairs = TRUE, outputFile = "TNWRRBT_GP_Final_OutDiffPW.txt", dememorization = 10000, batches = 100, iterations = 5000)



# 6.9-FIGURE 5: Pairwise Fst values and G-test for intra-drainage aggregation pairs.-------------------------------------------
# Compute Weir and Cockerham pairwise Fst for all aggregation pairs, round to four digits. 
WCpw <- pairwise.WCfst(HF) %>% 
  round(4)
# PW Fst matrix from hierfstat converted to column data.
WCpw_LTri <- WCpw
WCpw_LTri[lower.tri(WCpw_LTri)] <- NA
WCpw_LTri_melt <- as_tibble(melt(WCpw_LTri)) %>% 
  filter(Var1 == Var2 | !(value == "NA")) %>% 
  mutate(value = replace(value, value < 0, 0)) %>% 
  mutate(FstCat = cut(value,breaks=c(-0.05,0.05,0.10,0.20,0.30,0.40, max(value, na.rm = T)), 
                                      labels=c("0-0.05","0.05-0.10","0.10-0.20","0.20-0.30","0.30-0.04", ">0.40"))) %>% 
  mutate(FstCat = factor(FstCat,levels = rev(levels(FstCat))))

# PW Fst intradrainage
PW_Fst_IntraDrainage <- filter(WCpw_LTri_melt, str_sub(Var1, start = 1L, end = 2L) == str_sub(Var2, start = 1L, end = 2L), !value == "NA") %>% 
  mutate(Drainage = str_sub(Var1, start = 1L, end = 2L)) %>% 
  unite(PW, c("Var1", "Var2"), sep = "/")

# PW G-test results from GenePop. The results text file was manually edited in a text editor to create the *_sub file.
PW_Gtest <- read_table("TNWRRBT_GP_Final_OutDiffPW_sub.txt", skip = 5, col_names = F) %>% 
  select(X1, X3, X5) %>% 
  rename(Pop1 = X1, Pop2 = X3, PW_Pval = X5) %>% 
  mutate(PW_Pval = recode(PW_Pval, "Highly sign." = "0.000000")) %>% 
  mutate(PW_Pval = str_remove(PW_Pval, "<")) %>% 
  mutate(PW_Pval = as.double(PW_Pval))

# PW G-test intradrainage
PW_Gtest_IntraDrainage <- filter(PW_Gtest, str_sub(Pop1, start = 1L, end = 2L) == str_sub(Pop2, start = 1L, end = 2L), !PW_Pval == "NA") %>% 
  unite(PW, c("Pop1", "Pop2"), sep = "/")

# Combine intradrainage Fst and G-test results 
PW_IntraDrainage <-  left_join(PW_Fst_IntraDrainage, PW_Gtest_IntraDrainage,by = "PW") %>% 
  rename(Fst = value, PW_gtestPval = PW_Pval, AggPair = PW) %>% 
  select(Drainage, AggPair, Fst, PW_gtestPval) %>% 
  mutate(Gtest = ifelse(PW_gtestPval <= 0.05, PW_gtestPval, "ns")) %>% 
  mutate(GtestSig = ifelse(Gtest == "ns", "ns", ifelse(Gtest <= 0.001, "***", ifelse(Gtest <= 0.01, "**", "*")))) %>% 
  mutate(BonF = c("b", rep(c(""), 4), "b", "", "b", rep(c(""), 6), "b", rep(c(""), 2), "b", rep(c(""), 6), "b", "b", rep(c(""), 3), "b", "b", "")) %>% 
  mutate(GtestBonF = paste(GtestSig, BonF, sep = "")) %>% 
  mutate_if(is.double, round, 3) %>% 
  mutate(Drainage = factor(Drainage, unique(Drainage)), AggPair = factor(AggPair, unique(AggPair))) %>% 
#  mutate(Fst = sprintf('%.3f', Fst), PW_gtestPval = sprintf('%.3f', PW_gtestPval)) %>% 
  select(Drainage, AggPair, Fst, Gtest, GtestSig, BonF, GtestBonF) %>% 
  rename(Drg = Drainage)

PW_IntraDrainageMeanFst <- group_by(PW_IntraDrainage, Drg) %>% summarize(mFst_a = mean(Fst), PWtests = n(), 
                                                      Aggs = (1+sqrt(1-(4*-2*n())))/2, GtestSig = sum(!Gtest == "ns")) %>% 
  rename(N_a = Aggs) %>% mutate(SigG_a = paste(GtestSig, "(", PWtests, ")", sep = ""))  %>% select(Drg, N_a, mFst_a, SigG_a) %>% 
  mutate(Drg = factor(Drg, levels = WatershedFactors))

saveRDS(PW_IntraDrainageMeanFst, "PW_IntraDrainageMeanFst.rds")
saveRDS(PW_IntraDrainage, "PW_IntraDrainage.rds")

# 6.9.1-Build FIGURE 5-----------------------------------------------------------------------------------------------------------

PW_IntraDrainage <- readRDS("PW_IntraDrainage.rds") %>% 
  separate(AggPair, c("P1", "P2"), sep="/") %>% 
  mutate(Drg = factor(Drg, levels = WatershedFactors)) %>% 
  mutate(Pop = "ByLoca") %>% 
  filter(!(Drg == "SN")) %>% 
  mutate(Pair = paste(Drg, paste(str_sub(P1, start = -1L, end = -1L), str_sub(P2, start = -1L, end = -1L), sep = "_"))) %>% 
  mutate(Pair = factor(Pair, unique(Pair))) %>% 
  mutate(Fst = replace(Fst, Fst == 0, 0.0005))

Fig5plot <- ggplot(PW_IntraDrainage, aes(x = Pair, y = Fst, fill = "black", color = "black")) +
  #Data visualization  
  geom_bar(stat = "identity", width = 0.50) +
  #layout  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 10), axis.text.x = element_text(size = 6, color = "black"), 
        axis.title.x = element_text(size = 10), legend.position = "none") + 
  #overrides 
  scale_x_discrete(limits = rev(levels(PW_IntraDrainage$Pair))) +
  coord_flip() +  
  ylab("PW_Fst") +
  scale_color_manual(values = c("black")) +
  scale_fill_manual(values = c("black")) +  
  geom_text(aes(label = ifelse(GtestBonF == "ns", "", GtestBonF)), vjust = 0.75, hjust = -0.25, size = 2) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0.00,0.08,0.02), limits = c(0.00, 0.075)) +
  annotate("text", x = 18, y = 0.07, label = "Gtest   \n* P<0.050\n** P<0.010\n*** P<0.001\nb  SigBonf", hjust = 1, size = 2)

# 6.9.2-Print FIGURE 5 to directory-----------------------------------------------------------------------------------------
ggsave(filename = "Fig5_TNWRRBTPaper.tiff", plot = Fig5plot, device = "tiff",
       height=3.5, width = 3.5, units="in", dpi=600)

ggsave(filename = "Fig5_TNWRRBTPowerPoint.tiff", plot = Fig5plot, device = "tiff",
       height=3.5, width = 3.5, units="in", dpi=600)

# 6.9.3-Optional TABLE for PW intra-drainage Fst and Gtest-------------------------------------------------------------------------------------------------

PW_IntraDrainageMeanFst <- readRDS("PW_IntraDrainageMeanFst.rds")
PW_IntraDrainageMeanFst_dapc <- readRDS("PW_IntraDrainageMeanFst_dapc.rds")
KrangeJoin <- readRDS("KrangeJoin.rds")

Table4 <- left_join(KrangeJoin, PW_IntraDrainageMeanFst) %>% 
  left_join(PW_IntraDrainageMeanFst_dapc) %>% 
  mutate(N_d = dapc_K) %>% 
  mutate(N_d_range = K_range) %>% 
  select(Drg, N_a, mFst_a, SigG_a, N_d, N_d_range, mFst_d, SigG_d) %>%
  mutate(N_a = replace_na(N_a, 1)) %>% 
  mutate(mFst_a = sprintf('%.3f', mFst_a), mFst_d = sprintf('%.3f', mFst_d)) 

write_csv(Table4,"Table4_TNWRRBT.txt",quote_escape = F, col_names=T, na = "--")



# 6.10-Use GeneClass outside of R to do intra-drainage individual assignment-----------------------------------------
# Modify input file for conversion to genepop file using the function Tib2Genepop
input2GP <- pivot_longer(inputFinal, c("Loc1":"Loc55"), names_to = "Loci", values_to = "Genotype") %>% 
  mutate(Genotype = str_replace_all(Genotype, "/", "")) %>% 
  pivot_wider(names_from = "Loci", values_from = "Genotype") %>% 
  select(Ind, RiverCode, Loc1:Loc55) %>% 
  rename(IndList = Ind, PopList = RiverCode)
# add a filter for selecting a subset of the data (e.g., all pops in PopList with KN)  
filter(str_detect(PopList, "KN"))

# Input and two output files for function Tib2Genepop
# Use output files as input for Geneclass (ind or pop level)
Tib2Genepop(input2GP, "TNWRRBT_GP_Ind_Final.txt", "TNWRRBT_GP_Pop_Final.txt")


# 6.11-TABLE 3: Intra-drainage individual assignment results from geneclass-------------------------------------------------------
# collating and combining assignment results from GeneClass for each drainage.
AssignmentDataDir <- "./AssignmentFiles/"
# The object "files" is a character vector of the file names in "InputFiles"
files = dir(path = AssignmentDataDir, pattern="*.csv")
# Define input
# Read the output data files from the "AssignmentFiles" folder in the TogiakRainbowTrout R-project and create tibble "AssignmentResults"
AssignmentResults <-  paste(AssignmentDataDir,files,sep="") %>% 
  # Use the map function in purrr to read the output data files from geneclass as lists. The first 12 rows are skipped.
  map(read_delim, delim = ";", skip = 12) %>%
  map(select, 1:3, 5:6) %>% 
  map(rename, AssignedSample = 'Assigned sample', AssignedFirst = '1', Score1 = '%', AssignedSecond = '2', Score2 = '%_1') %>% 
  map(mutate, AssignedSample = str_replace_all(AssignedSample, "/", "")) %>%   
  # Use the reduce function in purrr to combine the lists (data files) by row using rbind.
  reduce(rbind) %>% 
  mutate(Drainage = str_sub(AssignedSample, start = 1L, end = 2L)) %>% 
  select(Drainage, everything()) %>% 
  mutate(AssignedSample = factor(AssignedSample, levels = RiverFactors), Drainage = factor(Drainage, levels = WatershedFactors)) %>% 
  # Add column that labels individuals as unassigned if likelihood score is < 90, otherwise label as assigned_home or assigned_other.
  mutate(Assigned = ifelse(Score1 < 90, "unassigned", ifelse(AssignedFirst == AssignedSample, "assigned_home", "assigned_other")))

# Summarize assignment results for each aggregation by watershed as assigned_home, assigned_other, unassigned.
AssignmentSum <- AssignmentResults %>% 
  group_by(Drainage, AssignedSample, Assigned) %>% 
  tally() %>% 
  mutate(freq = n/sum(n)) %>%
  mutate(freq = round(freq, 2)) %>% 
  filter(!(Drainage == "SN"))

Table3 <- ungroup(AssignmentSum) %>% select(Drainage, AssignedSample, Assigned, freq) %>% 
  pivot_wider(names_from = Assigned, values_from = freq) %>% 
  rename(Ws = Drainage, Ag = AssignedSample, pHome = assigned_home, pOther = assigned_other, pUn = unassigned) %>% 
  select(Ag, everything())

# 6.11.1-Write TABLE 3 to directory-----------------------------------------------------------------------------------------
# Replace NA with zero
Table3[is.na(Table3)] <- 0
# Write table
write_delim(Table3, "Table3_TNWRRBT.txt", delim = "\t", quote_escape = F, col_names=T)


# Optional plot assignment results by drainage.
AssignmentSum_Plot1 <- ggplot(filter(AssignmentSum, Drainage %in% c("KN", "GN", "TG", "UN")), aes(x = AssignedSample, y = freq, fill = Assigned)) +
  geom_col(color = "black", position = position_fill(reverse = T)) +
  theme_bw() + 
  facet_grid(~Drainage, scales = "free_x", space = "free") +
  scale_fill_manual(values = c("black", "gray", "white")) +   
  theme(legend.position = "none", legend.title = element_blank(), axis.title = element_blank(), 
        axis.text = element_text(size = 6, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.x =element_text(size = 6, margin = margin(.1, 0, .1, 0, "cm")))

AssignmentSum_Plot2 <- ggplot(filter(AssignmentSum, Drainage %in% c("IG", "WD")), aes(x = AssignedSample, y = freq, fill = Assigned)) +
  geom_col(color = "black", position = position_fill(reverse = T)) +
  theme_bw() + 
  facet_grid(~Drainage, scales = "free_x", space = "free") +
  scale_fill_manual(values = c("black", "gray", "white")) +   
  theme(legend.position = c(0.5, -0.6), legend.title = element_blank(), legend.text = element_text(size = 6), legend.direction = "horizontal", 
        legend.key.size = unit(0.05, "cm"), legend.key.width = unit(0.15, "cm"), axis.title = element_blank(), 
        axis.text = element_text(size = 6, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.x =element_text(size = 6, margin = margin(.1, 0, .1, 0, "cm")))

Fig6Plot <- ggarrange(AssignmentSum_Plot1, AssignmentSum_Plot2, ncol = 1, nrow = 3, 
                      heights = c(6, 6, 2), common.legend = FALSE)  %>% 
  annotate_figure(left = text_grob("Proportion", rot = 90, size = 10, hjust = 0.20), bottom = text_grob("Aggregation Code", vjust = -4.0, size = 10))

# Optional print assignement plot to directory
ggsave(filename = "Fig6_TNWRRBTPaper.tiff", plot = Fig6Plot, device = "tiff",
       height=3.5, width = 3.5, units="in", dpi=600)



# 6.12-FIGURE 6: Results of DAPC analysis to estimate number of aggregations in each drainage.---------

# 6.12.1-Infer the number of aggregations in each drainage using dapc and find.clusters-------------------------------------------------- 
#KW
inputKW <-  filter(inputFinal, WatershedCode %in% c("KW")) %>%  
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputKW, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputKW$Ind, pop=inputKW$RiverCode, NA.char="00/00", type="codom") 
KWGrp <- find.clusters(AD, max.n.clust = 5)
inputKW_dapc <- mutate(inputKW, dapcCluster = KWGrp$grp) %>% 
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything())

#KN
inputKN <-  filter(inputFinal, WatershedCode %in% c("KN")) %>%  
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputKN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputKN$Ind, pop=inputKN$RiverCode, NA.char="00/00", type="codom") 
KNGrp <- find.clusters(AD, max.n.clust = 5)
inputKN_dapc <- mutate(inputKN, dapcCluster = KNGrp$grp) %>% 
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything())

#AR
inputAR <-  filter(inputFinal, WatershedCode %in% c("AR")) %>%  
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputAR, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputAR$Ind, pop=inputAR$RiverCode, NA.char="00/00", type="codom") 
ARGrp <- find.clusters(AD, max.n.clust = 5)
inputAR_dapc <- mutate(inputAR, dapcCluster = ARGrp$grp) %>% 
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything())

# GN
inputGN <- filter(inputFinal, WatershedCode %in% c("GN")) %>% 
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputGN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputGN$Ind, pop=inputGN$RiverCode, NA.char="00/00", type="codom")
GNGrp <- find.clusters(AD, max.n.clust = 5)
inputGN_dapc <- mutate(inputGN, dapcCluster = GNGrp$grp) %>% 
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything())

#OS
inputOS <-  filter(inputFinal, WatershedCode %in% c("OS")) %>% 
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputOS, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputOS$Ind, pop=inputOS$RiverCode, NA.char="00/00", type="codom")
OSGrp <- find.clusters(AD, max.n.clust = 5)
inputOS_dapc <- mutate(inputOS, dapcCluster = OSGrp$grp) %>% 
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything())

# TG
inputTG <-  filter(inputFinal, WatershedCode %in% c("TG")) %>% 
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputTG, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputTG$Ind, pop=inputTG$RiverCode, NA.char="00/00", type="codom")
TGGrp <- find.clusters(AD, max.n.clust = 5)
inputTG_dapc <- mutate(inputTG, dapcCluster = TGGrp$grp) %>%
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything()) 

# UN
inputUN <-  filter(inputFinal, WatershedCode %in% c("UN")) %>% 
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputUN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputUN$Ind, pop=inputUN$RiverCode, NA.char="00/00", type="codom")
UNGrp <- find.clusters(AD, max.n.clust = 5)
inputUN_dapc <- mutate(inputUN, dapcCluster = UNGrp$grp) %>%
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything())  

#IG
inputIG <-  filter(inputFinal, WatershedCode %in% c("IG")) %>% 
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputIG, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputIG$Ind, pop=inputIG$RiverCode, NA.char="00/00", type="codom")
IGGrp <- find.clusters(AD, max.n.clust = 5)
inputIG_dapc <- mutate(inputIG, dapcCluster = IGGrp$grp) %>%
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything())  

#SN
inputSN <-  filter(inputFinal, WatershedCode %in% c("SN")) %>% 
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputSN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputSN$Ind, pop=inputSN$RiverCode, NA.char="00/00", type="codom")
SNGrp <- find.clusters(AD, max.n.clust = 5)
inputSN_dapc <- mutate(inputSN, dapcCluster = SNGrp$grp) %>%
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything()) 

#WD
inputWD <-  filter(inputFinal, WatershedCode %in% c("WD")) %>% 
  filter(!(Ind %in% FirstGenMig$Ind))
AD <- df2genind(select(inputWD, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputWD$Ind, pop=inputWD$RiverCode, NA.char="00/00", type="codom")
WDGrp <- find.clusters(AD, max.n.clust = 10)
inputWD_dapc <- mutate(inputWD, dapcCluster = WDGrp$grp) %>%
  mutate(dapcCluster = paste(WatershedCode, dapcCluster, "*", sep = "")) %>% 
  select(Ind:RiverNo, dapcCluster, everything()) 

#Summarize the dapc find.cluster results to show the inferred number of clusters and range (inferred cluster +/- 6 BIC)
newClustersList <- list(KWGrp$Kstat, KNGrp$Kstat, ARGrp$Kstat, GNGrp$Kstat, OSGrp$Kstat, TGGrp$Kstat, UNGrp$Kstat, IGGrp$Kstat, SNGrp$Kstat, WDGrp$Kstat)
names(newClustersList) <- WatershedFactors

saveRDS(newClustersList, "newClustersList.rds")
newClustersList <- readRDS("newClustersList.rds")

Kselect <- tibble::enframe(sapply(newClustersList, which.min, USE.NAMES = TRUE), name = "Drg", value = "dapc_K") %>% 
  mutate(Drg = WatershedFactors)

Krange <- sapply(newClustersList, function(x) (x-min(x) < 6), USE.NAMES = TRUE)

KrangeNS <- sapply(Krange, function(x) which(x, arr.ind = TRUE), USE.NAMES = TRUE)

KrangeNSmin <- tibble::enframe(sapply(KrangeNS, min, USE.NAMES = TRUE), name = "Drg", value = "min")
KrangeNSmax <- tibble::enframe(sapply(KrangeNS, max, USE.NAMES = TRUE), name = "Drg", value = "max")

KrangeJoin <- left_join(Kselect, KrangeNSmin) %>% left_join(KrangeNSmax) %>% unite("K_range", min:max, sep = "-") %>% 
  mutate(Drg = factor(Drg, levels = WatershedFactors)) 

saveRDS(KrangeJoin, "KrangeJoin.rds")

# 6.12.2-Add dapcClusters to inputFinal and compute PW Fst and Gtest and do HWP test for intra-drainage dapcClusters--------------------

inputFinal_dapcClusters <- bind_rows(inputKW_dapc, inputKN_dapc, inputAR_dapc, inputGN_dapc, inputOS_dapc, 
                                     inputTG_dapc, inputUN_dapc, inputIG_dapc, inputSN_dapc, inputWD_dapc)
dapcFactors <- c("KW1*", "KN1*", "KN2*", "AR1*", "GN1*", "GN2*", "GN3*", "OS1*", "OS2*", "TG1*", "TG2*", "TG3*", "UN1*", "UN2*", 
                 "UN3*", "IG1*", "IG2*", "SN1*", "SN2*", "SN3*", "SN4*", "WD1*", "WD2*", "WD3*", "WD4*", "WD5*")
inputFinal_dapcClusters <- inputFinal_dapcClusters %>% 
  arrange(WatershedCode, dapcCluster) %>% 
  mutate(dapcCluster = factor(inputFinal_dapcClusters$dapcCluster, levels = dapcFactors))

saveRDS(inputFinal_dapcClusters, "inputFinal_dapcClusters.rds")
inputFinal_dapcClusters <- readRDS("inputFinal_dapcClusters.rds")

# Compute Fst and do G-test for inferred clusters
PopList <- inputFinal_dapcClusters$dapcCluster
# Vector of individuals (for adegenet input)
IndList <- inputFinal_dapcClusters$Ind
# Vector of loci
LociList <- names(select(inputFinal_dapcClusters,-(Ind:dapcCluster)))
# Use 'select' to remove columns without genotypes.
input4AD <- select(inputFinal_dapcClusters, -(Ind:dapcCluster)) 
# Convert to adegenet input
AD <- df2genind(input4AD,sep = '/', ncode=4, ind.names=IndList, pop=PopList, NA.char="00/00", type="codom") 
# Convert to hierfsat input
HF<- genind2hierfstat(AD) 

input2GP <- pivot_longer(inputFinal_dapcClusters, c("Loc1":"Loc55"), names_to = "Loci", values_to = "Genotype") %>% 
  mutate(Genotype = str_replace_all(Genotype, "/", "")) %>% 
  pivot_wider(names_from = "Loci", values_from = "Genotype") %>% 
  select(Ind, dapcCluster, Loc1:Loc55) %>% 
  rename(IndList = Ind, PopList = dapcCluster) %>% 
  arrange(PopList)

# Input and two output files for function Tib2Genepop
Tib2Genepop(input2GP, "TNWRRBT_GPdapc_Ind_Final.txt", "TNWRRBT_GPdapc_Pop_Final.txt")

# Use GenePop to test for genotypic frequency heterogeniety between pairs of clusters
locinfile <-  "TNWRRBT_GPdapc_Pop_Final.txt"
test_diff(locinfile, genic = FALSE, pairs = TRUE, outputFile = "TNWRRBT_GPdapc_Final_OutDiffPW.txt", dememorization = 10000, batches = 100, iterations = 5000)
# Test HWP (Hardy_Weinberg proportions) for all loci and clusters
test_HW(locinfile, outputFile = "TNWRRBT_GPdapc_HW.txt", enumeration = TRUE, dememorization = 10000, batches = 100, iterations = 5000)

# Use hierFstat to compute Ho, Hs, Fis, Fst for all loci and clusters
B <-  basic.stats(HF) 
PopFis_mean <- round(colMeans(B$Fis,na.rm=TRUE),3)
WC <- wc(HF)

# Compute PW Fst
WCpw <- pairwise.WCfst(HF) %>% 
  round(4)
# Convert PW Fst into table
WCpw_LTri <- WCpw
WCpw_LTri[lower.tri(WCpw_LTri)] <- NA
WCpw_LTri_melt <- as_tibble(melt(WCpw_LTri)) %>% 
  filter(Var1 == Var2 | !(value == "NA")) %>% 
  mutate(value = replace(value, value < 0, 0)) %>% 
  mutate(FstCat = cut(value,breaks=c(-0.05,0.05,0.10,0.20,0.30,0.40, max(value, na.rm = T)), 
                      labels=c("0-0.05","0.05-0.10","0.10-0.20","0.20-0.30","0.30-0.04", ">0.40"))) %>% 
  mutate(FstCat = factor(FstCat,levels = rev(levels(FstCat))))

# PW Fst intradrainage
PW_Fst_IntraDrainage_dapc <- filter(WCpw_LTri_melt, str_sub(Var1, start = 1L, end = 2L) == str_sub(Var2, start = 1L, end = 2L), !value == "NA") %>% 
  mutate(Drainage = str_sub(Var1, start = 1L, end = 2L)) %>% 
  unite(PW, c("Var1", "Var2"), sep = "/")

# PW G-test results from GenePop. The results text file was manually edited in a text editor to create the *_sub file.
PW_Gtest_dapc <- read_table("TNWRRBT_GPdapc_Final_OutDiffPW_sub.txt", skip = 5, col_names = F) %>% 
  select(X1, X3, X5) %>% 
  rename(Pop1 = X1, Pop2 = X3, PW_Pval = X5) %>% 
  mutate(PW_Pval = recode(PW_Pval, "Highly sign." = "0.000000")) %>% 
  mutate(PW_Pval = str_remove(PW_Pval, "<")) %>% 
  mutate(PW_Pval = as.double(PW_Pval))

# PW G-test intradrainage
PW_Gtest_IntraDrainage_dapc <- filter(PW_Gtest_dapc, str_sub(Pop1, start = 1L, end = 2L) == str_sub(Pop2, start = 1L, end = 2L), !PW_Pval == "NA") %>% 
  unite(PW, c("Pop1", "Pop2"), sep = "/")

# Combine intradrainage Fst and G-test results 
PW_IntraDrainage_dapc <-  left_join(PW_Fst_IntraDrainage_dapc, PW_Gtest_IntraDrainage_dapc,by = "PW") %>% 
  rename(Fst = value, PW_gtestPval = PW_Pval, AggPair = PW) %>% 
  select(Drainage, AggPair, Fst, PW_gtestPval) %>% 
  mutate(Gtest = ifelse(PW_gtestPval <= 0.001, "***", ifelse(PW_gtestPval <= 0.01, "**", ifelse(PW_gtestPval <= 0.05, "*", "ns")))) %>% 
  mutate(Gtest = ifelse(PW_gtestPval <= 0.05, "P<0.05", "ns")) %>% 
  #  mutate_if(is.double, round, 3) %>% 
  mutate(Drainage = factor(Drainage, unique(Drainage)), AggPair = factor(AggPair, unique(AggPair))) %>% 
  #  mutate(Fst = sprintf('%.3f', Fst), PW_gtestPval = sprintf('%.3f', PW_gtestPval)) %>% 
  select(Drainage, AggPair, Fst, Gtest) %>% 
  rename(Drg = Drainage)

PW_IntraDrainageMeanFst_dapc <- group_by(PW_IntraDrainage_dapc, Drg) %>% summarize(mFst_d = mean(Fst), PWtests = n(), 
                                                                                   Aggs = (1+sqrt(1-(4*-2*n())))/2, GtestSig = sum(!Gtest == "ns")) %>% 
  rename(N_d = Aggs) %>% mutate(SigG_d = paste(GtestSig, "(", PWtests, ")", sep = ""))  %>% select(Drg, N_d, mFst_d, SigG_d) %>% 
  mutate(Drg = factor(Drg, levels = WatershedFactors))

saveRDS(PW_IntraDrainageMeanFst_dapc, "PW_IntraDrainageMeanFst_dapc.rds")
saveRDS(PW_IntraDrainage_dapc, "PW_IntraDrainage_dapc.rds")

# 6.12.3-Build FIGURE 6--------------------

PW_IntraDrainage_dapc <- readRDS("PW_IntraDrainage_dapc.rds") %>% 
  separate(AggPair, c("P1", "P2"), sep="/") %>% 
  mutate(Drg = factor(Drg, levels = WatershedFactors)) %>% 
  mutate(Pop = "DAPC") %>% 
  filter(!(Drg == "SN")) %>% 
  add_row(Drg = "KW", Pop = "DAPC") %>% 
  add_row(Drg = "AR", Pop = "DAPC")

PW_IntraDrainage_All <- bind_rows(PW_IntraDrainage, PW_IntraDrainage_dapc) %>% 
  mutate(Pop = str_replace_all(Pop, c("ByLoca" = "Aggregation"))) %>% 
  mutate(Pop = factor(Pop, levels = c("DAPC", "Aggregation")))
  

Drainage <- c(rep(c("KW", "KN", "AR", "GN", "OS", "TG", "UN", "IG", "WD"), 2))
N_Pops <- c(1,2,1,4,1,2,2,2,7,1,2,1,3,2,3,3,2,5)
ByLocaByDAPC <- tibble(Drainage, N_Pops) %>% 
  mutate(Drainage = factor(Drainage, levels = c("KW", "KN", "AR", "GN", "OS", "TG", "UN", "IG", "SN", "WD"))) %>% 
  mutate(Pop = c(rep(c("Aggregation"), 9), c(rep(c("DAPC"), 9)))) %>% 
  mutate(Pop = factor(Pop, levels = c("DAPC", "Aggregation"))) %>% 
  mutate(Low = c(rep(c(NA), 9), 1,1,1,1,1,1,1,1,2)) %>% 
  mutate(High = c(rep(c(NA), 9), 5,5,5,5,5,5,5,5,9))

Fig6Plot1 <- ggplot(data = ByLocaByDAPC, aes(x = Drainage, y = N_Pops, color = Pop, fill = Pop)) +
  #Plot type  
  geom_point(size = 2, shape = 22, position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), aes(ymin = Low, ymax = High), width = 0.2, size = 0.2) +
  #layout
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 10), axis.text.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 6), axis.line = element_line(colour = "black"), 
        legend.position = c(0.3, 0.9), legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.01, "lines"), legend.title = element_blank(), 
        legend.text = element_text(size = 5), legend.background = element_blank(), plot.margin = unit(c(0,0.25,0,0.5), "cm")) +
  #overrides  
  scale_color_manual(values = c("Aggregation" = "black", "DAPC" = "black")) +
  scale_fill_manual(values = c("Aggregation" = "white", "DAPC" = "black")) +
  scale_shape_manual(values = c("Aggregation" = 24, "DAPC" = 24)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,10,1), limits = c(0,9.5)) +
  guides(linetype=FALSE,color=FALSE) +
  annotate("text", x = 1, y = 8.5, label = "A", size = 3)

Fig6plot2 <- ggplot(PW_IntraDrainage_All, aes(x = Drg, y = Fst, fill = Pop, color = Pop)) +
  #Data visualization  
  geom_point(size = 2, shape = 21, position = position_dodge(width = 0.4)) +
  #layout  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 10), axis.text.x = element_text(size = 6, color = "black"), 
        axis.title.x = element_blank(), legend.key = element_rect(size = 0.5), axis.line = element_line(colour = "black"), 
        legend.key.size = unit(0.01, "lines"), legend.position = c(0.3, 0.9), legend.title = element_blank(), 
        legend.text = element_text(size = 5), legend.background = element_blank()) +
  #overrides 
  ylab("PW_Fst") +
  scale_color_manual(values = c("Aggregation" = "black", "DAPC" = "black")) +
  scale_fill_manual(values = c("Aggregation" = "white", "DAPC" = "black")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(-0.02,0.11,0.02), limits = c(-0.01, 0.12)) +
  guides(linetype=FALSE,color=FALSE) +
  annotate("text", x = 1, y = 0.105, label = "B", size = 3)

Fig6plot <- ggarrange(Fig6Plot1, Fig6plot2, ncol = 1, nrow = 3, 
                      heights = c(6, 6.5, 0.5), common.legend = FALSE) %>% 
  annotate_figure(bottom = text_grob("Drainage", vjust = -1.0, size = 10))

# 6.12.4-Print FIGURE 6 to directory-----------------------------------------------------------------------------------------
ggsave(filename = "Fig6_TNWRRBTPaper.tiff", plot = Fig6plot, device = "tiff",
       height=3.5, width = 3.5, units="in", dpi=600)



# 6.13-Surplus code---------------------------------------------------------------------------------------------------------------

# This code works but was not used for the paper.

#dapc plots-------
dapc2_MemProbPlot2 <- ggplot(filter(MemProbWooSna, RiverCode %in% RiverFactors[3:9]), aes(x = Ind, y = PostProb, fill = Group, color = Group)) + 
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.001)) +
  scale_fill_brewer(palette = "Paired", aesthetics = c("color", "fill")) +
  #geom_text(data = filter(test2, home %in% WatershedFactors[1:4]), 
  #aes(x = Ind, label = "*"), nudge_y = 0.01, vjust = 0.75, angle = 90, color = "black", size = 7) +
  facet_grid(WatershedCode ~ RiverCode, scales = "free_x", space = "free", switch = "y") +
  labs(y = "Membership Probability", x = "Individuals", fill = "River", color = "River")


dapc2_plots1to2 <- ggarrange(dapc2_MemProbPlot1, dapc2_MemProbPlot2, ncol = 2, nrow = 1, 
                             widths = c(1,3), common.legend = TRUE, legend = "right")  %>% 
  annotate_figure(left = text_grob("Membership Probability", rot = 90), bottom = text_grob("Individuals"))


grp <- find.clusters(AD, max.n.clust = 4)


dapc1DF = as.data.frame(dapc1$ind.coord[,1:2])
dapc1DF$Group = dapc1$grp
dapc1GroupTib <- as_tibble(dapc1$grp.coord[,1:2], rownames = "LakeGroup")
myCol = c('firebrick1', 'Firebrick2', 'firebrick3', 'green', 'green2', 'green4', 'darkgreen', 'dodgerblue', 'deepskyblue', 'darksalmon')
GrpLab = paste(c(1:10), '-', unique(dapc1DF$Group), sep = '')

dapc1_plot <- ggplot(dapc1DF, aes(x = LD1, y = LD2, color = Group, shape = Group)) + 
  geom_point() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.01, "lines"), legend.position = c(0.95, 0.70), legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5)) +
  scale_color_manual(values = myCol, labels = GrpLab) +
  scale_shape_manual(values = c(rep(1,3), rep(2,4), rep(5,2), 2), labels = GrpLab) +
  scale_y_continuous(limits = c(-5,7)) +
  scale_x_continuous(limits = c(-6,6)) +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  geom_hline(yintercept = 0, linetype = 'longdash') +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  #annotate('label', x = dapc1$grp.coord[,1], y = dapc1$grp.coord[,2], label = c(1:9), cex = 2, fill = myCol, color = "black") +
  #stat_ellipse(type = "t", show.legend = F, size = 0.75) +
  annotate('label', x = dapc1$grp.coord[,1], y = dapc1$grp.coord[,2], label = c(1:10), cex = 4, fill = myCol, color = "white")

loadingplot(dapc1$var.contr, axis = 2, thres = 0.04, lab.jitter = 1)

dapc1Load1 <- as.tibble(dapc1$var.contr, rownames = "LocAllele")

#Figure 4 for manuscript

dapc1_plot_PC1load <- ggplot(dapc1Load1, aes(x = LocAllele, y = LD1)) +
  geom_bar(stat = "identity", width = 0.25, color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.title.x = element_text(size = 8), axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8)) +
  geom_text_repel(data = filter(dapc1Load1, LD1 > 0.10), 
                  aes(LocAllele, LD1, label = LocAllele), size = 2) +
  scale_y_continuous(name = "PC1 Loadings", limits = c(0,0.15), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(limits = dapc1Load1$LocAllele)

dapc1_plot_PC2load <- ggplot(dapc1Load1, aes(x = LocAllele, y = LD2)) +
  geom_bar(stat = "identity", width = 0.25, color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.title.x = element_text(size = 8), axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8)) +
  geom_text_repel(data = filter(dapc1Load1, LD2 > 0.10), 
                  aes(LocAllele, LD2, label = LocAllele), size = 2) +
  scale_y_continuous(name = "PC2 Loadings", limits = c(0,0.15)) +
  scale_x_discrete(limits = dapc1Load1$LocAllele)

tiff('Figure4_FrazerPaper.tiff', height = 4.6, width = 6.8, units = 'in', compression = 'none', res = 600) #creates a high resolution TIFF file for publication
ggarrange(dapc1_plot, dapc1_plot_PC1load, dapc1_plot_PC2load, ncol = 1, nrow = 3, heights = c(2,1,1))
dev.off()


# Used to make heat map of pairwise FSt values------
WCpw_IntraDrainage <- filter(WCpw_LTri_melt, str_sub(Var1, start = 1L, end = 2L) == str_sub(Var2, start = 1L, end = 2L), !value == "NA") %>% 
  mutate(Drainage = str_sub(Var1, start = 1L, end = 2L)) %>% 
  unite(PW, c("Var1", "Var2"), sep = "_") %>% 
  select(Drainage, PW, value)

# Heatmap plot------
Fig2Plot <- ggplot(WCpw_LTri_melt, aes(x=Var1, y=Var2, fill=FstCat)) + 
  geom_tile(color = "black") + 
  theme_bw() + 
  coord_fixed() +
  #  scale_fill_gradientn(colors = rev(grDevices::heat.colors(4)), values = c(0, 0.05, 0.15, 0.30, 0.50)) +
  #  scale_fill_continuous(low = "white", high = "red3", na.value = "white") +
  scale_fill_manual(breaks = levels(WCpw_LTri_melt$FstCat), values = c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4"), na.value = "white") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 4, hjust = 0), 
        axis.text.y = element_text(angle = 0, vjust = 0.5, size = 4, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.25), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.25), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        legend.position = c(0.9, 0.8), 
        legend.text = element_text(size = 4, hjust = 0.5), 
        legend.title = element_text(size = 6),  
        legend.key.size = unit(0.2, "cm"),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(fill=guide_legend(title="Fst"))+
  scale_x_discrete(position = 'top') +
  scale_y_discrete(limits = rev(levels(inputFinal$RiverCode))) +
  annotate("text", x = 1:24, y = 24:1, label = "-", size = 2) +
  annotate("text", x = PW_Gtest_NS$Pop1, y = PW_Gtest_NS$Pop2, label = "NS", size = 1)
# Print plot to directory
tiff('Fig2_TNWRRBTPaper.tiff', height = 3.6, width = 3.6, units = 'in', compression = 'none', res = 600) #creates a high resolution TIFF file for publication
Fig2Plot
dev.off()

# Plot of intradraingage pairwise Fst values--------------
ggplot(WCpw_IntraDrainage, aes(x = AggPair, y = Fst)) +
  theme_bw() +
  #coord_flip() +
  #scale_x_discrete(limits = rev(levels(WCpw_IntraDrainage$PW))) +
  #geom_jitter(width = 0.1, size = 4) +
  geom_point(size = 4) +
  facet_grid(~Gtest, scale = "free_x", space = "free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12), axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), axis.title.x = element_text(size = 12), legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.01, "lines"), legend.position = "none", legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5)) +
  ylab("Fst") +
  xlab("Aggregation Pairs") +
  scale_color_manual(values = c("black")) +
  scale_shape_manual(values = 19)

# subset of PW G-tests that were NS (>0.05)
PW_Gtest_NS <- PW_Gtest %>%
  filter(PW_Pval >= 0.05)

# For insetting one plot into another plot.  Used for Kodiak sockeye. Leave here for possible future use.--------
# Zoom in GD binomial plot to show locus pairs with GD in many collections 
 GDBinomTestPlotZoom <- GDBinomTestPlot + coord_cartesian(ylim = c(0,40)) + 
   ggtitle("Locus Pairs") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.3,0.8),
         axis.title = element_text(size=14), axis.text = element_text(size=10), plot.margin=unit(c(0.5,10,0.5,1.2),"cm"),
         legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=10), 
         plot.title = element_text(size = 14, vjust = -10, hjust = 0.5)) +
   annotate("text", x = c(24,30), y = 3, label = c("One_MHC2_190 x\nOne_MHC2_251", "One_Tf_ex10_750 x\nOne_Tf_ex3_182"), size = 4)
# Inset overall plot within zoom plot
 GDPlotInsetZoom <- ggdraw() +
   draw_plot(GDBinomTestPlotZoom) +
   draw_plot(GDBinomTestPlot, x = 0.45, y = 0.45, width = 0.5, height = 0.4)

 
 # Use DAPC in adegenet to estimate membership probability of individuals in collections within select watersheds.-------
 # KN
 inputKN <- filter(inputFinal, WatershedCode %in% c("KN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputKN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputKN$Ind, pop=inputKN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputKan$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 5:10).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputKan$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 5:10, n.rep = 1000, xval.plot = TRUE)
 KN_dapc <- dapc(AD, inputKN$RiverCode, n.pca = 8)
 KN_MemProb <- as_tibble(round(KN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputKN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputKN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 # GN
 inputGN <- filter(inputFinal, WatershedCode %in% c("GN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputGN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputGN$Ind, pop=inputGN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputGN$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 10:30).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputGN$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 10:30, n.rep = 1000, xval.plot = TRUE)
 GN_dapc <- dapc(AD, inputGN$RiverCode, n.pca = 16)
 GN_MemProb <- as_tibble(round(GN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputGN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputGN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 # TG
 inputTG <-  filter(inputFinal, WatershedCode %in% c("TG")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputTG, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputTG$Ind, pop=inputTG$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputTG$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 2:25).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputTG$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 2:25, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 2 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 TG_dapc <- dapc(AD, inputTG$RiverCode, n.pca = 2)
 TG_MemProb <- as_tibble(round(TG_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputTG$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputTG$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 # UN
 inputUN <-  filter(inputFinal, WatershedCode %in% c("UN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputUN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputUN$Ind, pop=inputUN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputUN$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 2:30).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputUN$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 2:30, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 2 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 UN_dapc <- dapc(AD, inputUN$RiverCode, n.pca = 2)
 UN_MemProb <- as_tibble(round(UN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputUN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputUN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 # IG
 inputIG <-  filter(inputFinal, WatershedCode %in% c("IG")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputIG, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputIG$Ind, pop=inputIG$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputIG$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 5:25).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputIG$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 5:25, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 5 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 IG_dapc <- dapc(AD, inputIG$RiverCode, n.pca = 5)
 # Table of posterior probability results for each individual for each collection.
 IG_MemProb <- as_tibble(round(IG_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputIG$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputIG$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 # SN
 inputSN <-  filter(inputFinal, WatershedCode %in% c("SN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputSN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputSN$Ind, pop=inputSN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputSN$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 5:20).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputSN$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 5:20, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 5 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 SN_dapc <- dapc(AD, inputSN$RiverCode, n.pca = 5)
 # Table of posterior probability results for each individual for each collection.
 SN_MemProb <- as_tibble(round(SN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputSN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputSN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 # WD
 inputWD <-  filter(inputFinal, WatershedCode %in% c("WD")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputWD, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputWD$Ind, pop=inputWD$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputWD$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 20:40).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputWD$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 20:40, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 30 PCs (determined from cross-validation). 6 discriminant functions were chosen.
 WD_dapc <- dapc(AD, inputWD$RiverCode, n.pca = 30)
 # Table of posterior probability results for each individual for each collection.
 WD_MemProb <- as_tibble(round(WD_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputWD$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputWD$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")

 
 #Temp for Osviak------------
 inputFinal_dapcClusters_OS <- filter(inputFinal_dapcClusters, WatershedCode == "OS") %>% 
   select(Ind, dapcCluster)
 
 MetaRBT_OS <- left_join(MetaRBT_OS, inputFinal_dapcClusters_OS) %>% 
   filter(!dapcCluster == "NA") %>% 
   mutate(Length = as.double(Length))
 
 OS_clusters <- MetaRBT_OS %>% 
   group_by(dapcCluster) %>% 
   summarize(Results = t.test(Length))
 
 summarise(ClusterMean = mean(Length))
 
 t.test(Length ~ dapcCluster, MetaRBT_OS)
 
 write.table(MetaRBT_OS,"LatLongTNWR_OS.txt",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",") 
 write_csv(MetaRBT_OS,"LatLongTNWR_OS.csv",quote_escape = F, col_names=T)
 
 
 #Create dapc plot for individuals identified as probable first generation migrants.-------------------
 FirstGenMigMemProb <- filter(MemProb, Ind %in% MemProbFirstGenMig$Ind) %>% 
   filter(PostProb != 0)
 dapc1_MemProbFirstGenPlot <- ggplot(FirstGenMigMemProb, aes(x = Ind, y = PostProb, fill = Group, color = Group)) +
   geom_col() +
   theme_bw() +
   theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
         axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12)) +
   scale_y_continuous(expand = c(0,0), limits = c(0,1.0)) +
   scale_fill_manual(values = WatershedCols[c(5:7,9,10)], aesthetics = c("color", "fill")) +
   facet_grid(. ~ WatershedCode, scales = "free_x", space = "free") +
   labs(x = "Individual", y = "Posterior Membership Probability ")
 

 #Alternative FIGURE 5 plot of intra-drainage PW Fst.------------------- 
 Fig5plot <- ggplot(PW_IntraDrainage, aes(x = Drg, y = Fst, fill = Gtest, color = Gtest)) +
   #Data visualization  
   geom_point(size = 2, shape = 21, position = position_dodge(width = 0.3)) +
   #layout  
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6, color = "black"),
         axis.title.y = element_text(size = 10), axis.text.x = element_text(size = 6, color = "black"), 
         axis.title.x = element_blank(), legend.key = element_rect(size = 0.75),
         legend.key.size = unit(0.5, "lines"), legend.position = c(0.85, 0.9), legend.title = element_text(size = 6), 
         legend.text = element_text(size = 5), legend.background = element_blank()) +
   #overrides 
   ylab("PW_Fst") +
   scale_color_manual(values = c("ns" = "black", "P<0.05" = "black")) +
   scale_fill_manual(values = c("ns" = "white", "P<0.05" = "black")) +
   scale_y_continuous(expand = c(0, 0), breaks = seq(-0.02,0.08,0.02), limits = c(-0.005, 0.07)) +
   annotate("text", x = 5.7, y = 0, label = "3", size = 2) +
   annotate("text", x = 5.7, y = 0.018, label = "2", size = 2) +
   annotate("text", x = 6.3, y = 0.037, label = "2", size = 2) +
   annotate("text", x = 6.3, y = 0.021, label = "2", size = 2) +
   annotate("text", x = 2.3, y = 0.015, label = "2", size = 2)
 
 
 
 #Alternative FIGURE 7 plot of DAPC intra-drainage PW Fst.-------------------  
 Fig7plot2 <- ggplot(PW_IntraDrainage_dapc, aes(x = Drg, y = Fst, fill = Gtest, color = Gtest)) +
   #Data visualization  
   geom_point(size = 2, shape = 21, position = position_dodge(width = 0.3)) +
   #layout  
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6, color = "black"),
         axis.title.y = element_text(size = 10), axis.text.x = element_text(size = 6, color = "black"), 
         axis.title.x = element_blank(), legend.key = element_rect(size = 0.5), axis.line = element_line(colour = "black"), 
         legend.key.size = unit(0.01, "lines"), legend.position = "none", legend.title = element_text(size = 6), 
         legend.text = element_text(size = 5)) +
   #overrides 
   ylab("PW_Fst") +
   scale_color_manual(values = c("ns" = "black", "P<0.05" = "black")) +
   scale_fill_manual(values = c("ns" = "white", "P<0.05" = "black")) +
   scale_y_continuous(expand = c(0, 0), breaks = seq(-0.02,0.11,0.02), limits = c(-0.01, 0.12)) +
   annotate("text", x = 1, y = 0.105, label = "B", size = 3)
 
 
 # Histogram of AMOVA results from Arelquin Fst vs Fct vs Fsc--------------------------------------------------------------------------------------
 F_Statistic <- c("Overall", "Among\nWatersheds", "Within\nWatersheds")
 Value <- c(0.363, 0.350, 0.020)
 FigZtable <- tibble(F_Statistic, Value) %>% 
   mutate(F_Statistic = factor(F_Statistic, levels = c("Overall", "Among\nWatersheds", "Within\nWatersheds")))
 
 FigZPlot <- ggplot(data = FigZtable, aes(x = F_Statistic, y = Value)) +
   #Plot type  
   geom_bar(stat = "identity", color = "black", fill = "gray60") +
   #layout
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
         axis.title.y = element_text(size = 10), axis.text = element_text(color = "black", size = 8), 
         axis.line = element_line(colour = "black"), panel.border = element_blank()) +
   #overrides  
   scale_y_continuous(expand = c(0,0), breaks = seq(0,0.5,0.1), limits = c(0,0.55)) +
   ylab("Population Divergence (Fst)") +
   ggtitle("Rainbow Trout")
 
 # Print histogram of AMOVA results to directory------------------------------------------------------------------------------------------
 ggsave(filename = "FigZ_TNWRRBTPaper.tiff", plot = FigZPlot, device = "tiff",
        height=3.5, width = 4.5, units="in", dpi=600)

 
 # Fst for Rainbow Trout and other salmonids in southwest Alaska----------------------------------------------------------- 

 #Rainbow Trout = 0.363
 #Lake Trout = 0.350 
 #Dolly Varden = 0.018 
 #Coho = 0.026
 #Sockeye = 0.026
 #Chinook = 0.008

 Species <- c("Rainbow\nTrout", "Lake\nTrout", "Dolly\nVarden", "Coho", "Sockeye", "Chinook")
 Value <- c(0.363, 0.350, 0.018, 0.026, 0.026, 0.008)
 FigUtable <- tibble(Species, Value) %>% 
   mutate(Species = factor(Species, levels = c("Rainbow\nTrout", "Lake\nTrout", "Dolly\nVarden", "Coho", "Sockeye", "Chinook")))
 
 FigUPlot <- ggplot(data = FigUtable, aes(x = Species, y = Value)) +
   #Plot type  
   geom_bar(stat = "identity", color = "black", fill = "gray60") +
   #layout
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
         axis.title.y = element_text(size = 10), axis.text = element_text(color = "black", size = 8), 
         axis.line = element_line(colour = "black"), panel.border = element_blank()) +
   #overrides  
   scale_y_continuous(expand = c(0,0), breaks = seq(0,0.5,0.1), limits = c(0,0.55)) +
   ylab("Population Divergence (Fst)") +
   ggtitle("Southwestern Alaska Salmonids") +
   annotate("text", x = 4.5, y = 0.075, label = c("anadromous"), size = 4)
 
 # Print histogram of Fst for other salmonid in southwest Alaska to directory------------------------------------------------------------------------------------------
 ggsave(filename = "FigU_TNWRRBTPaper.tiff", plot = FigUPlot, device = "tiff",
        height=3.5, width = 4.5, units="in", dpi=600)
 
 
 
 # Pairwise Fst among drainages------------------------------------------------------
 # Compute Weir and Cockerham pairwise Fst for all drainage pairs, round to four digits. 
 WCpw_watersheds <- pairwise.WCfst(HF_watershed) %>% 
   round(3)
 WCpw_watersheds[upper.tri(WCpw_watersheds)] <- NA
 Table2 <- as.data.frame(WCpw_watersheds) %>% rownames_to_column("Watershed") %>% as_tibble() %>% select(-(WD))
 write_csv(Table2,"Table2_TNWRRBT.txt",quote_escape = F, col_names=T, na = "")
 
 
 
 
 
 
 # Optional FIGURE: Use DAPC in adegenet to estimate membership probability of individuals in collections within select watersheds.---------
 # KN
 inputKN <- filter(inputFinal, WatershedCode %in% c("KN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputKN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputKN$Ind, pop=inputKN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputKan$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 5:10).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputKan$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 5:10, n.rep = 1000, xval.plot = TRUE)
 KN_dapc <- dapc(AD, inputKN$RiverCode, n.pca = 8)
 KN_MemProb <- as_tibble(round(KN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputKN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputKN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 
 # GN
 inputGN <- filter(inputFinal, WatershedCode %in% c("GN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputGN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputGN$Ind, pop=inputGN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputGN$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 10:30).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputGN$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 10:30, n.rep = 1000, xval.plot = TRUE)
 GN_dapc <- dapc(AD, inputGN$RiverCode, n.pca = 16)
 GN_MemProb <- as_tibble(round(GN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputGN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputGN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 
 # TG
 inputTG <-  filter(inputFinal, WatershedCode %in% c("TG")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputTG, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputTG$Ind, pop=inputTG$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputTG$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 2:25).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputTG$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 2:25, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 2 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 TG_dapc <- dapc(AD, inputTG$RiverCode, n.pca = 2)
 TG_MemProb <- as_tibble(round(TG_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputTG$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputTG$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 
 # UN
 inputUN <-  filter(inputFinal, WatershedCode %in% c("UN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputUN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputUN$Ind, pop=inputUN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputUN$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 2:30).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputUN$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 2:30, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 2 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 UN_dapc <- dapc(AD, inputUN$RiverCode, n.pca = 2)
 UN_MemProb <- as_tibble(round(UN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputUN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputUN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 
 # IG
 inputIG <-  filter(inputFinal, WatershedCode %in% c("IG")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputIG, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputIG$Ind, pop=inputIG$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputIG$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 5:25).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputIG$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 5:25, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 5 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 IG_dapc <- dapc(AD, inputIG$RiverCode, n.pca = 5)
 # Table of posterior probability results for each individual for each collection.
 IG_MemProb <- as_tibble(round(IG_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputIG$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputIG$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 
 # SN
 inputSN <-  filter(inputFinal, WatershedCode %in% c("SN")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputSN, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputSN$Ind, pop=inputSN$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputSN$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 5:20).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputSN$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 5:20, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 5 PCs (determined from cross-validation). 1 discriminant functions were chosen.
 SN_dapc <- dapc(AD, inputSN$RiverCode, n.pca = 5)
 # Table of posterior probability results for each individual for each collection.
 SN_MemProb <- as_tibble(round(SN_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputSN$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputSN$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 
 # WD
 inputWD <-  filter(inputFinal, WatershedCode %in% c("WD")) %>% 
   filter(!(Ind %in% FirstGenMig$Ind))
 AD <- df2genind(select(inputWD, -(Ind:RiverNo)), sep = '/', ncode=4, ind.names=inputWD$Ind, pop=inputWD$RiverCode, NA.char="00/00", type="codom")
 # Use cross-validation function in adegenet to identify the optimal number of PCs to retain for DAPC.
 # Initial cross-validation done to narrow range of the number of PCs evaluated.
 xval1 <- xvalDapc(tab(AD, NA.method = "mean"), inputWD$RiverCode, training.set = 0.9)
 # Final cross-validation evaluates narrow range of number of PCs (in this case 20:40).
 xval2 <- xvalDapc(tab(AD, NA.method = "mean"), inputWD$RiverCode, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = 20:40, n.rep = 1000, xval.plot = TRUE)
 # Perform dapc using 30 PCs (determined from cross-validation). 6 discriminant functions were chosen.
 WD_dapc <- dapc(AD, inputWD$RiverCode, n.pca = 30)
 # Table of posterior probability results for each individual for each collection.
 WD_MemProb <- as_tibble(round(WD_dapc$posterior, 4), rownames = "Ind") %>% 
   mutate(RiverCode = factor(inputWD$RiverCode, levels = RiverFactors)) %>%   
   mutate(WatershedCode = inputWD$WatershedCode) %>% 
   select(Ind, RiverCode, WatershedCode, everything()) %>%   
   pivot_longer(cols = -c("Ind", "RiverCode", "WatershedCode"), names_to = "Group", values_to = "PostProb")
 
 # Combine tables of posterior probabilities for each drainage
 AllAggr_MemProb <- bind_rows(KN_MemProb, GN_MemProb, TG_MemProb, UN_MemProb, IG_MemProb, SN_MemProb, WD_MemProb) %>% 
   mutate(Group = factor(Group, levels = c("KN1", "KN2", "GN1", "GN2", "GN3", "GN4", "TG1", "TG2", "UN1", "UN2", "IG1", "IG2", "SN1", "SN2", "WD1", 
                                           "WD2", "WD3", "WD4", "WD5", "WD6", "WD7"))) %>% 
   mutate(WatershedCode = factor(WatershedCode, levels = c("KN", "GN", "TG", "UN", "IG", "SN", "WD")))
 
 saveRDS(AllAggr_MemProb, "AllAggr_MemProb.rds")
 AllAggr_MemProb <- readRDS("AllAggr_MemProb.rds")
 
 # Filter maximum posterior probability for each individual and add a column for classification at 90% threshold assigned. < 90% classified as admixed.
 AllAggr_PostProb <- group_by(AllAggr_MemProb, RiverCode, Ind) %>% filter(PostProb == max(PostProb)) %>% 
   mutate(Classification = ifelse(PostProb < 0.90, "Admixed", ifelse(RiverCode == Group, "Home", "Other"))) %>% 
   mutate(Classification = factor(Classification, levels = c("Home", "Other", "Admixed")))
 # Tally the number and freq or individuals classified as home, other, or admixed for each aggregation
 AllAggr_ProbSum <- AllAggr_PostProb %>% 
   group_by(WatershedCode, RiverCode, Classification) %>% 
   tally() %>% 
   mutate(freq = n/sum(n)) %>% 
   rename(Drainage = WatershedCode, Aggregation = RiverCode)
 #Plot posterior probability summary. Two plots were used because facet_grid was used to allow for free space and scales (facet_wrap does not permit free space) 
 AllAggr_Plot1 <- ggplot(filter(AllAggr_ProbSum, Drainage %in% c("KN", "GN", "TG", "UN")), aes(x = Aggregation, y = freq, fill = Classification)) +
   geom_col(color = "black", position = position_fill(reverse = T)) +
   theme_bw() + 
   facet_grid(~Drainage, scales = "free_x", space = "free") +
   scale_fill_manual(values = c("black", "gray", "white")) +   
   theme(legend.position = "none", legend.title = element_blank(), axis.title = element_blank(), 
         axis.text = element_text(size = 5, color = "black"), strip.text = element_text(size = 6), 
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
 AllAggr_Plot2 <- ggplot(filter(AllAggr_ProbSum, Drainage %in% c("IG", "SN", "WD")), aes(x = Aggregation, y = freq, fill = Classification)) +
   geom_col(color = "black", position = position_fill(reverse = T)) +
   theme_bw() + 
   facet_grid(~Drainage, scales = "free_x", space = "free") +
   scale_fill_manual(values = c("black", "gray", "white")) +   
   theme(legend.position = c(0.5, -0.30), legend.title = element_blank(), legend.text = element_text(size = 7), legend.direction = "horizontal", 
         legend.key.size = unit(0.10, "cm"), legend.key.width = unit(0.25, "cm"), axis.title = element_blank(), 
         axis.text = element_text(size = 5, color = "black"), strip.text = element_text(size = 6),panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())
 
 # Figure X - Combine plots one and two.  
 FigXPlot <- ggarrange(AllAggr_Plot1, AllAggr_Plot2, ncol = 1, nrow = 3, 
                       heights = c(6, 6, 1), common.legend = FALSE)  %>% 
   annotate_figure(left = text_grob("Percent of Aggregation", rot = 90, size = 8), bottom = text_grob("Aggregation", vjust = -4.0, size = 8))
 
 #geom_col(position = position_dodge2(preserve = "single", width = 10), color = "black") +
 #labs(fill = "Assigned")
 
 
 # Optional: Print FIGURE to directory---------------------------------------------------------------------------------
 ggsave(filename = "FigX_TNWRRBTPaper.tiff", plot = Fig4Plot, device = "tiff",
        height=4.5, width = 3.5, units="in", dpi=600)
 
 
 
 # Optional: This Figure 5 was used in the original manuscript and was replaced with new Figure based on reviewer comment------------
 Fig5plot_old <- ggplot(PW_IntraDrainage, aes(x = Pair, y = Fst, fill = Gtest, color = Gtest)) +
   #Data visualization  
   geom_bar(stat = "identity", width = 0.50) +
   #layout  
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6, color = "black"),
         axis.title.y = element_text(size = 10), axis.text.x = element_text(size = 6, color = "black"), 
         axis.title.x = element_text(size = 10), legend.key = element_rect(size = 0.75),
         legend.key.size = unit(0.50, "lines"), legend.position = c(0.85, 0.9), legend.title = element_text(size = 6), 
         legend.text = element_text(size = 5), legend.background = element_blank()) +
   #overrides 
   scale_x_discrete(limits = rev(levels(PW_IntraDrainage$Pair))) +
   coord_flip() +  
   ylab("PW_Fst") +
   scale_color_manual(values = c("ns" = "black", "P<0.05" = "black")) +
   scale_fill_manual(values = c("ns" = "white", "P<0.05" = "black")) +
   scale_y_continuous(expand = c(0, 0), breaks = seq(0.00,0.08,0.02), limits = c(0.00, 0.07)) 
 