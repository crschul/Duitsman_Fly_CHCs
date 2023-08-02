################################################
# Final Results and Figures
# Duitsman: Fly CHC Analysis
# Corey Schultz - University of Georgia - 2023
################################################

# Set Working Directories
setwd("C:/RscriptsAndrew")
datawd <- "C:/RscriptsAndrew"
figwd <- "C:/RscriptsAndrew/NewFiguresAndrew"

# Packages used
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(vegan)
library(car)
library(Hotelling)

# Set scientific notation
options(scipen=999,digits = 4)
#############################################################

### Read in and format the data
fixed_peaks <- read.csv("C:/RscriptsAndrew/Fixed_Sexual_Dimorphism_Peaks.csv",header = TRUE,check.names=FALSE)

male_peaks <- filter(fixed_peaks, Sex == "Male")
female_peaks <- filter(fixed_peaks, Sex == "Female")
female_peaks$Run <-sub("^","Run ",female_peaks$Run)


# loop to find the average
column_names <- colnames(fixed_peaks)

# male
male_calc_withZero <- setNames(data.frame(matrix(0,ncol = length(fixed_peaks), nrow = length(unique(fixed_peaks$Species)))), column_names)
male_calc_withZero <- subset(male_calc_withZero, select=-c(Run))
male_calc_noZero <- setNames(data.frame(matrix(0,ncol = length(fixed_peaks), nrow = length(unique(fixed_peaks$Species)))), column_names)
male_calc_noZero <- subset(male_calc_noZero, select=-c(Run))

mspecies <- unique(male_peaks$Species)

for (i in 1:length(mspecies)) {
  species_df <- filter(male_peaks, Species == mspecies[i])
  species_df <- subset(species_df, select=-c(Run))
  
  male_calc_withZero[i,1:5] <- species_df[1,1:5]
  male_calc_noZero[i,1:5] <- species_df[1,1:5]
  print(i)
  
  for (j in 6:89){
    male_calc_withZero[i,j] <- (sum(species_df[,j])/length(species_df[,j]))
    male_calc_noZero[i,j] <- (sum(species_df[,j])/sum(species_df[,j] != 0)) 
  }
  
}


# female
female_calc_withZero <- setNames(data.frame(matrix(0,ncol = length(fixed_peaks), nrow = length(unique(fixed_peaks$Species)))), column_names)
female_calc_withZero <- subset(female_calc_withZero, select=-c(Run))
female_calc_noZero <- setNames(data.frame(matrix(0,ncol = length(fixed_peaks), nrow = length(unique(fixed_peaks$Species)))), column_names)
female_calc_noZero <- subset(female_calc_noZero, select=-c(Run))

fspecies <- unique(female_peaks$Species)

for (i in 1:length(fspecies)) {
  species_df <- filter(female_peaks, Species == mspecies[i])
  species_df <- subset(species_df, select=-c(Run))
  
  female_calc_withZero[i,1:5] <- species_df[1,1:5]
  female_calc_noZero[i,1:5] <- species_df[1,1:5]
  print(i)
  
  for (j in 6:89){
    female_calc_withZero[i,j] <- (sum(species_df[,j])/length(species_df[,j]))
    female_calc_noZero[i,j] <- (sum(species_df[,j])/sum(species_df[,j] != 0)) 
  }
  
}

# combine them
peak_means_withZero <- rbind(male_calc_withZero,female_calc_withZero)
peak_means_noZero <- rbind(male_calc_noZero,female_calc_noZero) %>% mutate_all(~ifelse(is.nan(.), 0, .))

the_names <- peak_means_noZero[c("Species", "Sex")]
the_names$label <- paste(the_names$Species,the_names$Sex, sep = "_")

euclidean_withZero <- as.data.frame(as.matrix(vegdist(peak_means_withZero[6:89], method = "euclidean")))
rownames(euclidean_withZero) <- the_names$label
colnames(euclidean_withZero) <- the_names$label

euclidean_noZero <- as.data.frame(as.matrix(vegdist(peak_means_noZero[6:89], method = "euclidean")))
rownames(euclidean_noZero) <- the_names$label
colnames(euclidean_noZero) <- the_names$label

# save these dataframes
# write.csv(peak_means_withZero,"peak_means_withZero.csv",row.names = FALSE)
# write.csv(peak_means_noZero,"peak_means_noZero.csv",row.names = FALSE)
# write.csv(euclidean_withZero,"euclidean_withZero.csv",row.names = FALSE)
# write.csv(euclidean_noZero,"euclidean_noZero.csv",row.names = FALSE)
# 
# euclidean_fixed_noMean <- as.data.frame(as.matrix(vegdist(fixed_peaks[7:90], method = "euclidean")))
# write.csv(euclidean_fixed_noMean,"euclidean_fixed_noMean.csv",row.names = FALSE)


######################################
### Analysis and Figures Generated ###
#
# 1. Heatmap of presence/absence
# 2. Permanova of CHCs and metadata
# 3. RDA to visualize CHC clusters
# 4. Mantel test of CHCs vs genetic distance
# 5. PCAs for Supplemental Informations
#######################################
# The Data Frames
peak_means_noZero
peak_means_withZero

male_nozero <- filter(peak_means_noZero, Sex == "Male")
female_nozero <- filter(peak_means_noZero, Sex == "Female")

male_withzero <- filter(peak_means_withZero, Sex == "Male")
female_withzero <- filter(peak_means_withZero, Sex == "Female")

true_means <- peak_means_withZero[,!names(peak_means_withZero) %in% 
                                    c('-8.34', '4.92','5.01', '5.07', '6.66', '6.75', '6.9')]

####### 1. Heatmap  ######################################

### all species

df <- peak_means_withZero %>% arrange(Phylogenetic.Group)

long_data <- peak_means_withZero %>% 
  gather(key="Peaks", value='value',-Species,-Sex, -Phylogenetic.Group,-Binary.Sex, -Dietary.State) %>%
  setNames(c("Species","Sex","Phylogenetic.Group","Binary.Sex","Dietary.State","Peaks", "value")) %>%
  mutate(Species=factor(Species)) %>%
  mutate(Sex=factor(Sex)) %>%
  mutate(Peaks=factor(Peaks)) %>%
  mutate(value=as.numeric(value))
long_data$Species_Fly <- paste(long_data$Species,long_data$Sex,sep = "_")

long_data["value"][long_data["value"] == 0] <- NA
long_data <- long_data %>% arrange(Phylogenetic.Group)

groups <- df$Phylogenetic.Group

groups <- sub("1","dark green",groups)
groups <- sub("2","purple",groups)
groups <- sub("3","springgreen",groups)
groups <- sub("4","dodgerblue",groups)
groups <- sub("5","firebrick",groups)

man_colors <- c("dodgerblue","dodgerblue",
                "dodgerblue","dodgerblue",
                "dodgerblue", "dodgerblue",
                "purple","purple",
                "purple","purple",
                "purple","purple",
                "purple","purple",
                "purple","purple",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "dark green","dark green",
                "springgreen", "springgreen",
                "springgreen","springgreen",
                "springgreen","springgreen",
                "springgreen","springgreen",
                "firebrick","firebrick")

species_order <- c("Tripunctata_Male","Tripunctata_Female",
                   "Cardini_Male","Cardini_Female",
                   "Acutilabella_Male", "Acutilabella_Female",
                   "Bizonata_Male","Bizonata_Female",
                   "Putrida_Male","Putrida_Female",
                   "Neotestacea_Male","Neotestacea_Female",
                   "Testacea_Male","Testacea_Female",
                   "Orientacea_Male","Orientacea_Female",
                   "Nigromaculata_Male","Nigromaculata_Female",
                   "Palustrus_Male","Palustrus_Female",
                   "Subpalustrus_Male","Subpalustrus_Female",
                   "Deflecta_Male","Deflecta_Female",
                   "Quinaria_Male","Quinaria_Female",
                   "Tenebrosa_Male","Tenebrosa_Female",
                   "Occidentalis_Male","Occidentalis_Female",
                   "Suboccidentalis_Male","Suboccidentalis_Female",
                   "Recens_Male","Recens_Female",
                   "Transversa_Male","Transversa_Female",
                   "Subquinaria_Male","Subquinaria_Female",
                   "Angularis_Male", "Angularis_Female",
                   "Falleni_Male","Falleni_Female",
                   "Innubila_Male","Innubila_Female",
                   "Phalerata_Male","Phalerata_Female",
                   "Immigrans_Male","Immigrans_Female")


p <- long_data %>%  ggplot(aes(x=Peaks, y=Species_Fly, fill=value))+
  geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10, face="bold")) +
  theme(axis.text.y = element_text(size = 8, color = rev(man_colors))) +
  scale_x_discrete(limits = names(peak_means_withZero)[6:ncol(peak_means_withZero)]) +
  scale_fill_gradient2(low = "tan",
                       mid = "orange",
                       high = "blue",
                       na.value="white") +
  ggtitle("Area Percent - All Species") + labs(fill="Area Percent") + scale_y_discrete(limits = rev(species_order)) 
p

ggsave("Heatmap_all_DEFENSE.png",
       p, device = "png", width = 16, height = 8, dpi = 700)

### Just Quinaria

df <- true_means %>% filter(Phylogenetic.Group == 1)

long_data <- df %>% 
  gather(key="Peaks", value='value',-Species,-Sex, -Phylogenetic.Group,-Binary.Sex, -Dietary.State) %>%
  setNames(c("Species","Sex","Phylogenetic.Group","Binary.Sex","Dietary.State","Peaks", "value")) %>%
  mutate(Species=factor(Species)) %>%
  mutate(Sex=factor(Sex)) %>%
  mutate(Peaks=factor(Peaks)) %>%
  mutate(value=as.numeric(value))
long_data$Species_Fly <- paste(long_data$Species,long_data$Sex,sep = "_")

long_data["value"][long_data["value"] == 0] <- NA
long_data <- long_data %>% arrange(Phylogenetic.Group)

species_order <- c(
  "Nigromaculata_Male","Nigromaculata_Female",
  "Palustrus_Male","Palustrus_Female",
  "Subpalustrus_Male","Subpalustrus_Female",
  "Deflecta_Male","Deflecta_Female",
  "Quinaria_Male","Quinaria_Female",
  "Tenebrosa_Male","Tenebrosa_Female",
  "Occidentalis_Male","Occidentalis_Female",
  "Suboccidentalis_Male","Suboccidentalis_Female",
  "Recens_Male","Recens_Female",
  "Transversa_Male","Transversa_Female",
  "Subquinaria_Male","Subquinaria_Female"
)


p <- long_data %>%  ggplot(aes(x=Peaks, y=Species_Fly, fill=value))+
  geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10, face="bold")) +
  theme(axis.text.y = element_text(size = 8, color = "dark green")) +
  scale_x_discrete(limits = names(true_means)[6:ncol(true_means)]) +
  scale_fill_gradient2(low = "tan",
                       mid = "orange",
                       high = "blue",
                       na.value="white") +
  ggtitle("Area Percent - True CHCs - Quinaria") + labs(fill="Area Percent") + scale_y_discrete(limits = species_order) 
p
ggsave("Heatmap_Quinaria_True_DEFENSE.png",
       p, device = "png", width = 16, height = 8, dpi = 700)

#############################


################################
##### 2. PERMANOVA
################################

## make df with pseudocount, clr, and then scale
combined_peaks <- rbind(male_peaks,female_peaks)

comb_pseudo <- combined_peaks
comb_pseudo[,7:90] <- comb_pseudo[,7:90] + .001

comb_clr <- comb_pseudo
comb_clr[,7:90] <- clr(comb_clr[,7:90])

comb_scale <- comb_clr
comb_scale[,7:90] <- scale(comb_scale[,7:90], center = FALSE, scale = TRUE)

comb_scale$Phylogenetic.Group <- as.factor(comb_scale$Phylogenetic.Group)
comb_scale$Dietary.State <- as.factor(comb_scale$Dietary.State)

df <- comb_scale
EucDist <- vegdist(df[,7:90], method="euclidean")


# Individual
adonis2(EucDist ~ df$Sex, by = "terms", parallel = 4)
adonis2(EucDist ~ df$Species, by = "terms", parallel = 4)
adonis2(EucDist ~ df$Phylogenetic.Group, by = "terms", parallel = 4)
adonis2(EucDist ~ df$Dietary.State, by = "terms", parallel = 4)
adonis2(EucDist ~ df$Run, by = "terms", parallel = 4) 

# Complex
adonis2(EucDist ~ df$Dietary.State/df$Species + df$Phylogenetic.Group/df$Species
        + df$Dietary.State/df$Sex + df$Phylogenetic.Group/df$Sex + df$Species*df$Sex,
        by = "terms")
adonis2(EucDist ~ df$Phylogenetic.Group/df$Species + df$Dietary.State/df$Species
        + df$Dietary.State/df$Sex + df$Phylogenetic.Group/df$Sex + df$Species*df$Sex,
        by = "terms")
adonis2(EucDist ~ df$Species*df$Sex + df$Phylogenetic.Group/df$Species + df$Dietary.State/df$Species
        + df$Dietary.State/df$Sex + df$Phylogenetic.Group/df$Sex,
        by = "terms")

# By Sex
male_pdf <- filter(comb_scale, Sex == "Male")
male_euc <- vegdist(male_pdf[,7:90], method="euclidean")

female_pdf <- filter(comb_scale, Sex == "Female")
female_euc <- vegdist(female_pdf[,7:90], method="euclidean")

adonis2(male_euc ~ male_pdf$Dietary.State + male_pdf$Phylogenetic.Group + male_pdf$Phylogenetic.Group/male_pdf$Species,
        by = "margin")
adonis2(female_euc ~ female_pdf$Dietary.State + female_pdf$Phylogenetic.Group + female_pdf$Phylogenetic.Group/female_pdf$Species,
        by = "terms")
#######################################################################

################################
##### 3. RDA Figures
################################

# DFs
combined_peaks
peak_means_withZero
setwd(figwd)

df <- comb_scale
rda1 <- rda(df[7:90] ~ df$Species*df$Sex )

f1 <- scores(rda1, display = "sites") %>% as.data.frame() %>% 
  cbind(df[,1:6])
f1$Phylogenetic.Group <- as.factor(f1$Phylogenetic.Group)

# This is how you use variance explained in th legends
rsum <- summary(rda1)
var1 <- rsum$cont$importance[2, "RDA1"]
var2 <- rsum$cont$importance[2, "RDA2"]
var1 <- format(round(var1, 3), nsmall = 3)
var2 <- format(round(var2, 3), nsmall = 3)

all_rda <- ggplot(data = f1, aes(x = RDA1, y = RDA2, shape = Sex, color = Phylogenetic.Group)) + geom_point(aes(size = 2)) +
  ggtitle("CHC Peaks RDA by Sex and Phylogenetic Group") + scale_color_manual(values=c("dark green","purple","springgreen","dodgerblue","firebrick")) +
  xlab(paste(var1," variance explained")) +
  ylab(paste(var2," variance explained"))

all_rda
ggsave("All_RDA.png",
       all_rda, device = "png", width = 8, height = 8, dpi = 700)

# male and female
setwd(figwd)
male_rda <- filter(f1,Sex == "Male")
female_rda <- filter(f1,Sex == "Female")


m_rda <- ggplot(data = male_rda, aes(x = RDA1, y = RDA2, color = Phylogenetic.Group)) + geom_point(aes(size = 2)) +
  ggtitle("CHC Peaks RDA by Phylogenetic Group - Males") + scale_color_manual(values=c("dark green","purple","springgreen","dodgerblue","firebrick"))
m_rda
ggsave("Male_RDA.png",
       m_rda, device = "png", width = 8, height = 8, dpi = 700)

f_rda <- ggplot(data = female_rda, aes(x = RDA1, y = RDA2, color = Phylogenetic.Group)) + geom_point(aes(size = 2)) +
  ggtitle("CHC Peaks RDA by Phylogenetic Group - Females") + scale_color_manual(values=c("dark green","purple","springgreen","dodgerblue","firebrick"))
f_rda
ggsave("Female_RDA.png",
       f_rda, device = "png", width = 8, height = 8, dpi = 700)

####################################### by means
means_pseudo <- peak_means_withZero
means_pseudo[,6:89] <- means_pseudo[,6:89] + .01

means_clr <- means_pseudo
means_clr[,6:89] <- clr(means_clr[,6:89])

means_scale <- means_clr
means_scale[,6:89] <- scale(means_scale[,6:89], center = FALSE, scale = TRUE)

df <- means_scale
rda1 <- rda(df[,6:89] ~ df$Species*df$Sex )

f1 <- scores(rda1, display = "sites") %>% as.data.frame() %>% 
  cbind(df[,1:5])
f1$Phylogenetic.Group <- as.factor(f1$Phylogenetic.Group)

rsum <- summary(rda1)
var1 <- rsum$cont$importance[2, "RDA1"]
var2 <- rsum$cont$importance[2, "RDA2"]
var1 <- format(round(var1, 3), nsmall = 3)
var2 <- format(round(var2, 3), nsmall = 3)

# mean data not clt
# all_mean_rda <- ggplot(data = f1, aes(x = RDA1, y = RDA2, shape = Sex, color = Phylogenetic.Group)) + geom_point(aes(size = 2)) +
#   ggtitle("CHC Peaks RDA by Sex and Phylogenetic Group") + scale_color_manual(values=c("dark green","purple","springgreen","dodgerblue","firebrick"))

all_mean_rda <- ggplot(data = f1, aes(x = RDA1, y = RDA2, shape = Sex, color = Species)) + geom_point(aes(size = 2)) +
  ggtitle("CHC Peaks RDA by Sex and Species") +
  xlab(paste(var1," variance explained")) +
  ylab(paste(var2," variance explained"))

all_mean_rda

ggsave("All_RDA_Means_Species.png",
       all_mean_rda, device = "png", width = 8, height = 8, dpi = 700)


####### subset by sex
df <- means_scale %>% filter(Sex == "Male")
rda1 <- rda(df[6:89] ~ df$Species )

f1 <- scores(rda1, display = "sites") %>% as.data.frame() %>% 
  cbind(df[,1:5])
f1$Phylogenetic.Group <- as.factor(f1$Phylogenetic.Group)
f1$Dietary.State <- as.factor(f1$Dietary.State)

# all_rda <- ggplot(data = f1, aes(x = RDA1, y = RDA2, color = Phylogenetic.Group)) + geom_point(aes(size = 2)) +
#   ggtitle("CHC Peaks RDA by Phylogenetic Group: Females") + scale_color_manual(values=c("dark green","purple","springgreen","dodgerblue","firebrick"))
rsum <- summary(rda1)
var1 <- rsum$cont$importance[2, "RDA1"]
var2 <- rsum$cont$importance[2, "RDA2"]
var1 <- format(round(var1, 3), nsmall = 3)
var2 <- format(round(var2, 3), nsmall = 3)

all_rda <- ggplot(data = f1, aes(x = RDA1, y = RDA2, shape = Sex, color = Species)) + geom_point(aes(size = 2)) +
  ggtitle("CHC Peaks RDA by Species: Males") +
  xlab(paste(var1," variance explained")) +
  ylab(paste(var2," variance explained"))
# geom_text(label = f1$Species,nudge_x = 0.15, nudge_y = 0.15, check_overlap = T)
all_rda

ggsave("Male_RDA_mean_Species.png",
       all_rda, device = "png", width = 8, height = 8, dpi = 700)

#############################################################


############################
#### 4. Mantel Test
############################
setwd("C:/RscriptsAndrew/Distance Matrix")

genetic_distance <- as.dist(read.csv("Test_Dist.csv", header = TRUE, sep = ","))
quin_dist <- as.dist(read.csv("Correct Matrix/NoMT_Distance.csv", header = TRUE, sep = ","))
# With or without all metabolites
peak_means_withZero
true_means

#means all
means_male <- filter(true_means, Sex == "Male") 
means_female <- filter(true_means, Sex == "Female")

euc_male <- vegdist(means_male[,6:82], method="euclidean")
euc_female <- vegdist(means_female[,6:82], method="euclidean")


mantel(euc_female, euc_male, method="pearson")
mantel(euc_male, genetic_distance, method="pearson")
mantel(euc_female, genetic_distance, method="pearson")

################################# just quin
quin_dist <- as.dist(read.csv("Correct Matrix/Genetic_Distance_Matrix_Quinaria_Group.csv", header = TRUE, sep = ","))
quin_dist <- as.dist(read.csv("Correct Matrix/NoMT_Quinaria_Group.csv", header = TRUE, sep = ","))

group1 <- filter(true_means, Phylogenetic.Group == "1")

means_male <- filter(group1, Sex == "Male") 
means_female <- filter(group1, Sex == "Female")

euc_male <- vegdist(means_male[,6:82], method="euclidean")
euc_female <- vegdist(means_female[,6:82], method="euclidean")

mantel(euc_female, euc_male, method="pearson")
mantel(euc_male, quin_dist, method="pearson")
mantel(euc_female, quin_dist, method="pearson")

#########################################################################

##############################
##### 5. Species Specific PCAs
##############################
setwd("C:/RscriptsAndrew/NewFiguresAndrew/True_clt_PCA_species")
df <- true_CHCs

for (s in unique(df$Species)) {
  print(s)
  
  subset_df <- filter(df, Species == s)
  rda1 <- rda(subset_df[7:82])
  
  f1 <- scores(rda1, display = "sites") %>% as.data.frame() %>% 
    cbind(subset_df[,1:6])
  f1$Phylogenetic.Group <- as.factor(f1$Phylogenetic.Group)
  
  
  species_pc <- ggplot(data = f1, aes(x = PC1, y = PC2, color = Sex)) + geom_point(aes(size = 2)) +
    ggtitle(paste("CHC Peaks PCA by Sex: ", s))
  ggsave(paste(s,"_PCA.png"),
         species_pc, device = "png", width = 8, height = 8, dpi = 700)
  
}

setwd("C:/RscriptsAndrew/NewFiguresAndrew/True_clt_PCA_phylo")
df <- true_CHCs

for (s in unique(df$Phylogenetic.Group)) {
  print(s)
  
  subset_df <- filter(df, Phylogenetic.Group == s)
  rda1 <- rda(subset_df[7:82])
  
  f1 <- scores(rda1, display = "sites") %>% as.data.frame() %>% 
    cbind(subset_df[,1:6])
  f1$Phylogenetic.Group <- as.factor(f1$Phylogenetic.Group)
  
  
  species_pc <- ggplot(data = f1, aes(x = PC1, y = PC2, color = Sex)) + geom_point(aes(size = 2)) +
    ggtitle(paste("CHC Peaks PCA by Sex: ", s))
  ggsave(paste("PhyloGroup_",s,"_PCA.png"),
         species_pc, device = "png", width = 8, height = 8, dpi = 700)
  
}

########################################################################

##########################################
##### 6. Check for true sexual dimporphism
###########################################
setwd("C:/RscriptsAndrew/NewFiguresAndrew/SexualDimorphism")
fixed_peaks

### All CHCs
long_sex <- fixed_peaks %>% 
  gather(key="Peaks", value='value',-Species,-Run, -Sex,-Phylogenetic.Group, -Binary.Sex, -Dietary.State) %>%
  setNames(c("Species","Run","Sex","Phylogenetic.Group","Binary.Sex","Dietary.State", "Peaks","value")) %>%
  mutate(Species=factor(Species)) %>%
  mutate(Sex=factor(Sex)) %>%
  mutate(Peaks=factor(Peaks)) %>%
  mutate(value=as.numeric(value))

# Are male and Female Flies Just Different? : Yes! - just the mean of the df though
t.test(value ~ Sex, data = long_sex, alternative = "two.sided")

anova(lm(value ~ Sex, data = long_sex))

mspecies <- unique(long_sex$Species)
npeaks <- unique(long_sex$Peaks)

# empty dfs
species_dimorphism_ttest_level2 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(species_dimorphism_ttest_level2) <- c('Species','P_value')

species_dimorphism_ttest_level3 <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(species_dimorphism_ttest_level3) <- c('Species','Peak','P_value')

species_dimorphism_anova_level2 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(species_dimorphism_anova_level2) <- c('Species','P_value')

species_dimorphism_anova_level3 <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(species_dimorphism_anova_level3) <- c('Species','Peak','P_value')

for (i in 1:length(mspecies)) {
  print("_______________________")
  print(toString(mspecies[i]))
  species_long <- filter(long_sex, Species == mspecies[i])
  
  t_res <- t.test(value ~ Sex, data = species_long, alternative = "two.sided")
  anv <- anova(lm(value ~ Sex, data = species_long))
  anv_2_pvalue <- anv$`Pr(>F)`[1]
  
  t2row <- c(toString(mspecies[i]),t_res$p.value)
  a2row <- c(toString(mspecies[i]),anv_2_pvalue)
  
  species_dimorphism_ttest_level2[nrow(species_dimorphism_ttest_level2) + 1, ] <- t2row
  species_dimorphism_anova_level2[nrow(species_dimorphism_anova_level2) + 1, ] <- a2row
  
  
  for (p in 1:length(npeaks)) {
    species_peak <- filter(species_long, Peaks == npeaks[p])
    
    t_res <- t.test(value ~ Sex, data = species_peak, alternative = "two.sided")
    anv <- anova(lm(value ~ Sex, data = species_peak))
    anv_3_pvalue <- anv$`Pr(>F)`[1]
    
    t3row <- c(toString(mspecies[i]), toString(species_peak[1,7]),t_res$p.value)
    a3row <- c(toString(mspecies[i]), toString(species_peak[1,7]),t_res$p.value)
    
    species_dimorphism_ttest_level3[nrow(species_dimorphism_ttest_level3) + 1, ] <- t3row
    species_dimorphism_anova_level3[nrow(species_dimorphism_anova_level3) + 1, ] <- a3row
    
  }
}

write.csv(species_dimorphism_anova_level2,file = "Level2_anova_TotalCHCs.csv")
write.csv(species_dimorphism_ttest_level2,file = "Level2_ttest_TotalCHCs.csv")



####### p adjust
species_dimorphism_ttest_level3$P_adjust <- p.adjust(species_dimorphism_ttest_level3$P_value, method = "BH")
species_dimorphism_ttest_level3 <- filter(species_dimorphism_ttest_level3, P_adjust <= .05)
unique(species_dimorphism_ttest_level3$Species)
write.csv(species_dimorphism_ttest_level3,file = "Level3_ttest_TotalCHCs.csv")

species_dimorphism_anova_level3$P_adjust <- p.adjust(species_dimorphism_anova_level3$P_value, method = "BH")
species_dimorphism_anova_level3 <- filter(species_dimorphism_anova_level3, P_adjust <= .05)
unique(species_dimorphism_anova_level3$Species)
write.csv(species_dimorphism_anova_level3,file = "Level3_anova_TotalCHCs.csv")




### True CHCs
true_fixed <- fixed_peaks[,!names(fixed_peaks) %in% 
                            c('-8.34', '4.92','5.01', '5.07', '6.66', '6.75', '6.9')]

long_sex <- true_fixed %>% 
  gather(key="Peaks", value='value',-Species,-Run, -Sex,-Phylogenetic.Group, -Binary.Sex, -Dietary.State) %>%
  setNames(c("Species","Run","Sex","Phylogenetic.Group","Binary.Sex","Dietary.State", "Peaks","value")) %>%
  mutate(Species=factor(Species)) %>%
  mutate(Sex=factor(Sex)) %>%
  mutate(Peaks=factor(Peaks)) %>%
  mutate(value=as.numeric(value))

# Are male and Female Flies Just Different? : Yes! - just the mean of the df though
t.test(value ~ Sex, data = long_sex, alternative = "two.sided")

anova(lm(value ~ Sex, data = long_sex))

mspecies <- unique(long_sex$Species)
npeaks <- unique(long_sex$Peaks)

# empty dfs
species_dimorphism_ttest_level2 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(species_dimorphism_ttest_level2) <- c('Species','P_value')

species_dimorphism_ttest_level3 <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(species_dimorphism_ttest_level3) <- c('Species','Peak','P_value')

species_dimorphism_anova_level2 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(species_dimorphism_anova_level2) <- c('Species','P_value')

species_dimorphism_anova_level3 <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(species_dimorphism_anova_level3) <- c('Species','Peak','P_value')

for (i in 1:length(mspecies)) {
  print("_______________________")
  print(toString(mspecies[i]))
  species_long <- filter(long_sex, Species == mspecies[i])
  
  t_res <- t.test(value ~ Sex, data = species_long, alternative = "two.sided")
  anv <- anova(lm(value ~ Sex, data = species_long))
  anv_2_pvalue <- anv$`Pr(>F)`[1]
  
  t2row <- c(toString(mspecies[i]),t_res$p.value)
  a2row <- c(toString(mspecies[i]),anv_2_pvalue)
  
  species_dimorphism_ttest_level2[nrow(species_dimorphism_ttest_level2) + 1, ] <- t2row
  species_dimorphism_anova_level2[nrow(species_dimorphism_anova_level2) + 1, ] <- a2row
  
  
  for (p in 1:length(npeaks)) {
    species_peak <- filter(species_long, Peaks == npeaks[p])
    
    t_res <- t.test(value ~ Sex, data = species_peak, alternative = "two.sided")
    anv <- anova(lm(value ~ Sex, data = species_peak))
    anv_3_pvalue <- anv$`Pr(>F)`[1]
    
    t3row <- c(toString(mspecies[i]), toString(species_peak[1,7]),t_res$p.value)
    a3row <- c(toString(mspecies[i]), toString(species_peak[1,7]),t_res$p.value)
    
    species_dimorphism_ttest_level3[nrow(species_dimorphism_ttest_level3) + 1, ] <- t3row
    species_dimorphism_anova_level3[nrow(species_dimorphism_anova_level3) + 1, ] <- a3row
    
  }
}

write.csv(species_dimorphism_anova_level2,file = "Level2_anova_TrueCHCs.csv")
write.csv(species_dimorphism_ttest_level2,file = "Level2_ttest_TrueCHCs.csv")



####### p adjust
species_dimorphism_ttest_level3$P_adjust <- p.adjust(species_dimorphism_ttest_level3$P_value, method = "BH")
species_dimorphism_ttest_level3 <- filter(species_dimorphism_ttest_level3, P_adjust <= .05)
unique(species_dimorphism_ttest_level3$Species)
write.csv(species_dimorphism_ttest_level3,file = "Level3_ttest_TrueCHCs.csv")

species_dimorphism_anova_level3$P_adjust <- p.adjust(species_dimorphism_anova_level3$P_value, method = "BH")
species_dimorphism_anova_level3 <- filter(species_dimorphism_anova_level3, P_adjust <= .05)
unique(species_dimorphism_anova_level3$Species)
write.csv(species_dimorphism_anova_level3,file = "Level3_anova_TrueCHCs.csv")


# True vs Total

#Welches T.Test and anova

# Level 1: All Species + All Peaks: Male vs Female
# Level 2: 1 Species + All Peaks: Male vs Female - Table
# Level 3: 1 Species + 1 Peak: Male vs Female - Table

###################################################


######################################
##### 7. Set up comet and broken stick
######################################
setwd(datawd)

means_clr

true_means_clr <- means_clr[,!names(means_clr) %in% 
                              c('-8.34', '4.92','5.01', '5.07', '6.66', '6.75', '6.9')]

cometdf <- filter(true_means_clr, Phylogenetic.Group == "1")

means_rda <- rda(cometdf[,6:82])
rda_scores <- summary(means_rda)
rda_scores$sites

pca_df <- cbind(cometdf$Species,cometdf$Sex)
colnames(pca_df) =  c('Species','Sex')
pca_df <- as.data.frame(cbind(pca_df,rda_scores$sites))

pca_male <- filter(pca_df, Sex == "Male")
pca_female <- filter(pca_df, Sex == "Female")

setwd(datawd)

write.csv(pca_male,file = "CoMet_Stuff/Q_Male_PCAs_True.csv")
write.csv(pca_female,file = "CoMet_Stuff/Q_Female_PCAs_True.csv")




################################ Broken Stick
bstick(means_rda)
screeplot(means_rda, bstick = TRUE, type = "lines")











