###################################################################
############# Juvenile wild morphometrics analyses ################
###################################################################

########################### Housekeeping ##############################
# Set working director
setwd("C:/Users/nanda/Desktop/Island/Morph Data/R code")

rm(list=ls()) 

#install packages
install.packages("stargazer")



# Load Libraries
library(ggfortify)
library(FactoMineR)
library(lattice)
library(vcd)
library(factoextra)
library(stats)
library(plyr)
library(dplyr)
library(ggcorrplot)
library(corrr)
library(zoom)            
library(geomorph)
library(abind) #combines multidimensional data frames like my tps files
library(ggplot2)
library(ggpubr)
library(plotrix)
library(vegan)
library(Morpho)
library(finalfit)
library(MASS)
library(rstan)
library(data.table)
library(sjPlot)
library(boot)
library(stargazer)
library(caret) 


########################## Preparation #################################

### Load  and prep all tps files
tps_1 <- readland.tps("warm.TPS", specID = "ID")
tps_2 <- readland.tps("cold.TPS", specID = "ID")


# Add ID where I forgot
#which(dimnames(all)[[3]] == "168")
#dimnames(all)[[3]][1:40] <- c(1:40) # if i want to rename the labels (if there are duplicate names) with numbers instead
#c(Zahl1:Zahl2) defines names for the renaming by dimnames; [1:40] = order and number of rows to be renamed
# Combine all tps fiels into one

all <- abind(tps_1,tps_2)

# save to csv
#write.csv(all,"C:/Users/nanda/Desktop/Morphometric photos/R code/test.csv")

# Remove space at the end
dimnames(all)
dimnames(all)[[3]] <- gsub(" ", "", dimnames(all)[[3]])
dimnames(all)

### Load and prep meta data
info <- read.delim2("CSHS1file.txt", head = TRUE, stringsAsFactors = TRUE, dec = ".")
info$SL <- (info$SL*10) # excel messed up the decimal numbers, so I had to convert from cm to mm for calculating BC_SL
info$BC_SL <- (((cs$weight*1000)/10)/cs$SL^3)

#write.csv(info, "info.csv", row.names = FALSE)

#compare missing data between 2 data sets; setdiff(larger set, smaller set)
#write.csv(all, "C:/Users/nanda/Desktop/Morphometric photos/R code/test.csv")
#IDS <- read.csv ("test.csv", header = TRUE)
#IDS <- unique(IDS)
#setdiff(info$names, IDS$ID)


# Substistute all "" for NA in landmarks to be remove column
info[info == ""] <- NA

#
boxplot (info$SL ~ info$sex:info$temp)
boxplot (info$SL ~ info$temp:info$sex)
boxplot (info$weight ~ info$temp:info$sex)
#


# Remove all id's which I couldn't landmark
info <- subset(info, !(ID %in% no_pics))

#
#convert values into dataframe, not needed but good to know
#all <- data.frame(all)
#


# Get list of fish id's with open mouth
open_mouth <- subset(all, mouth == "open")$ID

### Start with landmark stuff
# Remove landmarks 1-5 for fish with open mouth
all[c(1:2,28:29),,which(dimnames(all)[[3]] %in% "open mouth")] <- NA
#all[which(dimnames(all)[[3]] %in% "open mouth")] <- NA
head(all)
str(all)

#all[c(1:5),,which(dimnames(all)[[3]] %in% "open_mouth")] <- NA

# Remove all other landmarks that I couldn't place
for(i in 1:nrow(info)){
  if(is.na(info[i,"remove"]))
    next
  else{
    row_name <- which(dimnames(all)[[3]] == info[i,"ID"])
    for(l in 1:5){
      landmark <- strsplit(info[i,"remove"], split = ",")[[1]][l]
      if(is.na(landmark))
        next
      else
        all[as.numeric(landmark),,row_name] <- NA
    }
  }
}

### Estimate missing landmarks
all_final <- estimate.missing(all, method = "TPS") 

### Define sliding landmarks and save them
#sliders <- define.sliders(all, nsliders = 22) # not working, use a csv instead!
sliders <- read.csv("Sliders.csv", sep = ";",  head = TRUE, stringsAsFactors = F) 
#take seperators into account and adjust it with sep =
print (sliders)



########################## Head Shape Analyses #######################################

### Generalized Procrustes analysis
# Generalized Procrustean Analysis (GPA) is used in sensory data analysis prior to a 
# Preference Mapping to reduce the scale effects and to obtain a consensual configuration

gpa <- gpagen(all, curves = sliders)
plot(gpa)
plotOutliers(gpa$coords)
zm() # zoom in plot with zoom package, instructions below in console 

# Outliers starting with the one furthest away of upper quartile
outliers <- c("AS02664","AS02993","AS03261")

### Run analysis without outliers identified in procruster analysis
#info <- subset(info, !(ID %in% outliers))
all_final <- all[,,-which(dimnames(all)[[3]] %in% outliers)]


### Combine metadata with morph data
# Order metadata ID's by morph data
names <- dimnames(gpa$coords)[[3]]
info <- info[match(names, info$ID),]

#c(repX,Y) literally counts from row to row till the Y has been reach and labels the rows with X
# Build combined data frame with morph data, origin, length and diet
#all_info <- geomorph.data.frame(gpa, origin = info$origin, length = info$SL, fish_id = c(1:40), session = c(rep(1,20),rep(2,20)))
                  #comma in upper line and no ")" if there data needs to be included
all_info <- geomorph.data.frame(gpa, origin = info$origin, SL = info$SL, BC_SL = info$BC_SL, 
                                ID = info$ID, temp = info$temp, diet = info$F1_diet,
                                sex = info$sex, weight = info$weight)
                              #  adult_chironomidae = juv_diet$adult_chironomidae_percentage, 
                               # chironomidae = juv_diet$tot_chironomidae_percentage,
                                # small_CL = juv_diet$small_CL_percentage, 
                                # big_CL = juv_diet$big_CL_percentage,
                                # various = juv_diet$tot_various_percentage)

#inter observer
InterObserver <- procD.lm(coords ~ session, iter = 10000, SS.type = "III", data = all_info)
summary (InterObserver)
#it's also usable to compare between two treatment types
InterObserver <- procD.lm(coords ~ temp, iter = 10000, SS.type = "III", data = all_info)
summary (InterObserver)
#SS standard dev, MS, mean standard dev, PR high = unsignificant; session (as a treatment doesn't affect the data)

########################## Add Head Length #######################################

#linear_marks <- matrix(c(1,9,6,12,14,16), ncol=2, byrow=TRUE, dimnames = list(c("head_l","head_w", "eye_w"),c("start", "end")))
linear_marks <- matrix(c(1,20,4,47,23,25), ncol=2, byrow=TRUE, dimnames = list(c("head_l","head_w", "eye_w"),c("start", "end")))
length_measurements <- as.data.frame(interlmkdist(all, linear_marks))

# Order metadata ID's by morph data and add it to data set
names <- dimnames(all)[[3]]
info <- info[match(names, info$ID),]

#categorical data
length_measurements$ID <- info$ID

length_measurements$diet <- info$F1_diet
length_measurements$origin <- info$origin
length_measurements$sex <- info$sex
length_measurements$temp <- info$temp
#numeric data
length_measurements$SL <- info$SL
length_measurements$weight <- info$weight
length_measurements$BC_SL <- info$BC_SL

length_measurements$head_prop <- length_measurements$head_l/length_measurements$SL
length_measurements$eye_surf <- ((length_measurements$eye_w/2)^2)*3.14159265359

#define categorical data as factors (not IDs for obv. reasons)
length_measurements$origin <- factor(as.factor(length_measurements$origin), levels=c("hs","cs"))
length_measurements$diet <- factor(as.factor(length_measurements$diet), levels=c("bw","cl"))
length_measurements$sex <- factor(as.factor(length_measurements$sex), levels=c("m","f"))
length_measurements$temp <- factor(as.factor(length_measurements$temp), levels=c("warm","cold"))

#put all new defined variables into a new data frame
all_info_updated <- geomorph.data.frame(all_info, head_l = length_measurements$head_l, head_w = length_measurements$head_w,
                                        eye_w = length_measurements$eye_w, eye_surf = length_measurements$eye_surf, 
                                        head_prop = length_measurements$head_prop, SL = length_measurements$SL, 
                                        weight = length_measurements$weight, BC = length_measurements$BC_SL, 
                                        temp =length_measurements$temp)

#boxplots 
ggplot(length_measurements, aes(x = sex, y = SL)) +
  geom_boxplot()


########################## Linear Models ####################################### DOESN'T WORK
cor.test(all_info_updated$SL, all_info_updated$Csize) #Csize = average distance between landmark to middle; google it later 
                                                      # very highly correlated, one of these enough for linear model
cor.test(all_info_updated$SL, all_info_updated$head_l)# doesn't matter, SL used to calc. head_l
cor.test(all_info_updated$SL, all_info_updated$weight)# high correlation, don't use both in ln


cor.test(all_info_updated$head_l, all_info_updated$head_w)# very highly correlated, one of these enough for linear model
cor.test(all_info_updated$SL, all_info_updated$head_prop)# very highly (negatively) correlated...

#cor matrix
cor_matrix <- cor(data.frame(SL = info$SL, weight = info$weight, Csize = all_info_updated$Csize, 
                             head_l = all_info_updated$head_l, head_w = all_info_updated$head_w, 
                             eye_w = all_info_updated$eye_w, eye_surf = all_info_updated$eye_surf, 
                             BC = ((all_info_updated$weight*100)/(all_info_updated$SL*100^3))))

  ggcorrplot(cor_matrix)

#select specific columns of a data frame
  #df = select(df,c(column x, column y, column z)); c = column; e.g. column x: c(1)

#all correlated, don't use them in the same model, because they explain the same variation

#### Procrusts ANOVA origin & body length
# ANOVA including interaction: Interaction wasn't significant so I removed it
# Run length and origin separate as I already know they're correlated

#ANOVA with only 1 variable: use type I instead of III
  SL <- procD.lm(SL ~ sex * temp + diet, iter = 1000, SS.type = "III", data = all_info_updated)
  summary(SL)
  HP <- procD.lm(head_prop ~ sex * temp + diet, iter = 1000, SS.type = "III", data = all_info_updated) 
  summary (HP)
  

#linear model SL
lmSLO <- lm(SL ~ sex * temp + sex * diet + origin, length_measurements)
anova (lmSLO)
lmSL0<- aov(lm(SL ~ sex * temp + sex * diet + origin, length_measurements),
      SS.type='III')

lmSL1 <- lm(SL ~ sex * temp + diet, length_measurements)
summary (lmSL1)
lmSL2 <- lm(SL ~ sex * temp + sex * diet, length_measurements)#this one!
summary (lmSL2)
lmSL3 <- lm(SL ~ sex * diet + temp, length_measurements)
summary (lmSL3)
lmSL4 <- lm(SL ~ sex * temp * diet, length_measurements)
summary (lmSL4)

AIC(lmSL1,lmSL2,lmSL3,lmSL4, k = 2) # small aic value = better

lmHPO <- lm(head_prop ~ sex + temp + diet + origin, length_measurements)
anova (lmHPO)

lmHP1 <- lm(head_prop ~ sex + temp + diet, length_measurements)
summary(lmHP1)
lmHP2 <- lm(head_prop ~ sex * temp + diet, length_measurements)
summary(lmHP2)

AIC(lmHP1, lmHP2, k = 2) # small aic value = better

tab_model(lmSL2,
          digits = 4,
          show.df = TRUE,
          p.style = c("numeric_stars"),
          show.stat = TRUE,
          string.stat = "t")


tab_model(lmHP1,
          digits = 4,
          show.df = TRUE,
          p.style = c("numeric_stars"),
          string.stat = "t-value",
          show.stat = TRUE,
          string.stat = "t-value")



boxplot (length_measurements$head_prop ~ length_measurements$temp)  # head_prop --> larger head in comparison to the body length
boxplot (length_measurements$head_prop ~ length_measurements$sex:length_measurements$diet) #sticklebacks in poor nutrition environment: they invest more into head, so the body is proportionally smaller

ggplot(length_measurements, aes(x = diet, y = head_prop, col = temp)) + # model: head_prop ~ temp + sex + diet
  geom_boxplot() +
  facet_wrap(~sex) +
  xlab ("Temperature") + ylab ("Head length / Standard length") + 
  scale_color_manual(values = c("coral","blue"))
  

ggplot(length_measurements, aes(x = diet, y = SL, col = temp)) + # model: SL ~ temp + sex + diet
  geom_boxplot() +
  facet_wrap(~sex) +
  xlab ("Temperature") + ylab ("Standard length in mm") + 
  scale_color_manual(values = c("coral","blue"))

ggplot(length_measurements, aes(x = sex:diet, y = head_prop)) + # model: head_prop ~ temp + sex + diet
  geom_boxplot() +
  facet_wrap(~origin) +
  xlab ("Origin") + ylab ("Head length / Standard length") + 
  scale_color_manual(values = c("coral","blue"))

#test statistics - anova
anovaSL <- aov(SL ~ origin + temp * sex + diet * sex, data = length_measurements, SS.type = "III")
summary(anovaSL)

anovaHP <- aov(head_prop ~ origin + temp + sex + diet, data = length_measurements)
summary(anovaHP)


#  leaving all interactions out, since the interactions only had an effect size less than 1%
# so this might be the model we want to keep working on
lm_inter <- procD.lm(coords ~ sex + temp + diet, iter = 1000, SS.type = "II", data = all_info_updated) # type III since the interactions don't seem to be important
summary(lm_inter) # diet doesn't contribute that much into variations (Rsq quite low in comparison)

    boxplot (all_info_updated$SL ~ all_info_updated$sex:all_info_updated$temp)

### Principal and phylogenetically-aligned components analysis of shape data
pc <- gm.prcomp(all_info$coords)
plot(pc)
summary(pc) # Comp 1-3 contribute the most into the variance of the results (Comp 33% = 0.335 proportion of variance);
# cumulative variance = sum of all comp.variance contribution so far

########################## Pretty Plots #######################################

###################### PCA 1-3 in data frame ######################
### Get pc values into data frame
#pca_data <- data.frame(pc1 = pc$x[,1], pc2 = pc$x[,2], pc3 = pc$x[,3], origin = info$origin, ID = info$ID, SL = info$SL,
                       #BC_SL = info$BC_SL, diet = info$F1_diet, temp = info$temp, sex = info$sex, Csize = all_info_updated$Csize)
pca_data <- data.frame(origin = info$origin, diet = info$F1_diet, temp = info$temp, sex = info$sex, 
                       pc1 = pc$x[,1], pc2 = pc$x[,2], pc3 = pc$x[,3], 
                       Csize = all_info_updated$Csize, SL = info$SL)
head(pca_data) #what's inside my pca data now?
str(pca_data)

boxplot (all_info_updated$head_prop ~ all_info_updated$temp:all_info_updated$sex)
boxplot (all_info_updated$head_prop ~ all_info_updated$sex:all_info_updated$diet)

boxplot (all_info_updated$head_prop ~ all_info_updated$sex)
boxplot (all_info_updated$head_prop ~ all_info_updated$temp)
#
boxplot (all_info_updated$SL ~ all_info_updated$temp:all_info_updated$sex)
boxplot (all_info_updated$SL ~ all_info_updated$sex)
boxplot (all_info_updated$SL ~ all_info_updated$temp)


# Factorial Analysis of Mixed Data 
#pca_famd <- FAMD(pca_data, graph = FALSE)
#pca_mfa <- MFA(pca_data)
#pca_famd$eig
#pca_famd$quali.var
#df = select(df,c(column x, column y, column z)); c = column; e.g. column x: c(1)

#
#fviz_famd_ind(pca_famd, col.ind = "cos2",
#             gradient.cols = c("blue", "orange", "coral"),
#             repel = FALSE) # repel true = no overlaps between labels

### Create PC1 vs PC2 (pure and mean) Plots; distance on x-axis more meaningful than for the y-axis (PC1 contributes more into var)
# Adding in a zoomed in version of the mean plot

#temp, splitted by sex
PC1_2 <-   ggplot(pca_data, aes(x = pc1, y = pc2, col = temp)) + #pc1 and pc2 are the x and y coords from the coords
  geom_point() +
  scale_color_manual(values = c("blue","coral", "black")) +
  xlab("PC1 33.5%") + ylab("PC2 16.9%") +
  ggtitle("PC1 vs PC2") +
  stat_ellipse()+ facet_wrap(~sex:diet)# circles = means of each treatments; facet_wrap (~variable) splits it into categories
PC1_2 

    #### Calculate means and store them in a data frame   
pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                           mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                           mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                           sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                           sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                           sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                           temp = c("Cold", "Cold", "Warm", "Warm"),
                           sex = c("F", "M", "F", "M"))
    
    PC1_2_mean <- ggplot(pca_means_sd, aes(x = mean1, y = mean2, col = temp)) + 
      geom_point(size = 3) +
    #  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + # red dot for means but not necessary here
      scale_color_manual(values = c("blue","coral")) +
      xlab("PC1 33.5%") + ylab("PC2 16.9%") +
      ggtitle("Mean PC1 vs Mean PC2") +
      geom_errorbar(aes( ymin = mean2 - sd2, ymax = mean2 + sd2)) +
      geom_errorbar(aes(xmin = mean1 - sd1, xmax = mean1 + sd1)) +  facet_wrap(~ sex)
    PC1_2_mean
    #    #
      pg <- ggplot_build(PC1_2_mean) # get values from the plots 
      pg
    #theme(text = element_text(size = 25)) +
  
  #temp, splitted by diet and sex
  PC1_2 <-   ggplot(pca_data, aes(x = pc1, y = pc2, col = temp)) + 
    geom_point() +
    scale_color_manual(values = c("blue","coral")) +
    xlab("PC1%") + ylab("PC2") +
    ggtitle("PC1 vs PC2") +
    stat_ellipse()+ facet_wrap(~sex + diet)# circles = means of each treatments; facet_wrap (~variable) splits it into categories
  PC1_2 


#sex
PC1_2 <-   ggplot(pca_data, aes(x = pc1, y = pc2, col = sex)) + 
  geom_point() +
  scale_color_manual(values = c("blue","coral")) +
  xlab("PC1") + ylab("PC2") +
  ggtitle("PC1 vs PC2") +
  stat_ellipse()+ facet_wrap(~temp)# circles = means of each treatments; facet_wrap (~variable) splits it into categories
                                   # temp (warm) reduces variance and guides coordinates (morph) to a smaller size
PC1_2
    #### Calculate means and store them in a data frame   
pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                           mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                           mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                           sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                           sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                           sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                           temp = c("Cold", "Cold", "Warm", "Warm"),
                           sex = c("F", "M", "F", "M"))
    
    PC1_2_mean <- ggplot(pca_means_sd, aes(x = mean1, y = mean2, col = temp)) + 
      geom_point(size = 3) +
      #  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + # red dot for means but not necessary here
      scale_color_manual(values = c("blue","coral")) +
      xlab("PC1 33.5%") + ylab("PC2 16.9%") +
      ggtitle("Mean PC1 vs Mean PC2") +
      geom_errorbar(aes( ymin = mean2 - sd2, ymax = mean2 + sd2)) +
      geom_errorbar(aes(xmin = mean1 - sd1, xmax = mean1 + sd1)) 
    PC1_2_mean
    #
    pg <- ggplot_build(PC1_2_mean) # get values from the plots 
    pg
    #theme(text = element_text(size = 25)) +

#diet
PC1_2 <-   ggplot(pca_data, aes(x = pc1, y = pc2, col = diet)) + 
  geom_point() +
  scale_color_manual(values = c("blue","coral")) +
  xlab("PC1 33.5%") + ylab("PC2 16.9%") +
  ggtitle("PC1 vs PC2") +
  stat_ellipse()+ facet_wrap(~sex)# circles = means of each treatments; facet_wrap (~variable) splits it into categories
PC1_2 
    #### Calculate means and store them in a data frame   
pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                           mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                           mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                           sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                           sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                           sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                           temp = c("Cold", "Cold", "Warm", "Warm"),
                           sex = c("F", "M", "F", "M"))

    PC1_2_mean <- ggplot(pca_means_sd, aes(x = mean1, y = mean2, col = temp)) + 
      geom_point(size = 3) +
      #  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + # red dot for means but not necessary here
      scale_color_manual(values = c("blue","coral")) +
      xlab("PC1 33.5%") + ylab("PC2 16.9%") +
      ggtitle("Mean PC1 vs Mean PC2") +
      geom_errorbar(aes( ymin = mean2 - sd2, ymax = mean2 + sd2)) +
      geom_errorbar(aes(xmin = mean1 - sd1, xmax = mean1 + sd1)) 
    PC1_2_mean
      #
      pg <- ggplot_build(PC1_2_mean) # get values from the plots 
      pg


#origin
PC1_2 <-   ggplot(pca_data, aes(x = pc1, y = pc2, col = origin)) + 
  geom_point() +
  scale_color_manual(values = c("blue","coral")) +
  xlab("PC1 33.5%") + ylab("PC2 16.9%") +
  ggtitle("PC1 vs PC2") +
  stat_ellipse()+ facet_wrap(~diet)# circles = means of each treatments; facet_wrap (~variable) splits it into categories
PC1_2 
    #### Calculate means and store them in a data frame   
pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                           mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                           mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                           sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                           sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                           sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                           temp = c("Cold", "Cold", "Warm", "Warm"),
                           sex = c("F", "M", "F", "M"))
    
    PC1_2_mean <- ggplot(pca_means_sd, aes(x = mean1, y = mean2, col = temp)) + 
      geom_point(size = 3) +
      #  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + # red dot for means but not necessary here
      scale_color_manual(values = c("blue","coral")) +
      xlab("PC1 33.5%") + ylab("PC2 10.7%") +
      ggtitle("Mean PC1 vs Mean PC2") +
      geom_errorbar(aes( ymin = mean2 - sd2, ymax = mean2 + sd2)) +
      geom_errorbar(aes(xmin = mean1 - sd1, xmax = mean1 + sd1)) 
    PC1_2_mean
        #
        pg <- ggplot_build(PC1_2_mean) # get values from the plots 
        pg    

pc
#theme(text = element_text(size = 25))
        


#
#
### Create PC2 vs PC3 (pure and mean) Plots
#temp
PC2_3 <- ggplot(pca_data[,], aes(x = pc2, y = pc3, col = temp)) + 
  geom_point() +
  scale_color_manual(values = c("blue","coral")) +
  #scale_color_manual(values = c("black","coral", "coral","black","black")) + #north vs south basin
  xlab("PC2 16.9%") + ylab("PC3 10.7%") +
  ggtitle("PC2 vs PC3") +
  stat_ellipse()+ facet_wrap(~ sex:diet)
PC2_3 
#theme(text = element_text(size = 25)) 
    #### Calculate means and store them in a data frame   
pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                           mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                           mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                           sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                           sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                           sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                           temp = c("Cold", "Cold", "Warm", "Warm"),
                           sex = c("F", "M", "F", "M"))
    
    PC2_3_mean <- ggplot(pca_means_sd, aes(x = mean2, y = mean3, col = temp)) + 
      geom_point(size = 3) +
      scale_color_manual(values = c("blue","coral")) +
      xlab("PC2 16.9%") + ylab("PC3 10.7%") +
      ggtitle("Mean PC2 vs Mean PC3") +
      geom_errorbar(aes( ymin = mean3 - sd3, ymax = mean3 + sd3)) +
      geom_errorbar(aes(xmin = mean2 - sd2, xmax = mean2 + sd2))
    PC2_3_mean
        #
        pg <- ggplot_build(PC2_3_mean) # get values from the plots 
        pg    
#theme(text = element_text(size = 25))

#sex
        PC2_3 <- ggplot(pca_data[,], aes(x = pc2, y = pc3, col = sex)) + 
          geom_point() +
          scale_color_manual(values = c("blue","coral")) +
          #scale_color_manual(values = c("black","coral", "coral","black","black")) + #north vs south basin
          xlab("PC2 16.9%") + ylab("PC3 10.7%") +
          ggtitle("PC2 vs PC3") +
          stat_ellipse()+ facet_wrap(~ temp)
        PC2_3 
        #theme(text = element_text(size = 25)) 
        #### Calculate means and store them in a data frame   
        pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                             tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                             tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                             tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                                   mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                             tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                             tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                             tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                                   mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                             tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                             tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                             tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                                   sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                           tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                           tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                           tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                                   sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                           tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                           tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                           tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                                   sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                           tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                           tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                           tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                                   temp = c("Cold", "Cold", "Warm", "Warm"),
                                   sex = c("F", "M", "F", "M"))
        
        PC2_3_mean <- ggplot(pca_means_sd, aes(x = mean2, y = mean3, col = temp)) + 
          geom_point(size = 3) +
          scale_color_manual(values = c("blue","coral")) +
          xlab("PC2 16.9%") + ylab("PC3 10.7%") +
          ggtitle("Mean PC2 vs Mean PC3") +
          geom_errorbar(aes( ymin = mean3 - sd3, ymax = mean3 + sd3)) +
          geom_errorbar(aes(xmin = mean2 - sd2, xmax = mean2 + sd2))
        PC2_3_mean
        #
        pg <- ggplot_build(PC2_3_mean) # get values from the plots 
        pg    

#diet
PC2_3 <- ggplot(pca_data[,], aes(x = pc2, y = pc3, col = diet)) + 
  geom_point() +
  scale_color_manual(values = c("blue","coral")) +
  #scale_color_manual(values = c("black","coral", "coral","black","black")) + #north vs south basin
  xlab("PC2 16.9%") + ylab("PC3 10.7%") +
  ggtitle("PC2 vs PC3") +
  stat_ellipse()+ facet_wrap(~ temp)
PC2_3 
#theme(text = element_text(size = 25)) 
      #### Calculate means and store them in a data frame   
pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                           mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                           mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                           sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                           sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                           sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                           temp = c("Cold", "Cold", "Warm", "Warm"),
                           sex = c("F", "M", "F", "M"))

      PC2_3_mean <- ggplot(pca_means_sd, aes(x = mean2, y = mean3, col = temp)) + 
        geom_point(size = 3) +
      scale_color_manual(values = c("blue","coral")) +
      xlab("PC2 16.9%") + ylab("PC3 10.7%") +
      ggtitle("Mean PC2 vs Mean PC3") +
      geom_errorbar(aes( ymin = mean3 - sd3, ymax = mean3 + sd3)) +
      geom_errorbar(aes(xmin = mean2 - sd2, xmax = mean2 + sd2))
    PC2_3_mean
      #
      pg <- ggplot_build(PC2_3_mean) # get values from the plots 
      pg    
    #theme(text = element_text(size = 25))
#origin
      PC2_3 <- ggplot(pca_data[,], aes(x = pc2, y = pc3, col = diet)) + 
        geom_point() +
        scale_color_manual(values = c("blue","coral")) +
        #scale_color_manual(values = c("black","coral", "coral","black","black")) + #north vs south basin
        xlab("PC2 16.9%") + ylab("PC3 10.7%") +
        ggtitle("PC2 vs PC3") +
        stat_ellipse()+ facet_wrap(~ sex)
      PC2_3 
      #theme(text = element_text(size = 25)) 
          #### Calculate means and store them in a data frame   
      pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                           tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                           tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                           tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                                 mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                           tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                           tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                           tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                                 mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                           tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                           tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                           tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                                 sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                         tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                         tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                         tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                                 sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                         tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                         tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                         tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                                 sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                         tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                         tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                         tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                                 temp = c("Cold", "Cold", "Warm", "Warm"),
                                 sex = c("F", "M", "F", "M"))
      
      PC2_3_mean <- ggplot(pca_means_sd, aes(x = mean2, y = mean3, col = temp)) + 
        geom_point(size = 3) +
        scale_color_manual(values = c("blue","coral")) +
        xlab("PC2 16.9%") + ylab("PC3 10.7%") +
        ggtitle("Mean PC2 vs Mean PC3") +
        geom_errorbar(aes( ymin = mean3 - sd3, ymax = mean3 + sd3)) +
        geom_errorbar(aes(xmin = mean2 - sd2, xmax = mean2 + sd2))
      PC2_3_mean
              #
              pg <- ggplot_build(PC2_3_mean) # get values from the plots 
              pg    
### Create PC1 vs PC3 (pure and mean) Plots
#temp    
PC1_3 <- ggplot(pca_means_sd, aes(x = mean1, y = mean2, col = temp)) + 
  geom_point(size = 3) +
  scale_color_manual(values = c("blue","coral")) +
  #scale_color_manual(values = c("black","coral", "coral","black","black")) + #north vs south basin
  xlab("PC1 33.5%") + ylab("PC3 10.7%%") +
  ggtitle("PC1 vs PC3") +
  stat_ellipse() +  facet_wrap(~ sex)
PC1_3
    #### Calculate means and store them in a data frame   
pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                           mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                           mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                     tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                           sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                           sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                           sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                   tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                           temp = c("Cold", "Cold", "Warm", "Warm"),
                           sex = c("F", "M", "F", "M"))

    PC1_3_mean <- ggplot(pca_means_sd, aes(x = mean1, y = mean3, col = temp)) + 
      geom_point(size = 3) +
      scale_color_manual(values = c("blue","coral")) +
      xlab("PC1 33.5%") + ylab("PC3 10.7%%") +
      ggtitle("Mean PC1 vs Mean PC3") +
      geom_errorbar(aes( ymin = mean3 - sd3, ymax = mean3 + sd3)) +
      geom_errorbar(aes(xmin = mean1 - sd1, xmax = mean1 + sd1))
    PC1_3_mean
          #
          pg <- ggplot_build(PC1_3_mean) # get values from the plots 
          pg    

          
          #diet    
          PC1_3 <- ggplot(pca_data[,], aes(x = pc1, y = pc3, col = diet)) + 
            geom_point() +
            scale_color_manual(values = c("blue","coral")) +
            #scale_color_manual(values = c("black","coral", "coral","black","black")) + #north vs south basin
            xlab("PC1 33.5%") + ylab("PC3 10.7%%") +
            ggtitle("PC1 vs PC3") +
            stat_ellipse() +  facet_wrap(~ sex)
          PC1_3
          #### Calculate means and store them in a data frame   
          pca_means_sd <- data.frame(mean1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],mean)[1],
                                               tapply(pc$x[,1], info[, c("temp", "sex")],mean)[3],
                                               tapply(pc$x[,1], info[, c("temp", "sex")],mean)[2],
                                               tapply(pc$x[,1], info[, c("temp", "sex")],mean)[4]),
                                     mean2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],mean)[1],
                                               tapply(pc$x[,2], info[, c("temp", "sex")],mean)[3],
                                               tapply(pc$x[,2], info[, c("temp", "sex")],mean)[2],
                                               tapply(pc$x[,2], info[, c("temp", "sex")],mean)[4]),
                                     mean3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],mean)[1],
                                               tapply(pc$x[,3], info[, c("temp", "sex")],mean)[3],
                                               tapply(pc$x[,3], info[, c("temp", "sex")],mean)[2],
                                               tapply(pc$x[,3], info[, c("temp", "sex")],mean)[4]),
                                     sd1 = c(tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[1],
                                             tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[3],
                                             tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[2],
                                             tapply(pc$x[,1], info[, c("temp", "sex")],std.error)[4]),
                                     sd2 = c(tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[1],
                                             tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[3],
                                             tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[2],
                                             tapply(pc$x[,2], info[, c("temp", "sex")],std.error)[4]),
                                     sd3 = c(tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[1],
                                             tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[3],
                                             tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[2],
                                             tapply(pc$x[,3], info[, c("temp", "sex")],std.error)[4]),
                                     temp = c("Cold", "Cold", "Warm", "Warm"),
                                     sex = c("F", "M", "F", "M"))
          
          PC1_3_mean <- ggplot(pca_means_sd, aes(x = mean1, y = mean3, col = temp)) + 
            geom_point(size = 3) +
            scale_color_manual(values = c("blue","coral")) +
            xlab("PC1 33.5%") + ylab("PC3 10.7%%") +
            ggtitle("Mean PC1 vs Mean PC3") +
            geom_errorbar(aes( ymin = mean3 - sd3, ymax = mean3 + sd3)) +
            geom_errorbar(aes(xmin = mean1 - sd1, xmax = mean1 + sd1))
          PC1_3_mean
          #
          pg <- ggplot_build(PC1_3_mean) # get values from the plots 
          pg    
          

### Plot PCA's and mean PCA's
ggarrange(PC1_2, PC1_2_mean, PC2_3, PC2_3_mean, PC1_3, PC1_3_mean,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)

ggarrange(PC1_2, PC2_3,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
ggarrange(PC1_2,
          labels = c("A"),
          ncol = 1, nrow = 1)

###################### Plot Mean body shape & extremes
mean_morph <- mshape(gpa$coords)
par(mfrow = c(2, 2))# generally allows to plot smaller versions of plots next to each other, also boxplots 



#shape predictor for Aless' data
plot_links <- matrix(c(1,2,3,1,1,2,5,6,9,9,11,12,19,20,21,13,14,15,16,17,6,7,6,10,
                       2,3,19,5,4,4,6,9,10,11,12,22,20,21,22,14,15,16,17,13,7,8,10,11),
                     ncol = 2) #can also be done with define.links()
# my data body
plot_links <- matrix(c(1,3,4,5,6,7,8,9,11,12,13,15,16,47,37,29,28,4,20,17,19,22,23,24,25,
                       3,4,5,6,7,8,9,11,12,13,15,16,47,37,29,28,1,20,17,19,47,23,24,25,22),
                     ncol = 2) #can also be done with define.links()

pc1_pred <- shape.predictor(gpa$coords, pc$x[,1], Intercept = F, pred1 = min(pc$x[,1]), pred2 = max(pc$x[,1]))
plotRefToTarget(mean_morph, pc1_pred$pred1, links = plot_links, mag = 1.1)#lower end of pc1
plotRefToTarget(mean_morph, pc1_pred$pred2, links = plot_links, mag = 1.1)#upper end of pc1


pc2_pred <- shape.predictor(gpa$coords, pc$x[,2], Intercept = F, pred1 = min(pc$x[,2]), pred2 = max(pc$x[,2]))
plotRefToTarget(mean_morph, pc2_pred$pred1, links = plot_links, mag = 1.1)#lower end of pc2
plotRefToTarget(mean_morph, pc2_pred$pred2, links = plot_links, mag = 1.1)#upper end of pc2


###################### Plot Mean body shape against length and diet
mean_morph <- mshape(gpa$coords)
par(mfrow = c(1, 2))# generally allows to plot smaller versions of plots next to each other, also boxplots 

size_pred <- shape.predictor(gpa$coords, info$SL, Intercept = T, pred1 = min(info$SL), pred2 = max(info$SL))
plotRefToTarget(mean_morph, size_pred$pred1,links = plot_links, mag = 1)
plotRefToTarget(mean_morph, size_pred$pred2,links = plot_links, mag = 1)


# my data head
par(mfrow = c(2, 1)) # generally allows to plot smaller versions of plots next to each other, also boxplots 

head_links <- matrix(c(1,3,4,20,17,19,47,37,29,22,23,24,25,1,2,28,
                       3,4,20,17,19,47,37,29,1,23,24,25,22,2,28,1),
              ncol = 2) #can also be done with define.links()

head_size_pred <- shape.predictor(gpa$coords, all_info_updated$head_l, Intercept = T, pred1 = min(all_info_updated$head_l), pred2 = max(all_info_updated$head_l))
plotRefToTarget(mean_morph, size_pred$pred1,links = head_links, mag = 4)
plotRefToTarget(mean_morph, size_pred$pred2,links = head_links, mag = 4)

pc1_head <- shape.predictor(gpa$coords, pc$x[,1], Intercept = F, pred1 = min(pc$x[,1]), pred2 = max(pc$x[,1]))
plotRefToTarget(mean_morph, pc1_head$pred1, links = head_links, mag = 4)#lower end of pc1
plotRefToTarget(mean_morph, pc1_head$pred2, links = head_links, mag = 4)#upper end of pc1


pc2_head <- shape.predictor(gpa$coords, pc$x[,2], Intercept = F, pred1 = min(pc$x[,2]), pred2 = max(pc$x[,2]))
plotRefToTarget(mean_morph, pc2_head$pred1, links = head_links, mag = 0.2)#lower end of pc2
plotRefToTarget(mean_morph, pc2_head$pred2, links = head_links, mag = 0.2)#upper end of pc2

###################### Plot PCA's by SL
ggplot(pca_data, aes(x = SL, y = pc1 , col = SL)) + 
  geom_point(size = 2) +
  xlab("SL") + ylab("PC1") +
  ggtitle("PC1 vs SL") 
#SL
ggplot(pca_data, aes(x = SL, y = pc1, col = temp)) + 
  geom_point(size = 2) +
  xlab("SL") + ylab("PC1") +
  ggtitle("PC1 vs SL") +
  scale_color_manual(values = c("blue","coral")) +
  stat_ellipse()+ facet_wrap(~sex) #screw origin lmao
 #
