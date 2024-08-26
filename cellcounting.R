###############################################################
#
#Comparing cell counting data of volunteers and experts
#
###############################################################


setwd("Insert path of working directory")
Data <- read.table("Bestimmungsdaten_4.txt", header = T, sep = "\t")

head(Data)

# clear data --------------------------------------------------------

Data$fID <- as.factor(Data$ID)
summary(Data$fID)
nlevels(Data$fID)

Data$fD <- as.factor(Data$D)

Data$fTaxE <- as.factor(Data$TaxE)
Data$fTaxP <- as.factor(Data$TaxP)

Data$fGegTaxE <- as.factor(Data$GegTaxE)
Data$fGegTaxP <- as.factor(Data$GegTaxP)

Data$fSex <- as.factor(Data$Sex)
Data$fPersID <- as.factor(Data$PersID)
Data$fAge <- as.factor(Data$Age)
Data$fLevelKnow <- as.factor(Data$LevelKnow)
Data$fBackground <- as.factor(Data$Background)


# Delete data with missing values
Data2 <- subset(Data, fTaxP != "Leer")

# Delete data that could not be identified by experts
levels(Data2$fTaxP)
Data3 <- subset(Data2, fTaxP != "Unbestimmbar")
levels(factor(Data3$fTaxP))

Data4 <- Data3
Data4$fTaxP <- factor(Data4$fTaxP)
Data4$fTaxE <- factor(Data4$fTaxE)
summary(Data4$fTaxP)
summary(Data4$fTaxE)

# Summarise taxa
Data4[Data4$fTaxP == "Passaloecus eremita", ]$fTaxP <-
  "Passaloecus spp."
Data4[Data4$fTaxE == "Passaloecus eremita", ]$fTaxE <-
  "Passaloecus spp."

# Only taxa that occur > 10 times should remain in the data, as other data are not representative enough
Data4$fTaxP <- as.character(Data4$fTaxP)
Data4$fTaxE <- as.character(Data4$fTaxE)

Data4 <- subset(Data4, fTaxP != "Dipogon spp.")
Data4 <- subset(Data4, fTaxP != "Auplopus carbonarius")
Data4 <- subset(Data4, fTaxP != "Hoplitis claviventris")
Data4 <- subset(Data4, fTaxP != "Colletes daviesanus")
Data4 <- subset(Data4, fTaxP != "Pseudoanthidium nanum")
Data4 <- subset(Data4, fTaxP != "Megachile spp.")
Data4 <- subset(Data4, fTaxP != "Osmia brevicornis")

Data4 <- subset(Data4, fTaxE != "Auplopus carbonarius")
Data4 <- subset(Data4, fTaxE != "Hoplitis claviventris")
Data4 <- subset(Data4, fTaxE != "Pseudoanthidium nanum")
Data4 <- subset(Data4, fTaxE != "Colletes daviesanus")
Data4 <- subset(Data4, fTaxE != "Megachile spp.")
Data4 <- subset(Data4, fTaxE != "Osmia brevicornis")

Data4$fTaxP <- as.factor(Data4$fTaxP)
Data4$fTaxE <- as.factor(Data4$fTaxE)
summary(Data4$fTaxP)
summary(Data4$fTaxE)


# "Osmia" is a taxon identified on another level --> delete it
Data5 <- subset(Data4, fTaxE != "Osmia")


#Delete cavities which the volunteers have indicated as empty/only 0 or max. 1 started cell was present
Data5$er <- rep(0, nrow(Data5))
for (i in 1:nrow(Data5)) {
  if ((Data5$fTaxE[i] == "Leer") &&
      ((Data5$NZP[i] = 1) || (Data5$NZP[i] = 0))) {
    Data5$er[i] <- 1
  } else{
    Data5$er[i] <- 0
  }
}
Data5 <- Data5[!(Data5$er == 1) ,]


# Add required columns --------------------------------------------
# Summarise cavities duplicates
Data5$duplicate <- rep(0, nrow(Data5))
for (i in 1:(nrow(Data5) - 1)) {
  if ((Data5$ID[i] == Data5$ID[i + 1]) &&
      (Data5$BNr[i] == Data5$BNr[i + 1]) &&
      (Data5$NNr[i] == Data5$NNr[i + 1])) {
    Data5$duplicate[i + 1] <- 1
  } else{
    Data5$duplicate[i + 1] <- 0
  }
}
# Has at least one larva developed successfully?
for (i in 1:(nrow(Data5))) {
  if ((Data5$NerfZP[i] >= 1) || is.na(Data5$NerfZP[i])) {
    Data5$larvcontr[i] <- 1
  } else{
    Data5$larvcontr[i] <- 0
  }
}
#Calculate success rate of identification (1 = correct, 0 = incorrect compared to the result of experts)
for (i in 1:nrow(Data5)) {
  if (Data5$fTaxP[i] == Data5$fTaxE[i]) {
    Data5$successTax[i] <- 1
  } else{
    Data5$successTax[i] <- 0
  }
}
(nrow(subset(Data5, successTax == 1)) / nrow(Data5)) * 100

#Calculate successful counting (1 = correct, 0 = uncorrect compares to the result of experts)
for (i in 1:nrow(Data5)) {
  if (is.na(Data5$NerfZE[i]) ||
      is.na(Data5$NerfZP[i])) {
    Data5$successNerfZ[i] <- NA
  }
  else if (abs(Data5$NerfZE[i] - Data5$NerfZP[i]) <= 1) {
    Data5$successNerfZ[i] <- 1
  } else{
    Data5$successNerfZ[i] <- 0
  }
}
Data5 <- na.omit(Data5, successNerfZ)

# Total result of counting--------------------------------------------------------
library(dplyr)
summary(Data5) #3404
Dat <- Data5 %>% group_by(Data5$successTax)  %>% tally()
head(Dat) # 1 = 3215, 0= 189
3215 / 3404 * 100 #94,4%


# Model ------------------------------------------------------------------

#correlation?
pairs(Data5[c("fSex", "fAge", "fLevelKnow", "fBackground")])
# Background and Levelknow(Level of knowledge) have more or less the same information, only one is needed (Level of knowledge)


library(Matrix)
library(car)
library(glmmTMB)
library(performance)

#----Null model
M0 <-
  glmmTMB (successNerfZ ~ 1 + (1 |
                                 fPersID),
           data = Data5,
           family = "binomial")

#----Model
options(na.action = "na.omit")
M3 <-
  glmmTMB(
    successNerfZ ~  fTaxP + fD + larvcontr + fLevelKnow  + fSex + fAge + duplicate + (1 |
                                                                                        fPersID),
    data = Data5,
    family = "binomial"
  )
summary(M3)
Anova(M3)
step(M3)
#larvcontr is not significant

M4 <-
  glmmTMB(
    successNerfZ ~  fTaxP + fD + fLevelKnow  + fSex + fAge + duplicate + (1 |
                                                                            fPersID),
    data = Data5,
    family = "binomial"
  )
summary(M4)
Anova(M4)

check_model(M4) # Collinarity of Age, Levelknow und TaxP

#logical: first delete "Age"
M5 <-
  glmmTMB(
    successNerfZ ~  fTaxP + fD + fLevelKnow  + fSex + duplicate + (1 |
                                                                     fPersID),
    data = Data5,
    family = "binomial"
  )
summary(M5)
Anova(M5) # sex is not significant

M6 <-
  glmmTMB(
    successNerfZ ~  fTaxP + fD + fLevelKnow + duplicate + (1 |
                                                             fPersID),
    data = Data5,
    family = "binomial"
  )
summary(M6)
Anova(M6)
check_model(M6)
fixef(M6)
check_overdispersion(M6)

library(report)
report(M6)

#R?
library(pscl)
#calculate McFadden's R-squared for model
pR2(M6)['McFadden'] #0.13

library(multcompView)
library(multcomp)
library(emmeans)

#for Taxa
CIs = cld(emmeans(M6, c("fTaxP")),
          sort = FALSE,
          Letters = letters,
          type = "response")
CIs$.group = gsub(" ", "", CIs$.group, fixed = TRUE)
CIs

library(dbplyr)
library(forcats)
library(dplyr)
library(ggplot2)

CIs %>%
  mutate(
    fTaxP = fct_relevel(
      fTaxP,
      "Heriades spp.",
      "Trypoxylon spp.",
      "Chelostoma florisomne",
      "Chelostoma spp.",
      "Psenulus spp.",
      "Passaloecus spp.",
      "Hylaeus spp.",
      "Deuteragenia spp.",
      "Osmia bicornis/cornuta",
      "Osmia caerulescens",
      "Eumeninae",
      "Pemphredon spp.",
      "Hoplitis adunca"
    )
  ) %>%
  ggplot(aes(x = fTaxP, y = prob * 100)) +
  geom_point(shape = 16, size = 2) +
  geom_errorbar(aes(ymin = lower.CL * 100, ymax = upper.CL * 100), width =
                  0.1) +
  geom_text(aes(label = .group, y = upper.CL * 100 + 5)) +
  ylab("cell counting success rate (%)") +
  xlab("taxon") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(size=15, color="black"), 
        axis.title.x = element_text(size= 15),
        axis.title.y = element_text(size= 15), 
        axis.text.y = element_text(size=15, color="black",face="italic"))



# Model validation -------------------------------------------------------
library(mlbench)

###60/40 Split
splitmodel <- Data5

#store results in rows and randomly reorder the rows
set.seed(2108)
rows <- sample(nrow(splitmodel))
splitmodel <- splitmodel[rows,]

#identify proper row to split and store row numbers as split
split <- round(nrow(splitmodel) * 0.60)
split
#create train
train <- splitmodel[1:split,]
#save the last 40% as the test set
test <- splitmodel[(split + 1):nrow(splitmodel),]

#confirm size
nrow(train) / nrow(splitmodel)

#Confusion Matrix (predicted matrix vs reference)
model <-
  glmmTMB(
    successNerfZ ~ fTaxP + fD + fLevelKnow + duplicate + (1 |
                                                            fPersID),
    data = train,
    family = "binomial"
  )
p <- predict(model,  test, type = "response")

p_class <- ifelse(p > 0.5, "1", "0")
table(p_class)

#2 way frequency table
table(p_class, test[["successTax"]])

library(caret)
confusionMatrix(table(p_class, test[["successTax"]]))

library(pROC)
par(pty = "s")
roc(test$successTax, p, plot = TRUE) #AUC: 0.576
