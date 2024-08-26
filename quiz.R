###############################################################
#
# How is the identification success of volunteers and can it be improved through practice?
#
###############################################################

#####
#----clear data
setwd("Insert path of working directory")
Data <- read.table("DatenQuizzz_ber.txt", header = T, sep = "\t")

names(Data)
str(Data)

Data$Bild <- as.factor(Data$Bild)
Data$WB <- as.factor(Data$WB)
Data$ID <- as.factor(Data$ID)
Data$Pseudonym <- as.factor(Data$Pseudonym)
Data$Age <- as.factor(Data$Age)
Data$Sex <- as.factor(Data$Sex)
Data$informed <- as.factor(Data$informed)
Data$Taxon <- as.factor(Data$Taxon)
Data$Antagonist <- as.factor(Data$Antagonist)
Data$Taxon_corr <- as.factor(Data$Taxon_corr)
Data$Antagonist_corr <- as.factor(Data$Antagonist_corr)
Data$PersID <- as.factor(Data$PersID)
Data$week <- as.factor(Data$week)

head(Data)

# Descriptive statistic ---------------------------------------------------
#Calculate success rate of identification
for (i in 1:nrow(Data)) {
  if (Data$Taxon[i] == Data$Taxon_corr[i]) {
    Data$successTax[i] <- 0
  } else{
    Data$successTax[i] <- 1
  }
}
(nrow(subset(Data, successTax == 1)) / nrow(Data)) * 100

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)
library(forcats)

#-----Display rate of successful identification along time gradient (error rate1-error rate10) per person (fill=PersID)
successrate7 <-
  Data %>% group_by(Pseudonym, week, successTax)  %>% tally()

successrate7 <- successrate7 %>%
  group_by(Pseudonym, week) %>%
  mutate(sum_by_IDW = sum(n)) %>%
  mutate(proz = n / sum_by_IDW * 100)

successrate7 <- successrate7[successrate7$successTax != 1,]
successrate7$week <- as.integer(successrate7$week)
summary(successrate7$Pseudonym)

plot8 <-
  ggplot(successrate7, aes(x = week, y = proz)) +
  geom_point(show.legend = FALSE) +
  geom_smooth(
    method = lm ,
    formula=y ~ poly(x, 2),
    color = "red",
    fill = "#69b3a2",
    se = TRUE
  ) +
  #geom_path(group=successrate7$Pseudonym, show.legend = FALSE)+
  ylab("success rate of\n identification [%]") +
  xlab("week") +
  theme_classic () +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 20))
plot8


# Model------------------------------------------------------------------
library(dplyr)
library(Matrix)
library(car)
library(glmmTMB)
library(performance)

# Model week ------------------------------------------------------------

#Extract essential columns into new data set and throw out NA's
Data1 <- Data [, c(1, 6:9, 17)]
Data2 <- na.omit(Data1, successTax)
Data2 <- na.omit(Data2, week)
Data2 <- na.omit(Data2, Age)
Data2 <- na.omit(Data2, Sex)
Data2 <- na.omit(Data2, informed)
Data2 <- na.omit(Data2, PersID)

options(na.action = "na.omit")
M1 <-
  glmmTMB(successTax ~  week + Age + Sex + informed + (1 |
                                                         PersID),
          data = Data2,
          family = "binomial")
summary(M1)
step(M1)
Anova(M1) #Age und Sex are not significant


M2 <-
  glmmTMB(successTax ~  week + informed + (1 |
                                             PersID),
          data = Data2,
          family = "binomial")
#M3<- glmmTMB(successTax ~  week * informed + (1|PersID), data= Data, family= "binomial") #AIC 2116 --> interaction did not increase the AIC
Anova(M2)
check_overdispersion(M2)
check_model(M2)
fixef(M2)

#R?
library(pscl)
#calculate McFadden's R-squared for model
pR2(M2)['McFadden'] #0.11

library(report)
report(M2)

library(multcompView)
library(multcomp)
library(emmeans)

#for week:
CIs = cld(emmeans(M2, ~ week),
          type = "response",
          sort = TRUE,
          Letters = letters)
CIs$.group = gsub(" ", "", CIs$.group, fixed = TRUE)
CIs

ggplot(CIs, aes(x = week, y = prob * 100)) +
  geom_point(shape = 16, size = 2) +
  geom_errorbar(aes(ymin = lower.CL * 100, ymax = upper.CL * 100), width =
                  0.1) +
  geom_text(aes(label = .group, y = upper.CL * 100 + 5)) +
  ylab("success rate of identification (%)") +
  xlab("week") +
  theme_bw() +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15))


# Model vaidation -------------------------------------------------------
library(mlbench)
###60/40 Split
splitmodel <- Data2

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
  glmmTMB(successTax ~  week + informed + (1 |
                                             PersID),
          data = train,
          family = "binomial")
p <-
  predict(model,  test, type = "response", allow.new.levels = TRUE)

p_class <- ifelse(p > 0.5, "1", "0")
table(p_class)

#2 way frequency table
table(p_class, test[["successTax"]])

library(caret)
confusionMatrix(table(p_class, test[["successTax"]]))

library(pROC)
par(pty = "s")
roc(test$successTax, p, plot = TRUE)#AUC: 0.72