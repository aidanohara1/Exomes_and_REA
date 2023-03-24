library(Rmisc)
library(dplyr)

# read in
AGiles <- read.csv("AGiles.csv")
AGiles_OG <- read.csv("AGiles_originalOnly.csv")
# original reports are already isolated in Agiles_OG


# col names for ref
colnames(AGiles)[11] = "abcBuckets"
colnames(AGiles)[10] = "defBuckets"
colnames(AGiles_OG)[10] = "abcBuckets"
colnames(AGiles_OG)[9] = "defBuckets"


# correct positive case
AGiles$`Report.Outcome`[AGiles$`Report.Outcome` == "positive"] <- "Positive"

# factorize report outcome
AGiles_OG <- mutate(AGiles_OG, classification = ifelse(`Report.Outcome`=="Positive",1,
                                                       ifelse(`Report.Outcome`=="Uncertain",0,-1)))
# factorize report outcome
AGiles <- mutate(AGiles, classification = ifelse(`Report.Outcome`=="Positive",1,
                                                       ifelse(`Report.Outcome`=="Uncertain",0,-1)))

# depricated#
#AGiles <- mutate(AGiles, abcFactor = as.factor(abcBuckets))

write.csv(AGiles,"rEAExomes.csv")


# AGiles with only single latest exome entries
latestAgiles <- group_by(AGiles, `New.study.ID`) %>%
  slice(which.max(as.Date(`Report.Date`, '%M/%Y'))) %>%
  mutate(binReanaz = ifelse(`Original.Latest` != "Original", 1,0))

# not sure what this does
#removeOG <- filter(AGiles, `Initiator.of.reanalysis` != "N/A")


# tally of id entries per id, mainly for helper
idEntries <- group_by(AGiles, `New.study.ID`) %>%
  tally()

# sort out by reanalysis
# caveat... Should we run this model s.t. if an exome is reclassified its inital entry is ommitted...

# the agiles Ill use for analysis, 
# most important note is the ommision of Orgianal reports about exomes that have been reanalyzed.
# keeping multiple reanalysis still... There are tradeoffs
omitAgiles <- inner_join(AGiles, idEntries, by = "New.study.ID") %>%
  filter(!(n>1 & `Original.Latest` == "Original"))


# tiny eda
mycount <- count(AGiles_OG, `ABC.REA.buckets`)
barplot(sort(mycount$n))




#ignore/omit N/As?
### or ####
# treat N/A as its own classification?


# ANOVA
#https://webpower.psychstat.org/wiki/_media/grant/mai-zhang-2017.pdf

# NOT ALLOWED


# set up the null hypothesis
# the mean, ratio, rate? 
# of reclassification, and reanalysis of exomes of each REA classification
# as well as the mean classification outcome
# is the same

# we will reject the null hypothesis if its more than 1 mean off?

# use the f- stat statistic test

# Our decision criteria will be alpha = 0.05 

# calculate the degrees of freedom

# Compute the f stat

# MS between 
# divided by
# MS within groups

# models, mostly just sum of squares totals
reclassAnBad <- lm(data = omitAgiles, formula = Reclassified ~ abcBuckets)
reclassAn <- glm(data = omitAgiles, formula = Reclassified ~ abcBuckets, family = binomial)

reanazAn <- lm(data = omitAgiles, formula = Reanalyed ~ abcBuckets)

initClassAn <- lm(data = AGiles_OG, formula = classification ~ abcBuckets)
initClassAn
latestClassAn <- lm(data = latestAgiles, formula = classification ~ abcBuckets)
latestClassAn

oneReclass <- glm(data = latestAgiles, formula = Reanalyed ~ abcBuckets, family = binomial)

# contingency <- latestAgiles %>%
#   group_by(abcBuckets,Reanalyed) %>%
#   summarize()


install.packages("gmodels")
library(gmodels)
CrossTable(latestAgiles$abcBuckets, latestAgiles$Reanalyed, chisq = TRUE)

anova(reclassAn, test = 'F')
round(qf(0.05,7,11154,lower.tail = FALSE),2)

anova(latestClassAn, test = 'F')
round(qf(0.05,7,10913,lower.tail = FALSE),2)
latestClassAn

reclassAn
anova(reclassAn, test = 'Chisq')
round(qchisq(0.05,7,lower.tail = FALSE),2)

# you get an F
#passMyFTest <- function(){}

anova(oneReclass, test = 'Chisq')

### RECLASSSIFICATION plot ###
reclassMeans <- group_by(omitAgiles, abcBuckets) %>%
  dplyr::summarise(meanReclassification = mean(Reclassified),
            uci_perGroup = CI(Reclassified)['lower'],
            lci_perGroup = CI(Reclassified)['upper']) %>%
  arrange(abcBuckets)

reanazMeans <- group_by(omitAgiles, abcBuckets) %>%
  dplyr::summarise(meanReanaz = mean(Reanalyed),
                   uci_perGroup = CI(Reanalyed)['lower'],
                   lci_perGroup = CI(Reanalyed)['upper']) %>%
  arrange(abcBuckets)

classMeans <- group_by(omitAgiles, abcBuckets) %>%
  dplyr::summarise(meanClassification = mean(classification),
                   uci_perGroup = CI(classification)['lower'],
                   lci_perGroup = CI(classification)['upper']) %>%
  arrange(abcBuckets)

latClassMeans <- group_by(latestAgiles, abcBuckets) %>%
  dplyr::summarise(meanClassification = mean(classification),
                   uci_perGroup = CI(classification)['lower'],
                   lci_perGroup = CI(classification)['upper']) %>%
  arrange(abcBuckets)


library(ggplot2)
reclassGG <- ggplot(reclassMeans) +
  geom_bar( aes(x=abcBuckets, y=meanReclassification), stat = "identity", alpha = 0.05) +
  geom_errorbar( aes(x=abcBuckets, ymin=lci_perGroup, ymax=uci_perGroup, color = abcBuckets), alpha=0.9, size=1)+
  coord_flip()+
  theme_dark()

reanazGG <- ggplot(reanazMeans) +
  geom_bar( aes(x=abcBuckets, y=meanReanaz), stat = "identity", alpha = 0.05) +
  geom_errorbar( aes(x=abcBuckets, ymin=lci_perGroup, ymax=uci_perGroup, color = abcBuckets), alpha=0.9, size=1)+
  coord_flip()+
  theme_dark()

classGG <- ggplot(classMeans) +
  geom_bar( aes(x=abcBuckets, y=meanClassification), stat = "identity", alpha = 0.05) +
  geom_errorbar( aes(x=abcBuckets, ymin=lci_perGroup, ymax=uci_perGroup, color = abcBuckets), alpha=0.9, size=1)+
  ggtitle("Cumulative classifications") +
  coord_flip()+
  theme_dark()

latClassGG <- ggplot(latClassMeans) +
  geom_bar( aes(x=abcBuckets, y=meanClassification), stat = "identity", alpha = 0.05) +
  geom_errorbar( aes(x=abcBuckets, ymin=lci_perGroup, ymax=uci_perGroup, color = abcBuckets), alpha=0.9, size=1)+
  ggtitle("Up to Date classifications") +
  coord_flip()+
  theme_dark()

reclassGG
reanazGG
classGG
latClassGG



# TODO this week
 # pick a question that andrew giles wants answered that lets us ignore
 #  the problem with having multiple entries/reanalysis
 #  After that construct the model, including other pertinent covariates
 #  Try and be ready to work with dan by wednesday morning
 #  He'll present our model and chi-sq test to the profs
 #   PRESENT the updated results on Thursday

write.csv(latestAgiles, "latestAgiles.csv")

library(gmodels)

latest <- read.csv("latestAgiles.csv")

# Does the distribution of receiving at least one reanalysis differ by REA
oneReclass <- glm(data = latestAgiles, formula = Reanalyed ~ abcBuckets, family = binomial)

anova(oneReclass, test = 'Chisq')

CrossTable(latestAgiles$abcBuckets, latestAgiles$Reanalyed, chisq = TRUE)




# library(ggplot2)
# 
# AandB <- filter(AGiles, `ABC.REA.buckets` == "A" | `ABC.REA.buckets` == "B")
# Others <- filter(AGiles, !(`ABC.REA.buckets` == "A" | `ABC.REA.buckets` == "B"))
# 
# ggplot(AGiles, aes(fill=`Report.Outcome`, x=`Report.Outcome`)) + 
#   geom_histogram(stat = "count") +
#   ggtitle(" ") +
#   facet_wrap(~`ABC.REA.buckets`) +
#   theme_ipsum() +
#   theme(legend.position="none") +
#   xlab("")
# 


#https://www.statology.org/anova-unequal-sample-size/

#how many exomes are under each classification?

#Kruskal–Wallis one-way analysis of variance
#https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance



#https://www.youtube.com/watch?v=ZWoiWzn08tw


