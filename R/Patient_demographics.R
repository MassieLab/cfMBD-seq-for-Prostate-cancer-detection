setwd("~/OneDrive - University of Cambridge/Projects/MBD/")

install.packages("tidyverse")
library(data.table)
library(tidyverse)

MBD_patients <- fread("MBD_Patient_info.csv")

##compare age and PSA of benign vs (localised) vs metastatic patients in either test or training set
MBD_patients_train <- MBD_patients %>% filter(.,Test_or_train=="train")
MBD_patients_test <- MBD_patients %>% filter(.,Test_or_train=="test")

summary(aov(Age~Stage,data=MBD_patients_train))
summary(aov(PSA~Stage,data=MBD_patients_train))
summary(aov(Age~Stage,data=MBD_patients_test))
summary(aov(PSA~Stage,data=MBD_patients_test))

#Since there are 3 groups for the test set (ben, loc, mets), posthoc analysis:
TukeyHSD(aov(Age~Stage,data=MBD_patients_test))
TukeyHSD(aov(PSA~Stage,data=MBD_patients_test))

##compare age and PSA of test vs training set
MBD_patients_ben <- MBD_patients %>% filter(.,Stage=="benign")
MBD_patients_met <- MBD_patients %>% filter(.,Stage=="metastatic")

summary(aov(Age~Test_or_train,data=MBD_patients_ben))
summary(aov(PSA~Test_or_train,data=MBD_patients_ben))
summary(aov(Age~Test_or_train,data=MBD_patients_met))
summary(aov(PSA~Test_or_train,data=MBD_patients_met))



