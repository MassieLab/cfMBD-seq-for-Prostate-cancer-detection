library(GenomicRanges)
library(tidyverse)
library(data.table)
library(scales) 
library(ggsci)
library(RColorBrewer)
library("caret")
library(ggplot2)     
library(grid)
library(gridExtra)   
library(likert)

################################
########### Figures ############
################################

load("ICHOR_merged_DMRs_30_170123.RData")
load("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/Dataset/ICHOR_merged_DMRs_30_170123.RData")

dmrs_final_merged <- dmrs_final %>% group_by(seqnames,start,end) %>%dplyr::select(ends_with("_beta")) %>% summarise_all("mean",na.rm = TRUE)

dmrs_final %>% group_by(seqnames,start,end) %>%dplyr::select(ends_with("_beta"),ends_with("_beta_means")) %>% summarise_all("mean",na.rm = TRUE)

means <- dmrs_final %>% group_by(seqnames,start,end) %>% dplyr::select(ends_with("_beta_means")) %>% summarise_all("mean",na.rm = TRUE) %>% mutate(diff = BENIGN_beta_means-MET_beta_means,abs_diff=abs(diff))

mean_dmrs <- dmrs_final %>% group_by(seqnames,start,end) %>%dplyr::select(ends_with("_beta"),ends_with("_beta_means")) %>% summarise_all("mean",na.rm = TRUE) %>% mutate(diff = BENIGN_beta_means-MET_beta_means,abs_diff=abs(diff))

dmrs <- dmrs_final %>% group_by(seqnames,start,end) %>% dplyr::select(meta_train$regions,meta_test$regions,ends_with("_beta_means")) %>% summarise_all("mean",na.rm = TRUE) %>% mutate(diff = BENIGN_beta_means-MET_beta_means,abs_diff=abs(diff))


mean_dmrs %>% slice_max(n=10,order_by = abs_diff)
meta_test <- read.csv("test_samples.csv")
meta_test <- meta_test %>% filter(group != "CPG1" & group != "CPG5") %>% mutate(regions = case_when(str_detect(sample_name,"CMDL",negate = TRUE) ~ paste0("X",sample_name,"_beta"),
                                                                                                    str_detect(sample_name,"CMDL") ~ paste0(sample_name,"_beta")))

meta_test_CPG1 <- meta_test %>% filter(group == "CPG1") %>% mutate(regions = case_when(str_detect(sample_name,"CMDL",negate = TRUE) ~ paste0("X",sample_name,"_beta"),
                                                                                       str_detect(sample_name,"CMDL") ~ paste0(sample_name,"_beta")))

meta_test_CPG5 <- meta_test %>% filter(group == "CPG5") %>% mutate(regions = case_when(str_detect(sample_name,"CMDL",negate = TRUE) ~ paste0("X",sample_name,"_beta"),
                                                                                       str_detect(sample_name,"CMDL") ~ paste0(sample_name,"_beta")))

meta_train <- read.csv("train_samples.csv")
meta_train <- meta_train %>% mutate(regions = case_when(str_detect(sample_name,"CMDL",negate = TRUE) ~ paste0("X",sample_name,"_beta"),str_detect(sample_name,"CMDL") ~ paste0(sample_name,"_beta"))) %>% filter(regions != "X667_beta")

mutate_regions <- full_cohort_overlaps %>% mutate(regions = tolower(paste0("chr",seqnames, "_", start))) %>% 
  dplyr::select(regions, ends_with("_beta")) %>% dplyr::select(-c("X667_beta")) %>% na.omit()

DMRs <- read.csv("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/auc_GB_loop_ichor.csv",header = FALSE)

colnames(DMRs) <- c("n","fpr","tpr","auc","y_test","y_pred_test","y_pred")

as.data.frame(cbind(as.numeric(unlist(stringr::str_extract_all(DMRs$fpr, '\\d+([.,]\\d+)?'))),as.numeric(unlist(stringr::str_extract_all(DMRs$tpr, '\\d+([.,]\\d+)?'))),as.numeric(unlist(stringr::str_extract_all(DMRs$y_test, '\\d+([.,]\\d+)?'))),as.numeric(unlist(stringr::str_extract_all(DMRs$y_pred_test, '\\d+([.,]\\d+)?'))),as.numeric(unlist(stringr::str_extract_all(DMRs$y_pred, '\\d+([.,]\\d+)?')))))
colnames(c) <- c("fpr","tpr","group")

predictions <- as.data.frame(cbind(as.numeric(unlist(stringr::str_extract_all(DMRs$y_test, '\\d+([.,]\\d+\\D+\\d+)?'))),as.numeric(unlist(stringr::str_extract_all(DMRs$y_pred_test, '\\d+([.,]\\d+\\D+\\d+)?'))),as.numeric(unlist(stringr::str_extract_all(DMRs$y_pred, '\\d+([.,]\\d+\\D+\\d+)?')))))
colnames(predictions) <- c("y_test","y_pred_test","y_pred")

ichor_CM0035 <- read.delim("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/ICHOR/ichor_output_summary_CM0035.txt",sep = " ")
ichor_CM0017 <- read.delim("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/ICHOR/ichor_output_summary_CM0017.txt",sep = " ")

ichor_CM0035$sample_name <- gsub(".*_","",gsub(".r_1_val_1_grch38_bwa_sorted_dedup.bam.wig_ichor","",ichor_CM0035$Sample))
ichor_CM0017$sample_name <- str_extract(ichor_CM0017$Sample, "CMDL\\d{8}")

ichor <- rbind(ichor_CM0017,ichor_CM0035)

ichor_test <- right_join(ichor, meta_test, by = "sample_name")


ichor_pred <- cbind(predictions,ichor_test)

ichor_pred %>% ggplot(aes(x = Tumor_Fraction, y=y_pred, colour = y_test)) + geom_point()
log

test <- read.csv("/Users/tm739/Dropbox (Cambridge University)/Mac/Downloads/test_for_ML_withnames.csv")

ichor_pred <- cbind(test[,c(1,902)],predictions) %>% left_join(ichor_test,by = "regions")

cpg <- read.csv("/Users/tm739/Dropbox (Cambridge University)/Mac/Downloads/auc_GB_CPG15.csv",header = FALSE)

colnames(cpg) <- c("n","fpr","tpr","auc","y_test","y_pred_test","y_pred")

predictions_cpg1 <- as.data.frame(cbind(as.numeric(unlist(stringr::str_extract_all(cpg$y_pred_test[2], '\\d+([.,]\\d+)?'))),rep("CPG1",length(as.numeric(unlist(stringr::str_extract_all(cpg$y_pred_test[2], '\\d+([.,]\\d+)?')))))))
predictions_cpg5 <- as.data.frame(cbind(as.numeric(unlist(stringr::str_extract_all(cpg$y_pred_test[3], '\\d+([.,]\\d+)?'))),rep("CPG5",length(as.numeric(unlist(stringr::str_extract_all(cpg$y_pred_test[3], '\\d+([.,]\\d+)?')))))))

estimated_cpg1 <- ifelse(predictions_cpg1$V1>0.86,1,0)
estimated_cpg5 <- ifelse(predictions_cpg5$V1>0.86,1,0)
confusionMatrix(data = as.factor(estimated_cpg1), reference = factor(rep(1,20), levels = c(0,1)))
confusionMatrix(data = as.factor(estimated_cpg5), reference = factor(rep(1,20), levels = c(0,1)))


as.data.frame(cbind(predictions_cpg1$V1,predictions_cpg5$V1)) 

cpg1 <- as.data.frame(cbind(test_for_ML_CPG1_g$regions,predictions_cpg1$V1)) 
colnames(cpg1) <- c("regions","pred_score")

cpg5 <- as.data.frame(cbind(test_for_ML_CPG5_g$regions,predictions_cpg5$V1)) 
colnames(cpg5) <- c("regions","pred_score")

df_cpg1 <- cbind(cpg1,rep("CPG1",20))
df_cpg5 <- cbind(cpg5,rep("CPG5",20))
colnames(df_cpg1) <- c("regions", "pred_score","group")
colnames(df_cpg5) <- c("regions", "pred_score","group")
df <- ichor_pred %>% select(regions,y_pred,group.x)
colnames(df) <- c("regions", "pred_score","group")

df <- rbind(df,df_cpg1,df_cpg5)

df$pred_score <- as.numeric(df$pred_score)

df$group2 <- ifelse(df$group == "CPG1" | df$group == "CPG5","Localised",df$group)

cpg1_ben <- df %>% filter(group == "BENIGN" | group == "CPG1")
cpg5_ben <- df %>% filter(group == "BENIGN" | group == "CPG5")
met_ben <- df %>% filter(group == "BENIGN" | group == "MET")

local_ben <- df %>% filter(group == "BENIGN" | group2 == "Localised")

cpg1_ROC <- do.call(rbind.data.frame, roc(cpg1_ben$pred_score,as.factor(ifelse(cpg1_ben$group == "BENIGN",0,1)))) %>% 
  t() %>% as.data.frame() 

cpg5_ROC <- do.call(rbind.data.frame, roc(cpg5_ben$pred_score,as.factor(ifelse(cpg5_ben$group == "BENIGN",0,1)))) %>% 
  t() %>% as.data.frame() 

local_ROC <- do.call(rbind.data.frame, roc(local_ben$pred_score,as.factor(ifelse(local_ben$group == "BENIGN",0,1)))) %>% 
  t() %>% as.data.frame() 

met_ROC <- do.call(rbind.data.frame, roc(met_ben$pred_score,as.factor(ifelse(met_ben$group == "BENIGN",0,1)))) %>% 
  t() %>% as.data.frame() 

AUC_cpg1 <- auc(roc(cpg1_ben$pred_score,as.factor(ifelse(cpg1_ben$group == "BENIGN",0,1))))
AUC_cpg5 <- auc(roc(cpg5_ben$pred_score,as.factor(ifelse(cpg5_ben$group == "BENIGN",0,1))))
AUC_met <- auc(roc(met_ben$pred_score,as.factor(ifelse(met_ben$group == "BENIGN",0,1))))
AUC_local <- auc(roc(local_ben$pred_score,as.factor(ifelse(local_ben$group == "BENIGN",0,1))))


################################
########## Figure 3 ############
################################


################################
########## Figure 4 ############
################################

# ICHOR Prediction plot

ichor_pred %>% ggplot(aes(x = Tumor_Fraction, y=y_pred, colour = group.x)) + 
  geom_point(size = 4)  + #scale_x_continuous(trans=log10_trans()) + 
  xlab("ICHOR CNA tumour fraction") + ylab("Prediction probability") +
  geom_hline(yintercept = 0.86, linetype = "dashed") +
  geom_vline(xintercept = 0.15,linetype = "dashed") +
  annotate(geom="text", x=0.50, y= 0.80, label="Prediction Cutoff" , color="black", size = 8) +
  annotate(geom="text", x=0.18, y= 0.4, label="ICHOR tumour threshold", size = 8, angle = 90) +
  scale_colour_manual(values=group.colors) +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 20),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/ICHORvsPrediction.png",width = 8.6, height = 7, units = "in", dpi = 300)
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/ICHORvsPrediction.svg",width = 8.6, height = 7, units = "in", dpi = 300)

# Confusion Matrix
library(caret)
library("likert")
estimated <- ifelse(ichor_pred$y_pred>0.86,1,0)

cm <- confusionMatrix(data = as.factor(estimated), reference = as.factor(ichor_pred$y_test),positive = "1")
cm_d <- as.data.frame(cm$table) # extract the confusion matrix values as data.frame
cm_st <-data.frame(cm$overall) # confusion matrix statistics as data.frame
cm_st$cm.overall <- round(cm_st$cm.overall,2) # round the values
cm_d$diag <- cm_d$Prediction == cm_d$Reference # Get the Diagonal
cm_d$ndiag <- cm_d$Prediction != cm_d$Reference # Off Diagonal     
#cm_d[cm_d == 0] <- NA # Replace 0 with NA for white tiles
cm_d$Reference <-  reverse.levels(cm_d$Reference) # diagonal starts at top left
cm_d$ref_freq <- cm_d$Freq * ifelse(is.na(cm_d$diag),-1,1)

ggplot(data = cm_d, aes(x = Prediction , y =  Reference, fill = Freq))+
  geom_tile( data = cm_d,aes(fill = ref_freq)) +
  scale_fill_gradient2(guide = FALSE ,low="red3",high="lightblue", midpoint = 0,na.value = 'white') +
  geom_text(aes(label = Freq), color = 'black', size = 10) +
  scale_x_discrete(labels=c("0" = "Control", "1" = "Case"),position = "top") +
  scale_y_discrete(labels=c("0" = "Control", "1" = "Case")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 25),
        panel.border = element_blank(),
        plot.background = element_blank(),
        axis.line = element_blank(),
        axis.text=element_text(color="black", size = 25),
        axis.ticks = element_line(colour = "black", size = 1)
  )
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/ConfusinMatrix_ICHOR.svg", dpi = 300)

do.call(rbind.data.frame, AUC::sensitivity(ichor_pred$y_pred,as.factor(ichor_pred$y_test),perc.rank = FALSE)) %>% 
  t() %>% as.data.frame() %>% ggplot(aes(x= V1, y = V2)) + geom_line() +
  geom_vline(xintercept = 0.86,linetype = "dashed") +
  xlab("Cutoffs") + ylab("Sensitivity") +
  annotate(geom="text", x=0.82, y= 0.30, label="Prediction Cutoff" , color="black", size = 10, angle = 90) +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 30),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 

ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/Sensitivity.png", dpi = 300)
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/Sensitivity.svg", dpi = 300)

do.call(rbind.data.frame, AUC::specificity(ichor_pred$y_pred,as.factor(ichor_pred$y_test),perc.rank =FALSE)) %>% 
  t() %>% as.data.frame() %>% ggplot(aes(x= V1, y = V2)) + geom_line() +
  geom_vline(xintercept = 0.86,linetype = "dashed") +
  xlab("Cutoffs") + ylab("Specificity") +
  annotate(geom="text", x=0.82, y= 0.30, label="Prediction Cutoff" , color="black", size = 10, angle = 90) +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 30),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/Specificity.png", dpi = 300)
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/Specificity.svg", dpi = 300)


ggplot(met_ROC,aes(x= V2, y = V3)) + geom_line(size = 1, color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = "dotted", size = 1) +
  annotate(geom="text", x=0.90, y= 0.2, label=paste0("AUC = ",round(AUC_met,2)) , color="black", size = 8) +
  geom_rect(aes(xmin = 0.75, xmax = 1.05 , ymin = 0.15, ymax = 0.25),
            fill = "transparent", color = "black", size = 0.2) +
  xlab("False Positive Rate") + ylab("True Postive Rate") +
  scale_colour_npg() +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        axis.line = element_line(colour = 'black', size = 1.2),
        axis.ticks = element_line(colour = "black", size = 1.2),
        text = element_text(size = 20),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 

ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_AUROC.png", dpi = 300)


################################
########## Figure 5 ############
################################



group.colors.box <- c(BENIGN = "black", CPG1 = "yellow",CPG5 = "orange" , MET = "red")
df %>% ggplot(aes(x = group, y = pred_score, fill = group , group = group)) + 
  geom_boxplot(width=0.5,alpha = 0.5) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5,fill = "black", color = "black") +
  geom_hline(yintercept = 0.86, linetype = "dashed") + 
  xlab("Stage") + ylab("Prediction probability") +
  scale_fill_npg() +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"),
        text = element_text(size = 20), 
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_boxplots.png", dpi = 300)

my_comparisons <- list( c("BENIGN", "CPG1"), c("BENIGN", "CPG5"), c("BENIGN", "MET") )

df %>% ggplot(aes(x = group, y = pred_score, fill = group , group = group)) + 
  geom_boxplot(width=0.5,alpha = 0.5) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5,fill = "black", color = "black") +
  geom_hline(yintercept = 0.86, linetype = "dashed") + 
  xlab("Stage") + ylab("Prediction probability") +
  scale_fill_npg() +
  theme_classic() + 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.5) +
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"),
        text = element_text(size = 20), 
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_boxplots_withsignif.png", dpi = 300)
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_boxplots_withsignif.svg", dpi = 300)


group.colors.box <- c(BENIGN = "black", Localised = "yellow" , MET = "red")
my_comparisons <- list( c("BENIGN", "Localised"), c("Localised", "MET"), c("BENIGN", "MET") )
df <- df %>% mutate(stage = str_to_title(group2))

df %>% ggplot(aes(x = stage, y = pred_score, fill = stage , group = stage)) + 
  geom_boxplot(width=0.5,alpha = 0.5) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5,fill = "black", color = "black") +
  geom_hline(yintercept = 0.86, linetype = "dashed", size = 1) + 
  xlab("Stage") + ylab("Prediction probability") +
  scale_fill_npg() +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 20),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 

ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_boxplots_localised.png",width = 8.6, height = 7, units = "in", dpi = 300)


my_comparisons <- list( c("BENIGN", "Localised"), c("Localised", "MET"), c("BENIGN", "MET") )

df %>% ggplot(aes(x = group2, y = pred_score, fill = group2 , group = group2)) + 
  geom_boxplot(width=0.5,alpha = 0.5) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5,fill = "black", color = "black") +
  geom_hline(yintercept = 0.86, linetype = "dashed") + 
  xlab("Stage") + ylab("Prediction probability") +
  # scale_colour_manual(values=group.colors.box, aesthetics = "fill") +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #  stat_compare_means(aes(group = group2), label = "p.signif", label.y = 0.50,label.x = 1.5,method="wilcox.test", paired=F)+
  stat_compare_means(label.y = 1.5) +
  scale_fill_npg() +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 20),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 

ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_boxplots_localised_withsignif.png", dpi = 300)
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_boxplots_localised_withsignif.svg", dpi = 300)


a <- roc(df$pred_score,as.factor(ifelse(df$group == "BENIGN",0,1)))
AUC_all <- auc(a)


do.call(rbind.data.frame, AUC::roc(df$pred_score,as.factor(ifelse(df$group == "BENIGN",0,1)))) %>% 
  t() %>% as.data.frame() %>% ggplot(aes(x= V2, y = V3)) + geom_line() + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = "dotted") +
  annotate(geom="text", x=0.90, y=0.10, label=paste0("AUC = ",round(AUC_all,2)), color="black") +
  xlab("False Positive Rate") + ylab("True Postive Rate") +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 20),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_AUROC.png", dpi = 300)

show_col(pal_npg(palette = c("nrc"), alpha = 1)(3))

do.call(rbind.data.frame, AUC::roc(df$pred_score,as.factor(ifelse(df$group == "BENIGN",0,1)))) %>% 
  t() %>% as.data.frame() %>% ggplot(aes(x= V2, y = V3)) + geom_line( colour = "#E64B35FF",alpha = 1) +
  geom_line(data = local_ROC, colour = "#4DBBD5FF", alpha = 1) +
  geom_line(data = met_ROC, colour = "#00A087FF") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = "dotted",size = 1) +
  annotate(geom="text", x=0.90, y= 0.09, label=paste0("AUC Local= ",round(AUC_local,2)) , color="#4DBBD5FF", size = 8) +
  annotate(geom="text", x=0.90, y= 0.15, label=paste0("AUC All = ",round(AUC_all,2)) , color="#E64B35FF", size = 8) +
  annotate(geom="text", x=0.90, y= 0.21, label=paste0("AUC Met = ",round(AUC_met,2)) , color="#00A087FF", size = 8) +
  geom_rect(aes(xmin = 0.70, xmax = 1.1 , ymin = 0.05, ymax = 0.25),
            fill = "transparent", color = "black", size = 0.2) +
  xlab("False Positive Rate") + ylab("True Postive Rate") +
  theme_classic() + 
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 20),
        legend.position = "top", 
        legend.title = element_blank(),
        legend.justification = "center",
        axis.text.x=element_text(),
        axis.text.y = element_text()) 

ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_AUROC_local.png",width = 8.6, height = 7, units = "in", dpi = 300)
ggsave("/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/Ermira/Machine Learning/ICHOR_ML_V1/All_sample_AUROC_local.svg",width = 8.6, height = 7, units = "in", dpi = 300)
