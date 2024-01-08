#PLR analysis of PHOSP Olink Data 
# Felicity Liew 01/08/2023

#Libraries---- 
library(janitor)
library(tidyverse)
library(stringr)
library(glmnet)
library(writexl)
library(caret)
library(lme4)
library(broom)
library(clipr)
library(magrittr)
library('mlbench')
library(grid)
library(gridExtra)
library(data.table)
library(gridExtra)
library(grid)
library(scales)
library(ggpubr)
library(ipflasso)


#Read in data -----




#################################### Penalized logistic regression - Fatigue as example ###############################

olink_f = olink3_3m%>%
  filter(!is.na(Fatigue)) 

summary(olink_f$Fatigue)
olink_f = olink_f%>%
  mutate(any_comorb = NA)
olink_f$any_comorb[olink_f$neuro_comorb=='Yes'|olink_f$cp_comorb=='Yes'|olink_f$gi_comorb=='Yes'| 
                     olink_f$gen_fat_comorb=='Yes']='Yes'
olink_f$any_comorb[olink_f$comorb== 'None']='No'

olink_f$any_comorb=factor(olink_f$any_comorb, levels = c('No', 'Yes'))
olink_f = olink_f%>%
  filter(!is.na(any_comorb))

olink_f=olink_f%>%
  select(-cluster, -chronic_inf,- VASH1, - SULT2A1, - CCL28, -'HLA-DRA', - CTSC, 
         -BCL2L11, - RAB6A, -IL2RB, -IL33, -  ITGA6, -LILRB4 )%>% #remove variables with high missingness that were not sig in volcano plots
  select(-study_id, -PASC_4,-PASC_1, -PASC_2, -PASC_3,  -phosp_id, -SampleID, -neuro_comorb, -gi_comorb, -cp_comorb, 
         -comorb, -crf1a_resp_support_4levels, -psq_recovered, -gen_fat_comorb, -GI, -Neuro, -CP, -gad7_summary_2levels,
         -mocal_total_corrected_summary, -phq9_summary_2levels , -brain_fog, -Neuro,  -GI ,-affective, -aching_in_your_muscles_pai.x,
         -bpi_severity_summary, -pain )%>% 
  filter(!is.na(age_admission))%>%
  filter(!is.na(WHO_sev))%>%
  filter(!is.na(Fatigue))%>%
  filter(!is.na(WAS))%>%
  filter(!is.na(IL1A))%>%
  filter(!is.na(CCL22))


#split data to test and train
set.seed(123)
olink_f$crf1a_sex=ifelse(olink_f$crf1a_sex == "Female", 1, 0)
olink_f$WHO_sev=ifelse(olink_f$WHO_sev == "severe", 1, 0)

x = model.matrix(Fatigue~.,olink_f)[,-1]

#created outcome into numerical variable
y = ifelse(olink_f$Fatigue=="Fatigue", 1, 0)
cv.ridge=cv.glmnet(x,y, alpha = 0, family ="binomial")
plot(cv.ridge)
cv.ridge$lambda.min
coef(cv.ridge, cv.ridge$lambda.1se) 

model=glmnet(x,y, alpha = 0, family = "binomial", lamda = NULL)
## nested cv 
nestedcv_fat =cvr.glmnet(x,y, family = 'binomial', standardize = FALSE, nfolds = 5, ncv = 10, type.measure ='class' )
nestedcv_fat_auc =cvr.glmnet(x,y, family = 'binomial', standardize = FALSE, nfolds = 5, ncv = 10, type.measure ='auc' )
summary(nestedcv_fat$cvm)
summary(nestedcv_fat_auc$cvm)
fat_cvm = as.data.frame(nestedcv_fat$cvm)
fat_cvm_auc = as.data.frame(nestedcv_fat_auc$cvm)


#lasso method: Find the best lamda using cross-validation
cv.ridge = cv.glmnet(x,y,alpha = 0, family ="binomial")
#fit the final model on the training data
model = glmnet(x,y,alpha = 0, family = 'binomial', 
               lamda = cv.ridge$lambda.min)

fat_ridge=coef(cv.ridge, cv.ridge$lambda.min) 

fat_model_olink= fat_ridge%>%
  tidy()%>%
  write_clip()


#fatigue forest
beta_1 = fat_ridge
fat_or = as.data.frame(as.matrix(beta_1))
v = exp(fat_or$s1*1.44)
z = exp(fat_or$s1*0.56)
exp_beta_1_cp = exp(fat_or$s1)

fat_or = fat_or%>%
  mutate(OR=exp_beta_1_cp)%>%
  mutate(upper = v)%>%
  mutate(lower = z)


fat_or$predictor = row.names(fat_or)
fat_or = fat_or %>%
  filter(OR!=1)%>%
  filter(predictor!='(Intercept)')

fat_or$predictor[fat_or$predictor=='crf1a_sex' ]='Female'

#identify top quartiles to threshold number of mediators shown on plots
fat_or2 =fat_or %>%
  filter(OR> 1)%>%
  arrange(desc(s1))
quantile(fat_or2$OR, probs = seq(0.05,1, by = 0.05))
fat_or3 =fat_or %>%
  filter(OR< 1)%>%
  arrange(desc(s1))
quantile(fat_or3$OR, probs = seq(0.1,1, by = 0.1))

fat_or4= fat_or%>%
  filter(OR> 1.085 |OR<0.88)

fat_forest = ggplot(data =fat_or4,aes(y =reorder(predictor,OR, decreasing=F), xmax = upper, xmin = lower, x= OR))+
  geom_pointrange( show.legend = FALSE, shape = 20, fill = 'hotpink4', size = 1.1, color = 'hotpink4')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_classic()+labs(title = 'Fatigue (n=314) Vs Recovered (n=233)', x ='odds ratio', y = 'predictor')+
  theme(panel.background = element_rect (fill = NA, colour = "black"),panel.border = element_rect(colour='black', fill = NA),
        panel.grid.major.x=element_blank(),panel.grid.major.y=element_line())+
  theme(axis.text= element_text(size=13), axis.title =element_text(size=14))+
  theme(plot.title = element_text(size = 16));fat_forest




##### Make final box plot of nested cv accuracy----
cvm = cvm%>%
  rename(affective = 'nestedcv_aff$cvm' )%>%
  rename(cp = 'nestedcv_cp$cvm' )%>%
  rename(cognitive = 'nestedcv_b$cvm')%>%
  rename(fatigue ='nestedcv_fat$cvm')%>%
  rename(gi ='nestedcv_g$cvm' )

#make boxplot
head(cvm)
# make data long
cvm_long = cvm%>%
  pivot_longer(cols =affective:cognitive,
               names_to= 'Symptom',
               values_to = 'class')


p= ggplot(data=cvm_long, mapping=aes(x=Symptom, y=class ))
cvm_box=p+geom_boxplot(outlier.shape=NA,aes(color = Symptom))+ 
  theme(panel.background = element_rect (fill = 'white', colour = "black"))+
  scale_colour_manual(values = c('orange', 'turquoise4','tomato3','hotpink4', 'royalblue' ))+
  xlab('Symptom group') + coord_cartesian(ylim =c(0,1))+
  ylab('Classification error') + 
  theme(legend.position = 'none')+
  scale_x_discrete(limits=c('cp', 'fatigue','affective','cognitive', 'gi'), 
                   labels= c( 'Cardio-resp', 'Fatigue','Anxiety/Depression','Cognitive',  'GI'))+
  ggtitle("50 repeats 10-fold nested CV")+
  theme(axis.text= element_text(size=14), axis.title =element_text(size=16),
        axis.text.x= element_text(angle =30, hjust = 0.9))+
  theme(plot.title = element_text(size = 16, vjust = 2), panel.border = element_rect(colour='black', fill = NA));cvm_box

#auc curve
cvm_auc = cbind(aff_cvm_auc, cp_cvm_auc,fat_cvm_auc, b_auc,nestedcv_g_df_auc)
cvm_auc = cvm_auc%>%
  rename(affective = 'nestedcv_aff_auc$cvm' )
cvm_auc = cvm_auc%>%
  rename(cp = 'nestedcv_cp_auc$cvm' )%>%
  rename(cognitive = 'nestedcv_b_auc$cvm')%>%
  rename(fatigue ='nestedcv_fat_auc$cvm')%>%
  rename(gi ='nestedcv_g_auc$cvm' )
#make boxplot
head(cvm_auc)
# make data long
cvm_long_auc = cvm_auc%>%
  pivot_longer(cols =affective:gi,
               names_to= 'Symptom',
               values_to = 'AUC')
p= ggplot(data=cvm_long_auc, mapping=aes(x=Symptom, y=AUC ))
cvm_box_auc=p+geom_boxplot(outlier.shape=NA,aes(color = Symptom))+ 
  theme(panel.background = element_rect (fill = 'white', colour = "black"))+
  scale_colour_manual(values = c('orange', 'turquoise4','tomato3','hotpink4', 'royalblue' ))+
  xlab('Symptom group') + coord_cartesian(ylim =c(0,1))+
  ylab('AUC') + 
  theme(legend.position = 'none')+
  stat_compare_means( paired = FALSE, method = "kruskal.test",
                      label.y =0.5, step.increase = 0.26,label= 'p.signif',
                      tip.length = 0, vjust = 0.2)+
  scale_x_discrete(limits=c('cp', 'fatigue','affective','cognitive', 'gi'), 
                   labels= c( 'Cardio-resp', 'Fatigue','Anxiety/Depression','Cognitive',  'GI'))+
  ggtitle("50 repeats 10-fold nested CV")+
  theme(axis.text= element_text(size=14), axis.title =element_text(size=16),
        axis.text.x= element_text(angle =30, hjust = 0.9))+
  theme(plot.title = element_text(size = 16, vjust = 2), panel.border = element_rect(colour='black', fill = NA));cvm_box_auc
