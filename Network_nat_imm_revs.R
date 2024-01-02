#18/10/2023 Claudia Efstathiou Nat Imm revisions 
#Olink Network analysis 

#libraries
library(tidyverse)
library(qgraph)
library(bootnet)
library(psych)
library(Hmisc)
library(corrplot)
library(readxl)
library(ggcorrplot)
library(ragg)
library(ggpubr)

#read data 
olink <- readRDS("olink_final.RDS")
olink_time <- readRDS("olink_time3m_natimm.rds")

olink_t <- left_join(olink_time, olink, by= "study_id")%>% distinct(study_id, .keep_all = T)

#selecting only those collected 3 months (12weeks) after discharge 

olink_3m <- olink_t[olink_t$Time_to_sample >= 12,] %>% distinct(study_id, .keep_all = T)
                                                            

names_cp <- c("IL1R2", "MATN2", "COLEC12", "IL1RL2", "CLEC4G", "CD4",
           "DPP10", "PCDH1", "ISM1", "CD276", "ROBO1", "ANGPTL2","ANGPTL4",
            "IL3RA", "C1QA","CTSO", 
           "TGFA", "TPSAB1", "CSF3",
           "PTH1R", "AMN", "CD160", "IL18R1")

names_Fa <- c("ANGPTL2","CD276", "CD4","CTSO",
               "CLEC4G", "COLEC12", "IL1RL2","ROBO1",
              "CSF3", "DPP10", "ENAH", "ERBB3", "EIF5A",
             "IL1R2", "IL3RA", "ISM1",  "MATN2", "TNFRSF13B",
              "PCDH1", "SCG3", "SIT1", "TFF2", "TGFA", "TPSAB1",
              "CCL11", "AMN","CCL26", "CD160")

names_GI <- c("AMBN", "C1QA", "CD4", "IL6", "EIF5A",
              "CSF3", "DPP10", "ERBB3", "FLT3LG", "FST", 
               "IL1R2", "ISM1", "ITM2A","IL1RL2","IL3RA", "CCL22",
              "MATN2", "SCG3", "GALNT3", "CLEC4G")

names_Anx <- c("AMBN", "ANGPTL2","C1QA",
               "COLEC12","CLEC4C","CD276", "LGALS4",
               "ERBB3", "IL1RL2","TGFA","EIF5A","TFF2", 
              "IL1R2", "IL3RA", "MATN2")

names_cog <- c("C1QA", "CCL11", "CD4",
              "COLEC12","CTSO", "DPP10","LIFR","CXADR",
              "LGALS4", "MATN2", "NFASC", "SPON1", "TNFRSF11B")

Rec <- subset(olink, GI == "Rec")
CP <- subset(olink, CP == "CP")
Fa <- subset(olink, Fatigue == "Fatigue")
GI <- subset(olink, GI == "GI")
Anx <- subset(olink, affective == "affective")
Cog <-subset(olink, brain_fog == "brain_fog")

#remove non-numeric variable 
ONum <- olink[,c(names_cp)] 


x <- cor(ONum)
y<- rcorr(as.matrix(ONum))

corrplot(cor(ONum))
ggcorrplot(cor(ONum))


#Network graphs with model selection CP vs rec ------

CP_onum <- CP[,c(names_cp)] 


est_CP <-  estimateNetwork(
  CP_onum,
  default = "cor",
  #"IsingFit" for binary, also "pcor" or "cor"
  corMethod = "spearman",
  threshold = "fdr")
CP_mat <- est_CP[["graph"]]


w<-qgraph(CP_mat, layout = "spring", sampleSize = nrow(CP_onum), 
       title = "Cardiorespiratory",
       details = T, colFactor = 1.5, fade = T,  curveAll = T,
       curveDefault = 0.5, normalize = T, repulsion = 1, directed = F, border.width = 0.5,
       edge.width = 0.8, minimum=0.15, theme= "colorblind", labels = names_cp,label.scale.equal=F)



#Network graphs with model selection Fatigue vs rec ------
 
Fa_onum <- Fa[,c(names_Fa)] 


est_Fa <-  estimateNetwork(
  Fa_onum,
  default = "cor",
  #"IsingFit" for binary, also "pcor" or "cor"
  corMethod = "spearman",
  threshold = "fdr")

plot(est_Fa, layout = "spring", fade = T, repulsion = 2, edge.width = 0.65, negDashed = T,
     label.prop = 0.5)

mat_fa <- est_Fa[["graph"]]

r<- qgraph(mat_fa, layout = "spring", sampleSize = nrow(Fa_onum), 
       title = "Fatigue",
       details = T, colFactor = 1.5, fade = T,  curveAll = T,
       curveDefault = 0.5, normalize = T, repulsion = 1, directed = F, border.width = 0.5,
       edge.width = 0.8, minimum=0.15, theme= "colorblind", labels = names_Fa,label.scale.equal=F)

#Network graphs with model selection GI vs rec ------

GI_onum <- GI[,c(names_GI)] 


est_GI <-  estimateNetwork(
  GI_onum,
  default = "cor",
  #"IsingFit" for binary, also "pcor" or "cor"
  corMethod = "spearman",
  threshold = "fdr")


mat_GI <- est_GI[["graph"]]

y<- qgraph(mat_GI, layout = "spring", sampleSize = nrow(GI_onum), 
       title = "GI",
       details = T, colFactor = 1.5, fade = T,  curveAll = T,
       curveDefault = 0.5, normalize = T, repulsion = 1, directed = F, border.width = 0.5,
       edge.width = 0.8, minimum=0.15, theme= "colorblind", labels = names_GI,label.scale.equal=F)

#Network graphs with model selection Anx vs rec ------

Anx_onum <- Anx[,c(names_Anx)] 

est_Anx <-  estimateNetwork(
  Anx_onum,
  default = "cor",
  #"IsingFit" for binary, also "pcor" or "cor"
  corMethod = "spearman",
  threshold = "fdr")


mat_Anx <- est_Anx[["graph"]]


i<- qgraph(mat_Anx, layout = "spring", sampleSize = nrow(Anx_onum), 
       title = "Anxiety/Depression",
       details = T, colFactor = 1.5, fade = T,  curveAll = T,
       curveDefault = 0.5, normalize = T, repulsion = 1, directed = F, border.width = 0.5,
       edge.width = 0.8, minimum=0.15, theme= "colorblind", labels = names_Anx,
       label.cex = 1,label.scale.equal=F)

#Network graphs with model selection cog vs rec ------

Cog_onum <- Cog[,c(names_cog)] 


est_Cog <-  estimateNetwork(
  Cog_onum,
  default = "cor",
  #"IsingFit" for binary, also "pcor" or "cor"
  corMethod = "spearman",
  threshold = "fdr")

mat_Cog <- est_Cog[["graph"]]


p<- qgraph(mat_Cog, layout = "spring", sampleSize = nrow(Cog_onum), 
       title = "Cognative",
       details = T, colFactor = 1.5, fade = T,  curveAll = T,
       curveDefault = 0.5, normalize = T, repulsion = 1, directed = F, border.width = 0.5,
       edge.width = 0.8, minimum=0.15, theme= "colorblind", labels = names_cog)

pdf("sm_Fa_network_pos_revx.pdf")

# Creating a plot
plot(r)

# Closing the graphical device
dev.off() 
#centality models ------
#CP = w, Fa = r, GI= y, Anx = i, cog = p


centrality(w, alpha = 1, posfun = abs, all.shortest.paths = FALSE,
           weighted = TRUE, signed = TRUE, R2 = FALSE)

CP_cen_tab <-centralityTable(w, standardized = TRUE,  relative = FALSE, weighted =
                  TRUE, signed = TRUE)
Fa_cen_tab<-centralityTable(r, standardized = TRUE,  relative = FALSE, weighted =
                              TRUE, signed = TRUE)
GI_cen_tab<-centralityTable(y, standardized = TRUE,  relative = FALSE, weighted =
                              TRUE, signed = TRUE)
Anx_cen_tab<-centralityTable(i, standardized = TRUE,  relative = FALSE, weighted =
                              TRUE, signed = TRUE)
Cog_cen_tab<-centralityTable(p, standardized = TRUE,  relative = FALSE, weighted =
                               TRUE, signed = TRUE)

CP_cen_tab1 <- CP_cen_tab %>%subset(measure=="Strength")
Symptom <- rep("Cardiorespiratory", 23)
CP_cen_tab1 <-cbind(CP_cen_tab1, Symptom)

Fa_cen_tab1 <- Fa_cen_tab %>%subset(measure=="Strength")
Symptom <- rep("Fatigue", 32)
Fa_cen_tab1 <-cbind(Fa_cen_tab1, Symptom)

GI_cen_tab1 <- GI_cen_tab %>%subset(measure=="Strength")
Symptom <- rep("GI", 18)
GI_cen_tab1 <-cbind(GI_cen_tab1, Symptom)

Anx_cen_tab1 <- Anx_cen_tab %>%subset(measure=="Strength")
Symptom <- rep("Anxiety", 9)
Anx_cen_tab1 <-cbind(Anx_cen_tab1, Symptom)

Cog_cen_tab1 <- Cog_cen_tab %>%subset(measure=="Strength")
Symptom <- rep("Cognative", 12)
Cog_cen_tab1 <-cbind(Cog_cen_tab1, Symptom)


#clusteringPlot(w, scale = c("raw0", "raw", "z-scores", "relative"), signed = FALSE, theme_bw = TRUE, print = TRUE,
             # verbose = TRUE,orderBy = "default",
              # decreasing = FALSE)
b<-centrality_auto(w, weighted = TRUE, signed = TRUE)

n<-centralityTable(w, standardized = TRUE,  relative = FALSE, weighted =
                     TRUE, signed = TRUE)

#single graphs 
a<- CP_cen_tab1 %>%  mutate(node = fct_reorder(node, value)) %>% ggplot(aes(x=node, y=value))+
  geom_point(stat = "identity")+ coord_flip()+
  ggtitle("Cardiorespiratory Network Centrality") + xlab("Mediator")+
  geom_line(color= "tomato3", aes(group = measure),size = 1)+theme_bw();a

s<-Fa_cen_tab1 %>%  mutate(node = fct_reorder(node, value)) %>% ggplot(aes(x=node, y=value))+
  geom_point(stat = "identity")+ coord_flip()+
  ggtitle("Fatigue Network Centrality") + xlab("Mediator")+
  geom_line(color= "hotpink4", aes(group = measure),size = 1)+theme_bw();s

d<-GI_cen_tab1 %>%  mutate(node = fct_reorder(node, value)) %>% ggplot(aes(x=node, y=value))+
  geom_point(stat = "identity")+ coord_flip()+
  ggtitle("GI Network Centrality") + xlab("Mediator")+
  geom_line(color= "royalblue", aes(group = measure),size = 1)+theme_bw();d

f<-Anx_cen_tab1 %>%  mutate(node = fct_reorder(node, value)) %>% ggplot(aes(x=node, y=value))+
  geom_point(stat = "identity")+ coord_flip()+
  ggtitle("Anxiety/Depression Network Centrality") + xlab("Mediator")+
  geom_line(color= "orange", aes(group = measure), size = 1)+theme_bw();f

g<-Cog_cen_tab1 %>%  mutate(node = fct_reorder(node, value)) %>% ggplot(aes(x=node, y=value))+
  geom_point(stat = "identity")+ coord_flip()+
  ggtitle("Cognitive Network Centrality") + xlab("Mediator")+
  geom_line(color= "turquoise4", aes(group = measure), size = 1)+theme_bw();g

ggarrange(a,s,f,d,g, nrow=3, ncol=2)
