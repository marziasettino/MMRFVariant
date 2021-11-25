library(MMRFVariant)
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)
library(survminer)

library(formattable)




patient<-MMRF_CoMMpass_IA15_PER_PATIENT
trt<-MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP
variant.ann<-MMRF_CoMMpass_IA15a_All_Canonical_Variants
#(GRCh37)

#hg19



Listvariant<-c("rs2157615","rs61731685","rs370560636","rs116293337")
Listvariant1<-c("rs2157615")
Listvariant2<-c("rs61731685")
Listvariant3<-c("rs370560636")
Listvariant4<-c("rs116293337")

ListGene<-c("PRDM16", "AGO1","TNNI3K") 

ListSNPs<-c("rs10000250","rs10006630","rs10015833","rs10082511")

ListSNPs2<-c("rs755588843","rs200556051","rs745587729","rs2066497","rs760494041")


#surv<-MMRFVariant_SurvivalKM(patient,  
                         #     trt,
                          #    variant.ann,
                          #    Listvariant,
                           #   FilterBy="treatment", 
                           #   filename="Survival_Compact",
                           #   height=10,
                          #    width=20,
                           #   xlim = c(100,3000),
                          #    conf.range = FALSE,
                           #   color = c("Dark2"))






#-----------------------------------


summary.plot<-MMRFVariant_PlotbyEffect_Impact(variant.ann,topN=50,height=10, width=15, filenm="PlotbyEffectImpact")


variants.plot<-MMRFVariant_PlotVariantsbyGene(variant.ann,ListGene,filenm="PlotVariantsbyGene_heatmap")
variants.plot.all<-MMRFVariant_PlotVariantsbyGene(variant.ann,topN=20,filenm="PlotVariantsbyGene_heatmap_all")

#Compact plot
        surv<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     Listvariant,
                                     FilterBy="treatment",
                                     expand=FALSE,
                                     filename="KM_Plot_Compact",
                                     xlim = c(100,3000),
                                     conf.range = FALSE,
                                     color = c("Dark2"))

#Expanded plot
surv.summary<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     Listvariant,
                                     expand=TRUE,
                                     FilterBy="treatment", 
                                     filename="KM_Plot_Expanded",
                                     xlim = c(100,3000),
                                     conf.range = FALSE,
                                     color = c("Dark2"))





#---------------
surv<-MMRFVariant_SurvivalKM(patient,  
                             trt,
                             variant.ann,
                             Listvariant,
                             FilterBy="stage",
                             expand=FALSE,
                             filename="KM_Plot_Compact",
                             xlim = c(100,3000),
                             conf.range = FALSE,
                             color = c("Dark2"))

#Expanded plot
surv.summary<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     Listvariant,
                                     expand=TRUE,
                                     FilterBy="stage", 
                                     filename="KM_Plot_Expanded",
                                     xlim = c(100,3000),
                                     conf.range = FALSE,
                                     color = c("Dark2"))


#---------------
surv<-MMRFVariant_SurvivalKM(patient,  
                             trt,
                             variant.ann,
                             Listvariant,
                             FilterBy="bestresp",
                             expand=FALSE,
                             filename="KM_Plot_Compact",
                             xlim = c(100,3000),
                             conf.range = FALSE,
                             color = c("Dark2"))

#Expanded plot
surv.summary<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     Listvariant,
                                     expand=TRUE,
                                     FilterBy="bestresp", 
                                     filename="KM_Plot_Expanded",
                                     xlim = c(100,3000),
                                     conf.range = FALSE,
                                     color = c("Dark2"))

#---------------
surv<-MMRFVariant_SurvivalKM(patient,  
                             trt,
                             variant.ann,
                             Listvariant,
                             FilterBy="race",
                             expand=FALSE,
                             filename="KM_Plot_Compact",
                             xlim = c(100,3000),
                             conf.range = FALSE,
                             color = c("Dark2"))

#Expanded plot
surv.summary<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     Listvariant,
                                     expand=TRUE,
                                     FilterBy="race", 
                                     filename="KM_Plot_Expanded",
                                     xlim = c(100,3000),
                                     conf.range = FALSE,
                                     color = c("Dark2"))

#--------------------

impact.table<-MMRFVariant_getImpact(variant.ann,ListSNPs2)








