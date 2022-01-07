library(MMRFVariant)
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)
library(survminer)
library(survival)
library(formattable)




patient<-MMRF_CoMMpass_IA15_PER_PATIENT
trt<-MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP
variant.ann<-MMRF_CoMMpass_IA15a_All_Canonical_Variants
#(GRCh37)

#hg19



#step 1: heatmap of the N# of variants occurrence in genes of <ListGene>

ListGene<-c("KRAS", "NRAS","TP53","FAM46C","DIS3","BRAF")


variants.plot<-MMRFVariant_PlotVariantsbyGene(variant.ann,ListGene,height=20, width=30,topN=50,
                                              filenm="PlotVariantsbyGene_heatmap")
#step 2 - Get the SNPs found in genes of <ListGene>

ListSNPs.bycount<-MMRFVariant_GetVariantsbyGene(variant.ann,ListGene)
ListSNPs<-ListSNPs.bycount$dbSNP


#step 3.2 -  generates the impact table (ordered by ascending SIFT and descending Poliphen)
impact.table<-MMRFVariant_getImpact(variant.ann,ListSNPs)

   #For semplification purposes, we visualize a subset of columns and rows
   impact.table.sub<-dplyr::select(impact.table,dbSNP,Gene,REF,ALT,feature,Effect,
                     SIFT_Impact,Polyphen_Impact,Impact)

    head(unique(impact.table.sub),10)
    
    
#step 4 - Plot the Impact-Effect of SNPs obtained from step 3

plot.impact.effect<-MMRFVariant_PlotbyEffect_Impact(variant.ann,ListSNPs,topN=50,height=30, 
                                              width=15, filenm="PlotbyEffectImpact")



#-------------------Survival plots---------------

KRAS_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"KRAS")$dbSNP
#NRAS_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"NRAS")$dbSNP
#TP53_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"TP53")$dbSNP
#FAM46C_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"FAM46C")$dbSNP
#DIS3_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"DIS3")$dbSNP
#BRAF_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"BRAF")$dbSNP

#for a better visualization, we take into account only four SNPs 
KRAS_SNPs<-head(KRAS_SNPs,4)

KRAS_surv.treatment<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     KRAS_SNPs,
                                     FilterBy="Treatment", 
                                     filename="KM_Plot_KRAS_treatment",
                                     xlim = c(100,3000),
                                     height=22,
                                     width=12,
                                     conf.range = FALSE,
                                     color = c("Dark2"))




KRAS_surv.Effect<-MMRFVariant_SurvivalKM(patient,  
                                          trt,
                                          variant.ann,
                                          KRAS_SNPs,
                                          FilterBy="Effect", 
                                          filename="KM_Plot_KRAS_effect",
                                          xlim = c(100,3000),
                                          height=22,
                                          width=12,
                                          conf.range = FALSE,
                                          color = c("Dark2"))


KRAS_surv.Stage<-MMRFVariant_SurvivalKM(patient,  
                                          trt,
                                          variant.ann,
                                          KRAS_SNPs,
                                          FilterBy="Stage", 
                                          filename="KM_Plot_KRAS_stage",
                                          xlim = c(100,3000),
                                          height=22,
                                          width=12,
                                          conf.range = FALSE,
                                          color = c("Dark2"))


KRAS_surv.Bestresp<-MMRFVariant_SurvivalKM(patient,  
                                          trt,
                                          variant.ann,
                                          KRAS_SNPs,
                                          FilterBy="Bestresp", 
                                          filename="KM_Plot_KRAS_bestresp",
                                          xlim = c(100,3000),
                                          height=22,
                                          width=12,
                                          conf.range = FALSE,
                                          color = c("Dark2"))

KRAS_surv.Gender<-MMRFVariant_SurvivalKM(patient,  
                                          trt,
                                          variant.ann,
                                          KRAS_SNPs,
                                          FilterBy="Gender", 
                                          filename="KM_Plot_KRAS_gender",
                                          xlim = c(100,3000),
                                          height=22,
                                          width=12,
                                          conf.range = FALSE,
                                          color = c("Dark2"))



