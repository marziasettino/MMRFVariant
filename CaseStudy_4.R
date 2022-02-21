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
                                              filenm="PlotVarbyGene_heatmap")
#step 2 - Get the SNPs found in genes of <ListGene>

ListSNPs.bycount<-MMRFVariant_GetVariantsbyGene(variant.ann,ListGene)
ListSNPs<-ListSNPs.bycount$dbSNP


#step 3 -  Perform the impact table (ordered by ascending SIFT and descending Poliphen)
impact.table<-MMRFVariant_getImpact(variant.ann,ListSNPs)

   #For semplification purposes, we visualize a subset of columns and rows
   impact.table.sub<-dplyr::select(impact.table,dbSNP,Gene,REF,ALT,feature,Effect,
                     SIFT_Impact,Polyphen_Impact,Impact)

    head(unique(impact.table.sub),10)
    
    
#step 4 - Plot the Impact-Effect of SNPs obtained from step 3

plot.impact.effect<-MMRFVariant_PlotbyEffect_Impact(variant.ann,ListSNPs,topN=50,height=30, 
                                              width=15, filenm="PlotbyEffectImpact")



#-------------------Survival plots---------------

#KRAS_SNPs.tab<-MMRFVariant_GetVariantsbyGene(variant.ann,"KRAS")
#KRAS_SNPs<-KRAS_SNPs.tab[order(KRAS_SNPs.tab$count, decreasing = TRUE),]$dbSNP

NRAS_SNPs.tab<-MMRFVariant_GetVariantsbyGene(variant.ann,"NRAS")
NRAS_SNPs<-NRAS_SNPs.tab[order(NRAS_SNPs.tab$count, decreasing = TRUE),]$dbSNP


#TP53_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"TP53")$dbSNP
#FAM46C_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"FAM46C")$dbSNP
#DIS3_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"DIS3")$dbSNP
#BRAF_SNPs<-MMRFVariant_GetVariantsbyGene(variant.ann,"BRAF")$dbSNP

#We take into account only the top six SNPs by occurrence (see <KRAS_SNPs.tab>)
NRAS_SNPs<-head(NRAS_SNPs,5)


NRAS_surv.treatment<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     NRAS_SNPs,
                                     FilterBy="Treatment", 
                                     filename="KM_Plot_NRAS_treatment",
                                     xlim = c(100,3000),
                                     height=22,
                                     width=12,
                                     conf.range = FALSE,
                                     color = c("Dark2"))


# (*) SNPs are discarded if:
# a) only a group with respect to FilterBy parameter is found
# b) pvalue is>=0.05
#One descarded the SNPs, the resulting SNP dataset is: 

NRAS_SNPs.disc<-c("rs11554290","rs121913254","rs121913237", "rs121913255")


    
NRAS_surv.treatment.disc<-MMRFVariant_SurvivalKM(patient,  
                                            trt,
                                            variant.ann,
                                            NRAS_SNPs.disc,
                                            FilterBy="Treatment", 
                                            filename="KM_Plot_NRAS_treatment_disc",
                                            xlim = c(100,2000),
                                            height=15,
                                            width=12,
                                            conf.range = FALSE,
                                            color = c("Dark2"))






NRAS_surv.Effect<-MMRFVariant_SurvivalKM(patient,  #no significant results are found (all pvalue>0.05)
                                         trt,
                                         variant.ann,
                                         NRAS_SNPs,
                                         FilterBy="Effect", 
                                         filename="KM_Plot_NRAS_effect",
                                         xlim = c(100,100),
                                         height=22,
                                         width=12,
                                         conf.range = FALSE,
                                         color = c("Dark2"))


# see (*)
NRAS_surv.Stage<-MMRFVariant_SurvivalKM(patient,  
                                         trt,
                                         variant.ann,
                                         NRAS_SNPs,
                                         FilterBy="Stage", 
                                         filename="KM_Plot_NRAS_stage",
                                         xlim = c(100,3000),
                                         height=22,
                                         width=12,
                                         conf.range = FALSE,
                                         color = c("Dark2"))

NRAS_SNPs.disc<-c("rs121913254")



NRAS_surv.Stage.disc<-MMRFVariant_SurvivalKM(patient,  
                                        trt,
                                        variant.ann,
                                        NRAS_SNPs.disc,
                                        FilterBy="Stage", 
                                        filename="KM_Plot_NRAS_stage_disc",
                                        xlim = c(100,2000),
                                        height=22,
                                        width=12,
                                        conf.range = FALSE,
                                        color = c("Dark2"))


# see (*)
NRAS_surv.Bestresp<-MMRFVariant_SurvivalKM(patient,  
                                        trt,
                                        variant.ann,
                                        NRAS_SNPs,
                                        FilterBy="Bestresp", 
                                        filename="KM_Plot_NRAS_bestresp",
                                        xlim = c(100,3000),
                                        height=22,
                                        width=12,
                                        conf.range = FALSE,
                                        color = c("Dark2"))



NRAS_SNPs.disc<-c("rs11554290","rs121913254","rs121434595", "rs121913237")

NRAS_surv.Bestresp.disc<-MMRFVariant_SurvivalKM(patient,  
                                           trt,
                                           variant.ann,
                                           NRAS_SNPs.disc,
                                           FilterBy="Bestresp", 
                                           filename="KM_Plot_NRAS_bestresp_disc",
                                           xlim = c(100,3000),
                                           height=22,
                                           width=12,
                                           conf.range = FALSE,
                                           color = c("Dark2"))






# see (*)

NRAS_surv.Gender<-MMRFVariant_SurvivalKM(patient,  #All SNPs have pvalue<=0.05
                                           trt,
                                           variant.ann,
                                           NRAS_SNPs,
                                           FilterBy="Gender", 
                                           filename="KM_Plot_NRAS_gender",
                                           xlim = c(100,3000),
                                           height=22,
                                           width=12,
                                           conf.range = FALSE,
                                           color = c("Dark2"))




# see (*)

NRAS_surv.Biotype<-MMRFVariant_SurvivalKM(patient,  #All SNPs have have only a group with respect to FilterBy parameter 
                                         trt,
                                         variant.ann,
                                         NRAS_SNPs,
                                         FilterBy="Biotype", 
                                         filename="KM_Plot_NRAS_biotype",
                                         xlim = c(100,3000),
                                         height=22,
                                         width=12,
                                         conf.range = FALSE,
                                         color = c("Dark2"))



# see (*)
NRAS_surv.Ethnicity<-MMRFVariant_SurvivalKM(patient,  
                                          trt,
                                          variant.ann,
                                          NRAS_SNPs,
                                          FilterBy="Ethnicity", 
                                          filename="KM_Plot_NRAS_ethnicity",
                                          xlim = c(100,3000),
                                          height=22,
                                          width=12,
                                          conf.range = FALSE,
                                          color = c("Dark2"))

NRAS_SNPs.disc<-c("rs11554290","rs121913254","rs121913237")

NRAS_surv.Ethnicity<-MMRFVariant_SurvivalKM(patient,  
                                            trt,
                                            variant.ann,
                                            NRAS_SNPs.disc,
                                            FilterBy="Ethnicity", 
                                            filename="KM_Plot_NRAS_ethnicity_disc",
                                            xlim = c(100,3000),
                                            height=22,
                                            width=12,
                                            conf.range = FALSE,
                                            color = c("Dark2"))







