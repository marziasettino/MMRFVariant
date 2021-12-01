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



#Listvariant<-c("rs2157615","rs61731685","rs370560636","rs116293337")
#Listvariant1<-c("rs2157615")
#Listvariant2<-c("rs61731685")
#Listvariant3<-c("rs370560636")
#Listvariant4<-c("rs116293337")

ListGene2<-c("PRDM16", "AGO1","TNNI3K") 

#ListSNPs<-c("rs10000250","rs10006630","rs10015833","rs10082511")

#ListSNPs<-c("rs116293337","rs372554181","rs377287282","rs782606620","rs372344564", 
#            "rs370824753", "rs116293337", "rs782606620")

ListSNPs<-c("rs116293337", "rs782606620", "rs370165512", "rs372344564" ,"rs377287282",
"rs116293337", "rs782606620", "rs372554181" ,"rs781956502" ,"rs782752158",
"rs372124885" ,"rs78872461" , "rs768882275", "rs782130293" ,"rs368220735",
"rs117375787" ,"rs369885179", "rs117569183" ,"rs370165512" ,"rs78304725" ,
"rs371126656" ,"rs372344564" ,"rs372124885" ,"rs371353327" ,"rs371252738",
"rs182723052", "rs373312097", "rs371714049" ,"rs377287282", "rs376449862",
"rs369886334" ,"rs148494108", "rs373130593", "rs150237274" ,"rs78872461" ,
"rs116162069", "rs141950820", "rs369706433" ,"rs373120212", "rs372554181",
"rs144089698" ,"rs140030392" ,"rs370824753" ,"rs372724023", "rs61731685" ,
"rs142696826" ,"rs781956502" ,"rs782752158", "rs2157615",   "rs2157615")  




ListSNPs<-c("rs61731685","rs2157615","rs2157615","rs367953381","rs367953381",
 "rs370560636", "rs370560636", "rs73880715","rs116293337","rs116293337",
 "rs116293337", "rs116293337", "rs116293337", "rs116293337", "rs116293337",
 "rs116293337", "rs116293337", "rs116293337", "rs116293337", "rs116293337",
 "rs116293337", "rs116293337", "rs782606620", "rs782606620", "rs782606620",
 "rs782606620", "rs782606620", "rs782606620", "rs782606620", "rs782606620",
 "rs782606620", "rs782606620", "rs782606620", "rs782606620","rs782606620",
 "rs782606620", "rs11554290",  "rs11554290","rs370615998","rs370615998",
 "rs17851045",  "rs17851045" , "rs61733526",  "rs372344564","rs372344564",
 "rs372344564", "rs372344564", "rs372344564", "rs372344564","rs372344564")





#-----------------------------------



#step 1 - heatmap of the the top #N of SNPs ocurrence
variants.plot.all<-MMRFVariant_PlotVariantsbyGene(variant.ann,topN=100,
                                                  filenm="PlotVariantsbyGene_heatmap_all")
#step 2 - Get the top #N genes and the SNPs
ListSNPs<-unique(variants.plot.all$data$dbSNP)
ListGene<-unique(variants.plot.all$data$Gene)

#step 3 - Plot the occurrence #N of the SNPs identified in the genes list <ListGene> obtained from step 2
variants.plot<-MMRFVariant_PlotVariantsbyGene(variant.ann,ListGene,height=30, width=40,topN=50,
                                              filenm="PlotVariantsbyGene_heatmap")


#step 3 - Get the Impacts table of the SNPs in the ListSNPs obtained from step 2
variants.target<- variant.ann[variant.ann$ID %in% ListSNPs, ] 

 #step 3.1 -  excludes items with no impact values 
 variants.target<-subset(variants.target, variants.target$dbNSFP_SIFT_score != "." & 
                         dbNSFP_Polyphen2_HDIV_score != ".")

 #step 3.2 -  generates the impact table 
 impact.table<-MMRFVariant_getImpact(variants.target,ListSNPs)

#step 3 - Plot the Impact-Effect of SNPs obtained from step 3

summary.plot<-MMRFVariant_PlotbyEffect_Impact(variants.target,topN=50,height=30, 
                                              width=15, filenm="PlotbyEffectImpact")




#step 4
#Compact plot
        surv<-MMRFVariant_SurvivalKM(patient,  
                                     trt,
                                     variant.ann,
                                     ListSNPs,
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
                                     ListSNPs,
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
                             ListSNPs,
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
                                     ListSNPs,
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
                             ListSNPs,
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
                                     ListSNPs,
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
                             ListSNPs,
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
                                     ListSNPs,
                                     expand=TRUE,
                                     FilterBy="race", 
                                     filename="KM_Plot_Expanded",
                                     xlim = c(100,3000),
                                     conf.range = FALSE,
                                     color = c("Dark2"))

#--------------------










