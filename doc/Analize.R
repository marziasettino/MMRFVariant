## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  patient <- read.csv("~/MMRF_CoMMpass_IA15_PER_PATIENT")
#  trt <- read.csv("~/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP")
#  variant.ann <- read.csv("~/MMRF_CoMMpass_IA15a_All_Canonical_Variants")
#  
#  

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Treatment",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename2, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Treatment", out.width = '99%'----
knitr::include_graphics("imgs/KM_Surv_treatment.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Stage",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename3, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Stage", out.width = '50%'----
knitr::include_graphics("imgs/KM_Surv_stage.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Ethnicity",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename4, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Ethnicity", out.width = '99%'----
knitr::include_graphics("imgs/KM_Surv_ethnicity.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Bestresp",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename5, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Bestresp", out.width = '99%'----
knitr::include_graphics("imgs/KM_Surv_Bestresp.png")

