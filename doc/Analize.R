## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  patient <- read.csv("~/MMRF_CoMMpass_IA15_PER_PATIENT")
#  trt <- read.csv("~/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP")
#  variant.ann <- read.csv("~/MMRF_CoMMpass_IA15a_All_Canonical_Variants")
#  
#  

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  Listvariant<-c("rs2157615","rs61731685","rs370560636","rs116293337")
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                Listvariant,
#                                expand = TRUE,
#                                FilterBy="treatment",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename2, echo=FALSE, fig.cap="KM Survival curves ", out.width = '99%'----
knitr::include_graphics("imgs/KM_Survival_Expanded.png")

