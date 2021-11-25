## ----setup--------------------------------------------------------------------
library(MMRFVariant)

## ---- echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE---------------------
devtools::load_all(".")

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  variant.impact_effect<-MMRFVariant_PlotbyEffect_Impact(variant.ann,topN=50,
#                                                         height=10, width=15,
#                                                         filenm="PlotbyEffectImpact")
#  variant.impact_effect
#  

## ----figurename, echo=FALSE, fig.cap="Plot correlating Impact- Effect to variants", out.width = '90%'----
knitr::include_graphics("imgs/Plot_impact_Effect_Variants.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  ListGene<-c("PRDM16", "AGO1","TNNI3K")
#  
#  variantsbyGene.plot<-MMRFVariant_PlotVariantsbyGene(variant.ann,ListGene)
#  
#  variantsbyGene.plot
#  

## ----figurename2, echo=FALSE, fig.cap="Number of variants by gene in ListGene", out.width = '90%'----
knitr::include_graphics("imgs/PlotVariantsbyGene_heatmap.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  
#  
#  variantsbyGene.plot.all<-MMRFVariant_PlotVariantsbyGene(variant.ann)
#  
#  variantsbyGene.plot.all
#  

## ----fig, echo=FALSE, fig.cap="Number of variants by gene", out.width = '90%'----
knitr::include_graphics("imgs/PlotVariantsbyGene_heatmap_all.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  ListSNPs2<-c("rs755588843","rs200556051","rs745587729","rs2066497","rs760494041")
#  
#  impact.table<-MMRFVariant_getImpact(variant.ann,ListSNPs2)
#  
#  impact.table
#  

## ----fig4, echo=FALSE, fig.cap="ccc", out.width = '99%'-----------------------
knitr::include_graphics("imgs/Tab_getImpacts.png")

