## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning=FALSE, 
  message=FALSE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MMRFVariant)
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)
library(survminer)
library(survival)
library(formattable)

## ---- echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE---------------------
devtools::load_all(".")

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)

## ----figurename="workflow", echo=FALSE, fig.cap="workflow that describes graphically step by step the procedure carried out to perform this case of study ", out.width = '99%'----
knitr::include_graphics("imgs/workflow.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  ListGene<-c("KRAS", "NRAS","TP53","FAM46C","DIS3","BRAF")
#  
#  

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  variants.plot<-MMRFVariant_PlotVariantsbyGene(variant.ann,ListGene,height=20, width=30,topN=50,
#                                                filenm="PlotVariantsbyGene_heatmap")
#  variant.impact_effect
#  

## ----figurename="PlotVariantsbyGene_heatmap", echo=FALSE, fig.cap="Heatmap of the N# of variants occurrence in ListGene", out.width = '90%'----
knitr::include_graphics("imgs/PlotVariantsbyGene_heatmap.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  ListSNPs.bycount<-MMRFVariant_GetVariantsbyGene(variant.ann,ListGene)
#  ListSNPs<-ListSNPs.bycount$dbSNP
#  

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  impact.table<-MMRFVariant_getImpact(variant.ann,ListSNPs)
#  
#     #For semplification purposes, we visualize a subset of columns and rows
#     impact.table.sub<-dplyr::select(impact.table,dbSNP,Gene,REF,ALT,feature,Effect,
#                       SIFT_Impact,Polyphen_Impact,Impact)
#  
#      head(unique(impact.table.sub),10)

## ----figurename="ImpactTable", echo=FALSE, fig.cap="Impact table of each SNP in GeneList", out.width = '90%'----
knitr::include_graphics("imgs/ImpactTable.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  plot.impact.effect<-MMRFVariant_PlotbyEffect_Impact(variant.ann,ListSNPs,topN=50,height=30,
#                                                width=15, filenm="PlotbyEffectImpact")
#  
#  plot.impact.effect

## ----figurename="PlotEffectbyImpact", echo=FALSE, fig.cap="Impact table of each SNP in GeneList", out.width = '90%'----
knitr::include_graphics("imgs/PlotEffectbyImpact.png")

