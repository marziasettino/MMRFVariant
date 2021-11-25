## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning=FALSE, 
  message=FALSE,
  comment = "#>"
)

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
#  patient <- read.csv("~/MMRF_CoMMpass_IA15_PER_PATIENT")
#  trt <- read.csv("~/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP")
#  variant.ann <- read.csv("~/MMRF_CoMMpass_IA15a_All_Canonical_Variants")
#  
#  

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(variant.ann.example,options = list(scrollX = TRUE, keys = TRUE))


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(patient.example,options = list(scrollX = TRUE, keys = TRUE))


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(trt.example,options = list(scrollX = TRUE, keys = TRUE))


