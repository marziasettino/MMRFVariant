
#private
MMRFVariant_getVariantAnn<- function(variant.ann){
  names(variant.ann)[1]<-"public_id"
  variant.ann$public_id<-substr(variant.ann$public_id,1,9)
  
  
  names(variant.ann)[names(variant.ann) == 'ANN....EFFECT'] <- "Effect"
  names(variant.ann)[names(variant.ann) == 'ANN....GENE'] <- "Gene"
  names(variant.ann)[names(variant.ann) == 'ANN....BIOTYPE'] <- "Biotype"
  names(variant.ann)[names(variant.ann) == 'ANN....IMPACT'] <- "Impact"
  names(variant.ann)[names(variant.ann) == 'ID'] <- "dbSNP"
  names(variant.ann)[names(variant.ann) == 'ANN....FEATURE'] <- "feature_type"
  names(variant.ann)[names(variant.ann) == 'ANN....FEATUREID'] <- "feature"
  names(variant.ann)[names(variant.ann) == 'dbNSFP_SIFT_score'] <- "SIFT"
  names(variant.ann)[names(variant.ann) == 'dbNSFP_Polyphen2_HDIV_score'] <- "Polyphen2"
  
  
  
  
  return(variant.ann)
}



#private


sift_impact <- function(x) {
  if (x>=0 && x<=0.05) {
    "Intolerant"
  } else if (x>=0.051 && x<=0.1) {
    "Potentially Intolerant" 
  } else if (x>=0.101 && x<=0.2) {
    "Borderline"
  } else if (x>=0.201 && x<=1) {
    "Tolerant"
  } else {
    stop("Invalid `x` value")
  }
}

polyphen_impact <- function(x) {
  if (x>=0 && x<=0.999) {
    "Benign"
  } else if (x>=1 && x<=1.24) {
    "Borderline" 
  } else if (x>=1.25 && x<=1.49) {
    "Potentially damaging"
  } else if (x>=1.5 && x<=1.99) {
    "Possibly damaging"
  } else if (x>=2) {
    "probably damaging"  
  } else {
    stop("Invalid `x` value")
  }
}




#private
table_formatter <- function(df) {
  
  df<-dplyr::select(df,dbSNP,feature,min_SIFT,max_polyphen,
                    SIFT_Impact,Polyphen_Impact,Effect,Gene,
                    REF,ALT,Biotype,Impact,feature_type)
  
  
  
  
  
  df<-formattable(df, list(
    SIFT_Impact = formatter("span", 
                            style = ~style(display = "block",
                                           font.weight = "bold", 
                                           color = "white",
                                           "border-radius" = "4px",
                                           "padding-right" = "4px",
                                           "background-color" =  
                                             ifelse(SIFT_Impact == "Intolerant","red",
                                                    ifelse(SIFT_Impact == "Potentially Intolerant","orange",
                                                           ifelse(SIFT_Impact == "Borderline","blue",    
                                                                  ifelse(SIFT_Impact=="Tolerant", "lightblue",NA)))))),
    
    Polyphen_Impact = formatter("span", 
                                style = ~style(display = "block",
                                               font.weight = "bold", 
                                               color = "white",
                                               "border-radius" = "4px",
                                               "padding-right" = "4px",
                                               "background-color" =  
                                                 ifelse(Polyphen_Impact == "probably damaging","red",
                                                        ifelse(Polyphen_Impact == "Possibly damaging","orange",
                                                               ifelse(Polyphen_Impact == "Potentially damaging","green",
                                                                      ifelse(Polyphen_Impact == "Borderline","blue",
                                                                             ifelse(Polyphen_Impact=="Benign", "lightblue",NA))))))),
    
    Gene = formatter("span", style = ~style(font.weight = "bold", color =  "blue")),
    dbSNP = formatter("span", style = ~style(font.weight = "bold", color =  "blue"))
    
    
  ))
  
  return(df) 
}



#private


addImpacts <- function(df) {
  
  #df<- MMRFVariant_getVariantAnn(df)
  df<-dplyr::select(df,public_id,dbSNP,
                    Effect,Gene,REF,ALT,
                    Biotype,Impact,feature_type,
                    feature,SIFT,Polyphen2)
  
  #add column (Max damaging impact)
  df$min_SIFT<-""
  df$SIFT_Impact<-""
  
  df$max_polyphen<-""
  df$Polyphen_Impact<-""
  
  
  
  for(i in 1:nrow(df)){ 
    sift<-df$SIFT[i]
    polyphen<-df$Polyphen2[i]
    
    print(paste0("i: ", i))
    print(paste0("SIFT: ", sift))
    print(paste0("POLYPHEN: ", polyphen))
    
    
    result = tryCatch({
      sift<-unlist(strsplit(sift, ",")) %>% as.numeric(sift)
      min<-min(sift)
      df$min_SIFT[i]<-min
      df$SIFT_Impact[i]<-sift_impact(sift<-df$min_SIFT[i])
      
    }, warning = function(w) {
      df$min_SIFT[i]<-NA
    }, error = function(e) {
      df$min_SIFT[i]<-NA
    }
    ) 
    
    
    result = tryCatch({
      polyphen<-unlist(strsplit(polyphen, ",")) %>% as.numeric(polyphen)
      max<-max(polyphen)
      df$max_polyphen[i]<-max
      
      df$Polyphen_Impact[i]<-polyphen_impact(sift<-df$max_polyphen[i])
      
    }, warning = function(w) {
      df$max_polyphen[i]<-NA
    }, error = function(e) {
      df$max_polyphen[i]<-NA
    }
    )     
    
    
    
  }  #for  
  
  return(df) 
}

#------------------------------------------------



#' @title MMRFVariant_SurvivalKM_Summary
#' @description Creates a survival plot from MMRF-RG patient clinical data
#' using survival library. It uses the fields D_PT_deathdy, D_PT_PRIMARYREASON and D_PT_lstalive
#' columns for groups.
#' @param patient is the data.frame of the patient clinical data downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15_PER_PATIENT file) and imported into environment.
#' @param trt is the data.frame of the patient clinical data (i.e. treatment-response) downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP file) and imported into environment.
#' @param variant.ann is the data.frame of the annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15a_All_Canonical_Variants file) and imported into environment.
#' @param Listvariant is the list of the variants to analyze.
#' @param risk.table show or not the risk table
#' @param legend Legend title of the figure
#' @param xlim x axis limits e.g. xlim = c(0, 1000). Present narrower X axis, but not affect survival estimates.
#' @param main main title of the plot
#' @param labels labels of the plot
#' @param ylab y axis text of the plot
#' @param xlab x axis text of the plot
#' @param filename The name of the pdf file.
#' @param color Define the colors/Pallete for lines.
#' @param width Image width
#' @param height Image height
#' @param pvalue show p-value of log-rank test
#' @param conf.range  show confidence intervals for point estimates of survival curves.
#' @param dpi Figure quality
#' @import survminer
#' @import survival
#' @import gridExtra
#' @import ggplot2
#' @import stringr
#' @return Survival plot
#' @examples
#' 
#' patient <- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                    "MMRF_0002","MMRF_0003",
#'                                    "MMRF_0004","MMRF_0005",
#'                                    "MMRF_0006","MMRF_0007",
#'                                    "MMRF_0008","MMRF_0009"),
#'                       D_PT_PRIMARYREASON = c("Death",NA,NA,"Death","Death", 
#'                                              "Death","NA","NA","Death","Death"),
#'                       D_PT_deathdy =  c(NA,NA,2226,172,NA,NA,1328,681,786,NA),
#'                       D_PT_lstalive = c(250,356,NA,NA,1814,223,NA,NA,NA,1450),
#'                       DEMOG_GENDER = c(rep(1,5),rep(2,5)), #2 stands for female, 1 standas for male#'                       
#'                       DEMOG_ETHNICITY=c(rep("Hispanic or Latino",5),rep("Not Hispanic or Latino",5)),
#'                       D_PT_issstage_char=c(rep("Stage III",3),rep("Stage II",2),rep("Stage I",5))
#'  )
#' 
#' 
#' trt<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                               "MMRF_0002","MMRF_0003",
#'                               "MMRF_0004","MMRF_0005",
#'                               "MMRF_0006","MMRF_0007",
#'                               "MMRF_0008","MMRF_0009"),
#'                  trtclass=c(rep("Bortezomib-based",2),rep("IMIDs-based",5),rep("combined bortezomib/IMIDs-based",3)),                                                    
#'                  bestresp=c(rep("Partial Responsed",5),rep("Very Good Partial Response",5))               
#'                                    
#'  )
#' 
#' 
#' 
#'     surv<-MMRFVariant_SurvivalKM_Summary(patient,  
#                                          trt,
#                                          variant.ann,
#                                          Listvariant1,
#                                          FilterBy="treatment", 
#                                          filename=NULL,
#                                          xlim = c(100,3000),
#                                          conf.range = FALSE,
#                                          color = c("Dark2"))

MMRFVariant_SurvivalKM_Single <- function(
  patient,
  trt,
  variant.ann,
  list.variant,
  risk.table = FALSE,
  FilterBy = NULL,
  legend = "Legend",
  labels = NULL,
  xlim = NULL,
  main = "Kaplan-Meier Survival Curve",
  ylab = "Probability of survival",
  xlab = "Time since diagnosis (days)",
  filename = "survival.pdf",
  color = NULL,
  height = 5,
  width = 12,
  dpi = 300,
  pvalue = TRUE,
  conf.range = TRUE) {
  
  patient<-MMRFVariant_GetSamplesbyVariant(variant.ann,patient,list.variant)
  
  
  
  
  condition <- c("race","stage","treatment","bestresp","gender","dbSNP","effect","biotype")
  parameter <- c("DEMOG_ETHNICITY", "D_PT_issstage_char","trtclass","bestresp","DEMOG_GENDER","ID","ANN....EFFECT","ANN....BIOTYPE")
  tab.condition <- data.frame(condition, parameter)
  
  
  names(patient)[1]<-"public_id"
  
  patient<-dplyr::select(patient,public_id,D_PT_PRIMARYREASON,
                         D_PT_lstalive,D_PT_deathdy,
                         DEMOG_ETHNICITY,D_PT_issstage_char,DEMOG_GENDER
  )
  
  
  
  df.merge<-merge(x = patient, y = trt, by = "public_id", type=left)
  
  
  
  if (is.null(color)) {
    color <- rainbow(length(unique(patient[, FilterBy])))
  }
  
  group <- NULL
  
  
  
  
  if (is.null(FilterBy)) {
    stop("Please provide a value for FilterBy parameter")
  } else {
    par<-tab.condition[tab.condition$condition==FilterBy,]$parameter  #check tab.condition 
    FilterBy<-par
    
    if (length(unique(df.merge[, FilterBy])) == 1) {
      #  stop(
      paste0(
        "Only this group found:\n",
        unique(df.merge[, FilterBy])
      )
      #)
    }
  }#end
  
  
  
  notDead <- is.na(df.merge$D_PT_deathdy)
  
  if (any(notDead == TRUE)) {
    df.merge[notDead, "D_PT_deathdy"] <- df.merge[notDead, "D_PT_lstalive"]
  }
  
  #TRUE(DEAD)/FALSE (ALIVE)
  df.merge$s <- grepl("Death", df.merge$D_PT_PRIMARYREASON, ignore.case = TRUE)
  
  
  
  df.merge$type <- as.factor(df.merge[, FilterBy])
  df.merge <-  df.merge[, c("D_PT_deathdy", "s", "type")]
  #   formula 
  f.m <-formula(survival::Surv(as.numeric(df.merge$D_PT_deathdy), event = df.merge$s) ~ df.merge$type)
  fit <- do.call(survival::survfit, list(formula = f.m, data = df.merge))
  
  label.add.n <- function(x) {
    na.idx <- is.na(df.merge[, "D_PT_deathdy"])
    negative.idx <- df.merge[, "D_PT_deathdy"] < 0
    idx <- !(na.idx | negative.idx)
    return(paste0(x, " (n = ",
                  sum(df.merge[idx, "type"] == x), ")"))
  }
  
  if (is.null(labels)) {
    d <- survminer::surv_summary(fit, data = df.merge)
    order <-
      unname(sapply(levels(d$strata), function(x)
        unlist(str_split(x, "="))[2]))
    labels <- sapply(order, label.add.n)
  }
  
  
  
  
  if (length(xlim) == 1) {
    xlim <- c(0, xlim)
  }
  
  
  suppressWarnings({
    
    
    surv <- survminer::ggsurvplot( 
      fit,
      risk.table = risk.table,
      df.merge,
      pval = pvalue,
      conf.range = conf.range,
      xlim = xlim,
      main = main,
      xlab = xlab,
      legend.title = legend,
      palette =  color,
      legend = "right",
      legend.labs = levels(FilterBy)
      
    )
  })
  
  
  list.variant.str<-toString(list.variant)
  
  
  surv<-surv+ labs(title = paste("dbSNP Variant:",list.variant.str))
  
  
  
  if (!is.null(filename)) {
    
    filenm<-paste0(filename,".pdf")
    path<-file.path(getwd())
    path<-paste0(path,"/","ResultsPlot","/",filenm)
    
    ggsave(
      surv$plot,
      filename = path,
      # device = pdf,
      width = width,
      height = height,
      units = "in"
    )
    message(paste0("File saved as: ", path))
    
    
    if (risk.table) {
      g1 <- ggplotGrob(surv$plot)
      g2 <- ggplotGrob(surv$table)
      min_ncol <- min(ncol(g2), ncol(g1))
      g <-
        gridExtra::gtable_rbind(g1[, 1:min_ncol], g2[, 1:min_ncol], size = "last")
      ggsave(
        g,
        filename = filename,
        width = width,
        height = height,
        dpi = dpi
      )
    }
  } 
  
  return(surv)
  
}






