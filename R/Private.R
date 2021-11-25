
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
  
  df<-dplyr::select(df,public_id,dbSNP,Effect,Gene,
                    REF,ALT,Biotype,Impact,feature_type,feature,
                    min_SIFT,max_polyphen,SIFT_Impact,Polyphen_Impact)
  
  
  
  
  
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


