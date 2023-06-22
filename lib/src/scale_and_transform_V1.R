rm(list=ls())
#SCALE_&_TRANSFORM_NORMALIZED_GCT.R
#Document Control 
#Latest Update: 2018-11-06
#This script inputs a normalized set of gene count tables, sets the 3rd quartile value of each profile to 1000, then applies a log_2 transformation.
#The resulting file can be input into the OneRNA Gene Expression Viewer

ReQ_packages = c("matrixStats", "optparse")

for (pack in 1:length(ReQ_packages)) {
  if(ReQ_packages[pack] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_packages[pack], repos='http://cran.us.r-project.org')
  }
}

library(matrixStats)
library(optparse)

input_dir <- "test_normalized.csv"
#output_path <- "C:\\Users\\vn\\OneDrive\\Documents\\Projects\\fusion_genes\\out\\"
main <- function(input_dir, output_path){
  #Read in normalized data file, use first row as column names
  df = read.csv(input_dir, header = TRUE, as.is = TRUE, stringsAsFactors = FALSE, strip.white = TRUE, sep='\t', check.names=FALSE)
  print(paste("Reading in", input_dir, "and transforming data"))
  
  df = df[!duplicated(df[,1])&is.na(df[,1])==F,]
  rownames(df) <- df[,1]
  df = df[,-1]
  
  num_cols <- unlist(lapply(df, is.numeric))
  
  if (length(num_cols) != ncol(df)){
    stop("ERROR: Type Error. Please check that columns contain only numeric values")
  }
  
  #Convert data columns to class numeric in preparation for scalar multiplication
  df[, num_cols] <- as.numeric(as.matrix(df[, num_cols]))
  
  Q3.mean <- get_Q3_mean(df[,num_cols])
  
  #Set Q3 to 1000 and add 1 (to avoid problems with taking log)
  df[, num_cols] <-df[, num_cols]/Q3.mean*1000+1
  
  #Confirm that the transformation matches, value should be 1001
  #Q3.mean_confirm <- get_Q3_mean(df[,num_cols])
  
  #Apply Log2 transformation
  df[, num_cols] <- log2(df[, num_cols])
  
  #Replace "NA" with 0 in expression data
  df[, num_cols]<- replace(df[, num_cols], is.na(df[, num_cols]), 0)
  
  #Replace periods with dashes
  names(df)[num_cols] <- gsub("\\.","\\-",names(df)[num_cols])
  names(df)[num_cols] <- gsub("$","_log2",names(df)[num_cols])
  
  #Output CSV file
  print(paste("Finished! Writing output to", paste(output_path)))
  write.csv(df, file=paste(output_path), row.names = TRUE)
}

get_Q3_mean <- function(df){
  #Calculate third quartile value for each normalized expression profile
  #Convert to data to  matrix first
  df.matrix <- data.matrix(df)
  Quartile <- colQuantiles(df.matrix)
  Q3.mean <- mean(Quartile[,"75%"])
  return(Q3.mean)
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default=NULL,
              help = "OneRNA normalized gene count table"),
  make_option(c("-o", "--output"), type = "character", default=NULL,
              help = "Output file to write scaled and transformed gene count table to")
)

opt <- parse_args(OptionParser(option_list = option_list))

main(opt$input, opt$output)
