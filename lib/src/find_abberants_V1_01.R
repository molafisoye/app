rm(list=ls())
#Load requeired function libraries

ReQ_packages = c("ggplot2", "optparse", "reshape2", "plyr")

for (pack in 1:length(ReQ_packages)) {
  if(ReQ_packages[pack] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_packages[pack], repos='http://cran.us.r-project.org')
  }
}

suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(plyr))
suppressWarnings(library(optparse))


main <- function(gene_count_inputs_path, output_path, project_id){
  cases <- read.csv(gene_count_inputs_path, header = TRUE,stringsAsFactors = FALSE, check.names = FALSE)
  
  cases = cases[!duplicated(cases[,1])&is.na(cases[,1])==F,]
  rownames(cases) <- cases[,1]
  cases = cases[,-1]
  
  num_cols <- unlist(lapply(cases, is.numeric))
  
  # #Select genes of interest using an external file
  # #Read in gene panel list
  # OneRNA_panels.csv
  #genelist <- read.csv(genelist_path,  header = TRUE, as.is= TRUE, stringsAsFactors = FALSE, strip.white=TRUE)
  #cases.genelist <- cases[rownames(cases) %in% genelist[,1],num_cols]
  cases.genelist <- cases[,num_cols]
  
  #Reformat gene expression data for analysis
  cases.genelist <-t(cases.genelist)
  cases.genelist <-as.data.frame(cases.genelist)
  cases.genelist.melt <-melt(cases.genelist, measure.vars = 1:ncol(cases.genelist), variable.name="SYMBOL", value.name = "Expression")
  
  #COMPARE TO NORMAL TISSUES
  # Extract normal tissue data from all cases
  #breast
  norms <- cases.genelist[grepl('_NT_log2$|_NT$', rownames(cases.genelist)), ]
  tumors <- cases.genelist[!grepl('_NT_log2$|_NT$', rownames(cases.genelist)), ]

  #Convert data frame to long format for boxplot
  norms.m <- melt(norms, measure.vars = 1:ncol(norms), variable.name="SYMBOL", value.name = "Expression")
  tumor.m <-  melt(tumors, measure.vars = 1:ncol(norms), variable.name="SYMBOL", value.name = "Expression")
 
  #Compute boxplot stats for all genes in all"normals"
  A <- boxplot(Expression~SYMBOL, data=norms.m, plot = FALSE, range = 1.5)
  boxstats <- A$stats
  colnames(boxstats)<-A$names
  rownames(boxstats)<-c('Min','Q1','Median','Q3','Max')
  boxstats <-data.frame(t(boxstats))
  
  boxstats$IQR <-as.matrix(boxstats[1:nrow(boxstats),"Max"]-boxstats[1:nrow(boxstats), "Min"])
  boxstats$HIGH <- boxstats[1:nrow(boxstats), "Max"]
  boxstats$LOW <- boxstats[1:nrow(boxstats), "Min"]
  boxstats$VERY_HIGH <- boxstats[1:nrow(boxstats), "Median"]+2*boxstats[1:nrow(boxstats),"IQR"]
  boxstats$VERYLOW <- boxstats[1:nrow(boxstats), "Median"]-2*boxstats[1:nrow(boxstats),"IQR"]
  
  write.csv(boxstats, file=paste(output_path, "/", project_id,"_thresholds.csv", sep=""))
  
  #CASE SELECTION
  #Include all input cases in the analysis
  caselist  <- as.matrix(colnames(cases))
  
  # Load expression profiles for selected cases into a matrix
  exp_data <- as.matrix(t(cases.genelist[rownames(cases.genelist) %in% caselist,]))
  
  #Compile over and under expression status for all selected cases and genes
  num_genes <- ncol(cases.genelist)
  genes.over <-as.matrix(exp_data[1:num_genes, ]  > boxstats[1:num_genes, "Max"])
  genes.under <-as.matrix(exp_data[1:num_genes, ]  < boxstats[1:num_genes, "Min"])
  IQR <-as.matrix(boxstats[1:num_genes, "Max"]-boxstats[1:num_genes, "Min"])
  genes.high <-as.matrix(exp_data[1:num_genes, ]  > boxstats[1:num_genes, "Median"]+2*IQR[1:num_genes])
  genes.low <-as.matrix(exp_data[1:num_genes, ]  < boxstats[1:num_genes, "Median"]-2*IQR[1:num_genes])
  
  genes.states <- genes.high
  genes.states[] <- "NORMAL"
  genes.states[genes.over] <- "HIGH"
  genes.states[genes.under] <- "LOW"
  genes.states[genes.high] <- "VERY HIGH"
  genes.states[genes.low] <- "VERY LOW"
  
  ### melt gene states so that for each sample, we know aberrant expressed genes
  #genes.states.m <- melt(genes.states, varnames = c('SYMBOL', 'SAMPLE_ID'))
  #genes.stats.m.aberrant <- genes.states.m[genes.states.m$value != 'NORMAL',]
  
  # # #Save data to files
  write.csv(genes.states, paste(output_path,'/',project_id,"_select_gene_states.csv", sep=''))
  #write.csv(genes.stats.m.aberrant, paste(output_path,'/',project_id,"_aberrant_gene_states.csv", sep=''))
  
  #aberrant_counts <-cbind(t(rbind(colSums(genes.under, na.rm=T), colSums(genes.over, na.rm=T))))
  #aberrant_counts <-cbind(t(rbind(colSums(genes.low), colSums(genes.high))))
  #aberrant_counts <- as.data.frame(aberrant_counts)
  #aberrant_counts$NORMAL_SAMPLE <- 0
  #aberrant_counts[grepl('_NT_log2$|_NT$', rownames(aberrant_counts)),'NORMAL_SAMPLE'] <- 1
  #colnames(aberrant_counts) <- c("UNDER/LOW", "OVER/HIGH", "NORMAL")
  #write.csv(aberrant_counts, paste(output_path, '/',project_id,"_aberrant_counts.csv", sep = ''))
}

option_list <- list(
  make_option(c("-i", "--input_gct"), type = "character", default=NULL,
              help = "OneRNA normalized gene count table"),
  #make_option(c("-g", "--gene_list"), type = "character", default=NULL,
   #           help = "OneRNA genes of interest where the first column will be taken as the genes of interest. List must be greater than 1"),
  make_option(c("-o", "--output"), type = "character", default=NULL,
              help = "Output path"),
  make_option(c("-p", "--project_id"), type = "character", default=NULL,
              help = "project id to write files to")
)

opt <- parse_args(OptionParser(option_list = option_list))

main(opt$input_gct, opt$output, opt$project_id)
