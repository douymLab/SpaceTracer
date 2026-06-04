# Date: 2023/06/01
# Update: 2024/07/25
# Update: 2025/04/02
# Author: Qing Yang
# Work: Input the original rawdata and complete data preprocessing, which main purpose is extracting somatic posterior matrix, reads matrix and features dataframe.


##### Time #####
library(stringr)
start_time <- Sys.time()
print(str_c("Start time: ", start_time))


##### library
library(optparse)
option_list = list(
    make_option(c("-i", "--inputfile"), type = "character", default = NA, help = "the name of inputfile - posterior raw data & reads information"),
    make_option(c("-n", "--cellnum"), type = "integer", default = 22, help = "Cell number."),
    make_option(c("-o", "--outputpath"), type = "character", default = "data/", help = "The outputpath you can set"),
    make_option(c("-c", "--scid_file"), type = "character", default = "no", help = "Please enter the file containing the cellid used to calculate the genotype (col: sampleid) and the cellid consistent with the ground-truth tree (col: scid_basedTree). If you do not need to replace the cell order, please ignore this parameter directly."),
    make_option(c("-r", "--is_remove_cells"), type = "character", default = "yes", help = "Choice option (yes/no). If in order to validate whether the mutation can pass through the phylogenetic tree, it is recommended to set 'no' to perform the operation of not deleting; In order to build a phylogenetic tree from scratch, it is recommended to remove irrelevant cells."),
    make_option(c("-t", "--threshold"), type = "character", default = 0.9, help = "The shreshold of somatic posterior filteration."), 
    make_option(c("-s", "--indid"), type = "character", default = "UMB1465", help = "The individual id of your input.")
)
parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

library(purrr)
library(ggplot2)


##### Parameters
inputfile <- as.character(opt$inputfile)
cellnum <- as.numeric(as.character(opt$cellnum))
outputpath <- as.character(opt$outputpath)
dir.create(outputpath, recursive=TRUE)
scid_file <- as.character(opt$scid_file)
is_remove_cells <- as.character(opt$is_remove_cells)
valid_choices <- c("yes", "no")
if (!is_remove_cells %in% valid_choices) {
  stop("ParaError: --is_remove_cells must be 'yes' or 'no'")
}
cutoff <- as.numeric(as.character(opt$threshold))
indid <- as.character(opt$indid)

# inputfile <- "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_A549/PhyloSOLID_tree/mosaic_mutations/treeinput_100k/treeinput_100k_spot_c_2162.csv"
# cellnum <- 2162
# outputpath <- "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_A549/PhyloSOLID_tree/tree_valid"
# dir.create(outputpath, recursive=TRUE)
# scid_file <- "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Benchmark/data_A549/PhyloSOLID_tree/mosaic_mutations/treeinput_100k/treeinput_100k_scid_barcode.txt"
# is_remove_cells <- "no"
# valid_choices <- c("yes", "no")
# if (!is_remove_cells %in% valid_choices) {
#   stop("ParaError: --is_remove_cells must be 'yes' or 'no'")
# }
# cutoff <- 0.5
# indid <- "A549_100k"


# Display the parameters for verification
cat("Parameters:\n")
cat("Input file: ", inputfile, "\n")
cat("Cell number: ", cellnum, "\n")
cat("Output path: ", outputpath, "\n")
cat("SCID file: ", scid_file, "\n")
cat("Remove cells: ", is_remove_cells, "\n")
cat("Threshold: ", cutoff, "\n")
cat("Individual ID: ", indid, "\n")


##### Functions
get.mut_allele <- function(mut_type){
    str <- mut_type
    split <- unlist(strsplit(str, "[ ,()]"))
    gt_pairs <- split[split != ""]
    ref_gt <- gt_pairs[1:2]
    mut_gt <- gt_pairs[3:4]
    if(length(unique(mut_gt))==1){
        mut_allele <- unique(mut_gt)
    } else {
        mut_allele <- mut_gt[!mut_gt %in% ref_gt]
    }
    return(gsub("'", "", mut_allele))
}

get_list.base_info <- function(info_set){
    set_str <- gsub("\\{|\\}", "", info_set)
    key_value_pairs <- strsplit(set_str, ", ")[[1]]
    my_list <- list()
    for(pair in key_value_pairs){
      parts <- strsplit(pair, ": ")[[1]]
      key <- parts[1]
      value <- as.numeric(parts[2])
      my_list[[key]] <- value
    }
    return(my_list)
}

stat.allele_info <- function(info_set, mut_allele){
    mutant_idx <- grep(mut_allele, names(get_list.base_info(info_set)))
    mut_allele_count <- 0
    for(idx in mutant_idx){
        mut_allele_count <- mut_allele_count + get_list.base_info(info_set)[[names(get_list.base_info(info_set))[idx]]]
    }
    total_allele_count <- sum(unlist(get_list.base_info(info_set)))
    allele_info <- as.character(str_c(mut_allele_count, "/", total_allele_count))
    return(allele_info)
}

mean_non_zero <- function(x){
    average_non_zero <- mean(x[x != 0], na.rm = TRUE)
    return(average_non_zero)
}

likelihood_split <- function(str_posterior){
    list <- as.numeric(strsplit(strsplit(strsplit(str_posterior, "[", fixed=T)[[1]][2], "]", fixed=T)[[1]][1], ", ")[[1]])
    return(list)
}

reads_split <- function(str_reads){
    list <- as.character(strsplit(strsplit(strsplit(str_reads, "[", fixed=T)[[1]][2], "]", fixed=T)[[1]][1], ", ")[[1]])
    return(list)
}


##########################################################
##### Read in datafile
inputdata <- read.table(inputfile, header=FALSE, sep=","); dim(inputdata)
# [1]   142 13283
inputdata <- unique(inputdata); dim(inputdata)


##### ID reference (sc).
# You can select wether to replace the cell id.
if(scid_file=="no"){
    scid_data <- as.data.frame(cbind(scid_basedTree=paste0("sc" ,1:cellnum), scid=paste0("sc" ,1:cellnum)))
} else {
    scid_data <- as.data.frame(read.table(scid_file, header=TRUE))
    if(dim(scid_data)[2]==1){
        scid_data <- as.data.frame(scid_data[scid_data$scid_basedTree!="none",])
        colnames(scid_data) <- "scid_basedTree"
    } else {
        scid_data <- scid_data[scid_data$scid_basedTree!="none",]
    }
}
dim(scid_data)
# [1] 3319    1


##### Posterior data
posterior_extract_data <- inputdata[,c(1:(5+cellnum+1))]; dim(posterior_extract_data)
# [1]  142 3325
posterior_data <- as.data.frame(cbind(mutid=paste0(posterior_extract_data[,2], "_", posterior_extract_data[,3], "_", posterior_extract_data[,4], "_", posterior_extract_data[,5]), posterior_extract_data[,1:dim(posterior_extract_data)[2]])); dim(posterior_data)
# [1] 1000   27
rownames(posterior_data) <- 1:dim(posterior_data)[1]
colnames(posterior_data) <- c("mutid", "indid", "chr", "pos", "ref", "mut", "somatic_posterior_persite", scid_data$scid_basedTree); dim(posterior_data)
# [1]  142 3326
posterior_data[1:2, ]


##### Likelihood_mut
likelihood_mut_extract_data <- inputdata[,c(2:5,(7+cellnum*1):(7+cellnum*2-1))]; dim(likelihood_mut_extract_data)
# [1]  142 3323
likelihood_mut_data <- as.data.frame(cbind(mutid=paste0(likelihood_mut_extract_data[,1], "_", likelihood_mut_extract_data[,2], "_", likelihood_mut_extract_data[,3], "_", likelihood_mut_extract_data[,4]), likelihood_mut_extract_data[,5:(5+cellnum-1)])); dim(likelihood_mut_data)
# [1]  142 3320
colnames(likelihood_mut_data) <- c("mutid", paste0(scid_data$scid_basedTree, "_mut")); dim(likelihood_mut_data)
# [1]  142 3320
print("***************************likelihood_mut_data")
print(likelihood_mut_data[1:2, ])
likelihood_mut_data[1:2, ]


##### Likelihood_unmut
likelihood_unmut_extract_data <- inputdata[,c(2:5,(7+cellnum*2):(7+cellnum*3-1))]; dim(likelihood_unmut_extract_data)
# [1]  142 3323
likelihood_unmut_data <- as.data.frame(cbind(mutid=paste0(likelihood_unmut_extract_data[,1], "_", likelihood_unmut_extract_data[,2], "_", likelihood_unmut_extract_data[,3], "_", likelihood_unmut_extract_data[,4]), likelihood_unmut_extract_data[,5:(5+cellnum-1)])); dim(likelihood_unmut_data)
# [1]  142 3320
colnames(likelihood_unmut_data) <- c("mutid", paste0(scid_data$scid_basedTree, "_unmut")); dim(likelihood_unmut_data)
# [1]  142 3320
likelihood_unmut_data[1:2, ]


##### Reads info data => Allele_info (mutant allele count/total allele count)
reads_extract_data <- inputdata[,c(2:5,(7+cellnum*3):(7+cellnum*4))]; dim(reads_extract_data)
# [1]  142 3324
allele_data <- as.data.frame(cbind(mutid=paste0(reads_extract_data[,1], "_", reads_extract_data[,2], "_", reads_extract_data[,3], "_", reads_extract_data[,4]), reads_extract_data[,5:(5+cellnum)])); dim(allele_data)
# [1]  142 3321
colnames(allele_data) <- c("mutid", paste0(c("bulk", scid_data$scid_basedTree), "_AlleleStat")); dim(allele_data)
# [1]  142 3321
allele_data[1:2, ]


##### Merge posterior and read into a df
likelihoods_all_merged_data <- merge(likelihood_unmut_data, likelihood_mut_data, by="mutid"); dim(likelihoods_all_merged_data)
# [1]  142 6639
posterior_likelihoods_merged_data <- merge(posterior_data, likelihoods_all_merged_data, by="mutid"); dim(posterior_likelihoods_merged_data)
# [1]  142 9964
posterior_likelihoods_allele_merged_data <- merge(posterior_likelihoods_merged_data, allele_data, by="mutid"); dim(posterior_likelihoods_allele_merged_data)
# [1]   142 13284
all_merged_data <- posterior_likelihoods_allele_merged_data; dim(all_merged_data)
# [1]   142 13284
all_merged_data[all_merged_data == "0/0"] <- NA
all_merged_data[1,]
#              mutid   indid chr     start       end
# 1 Human_STR_273180 UMB1465  11 132566633 132566655
#                                       ref
# 1 (CAGAG, ACACACACACACACAGACACACA, TGAAG)
#                                         mut somatic_posterior_persite
# 1 (CAGAG, ACACACACACACACACAGACACACA, TGAAG)                         1
#                   sc8                sc4                sc3                 sc6
# 1 0.00461613160443611 0.0266195303674471 0.0542730468078042 0.00613897244843555
#                    sc1                 sc9 sc18 sc14 sc20 sc15
# 1 0.000733986273101432 0.00202139104113102    1    1    1    1
#                  sc13                sc7              sc10_2               sc5
# 1 0.00538142353909083 0.0017261808444772 0.00106395744009769 0.000830967584408
#                  sc12 sc19                sc11                 sc2 sc16_1
# 1 0.00186734980578539    1 0.00119335650485626 0.00225262155169897      1
#   sc16_2               sc10_1 sc17         sc8_unmut         sc4_unmut
# 1      1 0.000340830673138674    1 -4.51158127696236 -2.96457062187094
#           sc3_unmut        sc6_unmut         sc1_unmut         sc9_unmut
# 1 -1.96796634553353 -4.2229931990513 -6.46348636527646 -5.35100562583719
#          sc18_unmut        sc14_unmut        sc20_unmut        sc15_unmut
# 1 -49.9286692294588 -89.5551579017285 -3.19443436256129 -21.3261282314306
#          sc13_unmut         sc7_unmut      sc10_2_unmut         sc5_unmut
# 1 -4.38356424824702 -5.53482059988318 -6.00171372737831 -6.27616666274987
#          sc12_unmut        sc19_unmut        sc11_unmut        sc2_unmut
# 1 -5.43038133378075 -37.6674029836497 -5.91486898513716 -5.2434666370342
#        sc16_1_unmut      sc16_2_unmut      sc10_1_unmut        sc17_unmut
# 1 -35.6229815064425 -46.5480699086992 -7.28032610894618 -27.4039006510445
#         sc8_mut_by_BB      sc4_mut_by_BB       sc3_mut_by_BB
# 1 -0.0420220140110085 -0.269452917358134 -0.0140523371878356
#         sc6_mut_by_BB      sc1_mut_by_BB       sc9_mut_by_BB    sc18_mut_by_BB
# 1 -0.0400653458610075 -0.151212547924612 -0.0530718501975934 -160.206436211876
#      sc14_mut_by_BB    sc20_mut_by_BB    sc15_mut_by_BB      sc13_mut_by_BB
# 1 -408.404785976849 -80.0636756765025 -64.0900053867294 -0.0681700202366787
#         sc7_mut_by_BB    sc10_2_mut_by_BB       sc5_mut_by_BB
# 1 -0.0787165150195164 -0.0610305203485956 -0.0880903618735172
#        sc12_mut_by_BB    sc19_mut_by_BB      sc11_mut_by_BB      sc2_mut_by_BB
# 1 -0.0530275176804273 -195.589082849462 -0.0890898622066007 -0.054073350531677
#    sc16_1_mut_by_BB  sc16_2_mut_by_BB   sc10_1_mut_by_BB   sc17_mut_by_BB
# 1 -228.911834103525 -213.917408582092 -0.200554391893868 -674.84446427795
#       sc8_mut_by_prod    sc4_mut_by_prod     sc3_mut_by_prod
# 1 -0.0420220140110085 -0.269452917358134 -0.0140523371878356
#       sc6_mut_by_prod    sc1_mut_by_prod     sc9_mut_by_prod  sc18_mut_by_prod
# 1 -0.0400653458610075 -0.151212547924612 -0.0530718501975934 -160.206436211876
#    sc14_mut_by_prod  sc20_mut_by_prod  sc15_mut_by_prod    sc13_mut_by_prod
# 1 -408.404785976849 -80.0636756765025 -64.0900053867294 -0.0681700202366787
#       sc7_mut_by_prod  sc10_2_mut_by_prod     sc5_mut_by_prod
# 1 -0.0787165150195164 -0.0610305203485956 -0.0880903618735172
#      sc12_mut_by_prod  sc19_mut_by_prod    sc11_mut_by_prod    sc2_mut_by_prod
# 1 -0.0530275176804273 -195.589082849462 -0.0890898622066007 -0.054073350531677
#   sc16_1_mut_by_prod sc16_2_mut_by_prod sc10_1_mut_by_prod sc17_mut_by_prod
# 1  -228.911834103525  -213.917408582092 -0.200554391893868 -674.84446427795
#   sc8_AlleleStat sc4_AlleleStat sc3_AlleleStat sc6_AlleleStat sc1_AlleleStat
# 1           1/42           0/87           0/53          1/111           0/78
#   sc9_AlleleStat sc18_AlleleStat sc14_AlleleStat sc20_AlleleStat
# 1           0/67            0/48            0/44            0/79
#   sc15_AlleleStat sc13_AlleleStat sc7_AlleleStat sc10_2_AlleleStat
# 1            0/67           0/102           0/52              0/41
#   sc5_AlleleStat sc12_AlleleStat sc19_AlleleStat sc11_AlleleStat sc2_AlleleStat
# 1           0/34            0/12            0/78           0/106           0/65
#   sc16_1_AlleleStat sc16_2_AlleleStat sc10_1_AlleleStat sc17_AlleleStat
# 1             50/75            75/106              0/62            0/87


######################################################################
##### Filter cells that do not contain any mutations
if(is_remove_cells == "yes") {
  
  remained_scid <- c()
  removed_scid <- c()
  
  for(sc in scid_data$scid_basedTree) {
    sc_values <- all_merged_data[, sc]
    
    # 如果所有值都小于 cutoff 或者全是 NA，则移除
    if(all(sc_values < cutoff, na.rm = TRUE) || all(is.na(sc_values))) {
      removed_scid <- c(removed_scid, sc)
    } else {
      # 找出大于 cutoff 的位点
      site_pass <- which(!is.na(sc_values) & sc_values > cutoff)
      
      if(length(site_pass) == 0) {
        removed_scid <- c(removed_scid, sc)
      } else {
        # 提取 mutant allele count
        reads_list <- all_merged_data[site_pass, paste0(sc, "_AlleleStat")]
        mutant_dp <- sapply(reads_list, function(x) as.numeric(strsplit(x, "/")[[1]][1]))
        names(mutant_dp) <- rownames(all_merged_data)[site_pass]
        
        # 如果有 mutant allele count >=1，则保留
        if(any(mutant_dp >= 1, na.rm = TRUE)) {
          remained_scid <- c(remained_scid, sc)
        } else {
          removed_scid <- c(removed_scid, sc)
        }
      }
    }
  }
  
  # ===========================
  # 更新 merged_data，只保留存在的列
  # ===========================
  base_cols <- c("mutid", "indid", "chr", "pos", "ref", "mut", "somatic_posterior_persite")
  scid_cols <- remained_scid
  scid_unmut_cols <- paste0(remained_scid, "_unmut")
  scid_BB_cols <- paste0(remained_scid, "_mut_by_BB")
  scid_prod_cols <- paste0(remained_scid, "_mut_by_prod")
  scid_AlleleStat_cols <- paste0(remained_scid, "_AlleleStat")
  
  all_cols <- c(base_cols, scid_cols, scid_unmut_cols, scid_BB_cols, scid_prod_cols, scid_AlleleStat_cols)
  
  # 只保留实际存在的列
  cols_to_keep <- intersect(all_cols, colnames(all_merged_data))
  
  merged_data <- all_merged_data[, cols_to_keep]
  
  # 输出信息
  print(str_c(
    "Raw cell number is: ", length(scid_data$scid_basedTree), 
    "; After filtering cells that do not contain significant mutations, current cell number is: ", 
    length(remained_scid)
  ))
  
} else if(is_remove_cells == "no") {
  merged_data <- all_merged_data
  remained_scid <- scid_data$scid_basedTree
  print(str_c(
    "Raw cell number is: ", length(scid_data$scid_basedTree), 
    "; And no cells are removed!"
  ))
}

# 输出 merged_data 形状
print(str_c(
  "The shape of merged_data is : ", 
  dim(merged_data)[1], " rows X ", 
  dim(merged_data)[2], " columns"
))
# [1] "The shape of merged_data is : 142 rows X 13284 columns"


##### Calculate average mutant cells allele frequency
### Step1: bulk features
cov_in_pseudobulk <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[2]))))})
### Step2
avg_cov_percell <- apply(merged_data, 1, function(row) {mean_non_zero(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[2]))))})
### => step3: Mutant cells are defined as those that contain at least one mutant allele.
# Detect if mutant allele is present and calculate VAF cross all cells
avg_mutAF_percell <- apply(merged_data, 1, function(row) {mean_non_zero(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))))})
max_mutAF_percell <- apply(merged_data, 1, function(row) {max(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))))})
# The coverage of the max mutant AF cell
cov_in_maxmutAFcell <- apply(merged_data, 1, function(row) {
    alleles <- na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")]))
    if (length(alleles) == 0) {
        return(NA)  # 如果所有 AlleleStat 都是 NA，则返回 NA
    }
    max_ratio <- max(sapply(strsplit(alleles, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
    max_value_index <- which(sapply(strsplit(alleles, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]) == max_ratio), arr.ind = TRUE)[1]
    max_value_after_slash <- as.numeric(strsplit(alleles[max_value_index], "/")[[1]][2])
    return(max_value_after_slash)
})
# Count the maximum allele number cross all cells
total_mutAllele_num <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1]))))})
max_mutAllele_num <- apply(merged_data, 1, function(row) {max(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1]))))})
mutAllele_cellnum <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(strsplit(na.omit(unlist(row[paste0(remained_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])))>0)})
### => Step4: identify mutant cell by somatic_posterior_percell posterior more than a set value
# Count shared cell number by cell posterior value
moreCutoff_cellnum <- apply(merged_data, 1, function(row) {sum(as.numeric(sapply(na.omit(unlist(row[remained_scid])), function(x) as.numeric(x)>cutoff)))})
### => step5: Consider the number of cells with cutoff and with or without mutant allele
mutant_cellnum <- apply(merged_data, 1, function(row) {
    consider_scid <- remained_scid[!is.na(unlist(row[paste0(remained_scid, "_AlleleStat")])) & !is.na(unlist(row[remained_scid]))]
    sum(
        sapply(strsplit(na.omit(unlist(row[paste0(consider_scid, "_AlleleStat")])), "/"), function(x) as.numeric(x[1])) > 0 &
        sapply(na.omit(unlist(row[consider_scid])), function(x) as.numeric(x)) >= cutoff
    )
})
mutant_cell_fraction <- mutant_cellnum/length(remained_scid)
# out to dataframe
df_calAF_data <- as.data.frame(cbind(merged_data[,1:8], cov_in_pseudobulk, avg_cov_percell, avg_mutAF_percell, max_mutAF_percell, cov_in_maxmutAFcell, max_mutAllele_num, total_mutAllele_num, mutant_cell_fraction, mutant_cellnum, mutAllele_cellnum, moreCutoff_cellnum, merged_data[9:dim(merged_data)[2]])); dim(df_calAF_data)
# [1]   142 13294
filter_calAF_data <- df_calAF_data[df_calAF_data$somatic_posterior_persite>=cutoff, ]; dim(filter_calAF_data)
# [1]   142 13294
filter_calAF_data <- filter_calAF_data[filter_calAF_data$total_mutAllele_num!=0, ]; dim(filter_calAF_data)
# [1]   142 13294
# First output: Smmary_data (including only present in one cells)
write.table(filter_calAF_data, str_c(outputpath, "/Summary_data_", as.character(dim(filter_calAF_data)[1]), ".filtered_and_calculatedAF_and_mutantAnno.txt"), sep="\t", quote=FALSE, col.names=T, row.names=F)


##### Extract somatic_posterior_persite and posterior file from raw calculated results
final_data <- filter_calAF_data; dim(final_data)
# [1]   142 13294
# Following output: (except only present in one cells)

### 1. posterior
output.posterior_data <- as.data.frame(lapply(final_data[, remained_scid], as.numeric)); dim(output.posterior_data)
# [1]  142 3319
rownames(output.posterior_data) <- final_data$mutid
colnames(output.posterior_data) <- remained_scid
write.table(output.posterior_data, str_c(outputpath, "/data.posterior_matrix.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 2. likelihood_unmut
output.likelihood_unmut_data <- as.data.frame(lapply(final_data[, paste0(remained_scid, "_unmut")], as.numeric)); dim(output.likelihood_unmut_data)
# [1]  142 3319
rownames(output.likelihood_unmut_data) <- final_data$mutid
colnames(output.likelihood_unmut_data) <- remained_scid
write.table(output.likelihood_unmut_data, str_c(outputpath, "/data.likelihood_unmut_matrix.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 3. likelihood_mut_by_BB
print(final_data)
output.likelihood_mut_by_BB_data <- as.data.frame(lapply(final_data[, paste0(remained_scid, "_mut")], as.numeric)); dim(output.likelihood_mut_by_BB_data)
# [1]  142 3319
rownames(output.likelihood_mut_by_BB_data) <- final_data$mutid
colnames(output.likelihood_mut_by_BB_data) <- remained_scid
write.table(output.likelihood_mut_by_BB_data, str_c(outputpath, "/data.likelihood_mut_matrix.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 5. allele_info
output.allele_data <- final_data[, paste0(c("bulk", remained_scid), "_AlleleStat")]; dim(output.allele_data)
# [1]  142 3320
rownames(output.allele_data) <- final_data$mutid
colnames(output.allele_data) <- c("bulk", remained_scid)
write.table(output.allele_data, str_c(outputpath, "/data.allele_count.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 6. features: somatic_posterior_persite + average mutant cell AF + avg_cov_percell
output.features <- as.data.frame(lapply(final_data[, c("somatic_posterior_persite", "cov_in_pseudobulk", "avg_cov_percell", "avg_mutAF_percell", "max_mutAF_percell", "cov_in_maxmutAFcell", "max_mutAllele_num", "total_mutAllele_num", "mutant_cell_fraction", "mutant_cellnum", "mutAllele_cellnum", "moreCutoff_cellnum")], as.numeric)); dim(output.features)
# [1] 142  11
output.features <- as.data.frame(cbind(mutid=final_data[, "mutid"], output.features)); dim(output.features)
# [1] 142  12
output.features <- merge(final_data[, c("mutid", "indid", "chr", "pos", "ref", "mut")], output.features, by="mutid", all=FALSE); dim(output.features)
# [1] 142  18
rownames(output.features) <- output.features$mutid
output.features <- output.features[, -1]; dim(output.features)
# [1] 142  17
write.table(output.features, str_c(outputpath, "/features.preprocess_items.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

### 7.sequencing data
mutant_cellIndex <- apply(final_data, 1, function(row) {which(na.omit(sapply(strsplit(unlist(row[paste0(remained_scid, "_AlleleStat")]), "/"), function(x) as.numeric(x[1]))>0 & sapply(na.omit(unlist(row[remained_scid])), function(x) as.numeric(x))>cutoff))})
fa_df <- data.frame(cellid=remained_scid)
for(pos in 1:length(mutant_cellIndex)){
    fa_vec_pos <- rep(final_data[pos,"ref"], length(remained_scid))
    cellIndex <- mutant_cellIndex[[pos]]
    fa_vec_pos[cellIndex] <- final_data[pos,"mut"]
    fa_df <- as.data.frame(cbind(fa_df, fa_vec_pos))
}
fa_df_only <- fa_df[,(1:length(mutant_cellIndex))+1]; dim(fa_df_only)
# [1] 22 76
rownames(fa_df_only) <- fa_df[,1]
colnames(fa_df_only) <- final_data[,"mutid"]
fa_str_list <- apply(fa_df_only, 1, function(x) {paste(x, collapse="")})
output.fasta <- data.frame(content=str_c(length(remained_scid), " ", length(mutant_cellIndex)))
for(i in 1:length(fa_str_list)){
    content_i <- str_c(names(fa_str_list[i]), " ", fa_str_list[i])
    output.fasta <- as.data.frame(rbind(output.fasta, content_i))
}
write.table(output.fasta, str_c(outputpath, "/data.allcells_fasta.phy"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


######################### plot #########################
##### posterior percell distributuion
plot_df <- as.data.frame(unlist(output.posterior_data)); dim(plot_df)
colnames(plot_df) <- "value"
plot_df <- na.omit(plot_df)
plot_df_stat <- as.data.frame(table(plot_df)); dim(plot_df_stat)
# [1] 8629    2
ratio <- sum(plot_df_stat[as.numeric(as.character(plot_df_stat$value))>0.25&as.numeric(as.character(plot_df_stat$value))<0.75, "Freq"])/sum(plot_df_stat$Freq)
# [1] 0.08128342
ratio <- sprintf("%.2f%%", ratio * 100)
p1 <- ggplot(plot_df, aes(x=value)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.posterior_percellpersite") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15)) +
    geom_text(x=0.5, y=500, label=str_c("0.25~0.75: ", as.character(ratio)), hjust=0, vjust=-1, size=5)
ggsave(str_c(outputpath, "/Histplot.somatic_posterior_percell.pdf"), p1, width=6, height=4)

##### Generate plot_features dataframe.
plot_features <- as.data.frame(cbind(mutid=rownames(output.features), output.features)); dim(plot_features)
plot_features <- na.omit(plot_features)

##### somatic posterior persite distributuion
p2 <- ggplot(plot_features, aes(x=somatic_posterior_persite)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.somatic_posterior_persite") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.somatic_posterior_persite.pdf"), p2, width=6, height=4)

##### average coverage distribution percell
p3 <- ggplot(plot_features, aes(x=avg_cov_percell)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.avg_cov_percell") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.avg_cov_percell.pdf"), p3, width=6, height=4)

##### mutant AF distributuion
p4 <- ggplot(plot_features, aes(x=avg_mutAF_percell)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.avg_mutAF_percell") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.avg_mutAF_percell.pdf"), p4, width=6, height=4)
p5 <- ggplot(plot_features, aes(x=max_mutAF_percell)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.max_mutAF_percell") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.max_mutAF_percell.pdf"), p5, width=6, height=4)

##### mutant cell number distributuion
p6 <- ggplot(plot_features, aes(x=mutAllele_cellnum)) +
    geom_histogram(bins=50, fill="steelblue", color="white") +
    labs(x="Value", y="Frequency", title="Histplot.mutAllele_cellnum") +
    theme_bw() + 
    theme(plot.title=element_text(hjust=0.5, size=15))
ggsave(str_c(outputpath, "/Histplot.mutAllele_cellnum.pdf"), p6, width=6, height=4)

##### mutant cell number barplot
mutAllele_cellnum_count <- table(plot_features$mutAllele_cellnum)/length(plot_features$mutAllele_cellnum)
mutAllele_cellnum_count_order <- mutAllele_cellnum_count[order(as.numeric(names(mutAllele_cellnum_count)), decreasing=T)]
mutAllele_cellnum_rate <- sprintf("%.2f%%", mutAllele_cellnum_count_order * 100)
x_position=c()
for(mut_cn in (unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])){
    position_list <- plot_features[plot_features$mutAllele_cellnum==mut_cn, "mutid"]
    if(mut_cn==max(unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])){
        position <- position_list[length(position_list)]
    } else if(mut_cn==min(unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])){
        position <- position_list[1]
    } else {
        position <- position_list[ceiling(length(position_list)/2)]
    }
    x_position <- c(x_position, position)
}
p7 <- ggplot(plot_features, aes(x=reorder(mutid, -mutAllele_cellnum), y=mutAllele_cellnum)) +
    geom_bar(stat="identity", color="steelblue") +
    labs(x="Mutations", y="The number of mutant cells", title="Histplot.mutAllele_cellnum") +
    theme_bw() + 
    theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(axis.title = element_text(size = 10)) + 
    theme(plot.title=element_text(hjust=0.5, size=15)) + 
    annotate("text", x=x_position, y=(unique(plot_features$mutAllele_cellnum)[order(unique(plot_features$mutAllele_cellnum), decreasing=T)])+0.1, label=mutAllele_cellnum_rate, color="black", fontface="bold", size=3)
ggsave(str_c(outputpath, "/Barplot.mutAllele_cellnum.pdf"), p7, width=6, height=4)

# ##### bulk converage distribution
# p8 <- ggplot(plot_features, aes(x=bulk_cov)) +
#     geom_histogram(bins=50, fill="steelblue", color="white") +
#     labs(x="Value", y="Frequency", title="Histplot.bulk_cov") +
#     theme_bw() + 
#     theme(plot.title=element_text(hjust=0.5, size=15))
# ggsave(str_c(outputpath, "/Histplot.bulk_cov.pdf"), p8, width=6, height=4)

# ##### bulk VAF distribution
# p9 <- ggplot(plot_features, aes(x=bulk_VAF)) +
#     geom_histogram(bins=50, fill="steelblue", color="white") +
#     labs(x="Value", y="Frequency", title="Histplot.bulk_VAF") +
#     theme_bw() + 
#     theme(plot.title=element_text(hjust=0.5, size=15))
# ggsave(str_c(outputpath, "/Histplot.bulk_VAF.pdf"), p9, width=6, height=4)


##### Time #####
end_time <- Sys.time()
print(str_c("End time: ", end_time))
print(str_c("Program finished in ", as.character(round((end_time-start_time), 4)), " seconds"))

