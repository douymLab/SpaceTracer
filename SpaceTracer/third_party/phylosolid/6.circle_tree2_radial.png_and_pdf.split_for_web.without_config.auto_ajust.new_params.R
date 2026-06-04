# Date: 2023/06/01
# Update: 2024/12/09
# Author: Qing Yang, Yonghe Xia
# Work: Plot circos tree.


##### Time #####
library(stringr)
start_time <- Sys.time()
print(str_c("Start time: ", start_time))


##### environment and library
.libPaths(c("/home/douyanmeiLab/yangqing/R/x86_64-conda-linux-gnu-library/4.3", "/storage/douyanmeiLab/yangqing/anaconda/anaconda3/envs/pmg/lib/R/library"))

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  devtools::install_github("jokergoo/ComplexHeatmap")
}
library(ComplexHeatmap)
if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
  devtools::install_github("xiangpin/ggtreeExtra")
}
library(ggtreeExtra)
library(converTree)
library(ggnewscale)
library(phylogram)
library(circlize)
library(stringr)
library(ggplot2)
library(cowplot)
library(tibble)
library(treeio)
library(ggtree)
library(ggtext)
if (!requireNamespace("gsubfn", quietly = TRUE)) {
  install.packages('gsubfn')
}
library(gsubfn)
library(tidyr)
library(dplyr)
library(ape)
library(pheatmap)
library(Polychrome)
library(pals)
library(grid)
library(gtable)
library(gridExtra)
library(rlang)
library(paletteer)
library(patchwork)
library(ggforce)  # 为了绘制arc bar等


library(optparse)
option_list = list(
    make_option(c("--inputpath"), type="character", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Visualized/plot_circos/auto_adjust/phylo_P6_35muts_1234cells/",  help="The datapath including inputfile, cffile, features_file and total_flipping_file."),
    make_option(c("--outputpath"), type="character", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/Visualized/plot_circos/auto_adjust/phylo_P6_35muts_1234cells/plot_circos/",  help="The outputpath."),
    make_option(c("--annotation_file"), type="character", default="/storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/04.mutations/P6/tree_combine/combine_all_info_add_cellprop.txt",  help="The annotation file."),

    make_option(c("--target_mut"), type="character", default="no",  help="The target_mut"),

    make_option(c("--conflict_muts"), type="character", default=NULL,  help="The conflict_muts"),
    make_option(c("--selected_mutlist"), type="character", default="all",  help="The selected_mutlist"),
    make_option(c("--adding_mut_file"), type="character", default="no",  help="The adding_mut_file"),
    make_option(c("--manual_fp_file"), type="character", default="no",  help="The manual_fp_file"),
    make_option(c("--barcode_file"), type="character", default="no",  help="The barcode_file"),
    make_option(c("--barcode_name"), type="character", default="no",  help="The barcode_name")
    )
parseobj = OptionParser(option_list=option_list,
            usage="usage: %prog [options]",
            add_help_option=TRUE,
            description="Valid the converted tree and calculate the accuray of the phylogenetic tree.")
opt = parse_args(parseobj)
# config_file <- as.character(opt$config_file)

inputpath <- as.character(opt$inputpath)
outputpath <- as.character(opt$outputpath)
annotation_file <- as.character(opt$annotation_file)

target_mut <- as.character(opt$target_mut)
conflict_muts <- as.character(opt$conflict_muts)
selected_mutlist <- as.character(opt$selected_mutlist)
adding_mut_file <- as.character(opt$adding_mut_file)
manual_fp_file <- as.character(opt$manual_fp_file)
barcode_file <- as.character(opt$barcode_file)
barcode_name <- as.character(opt$barcode_name)


inputfile <- str_c(inputpath, "/final_cleaned_I_full_withNA3_for_circosPlot.txt")
cffile <- str_c(inputpath, "/final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix")
features_file <- str_c(inputpath, "/df_flipping_count_for_each_mut.txt")
total_flipping_file <- str_c(inputpath, "/df_total_flipping_count.txt")




###### manually setting and running
# inputpath <- "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_P6_rep1.80true/mutation_integrator/test_100k/"
# outputpath <- "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_P6_rep1.80true/mutation_integrator/test_100k/circos_plot"
# annotation_file <- "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/spatial/tree_rebuttal_2508/annofile/P6_combined_annotation_file.txt"
tip_label_offset=10
tip_label_size=2.5
tip_point_size=0.5
heatmap_circos_width=0.3
heatmap_circos_offset=0.04
flipping_point_size=0.2
plot_height=12
plot_width=18
# inputfile <- str_c(inputpath, "/final_cleaned_I_full_withNA3_for_circosPlot.txt")
# cffile <- str_c(inputpath, "/final_cleaned_M_lowR_basedPivots.filtered_sites_inferred.CFMatrix")
# features_file <- str_c(inputpath, "/df_flipping_count_for_each_mut.txt")
# total_flipping_file <- str_c(inputpath, "/df_total_flipping_count.txt")




##### help function ----------------------------------------------------------
create_unique_filename <- function(basename, extension) {
  counter <- 1
  unique_name <- paste0(basename, extension)

  # 检查文件是否存在，如果存在则添加一个数字后缀
  while(file.exists(unique_name)) {
    unique_name <- paste0(basename, "_", counter, extension)
    counter <- counter + 1
  }

  return(unique_name)
}

find_zero_mutation <- function(table){
  zero_mutation <- rownames(table)[rowSums(table) == 0]
  return(zero_mutation)
}

find_zero_barcode <- function(table){
  zero_barcode <- colnames(table)[colSums(table) == 0]
  return(zero_barcode)
}

# format tree with muts label and branch length = count(muts)
dat2tree <- function(dat) {

  count_muts <- function(dat) {
    tmp <- str_count(dat$label, fixed("|"))
    return(tmp+1)
  }

  mutation_count <- count_muts(dat)
  dat$branch.length <- mutation_count

  dat2=dat
  tree_ <- dat %>% ape::as.phylo()

  ## fix label info
  index_label <- function(node) {

    dat2[dat2$node==node,]$label

  }

  index_length <- function() {

    root_node = dat2[dat2$node==dat2$parent,]$parent
    dat3 <- dat2[dat2$node!=root_node,]
    dat3$branch.length

  }

  tree_[["edge.length"]] <- index_length()
  tree_[["tip.label"]] <- Map(index_label,tree_[["tip.label"]]) %>% unlist()
  tree_[["node.label"]] <- Map(index_label,tree_[["node.label"]]) %>% unlist()

  return(tree_)
}

## find subclone and clone of conflict-free matrix
subclone_finder <- function(cf_tableno0) {
  ## check cf_tableno0

  ## convert to long data
  cf_tableno0_long <- cf_tableno0 %>%
    gather(key = "mutation", value = "presence", -cellIDxmutID) %>%
    filter(presence == 1) %>%
    arrange(cellIDxmutID, mutation)

  # muts_factor factor variable
  muts_factor <- factor(cf_tableno0_long$mutation)

  # Calculate the count of each level
  level_counts <- table(muts_factor)

  # Order the levels by count
  ordered_levels <- names(sort(level_counts, decreasing = TRUE))

  # Reorder the factor
  muts_factor_reordered <- factor(muts_factor, levels = ordered_levels)

  # for convenience, add a column representing the mutations by number code
  cf_tableno0_long$mutation_num <- muts_factor_reordered %>%
    as.numeric()  %>%
    as.character()

  # Identifying subclones
  subclones <- cf_tableno0_long %>%
    group_by(cellIDxmutID) %>%
    summarize(subclone = paste(mutation, collapse = ","),
              subclone_num = paste(mutation_num, collapse = ",")) %>%
    group_by(subclone,subclone_num) %>%
    summarize(cells = paste(cellIDxmutID, collapse = ","))

  ## add stat of subclones
  subclones$subclone_size = subclones$subclone %>% str_count(",")+1

  # a dataset for testing
  # subclones <- data.frame(
  #   clone = c("clone1", "clone1", "clone2", "clone3", "clone3", "clone4", "clone5", "clone6", "clone7", "clone7"),
  #   subclone_num = c("4,8,11", "4,11", "4", "3,4,2", "4,2", "3,4,5,2,9", "4,6", "4,7", "4,1,10", "4,10")
  # )

  # 定义 is_superset 函数
  is_superset <- function(subclone_set, other_set) {
    all(other_set %in% subclone_set)
  }

  # 定义合并子集的函数
  merge_subsets <- function(clone_mapping) {
    to_remove <- rep(FALSE, length(clone_mapping))
    for (i in seq_along(clone_mapping)) {
      for (j in seq_along(clone_mapping)) {
        if (i != j && !to_remove[j] && is_superset(clone_mapping[[i]], clone_mapping[[j]])) {
          clone_mapping[[i]] <- unique(c(clone_mapping[[i]], clone_mapping[[j]]))
          to_remove[j] <- TRUE
        }
      }
    }
    return(clone_mapping[!to_remove])  # 移除被合并的元素
  }

  # 转换 subclone_int 为整数列表
  subclones$subclone_int <- lapply(subclones$subclone_num, function(x) sort(as.integer(unlist(strsplit(x, ",")))))

  # 创建克隆映射
  clone_mapping <- setNames(subclones$subclone_int, paste0("clone", seq_along(subclones$subclone_int)))

  # 合并子集
  clone_mapping <- merge_subsets(clone_mapping)

  # clone_mapping is the main clone finded
  # map to the subclones
  subclones$clone_list <- NA
  subclones$clone_list <- sapply(subclones$subclone_int, function(x) {
    matching_clone <- which(sapply(clone_mapping, function(y) is_superset(y, x)))
    if (length(matching_clone) == 0) {
      return(NA)
    } else {
      return(names(clone_mapping)[matching_clone])
    }
  })

  ## rename clone
  # Identify rows that appear multi clones by count ,
  subclones$count <- subclones$clone_list %>%
    as.character() %>%
    str_count(",") + 1

  # for rows contain multi clones, set as root clone
  # the top count is root clone
  # the second count is subroot clone
  # until the clone count is 1

  subclones$clone <- NA

  # find the clone count is 1
  subclones$clone[subclones$count == 1] <- unlist(subclones$clone_list[subclones$count == 1])

  # find the max one be root clone
  # if clone list contains all subclones, set as root clone
  contatins_all <- sapply(subclones$clone_list, function(x) {
    return(all(names(clone_mapping) %in% x))
  })
  subclones$clone[subclones$count != 1 & subclones$count == max(subclones$count) & contatins_all] <- "root_clone"

  # find the tiny one be subroot clone
  # for those clone list contains multiclone but not root clone
  # find is which clone is the subroot clone
  subroot_who <- sapply(subclones$clone_list[subclones$count != 1 & subclones$count < max(subclones$count)], function(x) {
    paste0("subroot:(",paste(x,collapse = ","),")")
  })
  subclones$clone[subclones$count != 1 & subclones$count < max(subclones$count)] <- subroot_who

  ## order subclones by their count and their relationship
  # subclones <- subclones[order(subclones$count),]
  ## is a subroot clone size is smaller and belong to another subroot clone
  ## they should be order together
  # subroot_list <- subclones$clone_list[subclones$count != 1 & subclones$count < max(subclones$count)]
  # names(subroot_list) <- subclones$count[subclones$count != 1 & subclones$count < max(subclones$count)]
  subroot_list <- subclones$clone_list
  names(subroot_list) <- rownames(subclones)

  # 创建包含关系矩阵
  contain_matrix <- matrix(FALSE, nrow = length(subroot_list), ncol = length(subroot_list))

  # 填充包含关系矩阵
  for (i in 1:length(subroot_list)) {
    for (j in 1:length(subroot_list)) {
      if (i != j && all(subroot_list[[i]] %in% subroot_list[[j]])) {
        contain_matrix[i, j] <- TRUE
      }
    }
  }

  # 计算每个元素的长度
  lengths <- sapply(subroot_list, length)

  # 排序逻辑：首先按是否被包含排序，然后按长度排序
  order_index <- order(-rowSums(contain_matrix), -lengths)

  # 应用排序
  sorted_list <- subroot_list[order_index]

  subclones <- subclones[names(sorted_list),]

  subclones$clone_order <- 1:nrow(subclones)

  # make a new name of clone
  clones <- subclones$clone

  # Extract numbers and handle special cases
  numbers <- gsub(".*?(\\d+)$", "\\1", clones)
  numbers <- as.numeric(numbers)  # Check that all numbers are numeric
  unique_numbers <- unique(na.omit(numbers))  # Exclude non-numeric parts like 'root_clone'

  # Map old numbers to new sequence
  new_numbers <- setNames(seq_along(unique_numbers), unique_numbers)

  match <- regexpr("\\d+$", clones, perl = TRUE)
  last_number <- regmatches(clones, match)

  subclones$clone[which(attr(match,"match.length") != -1)] = paste0("clone",new_numbers[last_number])

  # 对每个数字进行替换
  test = subclones$clone[subclones$count != 1 & subclones$count < max(subclones$count)]
  if (length(test)>0) {
    # Define a replacement function
    replacement_function <- function(x) paste0("clone", new_numbers[x])
    test2 <- gsubfn("clone(\\d+)", replacement_function, test)
    subclones$clone[subclones$count != 1 & subclones$count < max(subclones$count)] <- test2
  }


  # stat clone info
  subclones <- subclones %>%
    group_by(clone) %>%
    reframe(clone_size = sum(subclone_size),across(.cols = everything()))

  return(subclones)

}

remove_zero_barcode <- function() {
  # remove barcode with all zero mutations
  cffile_no0 <- paste0(cffile,".no0")
  # Read the file
  lines <- readLines(cffile)
  # Process the data
  filtered_lines <- suppressWarnings(lines[sapply(lines, function(line) sum(as.numeric(strsplit(line, "\t")[[1]][-1])) != 0 )])
  filtered_lines[1] <- lines[1]
  # Write the output
  writeLines(filtered_lines, cffile_no0)
  return(cffile_no0)
}

find_zero_barcode <- function(table){
  zero_barcode <- rownames(table)[rowSums(table) == 0]
  return(zero_barcode)
}

insert_newline <- function(text, n = 180) {
  # Splits the string into segments of the specified length
  split_text <- strsplit(text, NULL)[[1]]
  # Blocks the string to the specified length
  split_text <- split(split_text, ceiling(seq_along(split_text) / n))
  # Reassembles each block into a string and adds a newline
  result <- sapply(split_text, paste, collapse = "")
  # Combines the result into a complete string, concatenating segments with \n
  final_result <- paste(result, collapse = "\n")
  return(final_result)
}

# Check all paramters in input config file 
check_parameters <- function(params) {
  # Find all parameter names that are NULL or NA
  missing_params <- names(params)[sapply(params, function(x) is.null(x) || all(is.na(x)))]
  
  if (length(missing_params) > 0) {
    stop(paste("The following parameters are not assigned in the config_file: \n", paste(missing_params, collapse = ", ")))
  } else {
    message("All parameters have been assigned correctly in the config_file ! ")
  }
}


# sort！！！
sort_mutation_matrix <- function(df) {
  mat <- as.matrix(df)
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  
  # 在剩余的cells/mutations中找出最大的clone
  find_leader_clone <- function(available_items, by_row = TRUE) {
    if(length(available_items) == 0) return(NULL)
    
    if(by_row) {
      item_sums <- rowSums(mat[available_items, , drop = FALSE])
      leader_idx <- which.max(item_sums)
      leader <- available_items[leader_idx]
      members <- which(mat[leader, ] == 1)
    } else {
      item_sums <- colSums(mat[, available_items, drop = FALSE])
      leader_idx <- which.max(item_sums)
      leader <- available_items[leader_idx]
      members <- which(mat[, leader] == 1)
    }
    
    if(length(members) == 0) return(NULL)
    
    return(list(
      leader = leader,
      members = members
    ))
  }
  
  # 递归地对cells或mutations进行克隆划分
  recursive_clone_partition <- function(available_items, parent_members = NULL, by_row = TRUE) {
    if(length(available_items) == 0) return(NULL)
    
    if(!is.null(parent_members)) {
      if(by_row) {
        mat_subset <- mat[available_items, parent_members, drop = FALSE]
        row_sums <- rowSums(mat_subset)
        if(max(row_sums) == 0) return(NULL)
        
        leader_idx <- which.max(row_sums)
        leader <- available_items[leader_idx]
        members <- parent_members[which(mat_subset[leader_idx, ] == 1)]
      } else {
        mat_subset <- mat[parent_members, available_items, drop = FALSE]
        col_sums <- colSums(mat_subset)
        if(max(col_sums) == 0) return(NULL)
        
        leader_idx <- which.max(col_sums)
        leader <- available_items[leader_idx]
        members <- parent_members[which(mat_subset[, leader_idx] == 1)]
      }
    } else {
      leader_clone <- find_leader_clone(available_items, by_row)
      if(is.null(leader_clone)) return(NULL)
      
      leader <- leader_clone$leader
      members <- leader_clone$members
    }
    
    # 找出所有属于这个clone的items
    if(by_row) {
      # 对于cells，找出所有至少包含一个当前clone的mutation的cells
      clone_items <- available_items[rowSums(mat[available_items, members, drop = FALSE]) > 0]
      # 按照包含的mutation数量排序
      clone_items <- clone_items[order(rowSums(mat[clone_items, members, drop = FALSE]), 
                                     decreasing = TRUE)]
    } else {
      # 对于mutations，找出所有至少在一个当前clone的cell中出现的mutations
      clone_items <- available_items[colSums(mat[members, available_items, drop = FALSE]) > 0]
      # 按照出现的cell数量排序
      clone_items <- clone_items[order(colSums(mat[members, clone_items, drop = FALSE]), 
                                     decreasing = TRUE)]
    }
    
    # 创建当前clone的结构
    current_clone <- list(
      leader = leader,
      clone_items = clone_items,
      members = members
    )
    
    # 在完全不相关的items中寻找平行clone
    unused_items <- setdiff(available_items, clone_items)
    if(length(unused_items) > 0) {
      parallel_clone <- recursive_clone_partition(unused_items, by_row = by_row)
      if(!is.null(parallel_clone)) {
        current_clone$parallel <- parallel_clone
      }
    }
    
    return(current_clone)
  }
  
  # 从克隆结构中提取排序顺序
  extract_ordering <- function(clone_structure, by_row = TRUE) {
    ordered_items <- c()
    
    process_clone <- function(clone) {
      if(is.null(clone)) return()
      
      # 添加所有属于当前clone的items
      new_items <- setdiff(clone$clone_items, ordered_items)
      if(length(new_items) > 0) {
        ordered_items <<- c(ordered_items, new_items)
      }
      
      # 处理平行clone
      if(!is.null(clone$parallel)) {
        process_clone(clone$parallel)
      }
    }
    
    process_clone(clone_structure)
    return(ordered_items)
  }
  
  # 主排序流程
  # 1. 对cells进行递归克隆划分
  cell_clone_structure <- recursive_clone_partition(rownames(mat), by_row = TRUE)
  sorted_rows <- extract_ordering(cell_clone_structure)
  
  # 2. 对mutations进行递归克隆划分
  mut_clone_structure <- recursive_clone_partition(1:n_cols, by_row = FALSE)
  sorted_cols <- extract_ordering(mut_clone_structure)
  
  # 处理任何未分配的行和列
  remaining_rows <- setdiff(rownames(mat), sorted_rows)
  if(length(remaining_rows) > 0) {
    remaining_rows <- remaining_rows[order(
      -rowSums(mat[remaining_rows, , drop = FALSE])
    )]
    sorted_rows <- c(sorted_rows, remaining_rows)
  }
  
  remaining_cols <- setdiff(1:n_cols, sorted_cols)
  if(length(remaining_cols) > 0) {
    remaining_cols <- remaining_cols[order(
      -colSums(mat[, remaining_cols, drop = FALSE])
    )]
    sorted_cols <- c(sorted_cols, remaining_cols)
  }
  
  return(df[sorted_rows, colnames(df)[sorted_cols]])
}


get_used_colors <- function(p1_tree, other_colors) {
  # 获取所有用过的颜色
  tree_colors <- unique(p1_tree$data$fill)
  c(tree_colors, other_colors)
}

check_unique_color <- function(color, existing_colors) {
  if (color %in% existing_colors) {
    stop(paste("The fixed color", color, "is already in use! Please choose another color."))
  } else {
    message(paste("The fixed color", color, "is unique and will be used for the root node."))
  }
}

# 创建绘制 circos 的时候 fliiping 对应关系的 labels 对应表
format_flipping_label <- function(x) {
  x <- gsub("0", "non-mutant", x)
  x <- gsub("1", "mutant", x)
  x <- gsub("3", "missing", x)
  x <- gsub(">", "->", x)
  return(x)
}

trim_label <- function(label_str, target_mut, default_n = 2) {
  if (is.na(label_str) || label_str == "") return("")
  
  parts <- unlist(strsplit(label_str, "\\|"))
  
  if (target_mut %in% parts) {
    # 只显示 target_mut 这个mutation
    parts_trim <- target_mut
  } else {
    parts_trim <- head(parts, default_n)
  }
  
  paste(parts_trim, collapse = "\n")
}


# 递归克隆划分函数
recursive_clone_partition <- function(available_items, parent_members = NULL, by_row = TRUE) {
  if(length(available_items) == 0) return(NULL)
  
  if(!is.null(parent_members)) {
    if(by_row) {
      mat_subset <- mat[available_items, parent_members, drop = FALSE]
      row_sums <- rowSums(mat_subset)
      if(max(row_sums) == 0) return(NULL)
      
      leader_idx <- which.max(row_sums)
      leader <- available_items[leader_idx]
      members <- parent_members[which(mat_subset[leader_idx, ] == 1)]
    } else {
      mat_subset <- mat[parent_members, available_items, drop = FALSE]
      col_sums <- colSums(mat_subset)
      if(max(col_sums) == 0) return(NULL)
      
      leader_idx <- which.max(col_sums)
      leader <- available_items[leader_idx]
      members <- parent_members[which(mat_subset[, leader_idx] == 1)]
    }
  } else {
    # 初始寻找leader
    if(by_row) {
      item_sums <- rowSums(mat[available_items, , drop = FALSE])
      leader_idx <- which.max(item_sums)
      leader <- available_items[leader_idx]
      members <- which(mat[leader, ] == 1)
    } else {
      item_sums <- colSums(mat[, available_items, drop = FALSE])
      leader_idx <- which.max(item_sums)
      leader <- available_items[leader_idx]
      members <- which(mat[, leader] == 1)
    }
  }
  
  # 找出属于当前 clone 的 items
  if(by_row) {
    clone_items <- available_items[rowSums(mat[available_items, members, drop = FALSE]) > 0]
    clone_items <- clone_items[order(rowSums(mat[clone_items, members, drop = FALSE]), decreasing = TRUE)]
  } else {
    clone_items <- available_items[colSums(mat[members, available_items, drop = FALSE]) > 0]
    clone_items <- clone_items[order(colSums(mat[members, clone_items, drop = FALSE]), decreasing = TRUE)]
  }
  
  current_clone <- list(
    leader = leader,
    clone_items = clone_items,
    members = members
  )
  
  unused_items <- setdiff(available_items, clone_items)
  if(length(unused_items) > 0) {
    parallel_clone <- recursive_clone_partition(unused_items, by_row = by_row)
    if(!is.null(parallel_clone)) {
      current_clone$parallel <- parallel_clone
    }
  }
  
  return(current_clone)
}

# 递归提取 mutation 顺序
extract_mutation_order <- function(clone_structure) {
  ordered_mutations <- c()
  
  process_clone <- function(clone) {
    if (is.null(clone)) return()
    
    # 加入当前 clone 的 mutations
    new_items <- setdiff(clone$clone_items, ordered_mutations)
    if(length(new_items) > 0) {
      ordered_mutations <<- c(ordered_mutations, new_items)
    }
    
    # 递归平行 clone
    if (!is.null(clone$parallel)) {
      process_clone(clone$parallel)
    }
  }
  
  process_clone(clone_structure)
  return(ordered_mutations)
}


##### auto_adjust_plot_parameters ----------------------------------------------------------

predict_plot_params <- function(num_mutations, num_cells, num_anno_circos = 1) {
  # 手动训练集（去掉 N47）
  training <- data.frame(
    num_mutations = c(35, 29, 52, 11, 10, 65, 51),
    num_cells = c(1234, 1178, 1147, 256, 148, 2187, 1115),
    num_anno_circos = c(5, 2, 1, 1, 1, 2, 1),
    tip_label_offset = c(11, 7.5, 13.5, 2.5, 1.9, 2.6, 17),
    tip_label_size = c(1.5, 1.5, 1.5, 4, 4.5, 0.5, 1.2),
    tip_point_size = c(0.1, 0.3, 0.3, 2.5, 3, 0.005, 0.2),
    heatmap_circos_width = c(0.4, 0.3, 0.4, 0.4, 0.4, 0.4, 0.5),
    heatmap_circos_offset = c(0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04),
    flipping_point_size = c(0.5, 0.5, 0.4, 3, 3.5, 0.03, 0.4),
    plot_height = rep(12, 7),
    plot_width = rep(18, 7)
  )

  # 拟合函数（强制最小值为 0.01）
  fit_and_predict <- function(y) {
    model <- lm(y ~ num_mutations + num_cells + num_anno_circos, data = training)
    pred <- predict(model, newdata = data.frame(
      num_mutations = num_mutations,
      num_cells = num_cells,
      num_anno_circos = num_anno_circos
    ))
    return(max(pred, 0.01))  # 防止负值或零值
  }

  list(
    tip_label_offset = fit_and_predict(training$tip_label_offset),
    tip_label_size = fit_and_predict(training$tip_label_size),
    tip_point_size = fit_and_predict(training$tip_point_size),
    heatmap_circos_width = fit_and_predict(training$heatmap_circos_width),
    heatmap_circos_offset = 0.04,  # 固定值
    flipping_point_size = fit_and_predict(training$flipping_point_size),
    plot_height = 12,  # 固定值
    plot_width = 18    # 固定值
  )
}




##############################################################################
##############################################################################
##############################################################################

##### plot data preprocess ----------------------------------------------------------


### create prefix
pdf_lastfix <- paste0(format(Sys.time(), "%m%d_%H%M%S"), "_", substr(uuid::UUIDgenerate(), 1, 8))
# 输出示例：0708_021530_1a2b3c4d

# dir = fs::path_dir(inputfile)
dir = outputpath
dir.create(outputpath)


##### --- readtable ---
raw_cf_table <- read.table(cffile, header = TRUE, sep = "\t", row.names = 1)
### sort mut and cell in cf_table
cf_table <- sort_mutation_matrix(raw_cf_table)
all_mutid <- colnames(cf_table)
### sort input_table by order_cf_table
raw_input_table <- read.table(inputfile, header = TRUE, sep = "\t", row.names = 1)
input_table <- raw_input_table[match(rownames(cf_table), rownames(raw_input_table)), match(colnames(cf_table), colnames(raw_input_table))]

## Remove all 0 spots
outputree_dat_zero_barcode <- find_zero_barcode(cf_table)
# process cfflie
cffile_no0 <- remove_zero_barcode()
cf_tableno0 <- read.table(cffile_no0, header = TRUE, sep = "\t")
dim(cf_tableno0)
print(str_c("The spots drawn in circos are a total of ", as.character(dim(cf_tableno0)[1]), "."))


##### --- readtree ---
tree_dat <- cf2treedata(cffile_no0)
tree <- dat2tree(tree_dat)


##### --- spatial data ---
## import cell cluster and cell annotation information
mfile <- annotation_file
df_all_info <- read.table(mfile, header=TRUE, sep="\t"); dim(df_all_info)
# [1] 7038    8
# colnames(df_all_info)
# [1] "barcode"      "sample"       "cluster_info" "cell_type"    "tumor_score" 

### --- spatial data: cell cluster --- 
anno_color <- c()
if("cluster_info" %in% colnames(df_all_info)){
    # 构建 df_cluster
    df_cluster <- df_all_info[, c("barcode", "cluster_info")]
    df_cluster$cluster_info <- trimws(as.character(df_cluster$cluster_info))
    df_cluster <- df_cluster[!df_cluster$barcode %in% outputree_dat_zero_barcode, ]
    colnames(df_cluster) <- c("label", "cluster")
    # 设置 cluster 为 factor，排序也可以自定义
    cluster_levels <- sort(unique(df_cluster$cluster))  # 或自定义顺序
    df_cluster$cluster <- factor(df_cluster$cluster, levels = cluster_levels)
    # 设置颜色，确保与 cluster_levels 顺序对应
    set.seed(42)  # 保证颜色可复现（推荐加上）
    df_cluster_color <- rgb(runif(length(cluster_levels)), runif(length(cluster_levels)), runif(length(cluster_levels)))
    names(df_cluster_color) <- cluster_levels
    # 合并颜色
    anno_color <- c(anno_color, df_cluster_color)
}
### --- spatial data: cell type --- 
if("cell_type" %in% colnames(df_all_info)){
    # 设置 celltype 为 factor，并设定 levels 顺序
    df_celltype <- df_all_info[,c("barcode", "cell_type")]
    df_celltype$cell_type <- trimws(as.character(df_celltype$cell_type))
    df_celltype <- df_celltype[!df_celltype$barcode %in% outputree_dat_zero_barcode,]
    colnames(df_celltype) <- c("label","celltype")
    # 设置 celltype 为 factor，并设定 levels 顺序
    celltype_levels <- sort(unique(df_celltype$celltype))  # 或使用自定义顺序
    df_celltype$celltype <- factor(df_celltype$celltype, levels = celltype_levels)
    # 生成颜色，确保和 levels 顺序一一对应
    df_celltype_color <- paletteer_d("ggthemes::stata_s1rcolor")[seq_along(celltype_levels)]
    names(df_celltype_color) <- celltype_levels
    # 合并入 anno_color（如需要）
    # demo: for skin
    # df_celltype_color <- c("#F0E442", "#D55E00", "#56B4E9", "#009E73", "#5A14A5", "#C8C8C8", "#323232")
    # names(df_celltype_color) <- c('B_cell','Keratinocytes','Smooth_muscle_cells', 'T_cells', 'Fibroblasts','NK_cell','Epithelial_cells')
    anno_color <- c(anno_color, df_celltype_color)
}
### --- spatial data: sample --- 
if("sample" %in% colnames(df_all_info)){
    df_sample <- df_all_info[,c("barcode", "sample")]
    df_sample$sample <- trimws(as.character(df_sample$sample))
    df_sample <- df_sample[!df_sample$barcode %in% outputree_dat_zero_barcode,]
    colnames(df_sample) <- c("label","sample")
    # colors
    # list_sample <- df_sample$sample |> unique()
    # df_sample_color <- paletteer_d("RColorBrewer::Set3")[1:length(list_sample)]
    # names(df_sample_color) <- list_sample
    # demo: for skin
    df_sample_color <- c("#ebb16e", "#b2aad3")
    names(df_sample_color) <- c('P6_ST_vis_rep2','P6_ST_vis_rep1')
    anno_color <- c(anno_color, df_sample_color)
}
### --- spatial data: tumor score --- 
if("tumor_score" %in% colnames(df_all_info)){
    df_tumor <- df_all_info[,c("barcode", "tumor_score")]
    df_tumor$tumor_score <- as.numeric(as.character(df_tumor$tumor_score))
    min_value <- min(df_tumor$tumor_score[df_tumor$tumor_score != 0])
    max_value <- max(df_tumor$tumor_score)
    colnames(df_tumor) <- c("label","tumor_score")
    # colors
    df_tumor_color <- c("#00a3c4", "#ffebf6", "#ff64be")
    breaks <- c(0, min_value, max_value)
    anno_color <- c(anno_color, df_tumor_color)
}
### --- spatial data: tumor score --- 
if("B_cell_prop" %in% colnames(df_all_info)){
    df_Bcell <- df_all_info[,c("barcode", "B_cell_prop")]
    df_Bcell$B_cell_prop <- as.numeric(as.character(df_Bcell$B_cell_prop))
    min_value <- min(df_Bcell$B_cell_prop[df_Bcell$B_cell_prop != 0])
    max_value <- max(df_Bcell$B_cell_prop)
    median_value <- median(df_Bcell$B_cell_prop)
    colnames(df_Bcell) <- c("label","B_cell_prop")
    # colors
    df_Bcell_color <- c("#F4D166", "#DBE8B4", "#24693D")
    breaks <- c(min_value, median_value, max_value)
    anno_color <- c(anno_color, df_Bcell_color)
}



##### --- cells order ---
cf_tableno0_for_sorting_mat <- cf_tableno0[,2:dim(cf_tableno0)[2]]
rownames(cf_tableno0_for_sorting_mat) <- cf_tableno0[,1]

mat <- as.matrix(cf_tableno0_for_sorting_mat)

# 递归划分细胞克隆
cell_clone_structure <- recursive_clone_partition(rownames(mat), by_row = TRUE)
sorted_cells <- extract_mutation_order(cell_clone_structure)  # 细胞排序

# 递归划分 mutation 克隆
mut_clone_structure <- recursive_clone_partition(1:ncol(mat), by_row = FALSE)
sorted_mut_indices <- extract_mutation_order(mut_clone_structure)
sorted_mutations <- colnames(mat)[sorted_mut_indices]

# 用排序后的行列索引重新排列矩阵
sorted_mat <- mat[sorted_cells, sorted_mutations, drop = FALSE]

# 查看及保存
# print(head(sorted_mat[, 1:10]))
write.table(sorted_mat, file = str_c(outputpath, "/sorted_cf_matrix.txt"), sep = "\t", quote = FALSE)


# 生成要用于后续的数据框，形似 cf_tableno0
cf_tableno0_sorted <- data.frame(cellIDxmutID = rownames(sorted_mat), sorted_mat, stringsAsFactors = FALSE)
# 确保除第一列外，其他列全是 numeric
cf_tableno0_sorted[, -1] <- lapply(cf_tableno0_sorted[, -1], as.numeric)
rownames(cf_tableno0_sorted) <- rownames(cf_tableno0)


##### --- subclone order ---
subclones <- subclone_finder(cf_tableno0_sorted)
subclones$clone <- unlist(subclones$clone)
## sum cell row
cell_mutant_size <- data.frame(
    cells = cf_tableno0_sorted[,1],
    cell_mutant_size = rowSums(cf_tableno0_sorted[,-1])
)

## flat for each cells information
subclones_transformed <- separate_rows(subclones, cells, sep = ",")

## add cluster information
subclones_transformed <- subclones_transformed %>%
    left_join(
        if (exists("df_cluster")) df_cluster else df_celltype,
        by = c("cells" = "label")
  )

## add clone information to tree leaves (cells) df
leaf_df <- tree_dat[tree_dat$label %in% tree$tip.label,]
leaf_df <- leaf_df %>% left_join(subclones_transformed, by = c("label" = "cells"))

leaf_df$clone[is.na(leaf_df$clone)] <- "root_clone"

## target mutation as a clone
if(target_mut=="no"){
    print("There is no need to hightlight a target_mut.")
} else {
    target_clone <- cf_tableno0_sorted[cf_tableno0_sorted[[target_mut]] == 1, "cellIDxmutID"]
    leaf_df <- leaf_df %>% 
        mutate(target_clone = ifelse(label %in% target_clone, "clone_under_target_mut", "others"))
}


## order cells by clone
# 1. 设置行名
rownames(cf_tableno0_sorted) <- cf_tableno0_sorted$cellIDxmutID
# 2. 去掉第一列，只保留突变矩阵
mutation_mat <- cf_tableno0_sorted[, -1, drop = FALSE]
# 3. 确认 leaf_df$label 在 mutation 矩阵中是否都存在
missing_cells <- setdiff(leaf_df$label, rownames(mutation_mat))
if (length(missing_cells) > 0) {
  warning("这些细胞条码在突变矩阵中找不到，会导致计算出错：", paste(missing_cells, collapse = ", "))
}
# 4. 计算 mutation_count 并添加到 leaf_df
leaf_df$mutation_count <- rowSums(mutation_mat[leaf_df$label, , drop = FALSE])
# 5. 按多个条件排序 leaf_df
leaf_df <- leaf_df %>%
  arrange(subclone_size, clone_order) %>%
  group_by(clone_order) %>%
  arrange(desc(mutation_count), .by_group = TRUE) %>%
  ungroup()


## a preview of plot
clone_order_tree <- ape::rotateConstr(tree, leaf_df$label)

## sum mutation col
mutant_cell_size <- data.frame(
    mutant = colnames(cf_tableno0)[-1],
    mutant_cell_size = colSums(cf_tableno0[,-1])
)

merged_df <- leaf_df %>%
    left_join(tree_dat, by = c("parent" = "node"))

merged_df2 <- separate_rows(merged_df, label.y, sep = "\\|")

unique(merged_df2$label.y)

leaf_df_cell_size <- merged_df2 %>%
    left_join(mutant_cell_size, by = c("label.y" = "mutant"))

unique(leaf_df_cell_size$label.y)

## flat for each mutations information
subclones_transformed2 <- separate_rows(subclones, subclone, sep = ",")

# ## add cell size info
temp_cell_size <- unique(na.omit(leaf_df_cell_size[,c("subclone_num","mutant_cell_size","label.y")]))
subclones_transformed3 <- subclones_transformed2 %>%
    left_join(temp_cell_size, by = c("subclone_num" = "subclone_num"))

## more cell, more head
subclones_transformed4 <- subclones_transformed3 %>%
    arrange(desc(mutant_cell_size),clone_order)

# ## plot tree with heatmap of mutations
# ## reorder heatmap by clone
reordered_heatmap = cf_tableno0[match(leaf_df$label,cf_tableno0$cellIDxmutID),
                                colnames(cf_table)]

## flip ----
input_table_long <- input_table %>% as_tibble(rownames = "barcode") %>%
    tidyr::pivot_longer(-barcode,names_to = "mutation",values_to = "input")

cf_table_long <- cf_table %>% as_tibble(rownames = "barcode") %>%
    tidyr::pivot_longer(-barcode,names_to = "mutation",values_to = "cf")

## merge tibble of input and cf
all_long <- merge(input_table_long, cf_table_long)

## flipped type classification
all_long$flipped_type <- paste(all_long$input,all_long$cf,sep = ">")

## long to wide
all_wide <- all_long %>%
    select(-input, -cf) %>%
    tidyr::pivot_wider(values_from = flipped_type,
                       names_from = mutation,
                       values_fill = "not_flipped")%>%
    column_to_rownames("barcode")
# first column to row names
all_wide <- as.data.frame(all_wide)
all_mat  <- as.matrix(all_wide)
unique_types <- unique(all_long$flipped_type)
# # remove no fipped
# no_flipped_index <- strsplit(unique_types,">") %>%
#   lapply(function(x)var(x)==0) %>%
#   unlist()

# color pal ------------------------------------------------------------
# define color of ternary matrix
color_define = structure(c("#2c81be","#e17259","#d9d9d9"),
                         names = c("absence", "presence", "missing"))

color_number = structure(color_define,
                         names = c("0", "1", "3"))

# define before and after color of flipped type
# if start with 0, the before color is absence
# if start with 1, the before color is presence
# if start with 3, the before color is missing
# if end with 0, the after color is absence
# if end with 1, the after color is presence
# if end with 3, the after color is missing
unique_types <- unique(all_long$flipped_type)
# remove no fipped
no_flipped_index <- strsplit(unique_types,">") %>%
    lapply(function(x)var(x)==0) %>%
    unlist()
colors_before = dplyr::case_when(
    startsWith(unique_types, "0") ~ color_define["absence"],
    startsWith(unique_types, "1") ~ color_define["presence"],
    startsWith(unique_types, "3") ~ color_define["missing"]
)
names(colors_before) = unique_types

colors_after = dplyr::case_when(
    endsWith(unique_types, "0") ~ color_define["absence"],
    endsWith(unique_types, "1") ~ color_define["presence"],
    endsWith(unique_types, "3") ~ color_define["missing"]
)
names(colors_after) = unique_types

sub <- all_mat[!rownames(all_mat) %in% outputree_dat_zero_barcode,]
colnames_order = colnames(reordered_heatmap)[colnames(reordered_heatmap) %in% colnames(sub)]
sub <- sub[clone_order_tree$tip.label,colnames_order]

sub <- as.data.frame(sub)
sub$cellIDxmutID <- rownames(sub)

sub_long <- sub %>%
    pivot_longer(-cellIDxmutID, names_to = "mutants", values_to = "mut_count")

# pheatmap(reordered_heatmap,
#          cluster_rows = FALSE, cluster_cols = FALSE,
#          show_rownames = FALSE,
#          show_colnames = TRUE,
#          color = colorRampPalette(c("white", "black"))(100))

reordered_heatmap_long <- cf_tableno0 %>%
    pivot_longer(-cellIDxmutID, names_to = "mutants", values_to = "mut_count")

# clone color
## clone color
# clone_name <- unique(subclones$clone)
# clone_color <- c("#9d8982","#acb5b5","#A67F68","#c7a08f","#e0c288","#cdc4c3")
clone_order <- unique(subclones[order(subclones$clone_order,decreasing = TRUE),]$clone)
# root color
roots_index <- grep("root",clone_order)
colfunc <- colorRampPalette(c("#7b767a", "#a49b95", "#dad7d8"))
roots_color <- colfunc(length(roots_index))

# clone color
clone_index <- !grepl("root",clone_order)
clone_num <- sum(clone_index)
# if (clone_num <= 8) {
#   clone_color_o <- rev(RColorBrewer::brewer.pal(clone_num,"Pastel2"))
# } else {
#   clone_color_o <- rev(c(RColorBrewer::brewer.pal(clone_num-8,"Paired")[clone_num-8],
#                        RColorBrewer::brewer.pal(8,"Pastel2")))
# }
clone_color_o <- Polychrome::palette36.colors(clone_num)
clone_color <- vector(length = length(clone_order))
clone_color[roots_index] <- roots_color
clone_color[clone_index] <- clone_color_o
names(clone_color) <- clone_order

## subclone color == clone color
clone_subclone_names <- unique(subclones_transformed4[,c("clone","subclone_num","label.y","clone_order")])
clone_subclone_names <- clone_subclone_names %>% left_join(subclones[,c("subclone_num","subclone")],
                                                           by = c("subclone_num" = "subclone_num"))
clone_subclone_names <- clone_subclone_names[order(clone_subclone_names$clone_order),]
clone_subclone_names$color <- clone_color[clone_subclone_names$clone]

subclone_color <- clone_color[clone_subclone_names$clone]
names(subclone_color) <- clone_subclone_names$subclone

## subclone_mut_color
subclone_mut_color <- clone_color[clone_subclone_names$clone]
names(subclone_mut_color) <- clone_subclone_names$label.y

# write.table(clone_subclone_names, str_c(outputpath, "/clone_and_subclone_table.txt"), sep="\t", col.names=TRUE)

## mutatation color
# pivot_muts = unique(subclones_transformed4$label.y)
if(selected_mutlist=="all" | is.na(selected_mutlist)){
    pivot_muts = colnames(cf_table)
} else {
    pivot_muts = unlist(strsplit(selected_mutlist, ",\\s*"))
}
# color_muts = RColorBrewer::brewer.pal(length(pivot_muts),"Paired")
color_muts = pals::viridis(length(pivot_muts))
names(color_muts) = pivot_muts

sub_long$mutants_symbol <- paste("M",as.numeric(factor(sub_long$mutants,levels = pivot_muts)),sep = "")
for_sort_symbol <- unique(sub_long[,c("mutants","mutants_symbol")])
pivot_muts_symbol = for_sort_symbol$mutants_symbol
pivot_muts_sub = for_sort_symbol$mutants


##########################################################################
##########################################################################
##########################################################################

# #### plotting parameters ---------------------------------------------------------------

# # Perform adaptive and automated settings for the drawing parameters

# ### Get basic numbers
# num_cells <- length(clone_order_tree$tip.label)
# print(str_c("===== The phylogenetic tree will include ", num_cells, " cells ====="))

# num_mutations <- length(pivot_muts)  # 或 length(pivot_muts_symbol)
# print(str_c("===== The phylogenetic tree will include ", num_mutations, " mutations. ====="))

# target_cols <- c("cell_type", "cluster_info", "sample", "tumor_score", "B_cell_prop")
# num_anno_circos = sum(sapply(target_cols, function(x) sum(colnames(df_all_info) == x)))
# print(str_c("===== The number of annotation circles is ", num_anno_circos, " ====="))

# ### Predict plotting parameters automatically
# params <- predict_plot_params(num_mutations, num_cells, num_anno_circos)

# # Assign predicted values to plotting variables
# tip_label_offset <- params$tip_label_offset
# tip_label_size <- params$tip_label_size
# tip_point_size <- params$tip_point_size
# heatmap_circos_width <- params$heatmap_circos_width
# heatmap_circos_offset <- params$heatmap_circos_offset
# flipping_point_size <- params$flipping_point_size
# plot_height <- params$plot_height
# plot_width <- params$plot_width


# # ##### Setting plotting parameters manully
# tip_label_offset=7.5
# tip_label_size=2.5
# tip_point_size=0.5
# heatmap_circos_width=0.5
# heatmap_circos_offset=0.04
# flipping_point_size=1.2
# plot_height=12
# plot_width=18

# ### Display parameters (optional)
# print("===== Adaptive plotting parameters =====")
# print(params)




##########################################################################
##########################################################################
##########################################################################

##### plot ---------------------------------------------------------------
library(ggtreeExtra)


title_theme = element_text(
  size = 20,
  # face = "italic",
  colour = "#3c3c3c",
  angle = 0
)
clone_labels <- paste("(",clone_subclone_names$clone,") ", clone_subclone_names$subclone, sep = "")
new_clone_labels <- lapply(clone_labels, insert_newline, n = 180)
new_clone_labels <- unlist(new_clone_labels)

saveRDS(clone_order_tree, str_c(outputpath, "/CNVtree_data_clone_order_tree.rds"))


##### plot tree with clone label
num_cells <- length(clone_order_tree$tip.label)

# 复制一个新的 clone_order_tree，绘图的 mutation node 显示 id 的时候就最多显示两个
clone_order_tree_trim <- clone_order_tree
max_mutnum_on_one_node <- max(sapply(clone_order_tree$node.label, function(x) {
  if (!is.na(x)) {
    # 拆分字符串，并计算子元素数量
    return(length(strsplit(x, "\\|")[[1]]))
  } else {
    return(0)  # 如果标签是 NA，则返回 0
  }
}))
max_mutnum_on_one_node <- pmin(max_mutnum_on_one_node, 5)
clone_order_tree_trim$node.label <- sapply(
  clone_order_tree_trim$node.label, 
  function(x) trim_label(x, target_mut = target_mut, default_n = max_mutnum_on_one_node)
)


if(num_cells >= 72){

    if(target_mut=="no"){

        # 所有 label 都是白色
        No_target_clone_mut_color <- rep("white", length(clone_order_tree_trim$node.label))
        names(No_target_clone_mut_color) <- sapply(clone_order_tree_trim$node.label, function(x) gsub("\\|", "\n", x))

        p1_tree <- ggtree(clone_order_tree_trim, size = 0.1, 
                          branch.length = "none", 
                          ladderize = FALSE, 
                          layout = "radial", 
                          color = "black") +  # 将树的分支颜色设置为灰色
          geom_tippoint(color = "black", size = tip_point_size) +  # 将叶点颜色设置为灰色
          xlim(-1, NA) +
          geom_label2(aes(subset = !isTip,
                          label = gsub("\\|", "\\\n", label), 
                          fill = label),
                      color = "#3c3c3c", 
                      show.legend = FALSE, 
                      alpha = .6) +
          geom_tiplab(align = FALSE, 
                      offset = tip_label_offset,
                      linesize = 0, size = tip_label_size, 
                      show.legend = FALSE) +
          scale_fill_manual(values = No_target_clone_mut_color)  # 如果需要，可以修改 `No_target_clone_mut_color` 中的颜色

    } else {

        # print(str_c(target_mut, ": start! "))
        # # 创建颜色映射：label 包含 target_mut 的为红色，其他为白色
        # label_vec <- sapply(clone_order_tree_trim$node.label, function(x) gsub("\\|", "\n", x))
        # target_clone_mut_color <- ifelse(grepl(target_mut, clone_order_tree_trim$node.label), "red", "white")
        # names(target_clone_mut_color) <- label_vec

        # p1_tree <- ggtree(clone_order_tree_trim, size = 0.1, 
        #                   branch.length = "none", 
        #                   ladderize = FALSE,
        #                   layout = "radial", 
        #                   aes(color = target_clone)) %<+% leaf_df[,c("label", "target_clone")] +
        #     geom_tippoint(aes(color = target_clone), size = tip_point_size) +
        #     xlim(-1, NA) +
        #     geom_label2(aes(subset = !isTip,
        #                     label = gsub("\\|", "\\\n", label), 
        #                     fill = label),
        #                 color = "#3c3c3c", 
        #                 show.legend = FALSE, 
        #                 alpha = .6) +
        #     geom_tiplab(align = FALSE, 
        #                 offset = tip_label_offset,
        #                 linesize = 0, size = tip_label_size, 
        #                 show.legend = FALSE) +
        #     scale_fill_manual(values = target_clone_mut_color, guide = "none")

        print(str_c(target_mut, ": start! "))
        target_clone_mut_color <- c("red", rep("white", length(pivot_muts)-1))
        names(target_clone_mut_color) <- c(target_mut, setdiff(pivot_muts, target_mut))

        p1_tree <- ggtree(clone_order_tree, size = 0.1, 
                          branch.length = "none", 
                          ladderize = FALSE,
                          layout = "radial", 
                          aes(color = target_clone)) %<+% leaf_df[,c("label", "target_clone")] +
            geom_tippoint(aes(color = target_clone), size = tip_point_size) +
            xlim(-1, NA) +
            scale_color_manual(values = c("clone_under_target_mut" = "red", 
                                        "others" = "black"),
                               guide = "none") +
            geom_label2(aes(subset = !isTip,
                            label = gsub("\\|", "\\\n", label), 
                            fill = label),
                        color = "#3c3c3c", 
                        show.legend = FALSE, 
                        alpha = .6) +
            geom_tiplab(align = FALSE, 
                        offset = tip_label_offset,
                        linesize = 0, size = tip_label_size, 
                        show.legend = FALSE) +
            scale_fill_manual(values = target_clone_mut_color, guide = "none")

  }

} else {
    ### plot sector tree in case a unit would be too wide
    # Set the angular width of each cell in degrees
    angle_per_cell <- 5
    # Calculate the total Angle required
    total_angle <- angle_per_cell * num_cells
    # Make sure the total Angle does not exceed 360 degrees
    if (total_angle > 360) {
        total_angle <- 360
    }
    # Calculate open.angle (the Angle of the blank part of the circle)
    open_angle <- 360 - total_angle

    if(target_mut=="no"){

        # 所有 label 都是白色
        No_target_clone_mut_color <- rep("white", length(clone_order_tree_trim$node.label))
        names(No_target_clone_mut_color) <- sapply(clone_order_tree_trim$node.label, function(x) gsub("\\|", "\n", x))

        p1_tree <- ggtree(clone_order_tree_trim, size = 0.3, 
                          branch.length = "none", 
                          ladderize = FALSE,
                          layout = "fan",
                          color = "black",
                          open.angle = open_angle) +
            geom_tippoint(color = "black", size = tip_point_size) +
            xlim(-1, NA) +
            geom_label2(aes(subset = !isTip,
                            label = gsub("\\|", "\\\n", label), 
                            fill = label),
                        color = "#3c3c3c", 
                        show.legend = FALSE, 
                        alpha = .6) +
            geom_tiplab(align = FALSE, 
                        offset = tip_label_offset,
                        linesize = 0, size = tip_label_size, 
                        show.legend = FALSE) +
            scale_fill_manual(values = No_target_clone_mut_color)

    } else {

        # print(str_c(target_mut, ": start! "))
        # # 创建颜色映射：label 包含 target_mut 的为红色，其他为白色
        # label_vec <- sapply(clone_order_tree_trim$node.label, function(x) gsub("\\|", "\n", x))
        # target_clone_mut_color <- ifelse(grepl(target_mut, clone_order_tree_trim$node.label), "red", "white")
        # names(target_clone_mut_color) <- label_vec

        # p1_tree <- ggtree(clone_order_tree_trim, size = 0.3, 
        #                   branch.length = "none", 
        #                   ladderize = FALSE,
        #                   layout = "fan",
        #                   color = "black",
        #                   open.angle = open_angle) +
        #     geom_tippoint(color = "black", size = tip_point_size) +
        #     xlim(-1, NA) +
        #     geom_label2(aes(subset = !isTip,
        #                     label = gsub("\\|", "\\\n", label), 
        #                     fill = label),
        #                 color = "#3c3c3c", 
        #                 show.legend = FALSE, 
        #                 alpha = .6) +
        #     geom_tiplab(align = FALSE, 
        #                 offset = tip_label_offset,
        #                 linesize = 0, size = tip_label_size, 
        #                 show.legend = FALSE) +
        #     scale_fill_manual(values = target_clone_mut_color, guide = "none")

        print(str_c(target_mut, ": start! "))
        target_clone_mut_color <- c("red", rep("white", length(pivot_muts)-1))
        names(target_clone_mut_color) <- c(target_mut, setdiff(pivot_muts, target_mut))

        p1_tree <- ggtree(clone_order_tree, size = 0.1, 
                          branch.length = "none", 
                          ladderize = FALSE,
                          layout = "radial", 
                          aes(color = target_clone)) %<+% leaf_df[,c("label", "target_clone")] +
            geom_tippoint(aes(color = target_clone), size = tip_point_size) +
            xlim(-1, NA) +
            scale_color_manual(values = c("clone_under_target_mut" = "red", 
                                        "others" = "black"),
                               guide = "none") +
            geom_label2(aes(subset = !isTip,
                            label = gsub("\\|", "\\\n", label), 
                            fill = label),
                        color = "#3c3c3c", 
                        show.legend = FALSE, 
                        alpha = .6) +
            geom_tiplab(align = FALSE, 
                        offset = tip_label_offset,
                        linesize = 0, size = tip_label_size, 
                        show.legend = FALSE) +
            scale_fill_manual(values = target_clone_mut_color, guide = "none")


    }

}

# svg_filename_p1_tree <- str_c(paste0(outputpath,"/p1_tree.circle_tree_output_as_point.", pdf_lastfix),".svg")
# ggsave(file = svg_filename_p1_tree, width = plot_width*2.5, height = plot_height*2, plot = p1_tree)
# ggsave(file = str_c(paste0(outputpath,"/p1_tree.circle_tree_output_as_point.", pdf_lastfix),".pdf"), width = plot_width*2.5, height = plot_height*2, plot = p1_tree)


# #### plot cell annotation
p2_anno <- p1_tree
if ("cell_type" %in% colnames(df_all_info)) {
    p2_anno <- p2_anno + new_scale_fill() +
        geom_fruit(data = df_celltype,
                   geom = geom_tile, pwidth=0.2, width = 0.2,
                   mapping = aes(y = label, x = "celltype", fill = celltype),
                   offset = 0.03, color = NA) +
        scale_fill_manual(values = df_celltype_color,
                          guide = guide_legend(title = "Cell Type", 
                                               title.theme = title_theme,
                                               ncol = 3, title.position = "top",
                                               override.aes = list(shape = 18, size = 5),
                                               order = 1)) +
        theme(legend.direction = "horizontal",
              legend.box = "vertical",
              legend.box.just = "left",
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14))
} 
if ("cluster_info" %in% colnames(df_all_info)) {
    p2_anno <- p2_anno + new_scale_fill() +
        geom_fruit(data = df_cluster,
                   geom = geom_tile, pwidth=0.2, width = 0.2,
                   mapping = aes(y = label, x = "cluster", fill = cluster),
                   offset = 0.03, color = NA) +
        scale_fill_manual(values = df_cluster_color,
                          guide = guide_legend(title = "Cluster", 
                                               title.theme = title_theme,
                                               ncol = 6, title.position = "top",
                                               override.aes = list(shape = 18, size = 5),
                                               order = 2)) +
        theme(legend.direction = "horizontal",
              legend.box = "vertical",
              legend.box.just = "left",
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14))
}
if ("sample" %in% colnames(df_all_info)) {
    p2_anno <- p2_anno + new_scale_fill() +
        geom_fruit(data = df_sample,
                   geom = geom_tile, pwidth=0.2, width = 0.2,
                   mapping = aes(y = label, x = "sample", fill = sample),
                   offset = 0.03, color = NA) +
        scale_fill_manual(values = df_sample_color,
                          guide = guide_legend(title = "Sample", 
                                               title.theme = title_theme,
                                               ncol = 4, title.position = "top",
                                               override.aes = list(shape = 18, size = 5),
                                               order = 3)) +
        theme(legend.direction = "horizontal",
              legend.box = "vertical",
              legend.box.just = "left",
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14))
}
if ("tumor_score" %in% colnames(df_all_info)) {
    df_tumor$tumor_score <- as.numeric(df_tumor$tumor_score)
    p2_anno <- p2_anno + new_scale_fill() +
        geom_fruit(data = df_tumor,
                   geom = geom_tile, pwidth = 0.2, width = 0.2,
                   mapping = aes(y = label, x = "tumor_score", fill = tumor_score),
                   offset = 0.03, color = NA) +
        scale_fill_gradientn(colors = df_tumor_color,
                             values = scales::rescale(breaks),
                             guide = guide_colorbar(title = "Tumor Score", 
                                                    title.theme = title_theme,
                                                    title.position = "top",
                                                    barwidth = 10,  # 设置颜色条的宽度
                                                    barheight = 0.5,  # 设置颜色条的高度
                                                    order = 4)) +
        theme(legend.direction = "horizontal",
              legend.box = "vertical",
              legend.box.just = "left",
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14))
}
if("B_cell_prop" %in% colnames(df_all_info)){
    df_Bcell$B_cell_prop <- as.numeric(df_Bcell$B_cell_prop)
    p2_anno <- p2_anno + new_scale_fill() +
        geom_fruit(data = df_Bcell,
                   geom = geom_tile, pwidth = 0.2, width = 0.2,
                   mapping = aes(y = label, x = "B_cell_prop", fill = B_cell_prop),
                   offset = 0.03, color = NA) +
        scale_fill_gradientn(colors = df_Bcell_color,
                             values = scales::rescale(breaks),
                             guide = guide_colorbar(title = "B Cell Proportion", 
                                                    title.theme = title_theme,
                                                    title.position = "top",
                                                    barwidth = 10,  # 设置颜色条的宽度
                                                    barheight = 0.5,  # 设置颜色条的高度
                                                    order = 5)) +
        theme(legend.direction = "horizontal",
              legend.box = "vertical",
              legend.box.just = "left",
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14))
}

# svg_filename_p2_anno <- str_c(paste0(outputpath,"/p2_anno.circle_tree_output_as_point.", pdf_lastfix),".svg")
# ggsave(file = svg_filename_p2_anno, width = plot_width, height = plot_height, plot = p2_anno)
# ggsave(file = str_c(paste0(outputpath,"/p2_anno.circle_tree_output_as_point.", pdf_lastfix),".pdf"), width = plot_width, height = plot_height, plot = p2_anno)


##### plot p3_heatmap including heatmap
### the guide parameters of p3_heatmap
guide_filp_point <- guide_legend(title = "Genotyper (raw > inferred)", 
                                 title.position = "top",
                                 override.aes = list(size = 10), 
                                 ncol = 1,
                                 title.theme = title_theme)

## manual fp flipping
if(manual_fp_file=="no"){
    new_colors_before <- colors_before
    new_colors_after <- colors_after
} else {

    df_manual_fp <- read.table(manual_fp_file, sep='\t', header=TRUE); dim(df_manual_fp)
    # highlight unit
    sub_long_bk <- sub_long
    sub_long$mut_count <- ifelse(paste(sub_long$cellIDxmutID, sub_long$mutants) %in% 
                                 paste(df_manual_fp$cellID, df_manual_fp$mutID),
                                 str_c("manual1", sub(".*(>.)", "\\1", sub_long$mut_count)), 
                                 sub_long$mut_count)
    # sub_long[sub_long$cellIDxmutID=="TATTATGTTTGCCTGC-1",]
    # sub_long[sub_long$cellIDxmutID=="TATTATGTTTGCCTGC-1" & sub_long$mutants=="chr11_65426524_T_C", "mut_count"]
    new_colors_before <- c(colors_before, "manual1>0" = "#FFFF00", "manual1>1" = "#FFFF00") # 亮黄色
    new_colors_after <- c(colors_after, "manual1>0" = "#2c81be", "manual1>1" = "#e17259") # 亮黄色
    # remove no fipped
    unique_mut_counts <- names(new_colors_before)
    no_flipped_index <- sapply(unique_mut_counts, function(x) {
    if ("manual" %in% x) {
        FALSE
    } else {
        x_split <- strsplit(x, ">")[[1]]
        length(x_split) == 2 && x_split[1] == x_split[2]
    }
    })

}

labels_dict <- setNames(
  format_flipping_label(names(new_colors_before)),
  names(new_colors_before)
)

if(target_mut=="no"){

    p3_heatmap <- p2_anno + new_scale_fill()+
        geom_fruit_list(
            geom_fruit(data = sub_long,
                       geom = geom_tile, 
                       pwidth = heatmap_circos_width,
                       mapping = aes(y = cellIDxmutID, 
                                     x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                     fill = as.character(mut_count)),
                       offset = heatmap_circos_offset,
                       color = "white",
                       size = 0.2,
                       axis.params = list(axis = "x",
                                          line.size = NA,
                                          line.color = NA,
                                          text.size = 0,
                                          vjust = 0.5)),
            scale_fill_manual(values = new_colors_before, 
                              breaks = names(new_colors_before)[!no_flipped_index],
                              labels = labels_dict[names(new_colors_before)[!no_flipped_index]],
                              guide = guide_filp_point),
            new_scale_fill(),
            geom_fruit(data = sub_long,
                       geom = geom_point, 
                       pwidth = heatmap_circos_width,
                       mapping = aes(y = cellIDxmutID, 
                                     x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                     fill = as.character(mut_count)),
                       offset = heatmap_circos_offset, 
                       size = flipping_point_size,
                       shape = 21, 
                       color="transparent"),
            scale_fill_manual(values = new_colors_after, 
                              breaks = names(new_colors_after)[!no_flipped_index],
                              labels = labels_dict[names(new_colors_after)[!no_flipped_index]],
                              guide = guide_filp_point)
        )  +
        theme(legend.text = element_text(size=20),
              axis.title.x = element_blank(),
              axis.text.x = element_blank())

} else {

    ## highlight target muation and its conflict mutations
    if (length(conflict_muts) == 0 || is.na(conflict_muts)) {
        plot_conflict_muts <- NA
    } else {
        plot_conflict_muts <- unlist(strsplit(conflict_muts, ",\\s*"))
    }

    # highlight
    highlight_mutid <- sub_long %>%
    filter(mutants == gsub(":", ".", target_mut)) %>%
    pull(mutants_symbol) %>%
    unique() %>%
    as.character()

    # circos heatmap
    p3_heatmap <- p2_anno + new_scale_fill() +
        geom_fruit_list(
            geom_fruit(data = sub_long,
                       geom = geom_tile, 
                       pwidth = heatmap_circos_width,
                       mapping = aes(y = cellIDxmutID, 
                                     x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                     fill = as.character(mut_count)),
                       offset = heatmap_circos_offset,
                       color = ifelse(factor(sub_long$mutants_symbol, levels = pivot_muts_symbol) == highlight_mutid, "black", "white"),
                       size = ifelse(factor(sub_long$mutants_symbol, levels = pivot_muts_symbol) == highlight_mutid, 0.2, 0.2),
                       axis.params = list(axis = "x",
                                          line.size = NA,
                                          line.color = NA,
                                          text.size = 0,
                                          vjust = 0.5)),
            scale_fill_manual(values = new_colors_before, 
                              breaks = names(new_colors_before)[!no_flipped_index],
                              labels = labels_dict[names(new_colors_before)[!no_flipped_index]],
                              guide = guide_filp_point),
            new_scale_fill(),
            geom_fruit(data = sub_long,
                       geom = geom_point, 
                       pwidth = heatmap_circos_width,
                       mapping = aes(y = cellIDxmutID, 
                                     x = factor(mutants_symbol, levels = pivot_muts_symbol),
                                     fill = as.character(mut_count)),
                       offset = heatmap_circos_offset, 
                       size = flipping_point_size,
                       shape = 21, 
                       color = "transparent"),
            scale_fill_manual(values = new_colors_after, 
                              breaks = names(new_colors_after)[!no_flipped_index],
                              labels = labels_dict[names(new_colors_after)[!no_flipped_index]],
                              guide = guide_filp_point)
            ) +
        theme(legend.text = element_text(size=20),
              axis.title.x = element_blank(),
              axis.text.x = element_blank())

}

# svg_filename_p3_heatmap <- str_c(paste0(outputpath,"/p3_heatmap.circle_tree_output_as_point.", pdf_lastfix),".svg")
# ggsave(file = svg_filename_p3_heatmap, width = plot_width*2.5, height = plot_height*2, plot = p3_heatmap)
# ggsave(file = str_c(paste0(outputpath,"/p3_heatmap.circle_tree_output_as_point.", pdf_lastfix),".pdf"), width = plot_width*2.5, height = plot_height*2, plot = p3_heatmap)


##### plot Add mutation or barcode hotspots to the outermost ring

if (adding_mut_file=="no" && barcode_name=="no") {

    print("There is no need to add a circle line and the corresponding barcode to the plot.")

} else if (barcode_name=="no" && adding_mut_file!="no") {
##### Adding a circle heatmap for a site to be verified

    print("There is no need to plot corresponding barcode.")

    mutinfo_cols <- c("barcode_name", "cluster","in_tissue","pos_x","pos_y","depth","vaf","mutation_prob","mosaic_likelihood", "mutant_allele_num","has_mutant_allele","high_prob","high_likelihood","mutated")
    df_adding_mut <- read.table(adding_mut_file, sep = '\t', header = FALSE, col.names = mutinfo_cols); dim(df_adding_mut)
    # [1] 115  14
    df_mutantallele <- df_adding_mut[,c("barcode_name", "has_mutant_allele")]; dim(df_mutantallele)
    # [1] 115   2

    ### Get new circle data formated dataframe
    # Get new mut information
    last_element <- tail(strsplit(adding_mut_file, "/")[[1]], 1)
    adding_mutid <- strsplit(last_element, "\\.")[[1]][2]
    df_adding_bin <- as.data.frame(df_mutantallele$has_mutant_allele)
    colnames(df_adding_bin) <- adding_mutid
    rownames(df_adding_bin) <- df_mutantallele$barcode_name

    # fixed spot id
    fixed_spotid <- rownames(input_table); length(fixed_spotid)
    # [1] 431
    # 提取 df_adding_bin 中在 fixed_spotid 中存在的行
    intersect_spots <- intersect(fixed_spotid, rownames(df_adding_bin))
    df_intersect <- as.data.frame(df_adding_bin[intersect_spots, ])
    rownames(df_intersect) <- intersect_spots
    colnames(df_intersect) <- colnames(df_adding_bin)
    # 找出 fixed_spotid 中但不在 df_adding_bin 行名中的 spot id
    filled_spots <- setdiff(fixed_spotid, rownames(df_adding_bin))
    # 创建一个新数据框，行名为 filled_spots，列名与 df_adding_bin 相同，内容填充为 3
    df_filled <- data.frame(matrix(3, nrow = length(filled_spots), ncol = ncol(df_adding_bin)))
    rownames(df_filled) <- filled_spots
    colnames(df_filled) <- colnames(df_adding_bin)

    # 将 df_existing 和 df_filled 合并，形成新的 df_adding_bin
    df_adding_bin_updated <- rbind(df_intersect, df_filled)

    ### Generate new circle data for plotting
    adding_long <- df_adding_bin_updated %>% as_tibble(rownames = "barcode") %>%
        tidyr::pivot_longer(-barcode,names_to = "mutation",values_to = "input")

    ### plotting
    # color
    mut_types <- as.character(unique(df_adding_bin_updated[,1]))
    adding_colors = dplyr::case_when(
    endsWith(mut_types, "0") ~ color_define["absence"],
    endsWith(mut_types, "1") ~ color_define["presence"],
    endsWith(mut_types, "3") ~ color_define["missing"]
    )
    names(adding_colors) = mut_types
    # circos
    p4_adding <- p3_heatmap + new_scale_fill() +
        geom_fruit_list(
            geom_fruit(data = adding_long,
                     geom = geom_tile, 
                     pwidth = 0.2,
                     width = 0.2,
                     mapping = aes(y = barcode, 
                                   x = factor(mutation),
                                   fill = as.character(input)),
                     offset = 0.015,
                     color = "white",
                     size = 0.2,
                     axis.params = list(axis = "x",
                                        line.size = NA,
                                        line.color = NA,
                                        text.size = 0,
                                        vjust = 0.5,
                                        hjust = 0.2)),
            scale_fill_manual(values = adding_colors, 
                              breaks = names(adding_colors),
                              guide = guide_legend(title = "Mutation types", 
                                                   title.position = "top",
                                                   override.aes = list(size = 10), 
                                                   title.theme = title_theme))
          )

    # svg_filename_p4_adding <- str_c(paste0(outputpath,"/p4_adding.circle_tree_output_as_point.", pdf_lastfix),".svg")
    # svg(file = svg_filename_p4_adding, width = plot_width*3, height = plot_height*2)
    # p4_adding
    # dev.off()
    # ggsave(file = svg_filename_p1_tree, width = plot_width*2.5, height = plot_height*2, plot = p1_tree)
  
} else if (adding_mut_file=="no" && barcode_name!="no") {
##### Adding a circle heatmap for a cre-loxP_barcode

    print("There is no need to plot adding circle line.")

    df_barcode <- read.table(barcode_file, sep = '\t', header = TRUE); dim(df_barcode)
    colnames(df_barcode) <- c("barcode", "1", "149", "1B9", "9", "1BG", "3", "G", "I", "12I", "A", "123", "5", "C69", "12345", "129", "12G", "189", "1DCFE89", "7", "E4C", "IHG", "123D9", "12789", "12C", "167", "16789", "1HGFEDC", "369", "56CBAD9", "789", "123H9", "569", "C", "185", "1DC", "1H9", "589", "9DCH5", "349", "36I", "127", "12C45", "3BG", "E", "1B789", "34789", "389", "9FE", "A49", "C2IHG", "1B7", "CBA", "14569", "3DIHGFE", "9BA", "345", "12369", "14567", "1234789", "1236789", "1HGF9", "GF9", "56789", "56IHG", "1234I", "36789", "9H5", "CBA89", "E67", "C4789", "5FG8I", "EDCB9", "1HGFCB9", "347", "1BCD58I", "C89", "12389", "1456789", "1F9", "5HGF9", "1B589", "1D789", "ED789", "1FED789", "169", "749", "CBG85", "GD9", "3DI", "GF389", "12349", "12ED9", "34569", "EDC", "1234567", "58I", "129FC", "A89", "1D5", "I8GFEDC")
    #[1] 306 115

    ## Function to remove leading 'X' if present
    #remove_leading_x <- function(x) {
    #  ifelse(substr(x, 1, 1) == "X", substr(x, 2, nchar(x)), x)
    #}

    ## Apply the function to column names
    #colnames(df_barcode) <- remove_leading_x(colnames(df_barcode))

    df_majorbarcode <- df_barcode[,c("barcode", barcode_name)]; dim(df_majorbarcode)
    colnames(df_majorbarcode) <- c("barcode", "has_barcode")
    # [1] 306   2

    ### Get new circle data formated dataframe
    # Get barcode information
    df_adding_bin <- as.data.frame(df_majorbarcode$has_barcode)
    colnames(df_adding_bin) <- barcode_name
    rownames(df_adding_bin) <- df_majorbarcode$barcode

    # fixed spot id
    fixed_spotid <- rownames(input_table); length(fixed_spotid)
    # [1] 431
    # 提取 df_adding_bin 中在 fixed_spotid 中存在的行
    intersect_spots <- intersect(fixed_spotid, rownames(df_adding_bin))
    df_intersect <- as.data.frame(df_adding_bin[intersect_spots, ])
    rownames(df_intersect) <- intersect_spots
    colnames(df_intersect) <- colnames(df_adding_bin)
    # 找出 fixed_spotid 中但不在 df_adding_bin 行名中的 spot id
    filled_spots <- setdiff(fixed_spotid, rownames(df_adding_bin))
    # 创建一个新数据框，行名为 filled_spots，列名与 df_adding_bin 相同，内容填充为 3
    df_filled <- data.frame(matrix(3, nrow = length(filled_spots), ncol = ncol(df_adding_bin)))
    rownames(df_filled) <- filled_spots
    colnames(df_filled) <- colnames(df_adding_bin)

    # 将 df_existing 和 df_filled 合并，形成新的 df_adding_bin
    df_adding_bin_updated <- rbind(df_intersect, df_filled)
  
    ### Generate new circle data for plotting
    adding_long <- df_adding_bin_updated %>% as_tibble(rownames = "barcode") %>%
        tidyr::pivot_longer(-barcode,names_to = "creloxPbarcode",values_to = "input")
  
    ### plotting
    # color
    barcode_types <- as.character(unique(df_adding_bin_updated[,1]))
    adding_colors = dplyr::case_when(
        endsWith(barcode_types, "0") ~ "#EEEBDC",
        endsWith(barcode_types, "1") ~ "#514549",
    )
    names(adding_colors) = barcode_types
    # circos
    p5_barcode <- p3_heatmap + new_scale_fill() +
        geom_fruit_list(
            geom_fruit(data = adding_long,
                     geom = geom_tile, 
                     pwidth = 0.2,
                     width = 0.2,
                     mapping = aes(y = barcode, 
                                   x = factor(creloxPbarcode),
                                   fill = as.character(input)),
                     offset = 0.02,
                     color = "white",
                     size = 0.2,
                     axis.params = list(axis = "x",
                                        line.size = NA,
                                        line.color = NA,
                                        text.size = 0,
                                        vjust = 0.5,
                                        hjust = 0.2)),
            scale_fill_manual(values = adding_colors, 
                              breaks = names(adding_colors),
                              guide = guide_legend(title = "cre-loxP barcode", 
                                                   title.position = "top",
                                                   override.aes = list(size = 10), 
                                                   title.theme = title_theme))
        )

    # svg_filename_p5_barcode <- str_c(paste0(outputpath,"/p5_barcode.circle_tree_output_as_point.", pdf_lastfix),".svg")
    # svg(file = svg_filename_p5_barcode, width = plot_width*3, height = plot_height*2)
    # p5_barcode
    # dev.off()
    # ggsave(file = svg_filename_p1_tree, width = plot_width*2.5, height = plot_height*2, plot = p1_tree)
  
} else if (adding_mut_file!="no" && barcode_name!="no") {

    ##### adding
    mutinfo_cols <- c("barcode_name", "cluster","in_tissue","pos_x","pos_y","depth","vaf","mutation_prob","mosaic_likelihood", "mutant_allele_num","has_mutant_allele","high_prob","high_likelihood","mutated")
    df_adding_mut <- read.table(adding_mut_file, sep = '\t', header = FALSE, col.names = mutinfo_cols); dim(df_adding_mut)
    # [1] 115  14
    df_mutantallele <- df_adding_mut[,c("barcode_name", "has_mutant_allele")]; dim(df_mutantallele)
    # [1] 115   2

    ### Get new circle data formated dataframe
    # Get new mut information
    last_element <- tail(strsplit(adding_mut_file, "/")[[1]], 1)
    adding_mutid <- strsplit(last_element, "\\.")[[1]][2]
    df_adding_bin <- as.data.frame(df_mutantallele$has_mutant_allele)
    colnames(df_adding_bin) <- adding_mutid
    rownames(df_adding_bin) <- df_mutantallele$barcode_name

    # fixed spot id
    fixed_spotid <- rownames(input_table); length(fixed_spotid)
    # [1] 431
    # 提取 df_adding_bin 中在 fixed_spotid 中存在的行
    intersect_spots <- intersect(fixed_spotid, rownames(df_adding_bin))
    df_intersect <- as.data.frame(df_adding_bin[intersect_spots, ])
    rownames(df_intersect) <- intersect_spots
    colnames(df_intersect) <- colnames(df_adding_bin)
    # 找出 fixed_spotid 中但不在 df_adding_bin 行名中的 spot id
    filled_spots <- setdiff(fixed_spotid, rownames(df_adding_bin))
    # 创建一个新数据框，行名为 filled_spots，列名与 df_adding_bin 相同，内容填充为 3
    df_filled <- data.frame(matrix(3, nrow = length(filled_spots), ncol = ncol(df_adding_bin)))
    rownames(df_filled) <- filled_spots
    colnames(df_filled) <- colnames(df_adding_bin)

    # 将 df_existing 和 df_filled 合并，形成新的 df_adding_bin
    df_adding_bin_updated <- rbind(df_intersect, df_filled)

    ### Generate new circle data for plotting
    adding_long1 <- df_adding_bin_updated %>% as_tibble(rownames = "barcode") %>%
        tidyr::pivot_longer(-barcode,names_to = "mutation",values_to = "input")

    ### plotting
    # color
    mut_types <- as.character(unique(df_adding_bin_updated[,1]))
    adding_colors1 = dplyr::case_when(
        endsWith(mut_types, "0") ~ color_define["absence"],
        endsWith(mut_types, "1") ~ color_define["presence"],
        endsWith(mut_types, "3") ~ color_define["missing"]
    )
    names(adding_colors1) = mut_types

    ##### barcode
    df_barcode <- read.table(barcode_file, sep = '\t', header = TRUE); dim(df_barcode)
    colnames(df_barcode) <- c("barcode", "1", "149", "1B9", "9", "1BG", "3", "G", "I", "12I", "A", "123", "5", "C69", "12345", "129", "12G", "189", "1DCFE89", "7", "E4C", "IHG", "123D9", "12789", "12C", "167", "16789", "1HGFEDC", "369", "56CBAD9", "789", "123H9", "569", "C", "185", "1DC", "1H9", "589", "9DCH5", "349", "36I", "127", "12C45", "3BG", "E", "1B789", "34789", "389", "9FE", "A49", "C2IHG", "1B7", "CBA", "14569", "3DIHGFE", "9BA", "345", "12369", "14567", "1234789", "1236789", "1HGF9", "GF9", "56789", "56IHG", "1234I", "36789", "9H5", "CBA89", "E67", "C4789", "5FG8I", "EDCB9", "1HGFCB9", "347", "1BCD58I", "C89", "12389", "1456789", "1F9", "5HGF9", "1B589", "1D789", "ED789", "1FED789", "169", "749", "CBG85", "GD9", "3DI", "GF389", "12349", "12ED9", "34569", "EDC", "1234567", "58I", "129FC", "A89", "1D5", "I8GFEDC")
    #[1] 306 115

    ## Function to remove leading 'X' if present
    #remove_leading_x <- function(x) {
    #  ifelse(substr(x, 1, 1) == "X", substr(x, 2, nchar(x)), x)
    #}

    ## Apply the function to column names
    #colnames(df_barcode) <- remove_leading_x(colnames(df_barcode))

    df_majorbarcode <- df_barcode[,c("barcode", barcode_name)]; dim(df_majorbarcode)
    colnames(df_majorbarcode) <- c("barcode", "has_barcode")
    # [1] 306   2

    ### Get new circle data formated dataframe
    # Get barcode information
    df_adding_bin <- as.data.frame(df_majorbarcode$has_barcode)
    colnames(df_adding_bin) <- barcode_name
    rownames(df_adding_bin) <- df_majorbarcode$barcode

    # fixed spot id
    fixed_spotid <- rownames(input_table); length(fixed_spotid)
    # [1] 431
    # 提取 df_adding_bin 中在 fixed_spotid 中存在的行
    intersect_spots <- intersect(fixed_spotid, rownames(df_adding_bin))
    df_intersect <- as.data.frame(df_adding_bin[intersect_spots, ])
    rownames(df_intersect) <- intersect_spots
    colnames(df_intersect) <- colnames(df_adding_bin)
    # 找出 fixed_spotid 中但不在 df_adding_bin 行名中的 spot id
    filled_spots <- setdiff(fixed_spotid, rownames(df_adding_bin))
    # 创建一个新数据框，行名为 filled_spots，列名与 df_adding_bin 相同，内容填充为 3
    df_filled <- data.frame(matrix(3, nrow = length(filled_spots), ncol = ncol(df_adding_bin)))
    rownames(df_filled) <- filled_spots
    colnames(df_filled) <- colnames(df_adding_bin)

    # 将 df_existing 和 df_filled 合并，形成新的 df_adding_bin
    df_adding_bin_updated <- rbind(df_intersect, df_filled)

    ### Generate new circle data for plotting
    adding_long2 <- df_adding_bin_updated %>% as_tibble(rownames = "barcode") %>%
        tidyr::pivot_longer(-barcode,names_to = "creloxPbarcode",values_to = "input")

    ### plotting
    # color
    barcode_types <- as.character(unique(df_adding_bin_updated[,1]))
    adding_colors2 = dplyr::case_when(
        endsWith(barcode_types, "0") ~ "#EEEBDC",
        endsWith(barcode_types, "1") ~ "#514549",
    )
    names(adding_colors2) = barcode_types

    ##### plot adding and barcode
    # circos
    p6_merged <- p3_heatmap + new_scale_fill() +
        geom_fruit_list(
            geom_fruit(data = adding_long1,
                       geom = geom_tile, 
                       pwidth = 0.2,
                       width = 0.2,
                       mapping = aes(y = barcode, 
                                     x = factor(mutation),
                                     fill = as.character(input)),
                       offset = 0.015,
                       color = "white",
                       size = 0.2,
                       axis.params = list(axis = "x",
                                          line.size = NA,
                                          line.color = NA,
                                          text.size = 0,
                                          vjust = 0.5,
                                          hjust = 0.2)),
            scale_fill_manual(values = adding_colors1, 
                              breaks = names(adding_colors1),
                              guide = guide_legend(title = "Mutation types", 
                                                   title.position = "top",
                                                   override.aes = list(size = 10), 
                                                   title.theme = title_theme))
        ) + new_scale_fill() +
        geom_fruit_list(
            geom_fruit(data = adding_long2,
                       geom = geom_tile, 
                       pwidth = 0.2,
                       width = 0.2,
                       mapping = aes(y = barcode, 
                                     x = factor(creloxPbarcode),
                                     fill = as.character(input)),
                       offset = 0.02,
                       color = "white",
                       size = 0.2,
                       axis.params = list(axis = "x",
                                          line.size = NA,
                                          line.color = NA,
                                          text.size = 0,
                                          vjust = 0.5,
                                          hjust = 0.2)),
            scale_fill_manual(values = adding_colors2, 
                            breaks = names(adding_colors2),
                            guide = guide_legend(title = "cre-loxP barcode", 
                                                 title.position = "top",
                                                 override.aes = list(size = 10), 
                                                 title.theme = title_theme))
        )

    # svg_filename_p6_merged <- str_c(paste0(outputpath,"/p6_merged.circle_tree_output_as_point.", pdf_lastfix),".svg")
    # svg(file = svg_filename_p6_merged, width = plot_width*3, height = plot_height*2)
    # p6_merged
    # dev.off()
    # ggsave(file = svg_filename_p1_tree, width = plot_width*2.5, height = plot_height*2, plot = p1_tree)

}


##### legend mut
heatmap_mut_legend_colnum <- ceiling(length(pivot_muts)/15)
if(heatmap_mut_legend_colnum <= 2){
    heatmap_mut_legend_colnum = heatmap_mut_legend_colnum
} else {
    heatmap_mut_legend_colnum = 3
}
# write.table(for_sort_symbol, str_c(outputpath, "/mutation_list_in_order.for_sort_symbol.txt"), sep="\t")
p_mut <- ggplot(for_sort_symbol, aes(x=mutants, y=mutants_symbol, color=mutants)) +
    geom_text(aes(label = mutants_symbol)) +
    scale_color_manual(values = rep("#3c3c3c", length(pivot_muts_sub)),
                       breaks = pivot_muts_sub,
                       guide = guide_legend(title = "Mutations Symbol", 
                                            title.position = "top",
                                            override.aes = list(size = 8, label = for_sort_symbol$mutants_symbol),
                                            ncol = heatmap_mut_legend_colnum, 
                                            title.theme = title_theme)) +
    theme(legend.position = "right", 
          legend.text = element_text(size = 20), 
          legend.title = element_blank())

legend_mut <- get_legend(p_mut)

# svg(file = str_c(outputpath, "/legend_components.mut.svg"), width = 10, height = 5)
# print(plot_grid(legend_mut$grobs[[1]]))
# dev.off()
# ggsave(file = svg_filename_p1_tree, width = plot_width*2.5, height = plot_height*2, plot = p1_tree)


##### flipping count for all mutations
total_flipping <- read.table(total_flipping_file, header = TRUE, sep = "\t")
colnames(total_flipping) <- c("total_flipping_0_to_1", "total_flipping_1_to_0", "total_flipping_NA_to_0", "total_flipping_NA_to_1")
# Create the data frame for filpping count legend
data_total_flipping <- tibble(
    type = c("False negative", "False positive", "Fill NA with 0", "Fill NA with 1"),
    category = c("0>1", "1>0", "NA>0", "NA>1"),
    count = c(
        total_flipping$total_flipping_0_to_1 %||% 0,
        total_flipping$total_flipping_1_to_0 %||% 0,
        total_flipping$total_flipping_NA_to_0 %||% 0,
        total_flipping$total_flipping_NA_to_1 %||% 0
    )
)
data_total_flipping$labels <- paste(data_total_flipping$category, ":", data_total_flipping$count)
data_total_flipping$count <- as.factor(data_total_flipping$count)
# generate legend plot
p_total_flipping <- ggplot(data_total_flipping, aes(x=labels, y=category, color=labels)) + 
    geom_text(aes(label = labels)) +
    scale_color_manual(values = rep("#3c3c3c", length(unique(data_total_flipping$labels))),
                       breaks = unique(data_total_flipping$labels),
                       guide = guide_legend(title = "Flipping count of all mutations",
                                            title.position = "top",
                                            override.aes = list(size = 8, label = data_total_flipping$type),
                                            ncol = 1,
                                            title.theme = title_theme)) +
    theme(
        legend.position = "right",  # 设置 legend 位置
        legend.text = element_text(size = 20),  # 设置 legend 中的文本字体大小
        legend.title = element_text()  # 设置 legend 标题字体
    )
# extract legend
legend_total_flipping <- get_legend(p_total_flipping)
# # plot legend
# svg(file = str_c(outputpath, "/legend_components.total_flipping_count.svg"), width = 10, height = 5)
# print(plot_grid(legend_total_flipping$grobs[[1]]))
# dev.off()
# ggsave(file = svg_filename_p1_tree, width = plot_width*2.5, height = plot_height*2, plot = p1_tree)


##### flipping count for target_mut
if(target_mut=="no"){

    print("No target_mut, so no flipping_count.")

} else {

    df_features <- read.table(features_file, header = TRUE, sep = "\t", row.names = 1)
    flipping_columns <- grep("flipping", colnames(df_features), value = TRUE)
    df_filpping <- df_features[, flipping_columns]
    target_flipping <- df_filpping[target_mut, ]
    # Create the data frame for filpping count legend
    data_target_flipping <- tibble(
    type = c("False negative", "False positive", "Fill NA with 0", "Fill NA with 1"),
    category = c("0>1", "1>0", "NA>0", "NA>1"),
    count = c(
        target_flipping$tree_flipping_False_Negative %||% 0,
        target_flipping$tree_flipping_False_Positive %||% 0,
        target_flipping$tree_flipping_NA_to_0 %||% 0,
        target_flipping$tree_flipping_NA_to_1 %||% 0
    )
    )
    data_target_flipping$labels <- paste(data_target_flipping$category, ":", data_target_flipping$count)
    data_target_flipping$count <- as.factor(data_target_flipping$count)
    # generate legend plot
    p_target_flipping <- ggplot(data_target_flipping, aes(x=labels, y=category, color=labels)) + 
        geom_text(aes(label = labels)) +
        scale_color_manual(values = rep("#3c3c3c", length(unique(data_target_flipping$labels))),
                           breaks = unique(data_target_flipping$labels),
                           guide = guide_legend(title = str_c("Flipping count of ", target_mut),
                                                title.position = "top",
                                                override.aes = list(size = 8, label = data_target_flipping$type),
                                                ncol = 1,
                                                title.theme = title_theme)) +
        theme(
            legend.position = "right",  # 设置 legend 位置
            legend.text = element_text(size = 20),  # 设置 legend 中的文本字体大小
            legend.title = element_text()  # 设置 legend 标题字体
        )
    # extract legend
    legend_target_flipping <- get_legend(p_target_flipping)
    # # plot legend
    # svg(file = str_c(outputpath, "/legend_components.target_flipping_count.svg"), width = 10, height = 5)
    # print(plot_grid(legend_target_flipping$grobs[[1]]))
    # dev.off()
    # ggsave(file = svg_filename_p1_tree, width = plot_width*2.5, height = plot_height*2, plot = p1_tree)

}


##### Get ordered metadata for heatmap
### The order of cells/spots
# 提取绘制后的树的坐标数据
circostree <- p1_tree$data
# 选择叶节点（tiplab的标签）
tip_data <- circostree[circostree$isTip,]
# 获取每个叶节点的角度信息（radian）
tip_data$angle <- atan2(tip_data$y, tip_data$x)
# 按角度从小到大排序，逆时针排序
sorted_tips <- tip_data[order(tip_data$angle),]
# 获取排好序的 cell ids
tip_label_order <- sorted_tips$label

### The order of muts
# print(for_sort_symbol)
cols_to_sort <- for_sort_symbol$mutants

### ordered_metadata_for_heatmap
metadata_for_heatmap_temp <- cf_tableno0
# reorder columns
sorting_metadata_for_heatmap <- metadata_for_heatmap_temp[, c("cellIDxmutID", cols_to_sort)]
# reorder rows
sorting_metadata_for_heatmap$cellIDxmutID <- factor(sorting_metadata_for_heatmap$cellIDxmutID, levels = tip_label_order)
ordered_metadata_for_heatmap <- sorting_metadata_for_heatmap[order(sorting_metadata_for_heatmap$cellIDxmutID), ]
# add mutid_in_heatmap
ordered_metadata_for_heatmap[] <- lapply(ordered_metadata_for_heatmap, as.character)
new_row <- c("mutid_in_heatmap", for_sort_symbol$mutants_symbol)
ordered_metadata_for_heatmap_final <- rbind(new_row, ordered_metadata_for_heatmap)

# all(as.character(ordered_metadata_for_heatmap_final$cellIDxmutID) == c("mutid_in_heatmap", tip_label_order))
# all(as.character(colnames(ordered_metadata_for_heatmap_final)) == c("cellIDxmutID", cols_to_sort))
write.table(ordered_metadata_for_heatmap_final, str_c(outputpath, "/ordered_metadata_for_heatmap.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


##### Generate mutation id corresponding to the ordered tip.lables 
# 创建一个空的数据框用于存储结果
df_muts_corresponding_to_ordered_tiplabels <- data.frame(cellIDxmutID = ordered_metadata_for_heatmap_final$cellIDxmutID[2:nrow(ordered_metadata_for_heatmap_final)], 
                     mutation_in_heatmap = character(nrow(ordered_metadata_for_heatmap_final)-1))

# 对每一行进行操作
for (i in 2:nrow(ordered_metadata_for_heatmap_final)) {
    # 获取当前行的突变数据
    row_data <- ordered_metadata_for_heatmap_final[i, cols_to_sort]

    # 找到值为1的突变的列索引
    ones_indices <- which(row_data == 1)

    if (length(ones_indices) > 0) {
        # 找到最小的索引
        max_index <- max(ones_indices)

        # 获取对应的突变名字（cols_to_sort 中最小索引的列名）
        mutation_least_cells <- cols_to_sort[max_index]

        # 将结果保存到 df_muts_corresponding_to_ordered_tiplabels 数据框中
        df_muts_corresponding_to_ordered_tiplabels$mutation_in_heatmap[i-1] <- mutation_least_cells
    } else {
        # 如果没有1，则设置为NA或其他默认值
        df_muts_corresponding_to_ordered_tiplabels$mutation_in_heatmap[i-1] <- NA
    }
}

# 保存结果
write.table(df_muts_corresponding_to_ordered_tiplabels, str_c(outputpath, "/df_muts_corresponding_to_ordered_tiplabels_by_anticlockwise.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)




##### Finally plot -------------------------------------------------------
# Create the main plot without the legends
# adding_mut_file and barcode_name don't appear simultaneously
if (trimws(adding_mut_file)=="no" && trimws(barcode_name)=="no") {
    main_plot <- p3_heatmap + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    if (trimws(target_mut)=="no") {
        svg_filename <- str_c(paste0(dir,"/No_target.circle_tree_output_as_point.", pdf_lastfix),".svg")
    } else {
        svg_filename <- str_c(paste0(dir,"/target_", target_mut,".circle_tree_output_as_point.", pdf_lastfix),".svg")
    }
} else if (trimws(barcode_name)=="no" && trimws(adding_mut_file)!="no") {
    main_plot <- p4_adding + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    if (trimws(target_mut)=="no") {
        svg_filename <- str_c(paste0(dir,"/adding_", adding_mutid, ".No_target.circle_tree_output_as_point.", pdf_lastfix),".svg")
    } else {
        svg_filename <- str_c(paste0(dir,"/adding_", adding_mutid, ".target_", target_mut,".circle_tree_output_as_point.", pdf_lastfix),".svg")
    }
} else if (trimws(adding_mut_file)=="no" && trimws(barcode_name)!="no") {
    main_plot <- p5_barcode + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    if (trimws(target_mut)=="no") {
        svg_filename <- str_c(paste0(dir,"/barcode_", barcode_name, ".No_target.circle_tree_output_as_point.", pdf_lastfix),".svg")
    } else {
        svg_filename <- str_c(paste0(dir,"/barcode_", barcode_name, ".target_", target_mut,".circle_tree_output_as_point.", pdf_lastfix),".svg")
    }
} else if (trimws(adding_mut_file)!="no" && trimws(barcode_name)!="no") {
    main_plot <- p6_merged + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    if (trimws(target_mut)=="no") {
        svg_filename <- str_c(paste0(dir,"/barcode_", barcode_name, ".adding_", adding_mutid, ".No_target.circle_tree_output_as_point.", pdf_lastfix),".svg")
    } else {
        svg_filename <- str_c(paste0(dir,"/barcode_", barcode_name, ".adding_", adding_mutid, ".target_", target_mut,".circle_tree_output_as_point.", pdf_lastfix),".svg")
    }
}


##### Plot root node
fixed_root_color <- "#2B2B2B" # 指定固定的 root 颜色
# 获取现有颜色
if(target_mut=="no"){
    used_colors <- get_used_colors(
        main_plot,
        other_colors = c(
            No_target_clone_mut_color, # 目标克隆颜色
            color_muts,             # 突变颜色
            subclone_color,         # 子克隆颜色
            anno_color
        )
    )
} else {
    used_colors <- get_used_colors(
        main_plot,
        other_colors = c(
            target_clone_mut_color, # 目标克隆颜色
            color_muts,             # 突变颜色
            subclone_color,         # 子克隆颜色
            anno_color
        )
    )
}
# 检查固定颜色是否唯一
check_unique_color(fixed_root_color, used_colors)
# 绘制 Root 节点
na_node_data <- subset(main_plot$data, is.na(label))
if (nrow(na_node_data) > 0) {
    main_plot <- main_plot + new_scale_fill() + 
        geom_point(data = na_node_data, aes(x = x, y = y), 
                   inherit.aes = FALSE, 
                   color = fixed_root_color, size = 5)
}

# svg_filename_main_plot <- str_c(paste0(outputpath,"/main_plot.circle_tree_output_as_point.", pdf_lastfix),".svg")
# ggsave(file = svg_filename_main_plot, width = plot_width*2.5, height = plot_height*2, plot = main_plot)
# ggsave(file = str_c(paste0(outputpath,"/main_plot.circle_tree_output_as_point.", pdf_lastfix),".pdf"), width = plot_width*2.5, height = plot_height*2, plot = main_plot)


##### Draw circos and legend separately
### plot only circos main body
p_final <- main_plot + theme(legend.position = "none")

# ggsave(file = svg_filename, width = plot_width*2, height = plot_height*2, plot=p_final)
pdf_filename <- gsub('svg', 'pdf', svg_filename)
ggsave(file = pdf_filename, width = plot_width*2, height = plot_height*2, plot=p_final)

# ggsave(file = svg_filename, width = plot_width*2.5, height = plot_height*2.5, plot=p_final, limitsize = FALSE)
# pdf_filename <- gsub('svg', 'pdf', svg_filename)
# ggsave(file = pdf_filename, width = plot_width*2.5, height = plot_height*2.5, plot=p_final, limitsize = FALSE)


### plot legend for circos
p_legend_circos <- plot_grid(get_legend(main_plot))

svg_filename_legend_circos <- str_c(outputpath, "/legend_components.circos_annotation.svg")
pdf_filename_legend_circos <- gsub('svg', 'pdf', svg_filename_legend_circos)
# ggsave(file = svg_filename_legend_circos, width = plot_width, height = plot_height, plot = p_legend_circos)
ggsave(file = pdf_filename_legend_circos, width = plot_width, height = plot_height, plot = p_legend_circos)


### plot all the legends in an independent picture.
p_legend_total_flipping <- plot_grid(legend_total_flipping$grobs[[1]])

svg_filename_legend_filpping <- str_c(outputpath, "/legend_components.total_flipping_count.svg")
pdf_filename_legend_filpping <- gsub('svg', 'pdf', svg_filename_legend_filpping)
# ggsave(file = svg_filename_legend_filpping, width = plot_width, height = plot_height, plot = p_legend_total_flipping)
ggsave(file = pdf_filename_legend_filpping, width = plot_width, height = plot_height, plot = p_legend_total_flipping)



















# ############################################################
# ############################################################
# ############################################################

# # annovar_files = c(
# #   "/storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/04.mutations/P6_ST_vis_rep1/spatial_lineager/new_model_250921/new_hFDR/all_RF_filter_hFDR49_less_08_feature_no_AG.filter_poly.filter_PON.anno.variant_function", 
# #   "/storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/04.mutations/P6_ST_vis_rep2/spatial_lineager/new_model_250921/new_hFDR/all_RF_filter_hFDR49_less_08_feature_no_AG.filter_poly.filter_PON.anno.exonic_variant_function")


# # Mutation ID: chr1:39034563,T>A
# # Cell ID: CGTGCCGACATTTGT-1
# # Gene ID: NDUFS5
# # Genome region: UTR3
# # Function prediction: synonymous

# library(dplyr)
# library(tidyr)

# #' 处理遗传变异数据格式
# #'
# #' @param data 输入数据框
# #' @param identifier_col identifier列的名称
# #' @param region_col region_and_function列的名称
# #' @param gene_col gene_id列的名称
# #' @return 处理后的数据框
# process_variant_data <- function(data, identifier_col = "identifier", region_col = "region_and_function", gene_col = "gene_id") {
  
#   # 检查必要的列是否存在
#   required_cols <- c(identifier_col, region_col, gene_col)
#   missing_cols <- setdiff(required_cols, names(data))
#   if (length(missing_cols) > 0) {
#     stop("缺少必要的列: ", paste(missing_cols, collapse = ", "))
#   }
  
#   result <- data %>%
#     # 重命名列以便处理
#     rename(
#       identifier = !!sym(identifier_col),
#       region_function = !!sym(region_col),
#       gene = !!sym(gene_col)
#     ) %>%
#     # 拆分identifier列
#     separate(identifier, into = c("chr", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>%
#     # 创建新格式的字符串
#     mutate(new_identifier = paste0(chr, ":", pos, ",", ref, ">", alt)) %>%
#     # 处理region_and_function列
#     mutate(
#       region = ifelse(grepl("\\(", region_function),
#                      sub("\\(.*", "", region_function),
#                      region_function),
#       function_col = ifelse(grepl("\\(", region_function),
#                            gsub(".*\\(|\\)", "", region_function),
#                            "unknown")
#     ) %>%
#     # 恢复原始列名并选择输出列
#     rename(
#       !!identifier_col := "identifier",
#       !!region_col := "region_function", 
#       !!gene_col := "gene"
#     ) %>%
#     # 选择输出列，保留原始 identifier 列
#     select(identifier, new_identifier, !!gene_col, region, function_col)
  
#   return(result)
# }


# annovar_files = c("/storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/04.mutations/P6_ST_vis_rep1/spatial_lineager/new_model_250921/new_hFDR/all_RF_filter_hFDR49_less_08_feature_no_AG.filter_poly.filter_PON.anno_info.txt", "/storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/04.mutations/P6_ST_vis_rep2/spatial_lineager/new_model_250921/new_hFDR/all_RF_filter_hFDR49_less_08_feature_no_AG.filter_poly.filter_PON.anno_info.txt")

# # 初始化一个空的数据框来存储所有结果
# df_annovar_combined <- data.frame()

# # 循环处理每个文件
# for(annovar_file in annovar_files){
#   cat("正在处理文件:", annovar_file, "\n")
  
#   # 读取数据
#   df_annovar <- read.table(annovar_file, sep="\t", header=FALSE)
#   colnames(df_annovar) <- c("identifier", "region_and_function", "gene_id")
  
#   # 处理数据
#   df_annovar_split <- process_variant_data(df_annovar)
  
#   # 如果是第一个文件，直接赋值
#   if (nrow(df_annovar_combined) == 0) {
#     df_annovar_combined <- df_annovar_split
#   } else {
#     # 合并并去重：保留原本的行，新的行加到后面，然后去重
#     df_annovar_combined <- bind_rows(df_annovar_combined, df_annovar_split) %>%
#       distinct(new_identifier, .keep_all = TRUE)
#   }
  
#   cat("当前数据框行数:", nrow(df_annovar_combined), "\n")
# }
# colnames(df_annovar_combined) <- c("identifier", "mutation_id", "gene_id", "genome_region", "function_prediction")

# # 查看最终结果
# cat("最终去重后的数据框:\n")
# # print(df_annovar_combined)
# cat("总行数:", nrow(df_annovar_combined), "\n")
# write.table(df_annovar_combined, file = str_c(outputpath, "/identifier.anno_info.txt"), sep = "\t", quote = FALSE)
# # write.table(df_annovar_combined, file = str_c("/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_P6_merged.96true/mutation_integrator/phylo/circos_plot/identifier.anno_info.txt"), sep = "\t", quote = FALSE)




##### Time #####
end_time <- Sys.time()
print(str_c("End time: ", end_time))
print(str_c("Program finished in ", as.character(round((end_time-start_time), 4)), " seconds"))



