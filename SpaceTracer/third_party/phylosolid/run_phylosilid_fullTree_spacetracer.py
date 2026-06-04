#!/usr/bin/env python3

import sys
from pathlib import Path

# Add project root to Python path
# _project_root = str(Path(__file__).parent.parent)
# if _project_root not in sys.path:
#     sys.path.insert(0, _project_root)

# Date: 2025/09/16
# Update: 2025/10/13
# Author: Qing
# Work: Test scaffold builder


##### Time #####
import time
start_time = time.perf_counter()




################################################################################################
########################################## PhyloSOLID #########################################
################################################################################################
import os
import logging
import copy
import random
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import deepcopy

logger = logging.getLogger(__name__)

from src.data_loader import load_all
# from src.scrna_classifier import real_time_classifier_predict
from src.germline_filter import identify_germline_variants
from src.germline_filter import *
from src.scaffold_builder import build_scaffold_tree
from src.scaffold_builder import *
from src.mutation_integrator import *

# 设置所有随机种子
RANDOM_SEED = 42  # 可以任意指定
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# 如果有使用其他库也设置种子
try:
    import torch
    torch.manual_seed(RANDOM_SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(RANDOM_SEED)
        torch.cuda.manual_seed_all(RANDOM_SEED)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
except ImportError:
    pass

# ------------------------------
# 配置 logging
# ------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# ------------------------------
# 项目参数与文件路径设置
# ------------------------------

import multiprocessing as mp
import argparse
from argparse import ArgumentParser
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sampleid", default="sampleid", type=str, help="The sampleid you can set and check.")
parser.add_argument("-i", "--inputpath", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_100k.30true/data", type=str, help="The inputpath contains the preprocessing results from raw posterior-reads data.")
parser.add_argument("-o", "--outputpath", default="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_100k.30true/results", type=str, help="The outputpath you want to save results.")
parser.add_argument("-c", "--celltype_file", default=None, type=str, help="The celltype_file you should provide. If you can't generate this file, please set 'None'.")
parser.add_argument("--features_file", default=None, type=str, help="The features_file you should provide. If you can't generate this file, please set 'None'.")
parser.add_argument("--is_predict_germ", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to predict germline mutations.")
parser.add_argument("--is_detect_passtree_by_dp", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to run Dynamic programing step.")
parser.add_argument("--is_filter_quality", default="yes", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to filter mutations in scaffold steps by coverage quality.")
args = parser.parse_args()


# get parameters
sampleid = args.sampleid
inputpath = args.inputpath
outputpath = args.outputpath
celltype_file = args.celltype_file
features_file = args.features_file


outputpath_classifier = os.path.join(outputpath, "classifier_filter")
outputpath_germline = os.path.join(outputpath, "germline_filter")
outputpath_scaffold = os.path.join(outputpath, "scaffold_builder")
outputpath_full = os.path.join(outputpath, "mutation_integrator")
# 递归创建目录，如果目录已存在也不会报错
os.makedirs(outputpath_classifier, exist_ok=True)
os.makedirs(outputpath_germline, exist_ok=True)
os.makedirs(outputpath_scaffold, exist_ok=True)
os.makedirs(outputpath_full, exist_ok=True)

is_predict_germ = args.is_predict_germ  ##"no"  ## "yes" or "no"
is_detect_passtree_by_dp = args.is_detect_passtree_by_dp  ##"no"  ## "yes" or "no"
is_filter_quality = args.is_filter_quality  ##"no"  ## "yes" or "no"

# Display parameters for verification
logger.info(f"sampleid: {sampleid}")
logger.info(f"inputpath: {inputpath}")
logger.info(f"outputpath: {outputpath}")
logger.info(f"celltype_file: {celltype_file}")
logger.info(f"features_file: {features_file}")
logger.info(f"is_predict_germ: {is_predict_germ}")
logger.info(f"is_detect_passtree_by_dp: {is_detect_passtree_by_dp}")
logger.info(f"is_filter_quality: {is_filter_quality}")


# ------------------------------
# 参数设置
# ------------------------------

# tumor_cell_of_origin_mutations = ["chr1_46278500_C_T", "chr16_67944430_G_C", "chr17_29572953_G_A"]
# for i in tumor_cell_of_origin_mutations:
#     count_result = count_list(I_attached[i])
#     NA_prop = count_result['NA']/(count_result['NA']+count_result[0]+count_result[1])
#     print(count_result)
#     print(NA_prop)

# # {'NA': 1032, 0.0: 149, 1.0: 72}
# # {'NA': 1128, 0.0: 95, 1.0: 30}
# # {'NA': 1146, 0.0: 65, 1.0: 42}

# Default parameters matching Methods Section 3
SETTING_PARAMS = {
    "models_path": "phylosolid/models/scdna",
    
    # 1 data loader
    "p_thresh": 0.5,
    
    # 2 germline filter
    "mcf_cutoff": 0.05,
    "mcn_cutoff": 5,
    
    # Pairwise correlation criteria (Section 2.1)
    "pair_N11_min": 0,                 # minimum N11 for correlation
    "jaccard_thresh": 0.2,             # Jaccard threshold for strong correlation
    "jaccard_low": 0.1,                # lower Jaccard bound for parent-child
    "fraction_parent_child_thresh": 0.9, # fraction threshold for parent-child
    
    # 3.1 Initial filtration
    "posterior_threshold": 0.5,
    "maf_max_threshold": 0.3,
    "maf_mean_threshold": 0.1,
    
    # 3.2 Coverage-based filtration
    "na_prop_thresh_global": 0.95,      # p_NA(j) ≤ 0.9 (present in >10% cells)
    "cv_thresh": 6.0,                   # CV_j < 2
    
    # 3.3 Consensus correlation graph
    "consensus_runs": 100,             # number of randomized runs
    "consensus_clone_freq_thresh": 0.1, # minimum clone frequency threshold
    "resolution_of_graph": 1,
        
    # 3.4 Penalty-based placement
    "general_weight_NA": 0.001,        # NA imputation penalty parameter
    "fnfp_ratio": 0.1,                 # false negative/false positive ratio
    "phi": 1.0,                        # BIC complexity parameter
    
    # 4.1 Dynamic programming
    "pass_tree_cutoff": 0.9,
    "unpass_tree_cutoff": 0.1,
    
    # 4.2 fp_ratio and fn_ratio
    "fp_ratio_cutoff_across_tree": 0.2,
    "fn_ratio_cutoff_across_tree": 0.9,
    "fp_ratio_cutoff_within_subclone": 0.1,
    "fp_ratio_persite_cutoff": 0.1, 
    "fp_count_persite_cutoff": 0,
    
    "fp_ratio_per_mutation_cross_all_cells_cutoff": 0.2,
    "fp_count_per_mutation_cross_all_cells_cutoff": 5,
    "fp_ratio_per_cell_cross_all_muts_cutoff": 0.5,
    
    "intersection_vs_fn_flipping_ratio_cutoff": 0.2,
    "intersection_cell_count_on_mutation_cutoff": 5,
    "intersection_cell_ratio_on_mutation_cutoff": 0.2,
    "intersection_count_per_cells_cutoff": 1,
    "flipping_count_fn_per_cells_cutoff": 1,
    "flipping_to_1_count_per_cells_cutoff": 2
}

params = SETTING_PARAMS




# ------------------------------
# Step 1: Load data
# ------------------------------
logger.info("===== Step1: Loading data ...")
data = load_all(inputpath)
P_raw, V_raw, C_raw, A_raw = data["P"], data["V"], data["C"], data["A"]
df_features = data['features']
df_reads = data['df_reads']
I_raw = build_binary_I(P_raw, V_raw, C_raw, params["p_thresh"])
all_mutations = list(I_raw.columns)

P, V, C, A, I = P_raw.copy(), V_raw.copy(), C_raw.copy(), A_raw.copy(), I_raw.copy()

df_features_new, empty_mutations = update_features_matrix(I, df_reads, df_features, params["mcf_cutoff"])

logger.info(f"Loaded data: {len(P)} cells, {len(I_raw.columns)} mutations")

##### Output binary matrix with NA=0
I_raw_withNA0 = I_raw.replace({np.nan: 0}).astype(int)
I_raw_withNA0.to_csv(os.path.join(outputpath, "I_raw_withNA0.txt"), sep="\t")
I_raw_withNA0_T = I_raw_withNA0.T
I_raw_withNA0_T.to_csv(os.path.join(outputpath, "I_raw_withNA0_T.txt"), sep="\t")

I_raw_withNA3 = I_raw.replace({np.nan: 3}).astype(int)
I_raw_withNA3.to_csv(os.path.join(outputpath, "I_raw_withNA3.txt"), sep="\t")
I_raw_withNA3_T = I_raw_withNA3.T
I_raw_withNA3_T.to_csv(os.path.join(outputpath, "I_raw_withNA3_T.txt"), sep="\t")

##### Output for BSCITE bulk input
# 获取 pseudo bulk 数据
df_corrected = df_reads.copy()

# 分割数据并转换为数值
def split_and_sum(series):
    before_sum = 0
    after_sum = 0
    
    for value in series:
        if pd.notna(value) and '/' in str(value):
            parts = str(value).split('/')
            if len(parts) == 2:
                try:
                    before_sum += int(parts[0])
                    after_sum += int(parts[1])
                except ValueError:
                    continue
                    
    return f"{before_sum}/{after_sum}"

# 对每一列应用函数
for col in df_corrected.columns:
    df_corrected.loc['pseudo_bulk', col] = split_and_sum(df_corrected[col])

# 产生 bulk input file
output_bulk_data = []

# 选择三个样本：pseudo_bulk, bulk, 和一个随机细胞
sample1 = 'pseudo_bulk'
sample2 = 'bulk'

# 从剩下的细胞中随机选一个（排除pseudo_bulk和bulk）
remaining_cells = [idx for idx in df_corrected.index if idx not in ['pseudo_bulk', 'bulk']]
sample3 = random.choice(remaining_cells)

print(f"使用样本: {sample1}, {sample2}, {sample3}")

for i, snp_name in enumerate(df_corrected.columns):
    # 获取三个样本的数据
    sample1_value = df_corrected.loc[sample1, snp_name]
    sample2_value = df_corrected.loc[sample2, snp_name] 
    sample3_value = df_corrected.loc[sample3, snp_name]
    
    samples_mutant = []
    samples_reference = []
    
    # 处理每个样本的数据
    for value in [sample1_value, sample2_value, sample3_value]:
        if pd.isna(value):
            # 如果数据缺失，使用0/0
            samples_mutant.append(0)
            samples_reference.append(0)
        elif '/' in str(value):
            mutant, reference = str(value).split('/')
            samples_mutant.append(int(mutant))
            samples_reference.append(int(reference))
        else:
            # 如果格式不对，使用0/0
            samples_mutant.append(0)
            samples_reference.append(0)
    
    # 组合三个样本的数据
    mutant_counts = f"{samples_mutant[0]};{samples_mutant[1]};{samples_mutant[2]}"
    reference_counts = f"{samples_reference[0]};{samples_reference[1]};{samples_reference[2]}"
    
    # 获取染色体和位置信息
    chromosome, position = get_random_chromosome_position(snp_name)
    
    output_bulk_data.append({
        'ID': f'mut{i}',
        'Chromosome': chromosome,
        'Position': position,
        'MutantCount': mutant_counts,
        'ReferenceCount': reference_counts
    })

output_bulk_df = pd.DataFrame(output_bulk_data)
output_bulk_df.to_csv(os.path.join(outputpath, "input_BULK.txt"), sep='\t', index=False)




# ------------------------------
# Step 2: Classifier for identifying mosaic mutations
# ------------------------------
logger.info("===== Step2: Classifier ...")

models_path = params['models_path']

##### Load Celltype Data
if features_file is None or features_file == "None" or features_file == "no":
    candidate_mutations = list(I_raw.columns)
else:
    df_for_classifier = pd.read_csv(features_file, sep="\t")
    if 'uniqueid' in df_for_classifier.columns:
        df_for_classifier['mutation_id'] = df_for_classifier['uniqueid'].str.split('.').str[0]
    else:
        df_for_classifier['uniqueid'] = df_for_classifier['identifier']+"."+sampleid
        df_for_classifier['mutation_id'] = df_for_classifier['uniqueid'].str.split('.').str[0]
    
    results_of_classifier = real_time_classifier_predict(df_for_classifier, sampleid, outputpath_classifier)
    
    df_results_of_classifier = results_of_classifier.copy()
    df_candidate_mosaic = df_results_of_classifier.loc[df_results_of_classifier['predicted_label']=='mosaic',:]
    candidate_mutations = [i for i in list(set(list(df_candidate_mosaic['mutation_id']))) if i in list(I_raw.columns)]

logging.info(f"Predicted {len(candidate_mutations)} candidate mosaic mutations")


P_candidate = P[candidate_mutations].copy()
V_candidate = V[candidate_mutations].copy()
A_candidate = A[candidate_mutations].copy()
C_candidate = C[candidate_mutations].copy()
I_candidate = I[candidate_mutations].copy()
df_reads_candidate = df_reads[candidate_mutations].copy()




# ------------------------------
# Step 3. Germline filtering
# ------------------------------
logger.info("===== Step3: Predict germline mutations ...")

if is_predict_germ=="yes":
    logging.info("Running germline filtering ...")
    stats_df, germline_mutations = identify_germline_variants(
        P=P, V=V, C=C, df_reads=df_reads, df_features_new=df_features_new, 
        p_thresh=params["p_thresh"],
        mcf_cutoff=params["mcf_cutoff"],
        mcn_cutoff=params["mcn_cutoff"],
        outputpath=outputpath_germline,
        sampleid=sampleid
    )
else:
    germline_mutations = set()

predicted_germline_mutations = list(germline_mutations)
rescued_germline_mutations = [mut for mut in predicted_germline_mutations 
                             if mut in df_features_new.columns 
                             and df_features_new.loc['mutant_cell_fraction_detected', mut] < 0.5]

removed_germline_mutations = [i for i in predicted_germline_mutations if i not in rescued_germline_mutations]


logging.info(f"Identified {len(predicted_germline_mutations)} germline variants")
# 100k: {'chr4_78610227_T_G', 'chr8_123681204_C_T', 'chr19_46838237_G_A'}
# 10k: 
# P6_merged: {'chr2_187412173_G_C', 'chr16_67944430_G_C', 'chr14_106276454_C_G', 'chr20_45109473_A_G'}
# 151674: set()
# P4_merged: set()
# P6_merged:  {'chr20_45109473_A_G', 'chr19_4815978_C_A', 'chr2_187412173_G_C', 'chr14_106276454_C_G', 'chr16_67944430_G_C'}


# heatmap 中展示鉴定结果
plot_heatmap_with_germline_mutations(I, predicted_germline_mutations, os.path.join(outputpath_germline, sampleid+".heatmap_with_predicted_germline_mutations_and_histograms.pdf"))
# predicted_germline_mutations = set()

# 去掉 predicted_germline_mutations 中的突变（只保留 somatic）
# removed_artifact_mutations = ['chr17_41690997_C_T', 'chr15_55317839_C_G']
# removed_artifact_mutations = ['chr17_41690997_C_T', 'chr15_55317839_C_G', 'chr5_139341854_G_T', 'chr13_44433526_T_C', 'chr7_24698525_G_T', 'chr16_57231300_G_A']
removed_artifact_mutations = []
somatic_mutations = [i for i in all_mutations if i not in removed_germline_mutations and i not in removed_artifact_mutations]
P_somatic = P[somatic_mutations].copy()
V_somatic = V[somatic_mutations].copy()
A_somatic = A[somatic_mutations].copy()
C_somatic = C[somatic_mutations].copy()
I_somatic = I[somatic_mutations].copy()
df_reads_somatic = df_reads[somatic_mutations].copy()

I_somatic_withNA3 = I_somatic.replace({np.nan: 3}).astype(int)
I_somatic_withNA3.to_csv(os.path.join(outputpath_germline, "I_somatic_withNA3.txt"), sep="\t")
df_features_new = add_mutation_proportions_to_features(df_features_new, I_somatic)




# ------------------------------
# Step 4: Scaffold builder
# ------------------------------
logger.info("===== Step4: Construct scaffold tree ...")

##### Load Celltype Data
if celltype_file is None or celltype_file == "None":
    barcodes = df_reads_somatic.index.tolist()  # 获取所有条形码
    df_celltype = pd.DataFrame({
        "barcode": barcodes,
        "cell_type": ["default_type"] * len(barcodes)
    })
else:
    df_celltype = pd.read_csv(celltype_file, sep="\t")

df_celltype.to_csv(os.path.join(outputpath_scaffold, "df_celltype.txt"), sep="\t")
logger.info(f"Celltype data loaded: {df_celltype.shape[0]} cells")


##### Scaffold builder
logging.info("Running scaffold building ...")
immune_mutations = ['chr14_105707803_C_A', 'chr14_106276422_A_G', 'chr14_106276454_C_G', 'chr22_22881588_G_A', 'chr22_22895482_G_A', 'chr22_22895504_C_A', 'chr22_22901226_C_G', 'chr22_22901169_G_C', 'chr22_22895709_C_G']
results_of_scaffold = build_scaffold_tree(
    P_somatic = P_somatic, 
    V_somatic = V_somatic, 
    A_somatic = A_somatic, 
    C_somatic = C_somatic, 
    I_somatic = I_somatic,
    df_reads_somatic = df_reads_somatic,
    df_features_new = df_features_new,
    params = params,
    is_filter_quality = is_filter_quality,
    outputpath = outputpath_scaffold,
    sampleid = sampleid,
    df_celltype = df_celltype,
    immune_mutations = immune_mutations
)


T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count, final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold, backbone_mutations, mutation_group, spots_to_split, group_mutations, remained_mutations, high_cv_mutations = results_of_scaffold
# scaffold_mutations = [i for i in initial_scaffold_mutations if i not in remained_mutations_by_scaffold_building]
scaffold_mutations = list(M_scaffold.columns)
non_scaffold_mutations = [i for i in somatic_mutations if i not in scaffold_mutations]


mutation_clones_for_scaffold = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_scaffold, I_somatic)
df_barcode_clones_for_scaffold = assign_clone_labels(M_scaffold, mutation_clones_for_scaffold)

df_barcode_clones_for_scaffold.to_csv(os.path.join(outputpath_scaffold, "df_barcode_clones_from_phylo_tree.csv"), sep=',', index=False)


print_tree(T_scaffold)
# └─ ROOT
#   └─ chr17_7578893_G_T
#     └─ chr8_80170791_C_T|chr7_20381814_G_T
#       └─ chr4_122860726_G_T
#         └─ chr15_50499489_A_G
#           └─ chr16_497748_G_T
#             └─ chr10_73254076_G_T
#       └─ chr21_46323662_C_T
#         └─ chr10_118753098_G_T
#         └─ chr16_15464093_A_G
#   └─ chr6_34246136_G_C
#     └─ chr1_45618093_G_T
#       └─ chr2_171687790_A_T
#         └─ chr21_44856223_C_A
#           └─ chr11_9425206_C_A
#             └─ chr17_50662360_G_T
#               └─ chr7_87196234_G_T
#                 └─ chr2_171324260_C_T
#                   └─ chr18_9135639_C_G
#                   └─ chr12_51215297_G_T
#   └─ chr22_37967634_G_T
#     └─ chr4_74388653_G_T
#       └─ chr17_43766282_C_A
#     └─ chr12_102013001_C_A
#   └─ chrX_41232675_T_A
#     └─ chr1_31260040_C_A
#       └─ chr7_132879457_T_A
#         └─ chr10_104003258_C_T


logging.info(f"Identified {len(scaffold_mutations)} scaffold variants")


logging.info(f"The number of all_mutations is: {len(all_mutations)}")
logging.info(f"The number of predicted_germline_mutations is: {len(predicted_germline_mutations)}")
logging.info(f"The number of somatic_mutations is: {len(somatic_mutations)}")
logging.info(f"The number of scaffold_mutations is: {len(scaffold_mutations)}")
logging.info(f"The number of non_scaffold_mutations is: {len(non_scaffold_mutations)}")




# ------------------------------
# Step 5: DP pass tree & classification
# ------------------------------
logger.info("===== Step5: Run DP for passtree ...")

##### DP
if is_detect_passtree_by_dp=="yes":
    logging.info("Running dynamic programming ...")
    
    pass_tree_cutoff=params['pass_tree_cutoff']
    unpass_tree_cutoff=params['unpass_tree_cutoff']
    
    df_DP_results, passtree_mutations, onecell_mutations = run_dp_pass_tree(
        data = data, 
        df_features_new = df_features_new, 
        M_scaffold = M_scaffold, 
        outputpath_full = outputpath_full, 
        scaffold_mutations = scaffold_mutations,
        p_thresh=p_thresh,
        pass_tree_cutoff=pass_tree_cutoff,
        unpass_tree_cutoff=unpass_tree_cutoff,
        is_log_value_for_likelihoods=True
    )
    uppasstree_mutations = [i for i in all_mutations if i not in passtree_mutations]
else:
    logging.info("Skip running dynamic programming ...")
    passtree_mutations = all_mutations



# ##### LR
# logging.info("Running tree-based classification ...")












##### 准备建树用的数据

logging.info("Prepare the data for fuul-resolved tree building ...")

attached_mutations = [i for i in passtree_mutations if i not in scaffold_mutations and i not in removed_germline_mutations and i not in removed_artifact_mutations]
logging.info(f"The number of attached_mutations is: {len(attached_mutations)}")
# 2025-10-15 18:26:45,114 [INFO] The number of attached_mutations is: 53

I_attached_selected = I[scaffold_mutations + attached_mutations]
I_attached_selected_sorted = I_attached_selected[I_attached_selected.apply(lambda col: (col == 1).sum(), axis=0).sort_values(ascending=False).index]
I_attached_sorted_non_empty = I_attached_selected_sorted[I_attached_selected_sorted.eq(1).any(axis=1)]

P_attached_selected = P[scaffold_mutations + attached_mutations]
P_attached_selected_sorted = P_attached_selected[P_attached_selected.apply(lambda col: (col == 1).sum(), axis=0).sort_values(ascending=False).index]
P_attached_sorted_non_empty = P_attached_selected_sorted[P_attached_selected_sorted.eq(1).any(axis=1)]

I_attached_split, P_attached_split = split_spots_by_immune_mutations(spots_to_split, [i for i in immune_mutations if i in I_attached_sorted_non_empty.columns], I_attached_sorted_non_empty, P_attached_sorted_non_empty)
I_attached, sorting_stats_of_I_attached = reorder_columns_by_mutant_stats(
    I_attached_split, 
    df_features_new,
    min_cell_threshold=30,  # ≥30的作为高优先级组
    bin_size=5,             # 30以下每5个一组
    descending=True         # 从大到小排序
)
P_attached = P_attached_split[I_attached.columns]

IRank_mutations = I_attached.columns.tolist()
IRank_mutations_reversed = I_attached.columns[::-1].tolist()




# ------------------------------
# Step 6: full-resolved tree
# ------------------------------
logger.info("===== Step6: Construct full-resolved tree ...")

# T_scaffold
# M_scaffold

##### 往树上放 mutations
logger.info("Calculating penalties and refining placements ...")
T_current = copy.deepcopy(T_scaffold)

new_rows_for_current = I_attached.index.difference(M_scaffold.index)
new_data_for_current = pd.DataFrame(0, index=new_rows_for_current, columns=M_scaffold.columns)
M_current_each_mut = pd.concat([M_scaffold, new_data_for_current])
all_nodes_in_T_scaffold = T_scaffold.all_names_no_root()
M_current = merge_mutations(M_current_each_mut, all_nodes_in_T_scaffold)
M_current.insert(0, 'ROOT', 1)

# 计算罚分时的 NA 权重设置
ω_NA = params['general_weight_NA'] if params['general_weight_NA'] else 0.001
fnfp_ratio = params['fnfp_ratio']
φ = params['phi']

# Refining positions for each mutation
root_mutations = []




# ------------------------------
# Step 6.1: 第一次往 scaffold tree 上放 attached_mutations（包括 non_scaffold_mutations 和 rescued_germline_mutations）
# ------------------------------
logger.info("===== Step6.1: The first time to hang mutations on the scaffold tree ...")

##### 第一次往 scaffold tree 上放 mutations
sorted_attached_mutations = [i for i in I_attached.columns if i in attached_mutations]
external_mutations_of_attached_on_scaffold, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
    sorted_attached_mutations=sorted_attached_mutations,
    T_current=T_current,
    M_current=M_current,
    I_attached=I_attached,
    P_attached=P_attached,
    ω_NA=ω_NA,
    fnfp_ratio=fnfp_ratio,
    φ=φ,
    logger=logger,
    root_mutations=root_mutations  # 可选，如果已有根突变列表
)
logging.info(f"The number of external_mutations_of_attached_on_scaffold is: {len(external_mutations_of_attached_on_scaffold)}")

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)




# -----------------------------
# Step6.2: 处理 external_mutations_of_attached_on_scaffold, 也就是第一次在树上找不到 intersection 的 mut 再挂一次
# -----------------------------
logger.info("===== Step6.2: Handle 'external_mutations_of_attached_on_scaffold' after first hang mutations ...")

##### external_mutations_of_attached_on_scaffold
len(external_mutations_of_attached_on_scaffold)
# 11
print(external_mutations_of_attached_on_scaffold)
# ['chr11_47268814_C_A', 'chr22_22895482_G_A', 'chr2_55866672_T_A', 'chr15_64953013_C_T', 'chr22_34670933_A_T', 'chr9_114311431_A_T', 'chr1_5862839_A_T', 'chr12_6451607_T_C', 'chr14_105707803_C_A', 'chr16_69117490_C_T', 'chr22_50526071_G_A']

##### 主处理流程，把第一次挂树没挂上的 external_mutations 重新挂一遍
sorted_external_mutations_of_attached_on_scaffold = [i for i in I_attached.columns if i in external_mutations_of_attached_on_scaffold]
final_external_mutations_of_attached_on_scaffold, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
    sorted_attached_mutations=sorted_external_mutations_of_attached_on_scaffold,
    T_current=T_current,
    M_current=M_current,
    I_attached=I_attached,
    P_attached=P_attached,
    ω_NA=ω_NA,
    fnfp_ratio=fnfp_ratio,
    φ=φ,
    logger=logger,
    root_mutations=root_mutations  # 可选，如果已有根突变列表
)

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)




# -----------------------------
# Step 6.3 依旧没有在树上找到 intersection 的位点最后再在外部进行组合整体挂树
# -----------------------------
logger.info("===== Step6.3: Despite not finding the intersection point on the tree, we will finally combine and attach the entire structure externally ...")

###### 依旧未处理的再加到 ROOT 的新节点中

logging.info(f"The number of final_external_mutations_of_attached_on_scaffold is: {len(final_external_mutations_of_attached_on_scaffold)}")

external_mutations = final_external_mutations_of_attached_on_scaffold
logging.info(f"The number of external_mutations is: {len(external_mutations)}")


final_external_mutations = []
if len(external_mutations) > 0:
    
    sorted_external_mutations = [i for i in I_attached.columns if i in external_mutations]
    final_external_mutations, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)


##### 最后还没挂上的突变只能上 ROOT 了
logging.info(f"The number of final_external_mutations is: {len(final_external_mutations)}")

remained_mutations = []
if len(final_external_mutations)>0:
    
    subtree_groups = cluster_external_mutations_by_intersection(I_attached, final_external_mutations)
    
    logger.info("Processing remaining external mutations by building subtrees")
    
    remained_mutations, T_current, M_current, root_mutations = process_external_mutations_by_subtree_groups(
        subtree_groups=subtree_groups,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations
    )

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)




# -----------------------------
# Step 6.4 最后再处理一次 remained_mutations
# -----------------------------
logger.info("===== Step6.4: Finally, process the 'remained_mutations' one more time. ...")

logging.info(f"The number of remained_mutations is: {len(remained_mutations)}")

final_remained_mutations = []
if len(remained_mutations) > 0:
    
    sorted_remained_mutations = [i for i in I_attached.columns if i in remained_mutations]
    final_remained_mutations, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_remained_mutations,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

logging.info(f"The number of final_remained_mutations is: {len(final_remained_mutations)}")

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)


print(T_current)
# └─ ROOT
#   └─ chr17_7578893_G_T
#     └─ chr8_80170791_C_T|chr7_20381814_G_T
#       └─ chr4_122860726_G_T
#         └─ chr15_50499489_A_G
#           └─ chr16_497748_G_T
#             └─ chr10_73254076_G_T
#           └─ chr19_16085831_G_T
#       └─ chr21_46323662_C_T
#         └─ chr10_118753098_G_T
#         └─ chr16_15464093_A_G
#     └─ chr1_171541779_G_T
#   └─ chr6_34246136_G_C
#     └─ chr1_45618093_G_T
#       └─ chr2_171687790_A_T
#         └─ chr21_44856223_C_A
#           └─ chr11_9425206_C_A
#             └─ chr17_50662360_G_T
#               └─ chr7_87196234_G_T
#                 └─ chr2_171324260_C_T
#                   └─ chr18_9135639_C_G
#                   └─ chr12_51215297_G_T
#             └─ chr5_128188419_C_A
#   └─ chr22_37967634_G_T
#     └─ chr4_74388653_G_T
#       └─ chr17_43766282_C_A
#         └─ chr5_128171668_G_C
#           └─ chr11_83165726_G_C
#     └─ chr12_102013001_C_A
#   └─ chrX_41232675_T_A
#     └─ chr1_31260040_C_A
#       └─ chr7_132879457_T_A
#         └─ chr10_104003258_C_T
#         └─ chr7_128458208_G_T




# -----------------------------
# Step 6.5 在树中的每一个 subclone 内部计算每一个突变的 fp_ratio，将异常突变重新挂树
# -----------------------------
logger.info("===== Step6.5: FP_1st, process mutations by fpratio_within_subclone ...")

mutation_clones_for_subclone = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

##### 计算 T_current 中每一个 mutations 的 fp_ratio (within subclone)
T_checkpoint_fpratio_within_subclone = copy.deepcopy(T_current)
M_checkpoint_fpratio_within_subclone = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone = M_checkpoint_fpratio_within_subclone.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpratio_within_subclone = M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone, mutations_on_T_current_fpratio_within_subclone)

df_fp_ratio_fpratio_within_subclone, fp_mutations_dict_for_out_subclone_muts_fpratio_within_subclone, fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone = calculate_fp_ratios_within_subclone(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone, I_attached, mutation_clones_for_subclone)
df_fp_ratio_fpratio_within_subclone
#              identifier  ... ratio_fp_mutations_for_in_subclone_muts
# 0   chr10_104003258_C_T  ...                                0.000000
# 1   chr10_118753098_G_T  ...                                0.000000
# 2    chr10_73254076_G_T  ...                                0.000000
# 3    chr11_83165726_G_C  ...                                0.000000
# 4     chr11_9425206_C_A  ...                                0.000000
# 5   chr12_102013001_C_A  ...                                0.000000
# 6    chr15_50499489_A_G  ...                                0.066667
# 7    chr16_15464093_A_G  ...                                0.000000
# 8      chr16_497748_G_T  ...                                0.033333
# 9    chr17_43766282_C_A  ...                                0.000000
# 10   chr17_50662360_G_T  ...                                0.000000
# 11    chr17_7578893_G_T  ...                                0.133333
# 12    chr18_9135639_C_G  ...                                0.000000
# 13    chr1_31260040_C_A  ...                                0.000000
# 14    chr1_45618093_G_T  ...                                0.000000
# 15   chr21_44856223_C_A  ...                                0.000000
# 16   chr21_46323662_C_T  ...                                0.066667
# 17   chr22_37967634_G_T  ...                                0.000000
# 18   chr2_171324260_C_T  ...                                0.000000
# 19   chr2_171687790_A_T  ...                                0.000000
# 20   chr4_122860726_G_T  ...                                0.066667
# 21    chr4_74388653_G_T  ...                                0.000000
# 22   chr5_128171668_G_C  ...                                0.000000
# 23   chr5_128188419_C_A  ...                                0.000000
# 24    chr6_34246136_G_C  ...                                0.000000
# 25   chr7_128458208_G_T  ...                                0.000000
# 26   chr7_132879457_T_A  ...                                0.000000
# 27    chr7_87196234_G_T  ...                                0.000000
# 28    chr8_80170791_C_T  ...                                0.133333
# 29    chrX_41232675_T_A  ...                                0.000000

# [30 rows x 9 columns]


# 观察数据设置 fp_ratio_cutoff_within_subclone
# fp_ratio_cutoff_within_subclone = 0.1
fp_ratio_cutoff_within_subclone = params['fp_ratio_cutoff_within_subclone']

# 筛选导致 fp 高的位点
rehanged_mutations_by_fpratio_within_subclone = df_fp_ratio_fpratio_within_subclone[df_fp_ratio_fpratio_within_subclone['fp_ratio_within_subclone_for_in_subclone_muts'] >= fp_ratio_cutoff_within_subclone]['identifier'].tolist()
rehanged_mutations_by_fpratio_within_subclone_but_backbone = [i for i in rehanged_mutations_by_fpratio_within_subclone if i not in list(set(expanded_mutations_of_current_backbone_nodes+scaffold_mutations))]
print(len(rehanged_mutations_by_fpratio_within_subclone_but_backbone))
# 0
print(rehanged_mutations_by_fpratio_within_subclone_but_backbone)
# []

# 获取 fp_ratio>=0.1 的位点之间在树上的关系
ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(T_current, rehanged_mutations_by_fpratio_within_subclone_but_backbone)
earlist_mutation_of_ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone = list(ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone.keys())

# 找到 fp_ratio>=0.1 的位点导致发生 fp 的 mutations
filtered_fp_mutations_dict_by_fpratio_within_subclone = {mut: other_muts 
    for mut, other_muts in fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone.items() 
    if mut in rehanged_mutations_by_fpratio_within_subclone_but_backbone}

# 找到 fp_ratio>=0.1 的位点的 daughter mutations
nodes_rehanged_mutations_by_fpratio_within_subclone_but_backbone = list(set([find_mutation_column(mutation, M_current.columns) for mutation in rehanged_mutations_by_fpratio_within_subclone_but_backbone]))

node_dict = {node.name: node for node in T_current.traverse()}
daughters_nodes_dict_by_fpratio_within_subclone = {
    mutation: get_all_daughter_mutations(node_dict[mutation])
    for mutation in nodes_rehanged_mutations_by_fpratio_within_subclone_but_backbone
    if mutation in node_dict  # 确保突变存在于树中
}
daughters_mutations_dict_by_fpratio_within_subclone = [i for i in list(set([daughter for daughters_list in daughters_nodes_dict_by_fpratio_within_subclone.values() for daughter in daughters_list])) if i not in rehanged_mutations_by_fpratio_within_subclone_but_backbone]


# 需要重新处理的位点 pool 到一起并且去重分组
current_backbone_nodes_fpratio_within_subclone = get_first_level_backbone_nodes(T_current)
expanded_current_backbone_nodes_fpratio_within_subclone = [mutation for node in current_backbone_nodes_fpratio_within_subclone for mutation in node.split('|')]

daughters_to_leaf_mutations_fpratio_within_subclone = []
fp_mutations_fpratio_within_subclone = []
for idx, branch_mut in enumerate(ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone):
    daughter_list = ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone[branch_mut]
    
    daughters_to_leaf_mutations_fpratio_within_subclone = list({
        item 
        for key in daughter_list
        if key in daughters_mutations_dict_by_fpratio_within_subclone
        for item in [key] + daughters_mutations_dict_by_fpratio_within_subclone[key]
    })
    
    fp_mutations_fpratio_within_subclone_initail = list({
        item 
        for key in daughter_list
        if key in filtered_fp_mutations_dict_by_fpratio_within_subclone
        for item in [key] + filtered_fp_mutations_dict_by_fpratio_within_subclone[key]
    })
    fp_mutations_fpratio_within_subclone = [i for i in fp_mutations_fpratio_within_subclone_initail if i not in daughter_list and i not in daughters_to_leaf_mutations_fpratio_within_subclone and i not in expanded_current_backbone_nodes_fpratio_within_subclone]

# 在 cell of orgin 的 candidate 中加一个大于 5 的条件
daughters_to_leaf_mutations_fpratio_within_subclone_qc = [mut for mut in daughters_to_leaf_mutations_fpratio_within_subclone if df_features_new[mut]['mutant_cellnum'] > 5]
fp_mutations_fpratio_within_subclone_qc = set(fp_mutations_fpratio_within_subclone+[i for i in daughters_to_leaf_mutations_fpratio_within_subclone if i not in daughters_to_leaf_mutations_fpratio_within_subclone_qc])

# sorted by our sorting method
# likely_tumor_cell_of_origin_mutations = ['chr11_59636872_G_A']
# fp_mutations_fpratio_within_subclone_qc = [i for i in fp_mutations_fpratio_within_subclone_qc if i not in likely_tumor_cell_of_origin_mutations]
# daughters_to_leaf_mutations_fpratio_within_subclone_qc = daughters_to_leaf_mutations_fpratio_within_subclone_qc+likely_tumor_cell_of_origin_mutations
sorted_fp_mutations_fpratio_within_subclone = [i for i in I_attached.columns if i in fp_mutations_fpratio_within_subclone_qc]
sorted_daughters_to_leaf_mutations_fpratio_within_subclone = [i for i in I_attached.columns if i in daughters_to_leaf_mutations_fpratio_within_subclone_qc]
sorted_rehanged_mutations_all_fpratio_within_subclone = sorted_fp_mutations_fpratio_within_subclone + sorted_daughters_to_leaf_mutations_fpratio_within_subclone

# reversed_sorted_fp_mutations_fpratio_within_subclone = [i for i in IRank_mutations_reversed if i in fp_mutations_fpratio_within_subclone_qc]
# reversed_sorted_daughters_to_leaf_mutations_fpratio_within_subclone = [i for i in IRank_mutations_reversed if i in daughters_to_leaf_mutations_fpratio_within_subclone_qc]
# reversed_sorted_rehanged_mutations_all_fpratio_within_subclone = sorted_fp_mutations_fpratio_within_subclone + sorted_daughters_to_leaf_mutations_fpratio_within_subclone


##### 按照顺序依次处理 fp_mutations_fpratio_within_subclone 和 daughters_to_leaf_mutations_fpratio_within_subclone, 将他们从树上抽离再重挂（而且按理说它们应该是一定会挂在树上的）
external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone = []
external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone = []
if len(sorted_rehanged_mutations_all_fpratio_within_subclone) > 0:
    
    T_removed_fpratio_within_subclone, M_removed_fpratio_within_subclone = remove_mutations_from_tree_and_matrix(T_checkpoint_fpratio_within_subclone, M_checkpoint_fpratio_within_subclone, sorted_rehanged_mutations_all_fpratio_within_subclone)
    print_tree(T_removed_fpratio_within_subclone)
    logger.info(f"The shape of removed_tree to be refined is : {M_removed_fpratio_within_subclone.shape}")    
    
    T_current = copy.deepcopy(T_removed_fpratio_within_subclone)
    M_current = M_removed_fpratio_within_subclone.copy()

if len(sorted_fp_mutations_fpratio_within_subclone) > 0:
    # 首先重挂 sorted_fp_mutations_fpratio_within_subclone
    external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_fp_mutations_fpratio_within_subclone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

if len(sorted_daughters_to_leaf_mutations_fpratio_within_subclone) > 0:
    # 其次重挂 sorted_daughters_to_leaf_mutations_fpratio_within_subclone, 并且调高 fnfp_ratio
    external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_daughters_to_leaf_mutations_fpratio_within_subclone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)


##### 把按照顺序没有挂树的 external_mutations_fpratio_within_subclone 重新挂一遍

external_mutations_fpratio_within_subclone = list(set(external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone + external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone))

logging.info(f"The number of external_mutations_fpratio_within_subclone is: {len(external_mutations_fpratio_within_subclone)}")

final_external_mutations_fpratio_within_subclone = []
if len(external_mutations_fpratio_within_subclone) > 0:
    
    sorted_external_mutations_fpratio_within_subclone = [i for i in I_attached.columns if i in external_mutations_fpratio_within_subclone]
    final_external_mutations_fpratio_within_subclone, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations_fpratio_within_subclone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

logging.info(f"The number of final_external_mutations_fpratio_within_subclone is: {len(final_external_mutations_fpratio_within_subclone)}")

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)




##### 最后再算一次放进要保存的文件中以便查看
mutation_clones_for_subclone_v2 = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)

##### 计算 T_current 中每一个 mutations 的 fp_ratio (whitin subclone)
T_checkpoint_fpratio_within_subclone_v2 = copy.deepcopy(T_current)
M_checkpoint_fpratio_within_subclone_v2 = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2 = M_checkpoint_fpratio_within_subclone_v2.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpratio_within_subclone_v2 = M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2 = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2, mutations_on_T_current_fpratio_within_subclone_v2)

df_fp_ratio_fpratio_within_subclone_v2, fp_mutations_dict_for_out_subclone_muts_fpratio_within_subclone_v2, fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone_v2 = calculate_fp_ratios_within_subclone(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2, I_attached, mutation_clones_for_subclone_v2)
df_fp_ratio_fpratio_within_subclone_v2
#              identifier  ... ratio_fp_mutations_for_in_subclone_muts
# 0   chr10_104003258_C_T  ...                                0.000000
# 1   chr10_118753098_G_T  ...                                0.000000
# 2    chr10_73254076_G_T  ...                                0.000000
# 3    chr11_83165726_G_C  ...                                0.000000
# 4     chr11_9425206_C_A  ...                                0.000000
# 5   chr12_102013001_C_A  ...                                0.000000
# 6    chr15_50499489_A_G  ...                                0.066667
# 7    chr16_15464093_A_G  ...                                0.000000
# 8      chr16_497748_G_T  ...                                0.033333
# 9    chr17_43766282_C_A  ...                                0.000000
# 10   chr17_50662360_G_T  ...                                0.000000
# 11    chr17_7578893_G_T  ...                                0.133333
# 12    chr18_9135639_C_G  ...                                0.000000
# 13    chr1_31260040_C_A  ...                                0.000000
# 14    chr1_45618093_G_T  ...                                0.000000
# 15   chr21_44856223_C_A  ...                                0.000000
# 16   chr21_46323662_C_T  ...                                0.066667
# 17   chr22_37967634_G_T  ...                                0.000000
# 18   chr2_171324260_C_T  ...                                0.000000
# 19   chr2_171687790_A_T  ...                                0.000000
# 20   chr4_122860726_G_T  ...                                0.066667
# 21    chr4_74388653_G_T  ...                                0.000000
# 22   chr5_128171668_G_C  ...                                0.000000
# 23   chr5_128188419_C_A  ...                                0.000000
# 24    chr6_34246136_G_C  ...                                0.000000
# 25   chr7_128458208_G_T  ...                                0.000000
# 26   chr7_132879457_T_A  ...                                0.000000
# 27    chr7_87196234_G_T  ...                                0.000000
# 28    chr8_80170791_C_T  ...                                0.133333
# 29    chrX_41232675_T_A  ...                                0.000000

# [30 rows x 9 columns]

df_fp_ratio_fpratio_within_subclone_final = pd.merge(
    df_fp_ratio_fpratio_within_subclone, 
    df_fp_ratio_fpratio_within_subclone_v2, 
    on='identifier', 
    suffixes=('.1', '.2')
)




# -----------------------------
# Step 6.6 在整个树上计算 fp_ratio 和 fn_ratio，将异常突变重新挂树
# -----------------------------
logger.info("===== Step6.6: FP_2nd, process mutations by fpfnratio_across_tree ...")

current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

##### 计算 T_current 中每一个 mutations 的 fp_ratio 和 fn_ratio (across tree)
T_checkpoint_fpfnratio_across_tree = copy.deepcopy(T_current)
M_checkpoint_fpfnratio_across_tree = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree = M_checkpoint_fpfnratio_across_tree.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpfnratio_across_tree = M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree, mutations_on_T_current_fpfnratio_across_tree)

df_fp_ratio_and_fn_ratio_fpfnratio_across_tree, fp_mutations_dict_fpfnratio_across_tree = calculate_fp_fn_ratios_across_tree(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree, I_attached)
df_fp_ratio_and_fn_ratio_fpfnratio_across_tree
#              identifier  fp_ratio  fn_ratio
# 0     chr17_7578893_G_T  0.362069  0.760684
# 1     chr8_80170791_C_T  0.278161  0.526316
# 2    chr4_122860726_G_T  0.150575  0.333333
# 3    chr15_50499489_A_G  0.116092  0.611111
# 4      chr16_497748_G_T  0.091954  0.000000
# 5    chr10_73254076_G_T  0.000000  0.000000
# 6    chr21_46323662_C_T  0.093103  0.058824
# 7   chr10_118753098_G_T  0.006897  0.000000
# 8    chr16_15464093_A_G  0.000000  0.000000
# 9     chr6_34246136_G_C  0.293103  0.194915
# 10    chr1_45618093_G_T  0.114943  0.331776
# 11   chr2_171687790_A_T  0.097701  0.493976
# 12   chr21_44856223_C_A  0.097701  0.256410
# 13    chr11_9425206_C_A  0.045977  0.351351
# 14   chr17_50662360_G_T  0.045977  0.282609
# 15    chr7_87196234_G_T  0.022989  0.159091
# 16   chr2_171324260_C_T  0.022989  0.035714
# 17    chr18_9135639_C_G  0.000000  0.000000
# 18   chr22_37967634_G_T  0.034483  0.550000
# 19  chr12_102013001_C_A  0.000000  0.000000
# 20   chr17_43766282_C_A  0.000000  0.636364
# 21    chr4_74388653_G_T  0.000000  0.000000
# 22    chrX_41232675_T_A  0.000000  0.500000
# 23    chr1_31260040_C_A  0.000000  0.200000
# 24   chr7_132879457_T_A  0.000000  0.375000
# 25  chr10_104003258_C_T  0.000000  0.000000
# 26   chr5_128188419_C_A  0.000000  0.000000
# 27   chr5_128171668_G_C  0.000000  0.666667
# 28   chr7_128458208_G_T  0.000000  0.000000
# 29   chr11_83165726_G_C  0.000000  0.000000


# 观察数据设置 fp_ratio_cutoff_across_tree 和 fn_ratio_cutoff_across_tree
# fp_ratio_cutoff_across_tree = 0.1
# fn_ratio_cutoff_across_tree = 0.9
fp_ratio_cutoff_across_tree=params['fp_ratio_cutoff_across_tree']
fn_ratio_cutoff_across_tree=params['fn_ratio_cutoff_across_tree']


##### 筛选导致 fp 高的位点并处理
rehanged_fp_mutations_by_fpfnratio_across_tree = df_fp_ratio_and_fn_ratio_fpfnratio_across_tree[df_fp_ratio_and_fn_ratio_fpfnratio_across_tree['fp_ratio'] >= fp_ratio_cutoff_across_tree]['identifier'].tolist()
rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = [i for i in rehanged_fp_mutations_by_fpfnratio_across_tree if i not in list(set(expanded_mutations_of_current_backbone_nodes+scaffold_mutations))]
print(len(rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone))
# 0
print(rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone)
# []

# 获取 fp_ratio>=0.2 的位点之间在树上的关系
ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(T_current, rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone)
earlist_mutation_of_ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = list(ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone.keys())

# 找到 fp_ratio>=0.2 的位点导致发生 fp 的 mutations
filtered_fp_mutations_dict_by_fpfnratio_across_tree = {mut: other_muts 
    for mut, other_muts in fp_mutations_dict_fpfnratio_across_tree.items() 
    if mut in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone}

# 找到 fp_ratio>=0.2 的位点的 daughter mutations
nodes_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = list(set([find_mutation_column(mutation, M_current.columns) for mutation in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone]))

node_dict = {node.name: node for node in T_current.traverse()}
daughters_nodes_dict_by_fpfnratio_across_tree = {
    mutation: get_all_daughter_mutations(node_dict[mutation])
    for mutation in nodes_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
    if mutation in node_dict  # 确保突变存在于树中
}
daughters_mutations_dict_by_fpfnratio_across_tree = [i for i in list(set([daughter for daughters_list in daughters_nodes_dict_by_fpfnratio_across_tree.values() for daughter in daughters_list])) if i not in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone]


# 需要重新处理的位点 pool 到一起并且去重分组
daughters_to_leaf_mutations_fpfnratio_across_tree = []
fp_mutations_fpfnratio_across_tree = []
for idx, branch_mut in enumerate(ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone):
    daughter_list = ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone[branch_mut]
    
    daughters_to_leaf_mutations_fpfnratio_across_tree = list({
        item 
        for key in daughter_list
        if key in daughters_mutations_dict_by_fpfnratio_across_tree
        for item in [key] + daughters_mutations_dict_by_fpfnratio_across_tree[key]
    })
    
    fp_mutations_fpfnratio_across_tree_initail = list({
        item 
        for key in daughter_list
        if key in filtered_fp_mutations_dict_by_fpfnratio_across_tree
        for item in [key] + filtered_fp_mutations_dict_by_fpfnratio_across_tree[key]
    })
    fp_mutations_fpfnratio_across_tree = [i for i in fp_mutations_fpfnratio_across_tree_initail if i not in daughter_list and i not in daughters_to_leaf_mutations_fpfnratio_across_tree and i not in list(set(expanded_mutations_of_current_backbone_nodes+scaffold_mutations))]

# 在 cell of orgin 的 candidate 中加一个大于 5 的条件
daughters_to_leaf_mutations_fpfnratio_across_tree_qc = [mut for mut in daughters_to_leaf_mutations_fpfnratio_across_tree if df_features_new[mut]['mutant_cellnum'] > 5]
fp_mutations_fpfnratio_across_tree_qc = set(fp_mutations_fpfnratio_across_tree+[i for i in daughters_to_leaf_mutations_fpfnratio_across_tree if i not in daughters_to_leaf_mutations_fpfnratio_across_tree_qc])

# sorted by our sorting method
sorted_fp_mutations_fpfnratio_across_tree = [i for i in I_attached.columns if i in fp_mutations_fpfnratio_across_tree_qc]
sorted_daughters_to_leaf_mutations_fpfnratio_across_tree = [i for i in I_attached.columns if i in daughters_to_leaf_mutations_fpfnratio_across_tree_qc]
sorted_rehanged_mutations_all_fpfnratio_across_tree = sorted_fp_mutations_fpfnratio_across_tree + sorted_daughters_to_leaf_mutations_fpfnratio_across_tree

# reversed_sorted_fp_mutations_fpfnratio_across_tree = [i for i in IRank_mutations_reversed if i in fp_mutations_fpfnratio_across_tree_qc]
# reversed_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree = [i for i in IRank_mutations_reversed if i in daughters_to_leaf_mutations_fpfnratio_across_tree_qc]
# reversed_sorted_rehanged_mutations_all_fpfnratio_across_tree = sorted_fp_mutations_fpfnratio_across_tree + sorted_daughters_to_leaf_mutations_fpfnratio_across_tree


##### 按照顺序依次处理 fp_mutations_fpfnratio_across_tree 和 daughters_to_leaf_mutations_fpfnratio_across_tree, 将他们从树上抽离再重挂（而且按理说它们应该是一定会挂在树上的）
external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree = []
external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree = []
if len(sorted_rehanged_mutations_all_fpfnratio_across_tree) > 0:
    
    T_removed_fpfnratio_across_tree, M_removed_fpfnratio_across_tree = remove_mutations_from_tree_and_matrix(T_checkpoint_fpfnratio_across_tree, M_checkpoint_fpfnratio_across_tree, sorted_rehanged_mutations_all_fpfnratio_across_tree)
    print_tree(T_removed_fpfnratio_across_tree)
    logger.info(f"The shape of removed_tree to be refined is : {M_removed_fpfnratio_across_tree.shape}")    
    
    T_current = copy.deepcopy(T_removed_fpfnratio_across_tree)
    M_current = M_removed_fpfnratio_across_tree.copy()

if len(sorted_fp_mutations_fpfnratio_across_tree) > 0:
    # 首先重挂 sorted_fp_mutations_fpfnratio_across_tree
    external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_fp_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

if len(sorted_daughters_to_leaf_mutations_fpfnratio_across_tree) > 0:
    # 其次重挂 sorted_daughters_to_leaf_mutations_fpfnratio_across_tree, 并且调高 fnfp_ratio
    external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_daughters_to_leaf_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)


##### 筛选导致 fn 高的位点并处理
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

rehanged_fn_mutations_by_fpfnratio_across_tree = df_fp_ratio_and_fn_ratio_fpfnratio_across_tree[df_fp_ratio_and_fn_ratio_fpfnratio_across_tree['fn_ratio'] >= fn_ratio_cutoff_across_tree]['identifier'].tolist()
sorted_fn_mutations_fpfnratio_across_tree = [i for i in I_attached.columns if i in rehanged_fn_mutations_by_fpfnratio_across_tree]
len(sorted_fn_mutations_fpfnratio_across_tree)
# 0
sorted_fn_mutations_fpfnratio_across_tree
# []

external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree = []
if len(rehanged_fn_mutations_by_fpfnratio_across_tree) > 0:
    
    T_removed_fpfnratio_across_tree, M_removed_fpfnratio_across_tree = remove_mutations_from_tree_and_matrix(T_checkpoint_fpfnratio_across_tree, M_checkpoint_fpfnratio_across_tree, rehanged_fn_mutations_by_fpfnratio_across_tree)
    print_tree(T_removed_fpfnratio_across_tree)
    logger.info(f"The shape of removed_tree to be refined is : {M_removed_fpfnratio_across_tree.shape}")    
    
    T_current = copy.deepcopy(T_removed_fpfnratio_across_tree)
    M_current = M_removed_fpfnratio_across_tree.copy()

if len(sorted_fn_mutations_fpfnratio_across_tree) > 0:
    # 首先重挂 sorted_fn_mutations_fpfnratio_across_tree
    external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_fn_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )


##### 把按照顺序没有挂树的 external_mutations_fpfnratio_across_tree 重新挂一遍
external_mutations_fpfnratio_across_tree = list(set(external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree+external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree+external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree))
len(external_mutations_fpfnratio_across_tree)
# 0
external_mutations_fpfnratio_across_tree
# []

logging.info(f"The number of external_mutations_fpfnratio_across_tree is: {len(external_mutations_fpfnratio_across_tree)}")

final_external_mutations_fpfnratio_across_tree = []
if len(external_mutations_fpfnratio_across_tree) > 0:
    
    sorted_external_mutations_fpfnratio_across_tree = [i for i in I_attached.columns if i in external_mutations_fpfnratio_across_tree]
    final_external_mutations_fpfnratio_across_tree, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations_fpfnratio_across_tree,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

logging.info(f"The number of final_external_mutations_fpfnratio_across_tree is: {len(final_external_mutations_fpfnratio_across_tree)}")

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1184, 34)




##### 计算 T_current 中每一个 mutations 的 fp_ratio 和 fn_ratio (across tree)
T_checkpoint_fpfnratio_across_tree_v2 = copy.deepcopy(T_current)
M_checkpoint_fpfnratio_across_tree_v2 = M_current.copy()

M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2 = M_checkpoint_fpfnratio_across_tree_v2.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_fpfnratio_across_tree_v2 = M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2 = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, mutations_on_T_current_fpfnratio_across_tree_v2)

df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, fp_mutations_dict_fpfnratio_across_tree_v2 = calculate_fp_fn_ratios_across_tree(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, I_attached)
df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2
#              identifier  fp_ratio  fn_ratio
# 0     chr17_7578893_G_T  0.362069  0.760684
# 1     chr8_80170791_C_T  0.278161  0.526316
# 2    chr4_122860726_G_T  0.150575  0.333333
# 3    chr15_50499489_A_G  0.116092  0.611111
# 4      chr16_497748_G_T  0.091954  0.000000
# 5    chr10_73254076_G_T  0.000000  0.000000
# 6    chr21_46323662_C_T  0.093103  0.058824
# 7   chr10_118753098_G_T  0.006897  0.000000
# 8    chr16_15464093_A_G  0.000000  0.000000
# 9     chr6_34246136_G_C  0.293103  0.194915
# 10    chr1_45618093_G_T  0.114943  0.331776
# 11   chr2_171687790_A_T  0.097701  0.493976
# 12   chr21_44856223_C_A  0.097701  0.256410
# 13    chr11_9425206_C_A  0.045977  0.351351
# 14   chr17_50662360_G_T  0.045977  0.282609
# 15    chr7_87196234_G_T  0.022989  0.159091
# 16   chr2_171324260_C_T  0.022989  0.035714
# 17    chr18_9135639_C_G  0.000000  0.000000
# 18   chr22_37967634_G_T  0.034483  0.550000
# 19  chr12_102013001_C_A  0.000000  0.000000
# 20   chr17_43766282_C_A  0.000000  0.636364
# 21    chr4_74388653_G_T  0.000000  0.000000
# 22    chrX_41232675_T_A  0.000000  0.500000
# 23    chr1_31260040_C_A  0.000000  0.200000
# 24   chr7_132879457_T_A  0.000000  0.375000
# 25  chr10_104003258_C_T  0.000000  0.000000
# 26   chr5_128188419_C_A  0.000000  0.000000
# 27   chr5_128171668_G_C  0.000000  0.666667
# 28   chr7_128458208_G_T  0.000000  0.000000
# 29   chr11_83165726_G_C  0.000000  0.000000

df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_final = pd.merge(
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree, 
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, 
    on='identifier', 
    suffixes=('.1', '.2')
)


# 4 次全部合并并且保存
combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree = pd.merge(df_fp_ratio_fpratio_within_subclone_final,
                     df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_final,
                     on='identifier',  # 依据 'identifier' 这一列合并
                     how='left')       # 以 df_fp_ratio_fpratio_within_subclone_v2 为基准，缺失的列填充为 NaN

# combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree.to_csv(os.path.join(outputpath_full, "combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree.csv"), sep=",")




# -----------------------------
# Step 6.7 在树中的每一个 subclone 内部计算每一个突变自身的 fp_ratio，将异常突变重新挂树
# -----------------------------
logger.info("===== Step6.7: FP_3rd, process mutations by fp_ratio_persitefp ...")

mutation_clones_for_persitefp = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)
current_backbone_nodes = get_first_level_backbone_nodes(T_current)
expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]

##### 计算 T_current 中每一个 subclone 内部计算每一个突变自身的 fp_ratio
T_checkpoint_fp_ratio_persitefp = copy.deepcopy(T_current)
M_checkpoint_fp_ratio_persitefp = M_current.copy()

M_for_fp_ratio_persitefp = M_checkpoint_fp_ratio_persitefp.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_persitefp = M_for_fp_ratio_persitefp.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_persitefp = split_merged_columns(M_for_fp_ratio_persitefp, mutations_on_T_current_persitefp)

df_fp_ratio_persitefp = calculate_fp_ratios_persite_within_subclone(M_for_fp_ratio_persitefp, I_attached, mutation_clones_for_persitefp)
df_fp_ratio_persitefp
#              identifier subclone_representative  ...  total_ones_in_subclone  subclone_size
# 0     chr17_7578893_G_T       chr17_7578893_G_T  ...                     241              9
# 1     chr8_80170791_C_T       chr17_7578893_G_T  ...                      72              9
# 2    chr4_122860726_G_T       chr17_7578893_G_T  ...                      48              9
# 3    chr15_50499489_A_G       chr17_7578893_G_T  ...                      34              9
# 4      chr16_497748_G_T       chr17_7578893_G_T  ...                      33              9
# 5    chr10_73254076_G_T       chr17_7578893_G_T  ...                      11              9
# 6    chr21_46323662_C_T       chr17_7578893_G_T  ...                      34              9
# 7   chr10_118753098_G_T       chr17_7578893_G_T  ...                      17              9
# 8    chr16_15464093_A_G       chr17_7578893_G_T  ...                      15              9
# 9     chr6_34246136_G_C       chr6_34246136_G_C  ...                     421             10
# 10    chr1_45618093_G_T       chr6_34246136_G_C  ...                     152             10
# 11   chr2_171687790_A_T       chr6_34246136_G_C  ...                     118             10
# 12   chr21_44856223_C_A       chr6_34246136_G_C  ...                      60             10
# 13    chr11_9425206_C_A       chr6_34246136_G_C  ...                      50             10
# 14   chr17_50662360_G_T       chr6_34246136_G_C  ...                      37             10
# 15    chr7_87196234_G_T       chr6_34246136_G_C  ...                      37             10
# 16   chr2_171324260_C_T       chr6_34246136_G_C  ...                      36             10
# 17    chr18_9135639_C_G       chr6_34246136_G_C  ...                      18             10
# 18   chr5_128188419_C_A       chr6_34246136_G_C  ...                      14             10
# 19   chr22_37967634_G_T      chr22_37967634_G_T  ...                     129              6
# 20  chr12_102013001_C_A      chr22_37967634_G_T  ...                      13              6
# 21   chr17_43766282_C_A      chr22_37967634_G_T  ...                      13              6
# 22    chr4_74388653_G_T      chr22_37967634_G_T  ...                      11              6
# 23   chr5_128171668_G_C      chr22_37967634_G_T  ...                       8              6
# 24   chr11_83165726_G_C      chr22_37967634_G_T  ...                       7              6
# 25    chrX_41232675_T_A       chrX_41232675_T_A  ...                      34              5
# 26    chr1_31260040_C_A       chrX_41232675_T_A  ...                      22              5
# 27   chr7_132879457_T_A       chrX_41232675_T_A  ...                      15              5
# 28  chr10_104003258_C_T       chrX_41232675_T_A  ...                      10              5
# 29   chr7_128458208_G_T       chrX_41232675_T_A  ...                       9              5

# [30 rows x 6 columns]

# df_fp_ratio_persitefp.to_csv(os.path.join(outputpath_full, "df_fp_ratio_persitefp.csv"), sep=",")


# 观察数据设置 fp_ratio_cutoff_within_subclone
# fp_ratio_cutoff_within_subclone = 0.1
fp_count_persite_cutoff = params['fp_count_persite_cutoff']
fp_ratio_persite_cutoff = params['fp_ratio_persite_cutoff']

##### 筛选自身 fp 高的位点并处理
rehanged_mutations_by_persitefp = df_fp_ratio_persitefp[df_fp_ratio_persitefp['fp_ratio_persite'] >= fp_ratio_persite_cutoff]['identifier'].tolist()
rehanged_mutations_by_persitefp_but_backbone = [i for i in rehanged_mutations_by_persitefp if i not in list(set(expanded_mutations_of_current_backbone_nodes+scaffold_mutations))]
print(len(rehanged_mutations_by_persitefp_but_backbone))
# 0
print(rehanged_mutations_by_persitefp_but_backbone)
# []
for m in rehanged_mutations_by_persitefp_but_backbone:
    print(m)
    print(count_list(I_attached[m]))
    print(count_conditions(I_attached[m], M_current[[col for col in M_current.columns if m in col][0]]))

# chr15_50499489_A_G
# {0.0: 252, 'NA': 839, 1.0: 34}
# {'count_fp_1_0': 33, 'count_fn_0_1': 2, 'count_na_1': 11, 'count_na_0': 828}


sorted_rehanged_mutations_by_persitefp_but_backbone = [i for i in I_attached.columns if i in rehanged_mutations_by_persitefp_but_backbone]

external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone = []
if len(sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
    
    T_removed_fp_ratio_persitefp, M_removed_fp_ratio_persitefp = remove_mutations_from_tree_and_matrix(T_checkpoint_fp_ratio_persitefp, M_checkpoint_fp_ratio_persitefp, sorted_rehanged_mutations_by_persitefp_but_backbone)
    print_tree(T_removed_fp_ratio_persitefp)
    logger.info(f"The shape of removed_tree to be refined is : {M_removed_fp_ratio_persitefp.shape}")    
    
    T_current = copy.deepcopy(T_removed_fp_ratio_persitefp)
    M_current = M_removed_fp_ratio_persitefp.copy()

if len(sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
    # 首先重挂 sorted_rehanged_mutations_by_persitefp_but_backbone
    external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_rehanged_mutations_by_persitefp_but_backbone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )


##### 把按照顺序没有挂树的 external_mutations_fp_ratio_persitefp 重新挂一遍

final_external_mutations_fp_ratio_persitefp = []
if len(external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
    
    sorted_external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone = [i for i in I_attached.columns if i in external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone]
    final_external_mutations_fp_ratio_persitefp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

logging.info(f"The number of final_external_mutations_fp_ratio_persitefp is: {len(final_external_mutations_fp_ratio_persitefp)}")

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1178, 29)

for m in rehanged_mutations_by_persitefp_but_backbone:
    print(m)
    print(count_list(I_attached[m]))
    print(count_conditions(I_attached[m], M_current[[col for col in M_current.columns if m in col][0]]))




##### 再次计算 T_current 中每一个 subclone 内部计算每一个突变自身的 fp_ratio

mutation_clones_for_persitefp_v2 = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)

##### 计算 T_current 中每一个 subclone 内部计算每一个突变自身的 fp_ratio
T_checkpoint_fp_ratio_persitefp_v2 = copy.deepcopy(T_current)
M_checkpoint_fp_ratio_persitefp_v2 = M_current.copy()

M_for_fp_ratio_persitefp_v2 = M_checkpoint_fp_ratio_persitefp_v2.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_persitefp = M_for_fp_ratio_persitefp_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_persitefp_v2 = split_merged_columns(M_for_fp_ratio_persitefp_v2, mutations_on_T_current_persitefp)

df_fp_ratio_persitefp_v2 = calculate_fp_ratios_persite_within_subclone(M_for_fp_ratio_persitefp_v2, I_attached, mutation_clones_for_persitefp_v2)
df_fp_ratio_persitefp_v2
#              identifier subclone_representative  ...  total_ones_in_subclone  subclone_size
# 0     chr17_7578893_G_T       chr17_7578893_G_T  ...                     241              9
# 1     chr8_80170791_C_T       chr17_7578893_G_T  ...                      72              9
# 2    chr4_122860726_G_T       chr17_7578893_G_T  ...                      48              9
# 3    chr15_50499489_A_G       chr17_7578893_G_T  ...                      34              9
# 4      chr16_497748_G_T       chr17_7578893_G_T  ...                      33              9
# 5    chr10_73254076_G_T       chr17_7578893_G_T  ...                      11              9
# 6    chr21_46323662_C_T       chr17_7578893_G_T  ...                      34              9
# 7   chr10_118753098_G_T       chr17_7578893_G_T  ...                      17              9
# 8    chr16_15464093_A_G       chr17_7578893_G_T  ...                      15              9
# 9     chr6_34246136_G_C       chr6_34246136_G_C  ...                     421             10
# 10    chr1_45618093_G_T       chr6_34246136_G_C  ...                     152             10
# 11   chr2_171687790_A_T       chr6_34246136_G_C  ...                     118             10
# 12   chr21_44856223_C_A       chr6_34246136_G_C  ...                      60             10
# 13    chr11_9425206_C_A       chr6_34246136_G_C  ...                      50             10
# 14   chr17_50662360_G_T       chr6_34246136_G_C  ...                      37             10
# 15    chr7_87196234_G_T       chr6_34246136_G_C  ...                      37             10
# 16   chr2_171324260_C_T       chr6_34246136_G_C  ...                      36             10
# 17    chr18_9135639_C_G       chr6_34246136_G_C  ...                      18             10
# 18   chr5_128188419_C_A       chr6_34246136_G_C  ...                      14             10
# 19   chr22_37967634_G_T      chr22_37967634_G_T  ...                     129              6
# 20  chr12_102013001_C_A      chr22_37967634_G_T  ...                      13              6
# 21   chr17_43766282_C_A      chr22_37967634_G_T  ...                      13              6
# 22    chr4_74388653_G_T      chr22_37967634_G_T  ...                      11              6
# 23   chr5_128171668_G_C      chr22_37967634_G_T  ...                       8              6
# 24   chr11_83165726_G_C      chr22_37967634_G_T  ...                       7              6
# 25    chrX_41232675_T_A       chrX_41232675_T_A  ...                      34              5
# 26    chr1_31260040_C_A       chrX_41232675_T_A  ...                      22              5
# 27   chr7_132879457_T_A       chrX_41232675_T_A  ...                      15              5
# 28  chr10_104003258_C_T       chrX_41232675_T_A  ...                      10              5
# 29   chr7_128458208_G_T       chrX_41232675_T_A  ...                       9              5

# [30 rows x 6 columns]


df_fp_ratio_persitefp_final = pd.merge(
    df_fp_ratio_persitefp, 
    df_fp_ratio_persitefp_v2, 
    on='identifier', 
    suffixes=('.1', '.2')
)


# 6 次全部合并并且保存
final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio = pd.merge(combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree,
                     df_fp_ratio_persitefp_final,
                     on='identifier',  # 依据 'identifier' 这一列合并
                     how='left')       # 以 df_fp_ratio_fpratio_within_subclone_v2 为基准，缺失的列填充为 NaN

final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio
#              identifier  ... subclone_size.2_y
# 0   chr10_104003258_C_T  ...                 5
# 1   chr10_118753098_G_T  ...                 9
# 2    chr10_73254076_G_T  ...                 9
# 3    chr11_83165726_G_C  ...                 6
# 4     chr11_9425206_C_A  ...                10
# 5   chr12_102013001_C_A  ...                 6
# 6    chr15_50499489_A_G  ...                 9
# 7    chr16_15464093_A_G  ...                 9
# 8      chr16_497748_G_T  ...                 9
# 9    chr17_43766282_C_A  ...                 6
# 10   chr17_50662360_G_T  ...                10
# 11    chr17_7578893_G_T  ...                 9
# 12    chr18_9135639_C_G  ...                10
# 13    chr1_31260040_C_A  ...                 5
# 14    chr1_45618093_G_T  ...                10
# 15   chr21_44856223_C_A  ...                10
# 16   chr21_46323662_C_T  ...                 9
# 17   chr22_37967634_G_T  ...                 6
# 18   chr2_171324260_C_T  ...                10
# 19   chr2_171687790_A_T  ...                10
# 20   chr4_122860726_G_T  ...                 9
# 21    chr4_74388653_G_T  ...                 6
# 22   chr5_128171668_G_C  ...                 6
# 23   chr5_128188419_C_A  ...                10
# 24    chr6_34246136_G_C  ...                10
# 25   chr7_128458208_G_T  ...                 5
# 26   chr7_132879457_T_A  ...                 5
# 27    chr7_87196234_G_T  ...                10
# 28    chr8_80170791_C_T  ...                 9
# 29    chrX_41232675_T_A  ...                 5

# [30 rows x 31 columns]


final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio.to_csv(os.path.join(outputpath_full, "final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio.csv"), sep=",")




# -----------------------------
# Step 6.8 计算 parent muts 中的 intersection/FN_flip per mutations，最后找到一些要重挂的突变
# -----------------------------
logger.info("===== Step6.8: Calculate the intersection/FN_flip per mutation in the parent muts, and finally identify some mutations that need to be reattached ...")

T_checkpoint_outgroup = copy.deepcopy(T_current)
M_checkpoint_outgroup = M_current.copy()

M_for_fp_ratio_and_fn_ratio_outgroup = M_checkpoint_outgroup.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_outgroup = M_for_fp_ratio_and_fn_ratio_outgroup.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_outgroup = split_merged_columns(M_for_fp_ratio_and_fn_ratio_outgroup, mutations_on_T_current_outgroup)


df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation = calculate_intersection_and_inter_vs_fn_flipping_ratio_per_mutation(T_checkpoint_outgroup, M_checkpoint_outgroup, I_attached)
df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.shape
# (34, 10)

df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.to_csv(os.path.join(outputpath_full, "df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.csv"), sep=",")

##### 找到要重放的的突变，他们的子突变直接在树上重挂，他们自己本身用额外挂到树的 ROOT 上的方法重挂
# intersection_vs_fn_flipping_ratio_cutoff = 0.4
# intersection_cell_count_on_mutation_cutoff = 5
# intersection_cell_ratio_on_mutation_cutoff = 0.2
intersection_vs_fn_flipping_ratio_cutoff=params['intersection_vs_fn_flipping_ratio_cutoff']
intersection_cell_count_on_mutation_cutoff=params['intersection_cell_count_on_mutation_cutoff']
intersection_cell_ratio_on_mutation_cutoff=params['intersection_cell_ratio_on_mutation_cutoff']

outgroup_mutations = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation[(
    (df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['parent_retention_ratio'] <= intersection_vs_fn_flipping_ratio_cutoff) & 
    ((df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['intersection_cell_count_on_mutation'] <= intersection_cell_count_on_mutation_cutoff) | 
        (df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['intersection_cell_ratio_on_mutation'] <= intersection_cell_ratio_on_mutation_cutoff)))]['mutation'].tolist()

outgroup_mutations_but_backbone = [i for i in outgroup_mutations if i not in list(set(expanded_mutations_of_current_backbone_nodes))]
print(len(outgroup_mutations_but_backbone))
# 9
print(outgroup_mutations_but_backbone)
# ['chr8_80170791_C_T', 'chr4_122860726_G_T', 'chr15_50499489_A_G', 'chr17_43766282_C_A', 'chr5_128188419_C_A', 'chr7_20381814_G_T', 'chr5_128171668_G_C', 'chr11_83165726_G_C', 'chr1_171541779_G_T']

if outgroup_mutations_but_backbone:
    # 找到这些突变的子突变们
    nodes_outgroup_mutations_but_backbone = list(set([find_mutation_column(mutation, M_checkpoint_outgroup.columns) for mutation in outgroup_mutations_but_backbone]))
    
    node_dict = {node.name: node for node in T_checkpoint_outgroup.traverse()}
    daughter_nodes_of_outgroup_mutations_but_backbone = {
        mutation: get_all_daughter_mutations(node_dict[mutation])
        for mutation in nodes_outgroup_mutations_but_backbone
        if mutation in node_dict  # 确保突变存在于树中
    }
    daughter_mutations_of_outgroup_mutations_but_backbone = [i for i in list(set([daughter for daughters_list in daughter_nodes_of_outgroup_mutations_but_backbone.values() for daughter in daughters_list])) if i not in outgroup_mutations_but_backbone]
    
    sorted_outgroup_mutations_but_backbone = [i for i in I_attached.columns if i in outgroup_mutations_but_backbone]
    sorted_daughter_mutations_of_outgroup_mutations_but_backbone = [i for i in I_attached.columns if i in daughter_mutations_of_outgroup_mutations_but_backbone]
    sorted_rehanged_mutations_all_outgroup = sorted_outgroup_mutations_but_backbone + sorted_daughter_mutations_of_outgroup_mutations_but_backbone
    
    ##### 删掉以上的要放到 root 上的突变和受到影响的子突变
    if len(sorted_rehanged_mutations_all_outgroup) > 0:
        
        T_removed_outgroup, M_removed_outgroup = remove_mutations_from_tree_and_matrix(T_checkpoint_outgroup, M_checkpoint_outgroup, sorted_rehanged_mutations_all_outgroup)
        print_tree(T_removed_outgroup)
        M_removed_outgroup_modified = process_matrices_by_removed_some_mutations_from_tree(M_removed_outgroup, I_attached)[1]
        
        logger.info(f"The shape of removed_tree to be refined is : {M_removed_outgroup.shape}")    
        
        T_current = copy.deepcopy(T_removed_outgroup)
        M_current = M_removed_outgroup.copy()
    
    ##### 先挂应该位于主树上的那些
    external_mutations_by_sorted_daughter_mutations_of_outgroup_mutations_but_backbone = []
    if len(sorted_daughter_mutations_of_outgroup_mutations_but_backbone) > 0:
        # 首先重挂 sorted_daughter_mutations_of_outgroup_mutations_but_backbone
        external_mutations_by_sorted_daughter_mutations_of_outgroup_mutations_but_backbone, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_daughter_mutations_of_outgroup_mutations_but_backbone,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=ω_NA,
            fnfp_ratio=fnfp_ratio,
            φ=φ,
            logger=logger,
            root_mutations=root_mutations  # 可选，如果已有根突变列表
        )
    
    ##### 再挂要放到主树之外直接挂到 root 上的突变（就是那些没有找到 backbone mutation 的那些 clone 下的突变）
    # 先解决那些 conflict 的 cells 归属问题
    M_current_noROOT = M_current.drop(columns=['ROOT'], errors='ignore')
    M_current_split_and_noROOT = split_merged_columns(M_current_noROOT, M_current_noROOT.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist())
    
    I_attached_only_outgroup = I_attached[sorted_outgroup_mutations_but_backbone].copy()
    
    # 使用函数处理数据处理 conflict 的 cells
    M_current_split_and_noROOT_processed, I_attached_only_outgroup_processed = process_conflicting_cells_allto_outgroup(
        M_current_split_and_noROOT, 
        I_attached_only_outgroup
    )
    
    # 处理列名
    M_current_merged_and_noROOT = merge_duplicate_columns(M_current_split_and_noROOT_processed)
    source_cols = M_current_merged_and_noROOT.columns
    target_cols = M_current_noROOT.columns
    
    column_mapping = create_column_mapping(source_cols, target_cols)
    
    M_current_merged_and_noROOT_renamed = M_current_merged_and_noROOT.rename(columns=column_mapping)
    
    M_current = M_current_merged_and_noROOT_renamed.copy()
    M_current.insert(0, 'ROOT', 1)
    
    I_attached_for_group = I_attached.copy()
    I_attached_for_group.update(I_attached_only_outgroup_processed)
    
    # 重新挂树
    external_mutations_by_sorted_outgroup_mutations_but_backbone = []
    if len(sorted_outgroup_mutations_but_backbone)>0:
        
        subtree_groups = cluster_external_mutations_by_intersection(I_attached, sorted_outgroup_mutations_but_backbone)
        
        logger.info("Processing remaining external mutations by building subtrees")
        
        external_mutations_by_sorted_outgroup_mutations_but_backbone, T_current, M_current, root_mutations = process_external_mutations_by_subtree_groups(
            subtree_groups=subtree_groups,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached_for_group,
            P_attached=P_attached,
            ω_NA=ω_NA,
            fnfp_ratio=fnfp_ratio,
            φ=φ,
            logger=logger,
            root_mutations=root_mutations
        )
    
T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape




# -----------------------------
# Step 6.9 计算 parent muts 中的 intersection/FN_flip per mutations，最后找到一些找不到合适的突变的可以删除的 cells
# -----------------------------
logger.info("===== Step6.9: Calculate the intersection/FN_flip per mutation in the parent muts, and finally identify some cells that cannot find suitable mutations and can be deleted ...")

T_checkpoint_wireless_cells = copy.deepcopy(T_current)
M_checkpoint_wireless_cells = M_current.copy()

M_for_fp_ratio_and_fn_ratio_wireless_cells = M_checkpoint_wireless_cells.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_wireless_cells = M_for_fp_ratio_and_fn_ratio_wireless_cells.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_fp_ratio_and_fn_ratio_wireless_cells = split_merged_columns(M_for_fp_ratio_and_fn_ratio_wireless_cells, mutations_on_T_current_wireless_cells)


df_intersection_and_flipping_to_1_count_per_cell = calculate_intersection_and_flipping_to_1_count_per_cell(M_for_fp_ratio_and_fn_ratio_wireless_cells, I_attached)
df_intersection_and_flipping_to_1_count_per_cell
#                          intersection_count  ...  flipping_to_1_count
# cell                                         ...                     
# 100k_CTACGTCTCCACGCAG-1                   1  ...                    2
# 100k_ACGGCCAGTTGAGTTC-1                   1  ...                    2
# 100k_ATTGGTGCATTATCTC-1                   1  ...                    2
# 100k_GTAACTGGTGGGTATG-1                   1  ...                    2
# 100k_GTTTCTATCTTGGGTA-1                   3  ...                    1
# ...                                     ...  ...                  ...
# 100k_CATGGCGAGACTGTAA-1                   1  ...                    0
# 100k_CGGAGTCAGAACAACT-1                   1  ...                    0
# 100k_CAGATCAAGGTAAACT-1                   1  ...                    0
# 100k_CTCTGGTGTTAGATGA-1                   1  ...                    0
# 100k_CATGGCGGTAGATTAG-1                   1  ...                    0

# [1184 rows x 4 columns]

df_intersection_and_flipping_to_1_count_per_cell.to_csv(os.path.join(outputpath_full, "df_intersection_and_flipping_to_1_count_per_cell.csv"), sep=",")

##### 找到那些 intersection 少且 FN+NA->0 比较多的 cells，这种应该就是缺少正确突变的 cells
# intersection_count_per_cells_cutoff = 1
# flipping_count_fn_per_cells_cutoff = 1
# flipping_to_1_count_per_cells_cutoff = 2
intersection_count_per_cells_cutoff=params['intersection_count_per_cells_cutoff']
flipping_count_fn_per_cells_cutoff=params['flipping_count_fn_per_cells_cutoff']
flipping_to_1_count_per_cells_cutoff=params['flipping_to_1_count_per_cells_cutoff']

df_wireless_cells_filter = df_intersection_and_flipping_to_1_count_per_cell.loc[((df_intersection_and_flipping_to_1_count_per_cell['intersection_count']==intersection_count_per_cells_cutoff) & (df_intersection_and_flipping_to_1_count_per_cell['flipping_count_fn']>=flipping_count_fn_per_cells_cutoff) & (df_intersection_and_flipping_to_1_count_per_cell['flipping_to_1_count']>=flipping_to_1_count_per_cells_cutoff))]

identified_wireless_cells = list(df_wireless_cells_filter.index)

##### 再找到那些重挂到 root 上的突变们中的 doublet cells
conflicting_cells_as_doublets_from_parents_format_nested_list = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.loc[df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['mutation'].isin(outgroup_mutations), 'intersection_cells_on_mutation_parents']
conflicting_cells_as_doublets_from_parents_format_flat_list = sum(conflicting_cells_as_doublets_from_parents_format_nested_list.tolist(), [])

conflicting_cells_as_doublets_from_children_format_nested_list = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.loc[df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['mutation'].isin(outgroup_mutations), 'intersection_cells_on_mutation_children']
conflicting_cells_as_doublets_from_children_format_flat_list = sum(conflicting_cells_as_doublets_from_children_format_nested_list.tolist(), [])


##### 删掉有问题的 cells，修改 M_current 即可
to_be_removed_cells = list(set(identified_wireless_cells + conflicting_cells_as_doublets_from_parents_format_flat_list + conflicting_cells_as_doublets_from_children_format_flat_list))
logger.info(f"The number of identified doublet cells is : {len(to_be_removed_cells)}")

with open(os.path.join(outputpath_full, "likely_doublet_cells_removed_from_tree_by_fn_flipping.csv"), 'w') as f:
    for cell in to_be_removed_cells:
        f.write(cell + '\n')

M_current = M_current.drop(to_be_removed_cells, errors='ignore')

print(M_current.shape)
# (1021, 34)





# -----------------------------
# Step 6.10 计算 cross all cells 的 fp_ratio_per_mutation 和 cross all mutations 的 fp_ratio_per_cell，鉴定并去掉 artifact mutations 和 doublet cells
# -----------------------------
logger.info("===== Step6.10: Caculate fp_ratio_per_mutation_cross_all_cells and fp_ratio_per_cell_cross_all_muts ...")

##### 计算 T_current 中每一个 mutations 的 fp_ratio_per_mutation_cross_all_cells 和 fp_ratio_per_cell_cross_all_muts
T_checkpoint_artifact_and_doublet = copy.deepcopy(T_current)
M_checkpoint_artifact_and_doublet = M_current.copy()

M_for_artifact_and_doublet = M_checkpoint_artifact_and_doublet.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_current_artifact_and_doublet = M_for_artifact_and_doublet.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_for_artifact_and_doublet = split_merged_columns(M_for_artifact_and_doublet, mutations_on_T_current_artifact_and_doublet)

# 计算
df_fp_ratio_per_mutation_cross_all_cells, df_fp_ratio_per_cell_cross_all_muts, overall_metrics, fp_mutations_dict_cross_all_cells = calculate_comprehensive_fp_metrics(
    M_for_artifact_and_doublet, 
    I_attached.loc[M_checkpoint_artifact_and_doublet.index, M_for_artifact_and_doublet.columns]  # 同时筛选行和列
)

df_fp_ratio_per_mutation_cross_all_cells.to_csv(os.path.join(outputpath_full, "df_fp_ratio_per_mutation_cross_all_cells.csv"), sep=",")
df_fp_ratio_per_cell_cross_all_muts.to_csv(os.path.join(outputpath_full, "df_fp_ratio_per_cell_cross_all_muts.csv"), sep=",")

df_fp_ratio_per_mutation_cross_all_cells
#              identifier  fp_cells_count  ...  total_zeros_cells  total_coverage_cells
# 0     chr6_34246136_G_C               5  ...                614                  1040
# 1     chr17_7578893_G_T               4  ...                898                  1143
# 2     chr1_45618093_G_T               1  ...                395                   548
# 3    chr22_37967634_G_T               6  ...                925                  1060
# 4    chr2_171687790_A_T               2  ...                461                   581
# 5     chr8_80170791_C_T               2  ...                275                   349
# 6    chr21_44856223_C_A               0  ...                214                   274
# 7     chr11_9425206_C_A               1  ...                217                   268
# 8    chr4_122860726_G_T               1  ...                164                   213
# 9    chr17_50662360_G_T               1  ...                207                   245
# 10    chr7_87196234_G_T               0  ...                246                   283
# 11   chr2_171324260_C_T               0  ...                111                   147
# 12   chr21_46323662_C_T               1  ...                 57                    92
# 13   chr15_50499489_A_G               1  ...                249                   283
# 14    chrX_41232675_T_A               0  ...                342                   376
# 15     chr16_497748_G_T               1  ...                 73                   106
# 16    chr1_31260040_C_A               0  ...                311                   333
# 17  chr10_118753098_G_T               2  ...                174                   192
# 18    chr18_9135639_C_G               0  ...                 81                    99
# 19   chr16_15464093_A_G               2  ...                 52                    68
# 20  chr12_102013001_C_A               2  ...                329                   344
# 21   chr17_43766282_C_A               2  ...                467                   482
# 22   chr7_132879457_T_A               0  ...                272                   287
# 23   chr5_128188419_C_A               0  ...                330                   344
# 24    chr4_74388653_G_T               1  ...                157                   169
# 25  chr10_104003258_C_T               1  ...                224                   235
# 26   chr10_73254076_G_T               0  ...                 94                   105
# 27   chr5_128171668_G_C               2  ...                423                   433
# 28   chr7_128458208_G_T               0  ...                176                   185
# 29   chr11_83165726_G_C               1  ...                421                   429

# [30 rows x 6 columns]

df_fp_ratio_per_cell_cross_all_muts
#                       cell_id  fp_muts_count  ...  total_zeros_muts  total_coverage_muts
# 0     100k_TTTACTGTCTGATTCT-1              0  ...                 5                    6
# 1     100k_GCTCTGTCAGGTGCCT-1              0  ...                 9                   10
# 2     100k_CTCGAGGAGCGGCTTC-1              0  ...                 9                   10
# 3     100k_GATCTAGCACATGACT-1              0  ...                 7                    8
# 4     100k_CAGCTAATCAAACCAC-1              0  ...                 3                    4
# ...                       ...            ...  ...               ...                  ...
# 1162  100k_TGACGGCTCCCACTTG-1              0  ...                 9                   10
# 1163  100k_ACTTGTTAGACTACAA-1              0  ...                 9                   10
# 1164  100k_CAGCTGGAGACTCGGA-1              0  ...                 8                    9
# 1165  100k_TGAGAGGTCTTACCGC-1              0  ...                11                   12
# 1166  100k_GAGGTGAAGAACTCGG-1              0  ...                 7                   10

# [1167 rows x 6 columns]


# 观察数据设置 fp_ratio_per_mutation_cross_all_cells_cutoff 和 fp_ratio_per_cell_cross_all_muts_cutoff
# fp_ratio_per_mutation_cross_all_cells_cutoff = 0.1
# fp_ratio_per_cell_cross_all_muts_cutoff = 0.9
fp_ratio_per_mutation_cross_all_cells_cutoff=params['fp_ratio_per_mutation_cross_all_cells_cutoff']
fp_count_per_mutation_cross_all_cells_cutoff=params['fp_count_per_mutation_cross_all_cells_cutoff']
fp_ratio_per_cell_cross_all_muts_cutoff=params['fp_ratio_per_cell_cross_all_muts_cutoff']


##### 筛选导致 fp 高的位点并处理
# rehanged_fp_mutations_cross_all_cells = df_fp_ratio_per_mutation_cross_all_cells[df_fp_ratio_per_mutation_cross_all_cells['fp_cells_ratio_per_mutation'] >= fp_ratio_per_mutation_cross_all_cells_cutoff]['identifier'].tolist()
rehanged_fp_mutations_cross_all_cells = df_fp_ratio_per_mutation_cross_all_cells[(df_fp_ratio_per_mutation_cross_all_cells['fp_cells_ratio_per_mutation'] >= fp_ratio_per_mutation_cross_all_cells_cutoff) & (df_fp_ratio_per_mutation_cross_all_cells['fp_cells_count'] >= fp_count_per_mutation_cross_all_cells_cutoff)]['identifier'].tolist()
rehanged_fp_mutations_cross_all_cells_but_backbone = [i for i in rehanged_fp_mutations_cross_all_cells if i not in list(set(expanded_mutations_of_current_backbone_nodes+scaffold_mutations))]
print(len(rehanged_fp_mutations_cross_all_cells_but_backbone))
# 3
print(rehanged_fp_mutations_cross_all_cells_but_backbone)
# ['chr3_197511839_C_T', 'chr11_105044420_G_A', 'chr20_58360500_A_T']

# 获取 fp_ratio>=0.4 的位点之间在树上的关系
ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(T_current, rehanged_fp_mutations_cross_all_cells_but_backbone)
# ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist_relaxed(T_current, rehanged_fp_mutations_cross_all_cells_but_backbone)
earlist_mutation_of_ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone = list(ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone.keys())

# 找到 fp_ratio>=0.4 的位点导致发生 fp 的 mutations
filtered_fp_mutations_dict_cross_all_cells = {mut: other_muts 
    for mut, other_muts in fp_mutations_dict_cross_all_cells.items() 
    if mut in rehanged_fp_mutations_cross_all_cells_but_backbone}

# 找到 fp_ratio>=0.4 的位点的 daughter mutations
nodes_rehanged_fp_mutations_cross_all_cells_but_backbone = list(set([find_mutation_column(mutation, M_current.columns) for mutation in rehanged_fp_mutations_cross_all_cells_but_backbone]))

node_dict = {node.name: node for node in T_current.traverse()}
daughters_nodes_dict_cross_all_cells = {
    mutation: get_all_daughter_mutations(node_dict[mutation])
    for mutation in nodes_rehanged_fp_mutations_cross_all_cells_but_backbone
    if mutation in node_dict  # 确保突变存在于树中
}
daughters_mutations_dict_cross_all_cells = [i for i in list(set([daughter for daughters_list in daughters_nodes_dict_cross_all_cells.values() for daughter in daughters_list])) if i not in rehanged_fp_mutations_cross_all_cells_but_backbone]


# 需要重新处理的位点 pool 到一起并且去重分组
daughters_to_leaf_mutations_cross_all_cells = []
fp_mutations_cross_all_cells = []
for idx, branch_mut in enumerate(ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone):
    daughter_list = ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone[branch_mut]
    
    daughters_to_leaf_mutations_cross_all_cells = list({
        item 
        for key in daughter_list
        if key in daughters_mutations_dict_cross_all_cells
        for item in [key] + daughters_mutations_dict_cross_all_cells[key]
    })
    
    fp_mutations_cross_all_cells_initail = list({
        item 
        for key in daughter_list
        if key in filtered_fp_mutations_dict_cross_all_cells
        for item in [key] + filtered_fp_mutations_dict_cross_all_cells[key]
    })
    fp_mutations_cross_all_cells = [i for i in fp_mutations_cross_all_cells_initail if i not in daughter_list and i not in daughters_to_leaf_mutations_cross_all_cells and i not in list(set(expanded_mutations_of_current_backbone_nodes+scaffold_mutations))]

# 在 cell of orgin 的 candidate 中加一个大于 5 的条件
daughters_to_leaf_mutations_cross_all_cells_qc = [mut for mut in daughters_to_leaf_mutations_cross_all_cells if df_features_new[mut]['mutant_cellnum'] > 5]
fp_mutations_cross_all_cells_qc = set(fp_mutations_cross_all_cells+[i for i in daughters_to_leaf_mutations_cross_all_cells if i not in daughters_to_leaf_mutations_cross_all_cells_qc])

# sorted by our sorting method
sorted_fp_mutations_cross_all_cells = [i for i in I_attached.columns if i in fp_mutations_cross_all_cells_qc]
sorted_daughters_to_leaf_mutations_cross_all_cells = [i for i in I_attached.columns if i in daughters_to_leaf_mutations_cross_all_cells_qc]
sorted_rehanged_mutations_all_cross_all_cells = sorted_fp_mutations_cross_all_cells + sorted_daughters_to_leaf_mutations_cross_all_cells

to_be_removed_mutations_by_fp_mutations_cross_all_cells = rehanged_fp_mutations_cross_all_cells_but_backbone
sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells = [i for i in sorted_rehanged_mutations_all_cross_all_cells if i not in rehanged_fp_mutations_cross_all_cells_but_backbone]
remove_mutations_for_rebuild = list(set(to_be_removed_mutations_by_fp_mutations_cross_all_cells+sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells))


##### 按照顺序依次处理 fp_mutations_cross_all_cells 和 daughters_to_leaf_mutations_cross_all_cells, 将他们从树上抽离再重挂（而且按理说它们应该是一定会挂在树上的）
external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells = []
if len(remove_mutations_for_rebuild) > 0:
    
    T_removed_cross_all_cells, M_removed_cross_all_cells = remove_mutations_from_tree_and_matrix(T_checkpoint_artifact_and_doublet, M_checkpoint_artifact_and_doublet, remove_mutations_for_rebuild)
    print_tree(T_removed_cross_all_cells)
    logger.info(f"The shape of removed_tree to be refined is : {M_removed_cross_all_cells.shape}")    
    
    T_current = copy.deepcopy(T_removed_cross_all_cells)
    M_current = M_removed_cross_all_cells.copy()

if len(sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells) > 0:
    # 首先重挂 sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells
    external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
        sorted_attached_mutations=sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells,
        T_current=T_current,
        M_current=M_current,
        I_attached=I_attached,
        P_attached=P_attached,
        ω_NA=ω_NA,
        fnfp_ratio=fnfp_ratio,
        φ=φ,
        logger=logger,
        root_mutations=root_mutations  # 可选，如果已有根突变列表
    )

T_test = copy.deepcopy(T_current)
M_test = M_current.copy()
M_test = M_test.drop(columns=['ROOT'], errors='ignore')
mutations_on_T_test = M_test.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
M_test = split_merged_columns(M_test, mutations_on_T_test)
final_cleaned_M_test = M_test.loc[(M_test != 0).any(axis=1)]  # 移除全0行
final_cleaned_M_test.shape
# (1253, 93)


##### 根据 fp ratio per cell 确定 doublet cells 并且删除
identified_doublet_cells = df_fp_ratio_per_cell_cross_all_muts.loc[(df_fp_ratio_per_cell_cross_all_muts['fp_muts_ratio_per_cell']>=fp_ratio_per_cell_cross_all_muts_cutoff)]['cell_id']

# 删掉有问题的 cells，修改 M_current 即可
to_be_removed_cells = list(identified_doublet_cells)
logger.info(f"The number of identified doublet cells by fp_ratio is : {len(to_be_removed_cells)}")

with open(os.path.join(outputpath_full, "likely_doublet_cells_removed_from_tree_by_fp_ratio.csv"), 'w') as f:
    for cell in to_be_removed_cells:
        f.write(cell + '\n')

M_current = M_current.drop(to_be_removed_cells, errors='ignore')

print(M_current.shape)
# (1021, 34)










# -----------------------------
# Step 7. 建树完成，初始处理
# -----------------------------
logger.info("===== Step7: The establishment is complete, and the processing and save results ...")

##### 结束建树，后续处理加输出
M_current_filtered = M_current.drop(columns=['ROOT'], errors='ignore')

for mut_on_root in root_mutations:
    M_current_filtered.insert(0, mut_on_root, 1)

mutations_on_T_current = M_current_filtered.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()

# Final scaffold tree and matrix
M_full_initial = M_current_filtered.copy()
T_full = copy.deepcopy(T_current)
M_full = split_merged_columns(M_current_filtered, mutations_on_T_current)

logger.info("Final scaffold tree:")
print_tree(T_full)
# └─ ROOT
#   └─ chr17_7578893_G_T
#     └─ chr21_46323662_C_T|chr16_497748_G_T
#       └─ chr10_118753098_G_T
#       └─ chr16_15464093_A_G
#         └─ chr19_16085831_G_T
#       └─ chr10_73254076_G_T
#   └─ chr6_34246136_G_C
#     └─ chr1_45618093_G_T
#       └─ chr2_171687790_A_T
#         └─ chr21_44856223_C_A
#           └─ chr11_9425206_C_A
#             └─ chr17_50662360_G_T
#               └─ chr7_87196234_G_T
#                 └─ chr2_171324260_C_T
#                   └─ chr18_9135639_C_G
#                   └─ chr12_51215297_G_T
#   └─ chr22_37967634_G_T
#     └─ chr12_102013001_C_A
#     └─ chr4_74388653_G_T
#   └─ chrX_41232675_T_A
#     └─ chr1_31260040_C_A
#       └─ chr7_132879457_T_A
#         └─ chr10_104003258_C_T
#         └─ chr7_128458208_G_T
#   └─ chr8_80170791_C_T
#     └─ chr4_122860726_G_T
#       └─ chr15_50499489_A_G
#       └─ chr7_20381814_G_T
#   └─ chr17_43766282_C_A
#     └─ chr5_128171668_G_C
#       └─ chr11_83165726_G_C
#   └─ chr5_128188419_C_A
#   └─ chr1_171541779_G_T

logger.info(f"Final scaffold matrix shape: {M_full.shape}")


# ------------------------------
# Step 7.1: Output Results
# ------------------------------
logger.info("===== Step7.1: Outputting results ...")

# Create output directory
phylo_dir = os.path.join(outputpath_full, "phylo")
if not os.path.exists(phylo_dir):
    os.makedirs(phylo_dir)

##### Prepare binary matrix with NA=3
I_full_withNA3 = I_attached.replace({np.nan: 3}).astype(int)
I_full_withNA3.to_csv(os.path.join(phylo_dir, "I_full_withNA3.txt"), sep="\t")
WriteTfile(os.path.join(phylo_dir, "M_full_basedPivots.filtered_sites_inferred"), M_full, M_full.index.tolist(), M_full.columns.tolist(), judge="yes")

# Clean M_full: remove all zeros columns(muts) or rows(cells)
final_cleaned_M_full = M_full.loc[:, (M_full != 0).any(axis=0)]  # 移除全0列
final_cleaned_M_full = final_cleaned_M_full.loc[(final_cleaned_M_full != 0).any(axis=1)]  # 移除全0行

# 获取保留的行列名
kept_rows = final_cleaned_M_full.index
kept_cols = final_cleaned_M_full.columns

# 从 I_full_withNA3 提取
final_cleaned_I_full_withNA3 = I_full_withNA3.loc[kept_rows, kept_cols]

WriteTfile(os.path.join(phylo_dir, "final_cleaned_M_full_basedPivots.filtered_sites_inferred"), 
           final_cleaned_M_full, final_cleaned_M_full.index.tolist(), final_cleaned_M_full.columns.tolist(), judge="yes")
final_cleaned_I_full_withNA3.to_csv(os.path.join(phylo_dir, "final_cleaned_I_full_withNA3_for_circosPlot.txt"), sep="\t")


# ##### Output binary matrix with NA=0
# I_full_withNA0 = I_attached.replace({np.nan: 0}).astype(int)
# final_cleaned_I_full_withNA0 = I_full_withNA0.loc[kept_rows, kept_cols]
# final_cleaned_I_full_withNA0.to_csv(os.path.join(phylo_dir, "final_cleaned_I_full_withNA0_for_other_methods.txt"), sep="\t")


# ------------------------------
# Step 7.2: Identify Flipping Spots
# ------------------------------
logger.info("===== Step7.2: Identifying flipping spots ...")

df_bin_withNA3_for_flipping = final_cleaned_I_full_withNA3.copy()
df_phylogeny = final_cleaned_M_full.copy()

# get false_negative_flipping spots
false_negative_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=0, condition_phylogeny=1)
)

# get NAto1 spots
NAto1_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=1)
)

# get false_positive_flipping spots
false_positive_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=1, condition_phylogeny=0)
)

# get NAto0 spots
NAto0_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=0)
)

# process na list
if false_negative_flipping_spots.empty:
    false_negative_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if false_positive_flipping_spots.empty:
    false_positive_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if NAto1_flipping_spots.empty:
    NAto1_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if NAto0_flipping_spots.empty:
    NAto0_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

# get flipping_spots dataframe
df_flipping_spots = pd.DataFrame({
    'Mutation': df_bin_withNA3_for_flipping.columns,
    'false_negative_flipping_spots': [', '.join(false_negative_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'false_positive_flipping_spots': [', '.join(false_positive_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'NAto1_flipping_spots': [', '.join(NAto1_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'NAto0_flipping_spots': [', '.join(NAto0_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns]
})
df_flipping_spots.to_csv(os.path.join(phylo_dir, "df_flipping_spots.txt"), sep="\t", index=False)


# ------------------------------
# Step 7.3: Calculate Total Flipping Counts
# ------------------------------
logger.info("===== Step7.3: Calculating total flipping counts ...")

total_FN_flipping = ((df_bin_withNA3_for_flipping == 0) & (df_phylogeny == 1)).sum().sum()
total_FP_flipping = ((df_bin_withNA3_for_flipping == 1) & (df_phylogeny == 0)).sum().sum()
total_NAto0 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 0)).sum().sum()
total_NAto1 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 1)).sum().sum()

logger.info(f"Total False Negative flipping count: {total_FN_flipping}")
logger.info(f"Total False Positive flipping count: {total_FP_flipping}")
logger.info(f"Total NA to 0 flipping count: {total_NAto0}")
logger.info(f"Total NA to 1 flipping count: {total_NAto1}")

df_total_flipping_count = pd.DataFrame({
    'total_flipping_False_Negative': [total_FN_flipping],
    'total_flipping_False_Positive': [total_FP_flipping],
    'total_flipping_NA_to_0': [total_NAto0],
    'total_flipping_NA_to_1': [total_NAto1]
})
df_total_flipping_count.to_csv(os.path.join(phylo_dir, "df_total_flipping_count.txt"), sep="\t", index=False)

logger.info("Full-resolved tree building completed successfully!")


df_flip_counts_tree = calculate_flip_counts_per_site(df_bin_withNA3_for_flipping, df_phylogeny)
df_flip_counts_tree.to_csv(os.path.join(phylo_dir, "df_flipping_count_for_each_mut.txt"), sep="\t", index=True)


print(df_total_flipping_count.iloc[0,:])
# total_flipping_False_Negative      203
# total_flipping_False_Positive       66
# total_flipping_NA_to_0           25211
# total_flipping_NA_to_1             719
# Name: 0, dtype: int64

logger.info(f"The shape of final_cleaned_M_full.shape: {final_cleaned_M_full.shape}")
print(final_cleaned_M_full.shape)
# (1061, 34)





# ------------------------------
# Step 7.4: Tree format as text and get clone info
# ------------------------------
logger.info("===== Step7.4: Tree format as text and get clone info ...")

##### 树 TreeNode 格式保存和读取

# 将 TreeNode 转换为字典格式并且保存为 JSON 文件
tree_dict = tree_to_dict(T_full)

with open(os.path.join(phylo_dir, 'final_cleaned_tree_node.json'), 'w') as f:
    json.dump(tree_dict, f, indent=4)

# 直接返回字符串文本格式并且保存
T_full.save_to_file(os.path.join(phylo_dir, 'final_cleaned_tree_node.txt'))



##### 根据树 TreeNode 格式的 mutation clone 划分 barcode clone

mutation_clones = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_full, I_attached)
df_barcode_clones = assign_clone_labels(M_full, mutation_clones)

df_barcode_clones.to_csv(os.path.join(phylo_dir, "df_barcode_clones_from_phylo_tree.csv"), sep=',', index=False)














# ------------------------------
# End of Process
# ------------------------------
end_time = time.perf_counter


##### Time #####
finish_time = time.perf_counter()
print("Program finished in {:.4f} seconds".format(finish_time-start_time))




