#!/usr/bin/env python3

###################################################################################################
######################### scRNA放松限制分类器核心函数 #######################
###################################################################################################

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
import joblib
from collections import Counter
from pathlib import Path
from src.reproducibility import set_seed, deterministic_permutation


# 定义特征列
SELECTED_FEATURES = [
    'falt', 
    'mutant_cell_fraction', 
    'vaf_in_pseudobulk',
    'vaf_mutant_avg', 
    'right_ref_querypos_num_remove_clip',
    'alt_UMI_avg_consistence_remove_single_read',
    'alt_mismatches_mean',
    'alt_consistence_soft_prop', 
    'alt_dp_mean_diff',
    'norm_alt_dp_in_pseudobulk', 
    'mismatches_p_adj',
    'sig_pvalue', 
    'alt_read_number_perUMI_median',
    'baseq_p_adj'
]

def build_relaxed_classifier_excluding_sample(df_training, exclude_sample_id):
    """
    构建放松限制的分类器（排除特定样本）
    """
    print(f"\n=== 构建放松限制的分类器 (排除样本 {exclude_sample_id}) ===")
    
    # 排除目标样本
    df_train = df_training[df_training['sampleid'] != exclude_sample_id].copy()
    print(f"训练数据大小: {df_train.shape}")
    print(f"训练数据类别分布: {Counter(df_train['label2'])}")
    
    # 准备特征和标签
    X_train = df_train[SELECTED_FEATURES].copy()
    y_train = df_train['label2'].copy()
    
    # 创建预处理管道
    pipeline = make_pipeline(
        SimpleImputer(strategy='median'),
        StandardScaler()
    )
    
    # 数据预处理和训练
    X_train_processed = pipeline.fit_transform(X_train)
    
    # 使用放松的参数设置
    rf_model = RandomForestClassifier(
        n_estimators=1000,
        random_state=42,
        class_weight={'mosaic': 1.0, 'artifact': 0.8},
        max_depth=15,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features='sqrt',
        bootstrap=True
    )
    rf_model.fit(X_train_processed, y_train)
    
    print(f"模型训练完成，类别: {rf_model.classes_}")
    return rf_model, pipeline

def predict_with_relaxed_threshold(model, pipeline, df_new, sample_id, mutation_ids=None):
    """
    使用放松的阈值进行预测
    """
    print(f"\n=== 为样本 {sample_id} 预测新突变 (放松阈值) ===")
    print(f"新数据形状: {df_new.shape}")
    
    # 检查必需的特征列
    missing_features = [feat for feat in SELECTED_FEATURES if feat not in df_new.columns]
    if missing_features:
        raise ValueError(f"缺少必需的特征列: {missing_features}")
    
    # 准备特征数据
    X_new = df_new[SELECTED_FEATURES].copy()
    
    # 使用相同的预处理管道
    X_new_processed = pipeline.transform(X_new)
    
    # 进行预测 - 使用概率而不是硬分类
    probabilities = model.predict_proba(X_new_processed)
    class_labels = model.classes_
    
    # 放松的预测逻辑
    predictions = []
    for i, prob_vector in enumerate(probabilities):
        mosaic_prob = prob_vector[list(class_labels).index('mosaic')] if 'mosaic' in class_labels else 0
        artifact_prob = prob_vector[list(class_labels).index('artifact')] if 'artifact' in class_labels else 0
        
        # 放松的决策规则
        if mosaic_prob > 0.4:  # 降低mosaic的阈值
            predictions.append('mosaic')
        elif artifact_prob > 0.6:  # 提高artifact的阈值
            predictions.append('artifact')
        else:
            if mosaic_prob > 0.2:  # 进一步降低mosaic的阈值
                predictions.append('mosaic')
            else:
                predictions.append(class_labels[np.argmax(prob_vector)])
    
    # 创建结果DataFrame
    results_df = pd.DataFrame({
        'sample_id': sample_id,
        'predicted_label': predictions
    })
    
    # 添加mutation_id列
    if mutation_ids is not None:
        results_df['mutation_id'] = mutation_ids
    elif 'mutation_id' in df_new.columns:
        results_df['mutation_id'] = df_new['mutation_id'].values
    elif 'identifier' in df_new.columns:
        results_df['mutation_id'] = df_new['identifier'].values
    else:
        # 如果没有mutation_id列，创建索引作为标识
        results_df['mutation_id'] = [f'mutation_{i+1}' for i in range(len(results_df))]
    
    # 添加概率分数
    for i, class_name in enumerate(class_labels):
        results_df[f'probability_{class_name}'] = probabilities[:, i]
    
    # 重新排列列的顺序，让mutation_id在前面
    cols = ['mutation_id', 'sample_id', 'predicted_label'] + [f'probability_{cls}' for cls in class_labels]
    results_df = results_df[cols]
    
    # 打印预测统计
    prediction_counts = pd.Series(predictions).value_counts().to_dict()
    print(f"放松阈值预测统计:")
    for label, count in prediction_counts.items():
        percentage = count / len(predictions) * 100
        print(f"  {label}: {count} 个位点 ({percentage:.1f}%)")
    
    return results_df

def real_time_classifier_predict(df_for_classifier_all, sampleid, outputpath):
    """
    scRNA放松限制分类器核心预测函数
    
    Args:
        df_for_classifier_all: 需要预测的数据框（包含所有特征）
        sampleid: 样本ID
        outputpath: 输出路径
        
    Returns:
        pd.DataFrame: 预测结果
    """
    print(f"\n{'='*60}")
    print(f"开始scRNA放松限制分类器预测 - 样本: {sampleid}")
    print(f"{'='*60}")
    
    # 确保输出目录存在
    os.makedirs(outputpath, exist_ok=True)
    
    # 1. 加载训练数据
    print("=== 加载训练数据 ===")
    script_dir = Path(__file__).parent
    features_file_labeled = script_dir / 'classifer' / 'scrna' / 'data_labeling_sampling.ratio_155_space.csv'
    df_training = pd.read_csv(features_file_labeled, sep="\t")
    print(f"训练数据形状: {df_training.shape}")
    print(f"训练数据类别分布: {Counter(df_training['label2'])}")
    
    # 2. 数据预处理
    print("=== 数据预处理 ===")
    df_features = df_for_classifier_all.copy()
    
    # 提取mutation_id（如果存在相关列）
    mutation_ids = None
    if 'mutation_id' in df_features.columns:
        mutation_ids = df_features['mutation_id'].values
    elif 'identifier' in df_features.columns:
        mutation_ids = df_features['identifier'].values
    
    # 数据清洗
    df_features_selected = df_features[SELECTED_FEATURES].copy()
    df_features_selected.replace('no', np.nan, inplace=True)
    
    for col in df_features_selected.columns:
        df_features_selected[col] = pd.to_numeric(df_features_selected[col], errors='coerce')
        median_value = df_features_selected[col].median(skipna=True)
        df_features_selected[col] = df_features_selected[col].fillna(median_value)
        finite_max = df_features_selected[col][np.isfinite(df_features_selected[col])].max()
        df_features_selected[col] = df_features_selected[col].replace([np.inf, -np.inf], finite_max)
    
    print(f"预处理后数据形状: {df_features_selected.shape}")
    
    # 3. 构建放松限制的分类器
    model, pipeline = build_relaxed_classifier_excluding_sample(df_training, sampleid)
    
    # 4. 使用放松的阈值进行预测
    results = predict_with_relaxed_threshold(model, pipeline, df_features_selected, sampleid, mutation_ids)
    
    # 5. 保存预测结果
    print("\n=== 保存预测结果 ===")
    
    # 保存所有位点的预测结果
    all_sites_file = os.path.join(outputpath, f"{sampleid}.feature_and_prediction.allsites.txt")
    results.to_csv(all_sites_file, index=False, sep="\t")
    print(f"所有位点预测结果保存到: {all_sites_file}")
    
    # 保存mosaic位点列表
    df_mosaic = results[results["predicted_label"] == "mosaic"]
    mosaic_list_file = os.path.join(outputpath, f"{sampleid}_mosaic_prediction.list.txt")
    df_mosaic['mutation_id'].to_csv(mosaic_list_file, index=False, header=False)
    print(f"mosaic位点列表保存到: {mosaic_list_file}")
    
    # 保存mosaic详细结果
    mosaic_detailed_file = os.path.join(outputpath, f"{sampleid}_mosaic_detailed_results.txt")
    df_mosaic.to_csv(mosaic_detailed_file, index=False, sep="\t")
    print(f"mosaic详细结果保存到: {mosaic_detailed_file}")
    
    # 保存统计信息
    mosaic_count_file = os.path.join(outputpath, f"{sampleid}_prediction_summary.txt")
    with open(mosaic_count_file, 'w') as f:
        f.write(f"样本 {sampleid} 预测结果汇总\n")
        f.write("=" * 40 + "\n")
        f.write(f"总位点数量: {len(results)}\n")
        f.write(f"mosaic位点数量: {len(df_mosaic)}\n")
        f.write(f"artifact位点数量: {len(results) - len(df_mosaic)}\n")
        f.write(f"mosaic比例: {len(df_mosaic)/len(results)*100:.1f}%\n")
        
        # 添加概率统计
        if 'probability_mosaic' in results.columns:
            f.write(f"\nMosaic概率统计:\n")
            f.write(f"  平均值: {results['probability_mosaic'].mean():.3f}\n")
            f.write(f"  中位数: {results['probability_mosaic'].median():.3f}\n")
            f.write(f"  最大值: {results['probability_mosaic'].max():.3f}\n")
            f.write(f"  最小值: {results['probability_mosaic'].min():.3f}\n")
    
    print(f"预测汇总保存到: {mosaic_count_file}")
    
    # 打印最终统计
    mosaic_count = len(df_mosaic)
    total_count = len(results)
    print(f"\n最终预测统计: {mosaic_count}/{total_count} 个位点被预测为mosaic ({mosaic_count/total_count*100:.1f}%)")
    
    print(f"\n{'='*60}")
    print(f"scRNA放松限制分类器预测完成 - 样本: {sampleid}")
    print(f"{'='*60}")
    
    return results

# # 使用示例
# if __name__ == "__main__":
#     # 您的数据框（包含mutation_id）
#     df_for_classifier_all = pd.DataFrame({
#         'mutation_id': [
#             'chr20_14326432_C_A', 'chr21_9001082_G_C', 'chr6_137868589_G_A',
#             'chr1_106340339_G_T', 'chr1_212476040_A_G', 'chr19_38404452_C_T'
#         ],
#         'falt': [0.1, 0.2, 0.05, 0.15, 0.08, 0.12],
#         'mutant_cell_fraction': [0.3, 0.4, 0.1, 0.35, 0.25, 0.28],
#         'vaf_in_pseudobulk': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012],
#         'vaf_mutant_avg': [0.15, 0.25, 0.08, 0.18, 0.12, 0.16],
#         'right_ref_querypos_num_remove_clip': [5, 8, 2, 6, 4, 5],
#         'alt_UMI_avg_consistence_remove_single_read': [0.1, 0.15, 0.05, 0.12, 0.08, 0.11],
#         'alt_mismatches_mean': [0.02, 0.03, 0.01, 0.025, 0.015, 0.022],
#         'alt_consistence_soft_prop': [0.05, 0.08, 0.02, 0.06, 0.04, 0.055],
#         'alt_dp_mean_diff': [0.12, 0.22, 0.06, 0.16, 0.1, 0.14],
#         'norm_alt_dp_in_pseudobulk': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012],
#         'mismatches_p_adj': [0.05, 0.08, 0.02, 0.06, 0.04, 0.055],
#         'sig_pvalue': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012],
#         'alt_read_number_perUMI_median': [2, 3, 1, 2, 2, 2],
#         'baseq_p_adj': [0.05, 0.08, 0.02, 0.06, 0.04, 0.055]
#     })
    
#     sampleid = "10k"
#     outputpath = "./scRNA_relaxed_predictions"
    
#     # 执行预测
#     results = real_time_classifier_predict(df_for_classifier_all, sampleid, outputpath)
    
#     print(f"\n预测结果预览:")
#     print(results.head())

