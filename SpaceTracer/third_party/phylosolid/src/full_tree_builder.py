#!/usr/bin/env python3
"""
full_tree_builder.py

This module encapsulates Steps 6-8 of the PhyloSOLID pipeline for building
a fully-resolved phylogenetic tree from a scaffold tree.

Sub-modules:
    - build_initial_tree(): Step 6 - Initial full-resolved tree construction
    - build_refined_tree(): Step 7 - Prune-and-regraft (PRG) refinement
    - tree_QC_and_filter(): Step 8 - Quality control and filtration
    - build_fully_resolved_tree(): Master function orchestrating all three steps
    - export_phylo_results_for_cv(): Export phylo results for each CV threshold

Input:  T_scaffold, M_scaffold, I_attached, P_attached, ...
Output: Refined tree (T_current, M_current) and QC statistics

Author: Qing
Date: 2026/07/22
Update: 2026/08/14 - Added export_phylo_results_for_cv function for parallel CV search
"""

import os
import copy
import json
import logging
import pandas as pd
import numpy as np
from typing import Tuple, Dict, Any, Optional, List

from src.mutation_integrator import *
from src.scaffold_builder import *


logger = logging.getLogger(__name__)


# ============================================================================
# Step 6: Initial full-resolved tree construction (Steps 6.1-6.4)
# ============================================================================

def build_initial_tree(
    T_scaffold: TreeNode,
    M_scaffold: pd.DataFrame,
    scaffold_mutations: List[str],
    I_attached: pd.DataFrame,
    P_attached: pd.DataFrame,
    attached_mutations: List[str],
    params: Dict[str, Any],
    all_conflict_mutations: List[str],
    logger_obj: Optional[logging.Logger] = None,
) -> Tuple[TreeNode, pd.DataFrame, List[str], List[str]]:
    """
    Step 6: Build initial fully-resolved tree by placing accessory mutations
    onto the scaffold tree using the Bayesian placement procedure.
    
    Sub-steps:
        6.1: First-pass integration of accessory mutations
        6.2: Second-pass integration of extraneous mutations
        6.3: Cluster-based integration of remaining extraneous mutations
        6.4: Final integration - root assignment of remained mutations
    
    Parameters
    ----------
    T_scaffold : TreeNode
        Scaffold tree
    M_scaffold : pd.DataFrame
        Scaffold mutation matrix
    scaffold_mutations : List[str]
        List of scaffold mutations
    I_attached : pd.DataFrame
        Mutation presence matrix
    P_attached : pd.DataFrame
        Posterior probability matrix
    attached_mutations : List[str]
        Accessory mutations to be placed
    params : Dict[str, Any]
        Parameters dictionary
    all_conflict_mutations : List[str]
        Existing conflict mutations (will be extended)
    logger_obj : logging.Logger, optional
        Logger instance
    
    Returns
    -------
    T_current : TreeNode
        Initial fully-resolved tree
    M_current : pd.DataFrame
        Initial mutation matrix
    root_mutations : List[str]
        Mutations assigned to ROOT
    all_conflict_mutations : List[str]
        Updated conflict mutations
    """
    if logger_obj is None:
        logger_obj = logging.getLogger(__name__)
    
    # Extract penalty parameters
    omega_NA = params.get('general_weight_NA', 0.001)
    fnfp_ratio = params.get('fnfp_ratio', 0.1)
    phi = params.get('phi', 1.0)
    
    # Prepare for mutation placement
    logger_obj.info("=" * 80)
    logger_obj.info("===== Step6: Initial full-resolved tree ...")
    logger_obj.info("PART I: INITIAL SCAFFOLD INTEGRATION")
    logger_obj.info("=" * 80)
    logger_obj.info("")
    
    logger_obj.info("Calculating penalties and refining placements ...")
    T_current = copy.deepcopy(T_scaffold)
    
    new_rows_for_current = I_attached.index.difference(M_scaffold.index)
    new_data_for_current = pd.DataFrame(0, index=new_rows_for_current, columns=M_scaffold.columns)
    M_current_each_mut = pd.concat([M_scaffold, new_data_for_current])
    all_nodes_in_T_scaffold = T_scaffold.all_names_no_root()
    M_current = merge_mutations(M_current_each_mut, all_nodes_in_T_scaffold)
    M_current.insert(0, 'ROOT', 1)
    
    root_mutations = []
    
    # --------------------------------------------------------------------------
    # Step 6.1: First-pass integration of accessory mutations
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 6.1: First-pass integration of accessory mutations")
    logger_obj.info("-" * 80)
    logger_obj.info("Integrating accessory mutations (M_accessory) onto the scaffold")
    logger_obj.info("tree (T_scaffold) using the Bayesian placement procedure.")
    logger_obj.info("")
    logger_obj.info("  M_scaffold: scaffold/backbone mutations (retained)")
    logger_obj.info("  M_accessory: accessory mutations to be placed")
    logger_obj.info("  M_extraneous: mutations with no clear intersection (flagged)")
    logger_obj.info("-" * 80)
    
    sorted_attached_mutations = [i for i in I_attached.columns if i in attached_mutations]
    
    if len(sorted_attached_mutations) > 0:
        external_mutations_of_attached_on_scaffold, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_attached_mutations,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    else:
        external_mutations_of_attached_on_scaffold = []
    
    logger_obj.info(f"  |M_scaffold|: {len(scaffold_mutations)}")
    logger_obj.info(f"  |M_accessory|: {len(sorted_attached_mutations)}")
    logger_obj.info(f"  |M_accessory successfully placed|: {len(sorted_attached_mutations) - len(external_mutations_of_attached_on_scaffold)}")
    logger_obj.info(f"  |M_extraneous|: {len(external_mutations_of_attached_on_scaffold)}")
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 6.2: Second-pass integration of extraneous mutations
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 6.2: Second-pass integration of extraneous mutations")
    logger_obj.info("-" * 80)
    logger_obj.info("Re-attempting integration of extraneous mutations (M_extraneous)")
    logger_obj.info("that failed to find clear intersections in the first pass.")
    logger_obj.info("-" * 80)
    
    if len(external_mutations_of_attached_on_scaffold) > 0:
        sorted_external_mutations = [i for i in I_attached.columns if i in external_mutations_of_attached_on_scaffold]
        final_external_mutations_of_attached_on_scaffold, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_external_mutations,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    else:
        final_external_mutations_of_attached_on_scaffold = []
    
    logger_obj.info(f"  |M_extraneous remaining after second pass|: {len(final_external_mutations_of_attached_on_scaffold)}")
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 6.3: Cluster-based integration of remaining extraneous mutations
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 6.3: Cluster-based integration of remaining extraneous mutations")
    logger_obj.info("-" * 80)
    logger_obj.info("Grouping remaining M_extraneous by intersection patterns and")
    logger_obj.info("integrating them as coherent phylogenetic clusters.")
    logger_obj.info("-" * 80)
    
    external_mutations = final_external_mutations_of_attached_on_scaffold
    logger_obj.info(f"The number of external_mutations is: {len(external_mutations)}")
    
    final_external_mutations = []
    if len(external_mutations) > 0:
        sorted_external_mutations = [i for i in I_attached.columns if i in external_mutations]
        final_external_mutations, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_external_mutations,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    logger_obj.info(f"The number of final_external_mutations is: {len(final_external_mutations)}")
    
    remained_mutations = []
    if len(final_external_mutations) > 0:
        subtree_groups = cluster_external_mutations_by_intersection(I_attached, final_external_mutations)
        logger_obj.info("Processing remaining external mutations by building subtrees")
        
        remained_mutations, conflict_mutations_temp, T_current, M_current, root_mutations = process_external_mutations_by_subtree_groups(
            subtree_groups=subtree_groups,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    logger_obj.info(f"  |M_extraneous resolved via clustering|: {len(final_external_mutations)}")
    logger_obj.info(f"  |M_remained|: {len(remained_mutations)}")
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 6.4: Final integration - root assignment of remained mutations
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 6.4: Final integration - root assignment of remained mutations")
    logger_obj.info("-" * 80)
    logger_obj.info("Performing final integration attempt for M_remained mutations.")
    logger_obj.info("Mutations that remain unresolved are assigned to ROOT as M_root.")
    logger_obj.info("-" * 80)
    
    final_remained_mutations = []
    if len(remained_mutations) > 0:
        sorted_remained_mutations = [i for i in I_attached.columns if i in remained_mutations]
        final_remained_mutations, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_remained_mutations,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    logger_obj.info(f"  |M_root|: {len(root_mutations)}")
    logger_obj.info(f"  |M_ambiguous|: {len(all_conflict_mutations)}")
    logger_obj.info("")
    
    # Step 6 Summary
    logger_obj.info("=" * 80)
    logger_obj.info("Step 6 COMPLETED: Initial fully-resolved tree constructed")
    logger_obj.info("-" * 80)
    
    logger_obj.info("Current tree structure:")
    logger_obj.info("-" * 40)
    print_tree_logger(T_current)
    logger_obj.info("-" * 40)
    
    mutations_on_tree = [m for m in M_current.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist() if m != 'ROOT']
    logger_obj.info(f"  Mutations on tree (excluding ROOT): {len(mutations_on_tree)}")
    logger_obj.info(f"  Cells: {M_current.shape[0]}")
    logger_obj.info(f"  Merged columns (including ROOT): {M_current.shape[1]}")
    logger_obj.info("=" * 80)
    logger_obj.info("")
    
    return T_current, M_current, root_mutations, all_conflict_mutations


# ============================================================================
# Step 7: Prune-and-regraft (PRG) refinement (Steps 7.1-7.4)
# ============================================================================

def build_refined_tree(
    T_current: TreeNode,
    M_current: pd.DataFrame,
    scaffold_mutations: List[str],
    I_attached: pd.DataFrame,
    P_attached: pd.DataFrame,
    df_features_new: pd.DataFrame,
    params: Dict[str, Any],
    all_conflict_mutations: List[str],
    root_mutations: List[str],
    outputpath_full: str,
    logger_obj: Optional[logging.Logger] = None,
) -> Tuple[TreeNode, pd.DataFrame, List[str], List[str], float, List[str]]:
    """
    Step 7: Prune-and-regraft (PRG) refinement.
    
    Identifies mutations with elevated discordance rates using four complementary
    metrics, prunes them from the tree, and re-integrates them using the
    Bayesian placement procedure.
    
    Metrics:
        7.1: r_FP^(clone) - Intra-clonal false positive discordance rate
        7.2: r_FP^(global), r_FN^(global) - Global discordance rates
        7.3: r_FP^(locus) - Locus-specific false positive discordance rate
        7.4: eta - Ancestral retention fraction
    
    Parameters
    ----------
    T_current : TreeNode
        Initial fully-resolved tree
    M_current : pd.DataFrame
        Initial mutation matrix
    scaffold_mutations : List[str]
        List of scaffold mutations
    I_attached : pd.DataFrame
        Mutation presence matrix
    P_attached : pd.DataFrame
        Posterior probability matrix
    df_features_new : pd.DataFrame
        Features matrix
    params : Dict[str, Any]
        Parameters dictionary
    all_conflict_mutations : List[str]
        Existing conflict mutations (will be extended)
    root_mutations : List[str]
        Existing root mutations (will be extended)
    outputpath_full : str
        Output directory path for saving QC files
    logger_obj : logging.Logger, optional
        Logger instance
    
    Returns
    -------
    T_current : TreeNode
        Refined tree after PRG
    M_current : pd.DataFrame
        Refined mutation matrix after PRG
    root_mutations : List[str]
        Updated root mutations
    all_conflict_mutations : List[str]
        Updated conflict mutations
    omega_before_qc : float
        Omega pre-QC (computed after PRG)
    outgroup_mutations : List[str]
        Orphaned mutations identified in Step 7.4
    """
    if logger_obj is None:
        logger_obj = logging.getLogger(__name__)
    
    # Extract parameters
    fnfp_ratio = params.get('fnfp_ratio', 0.1)
    
    logger_obj.info("=" * 80)
    logger_obj.info("Step 7: Prune-and-regraft (PRG) refinement")
    logger_obj.info("=" * 80)
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 7.1: Intra-clonal FP discordance pruning (M_FP^(clone))
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 7.1: Intra-clonal FP discordance pruning")
    logger_obj.info("-" * 80)
    logger_obj.info("  r_FP^(clone)(j) = delta_FP(j) / |C_j|")
    logger_obj.info("  where C_j is the set of cells in the clone where mutation j")
    logger_obj.info("  is phylogenetically placed.")
    logger_obj.info("")
    logger_obj.info("  M_FP^(clone): mutations flagged by intra-clonal FP discordance")
    logger_obj.info("  M_FP-daughter^(clone): daughter mutations co-flagged")
    logger_obj.info("  G_FP^(clone): clone-associated group containing M_FP^(clone)")
    logger_obj.info("-" * 80)
    
    # Recalculate backbone nodes
    current_backbone_nodes = get_first_level_backbone_nodes(T_current)
    expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]
    
    mutation_clones_for_subclone = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)
    
    # Compute intra-clonal FP discordance rates
    T_checkpoint_fpratio_within_subclone = copy.deepcopy(T_current)
    M_checkpoint_fpratio_within_subclone = M_current.copy()
    
    M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone = M_checkpoint_fpratio_within_subclone.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_fpratio_within_subclone = M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone, mutations_on_T_current_fpratio_within_subclone)
    
    df_fp_ratio_fpratio_within_subclone, fp_mutations_dict_for_out_subclone_muts_fpratio_within_subclone, fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone = calculate_fp_ratios_within_subclone(
        M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone, I_attached, mutation_clones_for_subclone
    )
    
    fp_ratio_cutoff_within_subclone = params.get('fp_ratio_cutoff_within_subclone', 0.1)
    
    # Identify M_FP^(clone)
    rehanged_mutations_by_fpratio_within_subclone = df_fp_ratio_fpratio_within_subclone[
        df_fp_ratio_fpratio_within_subclone['fp_ratio_within_subclone_for_in_subclone_muts'] >= fp_ratio_cutoff_within_subclone
    ]['identifier'].tolist()
    
    rehanged_mutations_by_fpratio_within_subclone_but_backbone = [
        i for i in rehanged_mutations_by_fpratio_within_subclone 
        if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
    ]
    
    logger_obj.info(f"  |M_FP^(clone)|: {len(rehanged_mutations_by_fpratio_within_subclone)}")
    logger_obj.info(f"  |M_FP^(clone) \\ M_backbone|: {len(rehanged_mutations_by_fpratio_within_subclone_but_backbone)}")
    
    # Identify M_FP-daughter^(clone)
    ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
        T_current, rehanged_mutations_by_fpratio_within_subclone_but_backbone
    )
    
    filtered_fp_mutations_dict_by_fpratio_within_subclone = {
        mut: other_muts 
        for mut, other_muts in fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone.items() 
        if mut in rehanged_mutations_by_fpratio_within_subclone_but_backbone
    }
    
    nodes_rehanged_mutations_by_fpratio_within_subclone_but_backbone = list(set([
        find_mutation_column(mutation, M_current.columns) 
        for mutation in rehanged_mutations_by_fpratio_within_subclone_but_backbone
    ]))
    
    node_dict = {node.name: node for node in T_current.traverse()}
    daughters_nodes_dict_by_fpratio_within_subclone = {
        mutation: get_all_daughter_mutations(node_dict[mutation])
        for mutation in nodes_rehanged_mutations_by_fpratio_within_subclone_but_backbone
        if mutation in node_dict
    }
    daughters_mutations_dict_by_fpratio_within_subclone = list(set([
        daughter for daughters_list in daughters_nodes_dict_by_fpratio_within_subclone.values() 
        for daughter in daughters_list
        if daughter not in rehanged_mutations_by_fpratio_within_subclone_but_backbone
    ]))
    
    logger_obj.info(f"  |M_FP-daughter^(clone)|: {len(daughters_mutations_dict_by_fpratio_within_subclone)}")
    logger_obj.info(f"  |G_FP^(clone)|: {len(ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone)} groups")
    
    # Pool mutations for PRG
    daughters_to_leaf_mutations_fpratio_within_subclone = []
    fp_mutations_fpratio_within_subclone = []
    
    for branch_mut in ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone:
        daughter_list = ordered_branch_groups_for_rehanged_mutations_by_fpratio_within_subclone_but_backbone[branch_mut]
        
        daughters_to_leaf_mutations_fpratio_within_subclone = list({
            item 
            for key in daughter_list
            if key in daughters_mutations_dict_by_fpratio_within_subclone
            for item in [key] + daughters_mutations_dict_by_fpratio_within_subclone[key]
        })
        
        fp_mutations_fpratio_within_subclone_init = list({
            item 
            for key in daughter_list
            if key in filtered_fp_mutations_dict_by_fpratio_within_subclone
            for item in [key] + filtered_fp_mutations_dict_by_fpratio_within_subclone[key]
        })
        fp_mutations_fpratio_within_subclone = [
            i for i in fp_mutations_fpratio_within_subclone_init 
            if i not in daughter_list 
            and i not in daughters_to_leaf_mutations_fpratio_within_subclone 
            and i not in expanded_mutations_of_current_backbone_nodes
        ]
    
    # QC filter: require mutant cell number > 5
    daughters_to_leaf_mutations_fpratio_within_subclone_qc = [
        mut for mut in daughters_to_leaf_mutations_fpratio_within_subclone 
        if df_features_new[mut]['mutant_cellnum'] > 5
    ]
    fp_mutations_fpratio_within_subclone_qc = set(
        fp_mutations_fpratio_within_subclone + 
        [i for i in daughters_to_leaf_mutations_fpratio_within_subclone 
         if i not in daughters_to_leaf_mutations_fpratio_within_subclone_qc]
    )
    
    sorted_fp_mutations_fpratio_within_subclone = [
        i for i in I_attached.columns 
        if i in fp_mutations_fpratio_within_subclone_qc
    ]
    sorted_daughters_to_leaf_mutations_fpratio_within_subclone = [
        i for i in I_attached.columns 
        if i in daughters_to_leaf_mutations_fpratio_within_subclone_qc
    ]
    sorted_rehanged_mutations_all_fpratio_within_subclone = (
        sorted_fp_mutations_fpratio_within_subclone + 
        sorted_daughters_to_leaf_mutations_fpratio_within_subclone
    )
    
    # Prune and regraft M_FP^(clone) and M_FP-daughter^(clone)
    external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone = []
    external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone = []
    
    if len(sorted_rehanged_mutations_all_fpratio_within_subclone) > 0:
        T_removed_fpratio_within_subclone, M_removed_fpratio_within_subclone = remove_mutations_from_tree_and_matrix(
            T_checkpoint_fpratio_within_subclone, M_checkpoint_fpratio_within_subclone, 
            sorted_rehanged_mutations_all_fpratio_within_subclone
        )
        logger_obj.info(f"  Tree after pruning: {M_removed_fpratio_within_subclone.shape[0]} cells, {M_removed_fpratio_within_subclone.shape[1]} mutations")
        T_current = copy.deepcopy(T_removed_fpratio_within_subclone)
        M_current = M_removed_fpratio_within_subclone.copy()
    
    omega_NA = params.get('general_weight_NA', 0.001)
    fnfp_ratio = params.get('fnfp_ratio', 0.1)
    phi = params.get('phi', 1.0)
    
    if len(sorted_fp_mutations_fpratio_within_subclone) > 0:
        external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_fp_mutations_fpratio_within_subclone,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    if len(sorted_daughters_to_leaf_mutations_fpratio_within_subclone) > 0:
        external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_daughters_to_leaf_mutations_fpratio_within_subclone,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    # Re-attempt any remaining external mutations
    external_mutations_fpratio_within_subclone = list(set(
        external_mutations_fpratio_within_subclone_by_sorted_fp_mutations_fpratio_within_subclone + 
        external_mutations_fpratio_within_subclone_by_sorted_daughters_to_leaf_mutations_fpratio_within_subclone
    ))
    
    if len(external_mutations_fpratio_within_subclone) > 0:
        logger_obj.info(f"  Re-attempting {len(external_mutations_fpratio_within_subclone)} external mutations...")
        sorted_external_mutations_fpratio_within_subclone = [
            i for i in I_attached.columns 
            if i in external_mutations_fpratio_within_subclone
        ]
        final_external_mutations_fpratio_within_subclone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_external_mutations_fpratio_within_subclone,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
        logger_obj.info(f"  Final external mutations: {len(final_external_mutations_fpratio_within_subclone)}")
    
    logger_obj.info("")
    
    # Recompute metrics after pruning
    mutation_clones_for_subclone_v2 = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)
    
    T_checkpoint_fpratio_within_subclone_v2 = copy.deepcopy(T_current)
    M_checkpoint_fpratio_within_subclone_v2 = M_current.copy()
    
    M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2 = M_checkpoint_fpratio_within_subclone_v2.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_fpratio_within_subclone_v2 = M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2 = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2, mutations_on_T_current_fpratio_within_subclone_v2)
    
    df_fp_ratio_fpratio_within_subclone_v2, fp_mutations_dict_for_out_subclone_muts_fpratio_within_subclone_v2, fp_mutations_dict_for_in_subclone_muts_fpratio_within_subclone_v2 = calculate_fp_ratios_within_subclone(
        M_for_fp_ratio_and_fn_ratio_fpratio_within_subclone_v2, I_attached, mutation_clones_for_subclone_v2
    )
    
    df_fp_ratio_fpratio_within_subclone_final = pd.merge(
        df_fp_ratio_fpratio_within_subclone, 
        df_fp_ratio_fpratio_within_subclone_v2, 
        on='identifier', 
        suffixes=('.1', '.2')
    )
    
    # --------------------------------------------------------------------------
    # Step 7.2: Global FP/FN discordance pruning (M_FP^(global), M_FN^(global))
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 7.2: Global FP/FN discordance pruning")
    logger_obj.info("-" * 80)
    logger_obj.info("  r_FP^(global)(j) = delta_FP(j) / N_total")
    logger_obj.info("  r_FN^(global)(j) = delta_FN(j) / N_total")
    logger_obj.info("")
    logger_obj.info("  M_FP^(global): mutations with r_FP^(global) >= threshold")
    logger_obj.info("  M_FN^(global): mutations with r_FN^(global) >= threshold")
    logger_obj.info("  M_FP-daughter^(global): daughter mutations co-flagged")
    logger_obj.info("  G_FP^(global): global FP-associated group")
    logger_obj.info("  G_FN^(global): global FN-associated group")
    logger_obj.info("-" * 80)
    
    # Recalculate backbone nodes
    current_backbone_nodes = get_first_level_backbone_nodes(T_current)
    expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]
    
    # Compute global FP/FN discordance rates
    T_checkpoint_fpfnratio_across_tree = copy.deepcopy(T_current)
    M_checkpoint_fpfnratio_across_tree = M_current.copy()
    
    M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree = M_checkpoint_fpfnratio_across_tree.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_fpfnratio_across_tree = M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree, mutations_on_T_current_fpfnratio_across_tree)
    
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree, fp_mutations_dict_fpfnratio_across_tree = calculate_fp_fn_ratios_across_tree(
        M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree, I_attached
    )
    
    fp_ratio_cutoff_across_tree = params.get('fp_ratio_cutoff_across_tree', 0.2)
    fn_ratio_cutoff_across_tree = params.get('fn_ratio_cutoff_across_tree', 0.9)
    
    # Identify M_FP^(global)
    rehanged_fp_mutations_by_fpfnratio_across_tree = df_fp_ratio_and_fn_ratio_fpfnratio_across_tree[
        df_fp_ratio_and_fn_ratio_fpfnratio_across_tree['fp_ratio'] >= fp_ratio_cutoff_across_tree
    ]['identifier'].tolist()
    
    rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = [
        i for i in rehanged_fp_mutations_by_fpfnratio_across_tree 
        if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
    ]
    
    logger_obj.info(f"  |M_FP^(global)|: {len(rehanged_fp_mutations_by_fpfnratio_across_tree)}")
    logger_obj.info(f"  |M_FP^(global) \\ M_backbone|: {len(rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone)}")
    
    # Identify G_FP^(global)
    ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
        T_current, rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
    )
    
    filtered_fp_mutations_dict_by_fpfnratio_across_tree = {
        mut: other_muts 
        for mut, other_muts in fp_mutations_dict_fpfnratio_across_tree.items() 
        if mut in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
    }
    
    nodes_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone = list(set([
        find_mutation_column(mutation, M_current.columns) 
        for mutation in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
    ]))
    
    node_dict = {node.name: node for node in T_current.traverse()}
    daughters_nodes_dict_by_fpfnratio_across_tree = {
        mutation: get_all_daughter_mutations(node_dict[mutation])
        for mutation in nodes_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
        if mutation in node_dict
    }
    
    daughters_mutations_dict_by_fpfnratio_across_tree = list(set([
        daughter for daughters_list in daughters_nodes_dict_by_fpfnratio_across_tree.values() 
        for daughter in daughters_list
        if daughter not in rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone
    ]))
    
    logger_obj.info(f"  |M_FP-daughter^(global)|: {len(daughters_mutations_dict_by_fpfnratio_across_tree)}")
    logger_obj.info(f"  |G_FP^(global)|: {len(ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone)} groups")
    
    # Pool mutations for PRG
    daughters_to_leaf_mutations_fpfnratio_across_tree = []
    fp_mutations_fpfnratio_across_tree = []
    
    for branch_mut in ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone:
        daughter_list = ordered_branch_groups_for_rehanged_fp_mutations_by_fpfnratio_across_tree_but_backbone[branch_mut]
        
        daughters_to_leaf_mutations_fpfnratio_across_tree = list({
            item 
            for key in daughter_list
            if key in daughters_mutations_dict_by_fpfnratio_across_tree
            for item in [key] + daughters_mutations_dict_by_fpfnratio_across_tree[key]
        })
        
        fp_mutations_fpfnratio_across_tree_init = list({
            item 
            for key in daughter_list
            if key in filtered_fp_mutations_dict_by_fpfnratio_across_tree
            for item in [key] + filtered_fp_mutations_dict_by_fpfnratio_across_tree[key]
        })
        fp_mutations_fpfnratio_across_tree = [
            i for i in fp_mutations_fpfnratio_across_tree_init 
            if i not in daughter_list 
            and i not in daughters_to_leaf_mutations_fpfnratio_across_tree 
            and i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
        ]
    
    daughters_to_leaf_mutations_fpfnratio_across_tree_qc = [
        mut for mut in daughters_to_leaf_mutations_fpfnratio_across_tree 
        if df_features_new[mut]['mutant_cellnum'] > 5
    ]
    fp_mutations_fpfnratio_across_tree_qc = set(
        fp_mutations_fpfnratio_across_tree + 
        [i for i in daughters_to_leaf_mutations_fpfnratio_across_tree 
         if i not in daughters_to_leaf_mutations_fpfnratio_across_tree_qc]
    )
    
    sorted_fp_mutations_fpfnratio_across_tree = [
        i for i in I_attached.columns 
        if i in fp_mutations_fpfnratio_across_tree_qc
    ]
    sorted_daughters_to_leaf_mutations_fpfnratio_across_tree = [
        i for i in I_attached.columns 
        if i in daughters_to_leaf_mutations_fpfnratio_across_tree_qc
    ]
    sorted_rehanged_mutations_all_fpfnratio_across_tree = (
        sorted_fp_mutations_fpfnratio_across_tree + 
        sorted_daughters_to_leaf_mutations_fpfnratio_across_tree
    )
    
    # Prune and regraft M_FP^(global) and M_FP-daughter^(global)
    external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree = []
    external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree = []
    
    if len(sorted_rehanged_mutations_all_fpfnratio_across_tree) > 0:
        T_removed_fpfnratio_across_tree, M_removed_fpfnratio_across_tree = remove_mutations_from_tree_and_matrix(
            T_checkpoint_fpfnratio_across_tree, M_checkpoint_fpfnratio_across_tree, 
            sorted_rehanged_mutations_all_fpfnratio_across_tree
        )
        logger_obj.info(f"  Tree after pruning: {M_removed_fpfnratio_across_tree.shape[0]} cells, {M_removed_fpfnratio_across_tree.shape[1]} mutations")
        T_current = copy.deepcopy(T_removed_fpfnratio_across_tree)
        M_current = M_removed_fpfnratio_across_tree.copy()
    
    if len(sorted_fp_mutations_fpfnratio_across_tree) > 0:
        external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_fp_mutations_fpfnratio_across_tree,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    if len(sorted_daughters_to_leaf_mutations_fpfnratio_across_tree) > 0:
        external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_daughters_to_leaf_mutations_fpfnratio_across_tree,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    # Identify M_FN^(global)
    rehanged_fn_mutations_by_fpfnratio_across_tree = df_fp_ratio_and_fn_ratio_fpfnratio_across_tree[
        df_fp_ratio_and_fn_ratio_fpfnratio_across_tree['fn_ratio'] >= fn_ratio_cutoff_across_tree
    ]['identifier'].tolist()
    
    sorted_fn_mutations_fpfnratio_across_tree = [
        i for i in I_attached.columns 
        if i in rehanged_fn_mutations_by_fpfnratio_across_tree
    ]
    
    logger_obj.info(f"  |M_FN^(global)|: {len(rehanged_fn_mutations_by_fpfnratio_across_tree)}")
    
    # Prune and regraft M_FN^(global)
    external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree = []
    
    if len(rehanged_fn_mutations_by_fpfnratio_across_tree) > 0:
        T_removed_fpfnratio_across_tree, M_removed_fpfnratio_across_tree = remove_mutations_from_tree_and_matrix(
            T_checkpoint_fpfnratio_across_tree, M_checkpoint_fpfnratio_across_tree, 
            rehanged_fn_mutations_by_fpfnratio_across_tree
        )
        logger_obj.info(f"  Tree after pruning: {M_removed_fpfnratio_across_tree.shape[0]} cells, {M_removed_fpfnratio_across_tree.shape[1]} mutations")
        T_current = copy.deepcopy(T_removed_fpfnratio_across_tree)
        M_current = M_removed_fpfnratio_across_tree.copy()
    
    if len(sorted_fn_mutations_fpfnratio_across_tree) > 0:
        external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_fn_mutations_fpfnratio_across_tree,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    # Re-attempt any remaining external mutations
    external_mutations_fpfnratio_across_tree = list(set(
        external_mutations_fpfnratio_across_tree_by_sorted_fp_mutations_fpfnratio_across_tree +
        external_mutations_fpfnratio_across_tree_by_sorted_daughters_to_leaf_mutations_fpfnratio_across_tree +
        external_mutations_fpfnratio_across_tree_by_sorted_fn_mutations_fpfnratio_across_tree
    ))
    
    if len(external_mutations_fpfnratio_across_tree) > 0:
        logger_obj.info(f"  Re-attempting {len(external_mutations_fpfnratio_across_tree)} external mutations...")
        sorted_external_mutations_fpfnratio_across_tree = [
            i for i in I_attached.columns 
            if i in external_mutations_fpfnratio_across_tree
        ]
        final_external_mutations_fpfnratio_across_tree, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_external_mutations_fpfnratio_across_tree,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
        logger_obj.info(f"  Final external mutations: {len(final_external_mutations_fpfnratio_across_tree)}")
    
    logger_obj.info("")
    
    # Recompute metrics after pruning
    T_checkpoint_fpfnratio_across_tree_v2 = copy.deepcopy(T_current)
    M_checkpoint_fpfnratio_across_tree_v2 = M_current.copy()
    
    M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2 = M_checkpoint_fpfnratio_across_tree_v2.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_fpfnratio_across_tree_v2 = M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2 = split_merged_columns(M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, mutations_on_T_current_fpfnratio_across_tree_v2)
    
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, fp_mutations_dict_fpfnratio_across_tree_v2 = calculate_fp_fn_ratios_across_tree(
        M_for_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, I_attached
    )
    
    df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_final = pd.merge(
        df_fp_ratio_and_fn_ratio_fpfnratio_across_tree, 
        df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_v2, 
        on='identifier', 
        suffixes=('.1', '.2')
    )
    
    # Merge with Step 7.1 results
    combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree = pd.merge(
        df_fp_ratio_fpratio_within_subclone_final,
        df_fp_ratio_and_fn_ratio_fpfnratio_across_tree_final,
        on='identifier',
        how='left'
    )
    
    # --------------------------------------------------------------------------
    # Step 7.3: Locus-specific FP discordance pruning (M_FP^(locus))
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 7.3: Locus-specific FP discordance pruning")
    logger_obj.info("-" * 80)
    logger_obj.info("  r_FP^(locus)(j) = delta_FP(j) / coverage_j")
    logger_obj.info("")
    logger_obj.info("  M_FP^(locus): mutations flagged by locus-specific FP discordance")
    logger_obj.info("  - These mutations exhibit site-specific bias")
    logger_obj.info("  - May indicate sequencing errors or mapping artifacts")
    logger_obj.info("-" * 80)
    
    # Recalculate backbone nodes
    current_backbone_nodes = get_first_level_backbone_nodes(T_current)
    expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]
    
    mutation_clones_for_persitefp = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)
    
    # Compute locus-specific FP discordance rates
    T_checkpoint_fp_ratio_persitefp = copy.deepcopy(T_current)
    M_checkpoint_fp_ratio_persitefp = M_current.copy()
    
    M_for_fp_ratio_persitefp = M_checkpoint_fp_ratio_persitefp.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_persitefp = M_for_fp_ratio_persitefp.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_persitefp = split_merged_columns(M_for_fp_ratio_persitefp, mutations_on_T_current_persitefp)
    
    df_fp_ratio_persitefp = calculate_fp_ratios_persite_within_subclone(M_for_fp_ratio_persitefp, I_attached, mutation_clones_for_persitefp)
    
    fp_ratio_persite_cutoff = params.get('fp_ratio_persite_cutoff', 0.1)
    
    # Identify M_FP^(locus)
    rehanged_mutations_by_persitefp = df_fp_ratio_persitefp[
        df_fp_ratio_persitefp['fp_ratio_persite'] >= fp_ratio_persite_cutoff
    ]['identifier'].tolist()
    
    rehanged_mutations_by_persitefp_but_backbone = [
        i for i in rehanged_mutations_by_persitefp 
        if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
    ]
    
    logger_obj.info(f"  |M_FP^(locus)|: {len(rehanged_mutations_by_persitefp)}")
    logger_obj.info(f"  |M_FP^(locus) \\ M_backbone|: {len(rehanged_mutations_by_persitefp_but_backbone)}")
    
    sorted_rehanged_mutations_by_persitefp_but_backbone = [
        i for i in I_attached.columns 
        if i in rehanged_mutations_by_persitefp_but_backbone
    ]
    
    # Prune and regraft M_FP^(locus)
    external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone = []
    
    if len(sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
        T_removed_fp_ratio_persitefp, M_removed_fp_ratio_persitefp = remove_mutations_from_tree_and_matrix(
            T_checkpoint_fp_ratio_persitefp, M_checkpoint_fp_ratio_persitefp, 
            sorted_rehanged_mutations_by_persitefp_but_backbone
        )
        logger_obj.info(f"  Tree after pruning: {M_removed_fp_ratio_persitefp.shape[0]} cells, {M_removed_fp_ratio_persitefp.shape[1]} mutations")
        T_current = copy.deepcopy(T_removed_fp_ratio_persitefp)
        M_current = M_removed_fp_ratio_persitefp.copy()
    
    if len(sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
        external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_rehanged_mutations_by_persitefp_but_backbone,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
    
    # Re-attempt any remaining external mutations
    if len(external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone) > 0:
        logger_obj.info(f"  Re-attempting {len(external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone)} external mutations...")
        sorted_external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone = [
            i for i in I_attached.columns 
            if i in external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone
        ]
        final_external_mutations_fp_ratio_persitefp, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
            sorted_attached_mutations=sorted_external_mutations_by_sorted_rehanged_mutations_by_persitefp_but_backbone,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
        all_conflict_mutations.extend(conflict_mutations_temp)
        logger_obj.info(f"  Final external mutations: {len(final_external_mutations_fp_ratio_persitefp)}")
    
    logger_obj.info("")
    
    # Recompute metrics after pruning
    mutation_clones_for_persitefp_v2 = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_current, I_attached)
    
    T_checkpoint_fp_ratio_persitefp_v2 = copy.deepcopy(T_current)
    M_checkpoint_fp_ratio_persitefp_v2 = M_current.copy()
    
    M_for_fp_ratio_persitefp_v2 = M_checkpoint_fp_ratio_persitefp_v2.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_persitefp = M_for_fp_ratio_persitefp_v2.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_persitefp_v2 = split_merged_columns(M_for_fp_ratio_persitefp_v2, mutations_on_T_current_persitefp)
    
    df_fp_ratio_persitefp_v2 = calculate_fp_ratios_persite_within_subclone(M_for_fp_ratio_persitefp_v2, I_attached, mutation_clones_for_persitefp_v2)
    
    df_fp_ratio_persitefp_final = pd.merge(
        df_fp_ratio_persitefp, 
        df_fp_ratio_persitefp_v2, 
        on='identifier', 
        suffixes=('.1', '.2')
    )
    
    final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio = pd.merge(
        combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree,
        df_fp_ratio_persitefp_final,
        on='identifier',
        how='left'
    )
    
    # Save combined metrics to outputpath_full
    final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio.to_csv(
        os.path.join(outputpath_full, "final_combined_df_fp_ratios_within_subclone_and_fpfn_ratios_across_tree_and_persite_fp_ratio.csv"), 
        sep=","
    )
    
    # --------------------------------------------------------------------------
    # Step 7.4: Ancestral retention-based orphaned mutation identification
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 7.4: Ancestral retention-based orphaned mutation identification")
    logger_obj.info("-" * 80)
    logger_obj.info("  eta(j, p(j)) = N_intersection / N_mutant(j)")
    logger_obj.info("")
    logger_obj.info("  M_orphan^(retention): mutations with eta < threshold")
    logger_obj.info("  M_orphan-progeny: daughter mutations of M_orphan^(retention)")
    logger_obj.info("")
    logger_obj.info("These mutations lack sufficient clonal parent support and will be")
    logger_obj.info("relocated to phylogenetically independent branches.")
    logger_obj.info("-" * 80)
    
    # Recalculate backbone nodes
    current_backbone_nodes = get_first_level_backbone_nodes(T_current)
    expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]
    
    # Compute ancestral retention fraction
    T_checkpoint_outgroup = copy.deepcopy(T_current)
    M_checkpoint_outgroup = M_current.copy()
    
    M_for_fp_ratio_and_fn_ratio_outgroup = M_checkpoint_outgroup.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_outgroup = M_for_fp_ratio_and_fn_ratio_outgroup.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_and_fn_ratio_outgroup = split_merged_columns(M_for_fp_ratio_and_fn_ratio_outgroup, mutations_on_T_current_outgroup)
    
    df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation = calculate_intersection_and_inter_vs_fn_flipping_ratio_per_mutation(
        T_checkpoint_outgroup, M_checkpoint_outgroup, I_attached
    )
    
    # Save to outputpath_full
    df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.to_csv(
        os.path.join(outputpath_full, "df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.csv"), 
        sep=","
    )
    
    intersection_vs_fn_flipping_ratio_cutoff = params.get('intersection_vs_fn_flipping_ratio_cutoff', 0.2)
    intersection_cell_count_on_mutation_cutoff = params.get('intersection_cell_count_on_mutation_cutoff', 5)
    intersection_cell_ratio_on_mutation_cutoff = params.get('intersection_cell_ratio_on_mutation_cutoff', 0.2)
    
    # Identify M_orphan^(retention)
    outgroup_mutations = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation[(
        (df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['parent_retention_ratio'] <= intersection_vs_fn_flipping_ratio_cutoff) & 
        ((df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['intersection_cell_count_on_mutation'] <= intersection_cell_count_on_mutation_cutoff) | 
         (df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['intersection_cell_ratio_on_mutation'] <= intersection_cell_ratio_on_mutation_cutoff))
    )]['mutation'].tolist()
    
    outgroup_mutations_but_backbone = [
        i for i in outgroup_mutations 
        if i not in list(set(expanded_mutations_of_current_backbone_nodes))
    ]
    
    logger_obj.info(f"  |M_orphan^(retention)|: {len(outgroup_mutations)}")
    logger_obj.info(f"  |M_orphan^(retention) \\ M_backbone|: {len(outgroup_mutations_but_backbone)}")
    
    # Identify M_orphan-progeny
    sorted_rehanged_mutations_all_outgroup = []
    
    if outgroup_mutations_but_backbone:
        nodes_outgroup_mutations_but_backbone = list(set([
            find_mutation_column(mutation, M_checkpoint_outgroup.columns) 
            for mutation in outgroup_mutations_but_backbone
        ]))
        
        node_dict = {node.name: node for node in T_checkpoint_outgroup.traverse()}
        daughter_nodes_of_outgroup_mutations_but_backbone = {
            mutation: get_all_daughter_mutations(node_dict[mutation])
            for mutation in nodes_outgroup_mutations_but_backbone
            if mutation in node_dict
        }
        daughter_mutations_of_outgroup_mutations_but_backbone = list(set([
            daughter for daughters_list in daughter_nodes_of_outgroup_mutations_but_backbone.values() 
            for daughter in daughters_list
            if daughter not in outgroup_mutations_but_backbone
        ]))
        
        sorted_outgroup_mutations_but_backbone = [
            i for i in I_attached.columns 
            if i in outgroup_mutations_but_backbone
        ]
        sorted_daughter_mutations_of_outgroup_mutations_but_backbone = [
            i for i in I_attached.columns 
            if i in daughter_mutations_of_outgroup_mutations_but_backbone
        ]
        sorted_rehanged_mutations_all_outgroup = (
            sorted_outgroup_mutations_but_backbone + 
            sorted_daughter_mutations_of_outgroup_mutations_but_backbone
        )
        
        logger_obj.info(f"  |M_orphan-progeny|: {len(daughter_mutations_of_outgroup_mutations_but_backbone)}")
    
    # Prune and regraft M_orphan^(retention) and M_orphan-progeny
    external_mutations_by_intersection_vs_fn = []
    
    if len(sorted_rehanged_mutations_all_outgroup) > 0:
        logger_obj.info(f"  Pruning and regrafting {len(sorted_rehanged_mutations_all_outgroup)} orphaned mutations...")
        
        for remove_mut_by_once in sorted_rehanged_mutations_all_outgroup:
            T_removed_outgroup, M_removed_outgroup = remove_mutations_from_tree_and_matrix(
                T_checkpoint_outgroup, M_checkpoint_outgroup, [remove_mut_by_once]
            )
            logger_obj.info(f"    Tree after pruning: {M_removed_outgroup.shape[0]} cells, {M_removed_outgroup.shape[1]} mutations")
            
            T_current = copy.deepcopy(T_removed_outgroup)
            M_current = M_removed_outgroup.copy()
            
            external_mutations_by_intersection_vs_fn_temp, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
                sorted_attached_mutations=[remove_mut_by_once],
                T_current=T_current,
                M_current=M_current,
                I_attached=I_attached,
                P_attached=P_attached,
                ω_NA=omega_NA,
                fnfp_ratio=fnfp_ratio,
                φ=phi,
                logger_obj=logger_obj,
                root_mutations=root_mutations
            )
            all_conflict_mutations.extend(conflict_mutations_temp)
            external_mutations_by_intersection_vs_fn.extend(external_mutations_by_intersection_vs_fn_temp)
    
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 7 Summary
    # --------------------------------------------------------------------------
    logger_obj.info("=" * 80)
    logger_obj.info("Step 7 COMPLETED: Tree refined by discordance-guided pruning")
    logger_obj.info("-" * 80)
    
    logger_obj.info("Refined tree structure:")
    logger_obj.info("-" * 40)
    print_tree_logger(T_current)
    logger_obj.info("-" * 40)
    
    mutations_on_tree = [m for m in M_current.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist() if m != 'ROOT']
    logger_obj.info(f"  Mutations on tree (excluding ROOT): {len(mutations_on_tree)}")
    logger_obj.info(f"  Cells: {M_current.shape[0]}")
    logger_obj.info(f"  Merged columns (including ROOT): {M_current.shape[1]}")
    logger_obj.info(f"  |M_ambiguous|: {len(all_conflict_mutations)}")
    logger_obj.info("=" * 80)
    logger_obj.info("")
    
    # ==========================================================================
    # Compute Omega pre-QC (Checkpoint)
    # ==========================================================================
    logger_obj.info("=" * 80)
    logger_obj.info("CHECKPOINT: Weighted Discordance Index (Omega) computation")
    logger_obj.info("-" * 80)
    logger_obj.info("Computing the weighted discordance index (Omega) to quantify the")
    logger_obj.info("aggregate phylogenetic inconsistency prior to quality control.")
    logger_obj.info("")
    logger_obj.info("  Omega = N_deltaFP + lambda * N_deltaFN,  where lambda = 0.1")
    logger_obj.info("")
    logger_obj.info("  - N_deltaFP: total false positive discordance count")
    logger_obj.info("  - N_deltaFN: total false negative discordance count")
    logger_obj.info("  - lambda: empirical FN/FP discordance weight")
    logger_obj.info("-" * 80)
    
    M_for_omega = M_current.drop(columns=['ROOT'], errors='ignore')
    mutations_on_tree_for_omega = M_for_omega.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_omega = split_merged_columns(M_for_omega, mutations_on_tree_for_omega)
    
    M_for_omega_clean = M_for_omega.loc[:, (M_for_omega != 0).any(axis=0)]
    M_for_omega_clean = M_for_omega_clean.loc[(M_for_omega_clean != 0).any(axis=1)]
    
    I_for_omega = I_attached.loc[M_for_omega_clean.index, M_for_omega_clean.columns].replace({np.nan: 3}).astype(int)
    
    N_deltaFP = ((I_for_omega == 1) & (M_for_omega_clean == 0)).sum().sum()
    N_deltaFN = ((I_for_omega == 0) & (M_for_omega_clean == 1)).sum().sum()
    N_NAto0 = ((I_for_omega == 3) & (M_for_omega_clean == 0)).sum().sum()
    N_NAto1 = ((I_for_omega == 3) & (M_for_omega_clean == 1)).sum().sum()
    
    omega_before_qc = N_deltaFP + params.get('fnfp_ratio', 0.1) * N_deltaFN
    
    logger_obj.info("")
    logger_obj.info("  ┌─────────────────────────────────────────────────────────────────────┐")
    logger_obj.info("  │              WEIGHTED DISCORDANCE INDEX (PRE-QC)                   │")
    logger_obj.info("  ├─────────────────────────────────────────────────────────────────────┤")
    logger_obj.info(f"  │  Weighted Discordance Index (Omega)        : {omega_before_qc:>10.4f}      │")
    logger_obj.info(f"  │    - delta_FP discordance                  : {N_deltaFP:>10}          │")
    logger_obj.info(f"  │    - delta_FN discordance                  : {N_deltaFN:>10}          │")
    logger_obj.info(f"  │    - NA->0 imputations                     : {N_NAto0:>10}          │")
    logger_obj.info(f"  │    - NA->1 imputations                     : {N_NAto1:>10}          │")
    logger_obj.info(f"  │    - FN/FP weight (lambda)                 : {params.get('fnfp_ratio', 0.1):>10.1f}      │")
    logger_obj.info("  ├─────────────────────────────────────────────────────────────────────┤")
    logger_obj.info(f"  │  Cells                                   : {M_for_omega_clean.shape[0]:>10}│")
    logger_obj.info(f"  │  Mutations                               : {M_for_omega_clean.shape[1]:>10}│")
    logger_obj.info("  └─────────────────────────────────────────────────────────────────────┘")
    
    logger_obj.info("")
    
    return T_current, M_current, root_mutations, all_conflict_mutations, omega_before_qc, outgroup_mutations


# ============================================================================
# Step 8: Quality control & filtration (Steps 8.1-8.4)
# ============================================================================

def tree_QC_and_filter(
    T_current: TreeNode,
    M_current: pd.DataFrame,
    scaffold_mutations: List[str],
    I_attached: pd.DataFrame,
    P_attached: pd.DataFrame,
    df_features_new: pd.DataFrame,
    params: Dict[str, Any],
    all_conflict_mutations: List[str],
    root_mutations: List[str],
    outgroup_mutations: List[str],
    remove_artifact_mutations: str,
    outputpath_full: str,
    logger_obj: Optional[logging.Logger] = None,
) -> Tuple[TreeNode, pd.DataFrame, List[str], List[str], List[str], List[str], List[str], List[str]]:
    """
    Step 8: Quality control and filtration.
    
    Identifies and removes:
        8.1: C_orphan - Phylogenetically orphaned cells
        8.2: M_artifact - Artifactual loci (persistently high global FP discordance)
        8.3: C_chimeric - Chimeric cells (doublets, excessive FP discordance)
        8.4: Final tree status summary (pre-output)
    
    Parameters
    ----------
    T_current : TreeNode
        Refined tree after PRG
    M_current : pd.DataFrame
        Refined mutation matrix after PRG
    scaffold_mutations : List[str]
        List of scaffold mutations
    I_attached : pd.DataFrame
        Mutation presence matrix
    df_features_new : pd.DataFrame
        Features matrix
    params : Dict[str, Any]
        Parameters dictionary
    all_conflict_mutations : List[str]
        Existing conflict mutations
    root_mutations : List[str]
        Existing root mutations
    outgroup_mutations : List[str]
        Orphaned mutations from Step 7.4
    remove_artifact_mutations : str
        'yes' or 'no'
    outputpath_full : str
        Output directory path for saving QC files
    logger_obj : logging.Logger, optional
        Logger instance    
    Returns
    -------
    T_current : TreeNode
        Final tree after QC
    M_current : pd.DataFrame
        Final mutation matrix after QC
    to_be_removed_cells : List[str]
        Cells removed as orphans
    identified_doublet_cells : List[str]
        Cells removed as doublets
    to_be_removed_mutations_by_fp_mutations_cross_all_cells : List[str]
        Mutations removed as artifacts
    final_remained_mutations : List[str]
        Final remained mutations
    final_conflict_mutations : List[str]
        Final conflict mutations
    all_conflict_mutations : List[str]
        Updated conflict mutations
    """
    if logger_obj is None:
        logger_obj = logging.getLogger(__name__)
    
    # Extract parameters
    fnfp_ratio = params.get('fnfp_ratio', 0.1)
    phi = params.get('phi', 1.0)
    omega_NA = params.get('general_weight_NA', 0.001)
    
    logger_obj.info("=" * 80)
    logger_obj.info("Step 8: Quality control & filtration")
    logger_obj.info("=" * 80)
    logger_obj.info("")
    
    logger_obj.info(f"  Remove artifact mutations: {remove_artifact_mutations}")
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 8.1: Identification of phylogenetically orphaned cells (C_orphan)
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 8.1: Identification of phylogenetically orphaned cells (C_orphan)")
    logger_obj.info("-" * 80)
    logger_obj.info("Identifying phylogenetically orphaned cells that lack sufficient")
    logger_obj.info("ancestral mutation support.")
    logger_obj.info("")
    logger_obj.info("  Criteria for C_orphan classification:")
    logger_obj.info("    - Intersection count <= 1 (minimal clonal support)")
    logger_obj.info("    - FN discordance >= 1 (persistent false negatives)")
    logger_obj.info("    - NA->1 imputations >= 2 (unreliable imputation)")
    logger_obj.info("-" * 80)
    
    # Recalculate backbone nodes
    current_backbone_nodes = get_first_level_backbone_nodes(T_current)
    expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]
    
    # Compute cell-level metrics
    T_checkpoint_wireless_cells = copy.deepcopy(T_current)
    M_checkpoint_wireless_cells = M_current.copy()
    
    M_for_fp_ratio_and_fn_ratio_wireless_cells = M_checkpoint_wireless_cells.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_wireless_cells = M_for_fp_ratio_and_fn_ratio_wireless_cells.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_fp_ratio_and_fn_ratio_wireless_cells = split_merged_columns(M_for_fp_ratio_and_fn_ratio_wireless_cells, mutations_on_T_current_wireless_cells)
    
    df_intersection_and_flipping_to_1_count_per_cell = calculate_intersection_and_flipping_to_1_count_per_cell(
        M_for_fp_ratio_and_fn_ratio_wireless_cells, I_attached
    )
    
    # Save to outputpath_full
    df_intersection_and_flipping_to_1_count_per_cell.to_csv(
        os.path.join(outputpath_full, "df_intersection_and_flipping_to_1_count_per_cell.csv"), 
        sep=","
    )
    
    # Identify C_orphan
    intersection_count_per_cells_cutoff = params.get('intersection_count_per_cells_cutoff', 1)
    flipping_count_fn_per_cells_cutoff = params.get('flipping_count_fn_per_cells_cutoff', 1)
    flipping_to_1_count_per_cells_cutoff = params.get('flipping_to_1_count_per_cells_cutoff', 2)
    
    df_wireless_cells_filter = df_intersection_and_flipping_to_1_count_per_cell.loc[
        ((df_intersection_and_flipping_to_1_count_per_cell['intersection_count'] == intersection_count_per_cells_cutoff) & 
         (df_intersection_and_flipping_to_1_count_per_cell['flipping_count_fn'] >= flipping_count_fn_per_cells_cutoff) & 
         (df_intersection_and_flipping_to_1_count_per_cell['flipping_to_1_count'] >= flipping_to_1_count_per_cells_cutoff))
    ]
    
    identified_wireless_cells = list(df_wireless_cells_filter.index)
    
    # Identify doublet candidates from orphaned mutations
    df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation = calculate_intersection_and_inter_vs_fn_flipping_ratio_per_mutation(
        T_checkpoint_wireless_cells, M_checkpoint_wireless_cells, I_attached
    )
    
    conflicting_cells_as_doublets_from_parents_format_nested_list = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.loc[
        df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['mutation'].isin(outgroup_mutations), 
        'intersection_cells_on_mutation_parents'
    ]
    conflicting_cells_as_doublets_from_parents_format_flat_list = sum(conflicting_cells_as_doublets_from_parents_format_nested_list.tolist(), [])
    
    conflicting_cells_as_doublets_from_children_format_nested_list = df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation.loc[
        df_intersection_and_inter_vs_fn_flipping_ratio_per_mutation['mutation'].isin(outgroup_mutations), 
        'intersection_cells_on_mutation_children'
    ]
    conflicting_cells_as_doublets_from_children_format_flat_list = sum(conflicting_cells_as_doublets_from_children_format_nested_list.tolist(), [])
    
    # Remove identified cells
    to_be_removed_cells = list(set(
        identified_wireless_cells + 
        conflicting_cells_as_doublets_from_parents_format_flat_list + 
        conflicting_cells_as_doublets_from_children_format_flat_list
    ))
    
    logger_obj.info(f"  |C_orphan| (wireless cells): {len(identified_wireless_cells)}")
    logger_obj.info(f"  |C_doublet_candidates|: {len(conflicting_cells_as_doublets_from_parents_format_flat_list)}")
    logger_obj.info(f"  |Total cells to remove|: {len(to_be_removed_cells)}")
    
    if len(to_be_removed_cells) > 0:
        with open(os.path.join(outputpath_full, "likely_doublet_cells_removed_from_tree_by_fn_flipping.csv"), 'w') as f:
            for cell in to_be_removed_cells:
                f.write(cell + '\n')
        
        M_current = M_current.drop(to_be_removed_cells, errors='ignore')
        logger_obj.info(f"  Removed {len(to_be_removed_cells)} cells from the tree")
    
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 8.2: Artifactual loci (M_artifact) and chimeric cells (C_chimeric)
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 8.2: Artifactual loci and chimeric cell detection")
    logger_obj.info("-" * 80)
    logger_obj.info("Identifying and permanently removing:")
    logger_obj.info("")
    logger_obj.info("  (1) Artifactual loci (M_artifact): mutations with persistently")
    logger_obj.info("      high global FP discordance (r_FP^(global) >= threshold)")
    logger_obj.info("  (2) Chimeric cells (C_chimeric): cells exhibiting excessive")
    logger_obj.info("      cross-mutation FP discordance (>50%)")
    logger_obj.info("")
    logger_obj.info("These represent the most severe data quality issues and are")
    logger_obj.info("permanently excluded from the final phylogeny.")
    logger_obj.info("-" * 80)
    
    # Recalculate backbone nodes
    current_backbone_nodes = get_first_level_backbone_nodes(T_current)
    expanded_mutations_of_current_backbone_nodes = [mutation for node in current_backbone_nodes for mutation in node.split('|')]
    
    to_be_removed_mutations_by_fp_mutations_cross_all_cells = []
    external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells = []
    
    # Compute comprehensive FP metrics
    T_checkpoint_artifact_and_doublet = copy.deepcopy(T_current)
    M_checkpoint_artifact_and_doublet = M_current.copy()
    
    M_for_artifact_and_doublet = M_checkpoint_artifact_and_doublet.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_current_artifact_and_doublet = M_for_artifact_and_doublet.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_for_artifact_and_doublet = split_merged_columns(M_for_artifact_and_doublet, mutations_on_T_current_artifact_and_doublet)
    
    df_fp_ratio_per_mutation_cross_all_cells, df_fp_ratio_per_cell_cross_all_muts, overall_metrics, fp_mutations_dict_cross_all_cells = calculate_comprehensive_fp_metrics(
        M_for_artifact_and_doublet, 
        I_attached.loc[M_checkpoint_artifact_and_doublet.index, M_for_artifact_and_doublet.columns]
    )
    
    # Save to outputpath_full
    df_fp_ratio_per_mutation_cross_all_cells.to_csv(
        os.path.join(outputpath_full, "df_fp_ratio_per_mutation_cross_all_cells.csv"), 
        sep=","
    )
    df_fp_ratio_per_cell_cross_all_muts.to_csv(
        os.path.join(outputpath_full, "df_fp_ratio_per_cell_cross_all_muts.csv"), 
        sep=","
    )
    
    fp_ratio_per_mutation_cross_all_cells_cutoff = params.get('fp_ratio_per_mutation_cross_all_cells_cutoff', 0.2)
    fp_count_per_mutation_cross_all_cells_cutoff = params.get('fp_count_per_mutation_cross_all_cells_cutoff', 5)
    fp_ratio_per_cell_cross_all_muts_cutoff = params.get('fp_ratio_per_cell_cross_all_muts_cutoff', 0.5)
    
    # Identify M_artifact (permanently removed)
    rehanged_fp_mutations_cross_all_cells = df_fp_ratio_per_mutation_cross_all_cells[
        (df_fp_ratio_per_mutation_cross_all_cells['fp_cells_ratio_per_mutation'] >= fp_ratio_per_mutation_cross_all_cells_cutoff) & 
        (df_fp_ratio_per_mutation_cross_all_cells['fp_cells_count'] >= fp_count_per_mutation_cross_all_cells_cutoff)
    ]['identifier'].tolist()
    
    rehanged_fp_mutations_cross_all_cells_but_backbone = [
        i for i in rehanged_fp_mutations_cross_all_cells 
        if i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
    ]
    
    artifact_candidates = rehanged_fp_mutations_cross_all_cells_but_backbone
    
    # Apply removal control
    if remove_artifact_mutations == "yes":
        to_be_removed_mutations_by_fp_mutations_cross_all_cells = artifact_candidates
        
        logger_obj.warning("=" * 80)
        logger_obj.warning(f"  PERMANENTLY REMOVED ARTIFACTUAL LOCI (M_artifact)")
        logger_obj.warning("-" * 80)
        logger_obj.warning(f"  |M_artifact|: {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)} total")
        for mut in to_be_removed_mutations_by_fp_mutations_cross_all_cells:
            logger_obj.warning(f"    - {mut}")
        logger_obj.warning("=" * 80)
        
        pd.Series(to_be_removed_mutations_by_fp_mutations_cross_all_cells).to_csv(
            os.path.join(outputpath_full, "artifact_mutations_permanently_removed.csv"),
            index=False, header=['mutation']
        )
        
        # Remove artifact mutations from the tree
        M_current = M_current.drop(columns=to_be_removed_mutations_by_fp_mutations_cross_all_cells, errors='ignore')
        logger_obj.info(f"  Removed {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)} artifact mutations from the tree")
        
        # Get branch groups and daughter mutations for reattachment
        if len(artifact_candidates) > 0:
            ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
                T_current, artifact_candidates
            )
            
            filtered_fp_mutations_dict_cross_all_cells = {
                mut: other_muts 
                for mut, other_muts in fp_mutations_dict_cross_all_cells.items() 
                if mut in artifact_candidates
            }
            
            nodes_rehanged_fp_mutations_cross_all_cells_but_backbone = list(set([
                find_mutation_column(mutation, M_current.columns) 
                for mutation in artifact_candidates
            ]))
            
            node_dict = {node.name: node for node in T_current.traverse()}
            daughters_nodes_dict_cross_all_cells = {
                mutation: get_all_daughter_mutations(node_dict[mutation])
                for mutation in nodes_rehanged_fp_mutations_cross_all_cells_but_backbone
                if mutation in node_dict
            }
            daughters_mutations_dict_cross_all_cells = list(set([
                daughter for daughters_list in daughters_nodes_dict_cross_all_cells.values() 
                for daughter in daughters_list
                if daughter not in artifact_candidates
            ]))
            
            daughters_to_leaf_mutations_cross_all_cells = []
            fp_mutations_cross_all_cells = []
            
            for branch_mut in ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone:
                daughter_list = ordered_branch_groups_for_rehanged_fp_mutations_cross_all_cells_but_backbone[branch_mut]
                
                daughters_to_leaf_mutations_cross_all_cells = list({
                    item 
                    for key in daughter_list
                    if key in daughters_mutations_dict_cross_all_cells
                    for item in [key] + daughters_mutations_dict_cross_all_cells[key]
                })
                
                fp_mutations_cross_all_cells_init = list({
                    item 
                    for key in daughter_list
                    if key in filtered_fp_mutations_dict_cross_all_cells
                    for item in [key] + filtered_fp_mutations_dict_cross_all_cells[key]
                })
                fp_mutations_cross_all_cells = [
                    i for i in fp_mutations_cross_all_cells_init 
                    if i not in daughter_list 
                    and i not in daughters_to_leaf_mutations_cross_all_cells 
                    and i not in list(set(expanded_mutations_of_current_backbone_nodes + scaffold_mutations))
                ]
            
            daughters_to_leaf_mutations_cross_all_cells_qc = [
                mut for mut in daughters_to_leaf_mutations_cross_all_cells 
                if df_features_new[mut]['mutant_cellnum'] > 5
            ]
            fp_mutations_cross_all_cells_qc = set(
                fp_mutations_cross_all_cells + 
                [i for i in daughters_to_leaf_mutations_cross_all_cells 
                 if i not in daughters_to_leaf_mutations_cross_all_cells_qc]
            )
            
            sorted_fp_mutations_cross_all_cells = [
                i for i in I_attached.columns 
                if i in fp_mutations_cross_all_cells_qc
            ]
            sorted_daughters_to_leaf_mutations_cross_all_cells = [
                i for i in I_attached.columns 
                if i in daughters_to_leaf_mutations_cross_all_cells_qc
            ]
            sorted_rehanged_mutations_all_cross_all_cells = (
                sorted_fp_mutations_cross_all_cells + 
                sorted_daughters_to_leaf_mutations_cross_all_cells
            )
            
            sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells = [
                i for i in sorted_rehanged_mutations_all_cross_all_cells 
                if i not in artifact_candidates
            ]
            remove_mutations_for_rebuild = list(set(
                artifact_candidates + 
                sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells
            ))
            
            # Prune and reattach
            if len(remove_mutations_for_rebuild) > 0:
                T_removed_cross_all_cells, M_removed_cross_all_cells = remove_mutations_from_tree_and_matrix(
                    T_checkpoint_artifact_and_doublet, M_checkpoint_artifact_and_doublet, 
                    remove_mutations_for_rebuild
                )
                logger_obj.info(f"  Tree after pruning: {M_removed_cross_all_cells.shape[0]} cells, {M_removed_cross_all_cells.shape[1]} mutations")
                T_current = copy.deepcopy(T_removed_cross_all_cells)
                M_current = M_removed_cross_all_cells.copy()
            
            if len(sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells) > 0:
                external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells, conflict_mutations_temp, T_current, M_current, root_mutations = attach_mutations_to_current_tree(
                    sorted_attached_mutations=sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells,
                    T_current=T_current,
                    M_current=M_current,
                    I_attached=I_attached,
                    P_attached=P_attached,
                    ω_NA=omega_NA,
                    fnfp_ratio=fnfp_ratio,
                    φ=phi,
                    logger_obj=logger_obj,
                    root_mutations=root_mutations
                )
                all_conflict_mutations.extend(conflict_mutations_temp)
    
    else:
        logger_obj.info(f"  Artifact mutation removal is DISABLED (remove_artifact_mutations=no)")
        logger_obj.info(f"  |M_artifact| candidates identified but NOT removed: {len(artifact_candidates)}")
        
        pd.Series(artifact_candidates).to_csv(
            os.path.join(outputpath_full, "artifact_mutations_identified_but_not_removed.csv"),
            index=False, header=['mutation']
        )
    
    # Identify C_chimeric (chimeric cells / doublets)
    identified_doublet_cells = df_fp_ratio_per_cell_cross_all_muts.loc[
        df_fp_ratio_per_cell_cross_all_muts['fp_muts_ratio_per_cell'] >= fp_ratio_per_cell_cross_all_muts_cutoff
    ]['cell_id'].tolist()
    
    logger_obj.info("")
    logger_obj.info(f"  |C_chimeric| (doublet cells): {len(identified_doublet_cells)}")
    logger_obj.info("  - These cells exhibit excessive cross-mutation FP discordance")
    logger_obj.info("  - Likely represent doublets or technical artifacts")
    
    if len(identified_doublet_cells) > 0:
        with open(os.path.join(outputpath_full, "likely_doublet_cells_removed_from_tree_by_fp_ratio.csv"), 'w') as f:
            for cell in identified_doublet_cells:
                f.write(cell + '\n')
        
        M_current = M_current.drop(identified_doublet_cells, errors='ignore')
        logger_obj.info(f"  Removed {len(identified_doublet_cells)} chimeric cells from the tree")
    
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 8.3: Attach conflict mutations to ROOT
    # --------------------------------------------------------------------------
    logger_obj.info("STEP 8.3: Attach conflict mutations to ROOT")
    logger_obj.info("-" * 80)
    logger_obj.info("Mutations with ambiguous placements (M_ambiguous) are attached")
    logger_obj.info("to the ROOT node as a final resolution step.")
    logger_obj.info("-" * 80)
    
    logger_obj.info(f"  |M_ambiguous| before resolution: {len(all_conflict_mutations)}")
    
    final_remained_mutations = []
    final_conflict_mutations = []
    
    if len(all_conflict_mutations) > 0:
        subtree_groups = cluster_external_mutations_by_intersection(I_attached, all_conflict_mutations)
        logger_obj.info(f"  Processing {len(subtree_groups)} conflict clusters")
        
        final_remained_mutations, final_conflict_mutations, T_current, M_current, root_mutations = process_external_mutations_by_subtree_groups(
            subtree_groups=subtree_groups,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            P_attached=P_attached,
            ω_NA=omega_NA,
            fnfp_ratio=fnfp_ratio,
            φ=phi,
            logger_obj=logger_obj,
            root_mutations=root_mutations
        )
    
    logger_obj.info(f"  |M_ambiguous| after resolution: {len(final_conflict_mutations)}")
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 8.4: Final tree status summary (pre-output)
    # --------------------------------------------------------------------------
    logger_obj.info("=" * 80)
    logger_obj.info("STEP 8.4: Final tree status summary (pre-output)")
    logger_obj.info("=" * 80)
    
    M_final = M_current.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_final = M_final.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_final_split = split_merged_columns(M_final, mutations_on_T_final)
    final_cleaned_M = M_final_split.loc[(M_final_split != 0).any(axis=1)]
    
    logger_obj.info("Final tree status (before output):")
    logger_obj.info(f"  Cells: {final_cleaned_M.shape[0]}")
    logger_obj.info(f"  Mutations: {final_cleaned_M.shape[1]}")
    logger_obj.info(f"  Shape: {final_cleaned_M.shape}")
    logger_obj.info("-" * 80)
    
    # Mutations not on tree
    summary_data = {
        'final_remained_mutations': len(final_remained_mutations),
        'final_conflict_mutations': len(final_conflict_mutations),
        'root_mutations': len(root_mutations),
        'artifact_mutations_permanently_removed': len(to_be_removed_mutations_by_fp_mutations_cross_all_cells),
    }
    
    try:
        external_count = len(set(
            external_mutations_by_intersection_vs_fn +
            external_mutations_cross_all_cells_by_sorted_rehanged_mutations_all_by_fp_mutations_cross_all_cells
        ))
        summary_data['all_external_mutations'] = external_count
    except NameError:
        summary_data['all_external_mutations'] = 0
    
    logger_obj.info("Mutations NOT on tree:")
    total = 0
    for name, count in summary_data.items():
        logger_obj.info(f"  {name}: {count}")
        total += count
    logger_obj.info(f"  Total: {total}")
    logger_obj.info("=" * 80)
    
    summary_file = os.path.join(outputpath_full, "final_tree_status_summary_pre_output.txt")
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("FINAL TREE STATUS SUMMARY (PRE-OUTPUT)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Cells: {final_cleaned_M.shape[0]}\n")
        f.write(f"Mutations: {final_cleaned_M.shape[1]}\n")
        f.write(f"Shape: {final_cleaned_M.shape}\n\n")
        f.write("-" * 60 + "\n")
        f.write("Mutations NOT on tree:\n")
        for name, count in summary_data.items():
            f.write(f"  {name}: {count}\n")
        f.write(f"  Total: {total}\n")
        f.write("=" * 60 + "\n")
    
    logger_obj.info(f"Saved to: {summary_file}")
    logger_obj.info("")
    
    # --------------------------------------------------------------------------
    # Step 8 Summary
    # --------------------------------------------------------------------------
    logger_obj.info("=" * 80)
    logger_obj.info("Step 8 COMPLETED: Quality control and filtration finished")
    logger_obj.info("-" * 80)
    
    logger_obj.info("Final tree structure after QC:")
    logger_obj.info("-" * 40)
    print_tree_logger(T_current)
    logger_obj.info("-" * 40)
    
    mutations_on_tree = [m for m in M_current.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist() if m != 'ROOT']
    logger_obj.info(f"  Mutations on tree (excluding ROOT): {len(mutations_on_tree)}")
    logger_obj.info(f"  Cells: {M_current.shape[0]}")
    logger_obj.info(f"  Merged columns (including ROOT): {M_current.shape[1]}")
    logger_obj.info(f"  |M_artifact| (permanently removed): {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)}")
    logger_obj.info(f"  |C_chimeric| (removed): {len(identified_doublet_cells)}")
    logger_obj.info(f"  |C_orphan| (removed): {len(to_be_removed_cells)}")
    logger_obj.info("=" * 80)
    logger_obj.info("")
    
    return (T_current, M_current, to_be_removed_cells, identified_doublet_cells,
            to_be_removed_mutations_by_fp_mutations_cross_all_cells,
            final_remained_mutations, final_conflict_mutations, all_conflict_mutations)


# ============================================================================
# Helper function: Export phylo results for each CV threshold
# ============================================================================

def export_phylo_results_for_cv(
    T_current: TreeNode,
    M_current: pd.DataFrame,
    I_attached: pd.DataFrame,
    root_mutations: List[str],
    all_conflict_mutations: List[str],
    outputpath_full: str,
    sampleid: str,
    cv_value: float,
    params: Dict[str, Any],
    logger_obj: Optional[logging.Logger] = None,
) -> bool:
    """
    Export phylo results (matrices, tree files, flipping statistics) for a specific CV threshold.
    
    This function is called within the parallel CV search to save complete phylo results
    for each CV value, enabling direct comparison across thresholds.
    
    Parameters
    ----------
    T_current : TreeNode
        Final tree
    M_current : pd.DataFrame
        Final mutation matrix
    I_attached : pd.DataFrame
        Original mutation matrix with NA values
    root_mutations : List[str]
        Root mutations
    all_conflict_mutations : List[str]
        Conflict mutations
    outputpath_full : str
        Base output directory for this CV run
    sampleid : str
        Sample ID
    cv_value : float
        CV threshold value
    params : Dict[str, Any]
        Parameters dictionary
    logger_obj : logging.Logger, optional
        Logger instance
    
    Returns
    -------
    bool : True if successful, False otherwise
    """
    if logger_obj is None:
        logger_obj = logging.getLogger(__name__)
    
    try:
        # Create phylo results directory
        phylo_dir = os.path.join(outputpath_full, "phylo_results")
        os.makedirs(phylo_dir, exist_ok=True)
        
        cv_label = str(cv_value).replace('.', '_')
        
        logger_obj.info(f"  Exporting phylo results for CV={cv_value} to: {phylo_dir}")
        
        # ---- Step 1: Prepare final data ----
        M_current_filtered = M_current.drop(columns=['ROOT'], errors='ignore')
        
        # Add root mutations back
        for mut_on_root in root_mutations:
            if mut_on_root not in M_current_filtered.columns:
                M_current_filtered.insert(0, mut_on_root, 1)
        
        # Expand merged columns
        mutations_on_T_current = M_current_filtered.columns.to_series().apply(
            lambda x: x.split("|")
        ).explode().unique().tolist()
        
        T_full = copy.deepcopy(T_current)
        M_full = split_merged_columns(M_current_filtered, mutations_on_T_current)
        
        # ---- Step 2: Export I matrix with NA=3 ----
        I_full_withNA3 = I_attached.replace({np.nan: 3}).astype(int)
        I_full_withNA3.to_csv(
            os.path.join(phylo_dir, f"I_full_withNA3_CV{cv_label}.txt"),
            sep="\t"
        )
        
        # ---- Step 3: Export M matrix ----
        WriteTfile(
            os.path.join(phylo_dir, f"M_full_basedPivots.filtered_sites_inferred_CV{cv_label}"),
            M_full,
            M_full.index.tolist(),
            M_full.columns.tolist(),
            judge="yes"
        )
        
        # ---- Step 4: Clean and export final matrices ----
        final_cleaned_M_full = M_full.loc[:, (M_full != 0).any(axis=0)]
        final_cleaned_M_full = final_cleaned_M_full.loc[(final_cleaned_M_full != 0).any(axis=1)]
        
        kept_rows = final_cleaned_M_full.index
        kept_cols = final_cleaned_M_full.columns
        
        final_cleaned_I_full_withNA3 = I_full_withNA3.loc[kept_rows, kept_cols]
        
        WriteTfile(
            os.path.join(phylo_dir, f"final_cleaned_M_full_basedPivots.filtered_sites_inferred_CV{cv_label}"),
            final_cleaned_M_full,
            final_cleaned_M_full.index.tolist(),
            final_cleaned_M_full.columns.tolist(),
            judge="yes"
        )
        final_cleaned_I_full_withNA3.to_csv(
            os.path.join(phylo_dir, f"final_cleaned_I_full_withNA3_for_circosPlot_CV{cv_label}.txt"),
            sep="\t"
        )
        
        # ---- Step 5: Identify flipping spots ----
        df_bin_withNA3_for_flipping = final_cleaned_I_full_withNA3.copy()
        df_phylogeny = final_cleaned_M_full.copy()
        
        # Compute flipping spots
        false_negative_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], 
                                          condition_in_bin=0, condition_phylogeny=1)
        )
        false_positive_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], 
                                          condition_in_bin=1, condition_phylogeny=0)
        )
        NAto1_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], 
                                          condition_in_bin=3, condition_phylogeny=1)
        )
        NAto0_flipping_spots = df_bin_withNA3_for_flipping.apply(
            lambda col: find_flipping_spots(col, df_phylogeny[col.name], 
                                          condition_in_bin=3, condition_phylogeny=0)
        )
        
        # Handle empty results
        if false_negative_flipping_spots.empty:
            false_negative_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        if false_positive_flipping_spots.empty:
            false_positive_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        if NAto1_flipping_spots.empty:
            NAto1_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        if NAto0_flipping_spots.empty:
            NAto0_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
        
        # Build flipping spots dataframe
        df_flipping_spots = pd.DataFrame({
            'Mutation': df_bin_withNA3_for_flipping.columns,
            'delta_FN_spots': [', '.join(false_negative_flipping_spots.get(col, [])) 
                             for col in df_bin_withNA3_for_flipping.columns],
            'delta_FP_spots': [', '.join(false_positive_flipping_spots.get(col, [])) 
                             for col in df_bin_withNA3_for_flipping.columns],
            'NA_to_1_spots': [', '.join(NAto1_flipping_spots.get(col, [])) 
                            for col in df_bin_withNA3_for_flipping.columns],
            'NA_to_0_spots': [', '.join(NAto0_flipping_spots.get(col, [])) 
                            for col in df_bin_withNA3_for_flipping.columns]
        })
        df_flipping_spots.to_csv(
            os.path.join(phylo_dir, f"df_flipping_spots_CV{cv_label}.txt"),
            sep="\t", index=False
        )
        
        # ---- Step 6: Calculate total flipping counts ----
        total_FN_flipping = ((df_bin_withNA3_for_flipping == 0) & (df_phylogeny == 1)).sum().sum()
        total_FP_flipping = ((df_bin_withNA3_for_flipping == 1) & (df_phylogeny == 0)).sum().sum()
        total_NAto0 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 0)).sum().sum()
        total_NAto1 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 1)).sum().sum()
        
        fnfp_ratio = params.get('fnfp_ratio', 0.1)
        omega_final = total_FP_flipping + fnfp_ratio * total_FN_flipping
        
        df_total_flipping_count = pd.DataFrame({
            'total_delta_FP': [total_FP_flipping],
            'total_delta_FN': [total_FN_flipping],
            'total_NA_to_0': [total_NAto0],
            'total_NA_to_1': [total_NAto1],
            'weighted_discordance_index_Omega': [omega_final]
        })
        df_total_flipping_count.to_csv(
            os.path.join(phylo_dir, f"df_total_flipping_count_CV{cv_label}.txt"),
            sep="\t", index=False
        )
        
        # Save per-site flipping counts
        df_flip_counts_tree = calculate_flip_counts_per_site(df_bin_withNA3_for_flipping, df_phylogeny)
        df_flip_counts_tree.to_csv(
            os.path.join(phylo_dir, f"df_flipping_count_for_each_mut_CV{cv_label}.txt"),
            sep="\t", index=True
        )
        
        # ---- Step 7: Tree format and clone information ----
        # Export tree as JSON
        tree_dict = tree_to_dict(T_full)
        with open(os.path.join(phylo_dir, f'final_cleaned_tree_node_CV{cv_label}.json'), 'w') as f:
            json.dump(tree_dict, f, indent=4)
        
        # Export tree as text
        T_full.save_to_file(
            os.path.join(phylo_dir, f'final_cleaned_tree_node_CV{cv_label}.txt')
        )
        
        # Assign clone labels to cells
        mutation_clones = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(
            T_full, I_attached
        )
        df_barcode_clones = assign_clone_labels(M_full, mutation_clones)
        df_barcode_clones.to_csv(
            os.path.join(phylo_dir, f"df_barcode_clones_from_phylo_tree_CV{cv_label}.csv"),
            sep=',', index=False
        )
        
        # ---- Step 8: Create summary file ----
        with open(os.path.join(phylo_dir, f"summary_CV{cv_label}.txt"), 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"PhyloSOLID Results Summary - CV={cv_value}\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Sample ID: {sampleid}\n")
            f.write(f"CV threshold: {cv_value}\n")
            f.write(f"Omega (final): {omega_final:.4f}\n")
            f.write(f"Total cells: {final_cleaned_M_full.shape[0]}\n")
            f.write(f"Total mutations: {final_cleaned_M_full.shape[1]}\n")
            f.write(f"Total FP flipping: {total_FP_flipping}\n")
            f.write(f"Total FN flipping: {total_FN_flipping}\n")
            f.write(f"FN/FP weight (lambda): {fnfp_ratio}\n")
            f.write("=" * 80 + "\n")
        
        logger_obj.info(f"  ✓ CV={cv_value} phylo results exported successfully")
        return True
        
    except Exception as e:
        logger_obj.error(f"  ✗ Failed to export phylo results for CV={cv_value}: {e}")
        import traceback
        traceback.print_exc()
        return False


# ============================================================================
# Master function: Orchestrate Steps 6-8
# ============================================================================

def build_fully_resolved_tree(
    T_scaffold: TreeNode,
    M_scaffold: pd.DataFrame,
    scaffold_mutations: List[str],
    I_attached: pd.DataFrame,
    P_attached: pd.DataFrame,
    df_features_new: pd.DataFrame,
    params: Dict[str, Any],
    outputpath_full: str,
    sampleid: str,
    attached_mutations: List[str],
    immune_mutations: List[str] = None,
    spots_to_split: List[str] = None,
    conflict_mutations: List[str] = None,
    remove_artifact_mutations: str = "yes",
    logger_obj: Optional[logging.Logger] = None,
    cv_value: Optional[float] = None,
    export_phylo: bool = False,
) -> Tuple[TreeNode, pd.DataFrame, List[str], List[str], float, List[str], List[str], List[str], List[str], List[str]]:
    """
    Master function orchestrating Steps 6-8 to build a fully-resolved tree.
    
    This function sequentially executes:
        1. build_initial_tree()      - Step 6: Initial tree construction
        2. build_refined_tree()      - Step 7: PRG refinement
        3. tree_QC_and_filter()      - Step 8: QC and filtration
    
    Parameters
    ----------
    T_scaffold : TreeNode
        Scaffold tree
    M_scaffold : pd.DataFrame
        Scaffold mutation matrix
    scaffold_mutations : List[str]
        List of scaffold mutations
    I_attached : pd.DataFrame
        Mutation presence matrix
    P_attached : pd.DataFrame
        Posterior probability matrix
    df_features_new : pd.DataFrame
        Features matrix
    params : Dict[str, Any]
        Parameters dictionary
    outputpath_full : str
        Output directory path
    sampleid : str
        Sample ID
    attached_mutations : List[str]
        Accessory mutations to be placed
    immune_mutations : List[str], optional
        Immune mutations
    spots_to_split : List[str], optional
        Spots to split
    conflict_mutations : List[str], optional
        Conflict mutations
    remove_artifact_mutations : str
        'yes' or 'no'
    logger_obj : logging.Logger, optional
        Logger instance
    cv_value : float, optional
        CV threshold value for file naming
    export_phylo : bool
        Whether to export phylo results for this CV
    
    Returns
    -------
    T_current : TreeNode
        Final refined tree after QC
    M_current : pd.DataFrame
        Final mutation matrix after QC
    root_mutations : List[str]
        Mutations assigned to ROOT
    all_conflict_mutations : List[str]
        Conflict mutations
    omega_before_qc : float
        Omega pre-QC
    to_be_removed_cells : List[str]
        Cells removed as orphans
    identified_doublet_cells : List[str]
        Cells removed as doublets
    to_be_removed_mutations_by_fp_mutations_cross_all_cells : List[str]
        Mutations removed as artifacts
    final_remained_mutations : List[str]
        Final remained mutations
    final_conflict_mutations : List[str]
        Final conflict mutations
    """
    if logger_obj is None:
        logger_obj = logging.getLogger(__name__)
    
    if immune_mutations is None:
        immune_mutations = []
    
    if spots_to_split is None:
        spots_to_split = []
    
    if conflict_mutations is None:
        conflict_mutations = []
    
    all_conflict_mutations = conflict_mutations.copy()
    
    # ============================================================
    # Step 6: Build initial tree
    # ============================================================
    logger_obj.info("=" * 80)
    logger_obj.info("BUILD FULLY RESOLVED TREE: Starting Step 6")
    logger_obj.info("=" * 80)
    
    T_current, M_current, root_mutations, all_conflict_mutations = build_initial_tree(
        T_scaffold=T_scaffold,
        M_scaffold=M_scaffold,
        scaffold_mutations=scaffold_mutations,
        I_attached=I_attached,
        P_attached=P_attached,
        attached_mutations=attached_mutations,
        params=params,
        all_conflict_mutations=all_conflict_mutations,
        logger_obj=logger_obj,
    )
    
    # ============================================================
    # Step 7: PRG refinement
    # ============================================================
    logger_obj.info("=" * 80)
    logger_obj.info("BUILD FULLY RESOLVED TREE: Starting Step 7")
    logger_obj.info("=" * 80)
    
    T_current, M_current, root_mutations, all_conflict_mutations, omega_before_qc, outgroup_mutations = build_refined_tree(
        T_current=T_current,
        M_current=M_current,
        scaffold_mutations=scaffold_mutations,
        I_attached=I_attached,
        P_attached=P_attached,
        df_features_new=df_features_new,
        params=params,
        all_conflict_mutations=all_conflict_mutations,
        root_mutations=root_mutations,
        outputpath_full=outputpath_full,
        logger_obj=logger_obj,
    )
    
    # ============================================================
    # Step 8: QC and filtration
    # ============================================================
    logger_obj.info("=" * 80)
    logger_obj.info("BUILD FULLY RESOLVED TREE: Starting Step 8")
    logger_obj.info("=" * 80)
    
    T_current, M_current, to_be_removed_cells, identified_doublet_cells, \
    to_be_removed_mutations_by_fp_mutations_cross_all_cells, \
    final_remained_mutations, final_conflict_mutations, all_conflict_mutations = tree_QC_and_filter(
        T_current=T_current,
        M_current=M_current,
        scaffold_mutations=scaffold_mutations,
        I_attached=I_attached,
        P_attached=P_attached,
        df_features_new=df_features_new,
        params=params,
        all_conflict_mutations=all_conflict_mutations,
        root_mutations=root_mutations,
        outgroup_mutations=outgroup_mutations,
        remove_artifact_mutations=remove_artifact_mutations,
        outputpath_full=outputpath_full,
        logger_obj=logger_obj,
    )
    
    # ============================================================
    # Final summary
    # ============================================================
    logger_obj.info("=" * 80)
    logger_obj.info("BUILD FULLY RESOLVED TREE: COMPLETED")
    logger_obj.info("-" * 80)
    logger_obj.info(f"  Final cells: {M_current.shape[0]}")
    logger_obj.info(f"  Final mutations: {M_current.shape[1]}")
    logger_obj.info(f"  Omega (pre-QC): {omega_before_qc:.4f}")
    logger_obj.info(f"  Cells removed (orphan): {len(to_be_removed_cells)}")
    logger_obj.info(f"  Cells removed (doublet): {len(identified_doublet_cells)}")
    logger_obj.info(f"  Mutations removed (artifact): {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells)}")
    logger_obj.info("=" * 80)
    logger_obj.info("")
    
    # ---- Export phylo results if requested ----
    if export_phylo and cv_value is not None:
        logger_obj.info("=" * 80)
        logger_obj.info("EXPORTING PHYLO RESULTS FOR THIS CV")
        logger_obj.info("=" * 80)
        
        # Store omega_before_qc in params for the summary
        params_with_omega = params.copy()
        params_with_omega['omega_before_qc'] = omega_before_qc
        
        export_success = export_phylo_results_for_cv(
            T_current=T_current,
            M_current=M_current,
            I_attached=I_attached,
            root_mutations=root_mutations,
            all_conflict_mutations=all_conflict_mutations,
            outputpath_full=outputpath_full,
            sampleid=sampleid,
            cv_value=cv_value,
            params=params_with_omega,
            logger_obj=logger_obj,
        )
        
        if export_success:
            logger_obj.info(f"  ✓ Phylo results exported to: {outputpath_full}/phylo_results/")
        logger_obj.info("=" * 80)
        logger_obj.info("")
    
    return (
        T_current,
        M_current,
        root_mutations,
        all_conflict_mutations,
        omega_before_qc,
        to_be_removed_cells,
        identified_doublet_cells,
        to_be_removed_mutations_by_fp_mutations_cross_all_cells,
        final_remained_mutations,
        final_conflict_mutations,
    )