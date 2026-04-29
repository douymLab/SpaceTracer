#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from SpaceTracer.cores.phasing import PhaseConfig,  run_phase_mode
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger

model_name = __name__
logger = get_logger(model_name)

class PhasingCandidateStep(BaseStep):

    def get_inputs(self, context):
        inputs = {
            "in_filter_bam": context.get("in_filter_bam"),
            "merged_germline_file": context.get("merged_germline_file"),
            "merged_ind_geno_filter_file": context.get("merged_ind_geno_filter_file")
        }

        if context.get("gtexGene"):
            inputs["gtexGene"] = context.get("gtexGene")
        if context.get("gencode"):
            inputs["gencode"] = context.get("gencode")

        return inputs

    def get_outputs(self, context):
        return {
            "phasing_result": os.path.join(self.work_dir,"phasing_results.txt"),
            "cluster_event_result": os.path.join(self.work_dir,"cluster_events.txt")
        }

    def get_step_config(self):
        return self.config.get("steps", {}).get("phasing", {})

    def _run(self, context):
        inputs=self.get_inputs(context)
        bam=inputs["in_filter_bam"]
        merged_germline_file=inputs["merged_germline_file"]
        merged_ind_geno_filter_file=inputs["merged_ind_geno_filter_file"]

        seq_type = self.config.get("sequence_type")
        fasta_file = self.config.get("genome_fasta")
        genome_details=self.config['genome_details']
        species=genome_details['species']
        gene_bed=self.config['gene_bed']
        bin_size=1 ## treat as bin1 level, if stereo-seq

        step_config=self.get_step_config()
        minprior=float(step_config["minprior"])
        alpha=float(step_config["alpha"])
        min_dp=int(step_config["min_dp"])
        min_total_dp=int(step_config["min_total_dp"])
        phasing_pad=int(step_config["phasing_pad"])
        merge_gap=int(step_config["merge_gap"])
        max_target=int(step_config["max_target"])
        seed=int(step_config["seed"])
        max_dist=10 # fixed

        autosomes=genome_details['chromosomes']['autosomes']

        outputs=self.get_outputs(context)
        out_phasing_file=outputs["phasing_result"]
        out_cluster_file=outputs["cluster_event_result"]

        phasing_chromosomes=autosomes
        thread = self.threads

        phase_config = PhaseConfig(
            fasta=fasta_file,
            bam=bam,
            germline=merged_germline_file,
            indgeno=merged_ind_geno_filter_file,
            seq_type=seq_type,
            bin_size=bin_size,
            minprior=minprior,
            phasing_chromosomes=phasing_chromosomes,
            thread=thread,
            species=species,
            gene_bed=gene_bed,
            min_dp=min_dp,
            min_total_dp=min_total_dp,
            out_phasing_file=out_phasing_file,
            out_cluster_file=out_cluster_file,
            alpha=alpha,
            max_dist=max_dist,
            phasing_pad=phasing_pad,
            merge_gap=merge_gap,
            max_target=max_target,
            seed=seed
        )
        outfile = run_phase_mode(phase_config)

        return outfile