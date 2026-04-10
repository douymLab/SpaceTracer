#!/usr/bin/env python3

from typing import Dict


def build_template_context(config: Dict) -> Dict:

    return {
        "config": config,
        "sample": "{sample}",
        "chrom": "{chrom}",
        "chunk": "{chunk}",
        "chunk_index": "{chunk}",
        "genome_details": config.get("genome_details"),
    }
