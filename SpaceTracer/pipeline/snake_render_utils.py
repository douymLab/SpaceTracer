#!/usr/bin/env python3
# not stable
from typing import Dict, List


def quote_path(path: str) -> str:
    return f'"{path}"'


def indent(text: str, n: int = 4) -> str:
    prefix = " " * n
    return "\n".join(prefix + line if line else line for line in text.splitlines())


def wildcard_flags_for_level(level: str) -> List[str]:
    if level == "sample":
        return ["--sample {wildcards.sample}"]
    elif level == "chrom":
        return ["--sample {wildcards.sample}", "--chrom {wildcards.chrom}"]
    elif level == "chunk":
        return ["--sample {wildcards.sample}", "--chunk {wildcards.chunk}"]
    else:
        raise ValueError(f"Unknown level: {level}")


def log_pattern(step_name: str, run_level: str) -> str:
    if run_level == "sample":
        return f'logs/{step_name}/{{sample}}.log'
    elif run_level == "chrom":
        return f'logs/{step_name}/{{sample}}_{{chrom}}.log'
    elif run_level == "chunk":
        return f'logs/{step_name}/{{sample}}_{{chunk}}.log'
    else:
        raise ValueError(f"Unknown run_level: {run_level}")


def flatten_output_dict(outputs: Dict[str, str]) -> List[str]:
  
    result = []
    for _, value in outputs.items():
        if isinstance(value, str):
            result.append(value)
    return result
