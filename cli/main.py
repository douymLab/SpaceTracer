#!/usr/bin/env python3
"""
Mutation Caller - Main CLI Entry Point
"""

import click
from SpaceTracer import __version__
from SpaceTracer.cli import run, split, call_batch, merge
from SpaceTracer.utils.logger import setup_logger

@click.group()
@click.version_option(version=__version__)
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.option('--debug', is_flag=True, help='Enable debug mode')
@click.pass_context
def cli(ctx, verbose, debug):
    """
    Mutation Caller - A tool for calling mutations from BAM files
    
    Examples:
    
        # Run full pipeline
        SpaceTracer run --bam sample.bam --genome hg38 --output results/
        
        # Split large regions
        SpaceTracer split --regions large.bed --chunk-size 50000
        
        # Batch processing
        SpaceTracer call-batch --manifest regions/manifest.json
        
        # Merge results
        SpaceTracer merge --input results/ --output final.vcf
    """
    # 设置日志
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['debug'] = debug
    
    log_level = 'DEBUG' if debug else ('INFO' if verbose else 'WARNING')
    setup_logger(log_level)

# 注册子命令
cli.add_command(run.run)
cli.add_command(call_batch.call_batch)
cli.add_command(merge.merge)

def main():
    """Main entry point"""
    cli(obj={})

if __name__ == '__main__':
    main()
