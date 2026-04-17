import os

from SpaceTracer.cores.phylosolid_runner import run_phylosolid_with_raw_args
from SpaceTracer.utils.phylosolid_wrapper import print_phylosolid_debug_info

from SpaceTracer.steps.base import BaseStep

class Phylogeny(BaseStep):
    def get_inputs(self,context):
        inputs={
            'final_vcf': context.get('final_vcf')
        }
        return inputs

    def get_outputs(self, context):
        sample_name="Sample"
        vcf_output_file = os.path.join(self.step_dir, sample_name + "_total_pred_truesites.vcf")
        vcf_pass_output_file = os.path.join(self.step_dir, sample_name + "_total_pred_truesites_PASS.vcf")
        outputs={
            'raw_pred_vcf':vcf_output_file,
            'final_vcf':vcf_pass_output_file
        }
        return outputs

def main():
    print_phylosolid_debug_info()

    phylosolid --workdir ./results scrna \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt

    result = run_phylosolid_with_raw_args(
        raw_args=["-h"
        ],
        capture_output=True,
        dry_run=False
    )

    print("return code: {}".format(result.returncode))

    if result.stdout:
        print("===== STDOUT =====")
        print(result.stdout)

    if result.stderr:
        print("===== STDERR =====")
        print(result.stderr)


if __name__ == "__main__":
    main()
