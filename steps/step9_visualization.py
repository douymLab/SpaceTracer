from typing import Dict

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger

model_name = __name__
logger = get_logger(model_name)

class Visualization(BaseStep):
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        inputs={
            'final_vcf':context.get("final_vcf")
        }
        return inputs
        
    def get_outputs(self, context: Dict) -> Dict[str, str]:
        outputs={}
        return outputs

    def _run(self,context):
        inputs=self.get_inputs(context)
        vcf_list=inputs['final_vcf']

        'input_mutation_list':context.get('input_mutation_list','') # this optional file allow the user provide their own mutation list
        
