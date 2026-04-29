import pandas as pd
from dataclasses import dataclass
from typing import Callable, List

@dataclass
class FilterStep:
    name: str
    func: Callable[[pd.DataFrame], pd.Series]  # return bool Series，True means pass filtration
    description: str = ""

class FilterPipeline:
    def __init__(self):
        self.steps: List[FilterStep] = []
    
    def add_filter(self, name: str, func: Callable, description: str = ""):
        self.steps.append(FilterStep(name=name, func=func, description=description))
        return self 
    
    def run(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        do filtration and mark each line
        
        Returns:
            df, with new cols:
            - filter_{name}: bool, True means pass filtration
            - filter_reason: filter reason, seperated by ";"
            - passed: bool, does this line pass all filtration
        """
        result_df = df.copy()
        
        filter_cols = []
        for step in self.steps:
            col = f"filter_{step.name}"
            result_df[col] = step.func(result_df)
            filter_cols.append(col)
        
        def get_filter_reason(row):
            failed = [
                step.name 
                for step, col in zip(self.steps, filter_cols) 
                if not row[col]
            ]
            return ';'.join(failed) if failed else 'PASS'
        
        result_df['filter_reason'] = result_df.apply(get_filter_reason, axis=1)
        result_df['passed'] = result_df['filter_reason'] == 'PASS'
        
        return result_df
    
    def summary(self, df: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for step in self.steps:
            col = f"filter_{step.name}"
            if col in df.columns:
                n_fail = (~df[col]).sum()
                n_pass = df[col].sum()
                rows.append({
                    'filter': step.name,
                    'description': step.description,
                    'n_pass': n_pass,
                    'n_fail': n_fail,
                    'pass_rate': f"{n_pass / len(df) * 100:.1f}%"
                })
        return pd.DataFrame(rows)


