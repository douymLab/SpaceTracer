import pysam
from collections import defaultdict
import pandas as pd

from SpaceTracer.utils.handle_UMI_combine import handle_seq_type

def get_cb_ub(bam_file,outname,run_type,bins):
    print(bam_file,outname,run_type,bins)
    bam = pysam.AlignmentFile(bam_file, "rb")
    cb_ub_stats = defaultdict(lambda: {'read_count': 0, 'ub_set': set()})

    for read in bam:
        barcode_name, UMI_name=handle_seq_type(read,run_type,bins)

        if barcode_name is not None and UMI_name is not None:
            cb_ub_stats[barcode_name]['read_count'] += 1
            cb_ub_stats[barcode_name]['ub_set'].add(UMI_name)

    bam.close()
    
    df = pd.DataFrame([
            {'barcode': cb, 'nUMI': len(stats["ub_set"]), 'nREAD': stats["read_count"]}
            for cb, stats in cb_ub_stats.items()
        ])

    df.to_csv(outname, sep='\t', index=False)
    return df





