import gzip
from pathlib import Path
import os
import pandas as pd
from SpaceTracer.utils.logger import get_logger
model_name=__name__
logger = get_logger(model_name)


class mappabilityFeatures:
    def __init__(self,bedgraph_path:str):
        self.bedgraph_path = bedgraph_path
        self.mappability_path=self._prepare_mappability()

    def _prepare_mappability(self):
        bedgraph_path = self.bedgraph_path
        input_path = Path(bedgraph_path)
        if not input_path.exists():
            raise FileNotFoundError(f"Mappability file not found: {bedgraph_path}")
        
        input_dir = input_path.parent
        base_name = input_path.stem
        if base_name.endswith('.bedgraph'):
            base_name = base_name.replace('.bedgraph', '')
        elif base_name.endswith('.bedgraph.gz'):
            base_name = base_name.replace('.bedgraph.gz', '')
        
        output_dir = input_dir / f'{base_name}_parquet'
        
        complete_marker = output_dir / '.complete'
        if output_dir.exists() and complete_marker.exists():
            parquet_files = list(output_dir.glob('*.parquet'))
            if len(parquet_files) > 0:
                return str(output_dir)
            else:
                logger.warning(f"Found .complete marker but no parquet files, will reprocess")
        
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            test_file = output_dir / '.write_test'
            test_file.touch()
            test_file.unlink()
        except (PermissionError, OSError) as e:
            raise PermissionError(
                f"No write permission in directory: {input_dir}"
                f"Error: {e}"
            )
        
        logger.info(f"Converting mappability file to parquet: {bedgraph_path}. This step will process only once.")
        
        is_gzipped = bedgraph_path.endswith('.gz')
        opener = gzip.open if is_gzipped else open
        open_mode = 'rt' if is_gzipped else 'r'
        
        CHUNK_SIZE = 1_000_000
        chrom_buffers = {}
        chrom_chunk_counts = {}  
        total_lines = 0
        
        try:
            with opener(bedgraph_path, open_mode) as f:
                for line in f:
                    if line.startswith(('track', '#', 'browser')) or not line.strip():
                        continue
                    
                    parts = line.strip().split()
                    if len(parts) != 4:
                        continue
                    
                    try:
                        chrom, start, end, score = parts
                        start = int(start)
                        end = int(end)
                        score = float(score)
                    except ValueError:
                        continue
                    
                    if chrom not in chrom_buffers:
                        chrom_buffers[chrom] = []
                        chrom_chunk_counts[chrom] = 0
                    
                    chrom_buffers[chrom].append({
                        'start': start,
                        'end': end,
                        'mappability': score
                    })
                    
                    total_lines += 1
                    
                    if len(chrom_buffers[chrom]) >= CHUNK_SIZE:
                        self._write_chunk(chrom, chrom_buffers[chrom], output_dir, chrom_chunk_counts[chrom])
                        chrom_chunk_counts[chrom] += 1
                        chrom_buffers[chrom] = []  # 清空缓冲区
                    
                    if total_lines % 10_000_000 == 0:
                        logger.debug(f"Processed {total_lines:,} lines, {len(chrom_buffers)} chromosomes...")
            
            for chrom, buffer in chrom_buffers.items():
                if buffer:
                    self._write_chunk(chrom, buffer, output_dir, chrom_chunk_counts[chrom])
                    chrom_chunk_counts[chrom] += 1
            
            saved_chroms = self._merge_chunks(chrom_chunk_counts, output_dir)
        
        except Exception as e:
            raise RuntimeError(f"Error processing mappability file: {e}")
        
        if total_lines == 0:
            raise ValueError(f"No valid data found in mappability file: {bedgraph_path}")
        
        if len(saved_chroms) == 0:
            raise RuntimeError(f"Failed to save any chromosome data to {output_dir}")
        
        try:
            with open(complete_marker, 'w') as f:
                f.write(f"Conversion completed at: {pd.Timestamp.now()}")
                f.write(f"Source file: {bedgraph_path}")
                f.write(f"Chromosomes: {', '.join(sorted(saved_chroms))}")
        except Exception as e:
            logger.warning(f"Could not create completion marker: {e}")
        
        logger.info(f"Mappability conversion complete: {output_dir}")
        logger.info(f"Saved {len(saved_chroms)} chromosome files")
        
        return str(output_dir)

    def _write_chunk(self, chrom, buffer, output_dir, chunk_id):
        if not buffer:
            return
        
        df = pd.DataFrame(buffer)
        chunk_file = output_dir / f'{chrom}.chunk_{chunk_id}.parquet'
        df.to_parquet(chunk_file, engine='pyarrow', compression='snappy', index=False)

    def _merge_chunks(self, chrom_chunk_counts, output_dir):
        saved_chroms = []
        
        for chrom, num_chunks in chrom_chunk_counts.items():
            if num_chunks == 0:
                continue
            
            dfs = []
            for i in range(num_chunks):
                chunk_file = output_dir / f'{chrom}.chunk_{i}.parquet'
                dfs.append(pd.read_parquet(chunk_file))
            
            df = pd.concat(dfs, ignore_index=True)
            df = df.sort_values('start').reset_index(drop=True)
            
            final_file = output_dir / f'{chrom}.parquet'
            df.to_parquet(final_file, engine='pyarrow', compression='snappy', index=False)
            
            for i in range(num_chunks):
                chunk_file = output_dir / f'{chrom}.chunk_{i}.parquet'
                chunk_file.unlink()
            
            saved_chroms.append(chrom)
            logger.info(f"Saved {chrom}: {len(df):,} intervals, "
                    f"mappability range: [{df['mappability'].min():.3f}, {df['mappability'].max():.3f}]")
        
        return saved_chroms
    
    def _load_mappability_for_chrom(self, chrom: str) -> pd.DataFrame:
        """load mappability data for one chromosome"""
        mappability_dir=self.mappability_path
        chrom_file = os.path.join(mappability_dir, f'{chrom}.parquet')
        
        if not os.path.exists(chrom_file):
            logger.warning(f"No mappability data for {chrom}")
            return pd.DataFrame(columns=['start', 'end', 'mappability'])
        
        return pd.read_parquet(chrom_file)
    
    def _query_mappability(self, mappability_df: pd.DataFrame, positions: list) -> list:
        """query mappability values for a list of positions"""
        import numpy as np
        
        if len(mappability_df) == 0:
            return [None] * len(positions)
        
        positions = np.array(positions)
        starts = mappability_df['start'].values
        ends = mappability_df['end'].values
        scores = mappability_df['mappability'].values
        
        results = np.zeros(len(positions))
        
        for i, pos in enumerate(positions):
            # searchsorted can find the index of positions quickly
            idx = np.searchsorted(starts, pos, side='right') - 1
            
            if idx >= 0 and idx < len(starts) and starts[idx] <= pos < ends[idx]:
                results[i] = scores[idx]
        
        return results.tolist()
    
  