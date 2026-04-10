import sqlite3


class ChunkIndexDB:
    def __init__(self, db_path):
        self.db_path = db_path

    def _connect(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def list_chunks(self):
        with self._connect() as conn:
            rows = conn.execute("""
                SELECT chunk_id
                FROM chunks
                ORDER BY chrom, chunk_idx
            """).fetchall()
        return [row["chunk_id"] for row in rows]

    def list_chroms(self):
        with self._connect() as conn:
            rows = conn.execute("""
                SELECT chrom
                FROM chroms
                ORDER BY chrom
            """).fetchall()
        return [row["chrom"] for row in rows]

    def get_chunk(self, chunk_id):
        with self._connect() as conn:
            row = conn.execute("""
                SELECT *
                FROM chunks
                WHERE chunk_id = ?
            """, (chunk_id,)).fetchone()
        return dict(row) if row else None

    def get_chunks_by_chrom(self, chrom):
        with self._connect() as conn:
            rows = conn.execute("""
                SELECT *
                FROM chunks
                WHERE chrom = ?
                ORDER BY chunk_idx
            """, (chrom,)).fetchall()
        return [dict(row) for row in rows]

    def get_chrom(self, chrom):
        with self._connect() as conn:
            row = conn.execute("""
                SELECT *
                FROM chroms
                WHERE chrom = ?
            """, (chrom,)).fetchone()
        return dict(row) if row else None