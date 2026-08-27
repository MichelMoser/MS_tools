import pyarrow.parquet as pq
import sys

pf = pq.ParquetFile(sys.argv[1])
print(pf.metadata.num_rows)
