import pyarrow as pa
import pyarrow.parquet as pq
import sys

test = pq.read_table(sys.argv[1])

for e in test: 
    print(e)



    
