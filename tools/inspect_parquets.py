import polars as pl

paths = [
    'Results/Covid/node_mutations.parquet',
    'Results/Influenza/node_mutations.parquet',
]

for p in paths:
    try:
        df = pl.read_parquet(p)
        print('PATH:', p)
        print('COLUMNS:', df.columns)
        print('DTYPES:', [str(t) for t in df.dtypes])
        print('EXAMPLE ROWS:')
        for r in df.head(5).to_dicts():
            print(r)
        print('\n')
    except Exception as e:
        print('ERROR reading', p, e)
