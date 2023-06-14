import pandas as pd

def load_data(paths_to_counts: str) -> list[str]:
    return [pd.read_csv(file, index_col=[0], sep='\t') for file in paths_to_counts]

def trim_data(list_of_dfs: list[pd.DataFrame]) -> list[pd.DataFrame]:
    return [df[df.columns[[-1]]] for df in list_of_dfs]

def capitalize_data(list_of_dfs: list[pd.DataFrame]) -> list[pd.DataFrame]:
    for i, df in enumerate(list_of_dfs):
        list_of_dfs[i].index = df.index.str.upper()
        list_of_dfs[i] = make_multi_index(df)
    return list_of_dfs

def join_data(list_of_dfs: list[pd.DataFrame]) -> pd.DataFrame:
    return pd.concat(list_of_dfs, join="inner", axis=1)

def get_experiment(column_header: str) -> str:
    exp_name = column_header.split('/')[-3]
    return exp_name

def make_multi_index(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = pd.MultiIndex.from_tuples([(get_experiment(df.columns[0]), df.columns[0])], names=[None, df.index.name])
    df.index.name = None
    return df

#df.to_csv(snakemake.output["merged_counts"], sep='\t', index=True)

files = snakemake.input["counts"]
#files = ["../../NMR_iPS/counts/SRR5985204.counts.txt", "../../NMR_iPS/counts/SRR5985215.counts.txt"]
dfs = load_data(files)
dfs = trim_data(dfs)
dfs = capitalize_data(dfs)
df = join_data(dfs)
file_out = snakemake.output["merged_counts"]
df.sort_index(ascending=True, inplace=True)
df.to_csv(file_out, sep='\t', index=True)
