import pandas as pd

# Work out extraction from a single example

base_dir = "/net/snowwhite/home/beckandy/research/phasing_clean/output/X_noRare/whatshap/"
# f_name = base_dir + "beagle" + "/eval_" + "1" + ".tsv"
# df = pd.read_table(f_name, sep="\t")
# print(df.head())
# print(df.columns)
# switch_flip = df['all_switchflips']
# # switch_flip has number of switches and flips with / delimiter
# # split switch_flip into two numbers
# switch_flip_split = switch_flip.str.split("/", expand=True)
# switches = switch_flip_split[0].astype(int)
# flips = switch_flip_split[1].astype(int)
# total_errors = switches + flips
# n_het = df['covered_variants'][0]

# write function to get dataframe row for each sample
def get_sample_errors(sample_id, method="beagle", base_dir=base_dir):
    f_name = base_dir + method + "/eval_" + str(sample_id) + ".tsv"
    df = pd.read_table(f_name, sep="\t")
    switch_flip = df['all_switchflips']
    switch_flip_split = switch_flip.str.split("/", expand=True)
    switches = switch_flip_split[0].astype(int)
    flips = switch_flip_split[1].astype(int)
    total_errors = switches + flips
    n_het = df['covered_variants'][0]
    return sample_id, method, switches.iloc[0], flips.iloc[0], total_errors.iloc[0], n_het

get_sample_errors(1, method="beagle")

# Initialize empty list to hold results
results = []

for sample_id in range(1, 1001):
    for method in ["beagle", "eagle", "shapeit"]:
        result = get_sample_errors(sample_id, method=method)
        results.append(result)

df_final = pd.DataFrame(results, columns=["sample_id", "method", "switches", "flips", "total_errors", "n_het"])
df_final.to_csv(base_dir + "X_noRare_errors_summary.tsv", sep="\t", index=False)