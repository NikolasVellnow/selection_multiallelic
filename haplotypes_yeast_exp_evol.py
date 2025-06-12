import pandas as pd
import matplotlib.pyplot as plt



haplotypes = ["DBVPG6044", "DBVPG6765", "Y12", "YPS128"]
haplo_files = ["Philips_et_al_data/DBVPG6044_hap_freqs_10kb.txt",
               "Philips_et_al_data/DBVPG6765_hap_freqs_10kb.txt",
               "Philips_et_al_data/Y12_hap_freqs_10kb.txt",
               "Philips_et_al_data/YPS128_hap_freqs_10kb.txt"]


costs = [0.0007, 0.0001, 0.0079, 0.0082, 0.0001, 0.0027, 0.0008, 0.0015, 0.0068, 0.0076, 0.0020, 0.0012, 0.0008, 0.0017, 0.0008, 0.0002]
norm_costs = [cost / 0.0429 for cost in costs]

# "category 1" SNP
#focal_rep = "3"
#focal_chrom = "C12"
#focal_pos = 436424

# "category 2" SNP
focal_rep = "3"
focal_chrom = "C16"
focal_pos = 457847

# "category 3" SNP
#focal_rep = "3"
#focal_chrom = "C16"
#focal_pos = 451847




for idx, haplotype in enumerate(haplotypes):
    df = pd.read_csv(haplo_files[idx], sep=" ")
    #print(df.info())
    #print(df.describe())
    #print(df.head())

    focal_row = df[(df["chr"] == focal_chrom) & (df["pos"] == focal_pos)]
    print(focal_row)
    col_names = list(df.columns)
    dict_list = []

    for col_name in col_names:
        
        value_dict = {}
        splitted_col_name = col_name.split("_")
        # only take values for focal replicate
        if col_name in ("chr", "pos"):
            continue

        if col_name == "hap_anc":
            value_dict["time_point"] = int(0)
            value_dict[f"{haplotype}_freq"] = focal_row[col_name].values[0]
            dict_list.append(value_dict)
        
        if splitted_col_name[1] == focal_rep:
            #print(col_name)
            value_dict["time_point"] = int(splitted_col_name[2])
            value_dict[f"{haplotype}_freq"] = focal_row[col_name].values[0]
            dict_list.append(value_dict)

    #print(dict_list)

    new_df = pd.DataFrame(dict_list)

    if idx == 0:
        haps_df = new_df
    else:
        haps_df = pd.merge(haps_df, new_df, on="time_point")



haps_df["sum"] = haps_df[["DBVPG6044_freq","DBVPG6765_freq","Y12_freq","YPS128_freq"]].sum(axis=1)
haps_df["gen"] = haps_df["time_point"]*30
print(haps_df)

haps_df.to_csv(f"hap_freqs_rep_{focal_rep}_chr_{focal_chrom}_pos_{focal_pos}.csv", index = False)

# Get a 'viridis' colormap with as many colors as there are lines
cmap = plt.get_cmap('viridis', len(haplotypes))

print(cmap(0))
print(cmap(1))
print(cmap(2))
print(cmap(3))

plt.plot(haps_df["gen"], haps_df["DBVPG6044_freq"], color = cmap(0))
plt.plot(haps_df["gen"], haps_df["DBVPG6765_freq"], color = cmap(1))
plt.plot(haps_df["gen"], haps_df["Y12_freq"], color = cmap(2))
plt.plot(haps_df["gen"], haps_df["YPS128_freq"], color = cmap(3))
plt.plot(haps_df["gen"][1:len(haps_df["gen"])], norm_costs, color = "black")
plt.ylim((0.0, 1.0))


plt.xlabel("Generation")
plt.ylabel("Haplotype frequency")
plt.title(f"Replicate: {focal_rep}, Chromosome: {focal_chrom}, position: {focal_pos}")
plt.show()



