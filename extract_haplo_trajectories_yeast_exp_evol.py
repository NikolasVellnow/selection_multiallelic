import pandas as pd
import matplotlib.pyplot as plt



haplotypes = ["DBVPG6044", "DBVPG6765", "Y12", "YPS128"]
haplo_files = ["Philips_et_al_data/DBVPG6044_hap_freqs_10kb.txt",
               "Philips_et_al_data/DBVPG6765_hap_freqs_10kb.txt",
               "Philips_et_al_data/Y12_hap_freqs_10kb.txt",
               "Philips_et_al_data/YPS128_hap_freqs_10kb.txt"]



# "category 1" SNP
#focal_rep = "3"
#focal_chrom = "C12"
#focal_pos = 436424

# "category 2" SNP
#focal_rep = "3"
#focal_chrom = "C16"
#focal_pos = 457847

# "category 3" SNP
focal_rep = "3"
#focal_chrom = "C16"
#focal_pos = 451847

main_df = pd.DataFrame()

df_start = pd.read_csv(haplo_files[0], sep=" ")
df_start = df_start.iloc[:, [0, 1]]
df_start['chr'] = df_start['chr'].astype(str)

#print(df_start.head())

#print("columns: ", df_start.columns)
#print("columns: ", df_start.dtypes)


#print(df_start['pos'])

for row in df_start.itertuples():
    print(row.chr)
    print(row.pos)
    focal_chrom = row.chr
    focal_pos = row.pos

    # loop for one 10 KB window / region; loops over four haplotypes
    for idx, haplotype in enumerate(haplotypes):
        df = pd.read_csv(haplo_files[idx], sep=" ")
        #print(df.info())
        #print(df.describe())
        #print(df.head())

        focal_row = df[(df["chr"] == focal_chrom) & (df["pos"] == focal_pos)]
        #print(focal_row)
        col_names = list(df.columns)
        dict_list = []      # Here we will save the data all timepoints of one haplotype

        for col_name in col_names:
        
            value_dict = {}
        
            splitted_col_name = col_name.split("_")

            if col_name == "chr":
                chrom = focal_row[col_name].values[0]
                continue

            if col_name == "pos":
                pos = focal_row[col_name].values[0]
                continue

            if col_name == "hap_anc":
                value_dict["rep"] = focal_rep
                value_dict["chr"] = chrom
                value_dict["pos"] = pos
                value_dict["time_point"] = int(0)
                value_dict[f"{haplotype}_freq"] = focal_row[col_name].values[0]
                dict_list.append(value_dict)
        
            # only take values for focal replicate
            if splitted_col_name[1] == focal_rep:
                #print(col_name)
                value_dict["rep"] = focal_rep
                value_dict["chr"] = chrom
                value_dict["pos"] = pos
                value_dict["time_point"] = int(splitted_col_name[2])
                value_dict[f"{haplotype}_freq"] = focal_row[col_name].values[0]
                dict_list.append(value_dict)

        #print("dict list")
        #print(dict_list)

        new_df = pd.DataFrame(dict_list)
        #print(new_df)

        if idx == 0:
            haps_df = new_df
        else:
            haps_df = pd.merge(haps_df, new_df, on=["rep", "chr", "pos", "time_point"])

    haps_df["sum"] = haps_df[["DBVPG6044_freq","DBVPG6765_freq","Y12_freq","YPS128_freq"]].sum(axis=1)
    haps_df["gen"] = haps_df["time_point"]*30
    #print(haps_df)
    main_df = pd.concat([main_df, haps_df], ignore_index=True)



main_df.to_csv(f"haplo_trajectories_rep_{focal_rep}.csv", index = False)



