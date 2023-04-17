# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 10:30:25 2023

@author: andrea.graziadei
"""

import os
from pyteomics import mzml, mass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import glob
from Bio import SeqIO

protein_name = "sp|P75591|NUSA_MYCPN"
fasta_file = glob.glob("*.fasta")[0]
peptides = pd.read_csv("peptides.txt", sep="\t")
cutoff_value = 0.66
missed_cleavages_allowed = 1


conditions = ['G1B1', 'G1B2', 'G1B3', 'G1B4',
       'G2B1', 'G2B2', 'G2B3', 'G2B4',
       'G3B1', 'G3B2', 'G3B3', 'G3B4']



# functions

def FilterPeptides(df, condition, cutoff):
    column_name = str("Intensity " +condition)
    df_filtered = df[df[column_name]>0]
    
    cutoff_number = int(round(cutoff * len(df_filtered)))
    
    df_filtered = df_filtered.nlargest(cutoff_number, column_name, keep="all")
    df_melt = pd.melt(df_filtered,
                      id_vars=[column_name, "y_value"],
                      value_vars=["Start position", "End position"])
    df_melt["Log2 Int"] = np.log2(df_melt[column_name])

    return df_melt
    

# get protein of interest
sequence_database = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
protein = sequence_database.get(protein_name)
sequence_length = len(protein)

# filter peptides
peptides = peptides[peptides["Reverse"]!="+"]
peptides = peptides[peptides["Potential contaminant"]!="+"]

peptides_filtered = peptides[peptides["Missed cleavages"]<=missed_cleavages_allowed]

peptides_filtered = peptides.filter(regex="Intensity .*")
peptides_filtered["Start position"] = peptides["Start position"]
peptides_filtered["End position"] = peptides["End position"]


peptides_sorted = peptides_filtered.sort_values(by="Start position", ignore_index=True)
peptides_sorted["y_value"] = peptides_sorted.index +1
peptides_sorted["y_value"] = peptides_sorted["y_value"].astype("category")



plt.clf()
fig, axs = plt.subplots(len(conditions),
                        layout='constrained',
                        figsize=(14,10), 
                        sharex=True)

for idx, condition in enumerate(conditions):
    band = FilterPeptides(peptides_sorted, condition, cutoff_value)
    sns.pointplot(data=band,
                  x="value",
                  y="y_value",
                  hue="Log2 Int",
                  ax=axs[idx],
                  palette=sns.color_palette("rocket_r"))
    axs[idx].get_legend().remove()
    axs[idx].set_xlim(0,sequence_length)
    axs[idx].set(yticklabels=[])
    axs[idx].set(ylabel=None)
    axs[idx].set(yticks=[])
    if idx==(len(conditions)-1):
        axs[idx].set_xlabel("residue")
    else:
        axs[idx].set(xlabel=None)
    axs[idx].set_title(condition)
    
plt.savefig("coverage.pdf", dpi=300)

