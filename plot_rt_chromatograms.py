# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 11:17:37 2023

@author: andrea.graziadei
"""
from pyteomics import mzml, mass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import glob


run_list = glob.glob("*mzML")

tolerance = 3 #tolerance in ppm
charge_state = 2
rt_min = 50 

mass_list = [985.522,
            1224.6189,
            990.5589,
            900.5524,
            843.4582,
            1389.6503,
            1171.5861,
            1545.7766,
            1114.6374,
            1600.8084,
            1488.7704,
            995.589,
            1144.5905,
            1358.7326,
            1572.8279]

hydrophobicity_factor = [7.56,
                        15.5,
                        15.52,
                        17.65,
                        19.15,
                        25.88,
                        25.24,
                        28.37,
                        32.18,
                        34.5,
                        34.96,
                        37.3,
                        40.42,
                        41.18,
                        46.66]


mass_df = pd.DataFrame({"mass":mass_list})
mass_df["HF"] = hydrophobicity_factor

def calc_mz(mass_value, charge_state):
    mz = (mass_value + (charge_state*mass.calculate_mass(formula="H")))/charge_state
    return mz

def mz_min(ref_mz, tolerance):
    ref_tolerance = (tolerance/1000000)*ref_mz
    ref_mz_lower = ref_mz - ref_tolerance
    return ref_mz_lower

def mz_max(ref_mz, tolerance):
    ref_tolerance = (tolerance/1000000)*ref_mz
    ref_mz_upper = ref_mz + ref_tolerance
    return ref_mz_upper

def AddPeak(spectrum, mz_value, mz_min_value, mz_max_value, mass_value, HF_value, run_value):
    rt = spectrum["scanList"]["scan"][0]["scan start time"]
    intensity_indices = np.where(np.logical_and(spectrum["m/z array"]>mz_min_value,
                                                spectrum["m/z array"]<mz_max_value))
    intensity_values = spectrum["intensity array"][intensity_indices]
    intensity_values = intensity_values[np.where(intensity_values >500000)]
    intensity_numbers = len(intensity_values)
    temp_d = {"mass":[mass_value]*intensity_numbers,
              "rt":[rt]*intensity_numbers,
              "mz":[mz_value]*intensity_numbers,
              "intensity":intensity_values.tolist(),
              "HF":[HF_value]*intensity_numbers,
              "run":[run_value]*intensity_numbers}
    temp_d = pd.DataFrame(temp_d)
    return temp_d

mass_df["charge state"] = charge_state
mass_df["mz"] = mass_df.apply(lambda x: calc_mz(x["mass"], charge_state), axis=1)
mass_df["mzmin"] = mass_df.apply(lambda x: mz_min(x["mz"], tolerance), axis=1)
mass_df["mzmax"] = mass_df.apply(lambda x: mz_max(x["mz"], tolerance), axis=1)


d = {'mass':[], 'rt': [], 'mz': [], 'intensity':[], "HF":[], "run":[]}
rt_intensity_mz = pd.DataFrame(d)

for run_name in run_list:
    spectra_list = []
    f = mzml.MzML(run_name)
    for spectrum in f:
        if spectrum["ms level"]==1:
            spectra_list.append(spectrum)
    for index, rt_peptide in mass_df.iterrows():
        rt_peptide["run"]=run_name
        print(rt_peptide["mass"], rt_peptide["run"])
        for ms1_spectrum in spectra_list:
            if rt_peptide["mass"]!=1572.8279:
                    temp_d = AddPeak(ms1_spectrum,
                                     rt_peptide["mz"],
                                     rt_peptide["mzmin"],
                                     rt_peptide["mzmax"],
                                     rt_peptide["mass"],
                                     rt_peptide["HF"],
                                     rt_peptide["run"])
                    if len(temp_d)<1:
                        continue
                    else:
                        rt_intensity_mz = pd.concat([rt_intensity_mz, temp_d], ignore_index=True)
            else:
                if ms1_spectrum["scanList"]["scan"][0]["scan start time"] > rt_min:
                        temp_d = AddPeak(ms1_spectrum,
                                         rt_peptide["mz"],
                                         rt_peptide["mzmin"],
                                         rt_peptide["mzmax"],
                                         rt_peptide["mass"],
                                         rt_peptide["HF"],
                                         rt_peptide["run"])
                        if len(temp_d)<1:
                            continue
                        else:
                            rt_intensity_mz = pd.concat([rt_intensity_mz, temp_d], ignore_index=True)
    f.reset()




rt_intensity_mz["mass_round"] = rt_intensity_mz["mass"].round(2)
rt_intensity_mz["mass_round"] = rt_intensity_mz["mass_round"].astype("category")
rt_intensity_mz["HF"] = rt_intensity_mz["HF"].astype("category")

d = {'mass':[], 'rt': [], 'mz': [], 'intensity':[], "HF":[], "mass_round":[]}
peak_rt = pd.DataFrame(d)

peak_rt = rt_intensity_mz.groupby(['mass_round', 'run'])['intensity'].max()


idx = rt_intensity_mz.groupby(['mass_round', 'run'])['intensity'].transform(max) == rt_intensity_mz['intensity']

peak_rt = rt_intensity_mz[idx]
peak_rt["mean_intensity"] = peak_rt.groupby(['mass_round'])['intensity'].transform('mean')
peak_rt["intensity_deviation"] = (1-(peak_rt["intensity"]/peak_rt["mean_intensity"]))*100
peak_rt["mean_rt"] = peak_rt.groupby(['mass_round'])['rt'].transform('mean')
peak_rt["rt_deviation"] = np.absolute(1-(peak_rt["rt"]/peak_rt["mean_rt"]))*100


plt.clf()

fig, axs = plt.subplots(3,
                        layout='constrained',
                        figsize=(8,10), 
                        sharex=False)

sns.lineplot(data=peak_rt, x="rt", y="HF", ax=axs[0], marker="o", hue="run")
axs[0].set_ylabel("hydrophobicity factor")
axs[0].set_xlabel("rt (minutes)")
axs[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

sns.lineplot(data=peak_rt, x="run", y="rt", ax=axs[1], marker="o", hue="mass_round")
axs[1].set_xlabel("run")
axs[1].set_ylabel("rt (minutes)")
axs[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
axs[1].set_xticklabels(list(peak_rt["run"].unique()), rotation = 90)

sns.lineplot(data=peak_rt, x="run", y="intensity", ax=axs[2],  marker="o", hue="mass_round")
axs[2].set_xlabel("run")
axs[2].set_ylabel("intensity")
axs[2].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
axs[2].set_xticklabels(list(peak_rt["run"].unique()), rotation = 90)

plt.savefig("peak_position.png", dpi=300)

plt.clf()

fig, axs = plt.subplots(len(run_list),
                        layout='constrained',
                        figsize=(8,14), 
                        sharex=True)

for idx, element in enumerate(run_list):
    axs_position = idx
    sns.lineplot(data=rt_intensity_mz[rt_intensity_mz["run"]==element],
                 x="rt",
                 y="intensity",
                 hue="mass_round",
                 ax=axs[axs_position])
    if axs_position==0:
        axs[axs_position].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    else:
        axs[axs_position].get_legend().remove()
    axs[axs_position].set_title(element)
    axs[axs_position].set(yticklabels=[])
    axs[axs_position].set(ylabel=None)
#plt.tight_layout()
plt.savefig("chromatogram.png", dpi=300)


plt.clf()

fig, axs = plt.subplots(len(list(peak_rt["mass_round"].unique())),
                        layout='constrained',
                        figsize=(8,14), 
                        sharex=True)
for idx, rt_peptide in enumerate(list(peak_rt["mass_round"].unique())):    
    sns.lineplot(data=peak_rt, x="run", y="rt_deviation", ax=axs[idx],  marker="o")
    axs[idx].set_xlabel("run")
    axs[idx].set_ylim(0,20)
    axs[idx].set_ylabel("rt %ge variance")
    axs[idx].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    axs[idx].set_xticklabels([])
    axs[idx].set_title(str(rt_peptide))
    if rt_peptide==list(peak_rt["mass_round"].unique())[-1]:
        axs[idx].set_xticklabels(list(peak_rt["run"].unique()), rotation = 90)
    

plt.suptitle("rt variation % from mean")
plt.savefig("rt_variation.png", dpi=300)


plt.clf()

fig, axs = plt.subplots(len(list(peak_rt["mass_round"].unique())),
                        layout='constrained',
                        figsize=(8,14), 
                        sharex=True,
                        sharey=True)
for idx, rt_peptide in enumerate(list(peak_rt["mass_round"].unique())):    
    sns.lineplot(data=peak_rt,
                 x="run",
                 y="intensity_deviation",
                 ax=axs[idx],
                 marker="o",
                 estimator="max")
    axs[idx].set_xlabel("run")
    axs[idx].set_ylabel(None)
    axs[idx].set_ylim(-20,20)
    axs[idx].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    axs[idx].set_xticklabels([])
    axs[idx].set_title(str(rt_peptide))
    if rt_peptide==list(peak_rt["mass_round"].unique())[-1]:
        axs[idx].set_xticklabels(list(peak_rt["run"].unique()), rotation = 90)


plt.suptitle("intensity variation % from mean")
plt.savefig("intensity_variation.png", dpi=300)
