# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import subprocess
import re
import argparse
import sys


parser = argparse.ArgumentParser(prog='Maxquant Report generation',
                    description='generates report from maxquant output files',
                    epilog='usage: report.py --infile E230203_05_QC_AG_100_HeLa200ng_OT_inj2.raw')

parser.add_argument("--infile", metavar="E230203_05_QC_AG_100_HeLa200ng_OT_inj2.raw",
                    required=True,
                    type=str,
                    help= str("the raw file at the top of the maxquant directory tree \n" + 
                              "the presence of the word OT or IT will determine ion trap or orbitrap\n" +
                              "the presence of the word FAIMS will classify it as a run with FAIMS\n" +
                              "the file name has to start with a letter for the instrument and a date in format YYMMDD followed by a _\n" +
                              "example file name E230203_05_QC_AG_100_HeLa200ng_OT_FAIMS_inj2.raw"))

parser.add_argument("--run_mspicture", metavar="",
                    required=False,
                    type=bool,
                    default=True,
                    help="whether to generate the RT vs m/z plot with Proteowizard mspicture. Set to True or False. Specify path to mspicture in the script")

parser.add_argument("--gradient_start", metavar="",
                    required=False,
                    type=float,
                    default=15,
                    help="start of gradient in minutes (to trim calculation of MS/MS id rate to only peptide-containing part of run")


parser.add_argument("--gradient_end", metavar="",
                    required=False,
                    type=float,
                    default=80,
                    help="start of gradient in minutes (to trim calculation of MS/MS id rate to only peptide-containing part of run")


args = parser.parse_args()


rawfile_name = args.infile

run_mspicture = args.run_mspicture


#start and end of peptide elution in minutes more or less to get msms id rate on gradient
pep_start = args.gradient_start
pep_end =args.gradient_end


#some config options

msPicture_command_path = str("C:\\Users\\andrea.graziadei\\AppData\\Local\\Apps\\ProteoWizard 3.0.22314.0cd8422 64-bit\\mspicture.exe")


#functions-----------------------------------------------------------
def checkFAIMS(raw_file):
    if "faims" in str(raw_file).lower():
        return "FAIMS"
    else:
        return "noFAIMS"

def checkDetector(raw_file):
    if "IT" in raw_file:
        return "IT"
    elif "OT" in raw_file:
        return "OT"
    else:
        return "other"


def AddType(df, from_other_source=False, other_source=None):
    if from_other_source==False:
        df["type"] = df["FAIMS"]+"_"+df["Detector"]
    if from_other_source==True:
        df["type"] = other_source["type"]
    return df

def CleanUpRawFile(df):
    df["Raw file"] = df["Raw file"].str.replace(".*\\\combined_","", regex=True).str.replace("raw.*", "raw")
    return df

def NumberOfProteins(df, min_peptides=0):
    '''extract number of proteins from ProteinGroups file'''
    #filter out reverse, contaminant, identified by site.
    df = df[df["Reverse"]!="+"]
    df = df[df["Potential contaminant"]!="+"]
    df = df[df["Only identified by site"]!="+"]
    df = df[df["Peptides"]>=min_peptides]
    number_of_proteins = len(df)
    return number_of_proteins

def QCplot(df, y_value, axs_value, title=None):
    '''plot time series of a particular value by detector and faims'''
    sns.lineplot(data=df, x="date",
                 y=y_value,
                 hue="type",
                 markers=True,
                 estimator="mean",
                 errorbar="sd",
                 style="type",
                 dashes=False,
                 ax=axs_value)
#    axs_value.tick_params(axis='x', rotation=90)
    axs_value.set_title(title)
    axs_value.legend(bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)


def MSPicture(raw_file, msPicture_command_path):
    raw_file_name = re.sub("combined_", "", raw_file)
    raw_file_name = re.sub("raw.*", "raw", raw_file_name)
    output_dir = str(raw_file_name+".png")
    mspicture_command = [msPicture_command_path,
                         "-o",
                         output_dir,
                         "--width",
                         "1000",
                         "--height",
                         "1000",
                         raw_file_name]
    subprocess.check_call(mspicture_command)
    run_image_path = glob.glob(str(output_dir+"\\*png"))
    return run_image_path


#generate report----------------------


raw_file_date = rawfile_name.split("_")[0][1:]
raw_file_date = datetime.strptime(raw_file_date, '%m%d%Y').date()

FAIMS = checkFAIMS(rawfile_name)

detector = checkDetector(rawfile_name)


figname = str(str(raw_file_date) +
              "_" + 
              FAIMS + 
              "_"+ 
              detector + 
              "_" + 
              rawfile_name +
              "_report.pdf")

if run_mspicture==True:
    #execute ProteoWizard MSPicture for nice mz vs RT picture of the run
    try:
        msimage_path_list = MSPicture(rawfile_name, msPicture_command_path)
        msimage_file = msimage_path_list[0]
    except subprocess.CalledProcessError:
        print("cannot run MSPicture for file "+ rawfile_name + ", skipping...")
        run_mspicture_iteration = False
try:
    msScan = pd.read_csv(str(".\\combined\\txt\\"+"msScans.txt"), delimiter="\t", low_memory=False)
    msScan["Cycle time rolling ave"] = msScan[["Cycle time"]].rolling(100).mean()
    ms_scan_gradient = msScan[(msScan["Retention time"] > pep_start) & (msScan["Retention time"] <pep_end)]
    msms_id_rate_gradient = ms_scan_gradient["MS/MS identification rate [%]"].mean().round(3)
    msScan_produced=True
except FileNotFoundError:
    msScan_produced=False
    print("no MS scan file was generated. Exiting...")
    sys.exit()
except:
    print(str("error at msScan file, Exiting..."))
    sys.exit()

#grab tables-----------
summary = pd.read_csv(str(".\\combined\\txt\\"+"summary.txt"), delimiter="\t", low_memory=False)
msms = pd.read_csv(str(".\\combined\\txt\\"+"msms.txt"), delimiter="\t", low_memory=False)

try:
    msms_Scans =  pd.read_csv(str(".\\combined\\txt\\" +"msmsScans.txt"), delimiter="\t", low_memory=False)
    msms_Scans = msms_Scans.fillna("-")
except FileNotFoundError:
    print(str("no msmsScans file, Exiting..."))
    sys.exit()
    
evidence = pd.read_csv(str(".\\combined\\txt\\"+"evidence.txt"), delimiter="\t", low_memory=False)
proteinGroups = pd.read_csv(str(".\\combined\\txt\\"+"proteinGroups.txt"), delimiter="\t", low_memory=False)

#summary parameters for table in report
ms1_spectra = summary["MS"][0]
ms2_spectra = summary["MS/MS"][0]

proteins_with_contaminants = len(proteinGroups[proteinGroups["Reverse"]!="+"])
contaminants = len(proteinGroups[proteinGroups["Potential contaminant"]=="+"])
proteins = NumberOfProteins(proteinGroups)
proteins_min_2_peptides = NumberOfProteins(proteinGroups, min_peptides=2)
try:
    quantified_proteins = len(proteinGroups[proteinGroups["Reverse"]!="+"][proteinGroups["iBAQ"]>0][proteinGroups["Potential contaminant"]!="+"][proteinGroups["Only identified by site"]!="+"])
except KeyError:
    quantified_proteins = "iBAQ not enabled"
msms_id_rate = summary["MS/MS identified [%]"][0]
isotope_patterns_detected = summary["Isotope patterns"][0]
isotope_atterns_sequenced_z_1 = summary["Isotope patterns sequenced (z>1)"][0]



summary_table = {"parameter" : ["ms1_spectra",
                              "ms2_spectra",
                              "proteins_with_contaminants",
                              "contaminants",
                              "proteins, no rev. cont. no only site",
                              "proteins_min_2_peptides",
                              "quantified_proteins",
                              "msms_id_rate [%]",
                              "msms_id_rate (on gradient) [%]",
                              "isotope_patterns_detected",
                              "isotope_atterns_sequenced_z_1"],
                 "number" : [ms1_spectra,
                             ms2_spectra,
                             proteins_with_contaminants,
                             contaminants,
                             proteins,
                             proteins_min_2_peptides,
                             quantified_proteins,
                             msms_id_rate,
                             msms_id_rate_gradient,
                             isotope_patterns_detected,
                             isotope_atterns_sequenced_z_1]}        
summary_table = pd.DataFrame.from_dict(summary_table)


if msScan_produced==True:
    plt.clf()
    fig = plt.figure(figsize=(3.5 * 4, 3.5 * 6))
    if run_mspicture==True:
        msimage = plt.imread(msimage_file)
        ap1 = plt.subplot2grid((9,3), (0,0), colspan=3)
        ap2 = plt.subplot2grid((9,3), (1,0), colspan=3)
        ap3 = plt.subplot2grid((9,3), (2,0), colspan=3)
        ap4 = plt.subplot2grid((9,3), (3,0), colspan=1)
        ap5 = plt.subplot2grid((9,3), (3,1), colspan=1)
        ap6 = plt.subplot2grid((9,3), (3,2), colspan=1)
        ap7 = plt.subplot2grid((9,3), (4,0), colspan=3)
        ap8 = plt.subplot2grid((9,3), (5,0), colspan=3)
        ap9 = plt.subplot2grid((9,3), (6,0), colspan=3, rowspan=3)
        # fig, axs = plt.subplots(nrows=3, ncols=3, layout='constrained',
        #                          figsize=(3.5 * 4, 3.5 * 6),
        #                          gridspec_kw={'height_ratios' : [1,1,1,1,1,1,1,1,3]})
    else:
        msimage = plt.imread(msimage_file)
        ap1 = plt.subplot2grid((6,3), (0,0), colspan=3)
        ap2 = plt.subplot2grid((6,3), (1,0), colspan=3)
        ap3 = plt.subplot2grid((6,3), (2,0), colspan=3)
        ap4 = plt.subplot2grid((6,3), (3,0), colspan=1)
        ap5 = plt.subplot2grid((6,3), (3,1), colspan=1)
        ap6 = plt.subplot2grid((6,3), (3,2), colspan=1)
        ap7 = plt.subplot2grid((6,3), (4,0), colspan=3)
        ap8 = plt.subplot2grid((6,3), (5,0), colspan=3)
    
    sns.lineplot(data=msScan, x="Retention time", y="Total ion current", ax=ap1)
    secondary = ap1.twinx()
    sns.lineplot(data=msScan, x="Retention time", y="MS/MS identification rate [%]", ax=secondary, color="orange")
    ap1.set_title("TIC & MSMS Id Rate")
    
    sns.lineplot(data=msScan, x="Retention time", y="MS/MS count", ax=ap2)
    ax2=ap2.twinx()
    sns.lineplot(data=msScan, x="Retention time", y="MS/MS / s", color="orange", ax=ax2)
    ap2.set_title("MSMS count & scan per second")
    
    
    sns.lineplot(data=msScan, x="Retention time", y="Cycle time", ax=ap3)
    ax3 = ap3.twinx()
    sns.lineplot(data=msScan, x="Retention time", y="Cycle time rolling ave", color="orange", ax=ax3)
    ap3.set_title("cycle times")
    
    sns.histplot(data=evidence, x="Uncalibrated mass error [ppm]", ax=ap4)
    ap4.axvline(x=evidence["Uncalibrated mass error [ppm]"].mean(), color="orange", lw=2.5)
    #axs[3].text(3, 500, str((evidence["Uncalibrated mass error [ppm]"].mean().round(3))+" ppm"))
    ap4.set_title("precursor mass error")
    ap4.set_xlim(-5,5)
    
    sns.histplot(data=msms, x="Intensity coverage", ax=ap5)
    ap5.set_title("Intensity coverage in msms")
    ap5.set_xlim(0,1)
    
    sns.histplot(data=msms, x="Isotope index", ax=ap6)
    ap6.set_title("Isotope index")
    
    sns.scatterplot(data=msms_Scans[msms_Scans["Identified"]=="-"],
                    x='Retention time',
                    y='m/z',
                    edgecolor=None,
                    color="gray",
                    ax=ap7)
    
    ax6=ap7.twinx()
    sns.scatterplot(data=msms_Scans[msms_Scans["Identified"]=="+"],
            x='Retention time',
            y="m/z",
            edgecolor=None,
            color="orange",
            ax=ax6)
    
    ap8.table(cellText=summary_table.values, colLabels=summary_table.columns, loc='center')
    ap8.axis("off")
    
    if run_mspicture==True:
        ap9.imshow(msimage)
    
    
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle(str("Report for" + 
                     rawfile_name + 
                     " , " +
                     detector + 
                     " " +
                     FAIMS +
                     " " +
                     str(raw_file_date)))

    plt.savefig(figname, dpi=200)
        