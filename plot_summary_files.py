import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import subprocess
import re


os.chdir("T:\\Processing_MQ\\") #path to top of directory of MQ runs

run_mspicture = True
msPicture_command_path = str("C:\\Users\\path\\to\\mspicture.exe")
specific_date = "2405"
specific_year = "E24"

#start and end of peptide elution in minutes more or less to get msms id rate on gradient
pep_start = 15
pep_end = 80

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

def CleanUpList(files_list,specific_year):
    files_list = pd.Series(files_list, name="Raw file")
    files_df = pd.DataFrame(files_list)
    files_df = files_df[files_df['Raw file'].str.contains(specific_year)]
    files_df = files_df[files_df['Raw file'].str.contains('200ng')]
    files_df["raw_date"] = files_df["Raw file"].str.split("_").str[0]
    files_df["date"] = pd.to_datetime(files_df["raw_date"], format='%y%m%d').dt.date
    files_df["FAIMS"] = files_df["Raw file"].apply(checkFAIMS)
    files_df["Detector"] = files_df["Raw file"].apply(checkDetector)
    files_df["run_name"] = files_df["Raw file"].str.replace(".*combined_","", regex=True).str.replace(".raw.*", ".raw", regex=True)
    files_df = files_df.reset_index()
    try:
        files_df = files_df.drop("index")
    except KeyError:
        pass
    return files_df


def CleanUpList_NoAmount(files_list):
    files_list = pd.Series(files_list, name="Raw file")
    files_df = pd.DataFrame(files_list)
#    files_df = files_df[files_df['Raw file'].str.contains('200ng')]
    files_df["raw_date"] = files_df["Raw file"].str.split("_").str[0]
    files_df["date"] = pd.to_datetime(files_df["raw_date"], format='%y%m%d').dt.date
    files_df["FAIMS"] = files_df["Raw file"].apply(checkFAIMS)
    files_df["Detector"] = files_df["Raw file"].apply(checkDetector)
    files_df["run_name"] = files_df["Raw file"].str.replace(".*combined_","", regex=True).str.replace(".raw.*", ".raw", regex=True)
    files_df = files_df.reset_index()
    try:
        files_df = files_df.drop("index")
    except KeyError:
        pass
    return files_df

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

#HISTORICAL SERIES---------------------------------------------------------------------------

#grab number of proteins from protein groups files-------------

proteinGroups_list = glob.glob("*HeLa*/*/txt/*proteinGroups.txt")

proteinGroups_df = CleanUpList(proteinGroups_list, specific_year)

proteins = []
proteins_2 = []
for element, rawfile in proteinGroups_df.iterrows():
    proteinGroups = pd.read_csv(rawfile["Raw file"], delimiter="\t")
    proteins.append(NumberOfProteins(proteinGroups, min_peptides=0))
    proteins_2.append(NumberOfProteins(proteinGroups, min_peptides=2))
proteinGroups_df["proteins"] = proteins
proteinGroups_df["proteins_min2"] = proteins_2



#grab summary files---------------
#needs MQ runs with a single raw file
summaryFiles_list =  glob.glob("*HeLa*/*/txt/*summary.txt")
summaryFiles_df = CleanUpList(summaryFiles_list, specific_year)

for element, rawfile in summaryFiles_df.iterrows():
    if element==0:
        summary = pd.read_csv(rawfile["Raw file"], delimiter="\t")
        if len(summary)==2:
            summary = summary[summary["Raw file"]!="Total"]
            summary = summary.drop(["Raw file"], axis=1)
            summary["Detector"] = rawfile["Detector"]
            summary["FAIMS"] = rawfile["FAIMS"]
            summary["date"] = rawfile["date"]
            summary["Raw file"] = rawfile["Raw file"]
        else:
            print("%s has more than 1 run" % rawfile)
    else:
        summary_line = pd.read_csv(rawfile["Raw file"], delimiter="\t")
        if len(summary_line)==2:
            summary_line = summary_line[summary_line["Raw file"]!="Total"]
            summary_line = summary_line.drop(["Raw file"], axis=1)
            summary_line["Detector"] = rawfile["Detector"]
            summary_line["FAIMS"] = rawfile["FAIMS"]
            summary_line["date"] = rawfile["date"]
            summary_line["Raw file"] = rawfile["Raw file"]
            summary = pd.concat([summary, summary_line], ignore_index=True)
        else:
            print("%s has more than 1 run" % rawfile)
            
#grab msms files-----------------
msmsFiles_list =  glob.glob("*HeLa*/*/txt/*msms.txt")
msmsFiles_df = CleanUpList(msmsFiles_list, specific_year)

for element, rawfile in msmsFiles_df.iterrows():
    if element==0:
        msms = pd.read_csv(rawfile["Raw file"], 
                            delimiter="\t", 
                            low_memory=False)
        msms = msms.drop(["Raw file"], axis=1)
        msms["date"] = rawfile["date"]
        msms["Raw file"] = rawfile["Raw file"]
        msms["FAIMS"] = rawfile["FAIMS"]
        msms["Detector"] = rawfile["Detector"]
    else:
        try:
            msms_line = pd.read_csv(rawfile["Raw file"], 
                                    delimiter="\t", 
                                    low_memory=False)
            msms_line = msms_line.drop(["Raw file"], axis=1)
            msms_line["date"] = rawfile["date"]
            msms_line["Raw file"] = rawfile["Raw file"]
            msms_line["FAIMS"] = rawfile["FAIMS"]
            msms_line["Detector"] = rawfile["Detector"]
            msms = pd.concat([msms, msms_line], ignore_index=True)
        except:
            print("Fail at\n")
            print(rawfile["Raw file"])
            continue
            
            
            
#grab evidence files-----------------
evidenceFiles_list =  glob.glob("*HeLa*/*/txt/*evidence.txt")
evidenceFiles_df = CleanUpList(evidenceFiles_list, specific_year)

for element, rawfile in evidenceFiles_df.iterrows():
    if element==0:
        evidence = pd.read_csv(rawfile["Raw file"], delimiter="\t", low_memory=False)
        evidence = evidence.drop(["Raw file"], axis=1)
        evidence["date"] = rawfile["date"]
        evidence["Raw file"] = rawfile["Raw file"]
        evidence["FAIMS"] = rawfile["FAIMS"]
        evidence["Detector"] = rawfile["Detector"]
    else:
        evidence_line = pd.read_csv(rawfile["Raw file"], delimiter="\t", low_memory=False)
        evidence_line = evidence_line.drop(["Raw file"], axis=1)
        evidence_line["date"] = rawfile["date"]
        evidence_line["Raw file"] = rawfile["Raw file"]
        evidence_line["FAIMS"] = rawfile["FAIMS"]
        evidence_line["Detector"] = rawfile["Detector"]
        evidence = pd.concat([evidence, evidence_line], ignore_index=True)


#grab msScans files
msScansFiles_list =  glob.glob("*HeLa*/*/txt/msScans.txt")
msScansFiles_df = CleanUpList(msScansFiles_list, specific_year)

for element, rawfile in msScansFiles_df.iterrows():
    if element==0:
        msScans = pd.read_csv(rawfile["Raw file"], delimiter="\t", low_memory=False)
        msScans = msScans.drop(["Raw file"], axis=1)
        msScans["date"] = rawfile["date"]
        msScans["Raw file"] = rawfile["Raw file"]
        msScans["FAIMS"] = rawfile["FAIMS"]
        msScans["Detector"] = rawfile["Detector"]
        msScans = msScans[msScans["Retention time"]>pep_start]
        msScans = msScans[msScans["Retention time"]<pep_end]
    else:
        msScans_line = pd.read_csv(rawfile["Raw file"], delimiter="\t", low_memory=False)
        msScans_line = msScans_line.drop(["Raw file"], axis=1)
        msScans_line["date"] = rawfile["date"]
        msScans_line["Raw file"] = rawfile["Raw file"]
        msScans_line["FAIMS"] = rawfile["FAIMS"]
        msScans_line["Detector"] = rawfile["Detector"]
        msScans_line = msScans_line[msScans_line["Retention time"]>pep_start]
        msScans_line = msScans_line[msScans_line["Retention time"]<pep_end]
        msScans = pd.concat([msScans, msScans_line], ignore_index=True)

    
#add identifiable attributes to data frames--------    
summary = AddType(summary)
proteinGroups_df = AddType(proteinGroups_df)
msms = AddType(msms)
evidence = AddType(evidence)
msScans = AddType(msScans)

summary = CleanUpRawFile(summary)
proteinGroups_df = CleanUpRawFile(proteinGroups_df)
msms = CleanUpRawFile(msms)
msScans = CleanUpRawFile(msScans)

evidence = CleanUpRawFile(evidence)



#merge summary and proteinGroups statistics----------------------
QC_df = pd.merge(right=proteinGroups_df,
          left=summary,
          on=["Detector", "FAIMS", "date", "Raw file", "type"], 
          how="outer")


#plots-------------------------
plt.clf()
fig, axs = plt.subplots(9, layout='constrained',
                          figsize=(3.5 * 4, 3.5 * 6))

#Standard line plots
QCplot(QC_df, "proteins", axs[0], "Identified proteins")
QCplot(QC_df, "proteins_min2", axs[1], "Identified proteins min. 2 peptides")
QCplot(QC_df, "MS", axs[2], "MS1 scans")
QCplot(QC_df, "MS/MS", axs[3], "MS2 scans")
QCplot(QC_df, "MS/MS submitted", axs[4], "MS2 submitted")
QCplot(QC_df, "MS/MS identified [%]", axs[5], "MS2 identification rate")
#violin plot
sns.violinplot(data=evidence, x="date",
              y="Uncalibrated mass error [ppm]",
              hue="type",
              style="type",
              errorbar="sd",
              estimator="median",
              inner=None,
              markers=True,
              dashes=False,
              ax=axs[6], rasterized=True)
axs[6].tick_params(axis='x', rotation=90)
axs[6].set_title("precursor mass error")
axs[6].legend(bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)

sns.lineplot(data=msms,x="date",
              y="Intensity coverage",
              hue="type",
              estimator="median",
              markers="o",
              errorbar="sd",
              err_style="bars",
              ax=axs[7], rasterized=True)
axs[7].set_title("MS2 median intensity coverage")
axs[7].legend(bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)

sns.lineplot(data=msScans,x="date",
             y="Total ion current",
             hue="type",
             estimator="mean",
             markers="o",
             errorbar="sd",
             err_style="bars",
             ax=axs[8],
             rasterized=True)
axs[8].set_title("MS1 mean intensity on gradient")
axs[8].legend(bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)
             

plt.tight_layout()
plt.savefig("QC_summary\\performance.pdf", dpi=300)


#PER-RUN REPORT------------------------------------------

run_lists = glob.glob("*HeLa*/*/txt/")
run_df = CleanUpList_NoAmount(run_lists)

#specify date
try:
    run_df = run_df[run_df["Raw file"].str.contains(specific_date)]

except:
    pass


for element, search_name in run_df.iterrows():
    if run_mspicture==True:
        run_mspicture_iteration = True
    else:
        run_mspicture_iteration = False
    
    search_name=AddType(search_name)
    
    figname = str("QC_summary\\" + str(search_name["date"]) +"_" + search_name["type"] + "_" +search_name["run_name"] +"_report.pdf")
    report_exists = len(glob.glob(figname))
    
    if report_exists>0:
        continue
    
    if run_mspicture==True:
        #execute ProteoWizard MSPicture for nice mz vs RT picture of the run
        try:
            msimage_path_list = MSPicture(search_name["Raw file"], msPicture_command_path)
            msimage_file = msimage_path_list[0]
            run_mspicture_iteration = True
        except subprocess.CalledProcessError:
            print("cannot run MSPicture for file "+ search_name["run_name"] + "because of subprocessError, skipping...")
            run_mspicture_iteration = False
        except NameError:
            print("cannot run MSPicture for file "+ search_name["run_name"] + " because NameError, skipping...")
            run_mspicture_iteration = False            
    try:
        msScan = pd.read_csv(str(search_name["Raw file"]+"msScans.txt"), delimiter="\t", low_memory=False)
        msScan = AddType(msScan, from_other_source=True, other_source=search_name)
        msScan["Cycle time rolling ave"] = msScan[["Cycle time"]].rolling(100).mean()
        ms_scan_gradient = msScan[(msScan["Retention time"] > pep_start) & (msScan["Retention time"] <pep_end)]
        msms_id_rate_gradient = ms_scan_gradient["MS/MS identification rate [%]"].mean().round(3)
        msScan_produced=True
    except FileNotFoundError:
        msScan_produced=False
        continue
    except:
        print(str("error at file " + search_name["run_name"] + ", skipping..."))
        continue
    #grab tables-----------
    summary = pd.read_csv(str(search_name["Raw file"]+"summary.txt"), delimiter="\t", low_memory=False)
    summary = AddType(summary, from_other_source=True, other_source=search_name)
    msms = pd.read_csv(str(search_name["Raw file"]+"msms.txt"), delimiter="\t", low_memory=False)
    msms = AddType(msms, from_other_source=True, other_source=search_name)
    try:
        msms_Scans =  pd.read_csv(str(search_name["Raw file"]+"msmsScans.txt"), delimiter="\t", low_memory=False)
        msms_Scans = msms_Scans.fillna("-")
        msms_Scans = AddType(msms_Scans, from_other_source=True, other_source=search_name)
    except FileNotFoundError:
        print(str("no msmsScans for file "+ search_name["run_name"] + ", skipping..."))
        continue
    evidence = pd.read_csv(str(search_name["Raw file"]+"evidence.txt"), delimiter="\t", low_memory=False)
    evidence = AddType(evidence, from_other_source=True, other_source=search_name)
    proteinGroups = pd.read_csv(str(search_name["Raw file"]+"proteinGroups.txt"), delimiter="\t", low_memory=False)
    proteinGroups = AddType(proteinGroups, from_other_source=True, other_source=search_name)
    
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
    
    #different summary parameters if msscans file is available-----
    if msScan_produced==False:
        summary_table = {"parameter" : ["ms1_spectra",
                                      "ms2_spectra",
                                      "proteins_with_contaminants",
                                      "contaminants",
                                      "proteins, no rev. cont. no only site",
                                      "proteins_min_2_peptides",
                                      "quantified_proteins",
                                      "msms_id_rate [%]",
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
                                      isotope_patterns_detected,
                                      isotope_atterns_sequenced_z_1]}
        
        summary_table = pd.DataFrame.from_dict(summary_table)
    elif msScan_produced==True:

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
        if run_mspicture==True and run_mspicture_iteration==True:
            msimage = plt.imread(msimage_file)
            ap1 = plt.subplot2grid((11,3), (0,0), colspan=3)
            ap2 = plt.subplot2grid((11,3), (1,0), colspan=3)
            ap3 = plt.subplot2grid((11,3), (2,0), colspan=3)
            ap4 = plt.subplot2grid((11,3), (3,0), colspan=1)
            ap5 = plt.subplot2grid((11,3), (3,1), colspan=1)
            ap6 = plt.subplot2grid((11,3), (3,2), colspan=1)
            ap7 = plt.subplot2grid((11,3), (4,0), colspan=3)
            ap8 = plt.subplot2grid((11,3), (5,0), colspan=3)
            ap9 = plt.subplot2grid((11,3), (6,0), colspan=3)
            ap10 = plt.subplot2grid((11,3), (7,0), colspan=3)
            ap11 = plt.subplot2grid((11,3), (8,0), colspan=3, rowspan=3)
            # fig, axs = plt.subplots(nrows=3, ncols=3, layout='constrained',
            #                          figsize=(3.5 * 4, 3.5 * 6),
            #                          gridspec_kw={'height_ratios' : [1,1,1,1,1,1,1,1,3]})
        else:
            ap1 = plt.subplot2grid((10,3), (0,0), colspan=3)
            ap2 = plt.subplot2grid((10,3), (1,0), colspan=3)
            ap3 = plt.subplot2grid((10,3), (2,0), colspan=3)
            ap4 = plt.subplot2grid((10,3), (3,0), colspan=1)
            ap5 = plt.subplot2grid((10,3), (3,1), colspan=1)
            ap6 = plt.subplot2grid((10,3), (3,2), colspan=1)
            ap7 = plt.subplot2grid((10,3), (4,0), colspan=3)
            ap8 = plt.subplot2grid((10,3), (5,0), colspan=3)
            ap9 = plt.subplot2grid((10,3), (6,0), colspan=3)
            ap10 = plt.subplot2grid((10,3), (7,0), colspan=3)
        
        ap1.table(cellText=summary_table.values, colLabels=summary_table.columns, loc='center')
        ap1.axis("off")
                
        sns.lineplot(data=msScan, x="Retention time", y="Total ion current", ax=ap2, rasterized=True)
        secondary = ap2.twinx()
        sns.lineplot(data=msScan, x="Retention time", y="MS/MS identification rate [%]", ax=secondary, color="orange")
        ap2.set_title("TIC & MSMS Id Rate")
        
        sns.lineplot(data=msScan, x="Retention time", y="MS/MS count", ax=ap3, rasterized=True)
        ax2=ap3.twinx()
        sns.lineplot(data=msScan, x="Retention time", y="MS/MS / s", color="orange", ax=ax2)
        ap3.set_title("MSMS count & scan per second")
        
        
        sns.histplot(data=evidence, x="Uncalibrated mass error [ppm]", ax=ap4, rasterized=True)
        ap4.axvline(x=evidence["Uncalibrated mass error [ppm]"].mean(), color="orange", lw=2.5)
        #axs[3].text(3, 500, str((evidence["Uncalibrated mass error [ppm]"].mean().round(3))+" ppm"))
        ap4.set_title("precursor mass error")
        ap4.set_xlim(-5,5)
        
        sns.histplot(data=msms, x="Intensity coverage", ax=ap5, rasterized=True)
        ap5.set_title("Intensity coverage in msms")
        ap5.set_xlim(0,1)
        
        sns.histplot(data=msms, x="Isotope index", ax=ap6, rasterized=True)
        ap6.set_title("Isotope index")
        ap6.set_xlim(-5,5)
        
        sns.scatterplot(data=msms_Scans[msms_Scans["Identified"]=="-"],
                        x='Retention time',
                        y='m/z',
                        edgecolor=None,
                        color="gray",
                        ax=ap7, rasterized=True)
        
        ax6=ap7.twinx()
        sns.scatterplot(data=msms_Scans[msms_Scans["Identified"]=="+"],
                x='Retention time',
                y="m/z",
                edgecolor=None,
                color="orange",
                ax=ax6, rasterized=True)
        
        
        sns.regplot(data=evidence,
                        x='Retention time',
                        y='Uncalibrated mass error [ppm]',
                        ax=ap8, scatter_kws={"rasterized":True, "edgecolor":None})
        ap8.axhline(y=evidence["Uncalibrated mass error [ppm]"].mean(), color="orange", lw=2.5)

        ap8.set_title("mass error vs retention time")
        sns.regplot(data=evidence,
                        x='m/z',
                        y='Uncalibrated mass error [ppm]',
                        ax=ap9, scatter_kws={"rasterized":True, "edgecolor":None})
        ap9.axhline(y=evidence["Uncalibrated mass error [ppm]"].mean(), color="orange", lw=2.5)

        ap9.set_title("mass error vs m/z")

        
        
        sns.lineplot(data=msScan, x="Retention time", y="Cycle time", ax=ap10, rasterized=True)
        ax10 = ap10.twinx()
        sns.lineplot(data=msScan, x="Retention time", y="Cycle time rolling ave", color="orange", ax=ax10, rasterized=True)
        ap10.set_title("cycle times")


        if run_mspicture==True and run_mspicture_iteration==True:
            ap11.imshow(msimage)
        
        
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.suptitle(str("Report for" + search_name["run_name"]+ " , "+ str(search_name["type"])+ " " + str(search_name["date"])))

        plt.savefig(figname, dpi=200)
        
