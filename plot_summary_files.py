import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os

os.chdir("S:\\Processing_MQ\\")

#functions-----------------------------------------------------------
def checkFAIMS(raw_file):
    if "faims" in str(raw_file).lower():
        return "FAIMS"
    else:
        return "no FAIMS"

def checkDetector(raw_file):
    if "IT" in raw_file:
        return "IT"
    elif "OT" in raw_file:
        return "OT"
    else:
        return "other"

def CleanUpList(files_list):
    files_list = pd.Series(files_list, name="Raw file")
    files_df = pd.DataFrame(files_list)


    files_df = files_df[files_df['Raw file'].str.contains('200ng')]
    files_df["raw_date"] = files_df["Raw file"].str.split("_").str[0]
    files_df["date"] = pd.to_datetime(files_df["raw_date"], format='%y%m%d').dt.date
    files_df["FAIMS"] = files_df["Raw file"].apply(checkFAIMS)
    files_df["Detector"] = files_df["Raw file"].apply(checkDetector)
    files_df = files_df.reset_index()
    try:
        files_df = files_df.drop("index")
    except KeyError:
        pass
    return files_df

def AddType(df):
    df["type"] = df["FAIMS"]+"_"+df["Detector"]
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
    axs_value.tick_params(axis='x', rotation=90)
    axs_value.set_title(title)
    axs_value.legend(bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)




#grab number of proteins from protein groups files-------------

proteinGroups_list = glob.glob("*HeLa*/*/txt/*proteinGroups.txt")

proteinGroups_df = CleanUpList(proteinGroups_list)

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
summaryFiles_df = CleanUpList(summaryFiles_list)

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
msmsFiles_df = CleanUpList(msmsFiles_list)

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
        msms_line = pd.read_csv(rawfile["Raw file"], 
                                delimiter="\t", 
                                low_memory=False)
        msms_line = msms_line.drop(["Raw file"], axis=1)
        msms_line["date"] = rawfile["date"]
        msms_line["Raw file"] = rawfile["Raw file"]
        msms_line["FAIMS"] = rawfile["FAIMS"]
        msms_line["Detector"] = rawfile["Detector"]
        msms = pd.concat([msms, msms_line], ignore_index=True)


#grab evidence files-----------------
evidenceFiles_list =  glob.glob("*HeLa*/*/txt/*evidence.txt")
evidenceFiles_df = CleanUpList(evidenceFiles_list)

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

   
#grab MS scans---------------------

#add identifiable attributes to data frames--------    
summary = AddType(summary)
proteinGroups_df = AddType(proteinGroups_df)
msms = AddType(msms)
evidence = AddType(evidence)

summary = CleanUpRawFile(summary)
proteinGroups_df = CleanUpRawFile(proteinGroups_df)
msms = CleanUpRawFile(msms)
evidence = CleanUpRawFile(evidence)



#merge summary and proteinGroups statistics----------------------
QC_df = pd.merge(right=proteinGroups_df,
         left=summary,
         on=["Detector", "FAIMS", "date", "Raw file", "type"], 
         how="outer")


#plots-------------------------
plt.clf()
fig, axs = plt.subplots(7, layout='constrained',
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
             ax=axs[6])
axs[6].tick_params(axis='x', rotation=90)
axs[6].set_title("precursor mass error")
axs[6].legend(bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)

plt.tight_layout()
plt.savefig("QC_summary\\performance.pdf", dpi=300)

# figure_list = [identifiedProteins, msmsSubmitted, msmsIDRate]


# fig, axs = plt.subplots(len(figure_list))
# for i, element in enumerate(figre_list):
    

# plt.close()
# plt.clf()
# sns.violinplot(data=evidence, x="date",
#              y="Uncalibrated mass error [ppm]",
#              hue="type",
#              style="type",
#              errorbar="sd",
#              estimator="median",
#              inner=None,
#              markers=True,
#              dashes=False)

# plt.xticks(rotation=90)
# plt.title("Uncalibrated precursor mass error")
# plt.legend(bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)
# plt.savefig("QC_summary\\performance.pdf", dpi=300)


