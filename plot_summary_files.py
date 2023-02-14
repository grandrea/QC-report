import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob

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
    files_df["date"] = pd.to_datetime(files_df["raw_date"], format='%y%m%d')
    files_df["FAIMS"] = files_df["Raw file"].apply(checkFAIMS)
    files_df["Detector"] = files_df["Raw file"].apply(checkDetector)
    files_df = files_df.reset_index()
    try:
        files_df = files_df.drop("index")
    except KeyError:
        pass
    return files_df

def NumberOfProteins(df):
    '''extract number of proteins from ProteinGroups file'''
    #filter out reverse, contaminant, identified by site.
    df = df[df["Reverse"]!="+"]
    df = df[df["Potential contaminant"]!="+"]
    df = df[df["Only identified by site"]!="+"]
    number_of_proteins = len(df)
    return number_of_proteins

#grab protein groups files-------------

proteinGroups_list = glob.glob("*HeLa*/*/txt/*proteinGroups.txt")

proteinGroups_df = CleanUpList(proteinGroups_list)

proteins = []
for element, rawfile in proteinGroups_df.iterrows():
    proteinGroups = pd.read_csv(rawfile["Raw file"], delimiter="\t")
    proteins.append(NumberOfProteins(proteinGroups))

proteinGroups_df["proteins"] = proteins


#grab summary files---------------
#needs MQ runs with a single raw file
summaryFiles_list =  glob.glob("*HeLa*/*/txt/*summary.txt")
summaryFiles_df = CleanUpList(summaryFiles_list)

for element, rawfile in summaryFiles_df.iterrows():
    if element==0:
        summary = pd.read_csv(rawfile["Raw file"], delimiter="\t")
        if len(summary)==2:
            summary = summary[summary["Raw file"]!="Total"]
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
precursor_error = []
precursor_sd= []
for element, rawfile in msmsFiles_df.iterrows():
    msms = pd.read_csv(rawfile["Raw file"], delimiter="\t")
    

            

summary["Raw file"] = summary["Raw file"].str.replace(".*\\\combined_","", regex=True).str.replace("raw.*", "raw")
proteinGroups_df["Raw file"] = proteinGroups_df["Raw file"].str.replace(".*\\\combined_","", regex=True).str.replace("raw.*", "raw")

#merge summary and proteinGroups statistics----------------------
QC_df = pd.merge(right=proteinGroups_df,
         left=summary,
         on=["Detector", "FAIMS", "date", "Raw file"], 
         how="outer")

QC_df["type"] = QC_df["FAIMS"]+"_"+QC_df["Detector"]

plt.close()
plt.clf()
sns.lineplot(data=QC_df, x="date",
             y="proteins",
             hue="type",
             markers=True,
             style="type",
             dashes=False)

plt.xticks(rotation=90)
plt.title("Identified proteins")

plt.close()
plt.clf()
sns.lineplot(data=QC_df, x="date",
             y="MS/MS submitted",
             hue="type",
             markers=True,
             style="type",
             dashes=False)

plt.xticks(rotation=90)
plt.title("MS2 submitted")

plt.close()
plt.clf()
sns.lineplot(data=QC_df, x="date",
             y="MS/MS identified [%]",
             hue="type",
             markers=True,
             style="type",
             dashes=False)

plt.xticks(rotation=90)
plt.title("MSMS identification rate")
