import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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


summary_table = pd.read_csv("QC_summary.csv")


#extract date from file names like  E230105_02_QC_AG_100_HeLa200ng_OT_FAIMS-45_inj1
summary_table["raw_date"] = summary_table["Raw file"].str.split("_").str[0].str[1:]
summary_table["date"] = pd.to_datetime(summary_table["raw_date"], format='%y%m%d')
summary_table.drop("raw_date", axis=1)

#annotate if exeriment is done with faims and orbitrap/iontrap
summary_table["FAIMS"] = summary_table["Raw file"].apply(checkFAIMS)
summary_table["Detector"] = summary_table["Raw file"].apply(checkDetector)

#keep only experiments at 200ng
summary_table = summary_table[summary_table['Raw file'].str.contains('200ng')]


plt.clf()
sns.lineplot(data=summary_table, x="date", y="MS/MS identified [%]", hue="Detector")
plt.show()
