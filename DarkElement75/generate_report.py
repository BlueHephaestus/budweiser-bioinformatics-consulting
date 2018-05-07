"""
Generates a report for the counts of all LINE_NMs and their associated PROFILEs
    for every SAMPLE_ID, as obtained from analyzing and parsing GenoResults.csv, 
    and as described both in "SOP for Data Conversion and Preparation of Purity Report.docx" 
    and through our direct discussions with Audrey from Budweiser.

WARNING: Should be run AFTER generate_key.py generates the key.csv file in the same directory.

-Blake Edwards / Dark Element
"""

import pandas as pd
import numpy as np
from base import *

#Only read in these columns
data = pd.read_csv("GenoResults.csv", usecols=["DNA_PLATE", "WELL", "SAMPLE_ID", "PLANT_ID", "LINE_NM", "MARKER_NM", "Call"])

#Remove all rows with empty entries
data = data.dropna()

#Get unique column values
line_nms = np.unique(data["LINE_NM"])
marker_nms = np.unique(data["MARKER_NM"])
sample_ids = np.unique(data["SAMPLE_ID"])

#Sort so that we group the same samples together, as well as order them.
data = data.sort_values(["SAMPLE_ID", "PLANT_ID", "MARKER_NM"])

#Add columns for each marker_nm value
for marker_nm in marker_nms:
    data[marker_nm] = ""

"""
We have groups of which share the same line_nm and plant_id in this table. 

For each of these groups, we collapse them so that each row's unique call value is put in the just-created columns
    for each marker_nm value. E.g. for rows that have 4:AA, 3:BB, 5:CC where 3, 4, 5 are marker_nm values
    and AA, BB, and CC are call values, we'd then combine these so that we had one result row like this:
    
    4  |  3  |  5
    AA    BB    CC
    
    Where 4, 3, and 5 are the new headings. So we're taking these call values from the related segments 
        and condensing theminto one row for each of these segments.
        
Though the segments may all be the same length as the number of marker_nm values, e.g. in the above example this
    length was 3, this isn't always the case, so we iterate through each row and know we're done with 
    the segment once we've filled out all the values. Once we're done with this segment, we increment the index
    for the destination row.
"""
dst = 0
marker_nm_loc = data.columns.get_loc("MARKER_NM")
call_loc = data.columns.get_loc("Call")
plant_id_loc = data.columns.get_loc("PLANT_ID")
well_loc = data.columns.get_loc("WELL")

#Loop through all individual rows
for src in range(len(data)):
    
    #If there are not any empty elements in the dst row, then we're done and can move on 
    if not np.any(data.iloc[dst] == ""):
        #New group, increment dst row index
        dst+=1
        
        #Set new row to match group info
        data.iloc[dst] = data.iloc[src]
        
    #Insert new call value in corresponding marker_nm column
    data.iat[dst, data.columns.get_loc(data.iat[src, marker_nm_loc])] = data.iat[src, call_loc]

#Remove all remaining original rows, i.e. the rows used for condensing after the last complete condensed row.
if not np.any(data.iloc[dst] == ""):
    #Include last row since it's complete
    data = data[:dst]
else:
    #Don't include the last row since it's incomplete
    data = data[:dst-1]
    
#MARKER_NM and Call rows are no longer necessary since they were used for condensing, so we remove these as well.
data = data.drop(["MARKER_NM", "Call"], axis=1)

#Create a new PROFILE column for the concatenation of the marker_nm columns.
data["PROFILE"] = data[[marker_nm for marker_nm in marker_nms]].apply(lambda x: ''.join(x), axis=1)

"""
We load our key.csv file from the previous generate_key.py script, made to map each LINE_NM
    value to each possible PROFILE it may have.
"""
try:
    profiles_key = pd.read_csv("key.csv")
except:
    print("YOU MUST FIRST RUN `python3 generate_key.py` TO GENERATE A key.csv FILE")
    
"""
We create a new mapping dictionary, to keep track of the number of times
    each LINE_NM maps to each of it's possible PROFILEs as described by the profiles_key
    for each SAMPLE_ID, as well as to keep track of the number of failures
    for each SAMPLE_ID.
    
We have a failure whenever the PROFILE is not valid, i.e. is not AA, BB, or AB,
    or whenever the PROFILE is not one of the profiles associated with the LINE_NM.

This mapping dictionary will be of the form:
    sample_id:line_nm:{profile1:count1, profile2:count2}
"""
#Mapping dictionary for all profiles and sample_ids
profiles = {}

#Add dictionaries for all unique sample_ids to avoid checks in the loop,
#And to add an entry to count failures.
for sample_id in sample_ids:
    profiles[sample_id] = {"FAIL":0}

#Additionally add dictionaries for all unique line_nms to avoid checks in the loop
for sample_id in sample_ids:
    for line_nm in np.unique(profiles_key["LINE_NM"]):
        profiles[sample_id][line_nm] = {}

#Add the line_nm:{profile1:count1, profile2:count2,...} mapping for every sample_id in the profiles dict
for i, row in profiles_key.iterrows():
    
    #Create profile only from the MARKER_NM values present in our data.
    profile = ""
    for marker_nm in marker_nms:
        profile+=row[marker_nm]
    
    #Insert completed profile with 0 initial count for every sample_id
    for sample_id in sample_ids:        
        profiles[sample_id][row["LINE_NM"]][profile] = 0    

#for each row in our data:
for i, row in data.iterrows():
    sample_id = row["SAMPLE_ID"]
    line_nm = row["LINE_NM"]
    profile = row["PROFILE"]
    #if the profile has an error:
    if not profile_is_valid(profile):
        #we increment the FAIL count in the profiles dict for this row's SAMPLE_ID and continue
        profiles[sample_id]["FAIL"]+=1
        continue
    
    #If the PROFILE is not one of the profiles associated with the LINE_NM for this row's SAMPLE_ID:
    if profile not in profiles[sample_id][line_nm]:
        #we increment the FAIL count in the profiles dict for this row's SAMPLE_ID and continue
        profiles[sample_id]["FAIL"]+=1
        continue
        
    #Otherwise, we increment the count in the profiles dict for this LINE_NM's PROFILE and SAMPLE_ID
    profiles[sample_id][line_nm][profile]+=1

#Get percentages of each count relative to the total counts for each sample_id, 
#    and replace each count with a (count, percentage) tuple
for sample_id in sample_ids:
    #Get total sum of all profile counts in this sample_id
    count_sum = profiles[sample_id]["FAIL"]
    for line_nm in line_nms:
        for profile in profiles[sample_id][line_nm].keys():
            count_sum += profiles[sample_id][line_nm][profile]
    
    #Use sum of counts to replace counts with (count, percentage) tuples
    for line_nm in line_nms:
        for profile in profiles[sample_id][line_nm].keys():
            profiles[sample_id][line_nm][profile] = (profiles[sample_id][line_nm][profile], 
                                                     profiles[sample_id][line_nm][profile]/count_sum*100)
    #Add fails
    profiles[sample_id]["FAIL"] = (profiles[sample_id]["FAIL"], profiles[sample_id]["FAIL"]/count_sum*100)
        
#Generate a report dataframe with our complete profiles dictionary
report = pd.DataFrame(columns=["SAMPLE_ID", "LINE_NM", "PROFILE", "COUNT", "PERCENTAGE"])
for sample_id in sample_ids:
    for line_nm in line_nms:
        for profile,profile_data in profiles[sample_id][line_nm].items():
            report = report.append({"SAMPLE_ID":sample_id, 
                           "LINE_NM":line_nm,
                           "PROFILE":profile, 
                           "COUNT":profile_data[0], 
                           "PERCENTAGE":profile_data[1]}, ignore_index=True)
    #Add fails
    report = report.append({"SAMPLE_ID":sample_id, 
                   "LINE_NM": "FAIL",
                   "PROFILE":"", 
                   "COUNT":profiles[sample_id]["FAIL"][0], 
                   "PERCENTAGE":profiles[sample_id]["FAIL"][1]}, ignore_index=True)
        
#Order the columns for data csv before writing
data_cols = ["SAMPLE_ID", "DNA_PLATE", "WELL", "PLANT_ID", "LINE_NM"]
data_cols.extend([marker_nm for marker_nm in marker_nms])
data_cols.append("PROFILE")

#With both our report dataframe and condensed data dataframe finished, we write both to csvs.
report.to_csv("report.csv", index=False)
data.to_csv("supplemental.csv", columns=data_cols, index=False)
