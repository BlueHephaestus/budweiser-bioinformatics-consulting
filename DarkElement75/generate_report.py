"""
Generates a report for the counts of all LINE_NMs and their associated PROFILEs
    for every SAMPLE_ID, as obtained from analyzing and parsing GenoResults.csv, 
    and as described both in "SOP for Data Conversion and Preparation of Purity Report.docx" 
    and through our direct discussions with Audrey from Budweiser.

-Blake Edwards / Dark Element
"""

import pandas as pd
import numpy as np
#For now we don't do anything with standards since we don't have a standards key file

#Only read in these columns
data = pd.read_csv("GenoResults.csv", usecols=["DNA_PLATE", "WELL", "SAMPLE_ID", "PLANT_ID", "LINE_NM", "MARKER_NM", "Call"])

#Remove all rows with empty entries
data = data.dropna()

#Get unique column values
marker_nms = np.sort(np.unique(data["MARKER_NM"]))
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
Create a mapping from PROFILE:LINE_NM using instances of each PROFILE for which PROFILE is
    only As and Bs, i.e. the PROFILE is complete and has no errors.
    
For every LINE_NM entry, get the associated PROFILE if it has no errors.

Since there are often several errors in PROFILE values, we have to find the true PROFILE for each LINE_NM
    value. We do this by counting the instances of each non-error PROFILE for each LINE_NM,
    and then obtaining the true profile via the PROFILE which occurs most often.
    
Since often times the counts will be similar to [893, 2, 1, 2, 1] due to most being correct, but
    there still being a few errors, in the above case 6 errors.
"""
#Mapping dictionary to find the true PROFILE for each LINE_NM value.
line_nms = {}

#Create these ahead of time to avoid checks in the loop
for line_nm in np.unique(data["LINE_NM"]):
    line_nms[line_nm] = {}

def profile_is_valid(profile):
    #Ensure profile is complete and has no errors, i.e. is only made up of A or B characters.
    valid = True
    for c in profile:
        if c not in ["A", "B"]:
            valid = False
            break
    return valid
            
for i, row in data.iterrows():
    profile = row["PROFILE"]
    
    #If it's not valid, move on to the next one.
    if not profile_is_valid(profile):
        continue

    #Otherwise, increment this profile's count in its respective line_nm entry in our line_nms dictionary,
    #setting its count to 1 if it doesn't already exist.
    line_nm = row["LINE_NM"]
    if profile in line_nms[line_nm]:
        line_nms[line_nm][profile] +=1
    else:
        line_nms[line_nm][profile] = 1
        

"""
Get the true profile for each line_nm, and insert the true profile into a new mapping dictionary,
    made to map from profile:line_nm for each sample_id value.
    
Since it will also be used to count the number of each true profile (and failures) for each sample_id,
    we also duplicate this mapping across each unique SAMPLE_ID and store a count value for each mapping as well, 
    initialized to 0.
    
So it will be of the form sample_id:profile:{LINE_NM:line_nm, COUNT:count}
"""
        
#Mapping dictionary for all profiles and sample_ids
profiles = {}

#Add dictionaries for all unique sample_ids to avoid checks in the loop,
#And to add an entry to count failures.

for sample_id in sample_ids:
    profiles[sample_id] = {"FAIL":{"COUNT":0}}
    
#Add true profile : line_nm, count mapping to every sample_id in the profiles dict.
for line_nm,profile_counts in line_nms.items():
    true_profile = max(profile_counts, key=(lambda key: profile_counts[key]))
    for sample_id in sample_ids:
        profiles[sample_id][true_profile] = {"LINE_NM":line_nm, "COUNT":0}

#for each row in our data:
for i, row in data.iterrows():
    sample_id = row["SAMPLE_ID"]
    profile = row["PROFILE"]
    #if the profile has an error:
    if not profile_is_valid(profile):
        #we increment the FAIL count in the profiles dict for this row's sample_id and continue
        profiles[sample_id]["FAIL"]["COUNT"]+=1
        continue
    
    #if the profile is not in our profiles dict for this row's sample_id:
    if profile not in profiles[sample_id]:
        #we increment the FAIL count in the profiles dict for this row's sample_id and continue
        profiles[sample_id]["FAIL"]["COUNT"]+=1
        continue
    
    #if the line_nm the profile maps to is not the same as the line_nm for this row:
    if profiles[sample_id][profile]["LINE_NM"] != row["LINE_NM"]:
        #we increment the FAIL count in the profiles dict for this row's sample_id and continue
        profiles[sample_id]["FAIL"]["COUNT"]+=1
        continue
        
    #Otherwise, we increment the count in the profiles dict for this profile and sample id
    profiles[sample_id][profile]["COUNT"]+=1
    
#Get percentages of each COUNT relative to the total counts for each sample_id, 
#    and add this entry alongside each COUNT
for sample_id in sample_ids:
    #Get sum of counts in this sample_id
    count_sum = 0
    for profile in profiles[sample_id]:
        count_sum += profiles[sample_id][profile]["COUNT"]
    
    #Use sum of counts to add PERCENTAGE attribute
    for profile in profiles[sample_id]:
        profiles[sample_id][profile]["PERCENTAGE"] = profiles[sample_id][profile]["COUNT"]/count_sum*100
        
#Generate a report dataframe with our complete profiles dictionary
report = pd.DataFrame(columns=["SAMPLE_ID", "LINE_NM", "PROFILE", "COUNT", "PERCENTAGE"])
for sample_id in sample_ids:
    for profile,profile_data in sorted(profiles[sample_id].items()):
        report = report.append({"SAMPLE_ID":sample_id, 
                       "LINE_NM":profile_data["LINE_NM"] if "LINE_NM" in profile_data else "FAIL", 
                       "PROFILE":profile, 
                       "COUNT":profile_data["COUNT"], 
                       "PERCENTAGE":profile_data["PERCENTAGE"]}, ignore_index=True)
        
#Order the columns for data csv before writing
data_cols = ["SAMPLE_ID", "DNA_PLATE", "WELL", "PLANT_ID", "LINE_NM"]
data_cols.extend([marker_nm for marker_nm in marker_nms])
data_cols.append("PROFILE")

#With both our report dataframe and condensed data dataframe finished, we write both to csvs.
report.to_csv("automated_report.csv", index=False)
data.to_csv("automated_supplemental_data.csv", columns=data_cols, index=False)
