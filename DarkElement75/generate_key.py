"""
Generate a master key file, mapping each LINE_NM to a PROFILE value based on majority rules 
    from the StandardsData.xlsx file. This resulting file just has a header, and each line in
    the file is a different LINE_NM to PROFILE mapping. There may be multiple duplicate LINE_NM keys, 
    that map to various different profiles. We determine the profiles by choosing all profiles 
    that LINE_NM maps to, for which it maps to them >= 25% of it's mappings.
    
-Blake Edwards / Dark Element
"""
import pandas as pd
import numpy as np
import sys
from base import *

#Read in file, with all columns except for PROJECT_NM, DNA_PLATE, PLANT_ID
#We keep rows with empty entries, since they may have some valid entries.
if len(sys.argv)>1:
    data = pd.read_excel(sys.argv[1])
else:
    data = pd.read_excel("StandardsData.xlsx")
data = data.drop(["PROJECT_NM", "DNA_PLATE", "PLANT_ID"], axis=1)
data = data.rename({"LINES":"LINE_NM"}, axis=1)

#Get lists of all marker_nm and line_nm values
marker_nms = list(data.columns[1:])
line_nms = np.unique(data["LINE_NM"])

"""
We now get the number of occurrences of each valid call value, 
    for each marker_nm value,
    for each line_nm value,
    throughout the entire file.
    
We use a multi-tiered dictionary of the form marker_nm -> line_nm -> call -> count,
    since using the form line_nm -> marker_nm -> call -> count would result in a large amount
    of duplicate values. 
    
We'll then use this to find all line_nm: profile mappings. 
"""
#Create these ahead of time to avoid checks in the loop, initialize our multi-tiered dictionary
profiles = {}
for marker_nm in marker_nms:
    profiles[marker_nm] = {}
    for line_nm in line_nms:
        profiles[marker_nm][line_nm] = {}

for i, row in data.iterrows():
    for marker_nm in marker_nms:
        #If this profile is not valid, move on to the next one.
        profile = row[marker_nm]
        
        if not profile_is_valid(profile):
            continue

        #Otherwise, increment this profile's count in it's respective line_nm and marker_nm entry 
        #in our profiles dictionary, setting its count to 1 if it doesn't already exist.
        line_nm = row["LINE_NM"]    
        if profile in profiles[marker_nm][line_nm]:
            profiles[marker_nm][line_nm][profile] += 1
        else:
            profiles[marker_nm][line_nm][profile] = 1

"""
Go through now complete profiles dictionary for each line_nm entry, and for every marker_nm
    entry, remove all profiles which make up < 25% of the total entries for this marker_nm for this line_nm.
"""
for marker_nm in marker_nms:
    for line_nm in line_nms:
        #Get sum number of profile occurences for this line_nm and marker_nm
        n = sum(profiles[marker_nm][line_nm].values())
        
        #Loop through profiles and remove every one that makes up <25% of the total entries
        for profile in list(profiles[marker_nm][line_nm].keys()):
            if profiles[marker_nm][line_nm][profile]/n < .25:
                del profiles[marker_nm][line_nm][profile]
        
"""
Use our now complete profiles dictionary to create an nx2 profile_map list for each pair mapping
    LINE_NM to PROFILE values. We additionally create multiple entries for a LINE_NM 
    in the case of multiple possible PROFILE combinations.
    
We use a list so that we can more efficiently append and sort, and since there may
    be both duplicate LINE_NMs and duplicate PROFILEs (but not both at once) throughout the list,
    and attempting to use a dictionary for this would be difficult.
"""
profile_map = []
for line_nm in line_nms:
    #To be used to create all possible profile combinations for this line_nm
    combinations = []
    for marker_nm in marker_nms:
        #Add a list of all profiles this line_nm maps to with this marker_nm
        combinations.append(list(profiles[marker_nm][line_nm].keys()))
    
    #Add all profile combinations for this line_nm to the main profile_map
    for combination in get_profile_combinations(combinations):
        profile_map.append([line_nm, combination])
        
"""
At this point profile_map is our completed list, mapping each line_nm to its profile combination(s). 
    The profile combinations are represented as lists of strings so that we can still reference marker_nm,
    and since some entries may be empty.

We use this to create our final key dataframe and CSV.
"""

#Initialize this ahead of time since we know the size
data = pd.DataFrame("", index=np.arange(len(profile_map)), columns=["LINE_NM"]+marker_nms)

#Fill with our values
for i, row in enumerate(profile_map):
    data.iat[i, 0] = row[0]
    for j, marker_nm in enumerate(row[1]):
        data.iat[i,j+1] = marker_nm

#Make CSV
data.to_csv("key.csv", index=False)
