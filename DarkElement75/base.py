def profile_is_valid(profile):
    #Ensure profile is complete and has no errors, i.e. is only made up of A or B characters.
    valid = True
    try:
        for c in profile:
            if c not in ["A", "B"]:
                valid = False
                break

        return valid

    except:
        #If profile is NaN, it's invalid
        return False


def get_profile_combinations(profiles):
    """
    Wrapper function for get_recursive_profile_combinations, so that 
        we can use it and don't have to use global variables, allowing us to have
        an encapsulated function call and return value. 
        
    This function will call get_recursive_profile_combinations with a
        newly initialized combinations list, and return the resulting
        complete combinations list.
    """
    combinations = []
    def get_recursive_profile_combinations(profiles, i=0, combination=[]):
        """
        Recursive function that makes use of the profiles list of lists given
            to create a list of lists of all profile string combinations in the combinations list.

        Will create each combination list through depth-first searching
            the tree resulting from the profiles list of lists. There may be quite a lot of combinations,
            since the number of total combinations is equal to the product of all non-zero sublist lengths
            in the profiles list. 
        """

        #Base Case
        if i == len(profiles):
            combinations.append(combination)
            return

        #Recursive Cases
        if len(profiles[i])==0:
            #Empty list case
            get_recursive_profile_combinations(profiles, i=i+1, combination=combination+[""])
        else:
            #Non-empty list case
            for profile in profiles[i]:
                get_recursive_profile_combinations(profiles, i=i+1, combination=combination+[profile])
    
    get_recursive_profile_combinations(profiles)
    return combinations
