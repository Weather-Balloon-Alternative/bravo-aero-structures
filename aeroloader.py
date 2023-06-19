import numpy as np
import pandas as pd

def loadaero(filepath,discard=True):
    """
    ARGS:
            filepath: string,           name or route to a specific .xlsx or .xlsm file.
    OUTPUTS:
            output: dictionary,         dictionary with all coefficients
    """

    #read the file to numpy
    df = pd.read_excel(filepath)
    data = df.to_numpy()

    if discard:
        #Discard all data until VSPAERO_Stab
        iddiscard = np.where(data == 'VSPAERO_Stab')[0][0]
        data = data[iddiscard::]
    

    
    #get a list of all names that are mentioned in the excel file
    namesraw = data[:,0]
    names = np.array(namesraw)

    
    #ensure that all elements in the dictionary will get a unique name
    for name in names:
        if len(np.where(names == name)[0]) > 1:

            #find the locations of the copies
            copy = np.where(names == name)[0]

            #change names to specific names
            for i in range(1,len(np.where(namesraw == name)[0])):
                names[copy[i]] = namesraw[copy[i]] + f'{i}'

    #create output dictionary
    output = {}

    #initialize creation of elements in ouput by making range with indicies
    for i in range(len(names)):
        # output[names[i]] = list(filter(lambda x: str(x) != 'nan', data[i,:]))

        # create for every element a list with all column data in case the data itself is not equal to not a number.
        output[names[i]] = [obj for obj in data[i,:] if str(obj) != 'nan']
    return output