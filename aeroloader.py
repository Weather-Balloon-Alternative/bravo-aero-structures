import numpy as np
import pandas as pd

def loadaero(filepath):
    df = pd.read_excel(filepath)
    data = df.to_numpy()
    namesraw = data[:,0]
    names = np.array(namesraw)
    for name in names:
        if len(np.where(names == name)[0]) > 1:
            copy = np.where(names == name)[0]
            for i in range(1,len(np.where(namesraw == name)[0])):
                names[copy[i]] = namesraw[copy[i]] + f'{i}'

    output = {}
    for i in range(len(names)):
        # output[names[i]] = list(filter(lambda x: str(x) != 'nan', data[i,:]))
        output[names[i]] = [obj for obj in data[i,:] if str(obj) != 'nan']
    return output


k = loadaero('Aero_output/glider_v7.xlsx')
print(k)