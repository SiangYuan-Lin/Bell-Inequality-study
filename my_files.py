import numpy as np # linear algebra
import pandas as pd
filepath = '../Updated signal features/npz files/'
def load_particle(file, index, unit):
    DF = pd.DataFrame()
    if index == None:
        f = np.load(filepath+file)
        for i in ['E', 'px', 'py', 'pz', 'm', 'pt']:
            print('Loaded '+ i +' from '+file+' collection ',f.files, 'rows: ', f[i].shape[0])
            DF[i] = f[i]*unit
        for i in ['eta', 'phi']:#, 'y']:
            print('Loaded '+ i +' from '+file+' collection ',f.files, 'rows: ', f[i].shape[0])
            DF[i] = f[i]
    else:
        f = np.load(filepath+file)
        for i in ['E', 'px', 'py', 'pz', 'm', 'pt']:
            print('Loaded '+ i +' from '+file+' collection ',f.files, 'rows: ', len(index))
            DF[i] = f[i][index]*unit
        for i in ['eta', 'phi']:#, 'y']:
            print('Loaded '+ i +' from '+file+' collection ',f.files, 'rows: ', len(index))
            DF[i] = f[i][index]
    return DF
def load_evt(file, label, index):
    if index == None:
        f = np.load(filepath+file)
        print('Loaded '+label+' from '+file+' collection ',f.files, 'rows: ', f[label].shape[0])
        return f[label]
    else:
        f = np.load(filepath+file)
        print('Loaded '+label+' from '+file+' collection ',f.files, 'rows: ', len(index))
        return f[label][index]

def load_p4(file,labels,index):
    labels = ['E','px','py','pz'] if labels == None else labels
    if index == None:
        DF = pd.DataFrame(columns= labels)
        f = np.load(filepath+file)
        DF[labels[0]] = f['E' ]
        DF[labels[1]] = f['px']
        DF[labels[2]] = f['py']
        DF[labels[3]] = f['pz']
        print('Loaded E, px, py, pz from '+file+' collection ',f.files, 'rows: ', DF.iloc[:,0].shape[0])
        return DF
    else:
        DF = pd.DataFrame(columns= labels)
        f = np.load(filepath+file)
        DF[labels[0]] = f['E' ][index]
        DF[labels[1]] = f['px'][index]
        DF[labels[2]] = f['py'][index]
        DF[labels[3]] = f['pz'][index]
        print('Loaded E, px, py, pz from '+file+' collection ',f.files, 'rows: ', len(index))
        return DF