import pandas as pd
import numpy as np
train = pd.read_csv('train.csv')
def doHist(data):
    h = np.zeros(600)
    for j in np.ceil(data.values).astype(int):
        h[j:] += 1
    h /= len(data)
    return h
hSystole = doHist(train.Systole)
hDiastole = doHist(train.Diastole)
sub = pd.read_csv('sample_submission_validate.csv', index_col='Id')
N = len(sub)//2
X = np.zeros((2*N,600))
for i in range(N):
    X[2*i,:] = hDiastole
    X[2*i+1,:] = hSystole
sub[sub.columns] = X
sub.to_csv('submission.csv')
