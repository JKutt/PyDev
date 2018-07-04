import numpy as np
import scipy
from scipy import sparse 

###############################################################################
# STACK FUNCTION
###############################################################################
def interpts(t,xt,time):
    cnt = -1
    for i in range(xt.size):
        cnt += 1
        if t[i] == time:
            return xt[i]
            break
        if cnt > 0:
            if t[i-1] < time < t[i]:
                interpData = (xt[i-1]+xt[i])/2
                return interpData
                break
        

# calculate stacked data
def stack(xt,freq,sampleRate):
    # determine time length of the recording
    rectime = xt.size/sampleRate
    #determine Period T
    T = 1/freq
    #determine how many periods are present int he signal
    numT = np.floor(rectime/T)
    # determine samples per period and half period
    sampPerT = T*sampleRate
    sampPerHalfT = sampPerT/2.0
    # now determine number of available stacks
    nStacks = np.floor(numT*2)

    tmp = np.ones(int(nStacks-4))
    tmp = tmp*4;
    # create the +ve/-ve
    for i in range(tmp.size):
        tmp[i] = tmp[i]*(pow((-1),(i+2)))

    # create full filter kernal
    f1 = np.zeros(tmp.size+4)
    f1[0] = 1
    f1[1] = -3
    f1[f1.size-2] = 3*(pow((-1),nStacks-2))
    f1[f1.size-1] = pow((-1),nStacks-1)

    for j in range(tmp.size):
        f1[j+2] = tmp[j]

    f1 = f1/(4.0*(nStacks-2))
    sizef2 = int(sampPerHalfT)
    f2 = np.zeros((sizef2,sizef2))
    for k in range(sizef2):
        f2[k,k] = 1
    
    filterK = sparse.kron(f1,f2).todense()
    print(filterK)
    print(filterK.shape)
    print(sampPerHalfT)

    # trim the data to proper length
    trimdata = np.zeros(int(filterK.shape[1]))
    for idx in range(int(filterK.shape[1])):
        trimdata[idx] = xt[idx]

    print(trimdata.size)
    print(filterK[599,1])
    # now do the stack calculation
    stackD = np.matmul(filterK,trimdata)

    stackedCheda = np.zeros(stackD.size)
    for idx in range(stackD.size):
        stackedCheda[idx] = stackD[0,idx]

    # now lets extract the value of the current
    # sum over 50% to 90% of the on-time
    sumStart = int(0.5*(sampPerHalfT/2))
    sumEnd = int(0.9 * (sampPerHalfT/2))
    amp = np.sum(stackedCheda[sumStart:sumEnd])/(sumEnd-sumStart)  
    #return the amplitude
    return stackedCheda
    
