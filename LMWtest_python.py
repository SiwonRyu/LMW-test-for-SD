import operator as op
import numpy as np
import pandas as pd
import math
import time
import xlrd

# Importing data part
data = pd.read_excel ('C:\data_temp\ex.xlsx') #for an earlier version of Excel, you may need to use the file extension of 'xls'
x1 = np.array(data)[:,0]
x2 = np.array(data)[:,1]
x3 = np.array(data)[:,2]
x4 = np.array(data)[:,3]

# Parameters
[alpha,n1,n2,ngrid,R,gen] = [0.05,100,100,40,100,2]

# Generate or import : when gen = 1, use following desing and when gen <> 1, real imported data(x1 and x2) will be used
if gen == 1 :
    # Design
    [mu1, sigma1, mu2, sigma2] = [0, 1, 2, 1]
    sample1 = mu1 + sigma1 * np.random.randn(n1,1,1,R)
    sample2 = mu2 + sigma2 * np.random.randn(n2,1,1,R)

    minx = min(sample1.min(3).min(0), sample2.min(3).min(0))[0]
    maxx = max(sample1.max(3).max(0), sample2.max(3).max(0))[0]
else :
    sample1 = x1[:,None,None,None] #This dimension adjustement is for fix dimension in following functions
    sample2 = x2[:,None,None,None]

    minx = min(sample1.max(0), sample2.max(0))[0][0]
    maxx = max(sample1.max(0), sample2.max(0))[0][0]
grid = np.linspace(minx,maxx, num=ngrid)
s = 1

b1 = 10
b2 = 10

def lmwtest(sample1,sample2,grid,s,b1,b2) :
    start_time = time.time()
    def operator(x, grid, s):
        return (np.apply_along_axis(op.lt, 1, x, grid) * np.apply_along_axis(op.sub, 1, x, grid) ** (
                s - 1) / math.factorial(s - 1))[:, :, 0]

    def ecdf(x, grid, s):
        return operator(x, grid, s).mean(0)[None]

    def stat(sample1, sample2, grid, s):
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        return (n1 * n2 / (n1 + n2)) ** (0.5) * (ecdf(sample1, grid, s) - ecdf(sample2, grid, s)).max(1)[:, None]

    def subsampling(sample, subsize, nsub):
        # nsub : # of subsamples
        # subindex : subsample size x 1 x nsub
        subindex = (np.arange(0, subsize, 1)[:, None] + np.arange(0, nsub, 1)[:, None].T)[:, None]

        # subsample : subsample size x 1 x nsub x 1
        subsample = sample[subindex][:, :, :, 0, 0]
        return subsample

    def substat(sample1, sample2, grid, s, b1, b2):
        # This part is to take into account different sample size (Linton(2005))
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        lbd = n2 / (n1 + n2)
        b = lbd * b1
        nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])
        subsample1 = subsampling(sample1, b1, nsub)
        subsample2 = subsampling(sample2, b2, nsub)
        return (b) ** (0.5) * (ecdf(subsample1, grid, s) - ecdf(subsample2, grid, s)).max(1)[:, None]

    lmw = stat(sample1, sample2, grid, s)
    pval = (substat(sample1, sample2, grid, s, b1, b2) > lmw).mean(2)[:, :, None]

    if s == 1:
        Torder = 'FSD'
    elif s == 2:
        Torder = 'SSD'
    elif s == 3:
        Torder = 'TSD'
    else:
        Torder = str(s) + 'th order SD'

    if R > 1 & gen == 1:
        print('LMW Test for Stochastic Dominance'
              '\n* H0 : sample1', Torder, 'sample2')
        print('* Design : \n\tsample1 ~ N(%3.1f' % mu1, ',%3.1f)' % sigma1,
              '\n\tsample2 ~ N(%3.1f' % mu2, ',%3.1f)\n' % sigma2)
        print('* #(sample1) \t\t = %6d' % sample1.shape[0],
              '\n* #(sample2) \t\t = %6d' % sample1.shape[0])
        print('* #(subsample1) \t = %6d' % b1,
              '\n* #(subsample2) \t = %6d' % b2)
        print('\n* SD order \t\t\t = %6d' % s,
              '\n* # of grid points \t = %6d' % ngrid,
              '\n* # of repetition \t = %6d\n' % R)
        print('* Rejection probabilities = %5.4f' % (pval < alpha).mean(3))


    else:
        print('LMW Test for Stochastic Dominance'
              '\n* H0 : sample1', Torder, 'sample2')
        if gen == 1:
            print('* Design : \n\tsample1 ~ N(%3.1f' % mu1, ',%3.1f)' % sigma1,
                  '\n\tsample2 ~ N(%3.1f' % mu2, ',%3.1f)' % sigma2)
        print('* #(sample1) \t\t = %6d' % sample1.shape[0],
              '\n* #(sample2) \t\t = %6d' % sample1.shape[0])
        print('* #(subsample1) \t = %6d' % b1,
              '\n* #(subsample2) \t = %6d' % b2)
        print('\n* SD order \t\t\t = %6d' % s,
              '\n* # of grid points \t = %6d\n' % ngrid)
        print('* LMW statistic \t\t = %5.4f' % lmw)
        print('* p-value(subsampling) \t = %5.4f' % pval)

    et = time.time() - start_time
    print('\n* Time elapsed : %5.2f Sec' % et)

lmwtest(sample1,sample2,grid,s,b1,b2)