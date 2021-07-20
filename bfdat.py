import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
import statistics
from scipy.optimize import leastsq

cat = pd.read_csv("C:\\Users\\Aman Arya\\Desktop\\astro\\open clusters\\catalouge\\1743new\\Trumpler_16.dat", delimiter='\s+')

def BF(data):
    ra = cat['RA_ICRS'].to_numpy()
    dec = cat['DE_ICRS'].to_numpy()
    d = 1/(cat['Plx'].to_numpy())
    print("number of stars:",len(d))
    def dist(a,b,c):
        d=[]
        m=[]
        for i in range(len(a)):
            d.append([])
            for j in range(len(b)):
                d[i].append(np.sqrt((a[j]-a[i])**2 + (b[j]-b[i])**2 + (c[j]-c[i])**2))
        m=np.matrix(d)
        m.sort()
        d1=m[:,1]
        d10=m[:,10]
        dd = np.array(d1/d10)
        dd = dd.astype('float64')
        dd = np.sort(dd, axis=None)
        return(dd)
    abc = dist(ra,dec,d)
    mean = statistics.mean(abc)
    sd = statistics.stdev(abc)
    #pdf = norm.pdf(abc, mean, sd)
    cdf = norm.cdf(abc, mean, sd)
    #cdff = norm.cdf(abc)
    #gauscdf = norm.pdf(cdf,np.mean(cdf),np.std(cdf))
    x=abc
    y=cdf
    def sinle_gaussian( x, params ):
        (c1, mu1, sigma1) = params
        res =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) ) 
        return res
    def sinle_gaussian_fit( params ):
        fit1 = sinle_gaussian( x, params )
        return (fit1 - y)
    
    # Least squares fit. Starting values found by inspection.
    #params for single gaussian fit 
        
    
    fit1 = leastsq( sinle_gaussian_fit, [0.15,0.2,0.1] )
    
    
    plt.plot( x, y, label='CDF',lw=3,c='r',alpha=0.8)
    plt.plot( x, sinle_gaussian( x, fit1[0] ), '--',label='Fit',c='g')
    plt.xlabel("d1/d10")
    plt.ylabel("CDF")
    plt.grid(True)
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()
    print(fit1[0])
    
    
    def double_gaussian( x, params ):
        (c1, mu1, sigma1, c2, mu2, sigma2) = params
        res =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) ) + c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )
        return res
    def double_gaussian_fit( params ):
        fit = double_gaussian( x, params )
        return (fit - y)
    # Least squares fit. Starting values found by inspection.
    # params for double gaussian fit 
    
    
    fit2 = np.array(leastsq( double_gaussian_fit, [0.15,0.2,0.1,0.6,0.4,0.2] ))
    
    plt.plot( x, y, label='CDF',lw=3,c='r',alpha=0.8)
    plt.plot( x, double_gaussian( x, fit2[0] ), '--',label='Fit',c='g')
    plt.xlabel("d1/d10")
    plt.ylabel("CDF")
    plt.grid(True)
    plt.legend(loc='best')
    #plt.xticks(np.arange(0,1.1,0.1))
    #plt.yticks(np.arange(0,0.2,0.02))
    plt.grid(True)
    #plt.axis([0,0.1,0,0.15])
    #plt.savefig("C:\\Users\\Aman Arya\\Desktop\\astro\\open clusters\\catalouge\\joshi1\\BF.png")
    plt.show()
    print(fit2[0])
    fit2 = abs(fit2)
    print('Binary Fraction is:',(fit2[0][3])/(fit2[0][0] + fit2[0][3]))

BF(cat)