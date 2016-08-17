# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 08:24:12 2016

@author: MVLSTECH
"""
from gaussfitter import gaussfit
#from spkfit import normdenoise
import  tifffile
import numpy as np
from matplotlib.pyplot import plot
from scipy import ndimage,signal,spatial
#from untitled7 import maskanddenoise
import pandas as pd
import easygui
import re
from skimage import restoration
from twodsparkdetection2 import gettimes,get_mask
from os import path
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

def getttp(a):
	i=a[::-1].argmax()
	ttp=0
	while (a[::-1][i]>a[::-1].max()*0.1) & (i<len(a)-1):
		i=i+1
		ttp=ttp+1
	return ttp

def ParseNumberList(s, item='', nmin=1, nmax=96, allowdup=False, sort=False, \
    debug=False):
    #Francis Burton has written this bit 
    if debug:
        print 'Input: ' + s
    numlist = []
    errstr = ''
    s = re.sub('^[^0-9]+', '', s)   # delete leading non-numbers
    s = re.sub('[^0-9]+$', '', s)   # delete trailing non-numbers
    s = re.sub('\s*-\s*', '-', s)   # delete whitespace in ranges
    s = re.sub('\s*to\s*', '-', s)  # ... and convert 'to' to hyphen
    s = re.sub('-+', '-', s)        # delete consecutive hyphens
    s = re.sub(',', ' ', s)         # convert commas to spaces
    s = re.sub(' +', ' ', s)        # compress multiple spaces to single
    if debug:
        print 'Numbers: ' + s
    
    for w in re.split('\s+', s):
        if '-' in w:
            r = re.split('-', w)
            if len(r) == 2:
                try:
                    n1 = int(r[0])
                except ValueError:
                    errstr = 'illegal first %snumber in range %s' % (item, w)
                    return ([], errstr)
                try:
                    n2 = int(r[1])
                except ValueError:
                    errstr = 'illegal second %snumber in range %s' % (item, w)
                    return ([], errstr)
            elif len(r) == 1:
                errstr = 'missing %snumber in range %s' % (item, w)
                return ([], errstr)
            else:
                errstr = 'extra %snumber in range %s' % (item, w)
                return ([], errstr)
            if n1 > n2:
                errstr = '%snumbers not ascending in range %s' % (item, w)
                n1, n2 = n2, n1
            if n1 < nmin:
                errstr = '%number %d less than min %d in range %s' \
                    % (item, n1, nmin, w)
            if n2 > nmax:
                errstr = '%number %d greater than maximum %d' \
                    % (item, n2, nmax, w)
            for i in xrange(n1, n2+1):
                numlist.append(i)    
        else:                       # single number
            try:
                n = int(w)
            except ValueError:
                errstr = 'illegal %snumber %s' % (item, w)
                return ([], errstr)
            if n < nmin:
                errstr = '%snumber %d less than min %d' % (item, n, nmin)
            if n > nmax:
                errstr = '%snumber %d greater than max %d' % (item, n, nmax)
            numlist.append(n)
            
        if not allowdup:
            if len(numlist) != len(set(numlist)):
                errstr = '%snumber list contains duplicate(s): %s' \
                    % (item, s)
    
        if sort:
            numlist = sorted(numlist)
            
    return (numlist, errstr)


def cutspark(im,pkpt,xysize=15,zsize=20):
   
    pkpt=np.array(pkpt,np.int64)
    halfsize=int(xysize//2  )
    halfzsize=zsize
    t=max([0,pkpt[1]-halfsize])
    b=min([im.shape[1],pkpt[1]+halfsize+1])
    l=max([0,pkpt[2]-halfsize])
    r=min([im.shape[2],pkpt[2]+halfsize+1])
    print t,b,l,r
    test=im[pkpt[0]-halfzsize:pkpt[0]+halfzsize+1,t:b,l:r]
    return test

def fitspk(im,pkpt,xysize=15,interval=2.5,pixsize=0.245,plotit=False,circle=False,excludebig=False):
    test=cutspark(im,pkpt,xysize=xysize,zsize=220/interval)
    tprof=signal.savgol_filter(test.mean(axis=1).mean(axis=1),7,2)
    z=np.arange(tprof.size)[tprof>tprof.max()/2]
    a=[]
    if plotit:
        fig, ax = plt.subplots()
    for i in z-1:
        com=gaussfit(test[i],circle=circle)[2:4]
        a.append(com)  
        if plotit:            
            plt.imshow(test[i], origin="lower")
            plt.scatter(*com[::-1])
            if plotit:
                ax.annotate(str(i), (com[1],com[0]))
    if plotit:
        plt.show()    
    x=np.array(a)[:,1]
    y=np.array(a)[:,0]
    x=x*pixsize
    y=y*pixsize
    #exclude very large movements. They are probably false
    if excludebig:
        include=np.arange(x.size-1)[np.abs(np.diff(x))<0.5]
    else:
        include=np.arange(x.size-1)
    x=x[include]
    y=y[include]
    z=z[include]
    return x,y,z
    
def fitspk2(im,pkpt,xysize=15,interval=2.5,pixsize=0.245,plotit=False,plotitg=False):
    #use top 80% as threshold
    #mask off rest of image
    #use centroid method to find middle of spark with subpixel precision
    test=cutspark(im,pkpt,xysize=xysize,zsize=220/interval)
    tprof=signal.savgol_filter(test.mean(axis=1).mean(axis=1),7,2)
    z=np.arange(tprof.size)[tprof>tprof.max()/2]
    if plotit or plotitg:
        fig, ax = plt.subplots()
    coms=[]
    gcoms=[]
    for i in z:
        hist, bins = np.histogram(test[i].ravel(), normed=True, bins=100)
        threshold = bins[np.cumsum(hist) * (bins[1] - bins[0]) > 0.8][0]
        mnorm2d = np.ma.masked_less(test[i],threshold)
        com = ndimage.measurements.center_of_mass(mnorm2d)
        mnorm2d = np.where(test[i]<threshold,threshold,test[i])
        gcom = gaussfit(mnorm2d,center=True)[2:4]
        if plotit:            
            if i==z[-1]:
                plt.imshow(test[i], origin="lower")
            plt.scatter(*com[::-1])
        if plotitg:            
            if i==z[-1]:
                plt.imshow(test[i], origin="lower")
            plt.scatter(*gcom)
        coms.append(com)
        gcoms.append(gcom)
    if plotit:
        for i,v in enumerate(coms):
            ax.annotate(str(i), (v[1],v[0]))
        plt.show()
    if plotitg:
        for i,v in enumerate(gcoms):
            ax.annotate(str(i), (v[0],v[1]))
        plt.show()
    y=np.array(coms)[:,0]*pixsize
    x=np.array(coms)[:,1]*pixsize
    gx=np.array(gcoms)[:,0]*pixsize
    gy=np.array(gcoms)[:,1]*pixsize
    return x,y,z,gx,gy

def ttp(d,pkpt,interval,zoomf=2,size=30):
    prof=cutspark(d,pkpt,zsize=size).mean(axis=1).mean(axis=1)
    
    noisefloor=prof[0:5].std()*3
    prof=ndimage.zoom(prof,zoomf)
    j=prof.argmax()
    while prof[j]>noisefloor:
        j=j-1
    return (prof.argmax()-j)/np.float32(zoomf)

def fdhm(d,pkpt,interval,zoomf=2,size=30):
    prof=cutspark(d,pkpt,zsize=size).mean(axis=1).mean(axis=1)
    baseline=prof[0:5].mean()
    prof=ndimage.zoom(prof,zoomf)
    ipk=prof.argmax()
    j=ipk
    hm=(prof[j]-baseline)/2.0
    while prof[j]>hm:
        j=j-1
    before =ipk-j
    j=ipk
    while prof[j]>hm:
        j=j+1
    after=j-ipk
        
    return (after+before+1)/np.float32(zoomf)

def strlisttoarray(ps):
    return (np.float32(ps.split('.')[0].split('[')[1]),np.float32(ps.split('.')[1]),np.float32(ps.split('.')[2]))

#   
def msd(xdata,ydata):
    r = np.sqrt(xdata**2 + ydata**2)
    diff = np.diff(r) #this calculates r(t + dt) - r(t)
    diff_sq = diff**2
    MSD = np.mean(diff_sq)
    return MSD
    
def sstrobrem(im,cellmask):
    smim=(im/ndimage.gaussian_filter(im,(0,5,5)))
    ncellmask=np.zeros(im.shape,np.bool)
    ncellmask[:]=cellmask
    ntest=np.where(cellmask==1,smim,np.NaN)
    corrim=im.copy()
    ntestmean=np.nanmean((ntest-1),axis=2)
    for i in range(im.shape[0]):
        corrim[i]=np.array([(-1*ntestmean[i]+1)]*im.shape[2]).transpose()*im[i]
    return corrim    

def getimangle_corr(im):
    l=[]
    a=ndimage.gaussian_filter(im,2)-ndimage.gaussian_filter(im,7)
    for i in np.arange(0,10,0.1):
        l.append(np.abs(ndimage.rotate(a,i).mean(axis=1)).sum())
    return np.array(l).argmax()*0.1

def getimangle(im):
    angles=[]
    for i in range(10):
        angles.append(getimangle_corr(im[i]))
    return (np.median(np.array(angles)))
    
if __name__ == '__main__':      
    #fn='C:/Users/MVLSTECH/Desktop/rattest3/Stream-6b.tif'
    #fn='f:/data/connor/sparks\\2msafter1hzcultured24hrs-4.tif'
    fn=easygui.fileopenbox()
    dffilename=(fn.lower()).replace('.tif','results.csv')
    df=pd.read_csv(dffilename)
    
    
    cellmaskname=(fn.lower()).replace('.tif','cellmask.tif.npy')
    
    if path.exists(cellmaskname):
        cellmask=np.load(cellmaskname)
    else:
        cellmask=get_mask(im.mean(axis=0))
    
    smname=(fn.lower()).replace('.tif','sm.tif.npy')
    if path.exists(smname):
        d=np.load(smname)
        l=np.array([cellmask]*d.shape[0])
    else:
        im=np.float32(tifffile.imread(fn))
        rotatev=getimangle(im)
        rotatev=np.float(easygui.enterbox('rotate by',default=str(rotatev)))
        im=ndimage.rotate(im,rotatev,axes=(2,1),reshape=False,mode='nearest')
        if rotatev!=0:
            im=sstrobrem(im,cellmask)
        immean=im.mean(axis=0)[cellmask==1].mean()
        imstd=im.std(axis=0)[cellmask==1].mean()
        l=np.array([cellmask]*im.shape[0])
        im=np.where(l==1,im,(np.random.randn(*(im.shape)))*imstd+immean)
        nim=im/im.mean(axis=0)
        m = np.array([ndimage.mean(nim[x],cellmask,1) for x in range(nim.shape[0])])
        #smooth it
        snim=signal.savgol_filter(m,35,2)
        #adjust baseline for variation in illumination intensity (temporal strobing effect)
        nim = nim/m[:,None,None]
        
        #adjust for change in baseling due to bleaching
        #nim = np.array([nim[x]*(1+(1-snim[x])) for x in range(nim.shape[0])])
        print 'smoothing'
        stdev=nim[l==1].std()*2    
        d=restoration.denoise_tv_chambolle(nim,stdev)-1.0
        
    
        
    
        
        
    
    
                              
    
    #includespk=[]
    
    #for ispk in range(df.index.max()):
    #    #check if spark is at the edge
    #    pos=np.int32(np.array((df['Peak Position'][ispk+1])[1:-1].split(),dtype=np.float32))
    #    if (cutspark(d,pos).shape[2]>8)and(cutspark(d,pos)[1]>8) :
    #        #if easygui.ynbox(msg="Is %.3g a spark"%np.float32(ispk),choices=("Yes","No"))==True:
    #        includespk.append(ispk)
    #    else:
    #        print 'spark too close to edge'
    includespk=ParseNumberList(easygui.enterbox(''))[0]
    interval=np.diff(np.array(gettimes(fn))).mean()*1000
    msdists=[]
    dists=[]
    areas=[]
    gmsdists=[]
    gdists=[]
    gareas=[]
    included=[]
    df2=pd.DataFrame(columns=['filename','spkno','time','x','y','z'])
    for i in range(len(includespk)):
        pos=np.int32(np.array((df['Peak Position'][includespk[i]])[1:-1].split(),dtype=np.float32))
        x,y,incl,gx,gy=fitspk2(im=d,pkpt=pos,interval=interval)
        pts=np.vstack((x,y)).transpose().tolist()
        gpts=np.vstack((gx,gy)).transpose().tolist()
        if len(pts)>3:
            included.append(includespk[i]+1)
            maxdist=spatial.distance.cdist(pts,pts).max()
            gmaxdist=spatial.distance.cdist(gpts,gpts).max()
            hull=ConvexHull(np.vstack((x,y)).transpose())
            ghull=ConvexHull(np.vstack((gx,gy)).transpose())
            areas.append(hull.area)
            gareas.append(ghull.area)
            msdist=msd(x,y)
            gmsdist=msd(gx,gy)
            print maxdist
            print gmaxdist
            dists.append(maxdist)
            gdists.append(gmaxdist)
            msdists.append(msdist)
            gmsdists.append(gmsdist)
            df.loc[i,'msdist']=msdist
            df.loc[i,'maxdist']=maxdist
            df.loc[i,'area']=hull.area
            df.loc[i,'gmsdist']=gmsdist
            df.loc[i,'gmaxdist']=gmaxdist
            df.loc[i,'garea']=ghull.area
            for p in range(len(pts)):
                df2.loc[len(df2)+1,'filename']=fn
                df2.loc[len(df2),'time']=incl[p]*interval
                df2.loc[len(df2),'spkno']=includespk[i]
                df2.loc[len(df2),'peakpos']=str(pos)
                df2.loc[len(df2),'x']=x[p]
                df2.loc[len(df2),'y']=y[p]
                df2.loc[len(df2),'gx']=gx[p]
                df2.loc[len(df2),'gy']=gy[p]
                df2.loc[len(df2),'z']=incl[p]
        else:
            'print spk %.3g is not long enough'%np.float32(includespk[i])
        
        
    
    if len(dists)>0:   
        dists=np.array(dists)
        #exclude crazy values
        dists[dists>2]=dists[dists<2].mean()
        df1=pd.DataFrame(columns=['filename','spkno','distance','msdist','areas','gdistance','gmsdist','gareas'])
        df1.areas=areas
        df1.spkno=included
        df1.msdist=msdists
        df1.distance=dists
        df1.areas=gareas
        df1.msdist=gmsdists
        df1.distance=gdists
        df1.to_csv((fn.lower()).replace('.tif','dists.csv'))      
        df2.to_csv((fn.lower()).replace('.tif','spkpos.csv'))
        df.to_csv(dffilename)
        #im=imread(fn)
        #d=normdenoise(im[0:300]-1)
        #pkpt=[188,7,30]
        #x,y=fitspk(d,pkpt)
        #plot(x,y)
    