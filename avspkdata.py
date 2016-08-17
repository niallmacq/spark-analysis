# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 12:32:19 2016

@author: MVLSTECH
"""

import pandas as pd
import numpy as np
import easygui
from PIL import Image
import re
from os import path

def getsec(im):
    # Grab time stamp from TIFF tag 306
    stamp = im.ifd.tagdata.get(306)
    # Search for everything between single quotes
    # Tag is a tuple: (2, '20100713 16:49:20.690\x00')
    #stamp = re.search('\'.*\'', str(stamp))
    # Strip single quotes
    #stamp = stamp.group(0).strip("'")
    # Break off end with null
    #stamp = stamp.split("\\")
    # Stringify first element of list
    s=stamp
    date, time = s.split(" ")
    hour, mins, s = time.split(":")
    s,ms=s.split('.')
    seconds = (360 * np.int(hour) + 60 * np.int(mins) +  np.int(s)+float(re.sub(u'\x00', '', ms))/1000) 
    return seconds
    
def gettimes(fn):
    im=Image.open(fn)
    ss=[]
    try:
        seconds=getsec(im)
        starttime=seconds
        im.seek(im.tell()+1) # seek to next image
        while 1:
            seconds=getsec(im)-starttime
            # Get found 
            #print stamp.group(0)
            ss.append(seconds) # blarg out seconds
            im.seek(im.tell()+1) # seek to next image
    except EOFError:
    	pass
    return ss
    
    
if __name__ == '__main__': 
    
    pixsize=0.245    
    
    fns=easygui.fileopenbox(multiple=True)
    for fn in fns:
        dffilename=(fn.lower()).replace('.tif','results.csv')
        df=pd.read_csv(dffilename)
        df1=pd.DataFrame(columns=['Filename','TotalTime','CellArea','SparkFreq','TTP','FDHM','FWHM'])
        
        
        nspks=len(df)
        
        
        cellmaskname=(fn.lower()).replace('.tif','cellmask.tif.npy')
        
        if path.exists(cellmaskname):
            cellmask=np.load(cellmaskname)
        else:
            print 'No mask file found for: '+fn
            
        #get interval    
        ss=np.float32(gettimes(fn))
        nframes=len(ss)+1
        interval=np.diff(ss).mean()*1000
        df1.loc[len(df1)+1,'Filename']=fn
        df1.loc[len(df1),'TotalTime']=interval*nframes
        df1.loc[len(df1),'CellArea']=cellmask.sum()*pixsize*pixsize
        df1.loc[len(df1),'SparkFreq']=(nspks/(interval*nframes/1000)/(cellmask.sum()*pixsize*pixsize))
        df1.loc[len(df1),'TTP']=df[df.TTP>interval].TTP.mean()
        df1.loc[len(df1),'FDHM']=df[df.fdhm>0].TTP.mean()
        df1.loc[len(df1),'FWHM']=df[(df.fwhm>1)&(df.fwhm<5)].fwhm.mean()
        
        
    print df1
    
    