"""
Io for UFXC detector
Works off a buffer in memory.
"""
import numpy as np
import time
import os


"""    Description of UFXC files.

    They are a sequence of uint32 values, each values is a
    packed frame,cts,pixel integer. Top 11 bits are frame,
    bits 17,16 are the 2 bit counter, and last 15 bits are pixel.

"""

class UFXCfile:
    '''The class for raw UFXC data files.
       The file is in 1 based numbering scheme (first frame is 0)
       This reads whole file into memory and parcels it out as
       asked for (or one each read).
    '''
    def __init__(self,filename,beg,end):
        '''UFXC initialization. Read whole file into memory.
            NOTE: At each record n, the file cursor points to record n+1,
            so ready to read next image.
        '''
        self.rows=128
        self.cols=256
        self.filename = filename
        self.recno = 1
        self.ifDK=False
        self.ifTHRESHOLD=False
        self.THRESHOLD=0.0
        self.ifNORM=False
        # some initialization stuff
        self.byts=2
        data=np.fromfile(filename,count=-1,dtype=np.uint32)
        no=data.size
        fr=np.int32(data>>21)
        for f in fr[:4096]:
            print(f)
        #handle overflows in frame number
        # pos=1+np.where((fr[1:]-fr[:-1])< 0)[0]
        # for p in pos:
        #     fr[p:] += 2048 #2**11
        # fr -= fr[0] #first frame number is random
        # self.cts=np.int32((data>>15)&0x3)
        # self.pix=np.int32(data&0x7fff)
        # del data #not used any more
        # self.fr_pos=np.r_[0,1+np.where(fr[1:]!=fr[:-1])[0],fr.size]
        # self.fr=fr
        # self.beg = beg
        # if(end==-1): end=fr[-1]+1;
        # self.number = end-beg+1

    def rdimg(self,n):
        '''raw read image n. Use rdframe(n) for processed image.
        '''
        p=self.pix[self.fr_pos[n-1]:self.fr_pos[n]]
        v=self.cts[self.fr_pos[n-1]:self.fr_pos[n]]
        img = np.zeros((self.rows*self.cols))
        np.put(img,p, v)
        img.shape = (self.cols,self.rows)
        return(img)
    def rdframe(self,n,tdk=True,tflag=True):
        ''' read processed images. Used by all analysis routines.
        '''
        img=self.rdimg(n).astype(dtype=np.float64) #newer numpy
        #img=self.rdimg(n).astype('float64') # older numpy
        if(self.ifDK):
                if(tdk): img -= self.DK
        if(self.ifTHRESHOLD):
                if(tflag): img *= (img>self.THRESHOLD)
        if(self.ifNORM): img *= self.NORM
        return(img)
    def rdrawframe(self,n):
        #note darks and thresholding down by compressing
        p=self.pix[self.fr_pos[n-1]:self.fr_pos[n]]
        v=np.float64(self.cts[self.fr_pos[n-1]:self.fr_pos[n]])
        if self.ifNORM: v *= self.NORM.ravel()[p]
        return(p,v)

#ldimgs
def ldimgs(FD,begframe=None,noframes=None):

    if(begframe==None):begframe=FD.beg
    if(noframes==None):noframes=FD.end-FD.beg+1
    ta = time.time()
    imgs=np.zeros([noframes, FD.row_end-FD.row_beg,FD.col_end-FD.col_beg])
    p=0
    for i in range(begframe,begframe+noframes):
        imgs[p,:,:]=FD.rdframe(i)
        p+=1
    tb = time.time()
    print("total time: {0}".format(tb-ta))
    return(imgs)


if __name__ == '__main__':
    import sys
    ufx = UFXCfile(sys.argv[1], 0, 100)