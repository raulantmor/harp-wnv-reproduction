"""Log1p weighted interval score; exact median/dispersion/over/under accounting."""
import numpy as np
from settings import QUANTILES

def wis(observed,quantiles,transform=True):
    q=np.asarray(quantiles,float);y=np.asarray(observed,float)
    if q.shape[-1]!=23 or np.any(np.diff(q,axis=-1)<0): raise ValueError('Expected 23 ordered quantiles')
    if transform:q=np.log1p(q);y=np.log1p(y)
    lower=q[...,:11];upper=q[...,-1:11:-1];median=q[...,11]
    alpha=2*QUANTILES[:11];norm=11.5
    dispersion=np.sum(alpha/2*(upper-lower),axis=-1)/norm
    over=(np.maximum(lower-y[...,None],0).sum(axis=-1)+.5*np.maximum(median-y,0))/norm
    under=(np.maximum(y[...,None]-upper,0).sum(axis=-1)+.5*np.maximum(y-median,0))/norm
    return dict(wis=dispersion+over+under,dispersion=dispersion,over=over,under=under)

def bootstrap(scores,totals,seed=0,n=10000):
    """scores: region x model x year; paired years; national region+year bootstrap."""
    rng=np.random.default_rng(seed);nr,nm,ny=scores.shape
    years=rng.integers(ny,size=(n,nr,ny))
    samples=np.empty((n,nr,nm))
    for r in range(nr):samples[:,r,:]=scores[r][:,years[:,r,:]].mean(axis=2).T
    ix=rng.integers(nr,size=(n,nr))
    selected=samples[np.arange(n)[:,None],ix]
    w=np.asarray(totals)[ix]
    national=(selected*w[:,:,None]).sum(axis=1)/w.sum(axis=1)[:,None]
    return np.concatenate([samples,national[:,None,:]],axis=1)
