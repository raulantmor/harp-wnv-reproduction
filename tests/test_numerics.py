import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.stats import nbinom
from model import log_posterior,log_gradient,design
from score import wis,bootstrap
from forecast import select
from settings import QUANTILES
import pandas as pd

def test_likelihood_matches_independent_scipy_parameterization():
    X=np.array([[1,-1],[1,0],[1,1.]])
    y=np.array([0,4,18]);t=np.array([2.,.4,np.log(3.)])
    mu=np.exp(X@t[:2]);phi=3.
    expected=nbinom.logpmf(y,phi,phi/(phi+mu)).sum()-.5*np.sum((t[:2]/2.5)**2)-phi+np.log(phi)
    assert np.isclose(log_posterior(t,X,y)[0],expected,atol=1e-12)

def test_analytic_gradient_agrees_with_central_differences():
    X=np.array([[1.,-1.2,.4],[1.,.3,-.8],[1.,2.,.1]])
    y=np.array([0.,7.,30.]);t=np.array([2.,.4,-.2,.6]);eps=1e-5
    eye=np.eye(4)*eps
    numerical=np.array([(log_posterior(t+e,X,y)[0]-log_posterior(t-e,X,y)[0])/(2*eps) for e in eye])
    np.testing.assert_allclose(log_gradient(t,X,y),numerical,atol=1e-7,rtol=1e-6)

def test_heldout_covariate_does_not_change_training_scaling():
    _,c,s=design(np.array([[1.],[2.],[3.]]))
    held,_,_=design(np.array([[100.]]),c,s)
    assert held[0,1]==98

def test_wis_quantile_loss_identity_and_components():
    rng=np.random.default_rng(12)
    for _ in range(10):
        q=np.sort(rng.integers(0,100,23));y=int(rng.integers(0,100));e=np.log1p(y)-np.log1p(q)
        independent=2*np.mean(np.maximum(QUANTILES*e,(QUANTILES-1)*e))
        s=wis(y,q)
        assert np.isclose(s['wis'],independent,atol=1e-12)
        assert np.isclose(s['wis'],s['dispersion']+s['over']+s['under'])
        assert all(v>=0 for v in s.values())

def test_paired_bootstrap_identical_models_exactly_tie():
    rng=np.random.default_rng(33);a=rng.uniform(size=(11,1,18));a=np.repeat(a,3,axis=1)
    b=bootstrap(a,np.arange(1,12),n=200)
    np.testing.assert_array_equal(b[:,:,0],b[:,:,1])

def test_secondary_predictor_from_different_group():
    x=pd.DataFrame({'tmean_01':np.arange(18), 'tdmean_01':np.arange(18)*2,'ppt_03':np.arange(18)+np.sin(np.arange(18))})
    pair,_=select(x,np.exp(np.arange(18)/10))
    assert pair[1]=='ppt_03'
