"""Parameter recovery, independent MCMC check, and analytic benchmark checks."""
import json
import numpy as np
import pandas as pd
import arviz as az
from scipy.stats import nbinom
from settings import RESULTS, SEED, QUANTILES
from model import fit, predict, log_posterior, design, parameter_summary
from score import wis

def metropolis(X,y,mode,cov,seed):
    rng=np.random.default_rng(seed);d=len(mode)
    theta=rng.multivariate_normal(mode,2*cov,size=4);lp=log_posterior(theta,X,y)
    L=np.linalg.cholesky(cov)*(2.38/np.sqrt(d))
    samples=[];accept=np.zeros(4)
    for i in range(14000):
        cand=theta+rng.normal(size=(4,d))@L.T
        new=log_posterior(cand,X,y);a=np.log(rng.uniform(size=4))<new-lp
        theta[a]=cand[a];lp[a]=new[a]
        if i>=4000:samples.append(theta.copy());accept+=a
    chain=np.transpose(samples,(1,0,2))
    data=az.from_dict(posterior={f'p{i}':chain[:,:,i] for i in range(d)})
    stats=az.summary(data)
    return chain,dict(rhat_max=float(az.rhat(data).to_array().values.max()),ess_bulk_min=float(az.ess(data,method='bulk').to_array().values.min()),acceptance_rate=(accept/10000).tolist())

def run():
    RESULTS.mkdir(exist_ok=True)
    rng=np.random.default_rng(SEED);x=rng.normal(size=(240,2))
    truth=np.array([3.,.55,-.4,4.])
    mu=np.exp(truth[0]+x@truth[1:3]);y=rng.negative_binomial(truth[-1],truth[-1]/(truth[-1]+mu))
    posterior=fit(x,y,SEED);summ=pd.DataFrame(parameter_summary(posterior));summ['truth']=truth
    summ['covered95']=(summ.lower<=truth)&(truth<=summ.upper)
    summ.to_csv(RESULTS/'synthetic_recovery.csv',index=False)
    X,_,_=design(x)
    chain,diag=metropolis(X,y,posterior['mode'],posterior['cov'],SEED+1)
    imean=posterior['weights']@posterior['theta']
    isd=np.sqrt(posterior['weights']@((posterior['theta']-imean)**2))
    diag['max_mean_difference_in_posterior_sd']=float(np.max(abs(chain.mean(axis=(0,1))-imean)/isd))
    diag['pass']=bool(diag['rhat_max']<=1.01 and diag['ess_bulk_min']>=400 and diag['max_mean_difference_in_posterior_sd']<.15)
    pred=predict(posterior,x,SEED+2,size=4000)
    ppc=dict(observed_mean=float(y.mean()),observed_variance=float(y.var()),observed_max=int(y.max()),
             replicate_mean_interval=np.quantile(pred.mean(axis=1),[.025,.5,.975]).tolist(),
             replicate_variance_interval=np.quantile(pred.var(axis=1),[.025,.5,.975]).tolist(),
             replicate_max_interval=np.quantile(pred.max(axis=1),[.025,.5,.975]).tolist())
    np.savez_compressed(RESULTS/'validation_draws.npz',chain=chain,observed=y,predictive=pred[:500])
    # Repeated simulations measure interval recovery rather than a single lucky fit.
    cover=[]
    for k in range(24):
        yy=rng.negative_binomial(truth[-1],truth[-1]/(truth[-1]+mu))
        s=pd.DataFrame(parameter_summary(fit(x,yy,SEED+100+k)))
        cover.append((s.lower.to_numpy()<=truth)&(truth<=s.upper.to_numpy()))
    # Intercept-only model is invariant to training-observation ordering.
    base=fit(np.empty((len(y),0)),y,SEED)
    rev=fit(np.empty((len(y),0)),y[::-1],SEED)
    delta=float(np.max(abs(base['mode']-rev['mode'])))
    identical=wis(5,np.full(23,5))['wis']==0
    # Independent quantile-loss identity for WIS = mean of 23 pinball losses * 2.
    q=nbinom.ppf(QUANTILES,3,.2);obs=17
    error=np.log1p(obs)-np.log1p(q)
    independent=2*np.mean(np.maximum(QUANTILES*error,(QUANTILES-1)*error))
    wis_delta=float(abs(wis(obs,q)['wis']-independent))
    result=dict(synthetic_importance=posterior['diagnostics'],mcmc=diag,posterior_predictive=ppc,
                recovery_24_simulations=np.mean(cover,axis=0).tolist(),benchmark_permutation_map_delta=delta,
                wis_identity_error=wis_delta,perfect_forecast_zero=bool(identical))
    result['pass']=bool(diag['pass'] and posterior['diagnostics']['numerical_pass'] and delta<1e-3 and wis_delta<1e-12 and identical)
    (RESULTS/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    if not result['pass']:raise RuntimeError('Numerical validation failed')

if __name__=='__main__':run()
