"""Full-data posterior predictive checks and observed-data diagnostic summaries."""
import json
import numpy as np
import pandas as pd
from settings import *
from forecast import inputs,select
from model import fit,predict,design
from validate import metropolis

def run():
    x,y=inputs();rows=[];diagnostics=[]
    for k,r in enumerate(ORDER):
        xx=x.loc[r].loc[YEARS];yy=y.loc[r,YEARS].to_numpy(float)
        selected,_=select(xx,yy)
        for j,ps in enumerate([[],selected[:1],selected]):
            a=xx[ps].to_numpy() if ps else np.empty((18,0))
            posterior=fit(a,yy,SEED+1000+k*10+j)
            pred=predict(posterior,a,SEED+2000+k*10+j,size=4000)
            for name,obs,rep in [('mean',yy.mean(),pred.mean(axis=1)),('variance',yy.var(),pred.var(axis=1)),('max',yy.max(),pred.max(axis=1))]:
                q=np.quantile(rep,[.025,.5,.975])
                rows.append(dict(region=r,model=['benchmark','available_uni','available_bi'][j],statistic=name,observed=obs,lower=q[0],median=q[1],upper=q[2],bayes_p=float(np.mean(rep>=obs))))
            diagnostics.append(dict(region=r,model=j,**posterior['diagnostics']))
    pd.DataFrame(rows).to_csv(RESULTS/'posterior_predictive_checks.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(RESULTS/'full_fit_diagnostics.csv',index=False)
    f=pd.read_csv(RESULTS/'forecasts.csv');d=pd.read_csv(RESULTS/'diagnostics.csv')
    # Challenge the numerical integrator on an observed heavy-weight-tail fold
    # and the influential Desert Southwest 2021 outbreak holdout.
    checks=[];chains={}
    stress=d[(d.region=='D14')&(d.year==2021)&(d.model=='published_bi')].iloc[0]
    for k,z in enumerate([d.loc[d.pareto_k.idxmax()],stress]):
        row=f[(f.region==z.region)&(f.year==z.year)&(f.model==z.model)].iloc[0]
        ps=[] if pd.isna(row.predictors) else row.predictors.split(',')
        xx=x.loc[z.region].loc[YEARS]
        a=xx[ps].to_numpy() if ps else np.empty((18,0))
        train=np.array(YEARS)!=z.year;yy=y.loc[z.region,YEARS].to_numpy(float)[train]
        posterior=fit(a[train],yy,SEED+3000+k)
        X,_,_=design(a[train])
        chain,check=metropolis(X,yy,posterior['mode'],posterior['cov'],SEED+4000+k)
        mean=posterior['weights']@posterior['theta']
        sd=np.sqrt(posterior['weights']@((posterior['theta']-mean)**2))
        delta=float(np.max(abs(chain.mean(axis=(0,1))-mean)/sd))
        checks.append(dict(region=z.region,year=int(z.year),model=z.model,
                           rhat_max=check['rhat_max'],ess_bulk_min=check['ess_bulk_min'],
                           max_mean_difference_in_posterior_sd=delta,
                           numerical_pass=bool(check['rhat_max']<=1.01 and check['ess_bulk_min']>=400 and delta<.15)))
        chains[f'{z.region}_{z.model}_{z.year}']=chain
    pd.DataFrame(checks).to_csv(RESULTS/'observed_mcmc_checks.csv',index=False)
    np.savez_compressed(RESULTS/'observed_mcmc_draws.npz',**chains)
    f.groupby(['region','model'])[['coverage50','coverage90','wis','dispersion','over','under']].mean().to_csv(RESULTS/'calibration.csv')
    biggest=f.loc[f.groupby(['region','model']).median_error.apply(lambda z:z.abs().idxmax())]
    biggest.to_csv(RESULTS/'largest_absolute_errors.csv',index=False)
    summary=dict(fits=len(d),failed_importance_checks=int((~d.numerical_pass).sum()),
                 min_importance_ess=float(d.importance_ess.min()),max_pareto_k=float(d.pareto_k.max()),
                 max_mcse_over_sd=float(d.max_mcse_over_sd.max()),map_optimizer_warnings=int((~d.map_optimizer_success).sum()),
                 observed_mcmc_checks_pass=all(z['numerical_pass'] for z in checks),
                 full_fit_checks_pass=all(z['numerical_pass'] for z in diagnostics))
    (RESULTS/'observed_validation.json').write_text(json.dumps(summary,indent=2))
    print(summary)
    if summary['failed_importance_checks'] or not summary['observed_mcmc_checks_pass'] or not summary['full_fit_checks_pass']:
        raise RuntimeError('Some fitted posteriors failed diagnostics')

if __name__=='__main__':run()
