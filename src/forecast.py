"""Fit all leave-one-year-out models, save quantiles, diagnostics and score samples."""
import argparse
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from settings import *
from model import fit,predict,parameter_summary
from score import wis,bootstrap

def inputs():
    c=pd.read_csv(PROCESSED/'regional_climate.csv')
    c['predictor']=c.variable+'_'+c.month.astype(str).str.zfill(2)
    x=c.pivot(index=['region','year'],columns='predictor',values='value')
    needed=[f'{v}_{m:02}' for v in ['tmean','tdmean','ppt','vpdmax','pdsi'] for m in MONTHS]
    if len(x)!=198 or x.reindex(columns=needed).isna().any().any():
        raise ValueError('Incomplete required PRISM/gridMET input matrix; finish downloads and processing')
    y=pd.read_csv(PROCESSED/'regional_cases.csv',index_col=0)
    y.columns=y.columns.astype(int)
    return x,y

def select(x,y):
    # Full-period screening is intentional: follows published retrospective design.
    if np.any(y<=0):raise ValueError('Exploratory log counts require positive regional annual counts')
    r={k:pearsonr(x[k],np.log(y)).statistic for k in x if x[k].notna().all()}
    ranked=sorted(r,key=lambda k:(-abs(r[k]),k))
    first=ranked[0];second=next(k for k in ranked if GROUP[k.split('_')[0]]!=GROUP[first.split('_')[0]])
    return [first,second],r

def run():
    RESULTS.mkdir(exist_ok=True);(RESULTS/'folds').mkdir(exist_ok=True)
    x,y=inputs(); corr=[];selection=[];forecasts=[];diagnostics=[];params=[];draws={}
    for region in ORDER:
        xr=x.loc[region].loc[YEARS];yr=y.loc[region,YEARS].to_numpy(float)
        chosen,r=select(xr,yr)
        for k,v in r.items():corr.append(dict(region=region,predictor=k,r=v))
        selection.append(dict(region=region,primary=chosen[0],secondary=chosen[1],
                              paper_primary=PUBLISHED[region][0],paper_secondary=PUBLISHED[region][1],
                              primary_match=chosen[0]==PUBLISHED[region][0],secondary_match=chosen[1]==PUBLISHED[region][1]))
        specs={'benchmark':[], 'available_uni':chosen[:1], 'available_bi':chosen}
        for n in [1,2]:
            ps=PUBLISHED[region][:n]
            if all(k in xr and xr[k].notna().all() for k in ps):specs['published_'+('uni' if n==1 else 'bi')]=ps
        for mi,(model,ps) in enumerate(specs.items()):
            xx=xr[ps].to_numpy() if ps else np.empty((18,0))
            for j,year in enumerate(YEARS):
                key=f'{region}_{model}_{year}'; path=RESULTS/'folds'/f'{key}.npz'
                seed=SEED+10000*ORDER.index(region)+100*mi+j
                payload=xx.tobytes()+yr.tobytes()+str(seed).encode()+','.join(ps).encode()
                for source in ['model.py','settings.py','forecast.py']:
                    payload+=(ROOT/'src'/source).read_bytes()
                signature=hashlib.sha256(payload).hexdigest()
                saved=np.load(path) if path.exists() else None
                if saved is not None and 'signature' in saved and str(saved['signature'])==signature:
                    pred=saved['pred'];diag=json.loads(str(saved['diagnostics']));summ=json.loads(str(saved['parameters']))
                else:
                    train=np.arange(18)!=j
                    posterior=fit(xx[train],yr[train],seed)
                    pred=predict(posterior,xx[j:j+1],seed+100000)[:,0]
                    diag=posterior['diagnostics'];summ=parameter_summary(posterior)
                    np.savez_compressed(path,pred=pred,diagnostics=json.dumps(diag),parameters=json.dumps(summ),signature=signature)
                q=np.quantile(pred,QUANTILES,method='inverted_cdf')
                base=dict(region=region,year=int(year),model=model)
                row=dict(**base,observed=int(yr[j]),predictors=','.join(ps),**{f'q{i:02}':float(v) for i,v in enumerate(q)})
                row.update({k:float(v) for k,v in wis(yr[j],q).items()})
                row['pit']=(np.mean(pred<yr[j])+np.mean(pred<=yr[j]))/2
                row['median_error']=q[11]-yr[j]
                row['coverage50']=bool(q[6]<=yr[j]<=q[16]);row['coverage90']=bool(q[2]<=yr[j]<=q[20])
                forecasts.append(row);diagnostics.append(dict(**base,**diag))
                params.extend([dict(**base,**p) for p in summ]);draws[key]=pred
            print(region,model,'complete',flush=True)
    pd.DataFrame(corr).to_csv(RESULTS/'correlations.csv',index=False)
    pd.DataFrame(selection).to_csv(RESULTS/'selected_predictors.csv',index=False)
    f=pd.DataFrame(forecasts);f.to_csv(RESULTS/'forecasts.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(RESULTS/'diagnostics.csv',index=False)
    pd.DataFrame(params).to_csv(RESULTS/'parameters.csv',index=False)
    np.savez_compressed(RESULTS/'predictive_draws.npz',**draws)
    # Use a common complete family for national aggregation; never mix regional models.
    models=['benchmark','available_uni','available_bi']
    scores=np.array([[f[(f.region==r)&(f.model==m)].sort_values('year').wis.to_numpy() for m in models] for r in ORDER])
    totals=y.loc[ORDER,YEARS].sum(axis=1).to_numpy()
    b=bootstrap(scores,totals,SEED)
    np.savez_compressed(RESULTS/'bootstrap.npz',scores=b,regions=ORDER+['National'],models=models)
    rows=[]
    for ri,r in enumerate(ORDER+['National']):
        for mi,m in enumerate(models[1:],1):
            diff=b[:,ri,0]-b[:,ri,mi];pct=100*diff/b[:,ri,0]
            q=np.quantile(pct,[.025,.25,.5,.75,.975]);dq=np.quantile(diff,[.025,.25,.5,.75,.975])
            rows.append(dict(region=r,model=m,lower=q[0],q25=q[1],median=q[2],q75=q[3],upper=q[4],
                             diff_lower=dq[0],diff_q25=dq[1],diff_median=dq[2],diff_q75=dq[3],diff_upper=dq[4],p=float(np.mean(diff<=0))))
    pd.DataFrame(rows).to_csv(RESULTS/'skill.csv',index=False)
    f.sort_values('wis',ascending=False).groupby(['region','model']).head(3).to_csv(RESULTS/'largest_errors.csv',index=False)

if __name__=='__main__':run()
