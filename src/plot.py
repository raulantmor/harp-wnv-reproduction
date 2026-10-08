"""All numerical panels read saved independent calculations; no published pixels."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/harp-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import box
from settings import *

OUT=ROOT/'figures'
COLORS={'benchmark':'#aaaaaa','available_uni':'#9dc9a1','available_bi':'#246a42','published_uni':'#9dc9a1','published_bi':'#246a42'}
LABELS={'benchmark':'Historical benchmark','available_uni':'One available predictor',
        'available_bi':'Two available predictors','published_uni':'One Table 1 predictor',
        'published_bi':'Two Table 1 predictors'}
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})

def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'));fig.savefig(OUT/(name+'.png'),dpi=150);plt.close(fig)

def map_panel(ax,gdf,column,norm,cmap):
    gdf.plot(column=column,ax=ax,cmap=cmap,norm=norm,edgecolor='white',linewidth=.15,
             missing_kwds={'color':'#dedede'})
    ax.set_axis_off()

def maps():
    counties=gpd.read_file(PROCESSED/'counties.geojson').to_crs(5070)
    domains=gpd.read_file(PROCESSED/'regions.geojson')
    # Some continental NEON domains also contain Alaska or Puerto Rico pieces.
    # Clip display geometry only; assignments and numerical reductions are unchanged.
    domains=domains[domains.domainID.le('D17')].to_crs(4326).clip(box(-125,24,-66,50)).to_crs(5070)
    regional=pd.read_csv(PROCESSED/'case_summary.csv').rename(columns={'domainID':'domainID'})
    cc=pd.read_csv(PROCESSED/'county_cases.csv',dtype={'GEOID':str}).set_index('GEOID')
    stats=pd.DataFrame({'mean':cc.mean(axis=1),'median':cc.median(axis=1)})
    rc=domains.merge(regional,on='domainID');ct=counties.merge(stats,on='GEOID')
    fig,axs=plt.subplots(2,2,figsize=(11,6))
    for ax,g,col,label,norm in zip(axs.flat,[rc,rc,ct,ct],['mean','median','mean','median'],
                                  ['(a) Regional mean','(b) Regional median','(c) County mean','(d) County median'],
                                  [Normalize(0,regional['mean'].max()),Normalize(0,regional['mean'].max()),LogNorm(.05,stats['mean'].max()),LogNorm(.05,stats['mean'].max())]):
        g=g.copy();g.loc[g[col]==0,col]=np.nan
        map_panel(ax,g,col,norm,'viridis');ax.set_title(label)
        if 'included' in g:
            g[~g.included].plot(ax=ax,facecolor='none',hatch='///',edgecolor='#555555',linewidth=.2)
        fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),ax=ax,shrink=.6,label='WNND cases / year')
    fig.suptitle('Independent CDC input-data reconstruction, 2005–2022')
    save(fig,'figure1')
    corr=pd.read_csv(RESULTS/'correlations.csv')
    for number,var in [(2,'pdsi'),(3,'tmean')]:
        fig,axs=plt.subplots(2,2,figsize=(10,6))
        for ax,m,letter in zip(axs.flat,[10,1,4,7],'abcd'):
            c=corr[corr.predictor==f'{var}_{m:02}']
            g=domains.merge(c,left_on='domainID',right_on='region',how='left')
            map_panel(ax,g,'r',Normalize(-1,1),'RdBu_r')
            ax.set_title(f'({letter}) '+{10:'Previous October',1:'January',4:'April',7:'July'}[m])
        fig.subplots_adjust(right=.87)
        ca=fig.add_axes([.90,.22,.02,.55]);fig.colorbar(ScalarMappable(norm=Normalize(-1,1),cmap='RdBu_r'),cax=ca,label='Pearson correlation with ln(WNND)')
        fig.suptitle('Independent '+('gridMET PDSI' if var=='pdsi' else 'PRISM temperature')+' correlations')
        save(fig,f'figure{number}')
    # Separate schematic, not counted as reproduced model results.
    fig,ax=plt.subplots(figsize=(8,4))
    domains.boundary.plot(ax=ax,color='#555',linewidth=.5)
    for _,row in domains[domains.domainID.isin(REGIONS)].iterrows():
        p=row.geometry.representative_point();ax.text(p.x,p.y,row.domainID,ha='center',fontsize=8)
    ax.set_axis_off();ax.set_title('Study-region locator (schematic; official NEON boundaries)')
    save(fig,'region_schematic')

def correlations():
    c=pd.read_csv(RESULTS/'correlations.csv');sel=pd.read_csv(RESULTS/'selected_predictors.csv').set_index('region')
    names=['Mean temperature','Mean dew point','Precipitation','Maximum VPD','PDSI','Soil moisture']
    for n,r in enumerate(REGIONS,1):
        table=c[c.region==r].set_index('predictor').r
        values=np.array([[table.get(f'{v}_{m:02}',np.nan) for m in MONTHS] for v in VARIABLES])
        fig,ax=plt.subplots(figsize=(8,3.5));img=ax.imshow(values,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
        for i in range(6):
            for j in range(10):
                v=values[i,j]
                if np.isnan(v):ax.text(j,i,'NA',ha='center',va='center',fontsize=7,color='#555')
                elif abs(v)<.32:ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor='white',edgecolor='#ddd',linewidth=.3))
                else:ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=7,color='white' if abs(v)>.65 else 'black')
        for ps,color in zip([sel.loc[r,'primary'],sel.loc[r,'secondary']],['black','#777']):
            v,m=ps.split('_');ax.add_patch(Rectangle((MONTHS.index(int(m))-.5,VARIABLES.index(v)-.5),1,1,fill=False,edgecolor=color,linewidth=2))
        ax.set_xticks(range(10),['Oct−1','Nov−1','Dec−1','Jan','Feb','Mar','Apr','May','Jun','Jul'])
        ax.set_yticks(range(6),names);ax.set_title(f'S{n}: {REGIONS[r]} — available-input screening')
        fig.colorbar(img,ax=ax,label='Pearson r');fig.tight_layout();save(fig,f'figureS{n}')

def timeseries():
    forecasts=pd.read_csv(RESULTS/'forecasts.csv')
    mapping={'D06':'figure4',**{r:f'figureS{n}' for n,r in enumerate([r for r in REGIONS if r!='D06'],12)}}
    statuses=[]
    for r,name in mapping.items():
        f=forecasts[forecasts.region==r];published=all(m in f.model.values for m in ['published_uni','published_bi'])
        models=['benchmark','published_uni','published_bi'] if published else ['benchmark','available_uni','available_bi']
        fig,ax=plt.subplots(figsize=(11,3.6))
        for k,m in enumerate(models):
            rows=f[f.model==m].sort_values('year');boxes=[]
            for _,z in rows.iterrows():boxes.append(dict(med=z.q11,q1=z.q06,q3=z.q16,whislo=z.q01,whishi=z.q21,fliers=[]))
            artists=ax.bxp(boxes,positions=np.arange(18)+(k-1)*.24,widths=.2,patch_artist=True,showfliers=False,
                           boxprops={'facecolor':COLORS[m],'edgecolor':COLORS[m]},medianprops={'color':'black'},
                           whiskerprops={'color':COLORS[m]},capprops={'color':COLORS[m]})
            ax.plot([],[],color=COLORS[m],linewidth=7,label=LABELS[m])
        obs=f[f.model=='benchmark'].sort_values('year');ax.scatter(np.arange(18),obs.observed,color='black',s=15,zorder=5,label='CDC observed')
        ax.set_xticks(np.arange(18),YEARS,rotation=45);ax.set_ylabel('WNND cases');ax.set_ylim(bottom=0)
        ax.legend(ncol=4,fontsize=8,loc='lower center',bbox_to_anchor=(.5,1.0),frameon=False)
        ax.set_title(REGIONS[r]+(' — independent fits using Table 1 predictors' if published else ' — partial: screening excludes unavailable soil moisture'),pad=28)
        fig.tight_layout();save(fig,name)
        statuses.append(dict(figure=name,region=r,published_predictors=published))
    pd.DataFrame(statuses).to_csv(RESULTS/'time_series_panels.csv',index=False)

def skill():
    data=pd.read_csv(RESULTS/'skill.csv')
    for name,prefix,ylabel in [('figure5','','WIS improvement relative to benchmark (%)'),('figureS22','diff_','Benchmark WIS minus climate WIS')]:
        fig,ax=plt.subplots(figsize=(11,4.5))
        for k,m in enumerate(['available_uni','available_bi']):
            d=data[data.model==m].set_index('region').loc[ORDER+['National']]
            boxes=[dict(med=z[prefix+'median'],q1=z[prefix+'q25'],q3=z[prefix+'q75'],whislo=z[prefix+'lower'],whishi=z[prefix+'upper'],fliers=[]) for _,z in d.iterrows()]
            ax.bxp(boxes,positions=np.arange(12)+(k-.5)*.3,widths=.25,patch_artist=True,showfliers=False,boxprops={'facecolor':COLORS[m]},medianprops={'color':'black'})
            ax.plot([],[],color=COLORS[m],linewidth=7,label=LABELS[m])
        ax.axhline(0,color='black',linewidth=.6);ax.set_xticks(np.arange(12),[REGIONS.get(r,r) for r in ORDER+['National']],rotation=45,ha='right')
        ax.set_ylabel(ylabel);ax.set_title('Partial reconstruction: available-input models; 10,000 paired bootstrap samples')
        ax.legend();fig.tight_layout();save(fig,name)
    forecasts=pd.read_csv(RESULTS/'forecasts.csv')
    for n,m in [(23,'available_uni'),(24,'available_bi')]:
        f=forecasts[forecasts.model==m]
        hist,edges=np.histogram(f.pit,bins=np.linspace(0,1,11))
        pd.DataFrame(dict(left=edges[:-1],right=edges[1:],count=hist,proportion=hist/len(f))).to_csv(RESULTS/f'pit_{m}.csv',index=False)
        fig,ax=plt.subplots(figsize=(6,3));ax.bar((edges[:-1]+edges[1:])/2,hist/len(f),width=.095,color=COLORS[m])
        ax.axhline(.1,linestyle='--',color='black');ax.set(xlabel='Mid-distribution PIT',ylabel='Proportion of 198 forecasts',title=f'Partial S{n}: '+LABELS[m].lower());save(fig,f'figureS{n}')

def validation():
    data=np.load(RESULTS/'validation_draws.npz');chain=data['chain']
    fig,axs=plt.subplots(1,2,figsize=(10,3))
    for i in range(4):axs[0].plot(chain[i,::20,1],alpha=.6,linewidth=.6)
    axs[0].set(xlabel='Retained draw / 20',ylabel='Standardized slope 1',title='Independent Metropolis chains (synthetic data)')
    obs=data['observed'];pp=data['predictive'];bins=np.arange(0,max(np.quantile(pp,.995),obs.max())+5,5)
    axs[1].hist(obs,bins=bins,density=True,histtype='step',color='black',label='Synthetic observations')
    axs[1].hist(pp.ravel(),bins=bins,density=True,alpha=.5,color='#357a50',label='Posterior predictive')
    axs[1].set(xlabel='Count',ylabel='Density',title='Synthetic posterior predictive check');axs[1].legend(fontsize=8)
    fig.tight_layout();save(fig,'validation')

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    maps();correlations();timeseries();skill();validation()
