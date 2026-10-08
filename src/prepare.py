"""County assignment, WNND sums, and independent raster-to-county reductions."""
import argparse
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from rasterio.features import rasterize
from rasterio.transform import from_origin
from settings import RAW, PROCESSED, REGIONS, YEARS, MONTHS

def geography():
    PROCESSED.mkdir(parents=True,exist_ok=True)
    counties=gpd.read_file('zip://'+str(RAW/'counties.zip')).to_crs(4326)
    counties=counties[~counties.STATEFP.isin(['02','15','60','66','69','72','78'])].copy()
    centers=pd.read_csv(RAW/'county_population_centers.txt',dtype={'STATEFP':str,'COUNTYFP':str})
    centers['GEOID']=centers.STATEFP+centers.COUNTYFP
    counties=counties.merge(centers[['GEOID','LATITUDE','LONGITUDE']],on='GEOID',validate='one_to_one')
    domains=gpd.read_file('zip://'+str(RAW/'neon.zip')+'!NEONDomains_2024/NEON_Domains.shp').to_crs(4326)
    points=gpd.GeoDataFrame(counties[['GEOID']].copy(),geometry=gpd.points_from_xy(counties.LONGITUDE,counties.LATITUDE),crs=4326)
    assignment=gpd.sjoin(points,domains[['domainID','geometry']],how='left',predicate='intersects')
    assert not assignment.GEOID.duplicated().any()
    missing=assignment.domainID.isna()
    if missing.any():
        nearest=gpd.sjoin_nearest(points.loc[missing].to_crs(5070),domains[['domainID','geometry']].to_crs(5070),distance_col='distance_m')
        nearest[['GEOID','domainID','distance_m']].to_csv(PROCESSED/'nearest_region_assignments.csv',index=False)
        assignment.loc[missing,'domainID']=nearest.domainID.values
    counties=counties.merge(assignment[['GEOID','domainID']],on='GEOID',validate='one_to_one').sort_values('GEOID').reset_index(drop=True)
    counties[['GEOID','NAME','domainID','LATITUDE','LONGITUDE']].to_csv(PROCESSED/'county_regions.csv',index=False)
    counties.geometry=counties.geometry.simplify(.002,preserve_topology=True)
    # Simplification is ONLY for display. Original polygons are read for raster means.
    counties[['GEOID','domainID','geometry']].to_file(PROCESSED/'counties.geojson',driver='GeoJSON')
    domains[['domainID','domainName','geometry']].to_file(PROCESSED/'regions.geojson',driver='GeoJSON')
    d=pd.read_csv(RAW/'cdc_county.csv',dtype={'County':str})
    ct=pd.read_csv(RAW/'cdc_ct_0.csv',dtype={'County':str})
    d=d[~((d.Year==2022)&d.County.str.startswith('09'))]
    d=pd.concat([d,ct],ignore_index=True)
    d=d[d.Year.isin(YEARS)].copy()
    d['County']=d.County.replace({'46113':'46102','51515':'51019'})
    d[(~d.County.isin(counties.GEOID))&(d['Neuroinvasive disease cases']>0)].to_csv(PROCESSED/'excluded_reported_cases.csv',index=False)
    d=d[d.County.isin(counties.GEOID)]
    assert d['Neuroinvasive disease cases'].notna().all(), 'Unresolved missing WNND counts'
    pivot=d.groupby(['County','Year'])['Neuroinvasive disease cases'].sum().unstack(fill_value=0)
    pivot=pivot.reindex(index=counties.GEOID,columns=YEARS,fill_value=0).fillna(0).astype(int)
    pivot.index.name='GEOID'; pivot.to_csv(PROCESSED/'county_cases.csv')
    region=pivot.join(counties.set_index('GEOID').domainID).groupby('domainID').sum()
    region.to_csv(PROCESSED/'regional_cases.csv')
    state=pd.read_csv(RAW/'cdc_state.csv')
    state=state[(state.Type=='Neuroinvasive disease cases')&state.Year.isin([str(v) for v in YEARS])&~state.State.isin(['AK','HI','PR'])]
    county_sum=pivot.sum();county_sum.index=county_sum.index.astype(str)
    state_sum=state.groupby('Year')['Reported Cases'].sum()
    audit=pd.DataFrame({'county_sum':county_sum,'state_sum':state_sum})
    audit['difference']=audit.county_sum-audit.state_sum
    audit.to_csv(PROCESSED/'case_state_audit.csv')
    summary=region.T.agg(['sum','mean','median']).T
    summary['included']=summary.index.isin(REGIONS)
    summary.to_csv(PROCESSED/'case_summary.csv')
    print(summary.to_string())
    threshold=set(summary.index[summary['sum']>=1000])
    (PROCESSED/'region_threshold_audit.json').write_text(json.dumps(dict(
        rule='at least 1000 cases in 2005-2022', qualifying=sorted(threshold),
        published=sorted(REGIONS),difference=sorted(threshold.symmetric_difference(REGIONS)),
        decision='Use the 11 explicitly named published regions for comparability'),indent=2))

def grid_operator(lats,lons,tag):
    """Area-weighted cell-center rasterization. Small counties use nearest land cell."""
    counties=gpd.read_file('zip://'+str(RAW/'counties.zip')).to_crs(4326)
    cross=pd.read_csv(PROCESSED/'county_regions.csv',dtype={'GEOID':str})
    counties=counties.set_index('GEOID').loc[cross.GEOID]
    assert np.all(np.diff(lats)<0) and np.all(np.diff(lons)>0)
    dx=np.diff(lons).mean();dy=-np.diff(lats).mean()
    transform=from_origin(lons[0]-dx/2,lats[0]+dy/2,dx,dy)
    labels=rasterize([(g,i+1) for i,g in enumerate(counties.geometry)],out_shape=(len(lats),len(lons)),transform=transform,dtype='int32')
    areas=np.broadcast_to(np.cos(np.deg2rad(lats))[:,None],labels.shape).ravel()
    ci=np.abs(cross.LONGITUDE.to_numpy()[:,None]-lons).argmin(axis=1)
    ri=np.abs(cross.LATITUDE.to_numpy()[:,None]-lats).argmin(axis=1)
    return labels.ravel(),areas,ri,ci,len(cross)

def reduce_grid(values,op):
    labels,areas,ri,ci,n=op
    values=np.asarray(values,dtype=float)
    a=values.ravel(); valid=np.isfinite(a)&(labels>0)
    sums=np.bincount(labels[valid],weights=(a*areas)[valid],minlength=n+1)[1:]
    weights=np.bincount(labels[valid],weights=areas[valid],minlength=n+1)[1:]
    out=np.divide(sums,weights,out=np.full(n,np.nan),where=weights>0)
    for k in np.where(~np.isfinite(out))[0]:
        # Find nearest valid grid cell to population center; needed on coastlines.
        rr,cc=ri[k],ci[k]
        for radius in [0,1,2,4,8,16]:
            r0,r1=max(0,rr-radius),min(values.shape[0],rr+radius+1)
            c0,c1=max(0,cc-radius),min(values.shape[1],cc+radius+1)
            r,c=np.where(np.isfinite(values[r0:r1,c0:c1]))
            if len(r):
                j=np.argmin((r+r0-rr)**2+(c+c0-cc)**2)
                out[k]=values[r[j]+r0,c[j]+c0]; break
    assert np.isfinite(out).all(), 'County climate coverage incomplete'
    return out, int((weights==0).sum())

def climate():
    cross=pd.read_csv(PROCESSED/'county_regions.csv',dtype={'GEOID':str})
    cases=pd.read_csv(PROCESSED/'county_cases.csv',dtype={'GEOID':str}).set_index('GEOID').loc[cross.GEOID]
    totals=cases.sum(axis=1).to_numpy()
    weights=np.log(np.maximum(totals,1))  # zero total excluded; total=1 has log weight zero
    rows=[];county_rows=[]; audit=[];ops={}
    files=sorted((RAW/'prism').glob('*.nc'))+sorted((RAW/'pdsi').glob('pdsi_????.nc'))+sorted((RAW/'nldas').glob('*.nc'))
    for k,f in enumerate(files):
        with xr.open_dataset(f,decode_times=False) as d:
            if f.parent.name=='prism':
                var,date=f.stem.split('_'); y,m=int(date[:4]),int(date[4:])
                lats=d.latitude.values;lons=d.longitude.values
                fields=[(y,m,d[var].values)];tag='prism'
            elif f.parent.name=='pdsi':
                var='pdsi';tag='pdsi';lats=d.lat.values;lons=d.lon.values
                import netCDF4
                times=netCDF4.num2date(d.day.values,d.day.units)
                a=d['daily_mean_palmer_drought_severity_index'].values
                fields=[]
                for m in MONTHS:
                    sel=np.array([t.month==m for t in times]); y=times[0].year
                    if sel.any():
                        den=np.isfinite(a[sel]).sum(axis=0)
                        mean=np.divide(np.nansum(a[sel],axis=0),den,out=np.full(a.shape[1:],np.nan),where=den>0)
                        fields.append((y,m,mean))
            else:
                var='soil';tag='nldas';lats=d.lat.values;lons=d.lon.values
                import netCDF4
                times=netCDF4.num2date(d.time.values,d.time.units)
                a=d['SoilM_0_10cm'].values
                fields=[(t.year,t.month,a[i]) for i,t in enumerate(times)]
            if lats[0]<lats[-1]:
                lats=lats[::-1];fields=[(y,m,a[::-1]) for y,m,a in fields]
            if tag not in ops: ops[tag]=grid_operator(lats,lons,tag)
            for y,m,a in fields:
                target=y+1 if m>=10 else y
                if target not in YEARS or m not in MONTHS: continue
                vals,nfallback=reduce_grid(a,ops[tag])
                audit.append(dict(file=f.name,year=y,month=m,nearest_cell_count=nfallback))
                for region in REGIONS:
                    mask=cross.domainID.eq(region).to_numpy()&(weights>0)
                    rows.append(dict(region=region,year=target,variable=var,month=m,value=np.average(vals[mask],weights=weights[mask])))
                county_rows.append(pd.DataFrame(dict(GEOID=cross.GEOID,year=y,month=m,variable=var,value=vals)))
        if k%50==0: print('Reduced',k,'/',len(files),flush=True)
    pd.DataFrame(rows).to_csv(PROCESSED/'regional_climate.csv',index=False)
    pd.concat(county_rows).to_csv(PROCESSED/'county_climate.csv.gz',index=False,compression='gzip')
    pd.DataFrame(audit).to_csv(PROCESSED/'raster_audit.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('step',choices=['geography','climate']);a=p.parse_args();globals()[a.step]()
