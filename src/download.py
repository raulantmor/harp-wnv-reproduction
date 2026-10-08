"""Download only independent public inputs and permitted scientific documents.

No author repository, institutional deposit, derived author data, or forecast
submission is accessed. Every response is hashed; mutable upstream data are
identified by content rather than by a guessed release number.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import threading
import time
from pathlib import Path
import requests
from settings import RAW, ROOT, MONTHS

LOCK = threading.Lock()
MANIFEST = ROOT / 'data/download_manifest.jsonl'

def fetch(url, filename, timeout=120):
    path = RAW / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return path
    rec = dict(url=url, file=str(path.relative_to(ROOT)),
               accessed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    try:
        # trust_env=False prevents accidental use of ~/.netrc credentials.
        session = requests.Session(); session.trust_env = False
        r = session.get(url, timeout=timeout)
        rec.update(status=r.status_code, final_url=r.url,
                   last_modified=r.headers.get('Last-Modified'), etag=r.headers.get('ETag'))
        r.raise_for_status()
        payload=r.content
        if filename.endswith('.nc') and not (payload.startswith(b'CDF') or payload.startswith(b'\x89HDF')):
            raise ValueError('Response is not a NetCDF file')
        rec.update(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
        path.write_bytes(payload)
    except Exception as e:
        rec['error']=str(e)
        raise
    finally:
        with LOCK:
            with MANIFEST.open('a') as f: f.write(json.dumps(rec)+'\n')
    return path

BASE = {
 'cdc_state.csv':'https://www.cdc.gov/wcms/vizdata/live/ncezid_dvbd/WNV/wnv_hum_historic%20-%20Historic%20Sta.csv',
 'cdc_county.csv':'https://www.cdc.gov/wcms/vizdata/live/ncezid_dvbd/WNV/wnv_hist_hum_nonhum_yearly.csv',
 'cdc-ct-config.json':'https://www.cdc.gov/west-nile-virus/data-maps/ct-2022.json',
 'county_population_centers.txt':'https://www2.census.gov/geo/docs/reference/cenpop2020/county/CenPop2020_Mean_CO.txt',
 'counties.zip':'https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_500k.zip',
 'neon.zip':'https://www.neonscience.org/sites/default/files/NEONDomains_2024.zip',
}

def base():
    for name,url in BASE.items(): fetch(url,name)
    d=json.loads((RAW/'cdc-ct-config.json').read_text())
    urls=set()
    def walk(x):
        if isinstance(x,dict):
            for v in x.values(): walk(v)
        elif isinstance(x,list):
            for v in x: walk(v)
        elif isinstance(x,str) and x.startswith('https://www.cdc.gov/') and '.csv' in x: urls.add(x)
    walk(d)
    for i,url in enumerate(sorted(urls)):
        fetch(url,'cdc_ct_%d.csv'%i)
        print('CT:',url,flush=True)

def prism():
    tasks=[]
    for y in range(2004,2023):
        for m in range(1,13):
            if m not in MONTHS or (y==2004 and m<10) or (y==2022 and m>7): continue
            for v in ['tmean','tdmean','ppt','vpdmax']:
                url=f'https://thredds.climate.ncsu.edu/thredds/fileServer/prism/monthly/{v}/{y}/prism_{v}_4km_{y}{m:02}.nc'
                tasks.append((url,f'prism/{v}_{y}{m:02}.nc'))
    failed=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        fs={pool.submit(fetch,u,n):(u,n) for u,n in tasks}
        for j,f in enumerate(concurrent.futures.as_completed(fs),1):
            try: f.result()
            except Exception as e:
                failed.append(fs[f]); print('FAILED',fs[f][1],str(e),flush=True)
            if j%40==0: print('PRISM',j,'/',len(tasks),flush=True)
    if failed: raise RuntimeError(f'{len(failed)} PRISM files failed; rerun to resume')

def pdsi(year=None):
    import netCDF4 as nc
    import numpy as np
    url='https://thredds.northwestknowledge.net/thredds/dodsC/agg_met_pdsi_1979_CurrentYear_CONUS.nc'
    source=nc.Dataset(url)
    dates=nc.num2date(source['day'][:],source['day'].units)
    var='daily_mean_palmer_drought_severity_index'
    for y in ([year] if year is not None else range(2004,2023)):
        path=RAW/f'pdsi/pdsi_{y}.nc';path.parent.mkdir(exist_ok=True)
        verified=path.with_suffix('.verified.json')
        if path.exists():
            digest=hashlib.sha256(path.read_bytes()).hexdigest()
            if verified.exists() and json.loads(verified.read_text()).get('sha256')==digest: continue
            with nc.Dataset(path) as check:
                annual=check[var][:]
                bad=any(np.std(grid.compressed())<.01 for grid in annual)
            if not bad:
                verified.write_text(json.dumps(dict(sha256=digest,check='Every pentad grid nonconstant')))
                continue
            path.rename(path.with_suffix('.rejected.nc'))
        ix=np.where(np.array([t.year==y for t in dates]))[0]
        a,b=int(ix[0]),int(ix[-1])+1
        chunks=[]
        for start in range(a,b,6):
            end=min(b,start+6)
            chunk=source[var][start:end,:,:]
            # Catch silent short/zero-filled DAP responses, observed during access.
            for i in range(len(chunk)):
                if np.std(chunk[i].compressed())<.01:
                    chunk[i]=source[var][start+i,:,:]
                if np.std(chunk[i].compressed())<.01:
                    raise ValueError(f'Constant PDSI grid at time index {start+i}')
            chunks.append(chunk)
        values=np.ma.concatenate(chunks)
        with nc.Dataset(str(path)+'.part','w') as out:
            for name,vals in [('day',source['day'][a:b]),('lat',source['lat'][:]),('lon',source['lon'][:])]:
                out.createDimension(name,len(vals));v=out.createVariable(name,'f8',(name,));v[:]=vals
                v.setncatts({k:source[name].getncattr(k) for k in source[name].ncattrs() if k!='_FillValue'})
            v=out.createVariable(var,'f4',('day','lat','lon'),zlib=True,fill_value=-9999.)
            v[:]=values;v.units='1';out.source=url;out.processing='OPeNDAP yearly subset; decoded scale and offset, no spatial averaging'
        Path(str(path)+'.part').rename(path)
        rec=dict(url=url,selection=f'{var}[{a}:1:{b-1}][0:1:584][0:1:1385]',
                 accessed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 file=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                 bytes=path.stat().st_size,processing='Decoded yearly OPeNDAP subset saved as float32 NetCDF4')
        with MANIFEST.open('a') as f:f.write(json.dumps(rec)+'\n')
        verified.write_text(json.dumps(dict(sha256=rec['sha256'],check='Every pentad grid nonconstant')))
        print('PDSI',y,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('dataset',choices=['base','prism','pdsi']);p.add_argument('--year',type=int); a=p.parse_args()
    RAW.mkdir(parents=True,exist_ok=True)
    if a.dataset=='pdsi':pdsi(a.year)
    else:globals()[a.dataset]()
