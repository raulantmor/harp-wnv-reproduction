"""Bounded parallel downloads; each process owns one NetCDF context/year."""
import concurrent.futures
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def year(y):
    subprocess.run([str(ROOT/'.venv/bin/python'),'src/download.py','pdsi','--year',str(y)],cwd=ROOT,check=True)
    return y

if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for result in pool.map(year,range(2004,2023)):print('Verified PDSI',result,flush=True)
