"""Constants transcribed from written specifications, never from deposited code/data."""
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw'
PROCESSED = ROOT / 'data/processed'
RESULTS = ROOT / 'results'
REGIONS = {'D01':'Northeast', 'D02':'Mid Atlantic', 'D03':'Southeast',
           'D05':'Great Lakes', 'D06':'Prairie Peninsula', 'D08':'Ozarks Complex',
           'D09':'Northern Plains', 'D10':'Central Plains', 'D11':'Southern Plains',
           'D14':'Desert Southwest', 'D17':'Pacific Southwest'}
ORDER = ['D17','D14','D09','D10','D11','D06','D08','D05','D03','D02','D01']
MONTHS = [10,11,12,1,2,3,4,5,6,7]
VARIABLES = ['tmean','tdmean','ppt','vpdmax','pdsi','soil']
GROUP = dict(tmean='air',tdmean='air',vpdmax='air',pdsi='soil',soil='soil',ppt='rain')
PUBLISHED = {'D17':['tmean_03','soil_06'], 'D14':['tmean_07','ppt_07'],
 'D09':['tmean_05','pdsi_03'], 'D10':['tdmean_06','pdsi_10'],
 'D11':['soil_10','ppt_03'], 'D06':['soil_03','tmean_05'],
 'D08':['pdsi_11','tdmean_06'], 'D05':['soil_03','tmean_01'],
 'D03':['soil_01','vpdmax_05'], 'D02':['tmean_05','ppt_07'],
 'D01':['soil_12','tmean_05']}
QUANTILES = np.array([.01,.025,.05] + list(np.arange(.1,.91,.05)) + [.95,.975,.99])
YEARS = np.arange(2005,2023)
SEED = 20261008
