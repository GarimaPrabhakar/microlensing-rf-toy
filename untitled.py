import numpy as np
from models import Exponential
from utils import adjacency_matrix as getam
import os
import pandas as pd
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("fn", type=str, help="Filename of Candidate.")
args = vars(parser.parse_args())

dir_data = str(args["fn"])

def pspl_model(mjd, f, u0, tE, t0):
    x = mjd
    t = x - np.min(x) - (np.max(x) - np.min(x))/2 # reformat the time to have the middle at 0
    u = (u0**2 + ((t - t0)/tE)**2)**0.5 # find u(t)
    A = (u + 2)/(u*(u + 4)**0.5) # find A(t)
    A = A/min(A) # normalize A(t)
    
    return mjd, A*f, A

def main():

    df = pd.read_csv(dir_data)
    df.columns = ['index', 'x', 'y', 'flux', 'deltx', 'delty', 'chi2_', 'dx', 'dy',
           'dflux', 'qf', 'rchi2', 'fracflux', 'fluxlbs', 'dfluxlbs', 'fwhm',
           'spread_model', 'dspread_model', 'fluxiso', 'xiso', 'yiso', 'sky',
           'BADIM FLAG', 'BADFIT FLAG', 'FAILED FIT', 'FIT CHI2',
           'APERTURE RADIUS', 'QUALITY FLAG', 'STDEV FLAG', 'RELERR FLAG',
           'ref_fn', 'forced_fn', 'fn_index', 's', 'seeing', 'magzp', 'band',
           'candidate_id', 'MJD', 'resim_']
    
    df = df.iloc[np.where(df.band=="g")]
    df = df.iloc[np.where((df["QUALITY FLAG"]>4))]
    df = df.iloc[np.where((df["QUALITY FLAG"]<6))]
    df = df.iloc[np.where((df["RELERR FLAG"]==0))]
    df = df.iloc[np.where((df["chi2_"]<120))]
    df = df.iloc[np.where((df["MAGZP"]!=0))]
    df = df.iloc[np.where((df["seeing"]!=0))]
    #df = df.iloc[np.where(((df["flux"]/df["dflux"])<0.2))]
    # dfn = dfn.iloc[np.where((dfn["seeing"]>1))]



    df["normflux"] = df["flux"]*(10**((df["MAGZP"]-29)/2.5))
    df["dnormflux"] = df["dflux"]*(10**((df["MAGZP"]-29)/2.5))
    df = df.iloc[np.where((df["dnormflux"]/df["normflux"])<0.2)]
    df = df.iloc[np.where((df["normflux"]>0))]
    # df = df.iloc[np.where(df["normflux"]<50000)]
    df = df.sort_values("MJD")
    df = df.reset_index()
    
    if len(df)<7:
        print("not enough g-band datapoints!")
        return None
    u0 = round(np.random.uniform(0.1, 1), 2)
    tE = round(np.random.uniform(5, 100))
    t0 = 0

    myt, myf, _ = pspl_model(np.array(df["MJD"]), np.array(df["normflux"]), u0, tE, t0)
    
    Af = myf/np.min(myf)
    
    _, _, am_ml = getam(numbins, Af, myt)
    
    with open('dip_cows_ns_result1_30nlay.csv','a') as csvfile:
         np.savetxt(np.flatten(am_ml),delimiter=',')

    
    
