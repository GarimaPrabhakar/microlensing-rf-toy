import argparse
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt


parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("fn", type=str, help="Filename of Candidate.")
parser.add_argument("numbins", type=int, help="Number of bins in adjacency matrix.")

args = vars(parser.parse_args())

def getam(numbins, f, x):
    # numbins = 4
    numbins = numbins + 1
    bins = np.linspace(min(f), max(f)+0.001, numbins)
    am = np.zeros((numbins, numbins))
    distm = np.zeros((numbins, numbins))
    countm = np.zeros((numbins, numbins))
    for counti, i in enumerate(f[:]):
        ibin = np.where(bins>i)[0][0]
        for countj, j in enumerate(f[counti:]):
            if counti==countj:
                continue
            jbin = np.where(bins>j)[0][0]
            countm[ibin-1, jbin-1] += 1
            distm[ibin-1, jbin-1] += abs(x[counti] - x[countj])
            am[ibin-1, jbin-1] += 1/abs(x[counti] - x[countj])

    return countm[:-1,:-1], distm[:-1,:-1], am[:-1,:-1]


def main():
    dir_data = str(args["fn"])
    numbins = int(args["numbins"])
    print(f"HELP ME PLEASE!! {dir_data}, {numbins}")
    df = pd.read_csv(dir_data)
    df.columns = ['index', 'x', 'y', 'flux', 'deltx', 'delty', 'chi2_', 'dx', 'dy',
           'dflux', 'qf', 'rchi2', 'fracflux', 'fluxlbs', 'dfluxlbs', 'fwhm',
           'spread_model', 'dspread_model', 'fluxiso', 'xiso', 'yiso', 'sky',
           'BADIM FLAG', 'BADFIT FLAG', 'FAILED FIT', 'FIT CHI2',
           'APERTURE RADIUS', 'QUALITY FLAG', 'STDEV FLAG', 'RELERR FLAG',
           'ref_fn', 'forced_fn', 'fn_index', 's', 'seeing', 'magzp', 'band',
           'candidate_id', 'MJD', 'resim_']
    
    
    df = df.iloc[np.where(df.band=="g")]
    df = df.iloc[np.where((df["QUALITY FLAG"]>3))]
    df = df.iloc[np.where((df["QUALITY FLAG"]<6))]
    df = df.iloc[np.where((df["RELERR FLAG"]==0))]
    df = df.iloc[np.where((df["chi2_"]<120))]
    df = df.iloc[np.where((df["flux"]>0))]
    df = df.iloc[np.where((df["magzp"]!=0))]
    df = df.iloc[np.where((df["seeing"]!=0))]
    #df = df.iloc[np.where(((df["flux"]/df["dflux"])<0.2))]
    # dfn = dfn.iloc[np.where((dfn["seeing"]>1))]
    df["normflux"] = df["flux"]*(10**((df["magzp"]-29)/2.5))
    df["dnormflux"] = df["dflux"]*(10**((df["magzp"]-29)/2.5))
    df = df.iloc[np.where((df["dnormflux"]/df["normflux"])<0.2)]
    df = df.sort_values("MJD")
    df = df.reset_index()
    
    if len(df)<10:
        print("not enough values")
        return None
    
    def pspl_model(mjd, f, u0, tE, t0):
        x = mjd
        t = x - np.min(x) - (np.max(x) - np.min(x))/2 # reformat the time to have the middle at 0
        u = (u0**2 + ((t - t0)/tE)**2)**0.5 # find u(t)
        A = (u**2 + 2)/(u*(u**2 + 4)**0.5) # find A(t)

        return mjd, A*f, A
    
    u0 = round(np.random.uniform(0.1, 1), 2)
    tE = round(np.random.uniform((max(df["MJD"]) - min(df["MJD"]))/5, 3*(max(df["MJD"]) - min(df["MJD"]))/2))
    x = df["MJD"]
    t = x - np.min(x) - (np.max(x) - np.min(x))/2 # reformat the time to have the middle at 0
    t0 = round(np.random.uniform(min(t) + 20, max(t) - 20))
    
    # MATRIX FOR INJECTED MICROLENSING OVER LC
    myt, myf, _ = pspl_model(df["MJD"], df["normflux"], u0, tE, t0)

    A = myf/np.median(myf)
    
    _, _, am_ml = getam(numbins, A, myt)
    
    mlmatrix = np.append(np.array([dir_data, "mlevent", tE, t0, u0]), am_ml.flatten())
    
    
    #MATRIX FOR ORIGINAL LC
    
    myt = df["MJD"]
    myf = df["normflux"]/np.median(df["normflux"])
    
    _, _, am_not_ml = getam(numbins, myf, myt)
    
    notmlmatrix = np.append(np.array([dir_data, "not_mlevent", 0, 0, 0]), am_not_ml.flatten())
    
    
    
    # Generate Random Noise Microlensing Event
    f = np.random.uniform(min(df["normflux"]), max(df["normflux"]), len(df["MJD"]))
    # f = f/np.median(f)
    
    #MATRIX FOR LC NOISE
    
    _, _, random_not_am_ml = getam(numbins, f/np.median(f), df["MJD"])
    
    random_not_mlmatrix = np.append(np.array([dir_data, "noise_not_mlevent", 0, 0, 0]), random_not_am_ml.flatten())
    
    #MATRIX FOR NOISE MICROLENSING EVENT
    
    myt, myf, _ = pspl_model(df["MJD"], f, u0, tE, t0)
    
    _, _, random_am_ml_orig = getam(numbins, myf/np.median(myf), myt)
    
    randommlmatrix_orig = np.append(np.array([dir_data, "noise_mlevent_orig", 0, 0, 0]), random_am_ml_orig.flatten())
    
    
    #MATRIX FOR ANOTHER RANDOM MICROLENSING EVENT
    
    u0 = round(np.random.uniform(0.1, 1), 2)
    tE = round(np.random.uniform((max(df["MJD"]) - min(df["MJD"]))/5, 3*(max(df["MJD"]) - min(df["MJD"]))/2))
    x = df["MJD"]
    t = x - np.min(x) - (np.max(x) - np.min(x))/2 # reformat the time to have the middle at 0
    t0 = round(np.random.uniform(min(t) + 20, max(t) - 20))
    
    myt, myf, _ = pspl_model(df["MJD"], f, u0, tE, t0)
    
    _, _, random_am_ml_new = getam(numbins, myf/np.median(myf), myt)
    
    randommlmatrix_new = np.append(np.array([dir_data, "noise_mlevent_new", tE, t0, u0]), random_am_ml_new.flatten())
    
    # Generate Low Noise Microlensing Event
    f = np.random.uniform(min(df["normflux"]), (max(df["normflux"]) - min(df["normflux"]))/5, len(df["MJD"]))
    # f = f/np.median(f)

    #MATRIX FOR LOW LC NOISE

    _, _, low_noise_not_ml_event = getam(numbins, f/np.median(f), df["MJD"])

    low_noise_not_mlevent = np.append(np.array([dir_data, "low_noise_not_mlevent", 0, 0, 0]), low_noise_not_ml_event.flatten())

    #MATRIX FOR LOW NOISE + INJECTED MICROLENSING EVENT

    myt, myf, _ = pspl_model(df["MJD"], f, u0, tE, t0)

    _, _, random_lownoise_am_ml_new = getam(numbins, myf/np.median(myf), myt)
    
    low_noise_randommlmatrix_new = np.append(np.array([dir_data, "low_noise_mlevent_new", tE, t0, u0]), random_lownoise_am_ml_new.flatten())
    
    
    with open(f'trainingset/microlensing.csv','a') as csvfile:
        np.savetxt(csvfile,np.reshape(np.array(mlmatrix), (len(mlmatrix), 1)).T,delimiter=',', fmt='%s')
        np.savetxt(csvfile,np.reshape(np.array(randommlmatrix_orig), (len(randommlmatrix_orig), 1)).T,delimiter=',', fmt='%s')
        np.savetxt(csvfile,np.reshape(np.array(randommlmatrix_new), (len(randommlmatrix_new), 1)).T,delimiter=',', fmt='%s')
        np.savetxt(csvfile,np.reshape(np.array(low_noise_randommlmatrix_new), (len(low_noise_randommlmatrix_new), 1)).T,delimiter=',', fmt='%s')
    
    with open(f'trainingset/not_microlensing.csv','a') as csvfile:
        np.savetxt(csvfile,np.reshape(np.array(random_not_mlmatrix), (len(random_not_mlmatrix), 1)).T,delimiter=',', fmt='%s')
        np.savetxt(csvfile,np.reshape(np.array(notmlmatrix), (len(notmlmatrix), 1)).T,delimiter=',', fmt='%s')
        np.savetxt(csvfile,np.reshape(np.array(low_noise_not_mlevent), (len(low_noise_not_mlevent), 1)).T,delimiter=',', fmt='%s')
    
    
    
if __name__ == "__main__":
    main()