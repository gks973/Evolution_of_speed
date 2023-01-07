import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("ggplot")

import numpy as np
import scipy.optimize
from scipy.signal import find_peaks

def monoExp(x, m, t,):
    return m * np.exp(t * x)


#k=1
#T=2
#c=1

k=1000 # Change to analyze other values of k (data files must exist)

E_plus_I = 0
pm=1


figdir="figures/"


Tfastest_peaks = []
Tfastest_fit = []

slopemax_fit = []
slopemax_peaks = []

gr_max = []

cs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10]

for c in cs:
    i=0
    Ts = []
    lambdas = []
    peakslopes = []
    cstr = ""
    if c<0:
        cstr=str(np.round(c,1))
    else:
        cstr=str(c)
    for T in np.arange(2,10,0.1):
        print("Data file:", "out_k" + str(k) + "_T" + str(round(T,1)) + "_c" + cstr +".dat")
        sir_out = pd.read_csv("out_k" + str(k) + "_T" + str(round(T,1)) + "_c" + cstr +".dat",sep=" ",header=2,names=["t","S", "E", "I","R"],index_col=False)
        print("S[end] =", np.array(sir_out["S"])[-1])
        myfilter = (sir_out["S"] > 0.99) * (sir_out["I"] > 0.5e-9)
        xs = sir_out["t"][myfilter]
        ys = sir_out["I"][myfilter] + E_plus_I*sir_out["E"][myfilter]
        xs = np.array(xs)
        ys = np.array(ys)
        peaks1, peaks2 = find_peaks(ys)
        peaks_y = []
        peaks_x = []
        for peak in peaks1:
            #print("Peak", ys[peak])
            #print("at", xs[peak])
            peaks_y.append(ys[peak])
            peaks_x.append(xs[peak])
        plt.scatter(peaks_x,peaks_y)
        lastidx=-1
        nlastidx=-2
        if len(peaks1)>=abs(nlastidx):
            m = (np.log(peaks_y[lastidx])-np.log(peaks_y[nlastidx]))/(peaks_x[lastidx]-peaks_x[nlastidx])
            print("Slope:", m, "for T=",T, "with peaks at", peaks_x[lastidx], "and", peaks_x[nlastidx] )
            if i>0:
                if (m/peakslopes[-1] > 1.5 or m/peakslopes[-1] < 0.66):  # Remove numerical artefacts
                    peakslopes.append(peakslopes[-1])
                else:    
                    peakslopes.append(m)
            else:
                peakslopes.append(m)
        i+=1
        
        # perform the fit
        p0 = (1e-10, 1e-2) # start with values near those we expect
        params, cv = scipy.optimize.curve_fit(monoExp, xs, ys, p0, maxfev=60000)
        m, t = params
        print(f"lambda for T={T}:", t)
        iline = plt.plot(xs,ys,linewidth=2, alpha=0.5, label=f"T={T}, slope={round(t,4)}")
        Ts.append(T)
        lambdas.append(t)



    plt.xlabel("Time",fontweight="bold")
    plt.ylabel("Fraction of population",fontweight="bold")
    #legend = plt.legend(title="Population",loc=5,bbox_to_anchor=(1,0.25))

    #legend = plt.legend(title="Population",loc="upper left")
    #frame = legend.get_frame()
    #frame.set_facecolor("white")
    #frame.set_linewidth(0)

    plt.yscale('log')
    plt.title(f"k={k}, c={c}")
    plt.savefig(figdir + "k_" + str(k) + "c_" + cstr + "_1.png", bbox_inches='tight', pad_inches=0, dpi=175)
    #plt.show()
    plt.close()
    plt.figure()
    plt.plot(Ts, lambdas)
    plt.title(f"k={k}, c={c}")
    plt.xlabel("T")
    plt.ylabel("slope")
    print("T_fastest (based on fitting):", Ts[np.argmax(lambdas)])
    Tfastest_fit.append(Ts[np.argmax(lambdas)])
    slopemax_fit.append(lambdas[np.argmax(lambdas)])
    if len(peakslopes)>0:
        slopemax_peaks.append(peakslopes[np.argmax(peakslopes)])
    else:
        slopemax_peaks.append(0)
    plt.title("T_fastest = " +str(Ts[np.argmax(lambdas)]))
    plt.savefig(figdir + "k_" + str(k) + "c_" + cstr + "_2.png", bbox_inches='tight', pad_inches=0, dpi=175)
    plt.close()

    if len(Ts)==len(peakslopes):
        plt.figure()
        plt.plot(Ts, peakslopes)
        plt.title(f"k={k}, c={c}")
        plt.xlabel("T")
        plt.ylabel("peak_slope")
        print("T_fastest (based on peaks):", Ts[np.argmax(peakslopes)])
        Tfastest_peaks.append(Ts[np.argmax(peakslopes)])
        plt.savefig(figdir + "k_" + str(k) + "c_" + cstr + "_3.png", bbox_inches='tight', pad_inches=0, dpi=175)
        plt.title("T_fastest = " + str(Ts[np.argmax(peakslopes)]))
        #plt.show()
        plt.close()
    else:
        Tfastest_peaks.append(0)

print("c", cs)
print("T_fastest(fit)", Tfastest_fit)
print("T_fastest(peaks)", Tfastest_peaks)

Tfastest_consensus = []

for i in range(len(Tfastest_peaks)):
    if Tfastest_peaks[i]>0.1:
        Tfastest_consensus.append(Tfastest_peaks[i])
        gr_max.append(slopemax_peaks[i])
    else:
        Tfastest_consensus.append(Tfastest_fit[i])
        gr_max.append(slopemax_fit[i])

print("T_fastest(consensus)", Tfastest_consensus)
print("max slopes(consensus)", gr_max)
