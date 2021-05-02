import numpy as np

def THDSI(cleanFFT, NoisyFFT, yAxis, binsize, threshold = 2, overlapFac=0,fs = None, TF = None):
    """
    Calculates the Total Harmonic Distortion for Speech Intelligibility (THDSI) value at every time step of the Short Time Fourier Transform (STFT) spectra
    
    v1.0
    
    Input
    -----
    * cleanFFT : 2darray  
    
        Where each row is the STFT spectra at a certain time (i.e. the typical np.fft.rfft() result that is amplitude corrected) for the clean signal.
        
    * NoisyFFT : 2darray  
    
        -Where each row is the STFT spectra at a certain time (i.e. the typical np.fft.rfft() result that is amplitude corrected) for the noisy signal.
        -Must be the same shape as cleanFFT and computed from the same STFT settings (e.g. binsize)
        
    * yAxis : 1Darray
    
        Array of the  center bin frequencies resulting from the cleanSTFT and noisySTFT computations
    
    * binsize : int
    
        The binsize used in the STFT calculations for cleanSTFT and noisy STFT
        
    * Threshold : float
    
        (Optional) The amplitude threshold used to determine if a valid fundamental frequency is found. Default = 2
        
    * OverlapFac : float
    
        (Optional) Value 0 - 1 representing the percent of overlap desired Example: 0.5 means 50% overlap. Default = 0
        
    * fs : int
    
        (Optional) Sampling rate of the original data used to make cleanSTFT and noisySTFT. Including fs will result in metrics being printed about the processing parameters. Default = None
        
    * TF : 1Darray
    
        (Optional): Transfer Function applied to the noisySTFT data. Must have length equal to the row size of noisySTFT. Default = None
    
    Output
    ------
    * THDSIvals : 1darray
    
        Contains the computed THDSI values at each time step. Returns NaN if no valid frequency is found or if no harmonics are found.
        
    * THDNSIvals : 1darray
    
        Contains the computed THDNSI values at each time step. Returns NaN if no valid frequency is found or if no harmonics are found.
        
    * harmonicVals : 1darray
    
        Contains the computed summation of all harmoincs at each time step. Returns NaN if no valid frequency is found or if no harmonics are found.
        
    * fundFreqs : 1darray
    
        Contains the estimated fundamental frequency determined from the cleanSTFT spectra. Returns NaN if no valid frequency is found or if no harmonics are found.
        
    """
        
    #Check input data quality
    assert (np.shape(cleanFFT) == np.shape(NoisyFFT)),"Clean and Noisy signals are not the exact same length"
   
    #Initialize variables
    THDSIvals = []
    THDNSIvals = []
    harmonicVals = []
    fundFreqs = []
    
    #Print some metrics if given fs
    if fs:
        print("Freq res: {:.2f} [Hz]".format(fs/binsize))
        print("FFT Frame size: {:.0f} [ms]".format(binsize/fs*1000))
        if overlapFac != 0:
            print("Increment between FFTs is {} [ms] - {}% overlap".format((binsize/fs*1000)*overlapFac,overlapFac*100))
        else:
            print("Increment between FFTs is {} [ms] - {}% overlap".format((binsize/fs*1000),overlapFac*100))

    
    ACF = 1/np.mean(np.hanning(binsize)) #1/mean
    ECF = 1/np.sqrt(np.mean(np.hanning(binsize)**2)) #1/rms
    for currCol, (cleanFFTData, noisyFFTData) in enumerate(zip(cleanFFT.T, NoisyFFT.T)):

        maxIdx = cleanFFTData.argmax()

        if np.any(TF):
            #Compute corrected STFT from CNT transfer function and apply it to the noisyData before THDSI calc
            noisyFFTData = noisyFFTData/TF
                        
        fund = noisyFFTData[maxIdx]
        if (fund > (threshold*np.mean(noisyFFTData))) & (yAxis[maxIdx]>20): #Threshold check     
            harmIdx = maxIdx
            harmonics = 0
            harmonicMultiplier = 2
            while (maxIdx * harmonicMultiplier) < len(noisyFFTData): #Loop through harmonics
                harmIdx = maxIdx * harmonicMultiplier
                harmonics += noisyFFTData[harmIdx]
                harmonicMultiplier += 1
            if (harmIdx != maxIdx): #Harmonics founds, data is good
                fundFreqs.append(yAxis[maxIdx])
                harmonicVals.append(harmonics)
                THDSIvals.append(harmonics/fund*100)
                THDNSIvals.append(((np.sum(noisyFFTData)-fund)/fund)*(ECF/ACF)*100) #Remember to convert from ACF to ECF when summing to an energy value
            else: #No harmoincs found so data no good. Write NaN
                fundFreqs.append(np.nan)
                harmonicVals.append(np.nan)
                THDSIvals.append(np.nan)
                THDNSIvals.append(np.nan)
        else: #Threshold not met, write NaN
            fundFreqs.append(np.nan)
            harmonicVals.append(np.nan)
            THDSIvals.append(np.nan)
            THDNSIvals.append(np.nan)

    return THDSIvals, THDNSIvals, harmonicVals, fundFreqs