import numpy as np
import xarray as xr
from scipy.signal import detrend
from math import floor, ceil, acos

def _tukey_window(M, alpha=0.05):
    """
    Create a Tukey window (cosine taper) of length M and taper fraction alpha.
    
    Parameters
    ----------
    M : int
        Length of the window
    alpha : float
        Fraction of the window length to taper (0 = rectangular, 1 = Hann)
    
    Returns
    -------
    w : np.ndarray
        Tukey window of length M
    """
    if alpha <= 0:
        return np.ones(M)
    elif alpha >= 1:
        return 0.5 * (1 - np.cos(2 * np.pi * np.arange(M) / (M - 1)))
    
    w = np.ones(M)
    edge = int(np.floor(alpha * (M - 1) / 2.0))
    
    # Rising cosine
    for n in range(edge):
        w[n] = 0.5 * (1 - np.cos(np.pi * n / edge))
    
    # Falling cosine
    for n in range(M - edge, M):
        w[n] = 0.5 * (1 - np.cos(np.pi * (M - n - 1) / edge))
    
    return w

def kf_filter(data_in, obsPerDay, tMin, tMax, kMin, kMax, hMin, hMax, waveName):
    """
    Python version of the NCL kf_filter (applies to (time, lon) input).
    
    Parameters
    ----------
    data_in : xarray.DataArray or 2D numpy array
        Input with shape (time, lon) or dims ("time","lon"). Missing values should be np.nan.
    obsPerDay : int
        number of timesteps per day (1 for daily)
    tMin, tMax : floats
        period cutoffs in days (positive). NCL uses jMin = round(timeDim/(tMax*obsPerDay)).
    kMin, kMax : int
        wavenumber cutoffs (NCL handles positive and negative as described).
    hMin, hMax : floats or None
        equivalent depths (m). If None (or np.nan) NCL treats as missing.
    waveName : str
        "Kelvin", "ER", "MRG"/"IG0", "IG1", "IG2" (case-insensitive)
    
    Returns
    -------
    same-type as input (xarray.DataArray or numpy ndarray) with filtered data
    """
    # detect xarray
    is_xr = isinstance(data_in, xr.DataArray)
    if is_xr:
        da = data_in
        coords = dict(time=da.coords['time'].values, lon=da.coords['lon'].values)
        attrs = da.attrs.copy()
        arr = da.values.astype(float)
    else:
        arr = np.array(data_in, dtype=float)
        coords = None
        attrs = {}

    # expecting shape (time, lon)
    if arr.ndim != 2:
        raise ValueError("kf_filter expects 2D array with shape (time, lon).")
    nt, nlon = arr.shape

    # decide about wrapFlag if coords available (NCL checks lon[0]+360 == lon[-1])
    wrapFlag = False
    if is_xr:
        lon_vals = coords['lon']
        if np.isfinite(lon_vals[0]) and np.isfinite(lon_vals[-1]):
            # approximate wrap detection used in NCL
            if np.isclose(lon_vals[0] + 360.0, lon_vals[-1]):
                wrapFlag = True

    # If wrap present, remove last duplicate lon for processing (like NCL)
    if wrapFlag:
        arr_proc = arr[:, 1:nlon]
        nlon_proc = nlon - 1
    else:
        arr_proc = arr.copy()
        nlon_proc = nlon

    # detrend (linear) along time for each lon
    arr_proc = detrend(arr_proc, axis=0, type='linear')

    # taper (NCL used taper(...,0.05,0)) -> use small Tukey-like window
    taper = _tukey_window(nt, alpha=0.10)
    arr_proc = arr_proc * taper[:, None]

    # NCL fft2df stores all zonal wavenumbers but only the non-negative
    # temporal frequencies.  Reproduce that layout as (frequency, wave).
    fft2 = np.fft.ifft(np.fft.rfft(arr_proc, axis=0), axis=1) * nlon_proc

    # dims for indexing
    kDim = nlon_proc
    freqDim = fft2.shape[0]

    # jMin, jMax (NCL: round(timeDim/(tMax*obsPerDay)); round(timeDim/(tMin*obsPerDay)))
    jMin = int(round((freqDim * 1.0) / (tMax * obsPerDay)))
    jMax = int(round((freqDim * 1.0) / (tMin * obsPerDay)))
    jMax = min(jMax, freqDim)

    # iMin/ iMax (NCL logic handling negative k via kDim + k)
    if kMin < 0:
        iMin = int(round(kDim + kMin))
        iMin = max(iMin, (kDim // 2))
    else:
        iMin = int(round(kMin))
        iMin = min(iMin, (kDim // 2))

    if kMax < 0:
        iMax = int(round(kDim + kMax))
        iMax = max(iMax, (kDim // 2))
    else:
        iMax = int(round(kMax))
        iMax = min(iMax, (kDim // 2))

    # Zero out frequencies outside jMin..jMax (NCL used inclusive/exclusive indexing; follow their logic)
    fft2[:jMin, :] = 0
    if jMax < (freqDim - 1):
        fft2[jMax+1:, :] = 0

    # Zero by wavenumber range
    if iMin < iMax:
        if iMin > 0:
            fft2[:, :iMin] = 0
        if iMax < (kDim - 1):
            fft2[:, iMax+1:] = 0
    else:
        # zero inside the wrapped range
        if (iMax + 1) <= (iMin - 1):
            fft2[:, iMax+1:iMin] = 0

    # dispersion masking constants (matching NCL)
    PI = acos(-1.0)
    beta = 2.28e-11
    # convert hMin/hMax: if missing use np.nan
    hMin_val = np.nan if hMin is None else float(hMin)
    hMax_val = np.nan if hMax is None else float(hMax)
    c_vals = np.array([np.nan, np.nan])
    if not np.isnan(hMin_val):
        c_vals[0] = np.sqrt(9.8 * hMin_val)
    if not np.isnan(hMax_val):
        c_vals[1] = np.sqrt(9.8 * hMax_val)

    # spc factor as in NCL
    spc = 24.0 * 3600.0 / (2.0 * PI * obsPerDay)

    # For each spatial wavenumber index i (0..kDim-1) compute freq curve and zero outside
    for i in range(kDim):
        # nondimensional k mapping as in NCL
        if i > (kDim / 2):
            k_val = (i - kDim) * 1.0 / (6.37e6)  # negative
        else:
            k_val = i * 1.0 / (6.37e6)           # positive

        # default jWave range
        jMinWave = 0
        jMaxWave = freqDim - 1

        wave = waveName.lower()
        freq_pair = np.array([np.nan, np.nan])

        # compute freq_pair consistent with NCL forms (freq in 1/s)
        if wave == "kelvin":
            # freq = k * c  (c array length 2)
            freq_pair = k_val * c_vals
        elif wave == "er":
            freq_tmp = np.full(2, np.nan)
            for m in range(2):
                c = c_vals[m]
                if np.isnan(c):
                    continue
                denom = (k_val**2 + 3.0 * beta / c)
                # protect denom == 0
                if denom == 0:
                    freq_tmp[m] = np.nan
                else:
                    freq_tmp[m] = -beta * k_val / denom
            freq_pair = freq_tmp
        elif wave in ("mrg", "ig0"):
            freq_tmp = np.full(2, np.nan)
            for m in range(2):
                c = c_vals[m]
                if np.isnan(c):
                    continue
                if k_val == 0:
                    freq_tmp[m] = np.sqrt(beta * c)
                else:
                    term = 1.0 + 4.0 * beta / (k_val**2 * c)
                    if term < 0:
                        freq_tmp[m] = np.nan
                    else:
                        if k_val > 0:
                            freq_tmp[m] = k_val * c * (0.5 + 0.5 * np.sqrt(term))
                        else:
                            freq_tmp[m] = k_val * c * (0.5 - 0.5 * np.sqrt(term))
            freq_pair = freq_tmp
        elif wave == "ig1":
            freq_tmp = np.full(2, np.nan)
            for m in range(2):
                c = c_vals[m]
                if np.isnan(c):
                    continue
                val = 3.0 * beta * c + (k_val**2) * (c**2)
                freq_tmp[m] = np.sqrt(val) if val >= 0 else np.nan
            freq_pair = freq_tmp
        elif wave == "ig2":
            freq_tmp = np.full(2, np.nan)
            for m in range(2):
                c = c_vals[m]
                if np.isnan(c):
                    continue
                val = 5.0 * beta * c + (k_val**2) * (c**2)
                freq_tmp[m] = np.sqrt(val) if val >= 0 else np.nan
            freq_pair = freq_tmp
        else:
            # unknown wave -> skip dispersion mask
            freq_pair = np.array([np.nan, np.nan])

        # Convert freq_pair (1/s) to j index as NCL: j = floor(freq * spc * timeDim)
        if np.isnan(freq_pair[0]):
            jMinWave = 0
        else:
            jMinWave = int(floor(freq_pair[0] * spc * freqDim))
        if np.isnan(freq_pair[1]):
            jMaxWave = freqDim - 1
        else:
            jMaxWave = int(ceil(freq_pair[1] * spc * freqDim))

        # clamp
        jMinWave = max(jMinWave, 0)
        jMaxWave = min(jMaxWave, freqDim - 1)

        # zero out fft bins outside this range for this i
        if jMinWave > 0:
            fft2[:jMinWave, i] = 0
        if jMaxWave < (freqDim - 1):
            fft2[jMaxWave+1:, i] = 0

    # NCL fft2db reconstructs a real field from the half-frequency spectrum.
    recon = np.fft.irfft(np.fft.fft(fft2, axis=1) / nlon_proc, n=nt, axis=0).real

    # if wrapFlag originally true, re-insert wrap column (lon0 = last col)
    if wrapFlag:
        out = np.empty_like(arr)
        out[:, 1:nlon] = recon
        out[:, 0] = out[:, nlon-1]
    else:
        out = recon if not is_xr else recon

    # wrap into DataArray if input was DataArray
    if is_xr:
        da_out = xr.DataArray(out, coords=[coords['time'], coords['lon']], dims=['time', 'lon'])
        # attach attributes similar to NCL
        da_out.attrs.update({
            'wavenumber': (kMin, kMax),
            'period': (tMin, tMax),
            'depth': (hMin_val, hMax_val),
            'waveName': waveName
        })
        return da_out
    else:
        return out

