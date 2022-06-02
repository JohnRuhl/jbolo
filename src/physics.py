# Mostly copied from bolo-calc's (https://github.com/KIPAC/bolo-calc)
# physics.py and noise.py
# which were inherited from BoloCalc (https://github.com/chill90/BoloCalc)
# combined because they're related.
#

# Everything in here should be in SI units.

import numpy as np
import pickle as pkl

# physical constants, etc.
h = 6.6261e-34
kB = 1.3806e-23
c = 299792458.0
mu0 = 1.256637e-6
ep0 = 8.854188e-12
Z0 = np.sqrt(mu0/ep0)
Tcmb = 2.725


def gaussian_spill_eff(freq, pixd, fnum, wf=3.0):
    """
    Pixel gaussian beam coupling efficiency given a frequency [Hz],
    pixel diameter [m], f-number, and beam waist factor

    Args:
    freq : frequencies [Hz]
    pixd : pixel size [m]
    fnum : f-number
    wf : waist factor. Defaults to 3.
    """
    #freq, pixd, fnum, wf = _check_inputs(freq, [pixd, fnum, wf])
    return 1. - np.exp(
        (-np.power(np.pi, 2)/2.) * np.power(
            (pixd / (wf * fnum * (c/freq))), 2))

def ruze_eff(freq, sigma):
    """
    Ruze efficiency given a frequency [Hz] and surface RMS roughness [m]

    Args:
    freq (float): frequencies [Hz]
    sigma (float): RMS surface roughness
    """
    #freq, sigma = _check_inputs(freq, [sigma])
    return np.exp(-np.power(4 * np.pi * sigma / (c / freq), 2.))


def ohmic_eff(freq, sigma):
    """
    Ohmic efficiency given a frequency [Hz] and conductivity [S/m]

    Args:
    freq (float): frequencies [Hz]
    sigma (float): conductivity [S/m]
    """
    #freq, sigma = _check_inputs(freq, [sigma])
    loss = 4. * np.sqrt(np.pi * freq * mu0 / sigma)
    loss = loss/Z0
    effic = 1-loss
    return effic

#def dielectric_loss(freq, thick, ind, ltan):
def loss_from_losstangent(freq, thick, ind, ltan):
    """
    The dielectric loss of a substrate given the frequency [Hz],
    substrate thickness [m], index of refraction, and loss tangent

    Args:
    freq (float): frequencies [Hz]
    thick (float): substrate thickness [m]
    ind (float): index of refraction
    ltan (float): loss tangent
    """
    #freq, thick, ind, ltan = _check_inputs(freq, [thick, ind, ltan])
    return 1. - np.exp(
        (-2. * np.pi * ind * ltan * thick) / (c/freq))


def n_occ(freq, temp):
    """
    Photon occupation number given a frequency [Hz] and
    blackbody temperature [K]

    freq (float): frequency [Hz]
    temp (float): blackbody temperature [K]
    """
    #freq, temp = _check_inputs(freq, [temp])
    fact = (h * freq)/(kB * temp)
    fact = np.where(fact > 100, 100, fact)
    #with np.errstate(divide='raise'):
    return 1. / (np.exp(fact) - 1.)

def a_omega(freq):
    """
    Throughput [m^2] for a diffraction-limited detector
    given the frequency [Hz]

    Args:
    freq (float): frequencies [Hz]
    """
    #freq = _check_inputs(freq)
    return lamb(freq)**2


def corr_facts(det_pitch, flamb_max=3.):
    """
    Calculate the Bose photon-noise correlation factor, for nearby pixels

    Args:
    det_pitch (float): detector pitch in f-lambda units
    flamb_max (float): the maximum detector pitch distance
       for which to calculate the correlation factor.
       Default is 3.

    Returns:
    aperture_factor :  the factor to apply to all elements but the stop
    stop_factor:  the factor to apply to the stop.

    Requires input files in the ApertureFuncs directory:
       coherentApertCorr.pkl
       coherentStopCorr.pkl
       # NOT REQUIRED incoherentApertCorr.pkl
       # NOT REQUIRED inocherentStopCorr.pkl
    """

    geo_fact = 6   #Geometric pitch factor, for HCP.  6 for temperature, less for NEQ, NEU

    # Only use the coherent apertures
    #  These are 100 element vectors.  The one starting with "p_" is spacing between horns
    #  in units of f-lambdas.  The other vector is the correlations (a complex number), which is 1
    #  for a completely correlated horn, zero for a completely uncorrelated horn.
    #  The correlation goes to 0.5 at around 0.5flambda, and is small by flambda.
    #  So, the idea here is to run over nearby horns in an array and see how correlated they are.  The
    #  effective number of independent horns is reduced by this correlation factor, to the extent that
    #  it is greater than zero.

    # Load the pre-calculated vectors
    p_c_apert,c_apert = pkl.load(open("ApertureFuncs/coherentApertCorr.pkl", "rb"),encoding="latin1")
    p_c_stop , c_stop  = pkl.load(open("ApertureFuncs/coherentStopCorr.pkl", "rb"),encoding="latin1")

    # Detector pitch vectors for calculated correlations
    det_p = p_c_apert

    # how many pixel spacings away we can go.
    # Not sure why this is restricted to an integer.
    ndets = int(round(flamb_max / (det_pitch), 0))

    # find the indices of the pixel spacings in the pre-calculated array that are closest to the actual detector
    # spacings we're going to consider.  We start with the nearest neighbor that is one pitch element away.
    inds1 = [np.argmin(np.abs(np.array(det_p) -
            det_pitch * (n + 1)))
            for n in range(ndets)]
    #
    # Now do the same but sqrt(3) further away considering the detectors at an zig-zag on a HCP array.
    inds2 = [np.argmin(np.abs(np.array(det_p) -
            det_pitch * (n + 1) * np.sqrt(3.)))
            for n in range(ndets)]
    #
    inds = np.sort(inds1 + inds2)

    aperture_factor = np.sqrt(1. + geo_fact*(np.sum([np.abs(c_apert[ind]) for ind in inds])))
    stop_factor     = np.sqrt(1. + geo_fact*(np.sum([np.abs( c_stop[ind]) for ind in inds])))

    return aperture_factor, stop_factor


def photon_NEP_single(P_nu, freqs):
    """
    Calculate photon NEPs [W/rtHz] for a single detector, don't worry about pixel-pixel correlations

    inputs:
      P_nu (numpy array): optical power spectral density, W/Hz.
      freqs (numpy array): frequencies of observation [Hz]

    outputs:
      nep : total photon nep
      nep_poisson : poisson contribution
      nep_bose : bose contribution

    """
    nep_bose    = np.sqrt(2.*np.trapz(P_nu**2,freqs))
    nep_poisson = np.sqrt(2.*h*np.trapz(freqs*P_nu,freqs))
    nep = np.sqrt(nep_bose**2 + nep_poisson**2)
    return nep, nep_poisson, nep_bose

def photon_NEP_single_v2(Pnu,nuvec):
    """
    Alternate version of single-detector photon noise calculation.
    This works for any number of modes or polarizations, because
    noises from those those add in the square.
    Uses photon occupation number.

    NEP_photon, NEP_poisson, NEP_bose, n_avg = photon_NEP_single_v2(Pnu,freqs)
    """
    n = Pnu/(h*nuvec)  # Photon occupation number.
    A = 2*(h**2)*(nuvec**2)  # The factor of 2 here is bandwidth

    NEP_poisson = np.sqrt(np.trapz(A*n,nuvec))
    NEP_bose = np.sqrt(np.trapz(A*n**2,nuvec))
    NEP_photon = np.sqrt(NEP_poisson**2 + NEP_bose**2)
    n_avg = np.mean(n)

    return NEP_photon, NEP_poisson, NEP_bose, n_avg


def photon_NEP_with_corrs(Pnu_apert, Pnu_stop, apert_factor, stop_factor, freqs):
    """
    Calculate photon NEP [W/rtHz] for a detector

    Args:
    popts (numpy vector): total power density at each frequency [W/Hz]
    freqs (numpy vector): frequencies of observation [Hz]
    """

    # JR:  need to work on this.
    Pnu_tot = Pnu_apert + Pnu_stop

    # Poisson noise doesn't care about the pixel-pixel correlations
    nep_corr_poisson = np.sqrt(2.*h*np.trapz(freqs*Pnu_tot,freqs))

    # Bose noise does care about them.
    # test:
    nep_corr_bose = np.sqrt(2.*np.trapz((Pnu_apert*apert_factor + Pnu_stop*stop_factor)**2,freqs))  # The factor of 2 here is wrong.

    nep_corr = np.sqrt(nep_corr_poisson**2 + nep_corr_bose**2)

    return nep_corr, nep_corr_poisson, nep_corr_bose


def bb_spec_rad(freq, temp, emiss=1.0):
    """
    Blackbody spectral radiance [W/(m^2 sr Hz)] given a frequency [Hz],
    blackbody temperature [K], and blackbody emissivity

    Args:
    freq (float): frequencies [Hz]
    temp (float): blackbody temperature [K]
    emiss (float): blackbody emissivity. Defaults to 1.
    """
    return (emiss * (2 * h * (freq**3) /
            (c**2)) * n_occ(freq, temp))


######### more purely bolometer-property things below here.

def Flink(beta, Tbath, Tc):
    """
    Link factor for the bolo to the bath

    Args:
    beta : thermal carrier index
    Tbath : bath temperature [K]
    Tc : TES transition temperature [K]
    """
    return (((beta + 1)/(2 * beta + 3)) * (1 - (Tbath / Tc)**(2 * beta + 3)) /
            (1 - (Tbath / Tc)**(beta + 1)))

def Gdynamic(psat, beta, Tbath, Tc):
    """
    Dynamic Thermal conductance between the bolo and the bath,
    Gdynamic = dP/dTbolo

    Args:
    psat : saturation power [W]
    beta : thermal carrier index
    Tbath  : bath temperature [K]
    Tc  : bolo transition temperature [K]
    """
    return (psat * (beta + 1) * (Tc**beta) /
            ((Tc**(beta + 1)) - (Tbath**(beta + 1))))

def NEP_phonon(flink, Gdynamic, Tc):   # bolo-calc called this "bolo_NEP"
    """
    Thremal carrier NEP [W/rtHz]

    Args:
    flink : link factor to the bolo bath (from Tc, Tbath, and beta)
    Gdynamic : thermal conduction between the bolo and the bath [W/K]
    Tc : bolo transition temperature [K]
    """
    return np.sqrt(4 * kB * flink * (Tc**2) * Gdynamic)


# Not used in jbolo;  pull apart the responsivity and read_NEP calcs.
def read_NEP(pelec, boloR, nei, sfact=1.):
    """
    Readout NEP [W/rtHz] for a voltage-biased bolo

    Args:
    pelec (float): bias power [W]
    boloR (float): bolometer resistance [Ohms]
    nei (float): noise equivalent current [A/rtHz]
    """
    responsivity = sfact / np.sqrt(boloR * pelec)
    return nei / responsivity

def dPdT(nu, T, effic, AOmega, Npol=1):
    """
    Change in power on the detector with change in CMB temperature [W/K]

    Args:
    freqs : observation frequencies [Hz]
    T : temperature of source we're varying
    effic : cumulative efficiency from source to detector
    AOmega : of detector (lambda^2 for a single-moded detector)
    Npol : number of polarizations
    """
    _x = (h*nu)/(kB*T)
    _prefac = Npol*h**2/(kB*c**2)
    _expx = np.exp(_x)
    _dBdT = _prefac*(nu**4 / T**2)*(_expx/(_expx-1)**2)

    return np.trapz(_dBdT*AOmega*effic, nu)
