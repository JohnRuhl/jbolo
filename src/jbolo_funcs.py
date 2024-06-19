import numpy as np
import yaml
import h5py as hp
import toml
import os
from copy import copy #np.copy breaks scalars

from .physics import *

# physical constants, etc.
h = 6.6261e-34
kB = 1.3806e-23
c = 299792458.0  # m/s
mu0 = 1.256637e-6
ep0 = 8.854188e-12
Z0 = np.sqrt(mu0/ep0)
Tcmb = 2.725

def get_atmos_from_hdf5(sim, nu, **kwargs):
    # Given: (site, elevation, and pwv), plus the frequencies of interest (in nu), where nu is in ghz
    #  read in (freq, tx, Tb) from hdf5 file in bolo-calc's format.
    #  and interpolate those to the frequencies in nu.
    #  returning tx,Tb at those frequencies.

    # If the path and filename is specified for the hdf5 file, use that.
    # If not, look in the jbolo directory, where it should also be if the
    # user put it there as the README suggests.
    if 'file' in sim['sources']['atmosphere'].keys():
        file_atmos = sim['sources']['atmosphere']['file']
    else:
        jbolo_path = os.environ.get("JBOLO_PATH", "./")
        file_atmos = os.path.join(jbolo_path,'atmos/atm_20201217.hdf5')

    site =       sim['sources']['atmosphere']['site']
    elev =   int(sim['sources']['atmosphere']['elevation'])
    
    # Override sim value of pwv if really wanted
    if 'pwv_value' in kwargs.keys():
        pwv = kwargs['pwv_value']
    else:
        pwv =    int(sim['sources']['atmosphere']['pwv'])

    atmos = hp.File(file_atmos,'r')

    # elevation can be used in key as is.  We need to read in the pwv values that are multiples of 100 microns,
    # just below and just above the requested pwv, then interpolate.  (first check if we're on a multiple of 100 microns).
    remainder = np.mod(pwv,100)
    if remainder == 0:
        # we don't need to interpolate
        keyname = str(pwv)+','+str(elev)
        _nu = atmos[site][keyname][0]
        _Tb = atmos[site][keyname][2]
        _tx = atmos[site][keyname][3]
    else:
        pwv1 = np.floor_divide(pwv,100)*100
        pwv2 = (1 + np.floor_divide(pwv,100))*100
        dpwv_frac = (pwv - pwv1)/100  # fractional distance requested pwv is above lower multiple of 100
        keyname1 = str(pwv1)+','+str(elev)
        keyname2 = str(pwv2)+','+str(elev)
        #
        _nu  = atmos[site][keyname1][0]
        #
        _Tb1 = atmos[site][keyname1][2]
        _Tb2 = atmos[site][keyname2][2]
        _Tb = _Tb1 + dpwv_frac*(_Tb2-_Tb1) # linear interpolation between the two.
        #
        _tx1 = atmos[site][keyname1][3]
        _tx2 = atmos[site][keyname2][3]
        _tx = _tx1 + dpwv_frac*(_tx2-_tx1) # linear interpolation between the two.

    # Interpolate to requested frequencies
    tx = np.interp(nu,_nu,_tx)
    Tb = np.interp(nu,_nu,_Tb)

    return(tx,Tb)
    
def get_atmos_from_textfile(sim, nu_ghz):
    filename = sim['sources']['atmosphere']['file']
    try:
        # File should have three columns, whitespace delimiters.
        _nu, _tx, _Tb = np.loadtxt(filename,comments='#',unpack=True)
    except:
        print('Error loading atmosphere text file.')
        
    # Interpolate to requested frequencies
    tx = np.interp(nu_ghz,_nu,_tx)
    Tb = np.interp(nu_ghz,_nu,_Tb)
    return(tx,Tb)


# Fill in the default values of all the optical element properties.
def fill_in_optics_defaults(sim):
    for elem in sim['optical_elements'].keys():
        for prop in sim['optics_defaults'].keys():
            if prop not in sim['optical_elements'][elem].keys():
                sim['optical_elements'][elem][prop] = sim['optics_defaults'][prop]


# Run optical calculations on all channels.
def run_optics(sim):
    # Note that there is one method difference from bolo_calc
    # Here we impose for each optical element that 1 = T + R + S + A
    # where T = transmission, R = reflection, S = scattering, A=absorption
    #  whereas bolocalc seems to have set
    #    T = (1-R)*(1-A)*(1-S)
    #
    #print('Running optics')

    # Go in and set all the optics components properties to the default
    # value if they've not been set explicitly already.
    fill_in_optics_defaults(sim)

    # Check whether dictionary keys exist for outputs and create them if needed.
    if 'outputs' not in sim.keys():
        sim['outputs'] = {}

    # Loop through all the channels, doing the optics calcs for each in turn.
    for ch in sim['channels'].keys():

        # create shortcut pointer to this channels' stuff, for use throughout
        sim_ch = sim['channels'][ch]
        chnum = sim_ch['chnum']

        # Make place for outputs for this channel if it doesn't exist yet.
        if (ch not in sim['outputs'].keys()):
            sim['outputs'][ch]={}
        # Create shortcut pointer to this channel's outputs, for use throughout
        sim_out_ch = sim['outputs'][ch]
        # place to store optics results
        sim_out_ch['optics'] = {}
        sim_out_ch['sources'] = {}

        # Optics calculations
        #
        # Create frequency vector, save it to dictionary
        nu_ghz = np.arange(sim_ch['nu_low'], sim_ch['nu_high'], sim['config']['dnu'])
        nu = nu_ghz*1e9  # in Hz
        sim_out_ch['nu'] = copy(nu_ghz)  # save copy for outputs

        # Decide how to handle A*Omega;  is it a fixed number  of modes,
        # or is it set (eg geometrically) to a value?
        if (sim['bolo_config']['AOmega_method'] == 'ModeCount'):
            AOmega = sim['bolo_config']['N_modes']*(c/nu)**2
        if (sim['bolo_config']['AOmega_method'] == 'Fixed'):
            AOmega = sim['bolo_config']['AOmega']*np.ones(nu)  # meters^2 * sr

        # Start at bolometer, then run through elements out to the cmb to get,
        # for each element and the total,
        #   P_optical_element, efficiency_element, cumulative_efficiency_to_just_past_that_element, each as f(nu)
        Pnu_total = np.zeros(len(nu))  # power per Hz
        effic_cumul = np.ones(len(nu))
        P_opt_cumul = 0  # integrated over photon frequency

        # This block of code handles the detector.
        sim_out_ch['optics']['detector'] = {}
        #
        # Flat bands, normalized by det_eff
        match sim_ch['band_response']['method']:
            case 'flat':
                if ('nu_lowedge' in sim_ch['band_response'].keys()):  # check if limits given
                    _nulow = sim_ch['band_response']['nu_lowedge']
                    _nuhigh = sim_ch['band_response']['nu_highedge']
                    band = np.where((nu_ghz>=_nulow) & (nu_ghz<=_nuhigh), 1.0, 1e-6)
                    effic = band*sim_ch['det_eff']
                    sim_out_ch['det_bandwidth'] = 1e9*(_nuhigh - _nulow)   # in Hz
                    sim_out_ch['det_bandcenter'] = 1e9*(_nuhigh+_nulow)/2  # in Hz
                else:
                    band = np.ones(len(nu))
                    effic = band*sim_ch['det_eff'] # just make a flat band covering all the integration range.
                    sim_out_ch['det_bandwidth'] = (nu[-1]-nu[0])   # in Hz
                    sim_out_ch['det_bandcenter'] = ((nu[-1]+nu[0])/2)
        #
        # Read band from file, normalized by det_eff.  File expected to be in GHz.
            case 'bandfile':
                nuband_in,band_in = np.loadtxt(sim_ch['band_response']['fname'], unpack=True)
                band = np.interp(nu_ghz,nuband_in,band_in,left=0, right=0)
                sim_out_ch['det_bandwidth'] = np.trapz(band, nu)/np.max(band)
                sim_out_ch['det_bandcenter'] = np.trapz(band*nu/np.max(band), nu)/sim_out_ch['det_bandwidth']
                effic = band*sim_ch['det_eff']  # shaped by band, so peak is probably the indicator of interest
        #
        # Specify bandshape in two numpy vectors already loaded into sim, nuband_in (frequency in GHz) and band_in.
        # Normalized by det_eff
            case 'band_vector':
                nuband_in = sim_ch['band_response']['nuband_in']
                band_in =   sim_ch['band_response']['band_in']
                band = np.interp(nu_ghz,nuband_in,band_in,left=0, right=0)
                sim_out_ch['det_bandwidth'] = np.trapz(band, nu)/np.max(band)
                sim_out_ch['det_bandcenter'] = np.trapz(band*nu/np.max(band), nu)/sim_out_ch['det_bandwidth']
                effic = band*sim_ch['det_eff']  # shaped by band, so peak is probably the indicator of interest
        #
        # Logistic model, normalized by det_eff, edges specified in GHz.
            case 'logistic':
                _nulow =  sim_ch['band_response']['nu_lowedge']
                _nuhigh = sim_ch['band_response']['nu_highedge']
                _a = sim_ch['band_response']['a']
                _n = sim_ch['band_response']['n']
                band = logistic_bandmodel(nu_ghz, _nulow, _nuhigh,_a,_n)
                sim_out_ch['det_bandwidth'] = np.trapz(band, nu)/np.max(band)
                sim_out_ch['det_bandcenter'] = np.trapz(band*nu, nu)/sim_out_ch['det_bandwidth']
                effic = band*sim_ch['det_eff']
                
            case 'AlphaPowerLaw':
                # This is really for blue leaks, some of which are thought to rise with frequency.
                # band = a + (nu/nu0)^b, with a maximum of 1
                # Note that this function does *not* use the detector efficiency in the yaml file, ie sim_ch['det_eff']
                _nulow =  sim_ch['band_response']['nu_lowedge']
                _nuhigh = sim_ch['band_response']['nu_highedge']
                _a = sim_ch['band_response']['a']
                _b = sim_ch['band_response']['b']
                _nu0 = sim_ch['band_response']['nu0']
                band = _a + (nu_ghz/_nu0)**_b
                band = np.where(band<1,band,1.0) 
                
            case _:
                print('Invalid detector response method.') 

        # Finally, if anywhere we have effic = too small, things will blow up later,
        # so make the minimum efficiency 1e-6.
        effic = np.where(effic>1e-6,effic,1e-6)
    
        sim_out_ch['optics']['detector']['band'] = copy(band)
        sim_out_ch['optics']['detector']['effic'] = copy(effic)
        sim_out_ch['optics']['detector']['effic_avg'] = sim_ch['det_eff']
        sim_out_ch['optics']['detector']['effic_cumul_avg'] = 1.0
        sim_out_ch['optics']['detector']['effic_cumul'] = copy(effic_cumul)  # all ones
        # Don't calculate a power of the detector on itself.
        sim_out_ch['optics']['detector']['P_opt'] = 0.0

        # for use later in optical element efficiency calculations
        detector_efficiency_integral = np.trapz(sim_out_ch['optics']['detector']['effic'], nu)
        np.trapz(effic*                         sim_out_ch['optics']['detector']['effic'],nu)/detector_efficiency_integral

        # Do update the cumulative efficiency, for use with the next element.
        effic_cumul *= effic

        # Now, work outward from detector (done above) to cmb.
        # We'll do all the instrument's optical elements first, in this block of code,
        # then move on to the "sources" like the atmosphere and cmb.
        # We start with the element nearest the detector.
        for elem in reversed(sim['optical_elements'].keys()):
            sim_elem = sim['optical_elements'][elem]
            sim_out_ch['optics'][elem] = {}

            # Find the absorption+ Ruze_scattering, using the method appropriate for the object type
            # For each element we use 1 = T + R + S + A
            # which is to say that the losses of each mechanism add.
            #obj_type = sim_elem['obj_type']
            match sim_elem['obj_type']:
                case 'Mirror':
                    ohmic_loss = 1 - ohmic_eff(nu,float(sim_elem['conductivity']))
                    ruze_loss = 1 - ruze_eff(nu,float(sim_elem['surface_rough']))
                    emiss = ohmic_loss + ruze_loss
                    
                case 'LossTangent':
                    emiss = loss_from_losstangent(nu,sim_elem['thickness'],sim_elem['index'],sim_elem['loss_tangent'])
                    
                case 'AlphaPowerLaw':
                    emiss = loss_from_alphapowerlaw(nu,sim_elem['thickness'],sim_elem['a'],sim_elem['b'])
                    
                case 'Bespoke':
                    if (type(sim_elem['absorption']) is list):
                        emiss = sim_elem['absorption'][chnum]*np.ones(len(nu))
                    elif isinstance(sim_elem['absorption'],str):
                        nuemiss_in,emiss_in = np.loadtxt(sim_elem['absorption'], unpack=True)
                        emiss = np.interp(nu_ghz,nuemiss_in,emiss_in,left=0, right=0)
                    else:
                        emiss = sim_elem['absorption']*np.ones(len(nu))
                        
                case 'ApertureStop':
                    pixd = sim_ch['horn_diameter']
                    fnum = sim['bolo_config']['f_number']
                    wf = sim['bolo_config']['waist_factor']
                    emiss = 1. - gaussian_spill_eff(nu, pixd, fnum, wf)
                    
                case _:
                    print(f"Invalid optical element type:{sim_elem['obj_type']}")

            # Calculate total efficiency, and P_optical given temperature, emissivity, etc.
            # All of the loss mechanisms EXCEPT REFLECTIONS are assumed to couple to the same temperature, that of this element.
            # Reflections are assumed to not couple to any tmeperature;  this appears to be what bolo-calc did, so after verifying we'll modify.
            tempdict = {}
            # Assign reflection, scattering and spillover values.
            # These can be from a list (one value per channel), or from a single value
            # or band file that applies to all channels.
            for item in ['reflection','scatter_frac','spillover']:
                if (type(sim_elem[item]) is list):  # if it's a list
                    tempdict[item] = sim_elem[item][chnum]
                elif isinstance(sim_elem[item],str): # if it's a band filename
                    nuitem_in,item_in = np.loadtxt(sim_elem[item], unpack=True)
                    tempdict[item] = np.interp(nu_ghz,nuitem_in,item_in,left=0, right=0)
                else:
                    tempdict[item] = sim_elem[item]
            # Total emissivity is the sum of all the losses from the various effects.
            # There are two easy ways to treat reflections:
            #  a) as in BoloCalc, don't count reflections in the emissivity, only count them in the loss of optical efficiency.  This
            #     is equivalent to all reflections terminating on zero Kelvin.
            #  b) count the reflections as part of the emissivity, which terminates them on the temperature of the object.
            #  We choose (a) below, but may modify it later.
            emiss_effective = emiss + tempdict['scatter_frac'] + tempdict['spillover'] #+ tempdict['reflection']
            effic = 1.0 - emiss_effective - tempdict['reflection']
            Inu = bb_spec_rad(nu, sim_elem['temperature'], emiss_effective)  # This has units of W/m^2/sr/Hz, and is coded for 2 polarizations.
            Pnu = Inu*AOmega*(sim['bolo_config']['N_polarizations']/2.0) # this is the power per Hz emitted by this element.
            sim_out_ch['optics'][elem]['Pnu'] = copy(Pnu)  # save the watts/Hz for this element.
            Pnu_total += Pnu*effic_cumul  # store how much of this element's Pnu gets to the bolometer, Watts/Hz.
            P_opt_elem = np.trapz(effic_cumul*Pnu,nu)  # Watts that get to the bolometer from this elmeent.
            P_opt_cumul += P_opt_elem   # add this element to the total.
            # Store the optical power absorbed by the bolo from this element, and the efficiency vs nu, of this element
            sim_out_ch['optics'][elem]['P_opt'] = copy(P_opt_elem)
            sim_out_ch['optics'][elem]['effic'] = copy(effic)
            
            # No longer in use:  sim_out_ch['optics'][elem]['effic_avg'] = np.mean(effic) # only make sense for flat bands.
            # Calculate the "average efficiency" for this optical element, which we take to be
            # eff_element_avg = [ integral( effic_element(nu)*effic_detector(nu)) ] / [ integral( effic_detector(nu)) ]
            sim_out_ch['optics'][elem]['effic_avg'] = \
                np.trapz(effic*sim_out_ch['optics']['detector']['effic'],nu)/detector_efficiency_integral
            
            # The cumulative efficiency we store at this element is *to* (not through) the relevant element.
            sim_out_ch['optics'][elem]['effic_cumul'] = copy(effic_cumul)
            sim_out_ch['optics'][elem]['effic_cumul_avg'] = \
                np.trapz(effic_cumul,nu)/(detector_efficiency_integral/sim_ch['det_eff'])

            # Update cumulative efficiency, for use with the next element.
            effic_cumul *= effic

        # Now that we've gotten through all the instrument elements, calculate
        # some instrument-only things.
        # instrument system (detector + optics) band center and bandwidth
        sim_out_ch['sys_bandwidth'] = np.trapz(effic_cumul,nu)/np.max(effic_cumul)
        sim_out_ch['sys_bandcenter'] = np.trapz(effic_cumul*nu/np.max(effic_cumul), nu)/sim_out_ch['sys_bandwidth']
        #
        # We want to store some properties of the instrument independent of the
        # detector efficiency.  To calculate these, we take the quantities we've
        # already calculated and divide by the detector efficiency.  Note that this
        # division will blow up if we're using a detector band model that has zeros in it,
        # so avoid doing that.
        #
        # Efficiency of all the detector + optical elements:
        sim_out_ch['inst_effic_total']=copy(effic_cumul)
        #
        # Efficiency of the optical elements (not detector)
        sim_out_ch['optics_effic_total']=effic_cumul/sim_out_ch['optics']['detector']['effic']
        sim_out_ch['optics_effic_total_avg']= \
            np.trapz(sim_out_ch['optics_effic_total']*sim_out_ch['optics']['detector']['effic'],nu)/detector_efficiency_integral
        sim_out_ch['inst_effic_avg'] =  sim_out_ch['optics_effic_total_avg'] * sim_ch['det_eff']



        #
        # Pnu of all the instrument things (optics) just before the detector
        sim_out_ch['optics_Pnu_total'] = Pnu_total/sim_out_ch['optics']['detector']['effic']

        # Now work through the "sources", from just after telescope to cmb.
        for src in reversed(sim['sources'].keys()):
            sim_src = sim['sources'][src]
            sim_out_ch['sources'][src]={}

            # atmosphere
            match src:
                case 'cmb':  # For historical reasons this is the same function as the graybody.
                    Inu = bb_spec_rad(nu,sim_src['T'],sim_src['emiss'])
                    Pnu = Inu*AOmega*(sim['bolo_config']['N_polarizations']/2.0)
                    effic = (1-sim_src['emiss'])*np.ones(len(nu))
                    
                case 'graybody':
                    Inu = bb_spec_rad(nu,sim_src['T'],sim_src['emiss'])
                    Pnu = Inu*AOmega*(sim['bolo_config']['N_polarizations']/2.0)
                    effic = (1-sim_src['emiss'])*np.ones(len(nu))
                    
                case 'atmosphere':
                    match sim_src['source_type']:
                        case 'hdf5':
                            # hdf5 file must be organized like the bolo-calc one.
                            tx_atmos, Tb_atmos = get_atmos_from_hdf5(sim,nu_ghz)
                    
                            # If we can, find dPdpwv.  This works for the hdf5 cube except for the McMurdo (balloon-borne) site.
                            if sim['sources']['atmosphere']['site'] != 'McMurdo':
                                # We have enough info to calculate dP_atmos/dpwv, which can be used later for g_pwv calc.
                                _pwv1 = 100*(sim_src['pwv']//100)
                                _pwv2 = _pwv1+100
                                _tx1,_Tb1 = get_atmos_from_hdf5(sim, nu_ghz, pwv_value=_pwv1)
                                _tx2,_Tb2 = get_atmos_from_hdf5(sim, nu_ghz, pwv_value=_pwv2)
                                _Inu1 = bb_spec_rad(nu, _Tb1, emiss=1.0)
                                _Inu2 = bb_spec_rad(nu, _Tb2, emiss=1.0)
                                dPatm = (_Inu2 - _Inu1)*AOmega*(sim['bolo_config']['N_polarizations']/2.0)
                                sim_out_ch['sources'][src]['dPdpwv'] = np.trapz(10*effic_cumul*dPatm,nu)

                    # Add option for file here, not done yet.
                        case 'textfile':
                            # Three column text file, (nu_ghz, transmission, T_blackbody)
                            tx_atmos, Tb_atmos = get_atmos_from_textfile(sim,nu_ghz)
                            
                        case _:
                            print(f"Invalid atmosphere source_type:{sim_src['source_type']}")

                    sim_out_ch['sources'][src]['Tb_atmos']=copy(Tb_atmos)
                    sim_out_ch['sources'][src]['tx_atmos']=copy(tx_atmos)
                    Inu = bb_spec_rad(nu, Tb_atmos, emiss=1.0)
                    Pnu = Inu*AOmega*(sim['bolo_config']['N_polarizations']/2.0)
                    effic = copy(tx_atmos)

                case _:    # default, no match
                    print(f'Invalid source:{src}, ignoring this source.')
                    continue # Go to next source
                
                # Add dust and synchrotron later.

            sim_out_ch['sources'][src]['Pnu'] = copy(Pnu)
            Pnu_total += Pnu*effic_cumul
            P_opt_elem = np.trapz(effic_cumul*Pnu,nu)
            P_opt_cumul += P_opt_elem
            # calculate and store the total optical power on the bolometer from this element.
            sim_out_ch['sources'][src]['P_opt'] = copy(P_opt_elem)
            sim_out_ch['sources'][src]['effic'] = copy(effic)
            sim_out_ch['sources'][src]['effic_cumul'] = copy(effic_cumul)

            #     
            sim_out_ch['sources'][src]['effic_avg'] = \
                np.trapz(sim_out_ch['optics']['detector']['effic']*sim_out_ch['sources'][src]['effic'],nu)/detector_efficiency_integral
            sim_out_ch['sources'][src]['effic_cumul_avg'] = \
                np.trapz(sim_out_ch['optics']['detector']['effic']*sim_out_ch['sources'][src]['effic_cumul'],nu)/detector_efficiency_integral

            # If we just did the atmosphere, calculate the band center and width to the celestial sources, ie above the atmos.
            if src == 'atmosphere':
                sim_out_ch['sky_bandwidth'] = np.trapz(effic_cumul,nu)/np.max(effic_cumul)
                sim_out_ch['sky_bandcenter'] = np.trapz(effic_cumul*nu/np.max(effic_cumul), nu)/sim_out_ch['sky_bandwidth']
                #
                # Find a measure of total system optical efficiency, including the atmosphere.  
                # 
                # We avoid doing more integrals because the shaped band is problematic for these cumulative efficiencies.
                sim_out_ch['effic_tot_avg'] = sim_out_ch['inst_effic_avg'] * sim_out_ch['sources'][src]['effic_avg']


            
            effic_cumul *= effic
        
        # report scalars of things summed over all elements, and over frequency
        sim['outputs'][ch]['P_opt'] = copy(P_opt_cumul)

        # report vectors of things summed over all elements, but not frequency
        sim['outputs'][ch]['Pnu_total'] = copy(Pnu_total)

        # calculate the correlation factor for that channel, due to horn size, f/#, and frequencies
        fnum = sim['bolo_config']['f_number']
        dpix = sim_ch['pixel_spacing']
        lam_mean = c/sim_out_ch['det_bandcenter']
        det_pitch_flam = dpix/(fnum*lam_mean)
        sim['outputs'][ch]['det_pitch_flam']=copy(det_pitch_flam)
        aperture_factor, stop_factor = corr_facts(det_pitch_flam, flamb_max = 4.)
        sim['outputs'][ch]['aperture_factor'] = copy(aperture_factor)
        sim['outputs'][ch]['stop_factor'] = copy(stop_factor)

        # Simple photon noise first, with no pixel-pixel correlations
        NEP_photonNC, NEP_photon_poissonNC, NEP_photon_boseNC = photon_NEP_single(Pnu_total,nu)
        sim_out_ch['NEP_photonNC'] = copy(NEP_photonNC)
        sim_out_ch['NEP_photon_poissonNC'] = copy(NEP_photon_poissonNC)
        sim_out_ch['NEP_photon_boseNC'] = copy(NEP_photon_boseNC)

        # Alternate version as a check
        NEP_photonNC, NEP_photon_poissonNC, NEP_photon_boseNC, n_avg = photon_NEP_single_v2(Pnu_total,nu)
        sim_out_ch['NEP_photonNC_v2'] = copy(NEP_photonNC)
        sim_out_ch['NEP_photon_poissonNC_v2'] = copy(NEP_photon_poissonNC)
        sim_out_ch['NEP_photon_boseNC_v2'] = copy(NEP_photon_boseNC)
        sim_out_ch['n_avg']=copy(n_avg)

        # Now do it with pixel-pixel correlations
        Pnu_stop = sim['outputs'][ch]['optics']['lyot']['Pnu']*sim_out_ch['optics']['lyot']['effic_cumul']
        Pnu_apert = Pnu_total - Pnu_stop
        NEP_photonC, NEP_photon_poissonC, NEP_photon_boseC = photon_NEP_with_corrs(Pnu_apert, Pnu_stop, aperture_factor, stop_factor, nu)
        sim_out_ch['NEP_photonC'] = copy(NEP_photonC)
        sim_out_ch['NEP_photon_poissonC'] = copy(NEP_photon_poissonC)
        sim_out_ch['NEP_photon_boseC'] = copy(NEP_photon_boseC)


        # Find conversion factor to NET, which is dP_opt/dT_cmb
        #   Pn_opt_cmb = eta_to_cmb*Bnu(Tcmb,nu)
        #   P_opt_cmb = trapz(Pn_opt_cmb,nu)
        #
        #   dPn_dT_cmb = eta_to_cmb*dPdT(Tcmb,nu)
        #   dP_opt/dT_cmb = trapz(dPn_dT_cmb)

        dpdt = dPdT(nu, Tcmb, sim_out_ch['sources']['cmb']['effic_cumul'], AOmega, sim['bolo_config']['N_polarizations'])
        sim_out_ch['dpdt'] = copy(dpdt)
        #
        dpdt_rj = dPdTrj(nu,  sim_out_ch['sources']['cmb']['effic_cumul'], AOmega, sim['bolo_config']['N_polarizations'])
        sim_out_ch['dpdt_rj'] = copy(dpdt_rj)
        #
        if 'dPdpwv' in sim_out_ch['sources']['atmosphere']:
            sim_out_ch['g_pwv'] = sim_out_ch['sources']['atmosphere']['dPdpwv']/sim_out_ch['dpdt']


    ################################


def run_bolos(sim):
    #print('Running bolos...')
    ################################
    #Bolometer calculations
    #
    #  - Psat
    #  - NEP_phonon
    #  - NEP_johnson == NEP_J_tot
    #  - NEP_total
    #  - NET total
    #  - NET_wafer, ignoriing nad including pixel photon noise correlations.
    ############################

    # Check whether dictionary keys exist for outputs and create them if needed.
    if 'outputs' not in sim.keys():
        sim['outputs'] = {}

    beta = sim['bolo_config']['beta']
    Tc =   sim['bolo_config']['T_c']
    Tbath =sim['bolo_config']['T_bath']
    R_bolo = sim['bolo_config']['R_bolo']
    if 'R_shunt' in sim['bolo_config'].keys():  # for johnson noise calc.
        R_shunt = sim['bolo_config']['R_shunt']
    else:
        R_shunt = R_bolo/100.  # use this as default if not specified, ie make it tiny.

    for ch in sim['channels'].keys():
        if (ch not in sim['outputs'].keys()):
            sim['outputs'][ch]={}
        sim_ch = sim['channels'][ch]
        sim_out_ch = sim['outputs'][ch]

        # calculate (if necessary) and save Psat
        if sim['bolo_config']['psat_method']=='specified':
            Psat = sim_ch['psat']
        if sim['bolo_config']['psat_method'] == 'from_optical_power':
            Psat = sim['bolo_config']['psat_factor']*sim_out_ch['P_opt']
        sim_out_ch['P_sat'] = copy(Psat)

        # Calculate and save electrical power, and from that the bolometer voltage
        P_elec = Psat - sim_out_ch['P_opt']
        if P_elec<0:
            print('Warning:  Detector saturated, channel: {0:s}'.format(ch))
        sim_out_ch['P_elec']= copy(P_elec)
        V_bolo = np.sqrt(P_elec*R_bolo)

        # Calculate NEP_phonon.
        # This also calculates Gdynamic, which may be needed for loop gain calc (next)
        Gdyn = Gdynamic(Psat, beta, Tbath, Tc)
        sim_out_ch['G_dynamic'] = copy(Gdyn)
        if 'F_link_method' in sim['bolo_config'].keys():
            if sim['bolo_config']['F_link_method']=='specified':
                sim_out_ch['F_link'] = sim_ch['F_link']
            if sim['bolo_config']['F_link_method']=='from_beta':
                sim_out_ch['F_link'] = Flink(beta, Tbath, Tc)
        else:
            sim_out_ch['F_link'] = Flink(beta, Tbath, Tc)
        sim_out_ch['NEP_phonon'] = NEP_phonon(sim_out_ch['F_link'], Gdyn, Tc)

        # Calculate and save the loop gain
        if 'loopgain_method' in sim['bolo_config'].keys():
            if sim['bolo_config']['loopgain_method']=='infinite':
                loopgain = 1000. # close enough to infinite for us.
            if sim['bolo_config']['loopgain_method']=='specified':
                loopgain = sim_ch['loopgain']
            if sim['bolo_config']['loopgain_method']=='from_alpha':
                loopgain = (sim_ch['alpha']*P_elec)/(Gdyn*Tc)
        else:
            loopgain = 1000. # set to near-infinite if user has not specified a method
        sim_out_ch['loopgain'] = copy(loopgain)

        # Calculate the magnitude of the current responsivity.  (We ignore the sign.)
        S_I = (loopgain/(loopgain+1))/V_bolo   # amps/watt
        sim_out_ch['S_I'] = copy(S_I)

        # Calculate NEP_johnson== NEP_J_tot, from Irwin and Hilton taking chsi(I) = 1.
        # Requires so many arguments I'm just coding it here.
        NEP2_J_bolo = 4*kB*Tc*P_elec/loopgain**2
        NEP2_J_shunt =(4*kB*Tbath*P_elec*R_shunt/R_bolo) * ((loopgain-1)/(loopgain))**2
        sim_out_ch['NEP_J_bolo']  = np.sqrt(NEP2_J_bolo)
        sim_out_ch['NEP_J_shunt'] = np.sqrt(NEP2_J_shunt)
        sim_out_ch['NEP_J_tot']   = np.sqrt(NEP2_J_bolo + NEP2_J_shunt)

        # Calculate NEP_readout, must come after everything else to use "fraction" method
        NEP_NC_allbutreadout2 = sim_out_ch['NEP_phonon']**2 + sim_out_ch['NEP_photonNC']**2 + sim_out_ch['NEP_J_tot']**2
        if (sim['readout']['method']=='fraction'):
            # readout noise leads to an increase in NEP_NC_total of X%
            # NEP_total**2 = (1+read_frac)**2*(NEP_photon**2 + NEP_phonon**2 + NEP_J_tot**2) = (NEP_photon**2 + NEP_phonon**2 + NEP_J_tot**2 + NEP_readout**2)
            # NEP_readout**2 = ((1+read_frac)**2 - 1) * (NEP_photon**2 + NEP_phonon**2 + NEP_J_tot**2)
            sim_out_ch['NEP_readout'] = np.sqrt(((1+sim_ch['read_frac'])**2 - 1.0)*NEP_NC_allbutreadout2)
            sim_out_ch['read_frac']=sim_ch['read_frac']
        if (sim['readout']['method']=='from_NEI'):
            sim_out_ch['NEP_readout'] = sim_ch['readout_NEI']/S_I  # convert NEI to NEP
            sim_out_ch['read_frac']= np.sqrt((NEP_NC_allbutreadout2 + sim_out_ch['NEP_readout']**2)/NEP_NC_allbutreadout2) - 1

        # Calculate NEP_total and NEP_dark (latter is without phonon noise)
        # Subscript "NC" means "Not using horn-horn photon noise correlations"
        # Subscript "C" means "using horn-horn photon noise correlations"
        # corr_factor is the ratio of those two, and depends on the relative size of various NEP contributions.
        sim_out_ch['NEP_NC_total'] = np.sqrt(sim_out_ch['NEP_readout']**2 + NEP_NC_allbutreadout2)  # ignore pixel correlations in bose noise
        sim_out_ch['NEP_C_total']  = np.sqrt(sim_out_ch['NEP_readout']**2 + sim_out_ch['NEP_phonon']**2 + sim_out_ch['NEP_photonC']**2 + sim_out_ch['NEP_J_tot']**2)   # include pixel correlations in bose noise
        sim_out_ch['corr_factor'] = sim_out_ch['NEP_C_total']/sim_out_ch['NEP_NC_total']
        #
        sim_out_ch['NEP_dark'] = np.sqrt(sim_out_ch['NEP_readout']**2 + sim_out_ch['NEP_phonon']**2 + sim_out_ch['NEP_J_tot']**2)

        # Convert NEPs [in Watts/sqrt(Hz)]to single-detector NETs [in K_cmb*sqrt(s)]
        net_conversion_factor    = 1/(sim_out_ch['dpdt']   *np.sqrt(2))
        netrj_conversion_factor  = 1/(sim_out_ch['dpdt_rj']*np.sqrt(2))
        sim_out_ch['NET_NC_total'] = sim_out_ch['NEP_NC_total']*net_conversion_factor
        sim_out_ch['NET_C_total'] = sim_out_ch['NEP_C_total']*net_conversion_factor
        sim_out_ch['NETrj_NC_total'] = sim_out_ch['NEP_NC_total']*netrj_conversion_factor
        sim_out_ch['NETrj_C_total'] = sim_out_ch['NEP_C_total']*netrj_conversion_factor


        # Make NET_wafer, taking correlation factor and yield and detector count into consideration.
        sim_out_ch['NET_C_wafer'] =  sim_out_ch['NET_C_total']/np.sqrt(sim['bolo_config']['yield']*sim_ch['num_det_per_wafer'])
        sim_out_ch['NET_NC_wafer'] = sim_out_ch['NET_NC_total']/np.sqrt(sim['bolo_config']['yield']*sim_ch['num_det_per_wafer'])
        sim_out_ch['NETrj_C_wafer'] =  sim_out_ch['NETrj_C_total']/np.sqrt(sim['bolo_config']['yield']*sim_ch['num_det_per_wafer'])
        sim_out_ch['NETrj_NC_wafer'] = sim_out_ch['NETrj_NC_total']/np.sqrt(sim['bolo_config']['yield']*sim_ch['num_det_per_wafer'])


##### Print a single optics channel's optical-chain info
def print_optics(sim,ch):
    print(ch)
    print('Element            Popt(pW)   Effic  Effic_cumul')

    popttotal = 0
    poptonly_total = 0
    for items in ['optics','sources']:
        for elem in sim['outputs'][ch][items].keys():
            _popt = sim['outputs'][ch][items][elem]['P_opt']
            popttotal += _popt
            if items == 'optics':
                poptonly_total += _popt
            _effic =       sim['outputs'][ch][items][elem]['effic_avg']
            _effic_cumul = sim['outputs'][ch][items][elem]['effic_cumul_avg']
            print('{0:15s}:  {1:8.4f}   {2:8.4f}  {3:8.4f}'.format(elem, _popt*1e12,_effic,_effic_cumul))

    print('P_opticsonly_total = {0:8.4e}'.format(poptonly_total))
    print('P_optical_total =  {0:8.4e}'.format(popttotal))


def print_detector(sim,ch):
    print(ch)
    sim_out_ch = sim['outputs'][ch]
    print('  P_sat:       {0:7.2f}'.format(1e12*sim_out_ch['P_sat']))
    print('  P_elec:      {0:7.2f}'.format(1e12*sim_out_ch['P_elec']))
    print('  F_link:      {0:7.2f}'.format(sim_out_ch['F_link']))
    print('  G_dynamic:   {0:7.2e}'.format(sim_out_ch['G_dynamic']))
    print('  Loop gain:   {0:7.2e}'.format(sim_out_ch['loopgain']))
    print('  NEP_phonon:  {0:7.2f}'.format(1e18*sim_out_ch['NEP_phonon']))
    print('  NEP_photon:  {0:7.2f}'.format(1e18*sim_out_ch['NEP_photonNC']))
    print('  NEP_J_total: {0:7.2f}'.format(1e18*sim_out_ch['NEP_J_total']))
    print('  NEP_readout: {0:7.2f}'.format(1e18*sim_out_ch['NEP_readout']))
    print('  NEP_total:   {0:7.2f}'.format(1e18*sim_out_ch['NEP_NC_total']))
    print('  NET_uncorr:  {0:7.2f}'.format(1e6*sim_out_ch['NET_NC_total']))
    print('  NET_NC_wafer:   {0:7.2f}'.format(1e6*sim_out_ch['NET_NC_wafer']))
    print('  With horn correlations:')
    print('  Corr_factor: {0:7.3f}'.format(sim_out_ch['corr_factor']))
    print('  NET_corr:    {0:7.2f}'.format(1e6*sim_out_ch['NET_C_total']))
    print('  NET_C_wafer:   {0:7.2f}'.format(1e6*sim_out_ch['NET_C_wafer']))

def print_one_line(sim,param,multiplier):
    print(param.rjust(22),': ',end='')
    for ch in sim['channels'].keys():
        print('{0:9.3f}'.format(sim['outputs'][ch][param]*multiplier),end='  ')
    print('')

def print_full_table(sim):
    param_list = {
        'det_bandcenter':1e-9,
        'det_bandwidth':1e-9,
        'sys_bandcenter':1e-9,
        'sys_bandwidth':1e-9,
        'sky_bandcenter':1e-9,
        'sky_bandwidth':1e-9,
        'optics_effic_total_avg':1.0,
        'inst_effic_avg':1.0,
        'effic_tot_avg':1.0,
        'P_opt':1e12,
        'n_avg':1.0,
        'P_elec':1e12,
        'P_sat':1e12,
        'loopgain':1.0,
        'F_link': 1.0,
        'G_dynamic': 1e12,
        'dpdt': 1e12,
        'dpdt_rj': 1e12,
        'NEP_readout':1e18,
        'NEP_J_tot':1e18,
        'NEP_phonon':1e18,
        'NEP_dark':1e18,
        'NEP_photonNC':1e18,
        'NEP_NC_total':1e18,
        'NETrj_NC_total':1e6,
        'NET_NC_total':1e6,
        'corr_factor':1.0,
        'NET_C_total':1e6,
        'NEP_photonC':1e18,
        'NET_NC_wafer':1e6,
        'NET_C_wafer':1e6}
    print((str(sim['version']['date'])+' : '+sim['version']['name']).ljust(16))
    print(' '.rjust(22),end='')
    for ch in sim['channels'].keys():
        print(ch.rjust(11),end='')
    print('')
    for param in param_list.keys():
        print_one_line(sim,param,param_list[param])


def print_lyot_efficiencies(sim):
    # Special case, go in and print 'lyot'  efficiencies on one line
    print('lyot effic_avg'.rjust(22),': ',end='')
    for ch in sim['channels'].keys():
        lyot_effic = sim['outputs'][ch]['optics']['lyot']['effic_avg']
        print('{0:9.3f}'.format(lyot_effic),end='  ')
    print('')

# A function that returns a detector (or other) passband that has soft
# edges...
def logistic_bandmodel(nuvec, nu_lowedge, nu_highedge,a,n):
    '''Returns a logistic-function band model,
       peaking at unity, given input central frequency
       and bandwidth in GHz.

       nuvec:  numpy vector of frequencies at which to return the function
       nu_lowedge:  50% power point of low edge of band.
       nu_highedge: 50% power point of high edge of band.
       a: prefactor used in logistic model
       n: exponent used in logistic model

       band = f1*f2, where f1 is a highpass, f2 is a lowpass
       k1 = a*(20/lowedge)**n
       k2 = a*(20/highedge)**m
       f1 =     1/(1+np.exp(-k1*(nu-nu_lowedge )))
       f2 = 1 - 1/(1+np.exp(-k2*(nu-nu_highedge)))

       Note that a = 2, n=0.7 gives "reasonable" band edge
       transitions for CMB-S4.
    '''
    k1 = a*(20/nu_lowedge)**n
    f1 = 1/(1+np.exp(-k1*(nuvec-nu_lowedge )))
    k2 = a*(20/nu_highedge)**n
    f2 = 1-1/(1+np.exp(-k2*(nuvec-nu_highedge)))
    f_band = f1*f2
    f_band = f_band/np.max(f_band) # normalize to unity
    return f_band


# this doesn't work, sim is too complicated now.
def sim_output_to_toml(mydict, filename):
    # Write to toml file
    with open(filename, "w") as toml_file:
        toml.dump(mydict, toml_file, encoder=toml.TomlNumpyEncoder())
