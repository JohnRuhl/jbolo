import numpy as np
import yaml
from physics import *
import h5py as hp
import toml

# physical constants, etc.
h = 6.6261e-34
kB = 1.3806e-23
c = 299792458.0  # m/s
mu0 = 1.256637e-6
ep0 = 8.854188e-12
Z0 = np.sqrt(mu0/ep0)
Tcmb = 2.725

def get_atmos_from_hdf5(sim, nu):
    # Given: (site, elevation, and pwv), plus the frequencies of interest (in nu), where nu is in ghz
    #  read in (freq, tx, Tb) from hdf5 file in bolo-calc's format.
    #  and interpolate those to the frequencies in nu.
    #  returning tx,Tb at those frequencies.

    file_atmos = sim['sources']['atmosphere']['file']
    site =       sim['sources']['atmosphere']['site']
    elev =   int(sim['sources']['atmosphere']['elevation'])
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

    fill_in_optics_defaults(sim)

    # Check whether dictionary keys exist for outputs and create them if needed.
    if 'outputs' not in sim.keys():
        sim['outputs'] = {}

    for ch in sim['channels'].keys():
        # shortcut pointer
        sim_ch = sim['channels'][ch]
        chnum = sim_ch['chnum']

        if (ch not in sim['outputs'].keys()):
            sim['outputs'][ch]={}

        # shorcut pointer
        sim_out_ch = sim['outputs'][ch]
        # place to store optics results
        sim_out_ch['optics'] = {}
        sim_out_ch['sources'] = {}

        # Optics calculations
        #
        # Create frequency vector, save it to dictionary
        nu_ghz = np.arange(sim_ch['nu_low'], sim_ch['nu_high'], sim['config']['dnu'])
        nu = nu_ghz*1e9
        sim_out_ch['nu'] = np.copy(nu_ghz)
        if (sim['bolo_config']['AOmega_method'] == 'ModeCount'):
            AOmega = sim['bolo_config']['N_modes']*(c/nu)**2
        if (sim['bolo_config']['AOmega_method'] == 'Fixed'):
            AOmega = sim['bolo_config']['AOmega']*np.ones(nu)

        # Start at bolometer, then run through elements cmb to get, for each element and the sum:
        # P_optical_element, efficiency_element, cumulative_efficiency_to_just_past_that_element, each as f(nu)
        Pnu_total = np.zeros(len(nu))  # power per Hz
        effic_cumul = np.ones(len(nu))
        P_opt_cumul = 0

        # start at detector, then go through the rest of the list
        sim_out_ch['optics']['detector'] = {}
        if (sim_ch['band_response']['method'] == 'flat'):
            effic = np.ones(len(nu))*sim_ch['det_eff']
            sim_out_ch['det_bandwidth'] = (nu[-1]-nu[0])   # in Hz
            sim_out_ch['det_bandcenter'] = ((nu[-1]+nu[0])/2)
        if (sim_ch['band_response']['method'] == 'bandfile'):
            nuband_in,band_in = np.loadtxt(sim_ch['band_response']['fname'], unpack=True)
            band = np.interp(nu_ghz,nuband_in,band_in,left=0, right=0)
            sim_out_ch['det_bandwidth'] = np.trapz(band, nu_ghz*1e9)/np.max(band)
            sim_out_ch['det_bandcenter'] = np.trapz(band*nu_ghz*1e9, nu_ghz*1e9)/sim_out_ch['det_bandwidth']
            effic = band*sim_ch['det_eff']  # shaped by band, so peak is probably the indicator of interest
        sim_out_ch['optics']['detector']['effic'] = np.copy(effic)
        sim_out_ch['optics']['detector']['effic_cumul'] = np.copy(effic_cumul)  # all ones
        # Don't calculate a power of the detector on itself.
        sim_out_ch['optics']['detector']['P_opt'] = 0.0
        # Do update the cumulative efficiency, for use with the next element.
        effic_cumul *= effic

        # work outward from detector (done above) to cmb, for this pass
        # Optical elements first.
        for elem in reversed(sim['optical_elements'].keys()):
            sim_elem = sim['optical_elements'][elem]
            sim_out_ch['optics'][elem] = {}

            # Find the absorption+ Ruze_scattering, using the method appropriate for the object type
            # For each element we use 1 = T + R + S + A
            # which is to say that the losses of each mechanism add.
            obj_type = sim_elem['obj_type']
            if obj_type == 'Mirror':
                ohmic_loss = 1 - ohmic_eff(nu,float(sim_elem['conductivity']))
                ruze_loss = 1 - ruze_eff(nu,float(sim_elem['surface_rough']))
                emiss = ohmic_loss + ruze_loss
            if obj_type == 'LossTangent':
                emiss = loss_from_losstangent(nu,sim_elem['thickness'],sim_elem['index'],sim_elem['loss_tangent'])
            if obj_type == 'Bespoke':
                if (type(sim_elem['absorption']) is list):
                    emiss = sim_elem['absorption'][chnum]*np.ones(len(nu))
                else:
                    emiss = sim_elem['absorption']*np.ones(len(nu))
            if obj_type == 'ApertureStop':
                pixd = sim_ch['horn_diameter']
                fnum = sim['bolo_config']['f_number']
                wf = sim['bolo_config']['waist_factor']
                emiss = 1. - gaussian_spill_eff(nu, pixd, fnum, wf)

            # Calculate total efficiency, and P_optical given temperature, emissivity, etc.
            # All of the loss mechanisms EXCEPT REFLECTIONS are assumed to couple to the same temperature, that of this element.
            # Reflections are assumed to not couple to any tmeperature;  this appears to be what bolo-calc did, so after verifying we'll modify.
            tempdict = {}
            for item in ['reflection','scatter_frac','spillover']:
                if (type(sim_elem[item]) is list):  # if it's a list
                    tempdict[item] = sim_elem[item][chnum]
                else:
                    tempdict[item] = sim_elem[item]
            emiss_effective = emiss + tempdict['scatter_frac'] + tempdict['spillover'] #+ tempdict['reflection']
            effic = 1.0 - emiss_effective - tempdict['reflection']

            Inu = bb_spec_rad(nu, sim_elem['temperature'], emiss_effective)  # This has units of W/m^2/sr/Hz
            Pnu = Inu*AOmega*(sim['bolo_config']['N_polarizations']/2.0) # this is the power per Hz emitted by this element.
            sim_out_ch['optics'][elem]['Pnu'] = np.copy(Pnu)
            Pnu_total += Pnu*effic_cumul  # store how much of this element's Pnu gets to the bolometer
            P_opt_elem = np.trapz(effic_cumul*Pnu,nu)
            P_opt_cumul += P_opt_elem
            # calculate and store the total optical power on the bolometer from this element.
            sim_out_ch['optics'][elem]['P_opt'] = np.copy(P_opt_elem)
            sim_out_ch['optics'][elem]['effic'] = np.copy(effic)
            # The efficiency we store is *to* (not through) the relevant element.
            sim_out_ch['optics'][elem]['effic_cumul'] = np.copy(effic_cumul)
            # Update cumulative efficiency, for use with the next element.
            effic_cumul *= effic

            # instrument system (detector + optics) band center and bandwidth
            sim_out_ch['sys_bandwidth'] = np.trapz(effic_cumul,nu)/np.max(effic_cumul)
            sim_out_ch['sys_bandcenter'] = np.trapz(effic_cumul*nu/np.max(effic_cumul), nu)/sim_out_ch['sys_bandwidth']


        # Now work through the "sources", from just after telescope to cmb.
        for src in reversed(sim['sources'].keys()):
            sim_src = sim['sources'][src]
            sim_out_ch['sources'][src]={}

            # probably the cmb, but you could use it for a room temperature load, say
            if (sim_src['source_type'] == 'blackbody'):
                Inu = bb_spec_rad(nu,sim_src['T'],sim_src['emiss'])
                Pnu = Inu*AOmega*(sim['bolo_config']['N_polarizations']/2.0)
                effic = (1-sim_src['emiss'])

            # atmosphere
            if (src == 'atmosphere'):
                if (sim_src['source_type'] == 'hdf5'):
                    # hdf5 file must be organized like the bolo-calc one.
                    tx_atmos, Tb_atmos = get_atmos_from_hdf5(sim,nu_ghz)
                # Add option for file here.
                Inu = bb_spec_rad(nu, Tb_atmos, emiss=1.0)
                Pnu = Inu*AOmega*(sim['bolo_config']['N_polarizations']/2.0)
                effic = np.copy(tx_atmos)

            # Add dust and synchrotron later.

            sim_out_ch['sources'][src]['Pnu'] = np.copy(Pnu)
            Pnu_total += Pnu*effic_cumul
            P_opt_elem = np.trapz(effic_cumul*Pnu,nu)
            P_opt_cumul += P_opt_elem
            # calculate and store the total optical power on the bolometer from this element.
            sim_out_ch['sources'][src]['P_opt'] = np.copy(P_opt_elem)
            sim_out_ch['sources'][src]['effic'] = np.copy(effic)
            sim_out_ch['sources'][src]['effic_cumul'] = np.copy(effic_cumul)
            effic_cumul *= effic

            # If we just did the atmosphere, calculate the band center and width to the celestial sources, ie above the atmos.
            if (sim_src['source_type'] == 'hdf5'):
                sim_out_ch['sky_bandwidth'] = np.trapz(effic_cumul,nu)/np.max(effic_cumul)
                sim_out_ch['sky_bandcenter'] = np.trapz(effic_cumul*nu/np.max(effic_cumul), nu)/sim_out_ch['sys_bandwidth']

        # report scalars of things summed over all elements, and over frequency
        sim['outputs'][ch]['P_opt'] = np.copy(P_opt_cumul)

        # report vectors of things summed over all elements, but not frequency
        sim['outputs'][ch]['Pnu_total'] = np.copy(Pnu_total)

        # calculate the correlation factor for that channel, due to horn size, f/#, and frequencies
        fnum = sim['bolo_config']['f_number']
        dpix = sim_ch['pixel_spacing']
        lam_mean = c/np.mean(nu)
        det_pitch_flam = dpix/(fnum*lam_mean)
        aperture_factor, stop_factor = corr_facts(det_pitch_flam, flamb_max = 3.)
        sim['outputs'][ch]['aperture_factor'] = np.copy(aperture_factor)
        sim['outputs'][ch]['stop_factor'] = np.copy(stop_factor)

        # Simple photon noise first, with no pixel-pixel correlations
        NEP_photonNC, NEP_photon_poissonNC, NEP_photon_boseNC = photon_NEP_single(Pnu_total,nu)
        sim_out_ch['NEP_photonNC'] = np.copy(NEP_photonNC)
        sim_out_ch['NEP_photon_poissonNC'] = np.copy(NEP_photon_poissonNC)
        sim_out_ch['NEP_photon_boseNC'] = np.copy(NEP_photon_boseNC)
        # Now do it with pixel-pixel correlations
        Pnu_stop = sim['outputs'][ch]['optics']['lyot']['Pnu']
        Pnu_apert = Pnu_total - Pnu_stop
        NEP_photonC, NEP_photon_poissonC, NEP_photon_boseC = photon_NEP_with_corrs(Pnu_apert, Pnu_stop, aperture_factor, stop_factor, nu)
        sim_out_ch['NEP_photonC'] = np.copy(NEP_photonC)
        sim_out_ch['NEP_photon_poissonC'] = np.copy(NEP_photon_poissonC)
        sim_out_ch['NEP_photon_boseC'] = np.copy(NEP_photon_boseC)


        # Find conversion factor to NET, which is dP_opt/dT_cmb
        #   Pn_opt_cmb = eta_to_cmb*Bnu(Tcmb,nu)
        #   P_opt_cmb = trapz(Pn_opt_cmb,nu)
        #
        #   dPn_dT_cmb = eta_to_cmb*dPdT(Tcmb,nu)
        #   dP_opt/dT_cmb = trapz(dPn_dT_cmb)

        dpdt = dPdT(nu, Tcmb, sim_out_ch['sources']['cmb']['effic_cumul'], AOmega, sim['bolo_config']['N_polarizations'])
        sim_out_ch['dpdt'] = np.copy(dpdt)



    ################################


def run_bolos(sim):
    #print('Running bolos...')
    ################################
    #Bolometer calculations
    #
    #  - Psat
    #  - NEP_phonon
    #  - NEP_johnson
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

    for ch in sim['channels'].keys():
        if (ch not in sim['outputs'].keys()):
            sim['outputs'][ch]={}
        sim_ch = sim['channels'][ch]
        sim_out_ch = sim['outputs'][ch]

        if sim['bolo_config']['psat_method']=='specified':
            Psat = sim_ch['psat']

        if sim['bolo_config']['psat_method'] == 'from_optical_power':
            Psat = sim['bolo_config']['psat_factor']*sim_out_ch['P_opt']

        # Calculate NEP_phonon
        sim_out_ch['P_sat'] = np.copy(Psat)
        sim_out_ch['G_dynamic'] = Gdynamic(Psat, beta, Tbath, Tc)
        sim_out_ch['F_link'] = Flink(beta, Tbath, Tc)
        sim_out_ch['NEP_phonon'] = NEP_phonon(sim_out_ch['F_link'], sim_out_ch['G_dynamic'], Tc)


        # Calculate NEP_johnson
        sim_out_ch['NEP_johnson'] = 0.0

        # Calculate NEP_readout, must come after everything else to use "fraction" method
        NEP_NC_allbutreadout2 = sim_out_ch['NEP_phonon']**2 + sim_out_ch['NEP_photonNC']**2 + sim_out_ch['NEP_johnson']**2
        if (sim['readout']['method']=='fraction'):
            # readout noise leads to an increase in NEP_NC_total of X%
            # NEP_total**2 = (1+read_frac)**2*(NEP_photon**2 + NEP_phonon**2 + NEP_johnson**2) = (NEP_photon**2 + NEP_phonon**2 + NEP_johnson**2 + NEP_readout**2)
            # NEP_readout**2 = ((1+read_frac)**2 - 1) * (NEP_photon**2 + NEP_phonon**2 + NEP_johnson**2)
            sim_out_ch['NEP_readout'] = np.sqrt(((1+sim_ch['read_frac'])**2 - 1.0)*NEP_NC_allbutreadout2)

        # Calculate NEP_total
        sim_out_ch['NEP_NC_total'] = np.sqrt(sim_out_ch['NEP_readout']**2 + NEP_NC_allbutreadout2)  # ignore pixel correlations in bose noise
        sim_out_ch['NEP_C_total']  = np.sqrt(sim_out_ch['NEP_readout']**2 + sim_out_ch['NEP_phonon']**2 + sim_out_ch['NEP_photonC']**2 + sim_out_ch['NEP_johnson']**2)   # include pixel correlations in bose noise
        sim_out_ch['corr_factor'] = sim_out_ch['NEP_C_total']/sim_out_ch['NEP_NC_total']

        # Convert NEPs to NETs
        net_conversion_factor = 1/(sim_out_ch['dpdt']*np.sqrt(2))
        sim_out_ch['NET_NC_total'] = sim_out_ch['NEP_NC_total']*net_conversion_factor
        sim_out_ch['NET_C_total'] = sim_out_ch['NEP_C_total']*net_conversion_factor

        # Make NET_wafer, taking correlation factor and yield and detector count into consideration.
        sim_out_ch['NET_C_wafer'] = sim_out_ch['NET_C_total']/np.sqrt(sim['bolo_config']['yield']*sim_ch['num_det_per_wafer'])
        sim_out_ch['NET_NC_wafer'] = sim_out_ch['NET_NC_total']/np.sqrt(sim['bolo_config']['yield']*sim_ch['num_det_per_wafer'])

    ################################
    # Save results to a toml file.
    #


##### Print a single optics channel's optical-chain info
def print_optics(sim,ch):
    print(ch)
    print('Element            Popt(pW)   Effic  Effic_cumul')

    popttotal = 0
    for items in ['optics','sources']:
        for elem in sim['outputs'][ch][items].keys():
            _popt = sim['outputs'][ch][items][elem]['P_opt']
            popttotal += _popt
            _effic = np.mean(sim['outputs'][ch][items][elem]['effic'])
            _effic_cumul = np.mean(sim['outputs'][ch][items][elem]['effic_cumul'])
            print('{0:15s}:  {1:8.4f}   {2:8.4f}  {3:8.4f}'.format(elem, _popt*1e12,_effic,_effic_cumul))

    print('P_optical_total =  {0:8.4e}'.format(popttotal))


def print_detector(sim,ch):
    print(ch)
    sim_out_ch = sim['outputs'][ch]
    print('  P_sat:       {0:7.2f}'.format(1e12*sim_out_ch['P_sat']))
    print('  F_link:      {0:7.2f}'.format(sim_out_ch['F_link']))
    print('  G_dynamic:   {0:7.2e}'.format(sim_out_ch['G_dynamic']))
    print('  NEP_phonon:  {0:7.2f}'.format(1e18*sim_out_ch['NEP_phonon']))
    print('  NET_uncorr:  {0:7.2f}'.format(1e6*sim_out_ch['NET_NC_total']))
    print('  NET_corr:    {0:7.2f}'.format(1e6*sim_out_ch['NET_C_total']))
    print('  Corr_factor: {0:7.3f}'.format(sim_out_ch['corr_factor']))
    print('  NET_NC_wafer:   {0:7.2f}'.format(1e6*sim_out_ch['NET_NC_wafer']))
    print('  NET_C_wafer:   {0:7.2f}'.format(1e6*sim_out_ch['NET_C_wafer']))

def print_one_line(sim,param,multiplier):
    print(param.rjust(15),': ',end='')
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
        'P_opt':1e12,
        'P_sat':1e12,
        'F_link': 1.0,
        'G_dynamic': 1e12,
        'dpdt': 1e12,
        'NEP_phonon':1e18,
        'NEP_photonNC':1e18,
        'NEP_photonC':1e18,
        'NEP_readout':1e18,
        'NEP_NC_total':1e18,
        'NET_NC_total':1e6,
        'corr_factor':1.0,
        'NET_C_total':1e6,
        'NET_C_wafer':1e6}
    print((str(sim['version']['date'])+' : '+sim['version']['name']).ljust(16))
    print(' '.rjust(15),end='')
    for ch in sim['channels'].keys():
        print(ch.rjust(11),end='')
    print('')
    for param in param_list.keys():
        print_one_line(sim,param,param_list[param])

def sim_output_to_toml(mydict, filename):
    # Write to toml file
    with open(filename, "w") as toml_file:
        toml.dump(mydict, toml_file, encoder=toml.TomlNumpyEncoder())
