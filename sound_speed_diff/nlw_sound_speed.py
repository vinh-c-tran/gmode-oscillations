"""
    sound speed difference and chemical potential partial derivatives
    for cmf models 

    06/19/2022

    trying out something new here 
""" 
import sys
sys.path.append('../equilibrium_structure/')
from nlw_classes import * 
import pandas as pd 

""" equilibrium sound speed methods """

def total_derivative(index, x_data, f_data):
    """ calculates derivative at a point using 5 point central differences for majority of points 
        note: this only works for evenly spaced points, ie, x[i+1] - x[i] = same value for all i
        also requires to be two points away from end points
    """
    
    h = x_data[2] - x_data[1]
    
    # treating the endpoints 
    if (index == 0):
        # use 3 step forward difference formula
        return 1/(2*h)*(-3*f_data[index] + 4*f_data[index + 1] - f_data[index + 2])
    elif (index == len(x_data) - 1):
        # use 3 step backwards difference formula 
        return 1/(2*h)*(f_data[index - 2] - 4 * f_data[index - 1] + 3 * f_data[index])
    elif (index == 1) or (index == len(x_data) - 2):
        # use 3 step central difference 
        return 1/(2*h)*(f_data[index + 1] - f_data[index - 1])
    
    else:
        # use 5 step central difference formula 
        term1 = f_data[index - 2]
        term2 = -8 * f_data[index - 1]
        term3 = 8 * f_data[index + 1]
        term4 = - f_data[index + 2]
        
        return 1/(12 * h) * (term1 + term2 + term3 + term4)


# def ce2(dataframe, eos):
#     """ calculate equilibrium sound speed and returns array of ce2 for entire nb range 
#         rather than a single value 
#     """
    
#     # extract tabulated pressure and energy density data 
#     pressure_array = np.array(dataframe.loc[:, 'Pressure (MeV/fm3)'])
#     energy_density_array = np.array(dataframe.loc[:, 'Energy Density (MeV/fm3)'])
#     nb_array = np.array(dataframe.loc[:, 'nB/n0']) * eos.n0 
    
#     # empty arrays for derivatives and sound speed
#     d_pressure_array = np.zeros(len(nb_array), dtype = 'np.float64')
#     d_energy_density_array = np.zeros(len(nb_array), dtype = 'np.float64')
    
#     # calculate derivatives at each nB step 
#     for index, nb in enumerate(nb_array):
#         d_pressure_array[index] = total_derivative(index, nb_array, pressure_array)
#         d_energy_density_array[index] = total_derivative(index, nb_array, energy_density_array)
    
#     # equilibrium sound speed array 
#     ce2_array = d_pressure_array / d_energy_density_array 
    
#     return ce2_array     

def ce2(baryon_list, meson_list, lepton_list, dataframe, eos):
    """ calculate equilibrium sound speed and returns array of ce2 for entire nb range 
        rather than a single value 

        updated 06/25/2022 to treat cusps better 
    """
    
    # extract tabulated pressure and energy density data 
    pressure_array = np.array(dataframe.loc[:, 'P (MeV/fm3)'])
    energy_density_array = np.array(dataframe.loc[:, 'E Dens (MeV/fm3)'])
    nb_array = np.array(dataframe.loc[:, 'nB/n0']) * eos.n0 
    
    # discretize the things 
    regions = discretize_endpoints(baryon_list, lepton_list, dataframe)

    # re-reset again 
    reset(baryon_list, meson_list, lepton_list)

    potential_baryons = potential_baryon_gen(baryon_list) 
    potential_leptons = potential_lepton_gen(lepton_list)
    
    # empty arrays for derivatives and sound speed
    d_pressure_array = np.zeros(len(nb_array))
    d_energy_density_array = np.zeros(len(nb_array))
    
    # calculate derivatives at each nB step 
    counter = 0 
    for index, nb in enumerate(nb_array):
        
        # check if new particles appear 
        if potential_baryons != []:
            for baryon in potential_baryons:
                if dataframe.loc[index, baryon.name + ' frac'] > 0.0: 
                    potential_baryons.remove(baryon)
                    
                    counter += 1 

        if potential_leptons != []: 
            for lepton in potential_leptons:
                if dataframe.loc[index, lepton.name + ' frac'] > 0.0: 
                    potential_leptons.remove(lepton)
                        
                    counter += 1 
        
        d_pressure_array[index] = total_derivative(index - regions[counter], nb_array[regions[counter]:regions[counter+1]], pressure_array[regions[counter]:regions[counter+1]])
        d_energy_density_array[index] = total_derivative(index - regions[counter], nb_array[regions[counter]:regions[counter+1]], energy_density_array[regions[counter]:regions[counter+1]])
    
    # equilibrium sound speed array 
    ce2_array = d_pressure_array / d_energy_density_array 
    
    return ce2_array


def cs2(equilibrium_sound_speed_array, sound_speed_diff_array):
    """ given ce2 and cs2 - ce2, solve for cs2 """
    
    return sound_speed_diff_array + equilibrium_sound_speed_array 





""" methods for calculating chemical potential partial derivatives """


def data_load_tier_one(index, baryon_list, meson_list, lepton_list, eos, dataframe):
    """ loads in data from dataframe """
    
    n0 = eos.n0  # saturation density 
    
    # load in initial data like kF, frac, and number density 
    particle_list = baryon_list + lepton_list
    
    for particle in particle_list:
        particle.n0 = n0 * hc**3
        particle.nb = dataframe.loc[index, 'nB/n0'] * n0 * hc**3 
        
        for column_name in dataframe.columns:
            # store baryon fermi momenta kF
            # by searching the dataframe for column with particle and kF
            if column_name.split()[0] == particle.name and column_name.split()[1] == 'kF':
                particle.kf = dataframe.loc[index, column_name]
            # store baryon number frac 
            elif column_name.split()[0] == particle.name and column_name.split()[1] == 'frac':
                particle.frac = dataframe.loc[index, column_name]
                particle.num_density = particle.frac * dataframe.loc[index, 'nB/n0'] * n0 * hc**3
                particle.frac_array = dataframe.loc[:, column_name]
    
    # load in meson field valeus 
    for meson in meson_list:
        for column_name in dataframe.columns:
            if column_name.split()[0] == meson.name:
                meson.field_value = dataframe.loc[index, column_name]
        

        


def data_load_tier_two(baryon_list, meson_list, lepton_list, eos):
    """ calculate and store tier two data 
        assuming tier one data has already been stored
    """
    
    # # assume sigma is first meson in list 
    # sig = meson_list[0] 
    
    for baryon in baryon_list:
        baryon.set_mass_eff(sigma.field_value)
        baryon.set_ef()
        baryon.set_log_fac()
        baryon.set_scalar_density()

        # other stuff 
        baryon.set_chem_pot() 
        baryon.set_d_kf()
    
    for meson in meson_list:
        if meson != sigma:
            meson.set_partial_nb(vector_meson_sys_gen, [omega.partial_nb, sigma.partial_nb], baryon_list, eos)
    
    for lepton in lepton_list:
        lepton.set_ef()
        lepton.set_log_fac()



def data_load_tier_three(baryon_list, meson_list, lepton_list, eos):
    """ tier three data: partial derivative of chemical potential """
    
    step_solve(baryon_list, eos)

    # can now solve for/set 
    # effective mass 
    for baryon in baryon_list:
        baryon.set_d_mass_eff()
        baryon.set_d_log_fac()
        #baryon.set_d_scalar_density() 

    for baryon in baryon_list:
        baryon.set_partial_chem_pot_two()
    
    for lepton in lepton_list: 
        lepton.set_partial_chem_pot()




def sound_speed_diff(baryon_list, meson_list, lepton_list, eos, dataframe, csv_name = 'sound_speed'):
    """ 
        calculates the sound speed difference given tabulated data 

        06/26/2022: Updated to include discrete region generator 
    """
    
    # loading in all the coupling constants 

    #sigma.set_sigma_self_coupling(eos)
    for baryon in baryon_list:
        baryon.init_eos_param(eos) 
    
    # nb array: note this is nB/n0 which is unitless 
    nb_array = np.array(dataframe.loc[:,'nB/n0']) 
    cs_ce_array = np.zeros(len(nb_array))

    reset(baryon_list, meson_list, lepton_list)

    # scan and gen the discrete regions for derivatives 
    regions = discretize_endpoints(baryon_list, lepton_list, dataframe) 

    # re-reset again 
    reset(baryon_list, meson_list, lepton_list) 

    current_baryons = [Neutron, Proton] 
    current_leptons = [electron] 

    potential_baryons = potential_baryon_gen(baryon_list) 
    potential_leptons = potential_lepton_gen(lepton_list)

    particle_list = current_baryons + current_leptons
    particle_list.remove(Neutron)
    particle_list.remove(electron)

    # dictionaries to store chemical potential partial derivatives 
    chem_pot_dict = {} 
    for particle in current_baryons + current_leptons:
        chem_pot_dict[particle] = np.zeros(len(nb_array)) 
    
    chem_pot_tilde_dict = {}
    for particle in particle_list:
        chem_pot_tilde_dict[particle] = np.zeros(len(nb_array)) 


    # meson field partial derivatives 
    meson_partial_nb_dict = {} 
    for meson in meson_list:
        meson_partial_nb_dict[meson] = np.zeros(len(nb_array)) 
    
    # baryon effective energy partial derivatives
    baryon_d_eff_dict = {} 
    for baryon in current_baryons:
        baryon_d_eff_dict[baryon] = np.zeros(len(nb_array)) 

    
    counter = 0 

    
    for index, nb_n0 in enumerate(nb_array):
        
        # check if new particles appear 
        if potential_baryons != []:
            for baryon in potential_baryons:
                if dataframe.loc[index, baryon.name + ' frac'] > 0.0: 
                    current_baryons.append(baryon) 
                    particle_list.append(baryon)
                    potential_baryons.remove(baryon)
                    chem_pot_dict[baryon] = np.zeros(len(nb_array)) 
                    chem_pot_tilde_dict[baryon] = np.zeros(len(nb_array)) 

                    counter += 1 

        if potential_leptons != []: 
            for lepton in potential_leptons:
                if dataframe.loc[index, lepton.name + ' frac'] > 0.0: 
                    current_leptons.append(lepton) 
                    particle_list.append(lepton)
                    potential_leptons.remove(lepton)
                    chem_pot_dict[lepton] = np.zeros(len(nb_array))
                    chem_pot_tilde_dict[lepton] = np.zeros(len(nb_array)) 

                    counter +=1 

        # check if we need to remove a particle 
        if index != 0:
            for particle in particle_list:
                if particle.type == 'Lepton':
                    if lepton_threshold(particle) != True:
                        current_leptons.remove(particle)
                        particle_list.remove(particle)      
        

        # load in basic data from eos dataframe 
        data_load_tier_one(index, current_baryons, meson_list, current_leptons, eos, dataframe)

        # calculate base quantities 
        data_load_tier_two(current_baryons, meson_list, current_leptons, eos)

        # perform chemical potential partial derivative calculations 
        data_load_tier_three(current_baryons, meson_list, current_leptons, eos)
        
        # store chemical potential partial derivatives in dictionary 
        for particle in chem_pot_dict:
            chem_pot_dict[particle][index] = particle.d_chem_pot 
        
        # store meson partial nb values 
        for meson in meson_partial_nb_dict:
            meson_partial_nb_dict[meson][index] = meson.partial_nb 
        
        # store derviative of effective energy fields 
        for baryon in baryon_d_eff_dict:
            baryon_d_eff_dict[baryon][index] = baryon.d_ef 

        # now partial derivatives should be calculated and stored 
        coeff = (nb_n0 * eos.n0)**2 / Neutron.chem_pot
        inner_sum = 0.0
        for particle in particle_list:
            if particle  == 'muon':
                particle.mu_tilde = electron.d_chem_pot - particle.d_chem_pot 
                chem_pot_tilde_dict[particle][index] = particle.mu_tilde 
            else: 
                particle.mu_tilde = Neutron.d_chem_pot - particle.charge * electron.d_chem_pot - particle.d_chem_pot
                chem_pot_tilde_dict[particle][index] = particle.mu_tilde 
            #particle.d_frac_dnb = total_derivative(index, nb_array * eos.n0, particle.frac_array)
            particle.d_frac_dnb = total_derivative(index - regions[counter], nb_array[regions[counter]: regions[counter + 1]] * eos.n0, np.array(particle.frac_array[regions[counter]: regions[counter+1]]))
            
            inner_sum += particle.mu_tilde * particle.d_frac_dnb 
        
        sound_speed_diff = coeff * inner_sum  * hc**3 
        
        cs_ce_array[index] = sound_speed_diff 
    
    # now we parse through the data 
    # create a datamatrix 
    num_rows = len(nb_array) 
    num_columns = 2 + len(chem_pot_dict) + len(chem_pot_tilde_dict) + len(meson_partial_nb_dict) + len(baryon_d_eff_dict)
    data_array = np.zeros((num_rows, num_columns), dtype = 'float') 

    data_array[:, 0] = nb_array 
    data_array[:, 1] = cs_ce_array

    for index, particle in enumerate(chem_pot_dict):
        data_array[:, index + 2] = chem_pot_dict[particle] 
    
    for index, particle in enumerate(chem_pot_tilde_dict):
        data_array[:, index + 2 + len(chem_pot_dict)] = chem_pot_tilde_dict[particle] 
    
    for index, particle in enumerate(meson_partial_nb_dict):
        data_array[:, index + 2 + len(chem_pot_dict) + len(chem_pot_tilde_dict)] = meson_partial_nb_dict[particle] 
    
    for index, particle in enumerate(baryon_d_eff_dict):
        data_array[:, index + 2 + len(chem_pot_dict) + len(chem_pot_tilde_dict) + len(meson_partial_nb_dict)] = baryon_d_eff_dict[particle] 
    
    data_frame = pd.DataFrame(data_array, columns = sound_speed_columns(chem_pot_dict, chem_pot_tilde_dict, meson_partial_nb_dict, baryon_d_eff_dict))
    data_frame.to_csv(csv_name + '.csv', float_format = '{:.12f}'.format) 

    
    return data_frame  




""" MISC """

def potential_baryon_gen(baryon_list):
    """ returns a list of the potential baryons """
    
    pot_baryon_list = [] 
    for baryon in baryon_list:
        if baryon.name != 'Neutron' and baryon.name != 'Proton':
            pot_baryon_list.append(baryon)
    return pot_baryon_list 
 

def potential_lepton_gen(lepton_list):
    """ returns a list of the potential leptons """
    
    pot_lepton_list = [] 
    for lepton in lepton_list:
        if lepton.name != 'electron':
            pot_lepton_list.append(lepton)
    return pot_lepton_list 


def reset(baryon_list, meson_list, lepton_list):
    """ resets attributes to zero """
    particle_list = baryon_list + meson_list + lepton_list 
    for particle in particle_list:
        particle.kf = 0.0 
        particle.ef = 0.0
        particle.chem_pot = 0.0 
        particle.frac = 0.0 
        particle.num_density = 0.0 

        if particle.type == 'Baryon':
            particle.mass_eff = 0.0
            particle.scalar_density = 0.0 
            particle.log_fac = 0.0 
        
        elif particle.type == "Meson":
            particle.field_value = 0.0 


def lepton_threshold(lepton):
    """ 
    checks to see if electron chemical potential is large enough
    to support entrance of other leptons, namely, muon
    """ 
    if electron.ef >= lepton.mass:
        return True
    else:
        return False 


def sound_speed_columns(chem_pot_dict, tilde_chem_pot_dict, meson_partial_nb_dict, baryon_d_eff_dict):
    """ generate column names """

    column_names = ['nB/n0', 'cs^2 - ce^2']

    for particle in chem_pot_dict:
        column_names.append(particle.name + " d chem pot") 
    
    for particle in tilde_chem_pot_dict:
        column_names.append(particle.name + " d tilde") 
    
    for meson in meson_partial_nb_dict:
        column_names.append(meson.name + " partial nb") 
    
    for baryon in baryon_d_eff_dict:
        column_names.append(baryon.name + " d eff") 
    

    return column_names  



def discretize_endpoints(baryon_list, lepton_list, dataframe):
    """ splits the entire nB interval into distinct regions whenever a new 
        particle populates

        necessary bc otherwise our method for taking derivatives would take derivatives 
        at cuspy points which is a big no-no 
    """
    nb_array = np.array(dataframe.loc[:, 'nB/n0'])
    
    end_points = [0] 
    counter = 0 

    potential_baryons = potential_baryon_gen(baryon_list) 
    potential_leptons = potential_lepton_gen(lepton_list)
    
    for index, nb in enumerate(nb_array):
        # check if new particles appear 
        if potential_baryons != []:
            for baryon in potential_baryons:
                if dataframe.loc[index, baryon.name + ' frac'] > 0.0: 
                    potential_baryons.remove(baryon)
                    end_points.append(index)
                    counter += 1 

        if potential_leptons != []: 
            for lepton in potential_leptons:
                if dataframe.loc[index, lepton.name + ' frac'] > 0.0: 
                    potential_leptons.remove(lepton)
                    end_points.append(index) 
                    counter += 1 
    
    end_points.append(len(dataframe))
        
    return end_points



""" vector meson contribution """

def d_omega_eom(d_omega, d_rho, baryon_list, eos):
    """ derivative of omega equation of motion """

    lhs = omega.mass**2 * d_omega 

    lhs += eos.xi / 2 * Neutron.g_omega**4 * omega.field_value**2 * d_omega 

    coeff1 = 2 * eos.lambda_omega * Neutron.g_rho**2 * Neutron.g_omega**2 
    lhs_term1 = 2 * rho.field_value * omega.field_value * d_rho 
    lhs_term2 = rho.field_value**2 * d_omega 

    lhs += coeff1 * (lhs_term1 + lhs_term2) 

    rhs = 0.0 

    for baryon in baryon_list:
        rhs += baryon.g_omega * baryon.frac 

    return lhs - rhs 


def d_rho_eom(d_omega, d_rho, baryon_list, eos):
    """ derivative of rho equation of motion """

    lhs = rho.mass**2 * d_rho 

    coeff1 = 2 * eos.lambda_omega * Neutron.g_rho**2 * Neutron.g_omega**2
    lhs_term1 = omega.field_value**2 * d_rho 
    lhs_term2 = 2 * rho.field_value * omega.field_value * d_omega 

    lhs += coeff1 * (lhs_term1 + lhs_term2) 

    rhs = 0.0 

    for baryon in baryon_list:
        rhs += baryon.g_rho * baryon.isospin * baryon.frac 

    return lhs - rhs 


def vector_meson_sys_gen(input_vec, baryon_list, eos):
    """ generates system of equations for system of vector mesons """ 

    # check if input vec is correct size
    if len(input_vec) != 2:
        print("Error") 
    
    # generate system of equations 
    d_omega_eqn = d_omega_eom(input_vec[0], input_vec[1], baryon_list, eos)
    d_rho_eqn = d_rho_eom(input_vec[0], input_vec[1], baryon_list, eos) 

    system = [d_omega_eqn, d_rho_eqn]
    return system 