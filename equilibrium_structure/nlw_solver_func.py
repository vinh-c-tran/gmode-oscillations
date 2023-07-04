from nlw_classes import * 
import pandas as pd
from scipy.optimize import root 

"""
    The functions in this file accomplishes the following tasks
        1. Sets up the system of non-linear equations that we have to solve in order to determine the 
            meson field values and baryon and lepton momenta (ie, fractions) as a function of baryon number density nB 
        2. Solves for the pressure and energy density as well 
        3. Returns the data in a .csv file with appropriately named columns

    The main function that would need to be called is just `full_solve`. Everything else is used to set up the problem. 

    In changing from model to model, the meson equations of motion would need to be updated by calculating them in the mean field approximation
    for the given Lagrangian. 

    In the following the functions make reference to `fermi_var_list`. This is a 2 x N array where the N is the number of baryons and leptons. 
    1. The first row corresponds to numeric values for the particle fermi momentum kF (initial guesses to pass into the non-linear solver). 
    2. The second row corresponds to the particle object which itself stores information such as mass and charge and so on. 

    Example. 
    if fermi_var_list = [0.3, 0.5]
                     [Neutron, Proton]
    then fermi_var_list[1][1].mass would return Proton.mass = 939.0 or so. 
    
"""

def sigma_eom(sigma_field, fermi_var_list, eos):
    """ sigma equation of motion """
    
    # left hand side 
    lhs = sigma.mass**2 * sigma_field
    lhs += eos.k / 2 * Neutron.g_sigma**3 * sigma_field**2 
    lhs += eos.l / 6 * Neutron.g_sigma**4 * sigma_field**3 
    
    # right hand side 
    rhs = 0.0
    
    for index, particle in enumerate(fermi_var_list[1]):
        if particle.type == 'Baryon':
            rhs += particle.g_sigma * particle.set_scalar_density(fermi_var_list[0][index], sigma_field)
    
    return lhs - rhs 


def omega_eom(omega_field, rho_field, fermi_var_list, eos):
    """ omega equation of motion in typical NL3 models """
    
    lhs = omega.mass**2 * omega_field
    lhs += eos.xi / (3 * 2) * Neutron.g_omega**4 * omega_field**3 
    lhs += 2 * eos.lambda_omega * Neutron.g_rho**2 * Neutron.g_omega**2 * rho_field**2 * omega_field

    rhs = 0.0 
    
    for index, particle in enumerate(fermi_var_list[1]):
        if particle.type == 'Baryon':
            rhs += particle.g_omega * fermi_var_list[0][index]**3 / (3 * np.pi**2)
        
    return lhs - rhs 


def rho_eom(omega_field, rho_field, fermi_var_list, eos):
    """ rho equation of motion """
    
    lhs = rho.mass**2 * rho_field 
    lhs += 2 * eos.lambda_omega * Neutron.g_rho**2 * Neutron.g_omega**2 * rho_field * omega_field**2 

    rhs = 0.0
    
    for index, particle in enumerate(fermi_var_list[1]):
        if particle.type == 'Baryon':
            rhs += particle.isospin * particle.g_rho * fermi_var_list[0][index]**3/(3* np.pi**2)
        
    return lhs - rhs 


def phi_eom(phi_field, fermi_var_list):
    """ rho equation of motion """
    
    lhs = phi.mass**2 * phi_field 
    rhs = 0.0
    
    for index, particle in enumerate(fermi_var_list[1]):
        if particle.type == 'Baryon':
            rhs += particle.g_phi * fermi_var_list[0][index]**3 / (3 * np.pi**2) 
        
    return lhs - rhs 


def charge_neutrality(fermi_var_list):
    """ generates charge neutrality equation """
    
    lhs, rhs = 0.0, 0.0 
    
    for index, particle in enumerate(fermi_var_list[1]):
        if particle.charge < 0.0:
            lhs += fermi_var_list[0][index]**3
        elif particle.charge > 0.0:
            rhs += fermi_var_list[0][index]**3
    
    return rhs - lhs



def baryon_num_conservation(fermi_var_list, nB):
    """ generates baryon conservation condition """
    
    lhs = 3 * pi**2 * nB
    rhs = 0.0 
    
    for index, particle in enumerate(fermi_var_list[1]):
        if particle.type == 'Baryon':
            rhs += fermi_var_list[0][index]**3
    
    return lhs - rhs 



def beta_equilibrium(fermi_var_list, meson_var_list):
    """ generates the (possibly more than 1) beta equilibrium conditions """
    
    beta_equilibrium_array = [] 
    
    # get the neutron and electron kf values explicitly 
    for index, particle in enumerate(fermi_var_list[1]):
        if particle.name == 'Neutron':
            neutron_kf = fermi_var_list[0][index] 
        elif particle.name == 'electron':
            electron_kf = fermi_var_list[0][index]
    
    
    # generate and append the other beta equilibrium conditions 
    for index, particle in enumerate(fermi_var_list[1]):
        
        if particle.type == 'Lepton' and particle.name != 'electron':
            eqn = particle.solver_chem_pot(fermi_var_list[0][index]) - electron.solver_chem_pot(electron_kf)
            beta_equilibrium_array.append(eqn)
            
            
        if particle.type == 'Baryon' and particle.name != 'Neutron':
            eqn = Neutron.solver_chem_pot(neutron_kf, meson_var_list) - particle.charge * electron.solver_chem_pot(electron_kf) - particle.solver_chem_pot(fermi_var_list[0][index], meson_var_list)
            beta_equilibrium_array.append(eqn)
    
    return beta_equilibrium_array 




def sys_eqn_new(input_vec, particle_list, meson_list, nB, eos):
    """ generate system of equations 
        new version taking particle list as argument with particles ordered in order of appearance
    """
    
    # check if input vec is correct size 
    if len(input_vec) != len(particle_list + meson_list):
        print("Error in input vector size" )

    # mesonic variable list 
    meson_var_list = [[],[]] 
    for i in range(len(meson_list)):
        meson_var_list[0].append(input_vec[i])
        meson_var_list[1].append(meson_list[i])
        
    #fermi momenta variable list 
    fermi_momenta_list = [[],[]]
    for i in range(len(particle_list)):
        fermi_momenta_list[0].append(input_vec[i + len(meson_list)])
        fermi_momenta_list[1].append(particle_list[i])
    
    
    # generate system of equations 
    sigma_eqn = sigma_eom(input_vec[0], fermi_momenta_list, eos)
    omega_eqn = omega_eom(input_vec[1], input_vec[2], fermi_momenta_list, eos)
    rho_eqn = rho_eom(input_vec[1], input_vec[2], fermi_momenta_list, eos)
    phi_eqn = phi_eom(input_vec[3], fermi_momenta_list)
                   
    charge_neu = charge_neutrality(fermi_momenta_list)
    baryon_cons = baryon_num_conservation(fermi_momenta_list, nB)
    beta_equil = beta_equilibrium(fermi_momenta_list, meson_var_list)
                   
                   
    # putting everything back together 
    system = [sigma_eqn, omega_eqn, rho_eqn, phi_eqn, charge_neu, baryon_cons] +  beta_equil 
    return system 


""" helper functions for implementing the full solve
    full solve means iterating through the entire range of nB values 
    so we need to implement a way for new particles to appear as densities get higher here 
        - threshold equations
        - manipulating baryon, lepton lists 
"""

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



def bare_chemical_potential(baryon, meson_list):
    """finds bare chemical potential for a baryon given a list of mesons
    assumes meson objects are filled with field values
    this code is ugly and is not easily generalizable... """
    
    bare_chem = 0.0 
    for meson in meson_list:
        if (meson == sigma):
            bare_chem += baryon.mass - baryon.g_sigma * meson.field_value
        elif (meson == omega):
            bare_chem += baryon.g_omega * meson.field_value
        elif (meson == rho):
            bare_chem += baryon.g_rho * baryon.isospin * meson.field_value
        elif (meson == phi):
            bare_chem += baryon.g_phi * meson.field_value
    return bare_chem 


def baryon_threshold(baryon, meson_list):
    """ checks to see if combination of neutron and electron chemical potential """ 
    
    neutron_chem_pot = Neutron.chem_pot 
    electron_chem_pot = electron.chem_pot
    
    if (neutron_chem_pot - baryon.charge * electron_chem_pot) > (bare_chemical_potential(baryon, meson_list) + 1.0):
        return True
    else:
        return False 
    

def lepton_threshold(lepton):
    """ 
    checks to see if electron chemical potential is large enough
    to support entrance of other leptons, namely, muon
    """ 
    if electron.ef >= lepton.mass:
        return True
    else:
        return False 




""" functions for calculating other things like pressure and energy """

# def U(sigma_field):
#     """ sigma field value """

#     if sigma.b == 0.0 or sigma.c == 0.0:
#         print("Sigma meson not initialized! Call init")
    
#     term1 = (1/3) * sigma.b * Neutron.mass * (Neutron.g_sigma * sigma_field)**3
#     term2 = (1/4) * sigma.c * (Neutron.g_sigma * sigma_field)**4 
    
#     return term1 + term2



def lepton_energy_density(lepton, kf_value):
    """ leptonic energy density """
    
    lepton.kf = kf_value 
    
    m = lepton.mass 
    kf = kf_value 
    e_eff = lepton.set_ef()
    l_fac = lepton.set_log_fac() 
    
    # term1 = kf * m**2 
    # term2 = - m**4 / e_eff * l_fac 
    # term3 = 2 * kf**3 
    
    # coeff = e_eff/8/np.pi**2 
    
    term1 = kf * e_eff**3 
    term2 = kf**3 * e_eff 
    term3 = - m**4 * l_fac

    coeff = 1/8/np.pi**2
    
    return coeff * (term1 + term2 + term3)


def baryon_energy_density(baryon, kf_value, sigma_value):
    """ gets energy density """
    
    baryon.kf = kf_value 
    
    kf = kf_value 
    m_eff = baryon.set_mass_eff(sigma_value)
    e_eff = baryon.set_ef()
    l_fac = baryon.set_log_fac()
    
    # term1 = kf * m_eff**2 
    # term2 = - m_eff**4 / e_eff * l_fac 
    # term3 = 2 * kf**3 
    
    # coeff = e_eff / 8 / np.pi**2 
    term1 = kf * e_eff**3 
    term2 = kf**3 * e_eff
    term3 = - m_eff**4 * l_fac 

    coeff = 1/8/np.pi**2
    
    return coeff * (term1 + term2 + term3)


def energy_density(index, data_dictionary, eos):
    """ gets the energy density in MeV/fm3 """
    
    result = 3 * eos.xi / 24 * Neutron.g_omega**4 * data_dictionary[omega][index]**4 
    result += eos.k / 6 * (Neutron.g_sigma * data_dictionary[sigma][index])**3
    result += eos.l / 24 * (Neutron.g_sigma * data_dictionary[sigma][index])**4 
    result += 3 * eos.lambda_omega * (Neutron.g_omega * Neutron.g_rho * data_dictionary[rho][index] * data_dictionary[omega][index])**2 
    
    for particle in data_dictionary:
        if particle.type == 'Meson':
            result += 1/2 * particle.mass**2 * data_dictionary[particle][index]**2 
        elif particle.type == "Baryon":
            result += baryon_energy_density(particle, data_dictionary[particle][index], data_dictionary[sigma][index])
        elif particle.type == "Lepton":
            result += lepton_energy_density(particle, data_dictionary[particle][index])
    
    return result/hc**3 



def pressure(index, data_dictionary, meson__list, energy_density):
    """ finds the pressure via legendre transformation 
        !! assumes that the energy density has already been calculated
        units MeV/fm3 
    """
    
    for particle in data_dictionary:
        
        # update meson field values
        if particle.type == "Meson":
            particle.field_value = data_dictionary[particle][index] 
        
    for particle in data_dictionary:
        if particle.type == "Baryon":

            # be careful and set kf again though it should be set 
            particle.kf = data_dictionary[particle][index]

            # set ef again!
            particle.set_mass_eff(sigma.field_value)
            particle.set_ef()
            
            # set the chemical potential 
            particle.set_chem_pot(meson_list = meson__list)
            
            # set the number density 
            particle.num_density = particle.kf**3 / 3 / np.pi**2 
        
        elif particle.type == "Lepton":
            particle.kf = data_dictionary[particle][index] 
            particle.set_ef()
            particle.num_density = particle.kf**3 / 3 / np.pi**2
    
    result = 0.0 
    for particle in data_dictionary:
        if particle.type != "Meson":
            result += particle.chem_pot * particle.num_density
    
    return (result - energy_density)/hc**3   



"""" Functions for generating fractions 
"""


def frac(fermi, nb):
    """ calculates fraction given fermi momentum """
    return fermi**3 /3 / np.pi**2 / nb 


def frac_array(kf_array, nb_array):
    """ calculates an array of frac values from a kf array and nb array """
    frac_array = np.zeros(len(nb_array))
    
    for index, kf in enumerate(kf_array):
        frac_array[index] = frac(kf, nb_array[index])
    
    return frac_array

def fraction_fill(data_dictionary, nb_MeV):
    """ generates fraction dictionary """
    frac_dictionary = {} 
    for particle in data_dictionary:
        if particle.type != "Meson":
            frac_dictionary[particle] = frac_array(data_dictionary[particle], nb_MeV)
    
    return frac_dictionary 




""" functions for generating the data frame 
    After getting the data, we then need to create column names to then export 
    the data as a dataframe 
"""

def column_name_gen(meson_list, data_dictionary):
    """ generates a list of column names """
    
    column_names = ['nB/n0'] 
    
    for meson in meson_list:
        column_names.append(meson.name + " " + "field (MeV)")
        
    for particle in data_dictionary:
        if particle.type != "Meson":
            column_names.append(particle.name + " " + "kF (MeV)")
        
    for particle in data_dictionary:
        if particle.type != "Meson":
            column_names.append(particle.name + " " + "frac") 
    
    for particle in data_dictionary:
        if particle.type != 'Meson':
            column_names.append(particle.name + " " + "chem pot") 
        
    column_names.append("P (MeV/fm3)")
    column_names.append("E Dens (MeV/fm3)")
    
    return column_names

""" resetter """

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


""" Solver """

def full_solve(eos, baryon_list, meson_list, lepton_list, initial_guess = [7.9, 4.50, -2.25, 0.01, 210.0, 45.0, 45.0], csv_name = 'data', meth = 'hybr'):
    """ full solve for entire density range and returns a dataframe with the relevant values """
    
    n0 = 0.153 
    hc = 197.32698
    
    # create list for baryon densities 
    nb_array = np.arange(0.27, 8.01, 0.01)       # this specifically is nb/n0 
    nb_MeV = nb_array * n0 * hc**3           # these are nb values in MeV   
    
    # reset particles in case of repeated calls
    reset(baryon_list, meson_list, lepton_list)

    # initialize baryons, mesons with coupling constants 
    #sigma.set_sigma_self_coupling(eos)
    sigma.mass = eos.sigma_mass 
    for baryon in baryon_list:
        baryon.init_eos_param(eos)
    
    # create the data array 
    row_size = len(nb_array) 
    column_size = len(meson_list) + 3 * len(baryon_list + lepton_list) + 3 
    data = np.zeros((row_size, column_size), dtype = np.float64)
    data[:,0] = nb_array 
    
    
    # create the initial system which is just npe
    #current_baryons = [Neutron, Proton] 
    #current_leptons = [electron]
    current_particles = [Neutron, Proton, electron]
    
    # create lists for potential particles
    potential_baryons = potential_baryon_gen(baryon_list)
    potential_leptons = potential_lepton_gen(lepton_list)
    
    # data list
    # append data here... 
    data_dictionary = {}
    for particle in meson_list + current_particles:
        data_dictionary[particle] = np.zeros(len(nb_MeV))
    
    input_vec_dict = {}
    for index, particle in enumerate(meson_list + current_particles):
        input_vec_dict[particle] = initial_guess[index]     
    
    chem_pot_dict = {} 
    for particle in current_particles:
        chem_pot_dict[particle] = np.zeros(len(nb_MeV)) 
    
    # pressure and energy density arrays 
    pressure_array = np.zeros(len(nb_MeV))
    energy_array = np.zeros(len(nb_MeV))
    
    # the initial guess
    x_guess = initial_guess
    

    # begin the solve 
    for index, nb in enumerate(nb_MeV):
        iter = 0

        # check if new particles enter
        if potential_baryons != []:
            for pot_baryon in potential_baryons:
                Bool = baryon_threshold(pot_baryon, meson_list)
                if Bool:
                    # current_baryons.append(pot_baryon)
                    current_particles.append(pot_baryon)
                    potential_baryons.remove(pot_baryon)
                    data_dictionary[pot_baryon] = np.zeros(len(nb_array))
                    chem_pot_dict[pot_baryon] = np.zeros(len(nb_array))
                    # x_guess = np.append(x_guess, 200.0 + 50.0*iter)

                    input_vec_dict[pot_baryon] = 250.0 
                    x_guess = np.append(x_guess, 250.0)
                    iter += 1
                    #break 
                    
                    
        if potential_leptons != []:
            for pot_lepton in potential_leptons:
                Bool = lepton_threshold(pot_lepton)
                if Bool:
                    # current_leptons.append(pot_lepton)
                    current_particles.append(pot_lepton)
                    potential_leptons.remove(pot_lepton)
                    data_dictionary[pot_lepton] = np.zeros(len(nb_array))
                    chem_pot_dict[pot_lepton] = np.zeros(len(nb_array))

                    input_vec_dict[pot_lepton] = 20.0 
                    x_guess = np.append(x_guess, 20.0)
        

        #  check if particles exit 
        if index != 0:
            for particle in current_particles:
                if particle.name == 'muon':
                    if lepton_threshold(particle) != True:
                        current_particles.remove(particle)
                        # try deleting the last element and seeing if it matters lmao 
                        x_guess = np.delete(x_guess, 7)
                        del input_vec_dict[particle]
        
        # do the solve!
        # uses scipy.optimize.root function 
        solution = root(sys_eqn_new, x_guess, method = meth, args = (current_particles, meson_list, nb, eos), tol = 1e-10)
        

        # fill solution into input_vec_dict 
        for ans_index, particle in enumerate(input_vec_dict):
            input_vec_dict[particle] = solution.x[ans_index]


        # add data to dictionary old
        # for column, particle in enumerate(data_dictionary):
        #     data_dictionary[particle][index] = solution.x[column]

        # add data to dictionary new (if removing particles)
        for particle in data_dictionary:
            if particle in input_vec_dict:
                data_dictionary[particle][index] = input_vec_dict[particle]
        
        # calculate energy density and pressure 
        energy_array[index] = energy_density(index, data_dictionary, eos)
        pressure_array[index] = pressure(index, data_dictionary, meson_list, energy_array[index] * hc**3)
        
        # calculate chemical potential 
        for particle in chem_pot_dict:
            if particle in input_vec_dict: 
                chem_pot_dict[particle][index] = particle.chem_pot 

        # update the guess 
        x_guess = solution.x 
        
    # if a particle species hasn't been populated... 
    # would likely only be a baryon so we only need this here 
    if potential_baryons != []:
        for pot_baryon in potential_baryons:
            data_dictionary[pot_baryon] = np.zeros(len(nb_array))
    
    if potential_leptons != []:
        for pot_lepton in potential_leptons:
            data_dictionary[pot_lepton] = np.zeros(len(nb_array))

    # fill in the fractions 
    fraction_dictionary = fraction_fill(data_dictionary, nb_MeV)

    
    # create a data matrix appended in the correct order... 
    num_rows = len(nb_MeV)
    num_columns = len(meson_list) + 3 * len(baryon_list + lepton_list) + 3 
    data_array = np.zeros((num_rows, num_columns), dtype = 'float') 
    
    data_array[:, 0] = nb_array 
    for index, particle in enumerate(data_dictionary):
        data_array[:, index + 1] = data_dictionary[particle]
    for index, particle in enumerate(fraction_dictionary):
        data_array[:, index + len(data_dictionary) + 1] = fraction_dictionary[particle]
    for index, particle in enumerate(chem_pot_dict):
        data_array[:, index + len(data_dictionary) + len(fraction_dictionary) + 1] = chem_pot_dict[particle]
    
    data_array[:,-1] = energy_array
    data_array[:,-2] = pressure_array 
    

    # create a dataframe 
    data_frame = pd.DataFrame(data_array, columns = column_name_gen(meson_list, data_dictionary))


    # write out dataframe to csv file 
    data_frame.to_csv(csv_name + '.csv', float_format = '{:.12f}'.format) 
    
    return data_frame 
    