"""
    Class declaration 
    April 18 2022
"""
import numpy as np
from scipy.optimize import root 
from scipy.linalg import solve as linalg_solve 


# from Root_Finder.solver_helper import Sigma_min, Sigma_plus 

hc = 197.32
n0 = 0.153
pi = np.pi


""" Defining all of the classes """


class eos():
    """ equation of state class """
    def __init__(self, name, n0 = 0.153, sigma_mass = 550.0,\
        g_sigma_N = 9.57, g_omega_N = 10.61, g_rho_N = 8.20, g_phi_N = 0.0,\
        k = 0.002947, l = -0.001070, xi = 0.0, lambda_omega = 0.0,\
        g_omega_l = 7.0733333, g_omega_sig = 7.0733333, g_omega_xi = 3.5366666,\
        g_rho_l = 0.0, g_rho_sig = 21.22, g_rho_xi = 10.61,\
        g_phi_l = -5.0016, g_phi_sig = -5.0016, g_phi_xi = -10.0032,\
        g_sigma_l = 0.0, g_sigma_xi = 0.0, g_sigma_sig = 0.0,\
        g_xi_l = 3.73, g_xi_sig = 0.0, g_xi_xi = 9.67):
        
        self.name = name 
        self.sigma_mass = sigma_mass 
        self.n0 = n0 
        self.g_sigma_N = g_sigma_N
        self.g_omega_N = g_omega_N
        self.g_rho_N = g_rho_N
        self.g_phi_N = g_phi_N
        self.g_xi_N = 0.0
        
        self.g_sigma_l = g_sigma_l
        self.g_omega_l = g_omega_l
        self.g_rho_l = g_rho_l
        self.g_phi_l = g_phi_l
        self.g_xi_l = g_xi_l
        
        self.g_sigma_sig = g_sigma_sig 
        self.g_omega_sig = g_omega_sig 
        self.g_phi_sig = g_phi_sig 
        self.g_rho_sig = g_rho_sig 
        self.g_xi_sig = g_xi_sig

        self.g_sigma_xi = g_sigma_xi 
        self.g_omega_xi = g_omega_xi 
        self.g_rho_xi = g_rho_xi 
        self.g_phi_xi = g_phi_xi 
        self.g_xi_xi = g_xi_xi 


        self.k = k
        self.l = l
        self.xi = xi 
        self.lambda_omega = lambda_omega 





class particle():
    """ base particle class """
    def __init__(self, name = 'NaN', mass = 0.0, charge = 0.0):
        
        self.name = name 
        self.mass = mass 
        self.charge = charge
        self.spin = 1/2
        
        # "base attributes" 
        self.num_density = 0.0 
        self.frac = 0.0
        self.kf = 0.0
        
        # "derived" attributes
        self.ef = 0.0 
        self.chem_pot = 0.0 
        
        self.d_kf = 0.0
        self.d_ef = 0.0 
        self.d_chem_pot = 0.0 
        self.d_chem_pot_tilde = 0.0 
        
        # coupling constants 
        self.g_sigma = 0.0 
        self.g_omega = 0.0 
        self.g_rho = 0.0 
        self.g_phi = 0.0 
        
    # methods for updating base attributes from dataframe
    # given a baryon density (index) 
    # assuming dataframe is named following this project's convention 
    
    def init_check(self):
        """ checks if particle has been instantized or not """

        if self.name == 'NaN':
            print("Name is 'NaN' check to see if baryon object has been instantized")
        if self.mass == 0.0:
            print("Mass is 0.0, check to see if baryon object has been instantized")
            
    
    def set_kf(self, index, data_frame):
        """  reads in the kf value from data frame """

        self.init_check(self)

        for column_name in data_frame.columns:
            if column_name.split()[0] == self.name and column_name.split()[1].lower() == 'kf':
                self.kf = data_frame.loc[index, column_name]
        
        return self.kf 

    
    def set_num_density(self, index, data_frame):
        """ reads in kf value and sets the number density """ 

        kf_val = self.set_kf(self, index, data_frame) 

        self.num_density = kf_val**3 / 3 / pi**2 

        return self.num_density 
    

    def set_frac(self, index, data_frame):
        """ reads in the fraction """
        
        num_dens = self.set_num_density(index, data_frame) 
        nB = data_frame.loc[index, 'nB/n0']
        self.frac = num_dens/nB

        return self.frac 



    def base_set(self, index, data_frame): 
        """ calls the three functions defined above """

        self.set_kf(index, data_frame) 
        self.set_num_density(index, data_frame) 
        self.set_frac(index, data_frame) 











class baryon(particle):
    """ baryon particle class """
    
    def __init__(self, name = 'NaN', mass = 0.0, charge = 0.0, kind = 'NaN', isospin = 0.0):
        super().__init__(name, mass, charge)
        self.kind = kind 
        self.type = 'Baryon'
        self.isospin = isospin 
    
        # other attributes
        self.mass_eff = 0.0 
        self.scalar_density = 0.0 
        self.log_fac = 0.0

        # derivative attributes for chemical potential partial derivatives
        self.d_mass_eff =  0.0 
        self.d_ef = 0.0 
        self.d_kf = 0.0 
        self.d_log_fac = 0.0 
        self.d_sigma = 0.0 

        # other things... 
        self.partial_mu_m = 0.0 

    
    
    """ baryon class methods 
        most of these are to calculate and set relevant quantities 
    """

    def init_eos_param(self, eos_object):
        """ loads in eos params """
        
        if self.kind == 'Nucleon':
            self.g_sigma = eos_object.g_sigma_N
            self.g_omega = eos_object.g_omega_N 
            self.g_rho = eos_object.g_rho_N 
            self.g_phi = eos_object.g_phi_N 
            self.g_xi = 0.0 

        elif self.kind == 'Hyperon':
            if self.name.split('_')[0] == 'Lambda':
                self.g_sigma = eos_object.g_sigma_l
                self.g_omega = eos_object.g_omega_l
                self.g_rho = eos_object.g_rho_l
                self.g_phi = eos_object.g_phi_l
                self.g_xi = eos_object.g_xi_l
            
            elif self.name.split('_')[0] == 'Sigma':
                self.g_sigma = eos_object.g_sigma_sig
                self.g_omega = eos_object.g_omega_sig
                self.g_rho = eos_object.g_rho_sig
                self.g_phi = eos_object.g_phi_sig
                self.g_xi = eos_object.g_xi_sig

            elif self.name.split('_')[0] == 'Xi':
                self.g_sigma = eos_object.g_sigma_xi
                self.g_omega = eos_object.g_omega_xi
                self.g_rho = eos_object.g_rho_xi
                self.g_phi = eos_object.g_phi_xi
                self.g_xi = eos_object.g_xi_xi 


        

    

    def set_mass_eff(self, sigma_field):
        """ sets effective mass """
        
        self.mass_eff = self.mass - self.g_sigma * sigma_field
        return self.mass_eff 
        
        
    def set_ef(self, kf = 'NaN', sigma_field = 0.0):
        """ sets ef for baryon """
        
        if kf == 'NaN':
            # used stored value for kf
            self.ef = np.sqrt(self.kf**2 + self.mass_eff**2)
            return self.ef
        else:
            # used specified kf 
            self.ef = np.sqrt(kf**2 + self.set_mass_eff(sigma_field)**2)
            return self.ef
    
    

    def set_log_fac(self, kf = 'NaN', sigma_field = 'NaN'):
        """ sets log factor """

        if kf == 'NaN' and sigma_field == 'NaN':
            # assume values are already stored and call from baryon
            numerator = (self.kf + self.ef)**2
            denominator = self.mass_eff**2
            self.log_fac = np.log(np.sqrt(numerator/denominator))
            return np.log(np.sqrt(numerator/denominator))
        
        elif kf != 'NaN' and sigma_field != 'NaN':
            # used specified values
            
            m_eff = self.set_mass_eff(sigma_field)
            e_f = self.set_ef(kf, sigma_field)
            
            numerator = (kf + e_f)**2
            denominator = m_eff**2
            self.log_fac = np.log(np.sqrt(numerator/denominator))
            return np.log(np.sqrt(numerator/denominator))
        
        else:
            print("self.set_log_fac Error: not all of kf or sigma_field are 'NaN's or full numeric values.")
    

    
    def set_scalar_density(self, kf = 'NaN', sigma_field = 'NaN'):
        """ sets and returns scalar density for given baryon """
        
        if kf == 'NaN' and sigma_field == 'NaN':
            # use stored values
            
            first_term = self.kf * self.ef 
            second_term = -self.mass_eff**2 * self.log_fac
            coeff = self.mass_eff/(2*np.pi**2)
            
            return coeff * (first_term + second_term)
        
        elif kf != 'NaN' and sigma_field != 'NaN':
            # use inputted values
            
            m_eff = self.set_mass_eff(sigma_field)
            e_f = self.set_ef(kf, sigma_field)
            l_fac = self.set_log_fac(kf, sigma_field)
            
            term1 = kf * e_f
            term2 = - m_eff**2 * l_fac
            
            coeff = m_eff/(2*np.pi**2)
            
            return coeff * (term1 + term2)

        else:
            print("self.scalar_density Error: not all of kf or sigma_field are 'NaN' or numeric values.")
            

    
    def set_chem_pot(self, kf = 'NaN', meson_list = []):
        """ sets chemical potential assuming meson objects have already been updated with new field values """
        
        if kf == 'NaN':
            # used stored values
            result = self.ef 

            # for meson in meson_list:
            #     if meson.name == 'omega':
            #         result += self.g_omega * meson.field_value
            #     elif meson.name == 'rho':
            #         result += self.isospin * self.g_rho * meson.field_value
            #     elif meson.name == 'phi':
            #         result += self.g_phi * meson.field_value 

            result += self.g_omega * omega.field_value 
            result += self.isospin * self.g_rho * rho.field_value 
            result += self.g_phi * phi.field_value 
            
            self.chem_pot = result 
            return self.chem_pot
            

        elif kf != 'NaN':
            # used inputted values 

            result = 0.0 

            for meson in meson_list:
                if meson.name == 'sigma':
                    sigma_field = meson.field_value 
                    e_f = self.set_ef(kf, sigma_field, xi.field_value)
                    result += e_f                
                if meson.name == 'omega':
                    result += self.g_omega * meson.field_value
                elif meson.name == 'rho':
                    result += self.isospin * self.g_omega * meson.field_value
                elif meson.name == 'phi':
                    result += self.g_phi * meson.field_value 
            
            self.chem_pot = result 
            return result
    


    def solver_chem_pot(self, kf, meson_var_list):
        """ calculates baryon chem pot in context of non linear solver
            similar but slightly different than the set chem pot function 
        """

        for index, particle in enumerate(meson_var_list[1]):
            if particle.name == 'sigma':
                sigma_field_val = meson_var_list[0][index]
        
        e_f = self.set_ef(kf, sigma_field_val) 

        # calculate remaining contribution 
        result = e_f 

        for index, particle in enumerate(meson_var_list[1]):
            if particle.name == 'omega':
                result += self.g_omega * meson_var_list[0][index] 
            elif particle.name == 'rho':
                result += self.g_rho * self.isospin * meson_var_list[0][index] 
            elif particle.name == 'phi':
                result += self.g_phi * meson_var_list[0][index]
        
        return result 
    


    def baryon_chem_pot(self, baryon_list, meson_list, eos):
        """ calculates dd-rmf baryon chemical potential 
            assuming baryon and meson objects have already been updated 
            There are three terms: E_eff, mu_meson, and Sigma r
            - E_eff is assumed already stored and updated in baryon object
        """
        
        # E_eff
        result = self.ef
        
        # mu_meson
        for meson in meson_list:
            if meson.name == 'omega':
                result += self.g_omega * meson.field_value
            elif meson.name == 'rho':
                result += self.isospin * self.g_rho * meson.field_value
            elif meson.name == 'phi':
                result += self.g_phi * meson.field_value 
        
        
        self.chem_pot = result
        return self.chem_pot 
    




    """ chemical potential partial derivative methods """

    def set_d_mass_eff(self, d_sigma = 'NaN'):
        """ 
            d_mass_eff in NL3 model, only sigma meson 
            updated 06/19/2022
        """

        if self.g_sigma == 0.0:
            print(" g sigma is zero ")

        if d_sigma == 'NaN':
            self.d_mass_eff = self.g_sigma * sigma.partial_nb 
            return self.d_mass_eff 

        else: 
            self.d_mass_eff = - self.g_sigma * d_sigma
            return self.d_mass_eff 


    def set_d_kf(self, frac = 'NaN', kf = 'NaN'):
        """ set dkf/dnb """

        if frac == "NaN" and kf == "NaN":
            # used stored values 
            self.d_kf = np.pi**2 * self.frac / self.kf**2
            return self.d_kf
        elif frac != 'NaN' and kf != 'NaN':
            self.d_kf = np.pi**2 * frac / kf**2
            return self.d_kf 
        else:
            print('Not all values initialized or specified ')
    


    def set_d_log_fac(self, d_sigma = 'NaN', d_ef = 'NaN', frac = 'NaN', kf = 'NaN', ef = 'NaN', sigma_field = 'NaN'):
        """ set d_log_fac / dnb """

        if d_sigma == 'NaN' and d_ef == 'NaN' and frac == 'NaN' and kf == 'NaN' and ef == 'NaN':
            # used stored values 
            term1 = (self.d_kf + self.d_ef) / (self.kf + self.ef)
            term2 = - self.d_mass_eff / self.mass_eff * self.mass_eff 
            self.d_log_fac = term1 + term2 
            return term1 + term2 

        else:
            # d_kf = self.d_kf(frac, kf) 
            # m_eff = self.set_mass_eff(sigma_field)
            # term1 = (d_kf + d_ef)/(kf + ef)
            # term2 = - 1 /  m_eff * self.set_d_mass_eff(d_sigma)

            term1 = (self.d_kf + d_ef) / (self.kf + self.ef) 
            term2 = - self.set_d_mass_eff(d_sigma)/self.mass_eff 

            self.d_log_fac = term1 + term2 

            return term1 + term2
    

    def set_partial_chem_pot_two(self):
        """ sets partial derivative of chemical potential 
            assuming everything has already been loaded in 

            06/18/2022

            assumes that partial derivatives of EFi already calculated and we already have 
            the meson field partial derivatives and found the derivative of the scalar densities 
        """
        
        result = self.d_ef 
        
        # next finding partial sigma 0 contribution 
        self.set_partial_mu_m()
        result += self.partial_mu_m
        
        
        self.d_chem_pot = result
        return self.d_chem_pot 


    def set_partial_mu_m(self):
        """ calculates mesonic contribution to partial derivative of chemical potential 
            for given baryon object 
        """
        
        result = self.g_omega * omega.partial_nb 
        result += self.g_rho * self.isospin * rho.partial_nb 
        result += self.g_phi * phi.partial_nb 
                
        self.partial_mu_m = result 
        
        return self.partial_mu_m





class lepton(particle):
    """ lepton particle class """
    
    def __init__(self, name = 'NaN', mass = 0.0, charge = 0.0):
        super().__init__(name, mass, charge)
        self.type = 'Lepton'
    
    
    # class methods 
    
    def set_ef(self, kf = 'NaN'):
        """ sets and returns ef for Lepton """
        
        if kf == 'NaN':
            self.ef = np.sqrt(self.kf**2 + self.mass**2)
            self.chem_pot = self.ef 
            return self.ef
        else:
            self.ef = np.sqrt(kf**2 + self.mass**2)
            self.chem_pot = self.ef 
            return self.ef 
    

    def solver_chem_pot(self, kf):
        """ gets the lepton chemical potential in context of non linear solver """
        return np.sqrt(kf**2 + self.mass**2)
    

    def set_log_fac(self, kf = 'NaN'):
        """ leptonic log factor """

        if kf == 'NaN':
            numerator = (self.kf + self.ef)**2
            denominator = self.mass**2 
            self.log_fac = np.log(np.sqrt(numerator/denominator))
            return self.log_fac 

        elif kf != 'NaN':
            e_f = self.set_ef(kf)
            
            numerator = (kf + e_f)**2
            denominator = self.mass**2
            self.log_fac = np.log(np.sqrt(numerator/denominator))
            return self.log_fac
    

    def set_partial_chem_pot(self):
        """ set partial derivative of chemical potential """
        numerator = np.pi**2 * self.frac 
        denominator = self.kf * self.ef 

        self.d_chem_pot = numerator/denominator 
        return self.d_chem_pot 


    

     
        
    






class meson(particle):
    """ meson particle class """

    type = 'Meson'
    field_value = 0.0 
    b = 0.0 
    c = 0.0 
    partial_nb = 0.0

    def set_sigma_self_coupling(self, eos):
        """ sets b and c """
        self.b = eos.b 
        self.c = eos.c 
    

    def set_partial_nb(self, baryon_list, eos):
        """ returns partial derivative of vector meson fields (ie, not the sigma or xi meson) wrt to nb """
        
        if self.name == 'sigma' or self.name == 'xi':
            print('Error: sigma eom is not generated by this function') 
            return 1 
        
        result = 0.0 
        
        if self.name == 'phi':
            for baryon in baryon_list:
                result += baryon.g_phi * baryon.frac 
            self.partial_nb = 1/self.mass**2 * result 
            return self.partial_nb 
        
        else: 
            self.partial_omega_rho(baryon_list, eos)
    


    def partial_omega_rho(self, baryon_list, eos):
        """ solves for partial omega partial rho in this 
            NL3 model with omega rho coupling 

            in this case the solution or x array is [partial omega/partial nb, partial rho/partial nb]
            in that order 
        """
        g_rho = Neutron.g_rho 
        g_omega = Neutron.g_omega 

        # first here we generate the matrices to solve the linear system of equations
        # Ax = b for the partial derivatives as shown in the thesis notes 

        A = np.zeros((2,2), dtype = np.float64) 
        b = np.zeros(2, dtype = np.float64)

        A[0][0] = omega.mass**2 + eos.xi / 2 * g_omega**2 * omega.field_value**2 + 2 * eos.lambda_omega * g_rho**2 * rho.field_value**2 
        A[0][1] = 4 * eos.lambda_omega * g_rho**2 * g_omega**2 * rho.field_value * omega.field_value 
        A[1][0] = A[0][1] 
        A[1][1] = rho.mass**2 + 2 * eos.lambda_omega * g_rho**2 * g_omega**2 * omega.field_value**2 

        for baryon in baryon_list:
            b[0] += baryon.g_omega * baryon.frac 
            b[1] += baryon.g_rho * baryon.isospin * baryon.frac 
        

        # solving the system of linear equations 
        solution = linalg_solve(A, b, assume_a = 'sym')

        omega.partial_nb = solution[0] 
        rho.partial_nb = solution[1]

        return solution 
    

    def set_partial_nb(self, sys_eqn, input_vec, baryon_list, eos):
        """ returns partial derivative of vector meson fields with respect to nB """

        solution = root(sys_eqn, input_vec, args = (baryon_list, eos), method = 'hybr', tol = 1e-9)

        if self.name == 'omega':
            self.partial_nb = solution.x[0] 
        elif self.name == 'rho':
            self.partial_nb = solution.x[1]   


        





                
        
        
        











""" Making a list of all the relevant particles and preloading them """ 


# baryons 
Neutron = baryon(name = 'Neutron', mass = 939.00, charge = 0.0, kind = 'Nucleon', isospin = -1/2)
Proton = baryon(name = 'Proton', mass = 939.00, charge = 1.0, kind = 'Nucleon', isospin = 1/2)
Lambda =  baryon(name = 'Lambda', mass = 1116.00, charge = 0.0, kind = 'Hyperon', isospin = 0.0)
Sigma_neu = baryon(name = 'Sigma_neu', mass = 1190.00, charge = 0.0, kind = 'Hyperon', isospin = 0.0)
Sigma_min = baryon(name = 'Sigma_min', mass = 1190.00, charge = - 1.0, kind = 'Hyperon', isospin = -1.0)
Sigma_plus = baryon(name = 'Sigma_plus', mass = 1190.00, charge = 1.0, kind = 'Hyperon', isospin = 1.0)
Xi_min = baryon(name = 'Xi_min', mass = 1315.00, charge = -1.0, kind = 'Hyperon', isospin = -1/2)
Xi_neu = baryon(name = 'Xi_neu', mass = 1315.00, charge = 0.0, kind = 'Hyperon', isospin = 1/2)



# leptons 
electron = lepton(name = 'electron', mass = 0.510, charge = -1.0)
muon = lepton(name = 'muon', mass = 105.65, charge = -1.0)


# mesons 
sigma = meson(name = 'sigma', mass = 492.730)
omega = meson(name = 'omega', mass = 782.5)
rho = meson(name = 'rho', mass = 763.0)
phi = meson(name = 'phi', mass = 1020.0)
xi = meson(name = 'xi', mass = 980.0)






""" Miscellaneous equations """

def d_sigma_eom(d_sigma, baryon_ef_list, eos):
    """ returns d sigma / dnB equation of motion 
        for dd-rmf model 
    """
    
    lhs, rhs = 0.0, 0.0 
    
    #lhs = d_sigma * (sigma.mass**2 + 2 * Neutron.g_sigma**3 * Neutron.mass * sigma.field_value + 3 * #Neutron.g_sigma**4 * sigma.field_value**2)

    lhs = d_sigma * sigma.mass**2
    lhs += d_sigma * eos.k * Neutron.g_sigma**3 * sigma.field_value
    lhs += d_sigma * eos.l / 2 * Neutron.g_sigma**4 * sigma.field_value**2 
    
    for index, baryon in enumerate(baryon_ef_list[1]):
        term1_coeff = baryon.g_sigma / (2 * np.pi**2) * baryon.set_d_mass_eff(d_sigma)
        #term1_a = baryon.scalar_density WRONG 
        term1_a = baryon.kf * baryon.ef - baryon.mass_eff**2 * baryon.log_fac 
        term1 = term1_coeff * term1_a
        
        term2_coeff = baryon.g_sigma * baryon.mass_eff / (2 * np.pi**2)
        term2_a = baryon.set_d_kf() * baryon.ef 
        term2_b = baryon.kf * baryon_ef_list[0][index]
        term2_c = -2 * baryon.mass_eff * baryon.log_fac * baryon.set_d_mass_eff(d_sigma)
        term2_d = - baryon.mass_eff**2 * baryon.set_d_log_fac(d_sigma, baryon_ef_list[0][index])
        
        term2 = term2_coeff * (term2_a + term2_b + term2_c + term2_d)
        
        rhs += term1 + term2
    
    return lhs - rhs


# def d_xi_eom(d_sigma, d_xi, baryon_ef_list):
#     """ returns d xi/ dnB equation of motion for NLW model """

#     lhs, rhs = 0.0, 0.0

#     lhs = d_xi * xi.mass**2 

#     for index, baryon in enumerate(baryon_ef_list[1]):
#         term1_coeff = baryon.g_xi / (np.pi**2) * baryon.set_d_mass_eff(d_sigma, d_xi)
#         term1_a = baryon.scalar_density 
#         term1 = term1_coeff * term1_a
        
#         term2_coeff = baryon.g_xi * baryon.mass_eff / (np.pi**2)
#         term2_a = baryon.set_d_kf() * baryon.ef 
#         term2_b = baryon.kf * baryon_ef_list[0][index]
#         term2_c = -2 * baryon.mass_eff * baryon.log_fac * baryon.set_d_mass_eff(d_sigma, d_xi)
#         term2_d = - baryon.mass_eff**2 * baryon.set_d_log_fac(d_sigma, d_xi, baryon_ef_list[0][index])
        
#         term2 = term2_coeff * (term2_a + term2_b + term2_c + term2_d)
        
#         rhs += term1 + term2
    
#     return lhs - rhs 


def baryon_d_ef(d_sigma, baryon_ef_list):
    """ returns the d_ef equations """
    
    baryon_d_ef_array = []
    
    for index, baryon in enumerate(baryon_ef_list[1]):
        lhs, rhs = 0.0, 0.0 
        
        lhs = baryon_ef_list[0][index]
        
        term1 = baryon.kf / baryon.ef * baryon.set_d_kf()
        term2 = baryon.mass_eff / baryon.ef * baryon.set_d_mass_eff(d_sigma)
        
        rhs = term1 + term2 
        
        baryon_d_ef_array.append(lhs - rhs)
        
    return baryon_d_ef_array


def equations_gen(input_vec, baryon_list, eos):

    """ generate the system of equations """

    # take input vec and link to objects 
    # first row is the input variable guess
    # second row is the corresponding baryon object 

    d_sigma = input_vec[0]
   #d_xi = input_vec[1] 

    # baryon EF variable list 
    baryon_ef_list = [[],[]]
    for i in range(len(baryon_list)):
        baryon_ef_list[0].append(input_vec[i + 1])
        baryon_ef_list[1].append(baryon_list[i])

    sigma_eqn = d_sigma_eom(d_sigma, baryon_ef_list, eos)
    #xi_eqn = d_xi_eom(d_sigma, d_xi, baryon_ef_list)

    baryon_ef = baryon_d_ef(d_sigma, baryon_ef_list)

    system = [sigma_eqn] + baryon_ef 

    return system 


def step_solve(baryon_list, eos):
    """ 
        solve for d_sigma and d_ef 
        using system of nonlinear equation 
    """

    initial_guess = [sigma.partial_nb] 
    
    for baryon in baryon_list:
        initial_guess.append(baryon.d_ef) 
    
    solutions = root(equations_gen, initial_guess, args = (baryon_list, eos), tol = 1e-10)

    sigma.partial_nb = solutions.x[0]

    for index, baryon in enumerate(baryon_list):
        baryon.d_ef = solutions.x[index + 1]
