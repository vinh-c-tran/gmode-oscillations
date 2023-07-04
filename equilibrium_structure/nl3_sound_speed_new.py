"""
    Trying to replicate the May sound speed stuff here 
"""

from nlw_classes import * 

"""
    Base derivative expressions 
"""

def d_m_eff(d_sigma, baryon):
    return - baryon.g_sigma * d_sigma 

def d_kf_dnb(baryon):
    return np.pi **2 * baryon.frac / baryon.kf**2 

def d_log_fac_dnb(d_sigma, d_ef, baryon):
    term1 = (baryon.d_kf + d_ef) / (baryon.kf + baryon.ef)
    term2 = - d_m_eff(d_sigma, baryon) / baryon.mass_eff 

    return term1 + term2 

def log_factor(baryon):
    numerator = baryon.kf + baryon.ef 
    denominator = baryon.mass_eff 

    return np.log(np.sqrt(numerator**2) / np.sqrt(denominator**2))


""" First equation """

def EF_deriv(d_sigma, baryon):
    """ gets expression for equation 1: dE_F/dnB """

    coeff = 1 / baryon.ef 
    term1 = baryon.kf * d_kf_dnb(baryon)
    term2 = baryon.mass_eff * d_m_eff(d_sigma, baryon)

    return coeff * (term1 + term2)

def equation_one(d_sigma, d_ef, baryon):
    return d_ef - EF_deriv(d_sigma, baryon)


""" Second Equation """

def U_second_deriv():
    """ d^2U / dsigma^2 """
    term1 = 2 * sigma.b * Neutron.mass * Neutron.g_sigma**3 * sigma.field_value 
    term2 = 3 * sigma.c * Neutron.g_sigma**4 * sigma.field_value 

    return term1 + term2 

def sigma_eom_lhs(d_sigma):

    return d_sigma * (sigma.num_mass**2 + U_second_deriv())


# rhs equations 


def first_term(d_sigma, baryon):

    coeff = 1/(2 * np.pi**2) * d_m_eff(d_sigma, baryon)  

    term1 = baryon.kf * baryon.ef 
    term2 = baryon.mass_eff**2 * baryon.log_fac 

    return coeff * (term1 + term2) 

 


def second_term(d_sigma, d_ef, baryon):

    coeff = 1/(2 * np.pi**2) * baryon.mass_eff 

    term1 = baryon.d_kf * baryon.ef 
    term2 = baryon.kf * d_ef 
    term3 = - 2 * baryon.mass_eff * d_m_eff(d_sigma, baryon) * log_factor(baryon)
    term4 = - baryon.mass_eff**2 * d_log_fac_dnb(d_sigma, d_ef, baryon)

    return coeff * (term1 + term2 + term3 + term4)



def rhs(d_sigma, d_ef, baryon_list):
    rhs = 0.0 

    for baryon in baryon_list: 
        rhs += baryon.g_sigma * (first_term(d_sigma, baryon) + second_term(d_sigma, d_ef, baryon))