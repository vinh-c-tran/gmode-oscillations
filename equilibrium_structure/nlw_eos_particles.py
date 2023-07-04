"""
    NLW equation of states and particles 

    This file declares instances of the various equations of state which fills in the
    relevant parameters for that specific model as well as the particles present 
"""

from nlw_classes import eos, particle, baryon, lepton, meson 


""" Initializing the various Equations of State / Models 

    Coupling constants are taken from literature. Big Apple was from Fattoyev et al 2020. 

    Nucleonic coupling constants are given, hyperonic ones are calculated using SU(6) symmetry 
"""

big_apple = eos(name = 'Big Apple', n0 = 0.155, sigma_mass = 492.730,\
        g_sigma_N = 9.669922, g_omega_N = 12.31160, g_rho_N = 14.161787, g_phi_N = 0.0,\
        k = 5.20326, l = -0.021739, xi = 0.00070, lambda_omega = 0.047471,\
        g_omega_l = 8.207733, g_omega_sig = 8.2077333, g_omega_xi = 4.1038666,\
        g_rho_l = 0.0, g_rho_sig = 28.323574, g_rho_xi = 14.161787,\
        g_phi_l = 5.80374389, g_phi_sig = 5.80374389, g_phi_xi = 11.60748779,\
        g_sigma_l = 5.765626, g_sigma_sig = 4.1313977, g_sigma_xi = 2.8555759,\
        g_xi_l = 0.0, g_xi_sig = 0.0, g_xi_xi = 0.0                 
        )


big_apple_new = eos(name = 'Big Apple', n0 = 0.155, sigma_mass = 492.730,\
        g_sigma_N = 9.669922, g_omega_N = 12.31160, g_rho_N = 14.161787, g_phi_N = 0.0,\
        k = 5.20326, l = -0.021739, xi = 0.00070, lambda_omega = 0.0,\
        g_omega_l = 8.207733, g_omega_sig = 8.2077333, g_omega_xi = 4.1038666,\
        g_rho_l = 0.0, g_rho_sig = 28.323574, g_rho_xi = 14.161787,\
        g_phi_l = 5.80374389, g_phi_sig = 5.80374389, g_phi_xi = 11.60748779,\
        g_sigma_l = 5.765626, g_sigma_sig = 4.1313977, g_sigma_xi = 2.8555759,\
        g_xi_l = 0.0, g_xi_sig = 0.0, g_xi_xi = 0.0                 
        )

iufsu = eos(name = 'IUFSU', n0 = 0.155, sigma_mass = 491.5,\
        g_sigma_N = 9.97129, g_omega_N = 13.03207, g_rho_N = 13.58999, g_phi_N = 0.0,\
        k = 3.38081, l = 0.000296, xi = 0.0300, lambda_omega = 0.046000,\
        g_omega_l = 8.68804666, g_omega_sig = 8.68804666, g_omega_xi = 4.34402333,\
        g_rho_l = 0.0, g_rho_sig = 27.17998, g_rho_xi = 13.58999,\
        g_phi_l = 6.1433767, g_phi_sig = 6.1433767, g_phi_xi = 12.2867534,\
        g_sigma_l = 6.0019784, g_sigma_sig = 4.3290279, g_sigma_xi = 2.9731067,\
        g_xi_l = 0.0, g_xi_sig = 0.0, g_xi_xi = 0.0  
        )

fsugarnet = eos(name = 'FSUGarnet', n0 = 0.153, sigma_mass = 496.939,\
            g_sigma_N = 10.50472, g_omega_N = 13.70017, g_rho_N = 13.88983, g_phi_N = 0.0,\
            k = 3.26018, l = -0.003551, xi = 0.02350, lambda_omega = 0.043377,\
            g_omega_l = 9.1334466, g_omega_sig = 9.1334466, g_omega_xi = 4.5667233,\
            g_rho_l = 0.0, g_rho_sig = 27.77966, g_rho_xi = 13.88983,\
            g_phi_l = 6.45832207, g_phi_sig = 6.45832207, g_phi_xi = 12.916644,\
            g_sigma_l = 6.366259, g_sigma_sig = 4.74118, g_sigma_xi = 3.156045,\
            g_xi_l = 0.0, g_xi_sig = 0.0, g_xi_xi = 0.0              
            )



""" Initializing the Particles """

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