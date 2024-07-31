import numpy as np
import core 

def z_m_plus_minus(m, s, z_e, z_s):
    z_m_plus = 2 * m * (s + 2 * z_e) + z_s
    z_m_minus = 2 * m * (s + 2 * z_e) - 2 * z_e - z_s
    return z_m_plus, z_m_minus

def phi(rho, z, mu_a, mu_sp, s, ni, ne = 1,  m_limit=100):
    D = 1/(3*mu_sp)
    phi_value = 0
    z_s = 1/mu_sp
    A = core.A_parameter(ne, ni)
    z_e = 2*A*D
    mu_eff = np.sqrt(3*mu_a*mu_sp)
    for m in range(-m_limit, m_limit + 1):
        z_m_plus, z_m_minus = z_m_plus_minus(m, s, z_e, z_s)
        
        sqrt_term_plus = np.sqrt(rho**2 + (z - z_m_plus)**2)
        sqrt_term_minus = np.sqrt(rho**2 + (z - z_m_minus)**2)
        
        term_plus = np.exp(-mu_eff * sqrt_term_plus) / sqrt_term_plus
        term_minus = np.exp(-mu_eff * sqrt_term_minus) / sqrt_term_minus
        if mu_eff * sqrt_term_plus > 100 or mu_eff * sqrt_term_minus >100:
            break
        phi_value += term_plus - term_minus
    
    return phi_value / (4 * np.pi * D)



result = phi(rho, z, mu_a, mu_sp, s, ni, ne = 1,  m_limit=100)
def thin_layer_absorber(th, mu_a)
    return np.exp(-mu_a*th)
final_fluence = result*thin_layer_absorber(th, mu_a)
