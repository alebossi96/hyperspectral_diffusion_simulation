import numpy as np
import matplotlib.pyplot as plt
import core
exp = np.exp
v = 29.9792/1.4 #cm/ns


A = core.A_parameter(1, 1.4)


def z_mean(t, mu_s):
    res = np.zeros((len(t),))
    D = 1/(3*mu_s)
    z_e = 2*A*D
    m_tot = 200
    _4Dvt = 4*D*v*t
    z_s = 1/mu_s
    for m in range(-m_tot, m_tot+1):
        if m == 0:
            continue
        
        M = -2*m*z_e - 2*m*2*D - z_s
        N = -2*m*2*D - (2*m-2)*z_e +z_s
        #TODO trova condizione per stop
        res += (D*t*v/m)*(exp(-M**2/_4Dvt)- exp(-N**2/_4Dvt))
    res = res/(- z_s*exp(-(z_s)**2/_4Dvt)- (2*z_e+z_s)*exp(-(2*z_e+z_s)**2/_4Dvt))
    return res/2
        
def fit_fract(function_values, x, t_values):
    """
    Find the time 't' when the height of the curve defined by diff_eq_slab is a fraction 'x' 
    of the maximum of the curve before the peak if x is negative, or after the peak if x is positive.
    
    :param diff_eq_slab: Function defining the curve
    :param x: Fraction of the maximum value to find
    :param t_values: time
    :return: Time 't' where the height of the curve is x times the maximum of the curve before/after the peak
    """
    # Find the peak
    peak_index = np.argmax(function_values)
    peak_value = function_values[peak_index]

    # Calculate the height to find
    height = abs(x) * peak_value

    # Find the index where the function reaches the desired height
    if x > 0:
        # After the peak
        target_indices = np.where(function_values[peak_index:] <= height)[0]
        if target_indices.size > 0:
            target_index = target_indices[0] + peak_index
        else:
            return None
    else:
        # Before the peak
        target_indices = np.where(function_values[:peak_index] <= height)[0]
        if target_indices.size > 0:
            target_index = target_indices[-1]
        else:
            return None

    # The time 't' where the height of the curve is x times the maximum of the curve before/after the peak
    return (t_values[target_index], target_index) 
    



def diff_eq_slab(t, rho, mu_s, mu_a = 0):
    """
    Assumo illuminazione su tutto piano
    CASO IPERSPETTRALE
    """
    D = 1/(3*mu_s)
    _4Dvt = 4*D*v*t
    z_e = 2*A*D
    z_s = 1/mu_s
    return (exp(-z_s**2/_4Dvt)-exp(-(4*z_e**2+z_s**2)/_4Dvt))*exp(-(rho**2+z_s**2)/_4Dvt -mu_a*v*t)/t**(1./2)

def z_average(mu_s, mu_a):
    t = np.arange(3, 4096)*3.05*1e-3
    z_mean_mu_s = z_mean(t, mu_s)
    z_mean_mu_s[:5] = 0 #correzione per early photons ad minchiam
    phi = diff_eq_slab(t, 0, mu_s = mu_s, mu_a = mu_a)
    return np.sum(z_mean_mu_s*phi/np.sum(phi))
    
for j in range(3,10, 2):
    print(j)
    z_averages = []
    for i in range(0,200, 10):
        z_avg = z_average(mu_s= j, mu_a= i/100)
        z_averages.append(z_avg)
    plt.plot(np.arange(0,200, 10)/100,z_averages)
plt.show()

