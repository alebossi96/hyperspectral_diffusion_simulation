import numpy as np
def A_parameter(ne, ni):
    """
    Computes parameter A; equation (27)

    Args:
    ne: external medium
    ni: diffusing medium

    Parameter A is computed using the solution presented in:
    Contini D, Martelli F, Zaccanti G.
    Photon migration through a turbid slab described by a model
    based diffusion approximation. I. Theory
    Applied Optics Vol 36, No 19, 1997, pp 4587-4599
    """

    # page 4590
    n = ni / ne

    if n > 1:
        # equations (30)
        t1 = 4 * (-1 - n**2 + 6*n**3 - 10*n**4 - 3*n**5 + 2*n**6 + 6*n**7 - 3*n**8 - (6*n**2 + 9*n**6) * np.sqrt(n**2 - 1)) / \
            (3 * n * (n**2 - 1)**2 * (n**2 + 1)**3)
        t2 = (-8 + 28*n**2 + 35*n**3 - 140*n**4 + 98*n**5 - 13*n**7 + 13*n * (n**2 - 1)**3 * np.sqrt(1 - (1 / n**2))) / \
            (105 * n**3 * (n**2 - 1)**2)
        t3 = 2 * n**3 * (3 + 2 * n**4) * np.log(
            ((n - np.sqrt(1 + n**2)) * (2 + n**2 + 2 * np.sqrt(1 + n**2)) * (n**2 + np.sqrt(n**4 - 1))) /
            (n**2 * (n + np.sqrt(1 + n**2)) * (-n**2 + np.sqrt(n**4 - 1)))
        ) / ((n**2 - 1)**2 * (n**2 + 1)**(7/2))
        t4 = ((1 + 6 * n**4 + n**8) * np.log((-1 + n) / (1 + n)) + 4 * (n**2 + n**6) * np.log((n**2 * (1 + n)) / (n - 1))) / \
            ((n**2 - 1)**2 * (n**2 + 1)**3)

        # equation (29)
        B = 1 + (3 / 2) * (2 * np.sqrt(1 - 1 / n**2) / 3 + t1 + t2 + ((1 + 6 * n**4 + n**8) * (1 - (n**2 - 1)**(3 / 2) / n**3)) / (3 * (n**4 - 1)**2) + t3)
        C = 1 - ((2 + 2 * n - 3 * n**2 + 7 * n**3 - 15 * n**4 - 19 * n**5 - 7 * n**6 + 3 * n**7 + 3 * n**8 + 3 * n**9) / (3 * n**2 * (n - 1) * (n + 1)**2 * (n**2 + 1)**2) + t4)
        A = B / C

    elif n == 1:
        # page 4591
        A = 1

    else:
        # equations (28)
        r1 = (-4 + n - 4 * n**2 + 25 * n**3 - 40 * n**4 - 6 * n**5 + 8 * n**6 + 30 * n**7 - 12 * n**8 + n**9 + n**11) / \
            (3 * n * (n**2 - 1)**2 * (n**2 + 1)**3)
        r2 = (2 * n**3 * (3 + 2 * n**4)) / ((n**2 - 1)**2 * (n**2 + 1)**(7 / 2)) * \
            np.log((n**2 * (n - np.sqrt(1 + n**2))) * (2 + n**2 + 2 * np.sqrt(1 + n**2)) / \
            ((n + np.sqrt(1 + n**2)) * (-2 + n**4 - 2 * np.sqrt(1 - n**4))))
        r3 = (4 * np.sqrt(1 - n**2) * (1 + 12 * n**4 + n**8)) / \
            (3 * n * (n**2 - 1)**2 * (n**2 + 1)**3)
        r4 = ((1 + 6 * n**4 + n**8) * np.log((1 - n) / (1 + n)) + 4 * (n**2 + n**6) * np.log((1 + n) / (n**2 * (1 - n)))) / \
            ((n**2 - 1)**2 * (n**2 + 1)**3)

        # equation (27)
        A = (1 + (3 / 2) * (8 * np.sqrt(1 - n**2) / (105 * n**3)) - (((n - 1)**2 * (8 + 32 * n + 52 * n**2 + 13 * n**3)) / (105 * n**3 * (1 + n)**2) + r1 + r2 + r3)) / \
            (1 - (-3 + 7 * n + 13 * n**2 + 9 * n**3 - 7 * n**4 + 3 * n**5 + n**6 + n**7) / (3 * (n - 1) * (n + 1)**2 * (n**2 + 1)**2) - r4)

    return A

