#
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt


def stiffness_matrix(E1, E2, v12, G12):
    """Returns stiffness matrix in material CS.

    Parameters
    ----------
    E1 : float
        Young's modulus in 1-direction.
    E2 : float
        Young's modulus in 2-direction.
    v12 : float
        Poisson's ratio.
    G12 : float
        Shear modulus.

    Returns
    -------
    C : NDArray(dtype=float, dim=2)
        Stiffness matrix (plane stress) in material coordinate
        system.

    """
    v21 = E2*v12/E1
    C = np.array([[E1/(1-v12*v21), v21*E1/(1-v12*v21), 0],
                  [v21*E1/(1-v12*v21), E2/(1-v12*v21), 0],
                  [0, 0, G12]])
    return C


def transformation_matrix(theta):
    """Returns transformation matrix.

    Parameters
    ----------
    theta : float
        Rotation angle in radians.

    Returns
    -------
    T : NDArray(dtype=float, dim=2)
        2D transformation matrix.

    """
    n = np.sin(theta)
    m = np.cos(theta)
    T = np.array([[m**2, n**2,     2*n*m],
                  [n**2, m**2,    -2*n*m],
                  [-m*n,  m*n, m**2-n**2]])
    return T


def Reuter():
    """Returns Reuter matrix."""
    R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    return R


def rotate_C(C, theta):
    """Returns stiffness matrix or matrices in ply CS.

    Parameters
    ----------
    C : NDArray(dtype=floats, dim=2)
        Stiffness matrix in material CS.
    theta : float or list
        Rotation angle or, in case of a laminate, list of rotation
        angles.

    Returns
    -------
    C_r : NDArray(dtype=float, dim=2) or list
        Stiffness matrix in ply CS or, in case of a theta is a list of
        angles, a list with stiffness matrices in ply CS.

    """
    R = Reuter()
    if isinstance(theta, float):
        T = transformation_matrix(theta)
        C_r = inv(T)@C@R@T@inv(R)
    elif isinstance(theta, list):
        C_r = [None]*len(theta)
        for i, phi in enumerate(theta):
            T = transformation_matrix(phi)
            C_r[i] = inv(T)@C@R@T@inv(R)
    else:
        raise TypeError("Please input a float or a list of floats.")
    return C_r


def rotate_alpha(alpha, theta):
    """Returns CTE's in ply CS.

    Parameters
    ----------
    alpha : NDArray(dtype=float, dim=1)
        CTEs in material CS.
    theta : float or list
        Rotation angle or, in case of a laminate, list of rotation
        angles (in radians).

    Returns
    -------
    alpha_r : NDArray(dype=float, dim=1) or list
        CTEs in ply CS or, in case of theta is a list of angles, a list
        with CTEs in ply CS for each ply.

    """
    R = Reuter()
    if isinstance(theta, float):
        T = transformation_matrix(theta)
        alpha_r = R@inv(T)@inv(R)@alpha
    elif isinstance(theta, list):
        alpha_r = [None]*len(theta)
        for i, phi in enumerate(theta):
            T = transformation_matrix(phi)
            alpha_r[i] = R@inv(T)@inv(R)@alpha
    else:
        raise TypeError("Please input a float or a list of floats.")
    return alpha_r


def ply_edges(h, N):
    """Returns location of ply edges.

    Parameters
    ----------
    h : float
        Ply thickness.
    N : int
        Number of plies.

    Returns
    -------
    z : NDArray(dtype=float, dim=1)
        Array with z-coordinates of the ply edges.

    """
    z = np.linspace(-N*h/2, N*h/2, N+1)
    return z


def ABD_matrix(C_r, z):
    """Returns ABD matrix.

    Parameters
    ----------
    C_r : list
        List of length N with stiffness matrices in ply CS for each ply.
    z : NDArray(dtype=float, dim=1)
        Location of the ply edges.

    Returns
    -------
    ABD : NDArray(dtype=float, dim=1)
        ABD matrix.

    """
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    for i in range(len(C_r)):
        A = A + C_r[i] * (z[i+1] - z[i])
        B = B + C_r[i] * (z[i+1]**2 - z[i]**2)/2
        D = D + C_r[i] * (z[i+1]**3 - z[i]**3)/3
    ABD = np.block([[A, B], [B, D]])
    return ABD


def thermal_force(C_r, alpha_r, z, deltaT):
    """Returns fictive thermal forces and moments.

    Parameters
    ----------
    C_r : list
        List of stiffness matrices for each ply.
    alpha_r : list
        List of CTE vectors for each ply.
    z : NDArray(dtype=float, dim=1)
        Array with z-coordinates of the ply edge locations.
    deltaT : float
        Temperature difference.

    Returns
    -------
    NM : NDArray(dtype=float, dim=1)
        Array with fictive in-plane loads and moments.

    """
    N = np.zeros(3)
    M = np.zeros(3)
    for i in range(len(C_r)):
        N = N + deltaT * C_r[i]@alpha_r[i] * (z[i+1]-z[i])
        M = M + 0.5*deltaT * C_r[i]@alpha_r[i] * (z[i+1]**2-z[i]**2)
    NM = np.block([N, M])
    return NM


def ply_stress(d, C_r, z, *argv):
    """Returns ply stresses in ply CS.

    Parameters
    ----------
    d : NDArray(dtype=float, dim=1)
        Laminate deformation vector with strains and curvatures.
    C_r : list
        List of stiffness matrices for each ply.
    z : NDArray(dtype=float, dim=1)
        Array with ply edge locations.
    alpha_r : list (optional)
        List with CTE vectors for each ply.
    deltaT : float (optional)
        Temperature difference.

    Returns
    -------
    stress_r : NDArray(dtype=float, dim=1)
        Array with stress state in ply CS at top and bottom of each ply.
        Column number i*2 gives the stress at the top of the i-th ply, while
        column i*2+1 gives the stress at the bottom of the i-th ply. Ply count
        starts at i=0.
    z_int : ndarray(2*N-2)
        Ply edge locations for each ply.

    """
    if len(argv) == 0:
        alpha_r = [np.array([0.0, 0.0, 0.0])]*len(C_r)
        deltaT = 0.0
    elif len(argv) == 2:
        alpha_r = argv[0]
        deltaT = argv[1]

    z_int = np.repeat(z, 2)[1:-1]
    stress_r = np.zeros((3, len(z_int)))
    eps0 = d[:3]
    kappa = d[3:]

    for i in range(len(C_r)):
        eps_top = eps0 + kappa*z_int[i*2] - alpha_r[i]*deltaT
        stress_r[:, i*2] = C_r[i]@eps_top
        eps_bot = eps0 + kappa*z_int[i*2+1] - alpha_r[i]*deltaT
        stress_r[:, i*2+1] = C_r[i]@eps_bot
    return stress_r, z_int


def rotate_stress_to_matCS(stress_r, theta):
    """Rotates stress from ply to material CS.

    Parameters
    ----------
    stress_r : NDArray(dtype=float, dim=2)
        Array of size (3, 2*N) with the stress in ply CS for each ply.
        Column 2*i and 2*i-1 provide the stress at the top and bottom
        of the i-th ply.
    theta : float or list
        Rotation angle or, in case of a laminate, list of rotation
        angles (in radians).

    Returns
    -------
    stress : NDArray(dtype=float, dim=2)
        Array with stresses in material CS.

    """
    stress = np.zeros(stress_r.shape)
    if isinstance(theta, float):
        T = transformation_matrix(theta)
        stress = T@stress_r
    elif isinstance(theta, list):
        for i, phi in enumerate(theta):
            T = transformation_matrix(phi)
            stress[:, i*2] = T@stress_r[:, i*2]
            stress[:, i*2+1] = T@stress_r[:, i*2+1]
    return stress


def CP_layup(N):
    """Returns layup for symmetric cross-ply laminate.

    Parameters
    ----------
    N : int
        Number of plies.

    Returns
    -------
    layup : list
        List with layup angles in radians.

    """
    if N % 4 == 0:
        unit = [0.0, np.pi/2]
        k = int(N/4)
        half = unit*k
        layup = half[:] + half[::-1]
    else:
        raise ValueError("The number of plies should be a multiple of 4.")
    return layup


def QI_layup(N):
    """Returns layup for symmetric quasi-isotropic laminate.

    Parameters
    ----------
    N : int
        Number of plies.

    Returns
    -------
    layup : list
        List with layup angles in radians.

    """
    if N % 8 == 0:
        unit = [np.pi/4, 0.0, -np.pi/4, np.pi/2]
        k = int(N/8)
        half = unit*k
        layup = half[:] + half[::-1]
    else:
        raise ValueError("The number of plies should be a multiple of 8.")
    return layup


def max_stress_criterion(stress, strength):
    """Returns stress components that exceed strength.

    Parameters
    ----------
    stress : NDArray(dtype=float, dim=1) or NDArray(dtype=float, dim=2)
        Array with the stress in material CS. Can also be an array with
        stress in material CS for each ply. In this case, column 2*i
        and 2*i-1 provide the stress at the top and bottom of the i-th ply.
    strength : NDArray(dtype=float, dim=2)
        Array of size (2, 3) with the strength values. The first column
        provides the compressive strengths in (1, 2, 6) direction, while
        the second column provides the tensile strengths.

    Returns
    -------
    failed : NDArray(dtype=int, dim=1) or NDArray(dtype=int, dim=2)
        Array of size similar to stress where -1 indicates a stress
        component exceeding the compressive strength and 1 indicating
        a stress component exceeding the tensile strength.

    """
    stress = np.expand_dims(stress, 1) if stress.ndim == 1 else stress
    failed = np.zeros(stress.shape, dtype=int)
    for i, s in enumerate(stress.T):
        failed[:, i] = failed[:, i] - 1*(s < strength[:, 0])
        failed[:, i] = failed[:, i] + 1*(s > strength[:, 1])
    return failed


def TsaiHill_criterion(stress, strength):
    """Returns stress components that exceed strength.

    Parameters
    ----------
    stress : NDArray(dtype=float, dim=1) or NDArray(dtype=float, dim=2)
        Array with the stress in material CS. Can also be an array with
        stress in material CS for each ply. In this case, column 2*i
        and 2*i-1 provide the stress at the top and bottom of the i-th ply.
    strength : NDArray(dtype=float, dim=2)
        Array of size (2, 3) with the strength values. The first column
        provides the compressive strengths in (1, 2, 6) direction, while
        the second column provides the tensile strengths.

    Returns
    -------
    failed : NDArray(dtype=int, dim=1)
        Array indicating whether (and where) failure occurred according
        to Tsai Hill.

    """
    stress = np.expand_dims(stress, 1) if stress.ndim == 1 else stress
    failed = np.zeros(stress.shape[1], dtype=int)
    for i, s in enumerate(stress.T):
        S1 = strength[0, 0] if s[0] < 0 else strength[0, 1]
        S2 = strength[1, 0] if s[1] < 0 else strength[1, 1]
        S6 = strength[2, 0]
        TH = s[0]**2/S1**2 - s[0]*s[1]/S1**2 + s[1]**2/S2**2 + s[2]**2/S6**2
        failed[i] = 1 if TH > 1 else 0
    return failed


def plot_stress(stress, z_int, i):
    """Plots stress distribution through thickness

    Parameters
    ----------
    stress : NDArray(dtype=float, dim=2)
        Array with stress for each ply. Column 2*i and 2*i-1 provide the
        stress at the top and bottom of the i-th ply.
    z_int : NDArray(dtype=float, dim=1)
        Z-coordinates of the ply-edges where the stresses are provided.
    i : int (0, 1, 2)
        Stress component to plot.

    """
    fig = plt.figure()
    ax = fig.gca()
    ax.plot([0, 0], [z_int.min(), z_int.max()], 'k:')
    ax.plot(stress[i], z_int)
    ax.set_xlabel("Stress, Pa")
    ax.set_ylabel("Z-coordinate, m")
    plt.show()
