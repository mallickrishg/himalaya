import numpy as np
import bemcs


def logisticfunction(x, k=1, x0=0):
    return 1 / (1 + np.exp(-k * (x - x0)))


def construct_topographymesh(height, dip, fault_x, nfault):
    # construct topography
    xmesh = np.concatenate(
        (-np.logspace(3, -0.5, 20), np.linspace(0, 95, 50), np.logspace(2, 3, 60))
    )
    ymesh = np.zeros_like(xmesh)
    ymesh[xmesh > 0] = height * logisticfunction(xmesh[xmesh > 0], k=0.1, x0=100)
    ntopo = len(xmesh) - 1
    x1 = xmesh[0:-1]
    x2 = xmesh[1:]
    y1 = ymesh[0:-1]
    y2 = ymesh[1:]

    # add fault
    # dip = 12  # in degrees
    xmesh = np.linspace(0, fault_x, nfault + 1)
    ymesh = -np.tan(np.deg2rad(dip)) * np.linspace(0, 100, nfault + 1)
    x1 = np.concatenate((x1, xmesh[0:-1]))
    x2 = np.concatenate((x2, xmesh[1:]))
    y1 = np.concatenate((y1, ymesh[0:-1]))
    y2 = np.concatenate((y2, ymesh[1:]))

    # create labels for each element
    labels = np.concatenate((np.repeat("topo", ntopo), np.repeat("fault", nfault)))
    # create vector of boundary conditions
    bctype = np.concatenate((np.repeat("t_local", ntopo), np.repeat("s_local", nfault)))
    bc_x = np.zeros(nfault + ntopo)
    bc_y = np.zeros(nfault + ntopo)
    bc_x[labels == "fault"] = 1

    # create bemcs data structure
    els = bemcs.initialize_els()
    els.x1 = x1
    els.y1 = y1
    els.x2 = x2
    els.y2 = y2
    bemcs.standardize_els_geometry(els)

    return els, labels, bctype, bc_x, bc_y


def setup_and_solve_bem(els, labels, bctype, bc_x, bc_y, mu=1, nu=0.25):
    n_els = len(els.x1)
    index_open, index_overlap, index_triple = bemcs.label_nodes(els)
    N_c = 2 * n_els  # central node equations
    N_o = 2 * len(index_open)  # open node equations
    N_i = 4 * len(index_overlap)  # overlapping node equations
    N_t = 6 * len(index_triple)  # triple junction equations

    Nequations = N_c + N_o + N_i + N_t
    Nunknowns = 6 * n_els

    # We will stack this with
    # equations for the element centers
    # equations at open nodes (RHS = 0)
    # equations at overlapping nodes (RHS = 0)
    # equations at triple junctions (RHS = 0)
    BC_c = np.zeros((N_c, 1))  # these are the only non-zero entries
    BC_o = np.zeros((N_o, 1))
    BC_i = np.zeros((N_i, 1))
    BC_t = np.zeros((N_t, 1))

    # apply BCs at central nodes
    BC_c[0::2, 0] = bc_x
    BC_c[1::2, 0] = bc_y

    # stack all the BCs into 1 big vector
    BCvector = np.vstack((BC_c, BC_o, BC_i, BC_t))

    # Design matrices (in x,y coordinates) for slip and slip gradients at each 3qn
    matrix_slip, matrix_slip_gradient = bemcs.get_matrices_slip_slip_gradient(
        els, reference="local"
    )

    # Patch center locations
    # (need to be shifted an infinitesimal amount in unit normal direction for displacement bcs)
    obs_xy = np.vstack((els.x_centers, els.y_centers)).T
    x_obs = (obs_xy[:, 0]).reshape(-1, 1)
    y_obs = (obs_xy[:, 1]).reshape(-1, 1)

    # Compute shear and tensile stress kernels evaluated ONLY at the center of each element
    kernels_s = bemcs.get_displacement_stress_kernel(x_obs, y_obs, els, mu, nu, "shear")
    kernels_n = bemcs.get_displacement_stress_kernel(
        x_obs, y_obs, els, mu, nu, "normal"
    )

    # Convert to traction kernels [Nobs x Ncoefficients]
    traction_kernels_s = bemcs.get_traction_kernels(els, kernels_s, flag="local")
    traction_kernels_n = bemcs.get_traction_kernels(els, kernels_n, flag="local")

    # Linear operator for central node BCs
    kerneleval_x = np.zeros((n_els, Nunknowns))
    kerneleval_y = np.zeros((n_els, Nunknowns))

    # x,y-kernels
    for i in np.unique(labels):
        index = np.where(labels == i)[0]
        for j in index:
            if bctype[j] == "u_global":
                for k in range(0, 3):
                    kerneleval_x[j, k::6] = kernels_s[3][j, k::3]
                    kerneleval_x[j, k + 3 :: 6] = kernels_n[3][j, k::3]
            elif bctype[j] == "t_local":
                for k in range(0, 3):
                    kerneleval_x[j, k::6] = traction_kernels_s[0][j, k::3]
                    kerneleval_x[j, k + 3 :: 6] = traction_kernels_n[0][j, k::3]
            elif bctype[j] == "s_local":
                kerneleval_x[j, :] = matrix_slip[2::6, :][j, :]
            else:
                raise ValueError("unrecognized boundary condition type")

            if bctype[j] == "u_global":
                for k in range(0, 3):
                    kerneleval_y[j, k::6] = kernels_s[4][j, k::3]
                    kerneleval_y[j, k + 3 :: 6] = kernels_n[4][j, k::3]
            elif bctype[j] == "t_local":
                for k in range(0, 3):
                    kerneleval_y[j, k::6] = traction_kernels_s[1][j, k::3]
                    kerneleval_y[j, k + 3 :: 6] = traction_kernels_n[1][j, k::3]
            elif bctype[j] == "s_local":
                kerneleval_y[j, :] = matrix_slip[3::6, :][j, :]
            else:
                raise ValueError("unrecognized boundary condition type")

    # Linear Operators for the appropriate boundary conditions
    matrix_system_c = np.zeros((N_c, Nunknowns))
    # populate matrix_system for central nodes
    matrix_system_c[0::2, :] = kerneleval_x
    matrix_system_c[1::2, :] = kerneleval_y

    matrix_system_o, matrix_system_i, matrix_system_t = bemcs.construct_smoothoperator(
        els, index_open, index_overlap, index_triple
    )

    # stack the matrices and create the full linear operator
    matrix_system = np.vstack(
        (matrix_system_c, matrix_system_o, matrix_system_i, matrix_system_t)
    )

    # compute quadratic node coefficients (in local (s,n) coordinates)
    quadratic_coefs = np.linalg.inv(matrix_system) @ BCvector
    print("Linear Operator Condition Number:", np.linalg.cond(matrix_system))

    # extract (s,n) components and store them in 2 separate vectors
    quadratic_coefs_s = np.zeros((3 * n_els, 1))
    quadratic_coefs_n = np.zeros((3 * n_els, 1))
    for i in range(n_els):
        quadratic_coefs_s[3 * i : 3 * (i + 1)] = quadratic_coefs[6 * i : 6 * i + 3]
        quadratic_coefs_n[3 * i : 3 * (i + 1)] = quadratic_coefs[
            6 * i + 3 : 6 * (i + 1)
        ]

    return quadratic_coefs, quadratic_coefs_s, quadratic_coefs_n
