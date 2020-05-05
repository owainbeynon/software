def switch_indices(model, A, B):
    '''
    Function for rearranging atomic indices in the structure.
    Returns a changed atoms object with retained calculator
    and ensures array of forces is rearranged accordingly.

    Parameters:

    model: Atoms object or string
        Structure requiring atom index rearrangement.
        If a string, e.g. 'name.traj', a file of this name will
        be read.
    A, B: integer
        Indices of atoms that need to be switched

    '''

    from ase import Atoms

    if isinstance(model, str):
        from ase.io import read
        model = read(model)

    # Index manipulation
    if not isinstance(A, int) and isinstance(B, int):
        raise ValueError
        print("Indices must be integers.")

    # Defining list of manipulated atoms
    t = [atom.index for atom in model]

    # Retrieve calculator information
    # TODO: Move this after the creation of the atoms object, so we reduce if statements.
    if model.get_calculator() is not None:
        # If it exists, forces array needs to be adjusted.
        prev_calc = model.get_calculator()
        prev_calc_results = prev_calc.results
        f = prev_calc_results["forces"]
    else:
        f = []

    t[A], t[B] = t[B], t[A]
    # Ensure function works if force information empty
    # TODO: Where is this information subsequently used?
    if not f == []:
        f[A], f[B] = f[B], f[A]

    # Generate a new model based on switched indices
    new_model = model[t]

    if model.get_calculator() is not None:
        new_model.set_calculator(prev_calc)
        # Trick calculator check_state by replacing atoms information
        # Can now use energy and forces as no changes in geometry detected
        prev_calc.atoms = new_model

    # User can interact with the new model
    return new_model

def check_interpolation(initial, final, n_max, interpolation="linear", verbose=True, save=True):
    '''
    Interpolates the provided geometries with n_max total images
    and checks whether any bond lengths are below sane defaults.
    Saves the interpolation in interpolation.traj

    Parameters:

    initial: Atoms object or string
        If a string, e.g. 'initial.traj', a file of this name will
        be read. Starting geometry for interpolation.
    final: Atoms object or string
        If a string, e.g. 'final.traj', a file of this name will
        be read. End point geometry for interpolation
    n_max: integer
        Desired total number of images for the interpolation
        including start and end point.
    interpolation: string
        "linear" or "idpp". First better for error identification, latter for
        use in NEB calculation
    verbose: boolean
        If verbose output of information is required
    save: boolean
        Whether to save the trajectory for transfer on to an NEB calculation
    '''

    from ase.neb import NEB
    from bonds import search_abnormal_bonds
    from ase.io.trajectory import Trajectory
    from ase.io import read

    # Pre-requirements
    if isinstance(initial, str):
        initial = read(initial)
    if isinstance(final, str):
        final = read(final)
    if not isinstance(n_max, int):
        raise ValueError
        print("Max number of images must be an integer.")

    # Make a band consisting of 10 images:
    images = [initial]
    images += [initial.copy() for i in range(n_max-2)]
    images += [final]
    neb = NEB(images, climb=True)
    # Interpolate linearly the potisions of the middle images:
    neb.interpolate(interpolation)

    #TODO: Tidy up this horrible mix of if statements.
    if save:
        t = Trajectory('interpolation.traj', 'w')

    flag = True
    for i in range(0, n_max):
        if verbose:
            print("Assessing image", str(i+1) + '.')
        updated_flag = search_abnormal_bonds(images[i], verbose)
        if save:
            t.write(images[i])
        if (not updated_flag):
            flag = updated_flag

    if save:
        t.close()

    return flag
