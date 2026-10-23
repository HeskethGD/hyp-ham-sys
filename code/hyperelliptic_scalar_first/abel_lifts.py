"""Shared integer period coordinates for compatible Abel lifts."""


def integer_lattice_coordinates(model, vector):
    c, g, periods = model.ctx, model.genus, model.periods
    lattice = c.matrix(g, 2*g)
    lattice[:, :g] = 2*periods.omega
    lattice[:, g:] = 2*periods.omega_prime
    real_system = c.matrix(2*g, 2*g)
    for i in range(g):
        for j in range(2*g):
            real_system[i, j] = lattice[i, j].real
            real_system[g+i, j] = lattice[i, j].imag
    target = c.matrix([z.real for z in vector]+[z.imag for z in vector])
    raw = c.lu_solve(real_system, target)
    integers = [int(c.nint(z)) for z in raw]
    residual = vector-lattice*c.matrix(integers)
    return integers, residual


