import numpy as np


class NewtonRaphsonSolver:
    def __init__(self, structure, tolerance=1e-8, max_iter=30):
        self.structure = structure
        self.tol = tolerance
        self.max_iter = max_iter

    def solve(self, P_ext):
        d = np.zeros(self.structure.neq)  # current displacements

        for iteration in range(self.max_iter):
            # Assemble internal forces and tangent stiffness
            f_int = np.zeros(self.structure.neq)
            K_t = np.zeros((self.structure.neq, self.structure.neq))

            for el in self.structure.elements.values():
                el.update_state(d)
                dofs = el.node_i.dofs + el.node_j.dofs  # list of 6 DOFs
                f_int[dofs] += el.get_internal_forces()  # 6‑element vector
                K_t[np.ix_(dofs, dofs)] += el.get_tangent_stiffness()  # 6×6 matrix

            # Unbalanced forces on free DOFs
            delta_U = P_ext[self.structure.free_dofs] - f_int[self.structure.free_dofs]

            # Check convergence (norm of unbalanced forces)
            if np.linalg.norm(delta_U) < self.tol * np.linalg.norm(
                P_ext[self.structure.free_dofs]
            ):
                print(f"Converged in {iteration+1} iterations")
                break

            # Solve for displacement increment
            K_ff = K_t[np.ix_(self.structure.free_dofs, self.structure.free_dofs)]
            delta_d_f = np.linalg.solve(K_ff, delta_U)

            # Update displacements
            d[self.structure.free_dofs] += delta_d_f

        else:
            print("Warning: did not converge within maximum iterations")

        return d
