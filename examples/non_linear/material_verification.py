"""
Verification example for material nonlinear static analysis of a cantilever beam.
Compares BilinearMaterial with CSVMaterial and shows progressive yielding behavior.
"""

import os
import numpy as np
import pandas as pd
from fem2d import Structure, Node, BeamElementMaterialNL
from fem2d.materials import BilinearMaterial, CSVMaterial
from fem2d.sections import FiberSection


def run_cantilever_analysis(material_model, load_steps=10):
    # Geometry
    L = 2.0
    b = 0.1
    h = 0.2

    # Set up nodes
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, L, 0.0)
    node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)  # cantilever fixity

    # Set up section
    section = FiberSection.from_rectangular_fibers(
        width=b, depth=h, material=material_model, n_fibers=20
    )

    # Element
    elem = BeamElementMaterialNL(1, node1, node2, section, n_int_pts=3)

    # Structure
    structure = Structure()
    structure.add_node(node1)
    structure.add_node(node2)
    structure.add_element(elem)

    # Apply external load at tip (exceeds elastic limit: plastic moment is b*h^2/4 * sigma_y = 1e-3 * 250e6 = 250 kN-m)
    # The load P = M / L => P = 250 kN / 2 = 125 kN.
    # We apply 100 kN, which is well into the plastic range.
    node2.set_load(fx=0.0, fy=-100e3, mz=0.0)

    # Solve incrementally
    structure.solve_nonlinear(n_steps=load_steps, tolerance=1e-8, max_iter=50)

    # Displacements
    disp = structure.disp[node2.dofs]
    return disp[1]  # vertical tip displacement


def main():
    print("=====================================================================")
    print("      Verification of Material Nonlinear Beam Bending Models")
    print("=====================================================================")

    # Material parameters
    E = 200e9       # Elastic Modulus (200 GPa)
    sigma_y = 250e6 # Yield stress (250 MPa)
    E_t = 2e9       # Tangent Modulus (2 GPa, 1% hardening)

    # Create temporary CSV file for custom CSVMaterial
    csv_path = "temp_bilinear_material.csv"
    data = {
        "strain": [0.0, sigma_y / E, 0.01],
        "stress": [0.0, sigma_y, sigma_y + E_t * (0.01 - sigma_y / E)]
    }
    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)

    print(f"Created temporary CSV envelope curve in '{csv_path}'.")

    # 1. Bilinear Material Model
    mat_bilinear = BilinearMaterial(E, sigma_y, E_t, beta=0.0)
    # 2. CSV Material Model (elastoplastic with memory)
    mat_csv = CSVMaterial(csv_path, x_col="strain", y_col="stress", path_dependent=True)

    print("\nRunning analysis with BilinearMaterial...")
    disp_bilinear = run_cantilever_analysis(mat_bilinear)

    print("\nRunning analysis with CSVMaterial...")
    disp_csv = run_cantilever_analysis(mat_csv)

    print("\nResults comparison:")
    print(f"BilinearMaterial Tip vertical displacement: {disp_bilinear:.6f} m")
    print(f"CSVMaterial Tip vertical displacement:      {disp_csv:.6f} m")
    print(f"Difference:                                 {abs(disp_bilinear - disp_csv):.3e} m")

    # Clean up temporary CSV file
    if os.path.exists(csv_path):
        os.remove(csv_path)


if __name__ == "__main__":
    main()
