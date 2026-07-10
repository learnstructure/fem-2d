"""
Tests for buckling analysis functionality.
"""

import pytest
from fem2d import Structure, Node, ElasticMaterial, BeamElement, TrussElementNL, buckling_analysis


def test_frame_triangle_simple():
    # Example frame with triangle geometry: three-node triangular truss/beam frame  
    # Buckling analysis on a simple rigid-frame configuration
    
    # 1. Create nodes
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, 4.0, 3.0)
    node3 = Node(3, 8.0, 0.0)

    # 2. Material properties
    E = 70e6  
    A = 645.6e-6
    I = 0.379389e-6  # m^4

    material = ElasticMaterial(E)

    # 3. Create beam elements
    elem1 = BeamElement(1, node1, node2, material, A, I)
    elem2 = BeamElement(2, node2, node3, material, A, I)
    elem3 = BeamElement(3, node3, node1, material, A, I)

    # 4. Build structure
    structure = Structure()
    structure.add_node(node1)
    structure.add_node(node2)
    structure.add_node(node3)
    structure.add_element(elem1)
    structure.add_element(elem2)
    structure.add_element(elem3)

    # 5. Apply supports
    node1.set_support(ux_fixed=True, uy_fixed=True)  # pinned
    node3.set_support(ux_fixed=False, uy_fixed=True)  # roller

    # 6. Apply external load
    node2.set_load(fx=0.0, fy=-10.0, mz=0.0)

    # 7. Run anaysis
    structure.solve()
    
    
    # 8. Perform buckling analysis 
    factors, modes = buckling_analysis(structure, num_modes=2)
    
    # The first buckling factor should be positive and around ~2.876 for this configuration  
    assert abs(factors[0] - 2.87630061) < 1e-4
    
    # Second mode has negative factor (indicating different deformation pattern) 
    assert abs(factors[1] - (-8.3502087)) < 1e-4

if __name__ == "__main__":
    test_frame_triangle_simple()  
