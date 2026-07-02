# import numpy as np

# from fem2d import SimpleFrame, buckling_analysis


# def test_truss_buckling():
#     frame = SimpleFrame()
#     ft = 12
#     frame.add_node(1, 0, 0)
#     frame.add_node(2, 60 * ft, 0)
#     frame.add_node(3, 30 * ft, 40 * ft)
#     frame.add_node(4, 30 * ft, 60 * ft)

#     E = 30e3
#     A = 3

#     frame.add_truss(1, 1, 3, E, A)
#     frame.add_truss(2, 2, 3, E, A)
#     frame.add_truss(3, 3, 4, E, A)

#     frame.add_support(1, [1, 1, 1])
#     frame.add_support(2, [1, 1, 1])
#     frame.add_support(4, [1, 1, 1])

#     frame.add_node_load(3, [5, -10, 0])
#     frame.solve()

#     factors, modes = buckling_analysis(frame.structure, num_modes=2)

#     assert factors.shape == (2,)
#     assert modes.shape == (frame.structure.neq, 2)
#     assert np.all(factors > 0)
