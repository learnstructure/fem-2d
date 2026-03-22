# Example 4.3: Logan D. L. (2022), A First Course in the Finite Element Method
# units in kips, inches

from fem2d import (
    SimpleFrame,
    Results,
)

frame = SimpleFrame()

frame.add_node(1, 0, 0)
frame.add_node(2, 3, 0)
frame.add_node(3, 6, 0)
frame.add_node(4, 6, -1)

E = 210e6
I = 2e-4
A = 100

frame.add_frame(1, 1, 2, E=E, A=A, I=I)
frame.add_frame(2, 2, 3, E=E, A=A, I=I)
frame.add_spring(3, 3, 4, 200)

frame.add_support(1, [1, 1, 1])
frame.add_support(2, [0, 1, 0])
frame.add_support(4, [1, 1, 1])

frame.add_node_load(3, [0, -50, 0])  # point load

frame.solve()

results = Results(frame)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("Reactions:\n", reactions)
print("Element Forces:\n", el_forces)

print(disp["uy"][2])  # result = -0.017441860465116272
print(el_forces["m_i"][1])  # result = 139.53488372093017
