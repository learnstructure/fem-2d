# Example 5.1: Logan D. L. (2017), A First Course in the Finite Element Method
# units in kips, inches
from fem2d import SimpleFrame
from fem2d.results import Results

frame = SimpleFrame()

frame.add_node(1, 0, 0)
frame.add_node(2, 0, 10 * 12)
frame.add_node(3, 10 * 12, 10 * 12)
frame.add_node(4, 10 * 12, 0)

E, A = 30e3, 10
I = 200

frame.add_frame(1, 1, 2, E, A, I)
frame.add_frame(2, 2, 3, E, A, I / 2)
frame.add_frame(3, 3, 4, E, A, I)

frame.add_support(1, [1, 1, 1])
frame.add_support(4, [1, 1, 1])

frame.add_node_load(2, [10, 0, 0])
frame.add_node_load(3, [0, 0, 5])

frame.solve()

results = Results(frame)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("Reactions:\n", reactions)
print("Element Forces:\n", el_forces)

print(disp["theta"][2])  # result = -0.0014859999862147104
print(el_forces["m_i"][2])  # result = 226.19833955969983
