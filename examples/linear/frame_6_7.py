# Example 6.7: Kassimali A. (2022), Matrix analysis of structures

# units in kN, m
from fem2d import SimpleFrame, DrawStructure
from fem2d.results import Results

frame = SimpleFrame()

frame.add_node(1, 0, 0)
frame.add_node(2, 0, 6)
frame.add_node(3, 0, 12)
frame.add_node(4, 9, 6)
frame.add_node(5, 9, 0)

E = 30e6
I = 4.8e-4
A = 75e-3

frame.add_frame(1, 1, 2, E, A, I)
frame.add_frame(2, 2, 3, E, A, I)
frame.add_frame(3, 3, 4, E, A, I)
frame.add_frame(4, 2, 4, E, A, I)
frame.add_frame(5, 4, 5, E, A, I)

frame.add_support(1, [1, 1, 1])
frame.add_support(5, [1, 1, 1])

frame.add_node_load(2, [80, 0, 0])  # point load
frame.add_node_load(3, [40, 0, 0])
frame.add_distributed_load(3, wy=12)  # udl

frame.solve()

results = Results(frame)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

# print("Node Displacements:\n", disp)
# print("Reactions:\n", reactions)
# print("Element Forces:\n", el_forces)

print(disp["theta"][2])  # result = 0.01789119564872733
print(el_forces["m_i"][2])  # result = -90.05939530594357

# drawer = DrawStructure(
#     frame.structure, scale=5
# )  # scale displacements by 100 for visibility
# drawer.draw()
