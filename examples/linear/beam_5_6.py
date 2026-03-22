# Example 5.6: Kassimali A. (2022), Matrix analysis of structures
# units in kN, m

from fem2d import SimpleFrame, DrawStructure
from fem2d.results import Results

frame = SimpleFrame()

frame.add_node(1, 0, 0)
frame.add_node(2, 4, 0)
frame.add_node(3, 7, 0)
frame.add_node(4, 11, 0)

E = 200e6
I = 108e-6
A = 10

frame.add_frame(1, 1, 2, E=E, A=A, I=I)
frame.add_frame(2, 2, 3, E=E, A=A, I=I)
frame.add_frame(3, 3, 4, E=E, A=A, I=I)

frame.add_support(1, [1, 1, 1])
frame.add_support(2, [0, 1, 0])
frame.add_support(3, [0, 1, 0])
frame.add_support(4, [1, 1, 1])

frame.add_element_point_load(1, py=-150, x=2)
frame.add_distributed_load(3, wy=-37.5)
frame.solve()

results = Results(frame)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("Reactions:\n", reactions)
print("Element Forces:\n", el_forces)

print(disp["theta"][1])  # result = 0.0019290123456790125
print(el_forces["fy_i"][1])  # result = 5.555555555555557

# drawer = DrawStructure(
#     frame.structure, scale=5
# )  # scale displacements by 100 for visibility
# drawer.draw()
