# Example 2.1: Logan D. L. (2017), A First Course in the Finite Element Method
# units in pound, inches

from fem2d import SimpleFrame, DrawStructure
from fem2d.results import Results

frame = SimpleFrame()

frame.add_node(1, 0, 0)
frame.add_node(2, 1, 0)
frame.add_node(3, 2, 0)
frame.add_node(4, 3, 0)

k1, k2, k3 = 1000, 2000, 3000

frame.add_spring(1, 1, 2, k1)
frame.add_spring(2, 2, 3, k2)
frame.add_spring(3, 3, 4, k3)

frame.add_support(1, [1, 1, 1])
frame.add_support(4, [1, 1, 1])

frame.add_node_load(3, [5000, 0, 0])

frame.solve()

results = Results(frame)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("Reactions:\n", reactions)
print("Element Forces:\n", el_forces)

print(disp["ux"][1])  # result = 0.9090909090909091
print(el_forces["fx_i"][1])  # result = -909.090909090909

# drawer = DrawStructure(frame.structure, scale=1)
# drawer.draw()
