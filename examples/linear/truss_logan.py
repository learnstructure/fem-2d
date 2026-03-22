# Problem: Logan D. L. (2017), A First Course in the Finite Element Method
# units in kips, inches
from fem2d import SimpleFrame
from fem2d.results import Results

frame = SimpleFrame()
ft = 12
frame.add_node(1, 0, 0)
frame.add_node(2, 60 * ft, 0)
frame.add_node(3, 30 * ft, 40 * ft)
frame.add_node(4, 30 * ft, 60 * ft)

E = 30e3
A = 3

frame.add_truss(1, 1, 3, E, A)
frame.add_truss(2, 2, 3, E, A)
frame.add_truss(3, 3, 4, E, A)

frame.add_support(1, [1, 1, 1])
frame.add_support(2, [1, 1, 1])
frame.add_support(4, [1, 1, 1])

frame.add_node_load(3, [5, -10, 0])

frame.solve()

results = Results(frame.structure)
print("Node Displacements:\n", results.node_displacements())
print("Reactions:\n", results.reactions())
print("Element Forces:\n", results.element_forces())

print(
    "Node 3 Displacement:\n", (results.node_displacements())["uy"][2]
)  # result is -0.01763668430335097
print(
    "Element 2 Force:\n", results.element_forces()["fx_i"][0]
)  # result is -2.0502645502645502
