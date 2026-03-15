# Example 5.2: Logan D. L. (2017), A First Course in the Finite Element Method
# units in kips, inches
from fem2d import SimpleFrame
from fem2d.results import Results
import math

frame = SimpleFrame()

# units in kips, in
frame.add_node(1, 0, 0)
frame.add_node(2, 30 * 12, 30 * 12)
frame.add_node(3, 70 * 12, 30 * 12)

E = 30 * 1000
I = 1000
A = 100

frame.add_frame(1, 1, 2, E=E, A=A, I=I)
frame.add_frame(2, 2, 3, E=E, A=A, I=I)

frame.add_support(1, [1, 1, 1])
frame.add_support(3, [1, 1, 1])

frame.add_distributed_load(2, wy=-1 / 12)  # udl
frame.solve()

results = Results(frame)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("Reactions:\n", reactions)
print("Element Forces:\n", el_forces)

print(disp["theta"][1])  # result = -0.0032917095717463736
print(el_forces["fy_i"][1])  # result = 17.39663896899728
