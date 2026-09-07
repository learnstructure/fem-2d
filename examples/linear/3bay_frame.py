# 3-bay frame example using the SimpleFrame high-level API
# units in kN, m
#
# Geometry: 3 bays x 6 m wide, single story, 3 m tall
#    1 --- 2 --- 3 --- 4 --- 5
#    |     |     |     |     |
#    |     |     |     |     |
#    6-----7-----8-----9-----10  (ground level - supports)
#
# Loading: 10 kN/m UDL (downward) on each of the 4 beams.
from fem2d import SimpleFrame, DrawStructure
from fem2d.results import Results

frame = SimpleFrame()

# --- Geometry ---
bay_width = 6.0   # m
height = 3.0      # m

# Roof-level nodes (top of columns)
frame.add_node(1, 0 * bay_width, height)
frame.add_node(2, 1 * bay_width, height)
frame.add_node(3, 2 * bay_width, height)
frame.add_node(4, 3 * bay_width, height)

# Ground-level nodes (supports)
frame.add_node(5, 0 * bay_width, 0.0)
frame.add_node(6, 1 * bay_width, 0.0)
frame.add_node(7, 2 * bay_width, 0.0)
frame.add_node(8, 3 * bay_width, 0.0)

# --- Material & section (same for all members) ---
E = 200e6   # kN/m^2
A = 75e-3   # m^2
I = 4.8e-4  # m^4

# --- Members: 4 columns + 4 beams ---
# Columns
frame.add_frame(1, 1, 5, E, A, I)
frame.add_frame(2, 2, 6, E, A, I)
frame.add_frame(3, 3, 7, E, A, I)
frame.add_frame(4, 4, 8, E, A, I)

# Beams (3 bays = 4 beams at roof level)
frame.add_frame(5, 1, 2, E, A, I)
frame.add_frame(6, 2, 3, E, A, I)
frame.add_frame(7, 3, 4, E, A, I)

# --- Supports: fixed at all four ground nodes ---
frame.add_support(5, [1, 1, 1])
frame.add_support(6, [1, 1, 1])
frame.add_support(7, [1, 1, 1])
frame.add_support(8, [1, 1, 1])

# --- Loading: 10 kN/m UDL on every beam (downward) ---
w = -10.0  # kN/m, negative wy means downward in local axis
frame.add_distributed_load(5, wy=w)
frame.add_distributed_load(6, wy=w)
frame.add_distributed_load(7, wy=w)

# --- Solve ---
frame.solve()

# --- Results ---
results = Results(frame)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("Reactions:\n", reactions)
print("Element Forces:\n", el_forces)

# Plot deformed shape
drawer = DrawStructure(frame.structure, scale=1000)
drawer.draw(show_deformed=True)
