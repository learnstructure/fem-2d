from fem2d import SimpleFrame

frame = SimpleFrame()
frame.add_node(1, 0, 0)
frame.add_node(2, 0, 10 * 12)
frame.add_node(3, 10 * 12, 10 * 12)
frame.add_node(4, 10 * 12, 0)

frame.add_frame(1, 1, 2, 30e3, 10, 200)
frame.add_frame(2, 2, 3, 30e3, 10, 100)
frame.add_frame(3, 3, 4, 30e3, 10, 200)

frame.add_support(1, [1, 1, 1])
frame.add_support(4, [1, 1, 1])

frame.add_node_load(2, [10, 0, 0])
frame.add_node_load(3, [0, 0, 5])

disp, reactions = frame.solve()
print(disp)
