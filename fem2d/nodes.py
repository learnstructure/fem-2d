class Node:
    def __init__(self, nid, x, y):
        self.id = nid
        self.x = x
        self.y = y
        self.support = [False, False, False]  # [ux, uy, rz] fixed?
        self.load = [0.0, 0.0, 0.0]  # nodal load [Fx, Fy, Mz]
        self.dofs = None  # will be set by Structure

    def set_support(self, ux_fixed=False, uy_fixed=False, rz_fixed=False):
        self.support = [ux_fixed, uy_fixed, rz_fixed]

    def set_load(self, fx=0.0, fy=0.0, mz=0.0):
        self.load = [fx, fy, mz]
