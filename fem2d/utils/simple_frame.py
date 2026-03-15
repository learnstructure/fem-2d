from fem2d.elements.truss import TrussElement
from fem2d.nodes import Node
from fem2d.materials import ElasticMaterial
from fem2d.sections import Section
from fem2d.elements import BeamElement
from fem2d.structure import Structure
from fem2d.loads import DistributedLoad


class SimpleFrame:
    def __init__(self):
        self.structure = Structure()
        self._node_counter = 0
        self._elem_counter = 0

    def add_node(self, id, x, y):
        node = Node(id, x, y)
        self.structure.add_node(node)
        return node

    def add_frame(self, id, node_i_id, node_j_id, E, A, I):
        node_i = self.structure.nodes[node_i_id]
        node_j = self.structure.nodes[node_j_id]
        material = ElasticMaterial(E)
        section = Section(A, I)
        elem = BeamElement(id, node_i, node_j, material, A, I)
        self.structure.add_element(elem)

    def add_truss(self, id, node_i_id, node_j_id, E, A):
        node_i = self.structure.nodes[node_i_id]
        node_j = self.structure.nodes[node_j_id]
        material = ElasticMaterial(E)
        # section = Section(A, I)
        elem = TrussElement(id, node_i, node_j, material, A)
        self.structure.add_element(elem)

    def add_support(self, node_id, fixity):
        node = self.structure.nodes[node_id]
        node.set_support(fixity[0], fixity[1], fixity[2])

    def add_node_load(self, node_id, load):
        node = self.structure.nodes[node_id]
        node.set_load(load[0], load[1], load[2])

    def add_distributed_load(self, element_id, wx=0.0, wy=0.0):
        """Add a uniformly distributed load to an element (local axes)."""
        element = self.structure.elements[element_id]
        load = DistributedLoad(element, wx, wy)
        self.structure.add_load(load)
        self.structure.add_load(load)

    def solve(self):
        self.structure.solve()
        return self.structure.disp, self.structure.reactions
