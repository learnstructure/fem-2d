"""
SimpleFrame module providing a high-level API for creating and solving 2D structural frames.
"""

from fem2d.elements.truss import TrussElement
from fem2d.nodes import Node
from fem2d.materials import ElasticMaterial
from fem2d.sections import Section
from fem2d.elements import BeamElement, SpringElement
from fem2d.structure import Structure
from fem2d.loads import DistributedLoad, ElementPointLoad


class SimpleFrame:
    """
    A high-level interface/wrapper for building and solving 2D frame structures.

    Simplifies the creation of nodes, frame members (beams), trusses, springs,
    supports, and load applications.

    Attributes
    ----------
    structure : Structure
        The underlying Structure object that holds the nodes, elements, and loads.
    """

    def __init__(self):
        """Initialize a SimpleFrame wrapper with an empty Structure."""
        self.structure = Structure()
        self._node_counter = 0
        self._elem_counter = 0

    def add_node(self, id, x, y):
        """
        Add a node to the frame.

        Parameters
        ----------
        id : int or str
            Unique identifier of the node.
        x : float
            x-coordinate of the node.
        y : float
            y-coordinate of the node.

        Returns
        -------
        Node
            The created and added Node object.
        """
        node = Node(id, x, y)
        self.structure.add_node(node)
        return node

    def add_frame(self, id, node_i_id, node_j_id, E, A, I):
        """
        Add an elastic beam/frame element to the structure.

        Parameters
        ----------
        id : int or str
            Unique identifier of the element.
        node_i_id : int or str
            ID of the start node.
        node_j_id : int or str
            ID of the end node.
        E : float
            Young's Modulus of the member material.
        A : float
            Cross-sectional area.
        I : float
            Moment of inertia of the member cross-section.
        """
        node_i = self.structure.nodes[node_i_id]
        node_j = self.structure.nodes[node_j_id]
        material = ElasticMaterial(E)
        section = Section(A, I)
        elem = BeamElement(id, node_i, node_j, material, A, I)
        self.structure.add_element(elem)

    def add_truss(self, id, node_i_id, node_j_id, E, A):
        """
        Add an elastic truss element to the structure.

        Parameters
        ----------
        id : int or str
            Unique identifier of the element.
        node_i_id : int or str
            ID of the start node.
        node_j_id : int or str
            ID of the end node.
        E : float
            Young's Modulus of the truss material.
        A : float
            Cross-sectional area.
        """
        node_i = self.structure.nodes[node_i_id]
        node_j = self.structure.nodes[node_j_id]
        material = ElasticMaterial(E)
        # section = Section(A, I)
        elem = TrussElement(id, node_i, node_j, material, A)
        self.structure.add_element(elem)

    def add_spring(self, id, node_i_id, node_j_id, stiffness):
        """
        Add a 2D spring element to the structure.

        Parameters
        ----------
        id : int or str
            Unique identifier of the element.
        node_i_id : int or str
            ID of the start node.
        node_j_id : int or str
            ID of the end node.
        stiffness : float
            Axial stiffness coefficient of the spring.
        """
        node_i = self.structure.nodes[node_i_id]
        node_j = self.structure.nodes[node_j_id]
        elem = SpringElement(id, node_i, node_j, stiffness)
        self.structure.add_element(elem)

    def add_support(self, node_id, fixity):
        """
        Apply boundary conditions (supports) to a node.

        Parameters
        ----------
        node_id : int or str
            ID of the node to constrain.
        fixity : list or tuple of bool
            Fixity flags [ux_fixed, uy_fixed, rz_fixed].
        """
        node = self.structure.nodes[node_id]
        node.set_support(fixity[0], fixity[1], fixity[2])

    def add_node_load(self, node_id, load):
        """
        Apply concentrated forces and moment to a node.

        Parameters
        ----------
        node_id : int or str
            ID of the node.
        load : list or tuple of float
            Nodal force and moment values [Fx, Fy, Mz].
        """
        node = self.structure.nodes[node_id]
        node.set_load(load[0], load[1], load[2])

    def add_distributed_load(self, element_id, wx=0.0, wy=0.0):
        """
        Add a uniformly distributed load to an element (local axes).

        Parameters
        ----------
        element_id : int or str
            ID of the element.
        wx : float, optional
            Axial load intensity per unit length. Defaults to 0.0.
        wy : float, optional
            Transverse load intensity per unit length. Defaults to 0.0.
        """
        element = self.structure.elements[element_id]
        load = DistributedLoad(element, wx, wy)
        self.structure.add_load(load)

    def add_element_point_load(self, element_id, px=0.0, py=0.0, mz=0.0, x=0.0):
        """
        Add a point load acting on an element at a specific distance from start node (local axes).

        Parameters
        ----------
        element_id : int or str
            ID of the element.
        px : float, optional
            Axial point force. Defaults to 0.0.
        py : float, optional
            Transverse point force. Defaults to 0.0.
        mz : float, optional
            Point moment. Defaults to 0.0.
        x : float, optional
            Distance from the start node along the element length. Defaults to 0.0.
        """
        element = self.structure.elements[element_id]
        load = ElementPointLoad(element, px, py, mz, x)
        self.structure.add_load(load)

    def solve(self):
        """
        Solve the linear structure.

        Returns
        -------
        disp : numpy.ndarray
            Global displacement vector.
        reactions : numpy.ndarray
            Global reaction forces vector at fixed degrees of freedom.
        """
        self.structure.solve()
        return self.structure.disp, self.structure.reactions
