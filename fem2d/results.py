"""
Results module for extracting and formatting analysis results into pandas DataFrames.
"""

import numpy as np
import pandas as pd


class Results:
    """
    Helper class to post-process and retrieve structural analysis results
    such as node displacements, reaction forces, and member end forces.

    Attributes
    ----------
    structure : Structure
        The analyzed Structure object from which results are retrieved.
    """

    def __init__(self, obj):
        """
        Initialize the Results post-processor.

        Parameters
        ----------
        obj : Structure or SimpleFrame
            The analyzed structure or the SimpleFrame wrapper containing the structure.
        """
        # If obj is a SimpleFrame, extract the underlying Structure
        if hasattr(obj, "structure") and hasattr(obj.structure, "nodes"):
            self.structure = obj.structure
        else:
            self.structure = obj

    def node_displacements(self):
        """
        Return displacements of all nodes in a DataFrame format.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing columns:
            - 'node': Node identifier
            - 'ux': Displacement in global x-direction
            - 'uy': Displacement in global y-direction
            - 'theta': Rotation about global z-axis
        """
        data = []
        for node in self.structure.nodes.values():
            disp = self.structure.disp[node.dofs]
            data.append(
                {"node": node.id, "ux": disp[0], "uy": disp[1], "theta": disp[2]}
            )
        return pd.DataFrame(data)

    def reactions(self):
        """
        Return support reaction forces in a DataFrame format.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing columns:
            - 'node': Node identifier
            - 'Fx': Reaction force in global x-direction
            - 'Fy': Reaction force in global y-direction
            - 'M': Reaction moment about global z-axis

        Raises
        ------
        ValueError
            If reactions have not been computed (run analysis first).
        """
        if self.structure.reactions is None:
            raise ValueError("Reactions not computed. Run analysis first.")
        data = []
        for node in self.structure.nodes.values():
            reac = self.structure.reactions[node.dofs]
            data.append({"node": node.id, "Fx": reac[0], "Fy": reac[1], "M": reac[2]})
        return pd.DataFrame(data)

    def element_forces(self):
        """
        Return local end forces for all elements in a DataFrame format.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing columns:
            - 'element': Element identifier
            - 'fx_i', 'fy_i', 'm_i': Axial force, shear force, and moment at end i
            - 'fx_j', 'fy_j', 'm_j': Axial force, shear force, and moment at end j
        """
        data = []
        for el in self.structure.elements.values():
            f_local = (
                el.get_local_forces()
            )  # each element knows how to compute its own forces
            data.append(
                {
                    "element": el.id,
                    "fx_i": f_local[0],
                    "fy_i": f_local[1],
                    "m_i": f_local[2],
                    "fx_j": f_local[3],
                    "fy_j": f_local[4],
                    "m_j": f_local[5],
                }
            )
        return pd.DataFrame(data)
