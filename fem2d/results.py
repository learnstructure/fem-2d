import numpy as np
import pandas as pd


class Results:
    def __init__(self, obj):
        """
        Parameters
        ----------
        obj : Structure or SimpleFrame
            The analysed structure (or the SimpleFrame wrapper).
        """
        # If obj is a SimpleFrame, extract the underlying Structure
        if hasattr(obj, "structure") and hasattr(obj.structure, "nodes"):
            self.structure = obj.structure
        else:
            self.structure = obj

    def node_displacements(self):
        """Return DataFrame with node displacements."""
        data = []
        for node in self.structure.nodes.values():
            disp = self.structure.disp[node.dofs]
            data.append(
                {"node": node.id, "ux": disp[0], "uy": disp[1], "theta": disp[2]}
            )
        return pd.DataFrame(data)

    def reactions(self):
        """Return DataFrame with reaction forces."""
        if self.structure.reactions is None:
            raise ValueError("Reactions not computed. Run analysis first.")
        data = []
        for node in self.structure.nodes.values():
            reac = self.structure.reactions[node.dofs]
            data.append({"node": node.id, "Fx": reac[0], "Fy": reac[1], "M": reac[2]})
        return pd.DataFrame(data)

    def element_forces(self):
        """Return DataFrame with element end forces in local coordinates."""
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
