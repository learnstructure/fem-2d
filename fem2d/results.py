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

    def _ensure_analysis_ready(self):
        """Ensure the structure has numbered DOFs and assembled matrices for reporting."""
        if not hasattr(self.structure, "neq") or getattr(self.structure, "neq", None) is None:
            self.structure.number_dofs()
        if getattr(self.structure, "K", None) is None:
            self.structure.assemble_stiffness()
        if getattr(self.structure, "F", None) is None:
            self.structure.assemble_loads()
        if not hasattr(self.structure, "free_dofs") or getattr(self.structure, "free_dofs", None) is None:
            self.structure.apply_boundary_conditions()
        if getattr(self.structure, "disp", None) is None:
            self.structure.solve()
        if getattr(self.structure, "reactions", None) is None:
            self.structure.solve()

    @staticmethod
    def _format_array(array):
        """Format a numpy array as a readable multiline string."""
        if array is None:
            return "Not computed"
        return np.array2string(np.asarray(array), precision=6, suppress_small=False)

    def create_report(self, print_report=True):
        """
        Create a text report summarizing the assembled matrices and results.

        The report is intended for classroom verification and includes element
        stiffness contributions in local and global coordinates, the assembled
        global system matrices, the reduced matrices for free DOFs, and the
        computed displacements and reactions.

        Parameters
        ----------
        print_report : bool, optional
            If True, print the report to stdout. Defaults to True.

        Returns
        -------
        str
            The formatted report as a multiline string.
        """
        self._ensure_analysis_ready()

        lines = []
        lines.append("Structural analysis report")
        lines.append("=" * 80)
        lines.append(f"Nodes: {len(self.structure.nodes)}")
        lines.append(f"Elements: {len(self.structure.elements)}")
        lines.append(f"Total DOFs: {self.structure.neq}")
        lines.append(f"Fixed DOFs: {self.structure.fixed_dofs}")
        lines.append(f"Free DOFs: {self.structure.free_dofs}")
        lines.append("")

        lines.append("Node support conditions")
        lines.append("-" * 80)
        for node in self.structure.nodes.values():
            lines.append(
                f"Node {node.id}: support={node.support}, dofs={node.dofs}, load={node.load}"
            )
        lines.append("")

        lines.append("Element reports")
        lines.append("-" * 80)
        for el in self.structure.elements.values():
            lines.append(f"Element {el.id}: nodes=({el.node_i.id}, {el.node_j.id})")
            lines.append(f"  DOFs: {el.node_i.dofs + el.node_j.dofs}")
            lines.append("  Element stiffness matrix in local coordinates")
            lines.append(self._format_array(el.local_stiffness()))
            lines.append("  Transformation matrix")
            lines.append(self._format_array(el.transformation_matrix()))
            lines.append("  Element stiffness matrix in global coordinates")
            lines.append(self._format_array(el.global_stiffness()))
            if hasattr(el, "equivalent_nodal_loads"):
                lines.append("  Equivalent nodal loads")
                lines.append(self._format_array(el.equivalent_nodal_loads()))
            lines.append("")

        lines.append("Structure stiffness matrix")
        lines.append("-" * 80)
        lines.append(self._format_array(self.structure.K))
        lines.append("")

        lines.append("Structure load vector")
        lines.append("-" * 80)
        lines.append(self._format_array(self.structure.F))
        lines.append("")

        lines.append("Reduced structure stiffness matrix")
        lines.append("-" * 80)
        k_ff = self.structure.K[np.ix_(self.structure.free_dofs, self.structure.free_dofs)]
        lines.append(self._format_array(k_ff))
        lines.append("")

        lines.append("Reduced structure load vector")
        lines.append("-" * 80)
        f_f = self.structure.F[self.structure.free_dofs]
        lines.append(self._format_array(f_f))
        lines.append("")

        lines.append("Displacements")
        lines.append("-" * 80)
        disp_df = self.node_displacements()
        lines.append(disp_df.to_string(index=False))
        lines.append("")

        lines.append("Reactions")
        lines.append("-" * 80)
        reaction_df = self.reactions()
        lines.append(reaction_df.to_string(index=False))
        lines.append("")

        lines.append("Element end forces")
        lines.append("-" * 80)
        lines.append(self.element_forces().to_string(index=False))

        report = "\n".join(lines)
        if print_report:
            print(report)
        return report
