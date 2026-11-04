"""
Adapter module bridging fem2d and struct_core.

Provides three conversion functions that keep fem2d and ``struct_core`` in
sync without one library hard-depending on the other:

- :func:`model_from_core` — convert a ``struct_core.StructuralModel`` (or
  ``struct_core.Project``) into a fully populated ``fem2d.Structure`` ready
  for analysis.
- :func:`model_to_core` — convert a ``fem2d.Structure`` (or a
  ``fem2d.SimpleFrame``) into a ``struct_core.StructuralModel``. The
  resulting model is suitable for re-import, JSON serialization, or feeding
  into downstream ``struct_core`` consumers.
- :func:`result_to_core` — convert the output of a fem2d analysis (a
  ``fem2d.Results`` post-processor, a ``fem2d.Structure``, or a
  ``fem2d.SimpleFrame``) into a ``struct_core.AnalysisResult`` carrying nodal
  displacements, support reactions, and element internal forces.

Naming convention
-----------------
Every function name follows ``<artefact>_<direction>_<target>``:

- ``model_from_core``  — model in, model out (struct_core → fem2d)
- ``model_to_core``    — model in, model out (fem2d → struct_core)
- ``result_to_core``   — results in, results out (fem2d → struct_core)

There is intentionally **no** ``result_from_core`` — round-tripping an
``AnalysisResult`` back into a ``fem2d.Results`` is out of scope.

Optional dependency
-------------------
``struct_core`` is imported lazily so fem2d remains usable without it
installed. The ``[ecosystem]`` extra pulls it in:

    pip install fem2d[ecosystem]
"""

from __future__ import annotations

import copy
from typing import Any, Dict, Optional, Union, TYPE_CHECKING

# --- struct_core -> fem2d direction: types are needed at module scope for
# the public function signatures. Imported eagerly so introspection tools and
# type checkers see them. Runtime cost is zero if user never calls the fn.
from struct_core import (
    BilinearMaterial as _ScBilinearMaterial,
    CircularSection as _ScCircularSection,
    ElasticMaterial as _ScElasticMaterial,
    GeneralSection as _ScGeneralSection,
    Project,
    RectangularSection as _ScRectangularSection,
    StructuralModel,
    Units,
    parse_units,
    validate_structural_model,
)

from fem2d.elements.beam import BeamElement as _FemBeamElement
from fem2d.elements.spring import SpringElement as _FemSpringElement
from fem2d.elements.truss import TrussElement as _FemTrussElement
from fem2d.loads import DistributedLoad, ElementPointLoad, PointLoad
from fem2d.materials.bilinear import BilinearMaterial as _FemBilinearMaterial
from fem2d.materials.elastic import ElasticMaterial as _FemElasticMaterial
from fem2d.nodes import Node as _FemNode
from fem2d.sections.section import Section as _FemSection
from fem2d.structure import Structure

try:
    from struct_load.engine.self_weight import calculate_self_weight as _calc_self_weight
    from struct_load.engine.mass_source import calculate_mass_source as _calc_mass_source
    _STRUCT_LOAD_AVAILABLE = True
except ImportError:
    _calc_self_weight = None
    _calc_mass_source = None
    _STRUCT_LOAD_AVAILABLE = False

if TYPE_CHECKING:  # pragma: no cover - typing only
    from struct_core import AnalysisResult


__all__ = ["model_from_core", "model_to_core", "result_to_core", "solve_combinations_and_envelope"]


# ============================================================================
# model_from_core: struct_core -> fem2d.Structure
# ============================================================================


def model_from_core(
    model_or_project: Union[StructuralModel, Project],
    load_case_id: Optional[Union[int, str, List[Union[int, str]], Tuple[Union[int, str], ...]]] = None,
    validate: bool = True,
    load_combination_id: Optional[str] = None,
    mass_source: Optional[Union[Any, str]] = None,
    include_self_weight: bool = True,
    units: Optional[Union[str, Units, dict]] = None,
    gravity_accel: Optional[float] = None,
    include_shear_deformation: Optional[bool] = None,
    flexural_stiffness_modifiers: Optional[Dict[Any, float]] = None,
) -> Structure:
    """
    Convert a struct_core ``StructuralModel`` (or ``Project``) into a fully
    populated ``fem2d.Structure`` ready for analysis.

    Parameters
    ----------
    model_or_project : StructuralModel or Project
        The input schema instance.
    load_case_id : int, str, or list of (int or str), optional
        Target load case ID or list of IDs. If None (default), all load cases
        are transferred and summed. Pass 'ALL' for all cases, or specific IDs like 'DEAD'.
    validate : bool, optional
        Whether to run structural validation rules before building.
        Defaults to ``True``.
    load_combination_id : str, optional
        The specific load combination ID to apply.
    mass_source : MassSource, dict, or str, optional
        Inertial mass source configuration. If provided, evaluates masses from dead/live load
        patterns and element self-masses and assigns them to nodes and frame elements.
    include_self_weight : bool, optional
        Whether to automatically compute and include element self-weight for any
        load case where `include_self_weight=True`. Defaults to True.
        Uses `struct_load.engine.self_weight.calculate_self_weight`.
    units : str, Units, or dict, optional
        Units specification (e.g. 'kN, m' or 'kip, in'). If None, inferred automatically
        from model_or_project.units or model_or_project.metadata.units. Defaults to 'kN, m'.
    gravity_accel : float, optional
        Gravity acceleration. Defaults to `active_units.gravity_acceleration`.
    include_shear_deformation : bool, optional
        Whether to consider shear deformation (Timoshenko beam theory) in beam elements.
        If None (default), reads `include_shear_deformation` from each `struct_core.BeamElement`.
        If True or False, overrides all beam elements in the model.
    flexural_stiffness_modifiers : dict, optional
        Multiplier on the moment of inertia per element ID (e.g. ACI 318 Table 6.6.3.1.1:
        0.70 for columns, 0.35 for beams). Area and mass are not modified.

    Returns
    -------
    Structure
        An initialized ``fem2d.Structure`` instance ready for analysis.
    """
    if isinstance(model_or_project, Project):
        model = model_or_project.model
    else:
        model = model_or_project

    if validate:
        validate_structural_model(model, raise_error=True)

    active_units = parse_units(
        units
        or getattr(getattr(model_or_project, "metadata", None), "units", None)
        or getattr(model_or_project, "units", None)
        or "kN, m"
    )

    structure = Structure(units=active_units)

    if gravity_accel is None:
        gravity_accel = active_units.gravity_acceleration

    mass_factor = active_units.consistent_mass_factor

    def _convert_rho(mat_schema: Any) -> float:
        # Dynamic mass density in consistent units: [F] * [T]^2 / [L]^4
        uw = getattr(mat_schema, "unit_weight", None)
        if uw is not None and float(uw) > 0.0:
            return float(uw) / gravity_accel
        raw_rho = getattr(mat_schema, "rho", 0.0)
        val = float(raw_rho or 0.0)
        if val <= 0.0:
            return 0.0
        return val * mass_factor

    def _convert_mass(raw_m: float) -> float:
        val = float(raw_m or 0.0)
        if val <= 0.0:
            return 0.0
        return val * mass_factor

    # 1. Convert Materials
    materials_map = {}
    for mat_schema in model.materials:
        rho_val = _convert_rho(mat_schema)
        uw_val = getattr(mat_schema, "unit_weight", None)
        if isinstance(mat_schema, _ScElasticMaterial):  # also concrete and rebar, elastic for analysis
            mat_obj = _FemElasticMaterial(E=mat_schema.E, rho=rho_val, nu=mat_schema.nu, unit_weight=uw_val)
        elif isinstance(mat_schema, _ScBilinearMaterial):
            mat_obj = _FemBilinearMaterial(E=mat_schema.E, sigma_y=mat_schema.fy, E_t=mat_schema.Et, rho=rho_val)
        else:
            raise TypeError(f"Material {mat_schema.id!r}: fem2d has no model for {type(mat_schema).__name__}")
        materials_map[mat_schema.id] = mat_obj

    # 2. Convert Sections (each with its material)
    sections_map = {}
    for sec_schema in model.sections:
        if sec_schema.material_id not in materials_map:
            raise ValueError(f"Section {sec_schema.id!r} references missing material {sec_schema.material_id!r}")
        material = materials_map[sec_schema.material_id]
        if isinstance(sec_schema, _ScRectangularSection):
            sec_obj = _FemSection.from_rectangle(width=sec_schema.b, depth=sec_schema.h, material=material)
        elif isinstance(sec_schema, _ScCircularSection):
            sec_obj = _FemSection.from_circle(diameter=sec_schema.d, material=material)
        elif isinstance(sec_schema, _ScGeneralSection):
            needs_shear = (
                include_shear_deformation is True
                or any(
                    getattr(e, "include_shear_deformation", False)
                    for e in model.elements
                    if getattr(e, "section_id", None) == sec_schema.id
                )
            )
            if needs_shear and sec_schema.As2 is None:
                raise ValueError(
                    f"Section {sec_schema.id!r} is a GeneralSection without shear area 'As2', but "
                    f"shear deformation was requested. Specify 'As2' on the section or disable shear deformation."
                )
            sec_obj = _FemSection(area=sec_schema.A, moi=sec_schema.I3, As=sec_schema.As2, material=material)
        else:
            raise TypeError(f"Section {sec_schema.id!r}: fem2d has no model for {type(sec_schema).__name__}")
        sections_map[sec_schema.id] = sec_obj

    # 3. Convert Nodes
    for n_schema in model.nodes:
        node_obj = _FemNode(nid=n_schema.id, x=n_schema.x, y=n_schema.y)
        if n_schema.mass > 0.0 or n_schema.inertia > 0.0:
            node_obj.set_mass(
                mass=_convert_mass(n_schema.mass),
                inertia=_convert_mass(n_schema.inertia),
            )
        structure.add_node(node_obj)

    # 4. Apply Supports
    for sup_schema in model.supports:
        if sup_schema.node_id in structure.nodes:
            node_obj = structure.nodes[sup_schema.node_id]
            node_obj.set_support(
                ux_fixed=sup_schema.ux,
                uy_fixed=sup_schema.uy,
                rz_fixed=sup_schema.rz,
            )
            kx = float(getattr(sup_schema, "kx", 0.0) or 0.0)
            ky = float(getattr(sup_schema, "ky", 0.0) or 0.0)
            krz = float(getattr(sup_schema, "krz", 0.0) or 0.0)
            kz = float(getattr(sup_schema, "kz", 0.0) or 0.0)
            krx = float(getattr(sup_schema, "krx", 0.0) or 0.0)
            kry = float(getattr(sup_schema, "kry", 0.0) or 0.0)
            if kz > 0.0 or krx > 0.0 or kry > 0.0:
                raise ValueError(
                    f"Support at node {sup_schema.node_id!r} specifies 3D spring components (kz={kz}, krx={krx}, kry={kry}) "
                    f"in a 2D model."
                )
            if kx > 0.0 or ky > 0.0 or krz > 0.0:
                node_obj.set_springs(kx=kx, ky=ky, krz=krz)

    # 5. Convert Elements
    for elem_schema in model.elements:
        if elem_schema.start_node not in structure.nodes or elem_schema.end_node not in structure.nodes:
            continue
        node_i = structure.nodes[elem_schema.start_node]
        node_j = structure.nodes[elem_schema.end_node]
        elem_type = elem_schema.type

        dx = float(node_j.x - node_i.x)
        dy = float(node_j.y - node_i.y)
        if (dx**2 + dy**2) < 1e-12 and elem_type != "spring":
            continue

        if elem_type == "truss":
            sec = sections_map[elem_schema.section_id]
            elem_obj = _FemTrussElement(
                eid=elem_schema.id,
                node_i=node_i,
                node_j=node_j,
                section=sec,
            )
        elif elem_type == "beam":
            sec = sections_map[elem_schema.section_id]
            elem_mods = getattr(elem_schema, "stiffness_modifiers", None)
            if elem_mods or (flexural_stiffness_modifiers and elem_schema.id in flexural_stiffness_modifiers):
                sec = copy.copy(sec)
                if elem_mods:
                    if "A" in elem_mods:
                        sec.A *= float(elem_mods["A"])
                    if "I" in elem_mods:
                        sec.I *= float(elem_mods["I"])
                    elif "I3" in elem_mods:
                        sec.I *= float(elem_mods["I3"])
                    if "As" in elem_mods and hasattr(sec, "As"):
                        sec.As *= float(elem_mods["As"])
                    elif "As2" in elem_mods and hasattr(sec, "As"):
                        sec.As *= float(elem_mods["As2"])
                if flexural_stiffness_modifiers and elem_schema.id in flexural_stiffness_modifiers:
                    sec.I *= float(flexural_stiffness_modifiers[elem_schema.id])
            elem_shear = (
                include_shear_deformation
                if include_shear_deformation is not None
                else getattr(elem_schema, "include_shear_deformation", False)
            )
            rel_i = getattr(elem_schema, "release_i", None)
            rel_j = getattr(elem_schema, "release_j", None)
            releases_i_2d = None
            releases_j_2d = None
            # Member local axes: [u1 (axial P), u2 (shear V2), r3 (moment M3)] mapped from 6-DOF [u1, u2, u3, r1, r2, r3]
            if rel_i is not None and any(rel_i):
                releases_i_2d = [bool(rel_i[0]), bool(rel_i[1]), bool(rel_i[5] if len(rel_i) >= 6 else rel_i[2])]
            if rel_j is not None and any(rel_j):
                releases_j_2d = [bool(rel_j[0]), bool(rel_j[1]), bool(rel_j[5] if len(rel_j) >= 6 else rel_j[2])]

            elem_obj = _FemBeamElement(
                eid=elem_schema.id,
                node_i=node_i,
                node_j=node_j,
                section=sec,
                geometric_nonlinear=elem_schema.is_nonlinear,
                include_shear_deformation=elem_shear,
                releases_i=releases_i_2d,
                releases_j=releases_j_2d,
            )
        elif elem_type == "spring":
            nonzero = [(name, val) for name, val in [("kx", elem_schema.kx), ("ky", elem_schema.ky), ("kz", elem_schema.kz)] if val > 0]
            if not nonzero:
                raise ValueError(f"Spring element '{elem_schema.id}' has zero stiffness (kx=ky=kz=0).")
            if len(nonzero) > 1:
                names = [n for n, _ in nonzero]
                raise ValueError(
                    f"fem2d SpringElement '{elem_schema.id}' specifies multiple directional stiffnesses ({', '.join(names)}). "
                    f"fem2d supports 1D axial spring elements along the member axis; define separate 1D spring elements for each direction."
                )
            stiffness = nonzero[0][1]
            elem_obj = _FemSpringElement(
                eid=elem_schema.id,
                node_i=node_i,
                node_j=node_j,
                stiffness=stiffness,
            )
        else:
            raise ValueError(f"Unsupported element type: '{elem_type}'")

        structure.add_element(elem_obj)

    # 6. Apply Loads (units and gravity_accel already resolved)

    # 5b. Evaluate Mass Source if specified
    ms_target = mass_source if mass_source is not None else getattr(model, "mass_source", None)
    if ms_target is not None:
        if _calc_mass_source is not None:
            ms_res = _calc_mass_source(
                model=model,
                mass_source=ms_target,
                gravity_accel=gravity_accel,
                units=active_units,
            )
            for nid, m_vec in ms_res.nodal_masses.items():
                if nid in structure.nodes:
                    m_val = m_vec[0] if isinstance(m_vec, (list, tuple)) else float(m_vec)
                    structure.nodes[nid].set_mass(mass=m_val)

            for eid, m_extra in ms_res.element_masses.items():
                if eid in structure.elements:
                    structure.elements[eid].mass += float(m_extra)

            structure.mass_source_result = ms_res
        else:
            import warnings
            warnings.warn(
                "A mass_source was specified, but 'struct_load' is not installed to evaluate it. "
                "Install struct_load to enable automated mass source calculation.",
                UserWarning,
                stacklevel=2,
            )

    sw_cache = {}

    def _to_fem_distributed_load(elem_obj, dl_schema, factor: float = 1.0) -> DistributedLoad:
        """Map a struct_core DistributedLoad onto fem2d's local (axial, transverse) components.

        struct_core stores global components (wx, wy) unless ``is_local`` is set, in which case
        (w1, w2) are already local. fem2d.DistributedLoad is always local, so global loads are
        rotated: w_axial = wx*c + wy*s, w_transverse = -wx*s + wy*c.
        """
        if getattr(dl_schema, "is_local", False):
            w_axial = float(dl_schema.w1 or 0.0)
            w_trans = float(dl_schema.w2 or 0.0)
            w1_e = getattr(dl_schema, "w1_end", None)
            w2_e = getattr(dl_schema, "w2_end", None)
            w_axial_end = float(w_axial if w1_e is None else w1_e)
            w_trans_end = float(w_trans if w2_e is None else w2_e)
        else:
            gx = float(dl_schema.wx or 0.0)
            gy = float(dl_schema.wy or 0.0)
            gx_e = getattr(dl_schema, "wx_end", None)
            gy_e = getattr(dl_schema, "wy_end", None)
            gx_end = float(gx if gx_e is None else gx_e)
            gy_end = float(gy if gy_e is None else gy_e)
            c = getattr(elem_obj, "cos", 1.0)
            s = getattr(elem_obj, "sin", 0.0)
            w_axial = gx * c + gy * s
            w_trans = -gx * s + gy * c
            w_axial_end = gx_end * c + gy_end * s
            w_trans_end = -gx_end * s + gy_end * c

        st = float(getattr(dl_schema, "start", 0.0) or 0.0)
        en = float(getattr(dl_schema, "end", 1.0) or 1.0)
        return DistributedLoad(
            element=elem_obj,
            wx=w_axial * factor,
            wy=w_trans * factor,
            wx_end=w_axial_end * factor,
            wy_end=w_trans_end * factor,
            start=st,
            end=en,
        )

    def _apply_case_self_weight_2d(case_id: str, factor: float = 1.0):
        if not include_self_weight:
            return
        if _calc_self_weight is None:
            import warnings
            warnings.warn(
                f"Load case '{case_id}' has include_self_weight=True, but 'struct_load' is not installed. "
                "Self-weight was not added to the analysis model. Install struct_load to enable automatic self-weight.",
                UserWarning,
                stacklevel=2,
            )
            return

        if case_id not in sw_cache:
            sw_cache[case_id] = _calc_self_weight(
                model=model,
                load_case_id=case_id,
                gravity_accel=gravity_accel,
                units=active_units,
                direction="-y",
            )
        sw_case = sw_cache[case_id]

        for dl in getattr(sw_case, "element_loads", []):
            eid = dl.element_id
            if eid in structure.elements:
                elem_obj = structure.elements[eid]
                w_mag = abs(float(getattr(dl, "wy", 0.0) or getattr(dl, "wx", 0.0) or 0.0))
                if w_mag <= 0.0:
                    continue
                if isinstance(elem_obj, _FemTrussElement):
                    total_w = w_mag * factor * elem_obj.length
                    structure.add_load(PointLoad(elem_obj.node_i, fy=-0.5 * total_w))
                    structure.add_load(PointLoad(elem_obj.node_j, fy=-0.5 * total_w))
                else:
                    structure.add_load(_to_fem_distributed_load(elem_obj, dl, factor))

    target_comb = None
    if load_combination_id is not None:
        target_comb = next(
            (c for c in getattr(model, "load_combinations", []) if c.id == load_combination_id), None
        )
    elif load_case_id is not None and getattr(model, "load_combinations", None):
        target_comb = next(
            (c for c in model.load_combinations if c.id == load_case_id), None
        )

    if target_comb is not None:
        load_cases_map = {lc.id: lc for lc in getattr(model, "load_cases", [])}
        for factor_obj in target_comb.factors:
            lc = load_cases_map.get(factor_obj.load_case_id)
            if not lc:
                raise KeyError(
                    f"Load combination '{target_comb.id}' references unknown load case '{factor_obj.load_case_id}'."
                )
            factor = factor_obj.factor
            for pl_schema in lc.point_loads:
                if pl_schema.node_id in structure.nodes:
                    node_obj = structure.nodes[pl_schema.node_id]
                    structure.add_load(
                        PointLoad(
                            node=node_obj,
                            fx=pl_schema.fx * factor,
                            fy=pl_schema.fy * factor,
                            mz=pl_schema.mz * factor,
                        )
                    )
            for dl_schema in lc.element_loads:
                if dl_schema.element_id in structure.elements:
                    elem_obj = structure.elements[dl_schema.element_id]
                    structure.add_load(_to_fem_distributed_load(elem_obj, dl_schema, factor))
            for epl_schema in getattr(lc, "element_point_loads", []):
                if epl_schema.element_id in structure.elements:
                    elem_obj = structure.elements[epl_schema.element_id]
                    pos = float(getattr(epl_schema, "position", 0.0)) * elem_obj.length
                    if getattr(epl_schema, "is_local", False):
                        px = float(getattr(epl_schema, "fx", 0.0)) * factor
                        py = float(getattr(epl_schema, "fy", 0.0)) * factor
                        mz = float(getattr(epl_schema, "mz", 0.0)) * factor
                    else:
                        gx = float(getattr(epl_schema, "fx", 0.0)) * factor
                        gy = float(getattr(epl_schema, "fy", 0.0)) * factor
                        mz = float(getattr(epl_schema, "mz", 0.0)) * factor
                        c = getattr(elem_obj, "cos", 1.0)
                        s = getattr(elem_obj, "sin", 0.0)
                        px = gx * c + gy * s
                        py = -gx * s + gy * c
                    structure.add_load(ElementPointLoad(element=elem_obj, px=px, py=py, mz=mz, x=pos))
            if getattr(lc, "include_self_weight", False):
                _apply_case_self_weight_2d(lc.id, factor=factor)
    else:
        cases_to_apply = []
        if model.load_cases:
            if load_case_id is None:
                cases_to_apply = model.load_cases
            elif isinstance(load_case_id, (list, tuple, set)):
                targets = {str(cid).upper() for cid in load_case_id}
                cases_to_apply = [
                    lc for lc in model.load_cases
                    if str(lc.id).upper() in targets or str(getattr(lc, "name", "")).upper() in targets
                ]
                if not cases_to_apply:
                    raise ValueError(
                        f"None of the requested load cases {load_case_id} found in model. Available: {[lc.id for lc in model.load_cases]}"
                    )
            else:
                cid_str = str(load_case_id).upper()
                if cid_str in ("ALL", "TOTAL", "ENVELOPE"):
                    cases_to_apply = model.load_cases
                else:
                    cases_to_apply = [
                        lc for lc in model.load_cases
                        if str(lc.id).upper() == cid_str or str(getattr(lc, "name", "")).upper() == cid_str
                    ]
                    if not cases_to_apply:
                        raise ValueError(
                            f"Requested load case '{load_case_id}' not found in model load cases: {[lc.id for lc in model.load_cases]}"
                        )

        for target_load_case in cases_to_apply:
            for pl_schema in target_load_case.point_loads:
                if pl_schema.node_id in structure.nodes:
                    node_obj = structure.nodes[pl_schema.node_id]
                    load_obj = PointLoad(
                        node=node_obj, fx=pl_schema.fx, fy=pl_schema.fy, mz=pl_schema.mz
                    )
                    structure.add_load(load_obj)

            for dl_schema in target_load_case.element_loads:
                if dl_schema.element_id in structure.elements:
                    elem_obj = structure.elements[dl_schema.element_id]
                    structure.add_load(_to_fem_distributed_load(elem_obj, dl_schema))
            for epl_schema in getattr(target_load_case, "element_point_loads", []):
                if epl_schema.element_id in structure.elements:
                    elem_obj = structure.elements[epl_schema.element_id]
                    pos = float(getattr(epl_schema, "position", 0.0)) * elem_obj.length
                    if getattr(epl_schema, "is_local", False):
                        px = float(getattr(epl_schema, "fx", 0.0))
                        py = float(getattr(epl_schema, "fy", 0.0))
                        mz = float(getattr(epl_schema, "mz", 0.0))
                    else:
                        gx = float(getattr(epl_schema, "fx", 0.0))
                        gy = float(getattr(epl_schema, "fy", 0.0))
                        mz = float(getattr(epl_schema, "mz", 0.0))
                        c = getattr(elem_obj, "cos", 1.0)
                        s = getattr(elem_obj, "sin", 0.0)
                        px = gx * c + gy * s
                        py = -gx * s + gy * c
                    structure.add_load(ElementPointLoad(element=elem_obj, px=px, py=py, mz=mz, x=pos))

            if getattr(target_load_case, "include_self_weight", False):
                _apply_case_self_weight_2d(target_load_case.id, factor=1.0)

    return structure


# ============================================================================
# model_to_core: fem2d.Structure -> struct_core.StructuralModel
# ============================================================================


def _resolve_structure(obj) -> Structure:
    """
    Accept any of ``Structure``, ``SimpleFrame``, or ``Results`` and return the
    underlying ``fem2d.Structure``. Mirrors the duck-typing used in
    ``fem2d.Results.__init__``.
    """
    # Results / SimpleFrame: expose ``.structure``
    if hasattr(obj, "structure") and isinstance(
        getattr(obj, "structure", None), Structure
    ):
        return obj.structure
    # Bare Structure: already the right type
    if isinstance(obj, Structure):
        return obj
    raise TypeError(
        f"expected a fem2d.Structure (or SimpleFrame / Results); "
        f"got {type(obj).__name__}."
    )


def _make_sc_distributed_load(*, element_id, w_axial: float, w_transverse: float):
    """Build a local-axis ``struct_core.DistributedLoad`` from fem2d's local (axial, transverse) load."""
    from struct_core import DistributedLoad as ScDistributedLoad

    return ScDistributedLoad(element_id=element_id, w1=w_axial, w2=w_transverse, is_local=True)


def model_to_core(
    obj,
    *,
    units: Optional[Union[str, Any, dict]] = None,
    default_load_case_id: str = "LIVE",
    validate: bool = True,
) -> Project:
    """
    Convert a ``fem2d.Structure`` (or ``SimpleFrame`` / ``Results`` wrapper)
    into a ``struct_core.Project``.

    Parameters
    ----------
    obj : fem2d.Structure, fem2d.SimpleFrame, or fem2d.Results
        Source fem2d model. ``SimpleFrame`` and ``Results`` are unwrapped
        automatically.
    units : str, Units, or dict, optional
        Units specification (e.g. 'kN, m'). If None, inferred from structure.units.
    default_load_case_id : str, optional
        Target load case ID to store existing structure loads into. Defaults to 'LIVE'.
    validate : bool, optional
        Run ``struct_core.validate_structural_model`` on the produced model.
        Defaults to ``True``.

    Returns
    -------
    struct_core.Project
        A project instance wrapping the structural model and units metadata.
    """
    try:
        from struct_core import (
            BeamElement as ScBeamElement,
            BilinearMaterial as ScBilinearMaterial,
            ElasticMaterial as ScElasticMaterial,
            GeneralSection as ScGeneralSection,
            LinearStaticAnalysis,
            LoadCase as ScLoadCase,
            Node as ScNode,
            PointLoad as ScPointLoad,
            SpringElement as ScSpringElement,
            Support as ScSupport,
            TrussElement as ScTrussElement,
            Units,
            parse_units,
        )
    except ImportError as exc:  # pragma: no cover - import guard
        raise ImportError(
            "struct_core is required for model_to_core(). "
            "Install with `pip install fem2d[ecosystem]` to enable schema export."
        ) from exc

    structure = _resolve_structure(obj)
    # Masses and mass densities: consistent units in fem2d, the model's mass unit in struct_core
    from struct_core import Metadata
    units_spec = units or getattr(structure, "units", None)
    target_units = (parse_units(units_spec) if isinstance(units_spec, (str, dict)) else units_spec) or Metadata().units
    mass_factor = target_units.consistent_mass_factor

    # ---- Materials (deduplicated by properties, not by Python object id) ---
    def _material_key(fem_mat) -> tuple:
        """Hashable key for material deduplication."""
        if isinstance(fem_mat, _FemBilinearMaterial):
            return ("bilinear", round(float(fem_mat.E), 12), round(float(fem_mat.E_t), 12),
                    round(float(fem_mat.sigma_y), 12), round(float(fem_mat.rho), 12))
        if not isinstance(fem_mat, _FemElasticMaterial):
            raise TypeError(f"struct_core has no schema for the fem2d material {type(fem_mat).__name__}")
        return ("elastic", round(float(fem_mat.E), 12), round(float(fem_mat.nu), 12), round(float(fem_mat.rho), 12),
                fem_mat.unit_weight)

    def _section_key(fem_sec) -> tuple:
        """Hashable key for section deduplication."""
        A = float(getattr(fem_sec, "A", getattr(fem_sec, "area", 0.0)))
        I = float(getattr(fem_sec, "I", getattr(fem_sec, "moi", 0.0)))
        As = getattr(fem_sec, "As", None)
        return (round(A, 12), round(I, 12), round(float(As), 12) if As is not None else None)

    # Springs have no physical section or material; they export their stiffness only
    framed = [elem for elem in structure.elements.values() if type(elem).__name__ != "SpringElement"]

    sc_materials: list = []
    materials_map: dict = {}  # material_key -> sc_id
    mat_counter = 0

    for elem in framed:
        fem_mat = elem.section.material
        key = _material_key(fem_mat)
        if key in materials_map:
            continue
        mat_counter += 1
        sc_id = f"M{mat_counter}"
        if isinstance(fem_mat, _FemBilinearMaterial):
            sc_materials.append(
                ScBilinearMaterial(
                    id=sc_id,
                    E=float(fem_mat.E),
                    Et=float(fem_mat.E_t),
                    fy=float(fem_mat.sigma_y),
                    rho=float(fem_mat.rho) / mass_factor,
                )
            )
        else:
            sc_materials.append(
                ScElasticMaterial(
                    id=sc_id,
                    E=float(fem_mat.E),
                    nu=float(fem_mat.nu),
                    rho=float(fem_mat.rho) / mass_factor,
                    unit_weight=fem_mat.unit_weight,
                )
            )
        materials_map[key] = sc_id

    # ---- Sections (deduplicated by A, I, As and material) -----------------
    sc_sections: list = []
    sections_map: dict = {}  # (section_key, material id) -> sc_id
    sec_counter = 0

    for elem in framed:
        fem_sec = elem.section
        mat_id = materials_map[_material_key(elem.section.material)]
        key = (_section_key(fem_sec), mat_id)
        if key in sections_map:
            continue
        sec_counter += 1
        sc_id = f"S{sec_counter}"
        A = float(getattr(fem_sec, "A", getattr(fem_sec, "area", 0.0)))
        I = float(getattr(fem_sec, "I", getattr(fem_sec, "moi", 0.0)))
        As = getattr(fem_sec, "As", None)
        as_float = float(As) if As is not None else None
        sc_sections.append(ScGeneralSection(id=sc_id, material_id=mat_id, A=A, I3=I, I2=I, J=2 * I if I > 0 else 1.0,
                                            As2=as_float, As3=as_float))
        sections_map[key] = sc_id

    # ---- Nodes ------------------------------------------------------------
    sc_nodes = []
    for nid, node in structure.nodes.items():
        sc_nodes.append(
            ScNode(
                id=nid,
                x=float(node.x),
                y=float(node.y),
                z=0.0,
                mass=float(getattr(node, "mass", 0.0)) / mass_factor,
                inertia=float(getattr(node, "inertia", 0.0)) / mass_factor,
            )
        )

    # ---- Elements ---------------------------------------------------------
    sc_elements = []
    for eid, elem in structure.elements.items():
        cls_name = type(elem).__name__
        if cls_name != "SpringElement":
            sec_id = sections_map[(_section_key(elem.section), materials_map[_material_key(elem.section.material)])]

        if cls_name in (
            "BeamElement",
            "BeamWithHingesElement",
            "BeamElementNL",
            "BeamElementMaterialNL",
        ):
            rel_i = list(getattr(elem, "releases_i", [False, False, getattr(elem, "hinge_i", False)]))
            rel_j = list(getattr(elem, "releases_j", [False, False, getattr(elem, "hinge_j", False)]))
            rel_i_6d = [bool(rel_i[0]), bool(rel_i[1]), False, False, False, bool(rel_i[2])] if len(rel_i) >= 3 else [False]*6
            rel_j_6d = [bool(rel_j[0]), bool(rel_j[1]), False, False, False, bool(rel_j[2])] if len(rel_j) >= 3 else [False]*6
            sc_elements.append(
                ScBeamElement(
                    id=eid,
                    start_node=elem.node_i.id,
                    end_node=elem.node_j.id,
                    section_id=sec_id,
                    release_i=rel_i_6d,
                    release_j=rel_j_6d,
                    is_nonlinear=bool(getattr(elem, "geometric_nonlinear", False)),
                    include_shear_deformation=bool(getattr(elem, "include_shear_deformation", False)),
                )
            )
        elif cls_name in ("TrussElement", "TrussElementNL"):
            sc_elements.append(
                ScTrussElement(
                    id=eid,
                    start_node=elem.node_i.id,
                    end_node=elem.node_j.id,
                    section_id=sec_id,
                )
            )
        elif cls_name == "SpringElement":
            sc_elements.append(
                ScSpringElement(
                    id=eid,
                    start_node=elem.node_i.id,
                    end_node=elem.node_j.id,
                    kx=float(getattr(elem, "k", 0.0)),
                )
            )
        else:  # pragma: no cover - defensive default
            raise TypeError(
                f"model_to_core(): unsupported element class '{cls_name}' "
                f"(id={eid})."
            )

    # ---- Supports ---------------------------------------------------------
    sc_supports = []
    for nid, node in structure.nodes.items():
        sup = getattr(node, "support", None) or [False, False, False]
        springs = getattr(node, "springs", None) or [0.0, 0.0, 0.0]
        if isinstance(sup, (list, tuple)) and len(sup) >= 3:
            ux_fixed, uy_fixed, rz_fixed = bool(sup[0]), bool(sup[1]), bool(sup[2])
        else:  # pragma: no cover
            ux_fixed = uy_fixed = rz_fixed = False
        kx = float(springs[0]) if len(springs) > 0 else 0.0
        ky = float(springs[1]) if len(springs) > 1 else 0.0
        krz = float(springs[2]) if len(springs) > 2 else 0.0
        if ux_fixed or uy_fixed or rz_fixed or kx > 0.0 or ky > 0.0 or krz > 0.0:
            sc_supports.append(
                ScSupport(node_id=nid, ux=ux_fixed, uy=uy_fixed, rz=rz_fixed, kx=kx, ky=ky, krz=krz)
            )

    # ---- Loads ------------------------------------------------------------
    sc_point_loads = []
    for nid, node in structure.nodes.items():
        load = getattr(node, "load", None)
        if load is None:
            continue
        fx = float(load[0]) if len(load) >= 1 else 0.0
        fy = float(load[1]) if len(load) >= 2 else 0.0
        mz = float(load[2]) if len(load) >= 3 else 0.0
        if fx == 0.0 and fy == 0.0 and mz == 0.0:
            continue
        sc_point_loads.append(ScPointLoad(node_id=nid, fx=fx, fy=fy, mz=mz))

    sc_element_loads = []
    for load_obj in getattr(structure, "loads", []) or []:
        if isinstance(load_obj, DistributedLoad):
            sc_element_loads.append(
                _make_sc_distributed_load(
                    element_id=load_obj.element.id,
                    w_axial=float(load_obj.wx),
                    w_transverse=float(load_obj.wy),
                )
            )

        elif isinstance(load_obj, PointLoad):
            nid = load_obj.node.id if hasattr(load_obj, "node") else None
            if nid is not None:
                fx = float(getattr(load_obj, "fx", 0.0))
                fy = float(getattr(load_obj, "fy", 0.0))
                mz = float(getattr(load_obj, "mz", 0.0))
                if fx == 0.0 and fy == 0.0 and mz == 0.0:
                    continue
                sc_point_loads.append(
                    ScPointLoad(node_id=nid, fx=fx, fy=fy, mz=mz)
                )

    sc_load_cases: list = []
    sc_analysis_cases: list = []
    if sc_point_loads or sc_element_loads:
        sc_load_cases.append(
            ScLoadCase(
                id=default_load_case_id,
                name="Live Load" if default_load_case_id == "LIVE" else default_load_case_id,
                point_loads=sc_point_loads,
                element_loads=sc_element_loads,
            )
        )
        sc_analysis_cases.append(
            LinearStaticAnalysis(id="Static", load_case_id=default_load_case_id)
        )

    model = StructuralModel(
        dimension="2D",
        vertical_axis="Y",
        nodes=sc_nodes,
        elements=sc_elements,
        materials=sc_materials,
        sections=sc_sections,
        supports=sc_supports,
        load_cases=sc_load_cases,
        analysis_cases=sc_analysis_cases,
    )

    if validate:
        validate_structural_model(model, raise_error=True)

    from struct_core import Metadata, Project
    return Project(model=model, metadata=Metadata(units=target_units))


def _build_node_displacement(ux: float, uy: float, theta: float):
    """Lazy import wrapper for ``struct_core.NodeDisplacement``."""
    from struct_core import NodeDisplacement

    return NodeDisplacement(ux=ux, uy=uy, uz=0.0, rx=0.0, ry=0.0, rz=theta)


def _build_node_reaction(fx: float, fy: float, mz: float):
    """Lazy import wrapper for ``struct_core.NodeReaction``."""
    from struct_core import NodeReaction

    return NodeReaction(fx=fx, fy=fy, fz=0.0, mx=0.0, my=0.0, mz=mz)


def _member_stations(length: float, f_i: tuple, loads: list) -> list:
    """Internal forces at STANDARD_STATIONS from the start-node end forces and the member's loads.

    ``f_i`` = (fx, fy, mz) are the local forces the start node exerts on the member. Equilibrium of
    the segment [0, x] with local member loads (uniform/trapezoidal/partial ``wx``, ``wy`` and concentrated
    ``px``, ``py``, ``mz`` at distance a) gives the actions on its cut face:
        P(x) = -(fx + Wx(x) + sum px)      V(x) = -(fy + Wy(x) + sum py)
        M(x) = -mz + fy x + Mw(x) + sum[py (x - a) - mz_p]
    where Wy(x) is the load resultant over [0, x] and Mw(x) its moment about the cut.
    At x = L these equal the end forces at node j.
    """
    from struct_core import BeamSectionForce
    from struct_core.analysis_result import STANDARD_STATIONS
    from fem2d.loads import ElementPointLoad, TriangularLoad

    fx, fy, mz = f_i
    dist_loads = []
    points = []
    for load in loads:
        if isinstance(load, DistributedLoad):
            dist_loads.append(load)
        elif isinstance(load, TriangularLoad):
            class _TriWrapper:
                start = 0.0
                end = 1.0
                wx = load.wx1
                wx_end = load.wx2
                wy = load.wy1
                wy_end = load.wy2
            dist_loads.append(_TriWrapper())
        elif isinstance(load, ElementPointLoad):
            points.append(load)
        else:
            raise TypeError(f"Unsupported member load {type(load).__name__} for internal force stations.")

    stations = []
    for s in STANDARD_STATIONS:
        x = s * length
        rx_total = 0.0
        ry_total = 0.0
        my_total = 0.0

        for dl in dist_loads:
            x1 = dl.start * length
            x2 = dl.end * length
            if x > x1 and x2 > x1:
                x_cut = min(x, x2)
                delta = x_cut - x1
                l_span = x2 - x1
                mx_slope = (dl.wx_end - dl.wx) / l_span
                my_slope = (dl.wy_end - dl.wy) / l_span

                rx_total += dl.wx * delta + 0.5 * mx_slope * delta**2
                ry = dl.wy * delta + 0.5 * my_slope * delta**2
                ry_total += ry
                my_total += (x - x1) * ry - (0.5 * dl.wy * delta**2 + (1.0 / 3.0) * my_slope * delta**3)

        acting = [p for p in points if p.x < x or (s == 1.0 and p.x <= length)]
        stations.append(BeamSectionForce(
            station=s,
            P=-(fx + rx_total + sum(p.px for p in acting)),
            V2=-(fy + ry_total + sum(p.py for p in acting)),
            M3=-mz + fy * x + my_total + sum(p.py * (x - p.x) - p.mz for p in acting),
        ))
    return stations


def result_to_core(
    obj,
    analysis_case_id: Union[int, str] = "Static",
    include_reactions: bool = True,
):
    """
    Convert the output of a fem2d analysis into a ``struct_core.AnalysisResult``.

    Parameters
    ----------
    obj : fem2d.Results, fem2d.Structure, or fem2d.SimpleFrame
        The analyzed structure. If a ``Results`` post-processor or a
        ``SimpleFrame`` is supplied, the underlying ``Structure`` is
        extracted automatically.
    analysis_case_id : int or str, optional
        Identifier written to ``AnalysisResult.analysis_case_id``.
        Defaults to ``"Static"``.
    include_reactions : bool, optional
        Whether to populate ``NodeResult.reaction`` for supported nodes.
        Defaults to ``True``. Has no effect when ``structure.reactions`` is
        ``None`` (i.e. the analysis has not produced reactions yet).

    Returns
    -------
    struct_core.AnalysisResult
        A Pydantic model containing nodal displacements, support reactions,
        and element internal forces with envelopes.

    Raises
    ------
    ImportError
        If ``struct_core`` is not installed. Install with
        ``pip install fem2d[ecosystem]``.
    """
    # Lazy import so fem2d remains usable without struct_core installed.
    try:
        from struct_core import AnalysisResult, ElementResult, NodeResult
    except ImportError as exc:  # pragma: no cover - import guard
        raise ImportError(
            "struct_core is required for result_to_core(). "
            "Install with `pip install fem2d[ecosystem]` to enable schema export."
        ) from exc

    structure = _resolve_structure(obj)

    # Lazy import Results locally to avoid circular import at module load.
    from fem2d.results import Results

    results = Results(structure)

    # --- Node results -------------------------------------------------------
    disp_df = results.node_displacements().set_index("node")
    has_reactions = (
        include_reactions and getattr(structure, "reactions", None) is not None
    )
    if has_reactions:
        reac_df = results.reactions().set_index("node")

    node_results = {}
    for nid in structure.nodes.keys():
        if nid not in disp_df.index:
            continue
        row = disp_df.loc[nid]
        displacement = _build_node_displacement(
            ux=float(row["ux"]),
            uy=float(row["uy"]),
            theta=float(row["theta"]),
        )
        reaction = None
        if has_reactions and nid in reac_df.index:
            rrow = reac_df.loc[nid]
            reaction = _build_node_reaction(
                fx=float(rrow["Fx"]),
                fy=float(rrow["Fy"]),
                mz=float(rrow["M"]),
            )
        node_results[nid] = NodeResult(
            node_id=nid,
            displacement=displacement,
            reaction=reaction,
        )

    # --- Element results ----------------------------------------------------
    loads_by_element: dict = {}
    for load in getattr(structure, "loads", []) or []:
        el = getattr(load, "element", None)
        if el is not None:
            loads_by_element.setdefault(el.id, []).append(load)

    force_df = results.element_forces()
    element_results = {}
    for _, row in force_df.iterrows():
        eid = row["element"]
        stations = _member_stations(
            structure.elements[eid].length,
            (float(row["fx_i"]), float(row["fy_i"]), float(row["m_i"])),
            loads_by_element.get(eid, []),
        )
        # Envelopes: the station value with the largest absolute magnitude.
        max_axial = max((s.P for s in stations), key=abs)
        max_shear = max((s.V2 for s in stations), key=abs)
        max_moment = max((s.M3 for s in stations), key=abs)

        element_results[eid] = ElementResult(
            element_id=eid,
            forces=stations,
            max_axial=max_axial,
            max_moment=max_moment,
            max_shear=max_shear,
        )

    return AnalysisResult(
        analysis_case_id=analysis_case_id,
        node_results=node_results,
        element_results=element_results,
    )


def solve_combinations_and_envelope(
    model_or_project: Union[StructuralModel, Project],
    combination_ids: Optional[list] = None,
    envelope_case_id: str = "Envelope",
):
    """
    Solve all load combinations for a StructuralModel or Project and compile the governing force envelope.

    Parameters
    ----------
    model_or_project : StructuralModel or Project
        Input schema model or project.
    combination_ids : list of str, optional
        Specific load combination IDs to solve. If None, all model.load_combinations are solved.
    envelope_case_id : str
        ID for the resulting enveloped AnalysisResult. Defaults to 'Envelope'.

    Returns
    -------
    Tuple[Dict[str, AnalysisResult], AnalysisResult]
        (combination_results_dict, envelope_result)
    """
    from struct_core.analysis_result import envelope_results

    if isinstance(model_or_project, Project):
        model = model_or_project.model
    else:
        model = model_or_project

    all_combs = getattr(model, "load_combinations", [])
    if combination_ids is not None:
        target_combs = [c for c in all_combs if c.id in combination_ids]
    else:
        target_combs = all_combs

    if not target_combs:
        struct = model_from_core(model)
        struct.solve()
        res = result_to_core(struct, analysis_case_id=envelope_case_id)
        if isinstance(model_or_project, Project):
            model_or_project.analysis_results[envelope_case_id] = res
        return {envelope_case_id: res}, res

    comb_results = {}
    for comb in target_combs:
        struct = model_from_core(model, load_combination_id=comb.id)
        struct.solve()
        res = result_to_core(struct, analysis_case_id=comb.id)
        comb_results[comb.id] = res

    # Governing envelope; per-combination station forces stay available to design
    envelope_analysis_res = envelope_results(list(comb_results.values()), output_case_id=envelope_case_id)

    if isinstance(model_or_project, Project):
        for cid, cres in comb_results.items():
            model_or_project.analysis_results[cid] = cres
        model_or_project.analysis_results[envelope_case_id] = envelope_analysis_res

    return comb_results, envelope_analysis_res