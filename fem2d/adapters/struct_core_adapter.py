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

from typing import Optional, Union, TYPE_CHECKING

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
    validate_structural_model,
)

from fem2d.elements.beam import BeamElement as _FemBeamElement
from fem2d.elements.spring import SpringElement as _FemSpringElement
from fem2d.elements.truss import TrussElement as _FemTrussElement
from fem2d.loads import DistributedLoad, PointLoad
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
    load_case_id: Optional[Union[int, str]] = None,
    validate: bool = True,
    load_combination_id: Optional[str] = None,
    mass_source: Optional[Union[Any, str]] = None,
    include_self_weight: bool = True,
    unit_system: Optional[str] = None,
    gravity_accel: Optional[float] = None,
    include_shear_deformation: Optional[bool] = None,
) -> Structure:
    """
    Convert a struct_core ``StructuralModel`` (or ``Project``) into a fully
    populated ``fem2d.Structure`` ready for analysis.

    Parameters
    ----------
    model_or_project : StructuralModel or Project
        The input schema instance.
    load_case_id : int or str, optional
        The specific load case ID to apply. If ``None`` and load cases exist,
        the first load case is applied by default.
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
    unit_system : str, optional
        'si' or 'us'. If model_or_project is a Project, inferred from `project.metadata.units`.
    gravity_accel : float, optional
        Gravity acceleration (defaults to 9.80665 for SI, 386.08858 for US).
    include_shear_deformation : bool, optional
        Whether to consider shear deformation (Timoshenko beam theory) in beam elements.
        If None (default), reads `include_shear_deformation` from each `struct_core.BeamElement`.
        If True or False, overrides all beam elements in the model.

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

    structure = Structure()

    # 1. Convert Materials
    materials_map = {}
    for mat_schema in model.materials:
        if isinstance(mat_schema, _ScElasticMaterial):
            mat_obj = _FemElasticMaterial(E=mat_schema.E, rho=mat_schema.rho)
        elif isinstance(mat_schema, _ScBilinearMaterial):
            mat_obj = _FemBilinearMaterial(
                E=mat_schema.E, Et=mat_schema.Et, fy=mat_schema.fy, rho=mat_schema.rho
            )
        else:
            mat_obj = _FemElasticMaterial(
                E=getattr(mat_schema, "E", 200e9),
                rho=getattr(mat_schema, "rho", 0.0),
            )
        materials_map[mat_schema.id] = mat_obj

    # 2. Convert Sections
    sections_map = {}
    for sec_schema in model.sections:
        if isinstance(sec_schema, _ScRectangularSection):
            sec_obj = _FemSection.from_rectangle(width=sec_schema.b, depth=sec_schema.h)
        elif isinstance(sec_schema, _ScCircularSection):
            sec_obj = _FemSection.from_circle(diameter=sec_schema.d)
        elif isinstance(sec_schema, _ScGeneralSection):
            as_val = getattr(sec_schema, "As2", getattr(sec_schema, "as2", getattr(sec_schema, "Ay", None)))
            sec_obj = _FemSection(area=sec_schema.A, moi=sec_schema.I3, As=as_val)
        else:
            sec_obj = _FemSection(
                area=getattr(sec_schema, "A", 1.0),
                moi=getattr(sec_schema, "I3", 1.0),
                As=getattr(sec_schema, "As2", None),
            )
        sections_map[sec_schema.id] = sec_obj

    # 3. Convert Nodes
    for n_schema in model.nodes:
        node_obj = _FemNode(nid=n_schema.id, x=n_schema.x, y=n_schema.y)
        if n_schema.mass > 0.0 or n_schema.inertia > 0.0:
            node_obj.set_mass(mass=n_schema.mass, inertia=n_schema.inertia)
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
            mat = materials_map[elem_schema.material_id]
            sec = sections_map[elem_schema.section_id]
            elem_obj = _FemTrussElement(
                eid=elem_schema.id,
                node_i=node_i,
                node_j=node_j,
                material=mat,
                section=sec,
            )
        elif elem_type == "beam":
            mat = materials_map[elem_schema.material_id]
            sec = sections_map[elem_schema.section_id]
            elem_shear = (
                include_shear_deformation
                if include_shear_deformation is not None
                else getattr(elem_schema, "include_shear_deformation", False)
            )
            elem_obj = _FemBeamElement(
                eid=elem_schema.id,
                node_i=node_i,
                node_j=node_j,
                material=mat,
                section=sec,
                geometric_nonlinear=elem_schema.is_nonlinear,
                include_shear_deformation=elem_shear,
            )
        elif elem_type == "spring":
            elem_obj = _FemSpringElement(
                eid=elem_schema.id,
                node_i=node_i,
                node_j=node_j,
                kx=elem_schema.kx,
                ky=elem_schema.ky,
                kz=elem_schema.kz,
            )
        else:
            raise ValueError(f"Unsupported element type: '{elem_type}'")

        structure.add_element(elem_obj)

    # 6. Apply Loads
    # Determine unit system and gravity acceleration for self-weight
    if unit_system is None:
        if isinstance(model_or_project, Project) and getattr(model_or_project, "metadata", None) and getattr(model_or_project.metadata, "units", None):
            units = model_or_project.metadata.units
            if getattr(units, "is_us_customary", False):
                unit_system = "us"
            elif getattr(units, "is_metric", False):
                unit_system = "si"
            else:
                length_u = str(getattr(units, "length", "m")).lower()
                force_u = str(getattr(units, "force", "kn")).lower()
                if length_u in ("in", "ft", "inch", "feet") or force_u in ("kip", "kips", "lb", "lbf", "psi", "ksi"):
                    unit_system = "us"
                else:
                    unit_system = "si"
            if gravity_accel is None:
                gravity_accel = getattr(units, "gravity_acceleration", None)
        else:
            unit_system = "si"
            if hasattr(model, "materials"):
                for m in model.materials:
                    uw = getattr(m, "unit_weight", None)
                    if uw is not None and 0 < uw < 1.0:
                        unit_system = "us"
                        break

    if gravity_accel is None:
        gravity_accel = 386.08858 if str(unit_system).lower() == "us" else 9.80665

    # 5b. Evaluate Mass Source if specified
    ms_target = mass_source if mass_source is not None else getattr(model, "mass_source", None)
    if ms_target is not None:
        if _calc_mass_source is not None:
            ms_res = _calc_mass_source(
                model=model,
                mass_source=ms_target,
                gravity_accel=gravity_accel,
                unit_system=unit_system,
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
                unit_system=unit_system,
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
                    c = getattr(elem_obj, "cos", 1.0)
                    s = getattr(elem_obj, "sin", 0.0)
                    wx_loc = -w_mag * factor * s
                    wy_loc = -w_mag * factor * c
                    structure.add_load(DistributedLoad(element=elem_obj, wx=wx_loc, wy=wy_loc))

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
                continue
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
                    structure.add_load(
                        DistributedLoad(
                            element=elem_obj,
                            wx=dl_schema.wx * factor,
                            wy=dl_schema.wy * factor,
                        )
                    )
            if getattr(lc, "include_self_weight", False):
                _apply_case_self_weight_2d(lc.id, factor=factor)
    else:
        cases_to_apply = []
        if model.load_cases:
            if load_case_id is not None and str(load_case_id).upper() not in ("ALL", "TOTAL", "ENVELOPE"):
                cases_to_apply = [lc for lc in model.load_cases if str(lc.id).upper() == str(load_case_id).upper()]
                if not cases_to_apply:
                    cases_to_apply = [lc for lc in model.load_cases if str(getattr(lc, "name", "")).upper() == str(load_case_id).upper()]
                if not cases_to_apply:
                    active_cases = [
                        lc for lc in model.load_cases
                        if lc.point_loads
                        or getattr(lc, "element_loads", [])
                        or getattr(lc, "distributed_loads", [])
                        or (include_self_weight and getattr(lc, "include_self_weight", False))
                    ]
                    cases_to_apply = active_cases if active_cases else [model.load_cases[0]]
            else:
                cases_to_apply = model.load_cases

        for target_load_case in cases_to_apply:
            for pl_schema in target_load_case.point_loads:
                if pl_schema.node_id in structure.nodes:
                    node_obj = structure.nodes[pl_schema.node_id]
                    load_obj = PointLoad(
                        node=node_obj, fx=pl_schema.fx, fy=pl_schema.fy, mz=pl_schema.mz
                    )
                    structure.add_load(load_obj)

            loads_list = list(getattr(target_load_case, "element_loads", []) or getattr(target_load_case, "distributed_loads", []))
            for dl_schema in loads_list:
                if dl_schema.element_id in structure.elements:
                    elem_obj = structure.elements[dl_schema.element_id]
                    dist_load_obj = DistributedLoad(
                        element=elem_obj,
                        wx=float(getattr(dl_schema, "wx", 0.0) or 0.0),
                        wy=float(getattr(dl_schema, "wy", 0.0) or 0.0),
                    )
                    structure.add_load(dist_load_obj)

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


def _make_sc_distributed_load(*, element_id, wx: float, wy: float):
    """Build a ``struct_core.DistributedLoad`` (lazy-imported)."""
    from struct_core import DistributedLoad as ScDistributedLoad

    return ScDistributedLoad(element_id=element_id, wx=wx, wy=wy, is_local=False)


def model_to_core(obj, *, validate: bool = True) -> StructuralModel:
    """
    Convert a ``fem2d.Structure`` (or ``SimpleFrame`` / ``Results`` wrapper)
    into a ``struct_core.StructuralModel``.

    Mapping
    -------
    - ``fem2d.Node`` → ``struct_core.Node`` (z=0 in the 2D world; mass and
      inertia preserved when set).
    - ``fem2d.ElasticMaterial`` → ``ElasticMaterial``;
      ``fem2d.BilinearMaterial`` → ``BilinearMaterial``. Other material
      classes are coerced to ``ElasticMaterial``.
    - ``fem2d.Section`` → ``struct_core.GeneralSection(A=A, Iz=I)`` so the
      area and second moment of area round-trip exactly.
    - ``BeamElement`` / ``BeamWithHingesElement`` / ``*NL`` → ``BeamElement``.
    - ``TrussElement`` / ``TrussElementNL`` → ``TrussElement``.
    - ``SpringElement`` → ``SpringElement``.
    - Node-level loads are gathered into one ``LoadCase`` named ``"LC1"``
      and a ``LinearStaticAnalysis`` named ``"Static"`` is added.

    Parameters
    ----------
    obj : fem2d.Structure, fem2d.SimpleFrame, or fem2d.Results
        Source fem2d model. ``SimpleFrame`` and ``Results`` are unwrapped
        automatically.
    validate : bool, optional
        Run ``struct_core.validate_structural_model`` on the produced model.
        Defaults to ``True``.

    Returns
    -------
    struct_core.StructuralModel
        A validated Pydantic model representing the fem2d input.

    Raises
    ------
    ImportError
        If ``struct_core`` is not installed. Install with
        ``pip install fem2d[ecosystem]``.
    TypeError
        If ``obj`` is not a fem2d ``Structure``, ``SimpleFrame``, or
        ``Results``.
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
        )
    except ImportError as exc:  # pragma: no cover - import guard
        raise ImportError(
            "struct_core is required for model_to_core(). "
            "Install with `pip install fem2d[ecosystem]` to enable schema export."
        ) from exc

    structure = _resolve_structure(obj)

    # ---- Materials (deduplicated by properties, not by Python object id) ---
    # Two different fem2d.ElasticMaterial instances with identical E/rho
    # collapse to a single struct_core material. This makes the M1/M2/...
    # assignment stable across repeated exports and means the schema
    # doesn't bloat when SimpleFrame.add_frame/E callers re-construct
    # the same material on every call.
    def _material_key(fem_mat) -> tuple:
        """Hashable key for material deduplication."""
        if isinstance(fem_mat, _FemBilinearMaterial):
            return (
                "bilinear",
                round(float(fem_mat.E), 12),
                round(float(getattr(fem_mat, "Et", fem_mat.E * 0.01)), 12),
                round(float(getattr(fem_mat, "fy", float("nan"))), 12),
                round(float(getattr(fem_mat, "rho", 0.0)), 12),
            )
        # Default: treat as elastic.
        return (
            "elastic",
            round(float(fem_mat.E), 12),
            round(float(getattr(fem_mat, "rho", 0.0)), 12),
        )

    def _section_key(fem_sec) -> tuple:
        """Hashable key for section deduplication."""
        A = float(getattr(fem_sec, "A", getattr(fem_sec, "area", 0.0)))
        I = float(getattr(fem_sec, "I", getattr(fem_sec, "moi", 0.0)))
        As = getattr(fem_sec, "As", None)
        return (round(A, 12), round(I, 12), round(float(As), 12) if As is not None else None)

    sc_materials: list = []
    materials_map: dict = {}  # material_key -> sc_id
    mat_counter = 0

    for elem in structure.elements.values():
        fem_mat = elem.material
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
                    Et=float(getattr(fem_mat, "Et", fem_mat.E * 0.01)),
                    fy=float(getattr(fem_mat, "fy", float("nan"))),
                    rho=float(getattr(fem_mat, "rho", 0.0)),
                )
            )
        else:
            sc_materials.append(
                ScElasticMaterial(
                    id=sc_id,
                    E=float(fem_mat.E),
                    rho=float(getattr(fem_mat, "rho", 0.0)),
                )
            )
        materials_map[key] = sc_id

    # ---- Sections (deduplicated by A, I, As) -------------------------------
    sc_sections: list = []
    sections_map: dict = {}  # section_key -> sc_id
    sec_counter = 0

    for elem in structure.elements.values():
        fem_sec = elem.section
        key = _section_key(fem_sec)
        if key in sections_map:
            continue
        sec_counter += 1
        sc_id = f"S{sec_counter}"
        A = float(getattr(fem_sec, "A", getattr(fem_sec, "area", 0.0)))
        I = float(getattr(fem_sec, "I", getattr(fem_sec, "moi", 0.0)))
        As = getattr(fem_sec, "As", None)
        as_float = float(As) if As is not None else None
        sc_sections.append(ScGeneralSection(id=sc_id, A=A, I3=I, I2=I, J=2 * I if I > 0 else 1.0, As2=as_float, As3=as_float))
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
                mass=float(getattr(node, "mass", 0.0)),
                inertia=float(getattr(node, "inertia", 0.0)),
            )
        )

    # ---- Elements ---------------------------------------------------------
    sc_elements = []
    for eid, elem in structure.elements.items():
        mat_id = materials_map[_material_key(elem.material)]
        sec_id = sections_map[_section_key(elem.section)]
        cls_name = type(elem).__name__

        if cls_name in (
            "BeamElement",
            "BeamWithHingesElement",
            "BeamElementNL",
            "BeamElementMaterialNL",
        ):
            rel_i = list(getattr(elem, "releases_i", [False, False, getattr(elem, "hinge_i", False)]))
            rel_j = list(getattr(elem, "releases_j", [False, False, getattr(elem, "hinge_j", False)]))
            sc_elements.append(
                ScBeamElement(
                    id=eid,
                    start_node=elem.node_i.id,
                    end_node=elem.node_j.id,
                    material_id=mat_id,
                    section_id=sec_id,
                    release_i=rel_i,
                    release_j=rel_j,
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
                    material_id=mat_id,
                    section_id=sec_id,
                )
            )
        elif cls_name == "SpringElement":
            sc_elements.append(
                ScSpringElement(
                    id=eid,
                    start_node=elem.node_i.id,
                    end_node=elem.node_j.id,
                    kx=float(getattr(elem, "kx", 0.0)),
                    ky=float(getattr(elem, "ky", 0.0)),
                    kz=float(getattr(elem, "kz", 0.0)),
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
        if isinstance(sup, (list, tuple)) and len(sup) >= 3:
            ux_fixed, uy_fixed, rz_fixed = bool(sup[0]), bool(sup[1]), bool(sup[2])
        else:  # pragma: no cover
            ux_fixed = uy_fixed = rz_fixed = False
        if ux_fixed or uy_fixed or rz_fixed:
            sc_supports.append(
                ScSupport(node_id=nid, ux=ux_fixed, uy=uy_fixed, rz=rz_fixed)
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
                    wx=float(getattr(load_obj, "wx", 0.0)),
                    wy=float(getattr(load_obj, "wy", 0.0)),
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
                id="LC1",
                name="LC1",
                point_loads=sc_point_loads,
                element_loads=sc_element_loads,
            )
        )
        sc_analysis_cases.append(
            LinearStaticAnalysis(id="Static", load_case_id="LC1")
        )

    model = StructuralModel(
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

    return model


def _build_node_displacement(ux: float, uy: float, theta: float):
    """Lazy import wrapper for ``struct_core.NodeDisplacement``."""
    from struct_core import NodeDisplacement

    return NodeDisplacement(ux=ux, uy=uy, uz=0.0, rx=0.0, ry=0.0, rz=theta)


def _build_node_reaction(fx: float, fy: float, mz: float):
    """Lazy import wrapper for ``struct_core.NodeReaction``."""
    from struct_core import NodeReaction

    return NodeReaction(fx=fx, fy=fy, fz=0.0, mx=0.0, my=0.0, mz=mz)


def _build_section_force(station: float, axial: float, shear: float, moment: float):
    """Lazy import wrapper for ``struct_core.BeamSectionForce``."""
    from struct_core import BeamSectionForce

    return BeamSectionForce(
        station=station,
        P=axial,
        V2=shear,
        V3=0.0,
        T=0.0,
        M2=0.0,
        M3=moment,
    )


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
    force_df = results.element_forces()
    element_results = {}
    for _, row in force_df.iterrows():
        eid = row["element"]
        # Two stations per element: end i (station 0.0) and end j (station 1.0).
        # fem2d reports local-coordinates forces; sign convention is preserved.
        s_i = _build_section_force(
            station=0.0,
            axial=float(row["fx_i"]),
            shear=float(row["fy_i"]),
            moment=float(row["m_i"]),
        )
        s_j = _build_section_force(
            station=1.0,
            axial=float(row["fx_j"]),
            shear=float(row["fy_j"]),
            moment=float(row["m_j"]),
        )
        stations = [s_i, s_j]
        # Envelopes: pick the station with the largest absolute magnitude.
        max_axial = s_i.P if abs(s_i.P) >= abs(s_j.P) else s_j.P
        max_shear = s_i.V2 if abs(s_i.V2) >= abs(s_j.V2) else s_j.V2
        max_moment = s_i.M3 if abs(s_i.M3) >= abs(s_j.M3) else s_j.M3

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
    from struct_core import AnalysisResult, ElementResult, NodeResult

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

    # Compute envelope across all combination results
    envelope_elements = {}
    first_res = next(iter(comb_results.values()))

    for eid in first_res.element_results.keys():
        all_elem_res = [cr.element_results[eid] for cr in comb_results.values() if eid in cr.element_results]
        if not all_elem_res:
            continue

        max_p = max(all_elem_res, key=lambda er: abs(er.max_axial or 0.0)).max_axial
        max_m = max(all_elem_res, key=lambda er: abs(er.max_moment or 0.0)).max_moment
        max_v = max(all_elem_res, key=lambda er: abs(er.max_shear or 0.0)).max_shear

        gov_elem = max(all_elem_res, key=lambda er: abs(er.max_moment or 0.0))
        envelope_elements[eid] = ElementResult(
            element_id=eid,
            forces=gov_elem.forces,
            max_axial=max_p,
            max_moment=max_m,
            max_shear=max_v,
        )

    envelope_nodes = {}
    for nid in first_res.node_results.keys():
        all_node_res = [cr.node_results[nid] for cr in comb_results.values() if nid in cr.node_results]
        gov_node = max(
            all_node_res,
            key=lambda nr: (nr.displacement.ux**2 + nr.displacement.uy**2) if nr.displacement else 0.0,
        )
        envelope_nodes[nid] = gov_node

    envelope_analysis_res = AnalysisResult(
        analysis_case_id=envelope_case_id,
        node_results=envelope_nodes,
        element_results=envelope_elements,
    )

    if isinstance(model_or_project, Project):
        for cid, cres in comb_results.items():
            model_or_project.analysis_results[cid] = cres
        model_or_project.analysis_results[envelope_case_id] = envelope_analysis_res

    return comb_results, envelope_analysis_res