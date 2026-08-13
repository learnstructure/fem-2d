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

if TYPE_CHECKING:  # pragma: no cover - typing only
    from struct_core import AnalysisResult


__all__ = ["model_from_core", "model_to_core", "result_to_core"]


# ============================================================================
# model_from_core: struct_core -> fem2d.Structure
# ============================================================================


def model_from_core(
    model_or_project: Union[StructuralModel, Project],
    load_case_id: Optional[Union[int, str]] = None,
    validate: bool = True,
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
            sec_obj = _FemSection(area=sec_schema.A, moi=sec_schema.Iz)
        else:
            sec_obj = _FemSection(
                area=getattr(sec_schema, "A", 1.0),
                moi=getattr(sec_schema, "Iz", 1.0),
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
        node_i = structure.nodes[elem_schema.start_node]
        node_j = structure.nodes[elem_schema.end_node]
        elem_type = elem_schema.type

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
            elem_obj = _FemBeamElement(
                eid=elem_schema.id,
                node_i=node_i,
                node_j=node_j,
                material=mat,
                section=sec,
                geometric_nonlinear=elem_schema.is_nonlinear,
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
    target_load_case = None
    if model.load_cases:
        if load_case_id is not None:
            target_load_case = next(
                (lc for lc in model.load_cases if lc.id == load_case_id), None
            )
        else:
            target_load_case = model.load_cases[0]

    if target_load_case is not None:
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
                dist_load_obj = DistributedLoad(
                    element=elem_obj, wx=dl_schema.wx, wy=dl_schema.wy
                )
                structure.add_load(dist_load_obj)

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
        return (round(A, 12), round(I, 12))

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

    # ---- Sections (deduplicated by A, I) ----------------------------------
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
        sc_sections.append(ScGeneralSection(id=sc_id, A=A, Iz=I))
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
            sc_elements.append(
                ScBeamElement(
                    id=eid,
                    start_node=elem.node_i.id,
                    end_node=elem.node_j.id,
                    material_id=mat_id,
                    section_id=sec_id,
                    is_nonlinear=bool(getattr(elem, "geometric_nonlinear", False)),
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
        Vy=shear,
        Vz=0.0,
        T=0.0,
        My=0.0,
        Mz=moment,
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
        max_shear = s_i.Vy if abs(s_i.Vy) >= abs(s_j.Vy) else s_j.Vy
        max_moment = s_i.Mz if abs(s_i.Mz) >= abs(s_j.Mz) else s_j.Mz

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