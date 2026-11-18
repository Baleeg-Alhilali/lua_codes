#!/usr/bin/env python3
"""Create a compact technical handoff for the Stokes-Darcy flow model."""

import base64
import collections
import os
import struct
import textwrap
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mpl-flow-handoff")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from PIL import Image
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image as RLImage, KeepTogether, PageBreak,
    PageTemplate, Paragraph, Spacer, Table, TableStyle,
)


ROOT = "/Users/alhilaba/Documents/Promesh_swc/Model"
TMP = os.path.join(ROOT, "tmp/pdfs")
OUT = os.path.join(ROOT, "output/pdf")
PDF_PATH = os.path.join(OUT, "stokes_darcy_flow_coupling_handoff.pdf")
MEMBRANE_VTU = os.path.join(
    ROOT, "Results_final/pressure_trace_darcy_fe_ref2_membrane_t0000.vtu"
)
INTER_VTU = os.path.join(
    ROOT, "Results_final/pressure_trace_darcy_fe_ref2_inter_t0000.vtu"
)
PRESSURE_FIG = os.path.join(TMP, "darcy_pressure_ref2.png")
SCHEMATIC_FIG = os.path.join(TMP, "coupling_schematic.png")

NAVY = colors.HexColor("#17324D")
BLUE = colors.HexColor("#2F6F9F")
PALE_BLUE = colors.HexColor("#EAF2F8")
PALE_ORANGE = colors.HexColor("#FFF1DF")
ORANGE = colors.HexColor("#D97925")
RED = colors.HexColor("#B53A3A")
GREEN = colors.HexColor("#2C7A5A")
GRAY = colors.HexColor("#5B6570")
LIGHT = colors.HexColor("#F5F7F9")


VTK_TYPES = {"Float32": "f", "Float64": "d", "Int32": "i", "Int8": "b"}


def read_array(node):
    raw = base64.b64decode("".join(node.itertext()).strip())
    size = struct.unpack("<I", raw[:4])[0]
    code = VTK_TYPES[node.attrib["type"]]
    item_size = struct.calcsize(code)
    vals = struct.unpack("<" + code * (size // item_size), raw[4:4 + size])
    ncomp = int(node.attrib.get("NumberOfComponents", "1"))
    arr = np.asarray(vals)
    return arr.reshape((-1, ncomp)) if ncomp > 1 else arr


def read_vtu(path, pressure_name="darcy_pressure_kinematic"):
    root = ET.parse(path).getroot()
    arrays = list(root.iter("DataArray"))
    points = read_array(next(a for a in arrays if "Name" not in a.attrib))
    named = {a.attrib["Name"]: read_array(a) for a in arrays if "Name" in a.attrib}
    conn = named["connectivity"].astype(np.int64)
    offsets = named["offsets"].astype(np.int64)
    types = named["types"].astype(np.int64)
    pressure = named[pressure_name]
    cells = []
    start = 0
    for stop, typ in zip(offsets, types):
        cells.append((int(typ), tuple(conn[start:stop])))
        start = stop
    return points, pressure, cells


def boundary_faces(cells):
    counter = collections.Counter()
    original = {}
    for typ, v in cells:
        if typ == 10 and len(v) == 4:
            a, b, c, d = v
            faces = ((a, b, c), (a, b, d), (a, c, d), (b, c, d))
        elif typ == 12 and len(v) == 8:
            a, b, c, d, e, f, g, h = v
            faces = ((a,b,c,d),(e,f,g,h),(a,b,f,e),(b,c,g,f),(c,d,h,g),(d,a,e,h))
        elif typ == 13 and len(v) == 6:
            a, b, c, d, e, f = v
            faces = ((a,b,c),(d,e,f),(a,b,e,d),(b,c,f,e),(c,a,d,f))
        else:
            continue
        for face in faces:
            key = tuple(sorted(face))
            counter[key] += 1
            original[key] = face
    return [original[k] for k, count in counter.items() if count == 1]


def triangulate(faces):
    triangles = []
    for face in faces:
        if len(face) == 3:
            triangles.append(face)
        else:
            for i in range(1, len(face) - 1):
                triangles.append((face[0], face[i], face[i + 1]))
    return np.asarray(triangles, dtype=np.int64)


def pca_projection(points):
    centered = points - points.mean(axis=0)
    _, _, vt = np.linalg.svd(centered[::max(1, len(centered)//10000)], full_matrices=False)
    axial = vt[0]
    vertical = vt[1]
    horizontal = np.cross(axial, vertical)
    horizontal /= np.linalg.norm(horizontal)
    vertical = np.cross(horizontal, axial)
    projected = np.column_stack((centered @ axial, centered @ vertical, centered @ horizontal))
    return projected


def create_pressure_figure():
    pts, p, cells = read_vtu(MEMBRANE_VTU)
    tris = triangulate(boundary_faces(cells))
    proj = pca_projection(pts)
    # Slight oblique projection shows axial bands and the cylindrical wall.
    x = proj[:, 0] + 0.18 * proj[:, 2]
    y = proj[:, 1] + 0.32 * proj[:, 2]
    xy = np.column_stack((x, y))
    # Painter ordering for the oblique view.
    depth = proj[tris, 2].mean(axis=1)
    order = np.argsort(depth)
    tris = tris[order]
    polys = xy[tris]
    values = p[tris].mean(axis=1)

    fig, ax = plt.subplots(figsize=(9.0, 4.15), dpi=220)
    coll = PolyCollection(polys, array=values, cmap="turbo", linewidths=0.03,
                          edgecolors=(0.08, 0.12, 0.25, 0.18))
    lo, hi = np.percentile(values, [1.0, 99.0])
    coll.set_clim(lo, hi)
    ax.add_collection(coll)
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_axis_off()
    ax.set_facecolor("#515a73")
    fig.patch.set_facecolor("#515a73")
    cbar = fig.colorbar(coll, ax=ax, fraction=0.028, pad=0.015)
    cbar.set_label("Darcy kinematic pressure [mm$^2$/s$^2$]", color="white", fontsize=8)
    cbar.ax.tick_params(colors="white", labelsize=7)
    cbar.outline.set_edgecolor("white")
    ax.text(0.01, 0.97, "Membrane pressure, numRefs = 2",
            transform=ax.transAxes, va="top", color="white", fontsize=10,
            weight="bold", bbox=dict(facecolor="#263248", alpha=0.72, edgecolor="none", pad=4))
    fig.tight_layout(pad=0.4)
    fig.savefig(PRESSURE_FIG, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return float(p.min()), float(p.max())


def create_schematic():
    fig, ax = plt.subplots(figsize=(8.4, 2.55), dpi=220)
    ax.set_aspect("equal")
    ax.set_xlim(-1.7, 6.7)
    ax.set_ylim(-1.55, 1.55)
    ax.axis("off")
    # Longitudinal cutaway: lumen, membrane and interstitium.
    ax.fill_between([0, 5.7], -0.55, 0.55, color="#B9D9F2", label="Lumen: Stokes FV1")
    ax.fill_between([0, 5.7], 0.55, 0.83, color="#F5B971", label="Membrane: Darcy FE")
    ax.fill_between([0, 5.7], -0.83, -0.55, color="#F5B971")
    ax.fill_between([0, 5.7], 0.83, 1.27, color="#D9E6D1", label="Inter: Darcy FE")
    ax.fill_between([0, 5.7], -1.27, -0.83, color="#D9E6D1")
    ax.plot([0, 5.7], [0.55, 0.55], color="#B85A20", lw=1.5)
    ax.plot([0, 5.7], [-0.55, -0.55], color="#B85A20", lw=1.5)
    ax.annotate("prescribed parabolic\nvelocity", xy=(0.05, 0), xytext=(-1.55, 0),
                arrowprops=dict(arrowstyle="->", color="#17324D", lw=1.5),
                ha="center", va="center", fontsize=8, color="#17324D")
    ax.annotate("zero normal-stress\noutlet", xy=(5.65, 0), xytext=(6.45, 0),
                arrowprops=dict(arrowstyle="->", color="#17324D", lw=1.5),
                ha="center", va="center", fontsize=8, color="#17324D")
    ax.annotate("Apical coupling", xy=(2.9, 0.57), xytext=(2.9, 1.48),
                arrowprops=dict(arrowstyle="-|>", color="#B53A3A", lw=1.3),
                ha="center", fontsize=8, weight="bold", color="#B53A3A")
    ax.text(2.85, 0.02, "Stokes", ha="center", va="center", fontsize=10, weight="bold", color="#17324D")
    ax.text(2.85, 0.70, "Membrane", ha="center", va="center", fontsize=8, color="#6A3B10")
    ax.text(2.85, 1.06, "Inter", ha="center", va="center", fontsize=8, color="#315634")
    ax.text(5.7, 1.39, "$p_D=p_{outer}$ on OuterWall", ha="right", fontsize=8, color="#315634")
    ax.text(0.05, -1.48,
            "Apical laws: pressure trace equality; merged normal conservation row; BJS tangential slip; Darcy-pressure traction on Stokes",
            fontsize=7.5, color="#333333")
    fig.tight_layout(pad=0.25)
    fig.savefig(SCHEMATIC_FIG, transparent=False, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D5DCE3"))
    canvas.line(18*mm, 13*mm, A4[0]-18*mm, 13*mm)
    canvas.setFont("Helvetica", 7.2)
    canvas.setFillColor(GRAY)
    canvas.drawString(18*mm, 8.5*mm, "Request for technical guidance | UG4 Stokes-Darcy coupling")
    canvas.drawRightString(A4[0]-18*mm, 8.5*mm, f"Page {doc.page}")
    canvas.restoreState()


def table(data, widths, header=True, font=7.3):
    wrapped = []
    for row_index, row in enumerate(data):
        wrapped_row = []
        for cell in row:
            color = "#FFFFFF" if header and row_index == 0 else "#20262C"
            weight = "Helvetica-Bold" if header and row_index == 0 else "Helvetica"
            wrapped_row.append(Paragraph(
                escape(str(cell)),
                ParagraphStyle(
                    f"table-{row_index}-{len(wrapped_row)}",
                    fontName=weight, fontSize=font, leading=font + 2.0,
                    textColor=colors.HexColor(color), spaceAfter=0,
                ),
            ))
        wrapped.append(wrapped_row)
    t = Table(wrapped, colWidths=widths, repeatRows=1 if header else 0, hAlign="LEFT")
    style = [
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("FONTNAME", (0,0), (-1,-1), "Helvetica"),
        ("FONTSIZE", (0,0), (-1,-1), font),
        ("LEADING", (0,0), (-1,-1), font + 2.0),
        ("GRID", (0,0), (-1,-1), 0.35, colors.HexColor("#CBD3DA")),
        ("LEFTPADDING", (0,0), (-1,-1), 4),
        ("RIGHTPADDING", (0,0), (-1,-1), 4),
        ("TOPPADDING", (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
    ]
    if header:
        style += [
            ("BACKGROUND", (0,0), (-1,0), NAVY),
            ("TEXTCOLOR", (0,0), (-1,0), colors.white),
            ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ]
        for row in range(1, len(data)):
            if row % 2 == 0:
                style.append(("BACKGROUND", (0,row), (-1,row), LIGHT))
    t.setStyle(TableStyle(style))
    return t


def make_pdf(pmin, pmax):
    styles = getSampleStyleSheet()
    title = ParagraphStyle("TitleX", parent=styles["Title"], fontName="Helvetica-Bold",
        fontSize=19, leading=22, textColor=NAVY, spaceAfter=5)
    subtitle = ParagraphStyle("Subtitle", parent=styles["Normal"], fontSize=9.2,
        leading=12, textColor=GRAY, spaceAfter=7)
    h1 = ParagraphStyle("H1X", parent=styles["Heading1"], fontName="Helvetica-Bold",
        fontSize=12.5, leading=15, textColor=NAVY, spaceBefore=3, spaceAfter=5)
    h2 = ParagraphStyle("H2X", parent=styles["Heading2"], fontName="Helvetica-Bold",
        fontSize=9.4, leading=11, textColor=BLUE, spaceBefore=3, spaceAfter=3)
    body = ParagraphStyle("BodyX", parent=styles["BodyText"], fontSize=8.25,
        leading=10.7, textColor=colors.HexColor("#252B31"), spaceAfter=4)
    small = ParagraphStyle("SmallX", parent=body, fontSize=7.3, leading=9.2)
    bullet = ParagraphStyle("BulletX", parent=body, leftIndent=10, firstLineIndent=-6,
        bulletIndent=2, spaceAfter=2.4)
    callout = ParagraphStyle("Callout", parent=body, fontSize=8.2, leading=10.6,
        leftIndent=7, rightIndent=7, borderColor=ORANGE, borderWidth=0.8,
        borderPadding=7, backColor=PALE_ORANGE, spaceBefore=4, spaceAfter=6)
    eq = ParagraphStyle("Eq", parent=body, alignment=TA_CENTER, fontName="Helvetica-Oblique",
        fontSize=9, leading=12, textColor=NAVY, spaceBefore=2, spaceAfter=4)
    caption = ParagraphStyle("Caption", parent=small, fontName="Helvetica-Oblique",
        textColor=GRAY, alignment=TA_LEFT, spaceBefore=2, spaceAfter=5)

    doc = BaseDocTemplate(PDF_PATH, pagesize=A4, rightMargin=17*mm, leftMargin=17*mm,
                          topMargin=15*mm, bottomMargin=17*mm,
                          title="Request for technical guidance on UG4 Stokes-Darcy coupling",
                          author="Stokes-Darcy model development notes")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="main")
    doc.addPageTemplates(PageTemplate(id="standard", frames=[frame], onPage=footer))
    story = []

    def consultation_story():
        content = [
            Paragraph("Request for Technical Guidance: UG4 Stokes-Darcy Coupling", title),
            Paragraph("Implementation rationale, numerical behavior, and questions about mathematical and discrete consistency", subtitle),
            Table([[Paragraph("CURRENT IMPLEMENTATION", small),
                    Paragraph("code/flow/working_flow_coupling_simple_long_tube_pressure_trace_darcy_fe_current_mesh.lua", small)],
                   [Paragraph("RESULTS DISCUSSED", small),
                    Paragraph("pressure_trace_darcy_fe_ref2_{lumen,membrane,inter}_t0000.vtu", small)]],
                  colWidths=[37*mm, 136*mm], style=[("BACKGROUND",(0,0),(0,-1),PALE_BLUE),
                  ("BOX",(0,0),(-1,-1),0.5,colors.HexColor("#B7C8D8")),
                  ("INNERGRID",(0,0),(-1,-1),0.3,colors.HexColor("#D6E0E8")),
                  ("VALIGN",(0,0),(-1,-1),"TOP"),("LEFTPADDING",(0,0),(-1,-1),5),
                  ("RIGHTPADDING",(0,0),(-1,-1),5),("TOPPADDING",(0,0),(-1,-1),4),
                  ("BOTTOMPADDING",(0,0),(-1,-1),4)]),
            Spacer(1, 5),
            RLImage(SCHEMATIC_FIG, width=173*mm, height=51*mm),
            Paragraph("Figure 1. Domains and coupling conditions represented in the current Lua implementation.", caption),
            Paragraph("1. Purpose and intended mathematical model", h1),
            Paragraph("We are developing a scalable monolithic model for flow from a tubular lumen into a porous membrane and surrounding interstitium. We are sharing this note to explain the choices made so far and to ask for guidance on whether the continuum interface model and its UG4 realization are mathematically consistent. We do <b>not</b> regard the present formulation as fully validated.", callout),
            Paragraph("The lumen is modeled by steady incompressible Stokes flow. Membrane and Inter are modeled by Darcy flow. Five P1 unknowns share one UG4 ApproximationSpace and one monolithic matrix: <i>(u,v,w,p<sub>S</sub>,p<sub>D,s</sub>)</i>. The porous pressure is stored as <i>p<sub>D,s</sub>=p<sub>D</sub>/S</i>, with <i>S=10<sup>4</sup></i>, and the Darcy diffusion coefficient is multiplied by the same factor. This was introduced to balance algebraic block magnitudes and is intended not to alter the physical solution.", body),
            Paragraph("q<sub>D</sub>=-(K/nu) grad(p<sub>D</sub>), &nbsp; div(q<sub>D</sub>)=0, &nbsp; div(u<sub>S</sub>)=0", eq),
            Paragraph("The intended Apical laws are normal mass conservation, Darcy pressure acting through the normal Stokes traction, and Beavers-Joseph-Saffman tangential slip. The code also imposes <i>p<sub>S</sub>=p<sub>D</sub></i>. Our reasoning was that, after neglecting the viscous normal-stress correction, <i>-n.sigma<sub>S</sub>n=p<sub>D</sub></i> approximately reduces to pressure equality. We need advice on whether using both an explicit Darcy-pressure traction and pressure equality duplicates or overconstrains the normal condition.", body),
            Paragraph("Boundary conditions used in the reported run", h2),
            table([
                ["Subset", "Current choice and reason"],
                ["Inlet", "Prescribed parabolic velocity; vmax = 3 mm/s and radius = 0.01 mm"],
                ["Outlet", "Stokes zero normal-stress outflow"],
                ["Apical", "Stokes flux closure, normal traction, pressure trace constraint, and BJS slip"],
                ["TripleJunction", "No-slip wall; included in the trace constraint to constrain endpoint values"],
                ["OuterWall", "Darcy pressure fixed at 28.7 Pa as an external pressure reference"],
                ["Basolateral", "Internal conforming Membrane/Inter interface; no external boundary condition"],
            ], [38*mm, 135*mm], font=7.05),
            PageBreak(),
            Paragraph("2. Why the present discretization and solver were chosen", h1),
            table([
                ["Choice", "Reason for choosing it", "Current uncertainty"],
                ["Two pressures", "Stokes pressure is an incompressibility multiplier; Darcy pressure solves an elliptic equation. A single field gave poor scaling/behavior.", "Whether a properly scaled shared trace is preferable."],
                ["Scaled pD", "Stokes and porous blocks differed strongly in magnitude; scaling improved algebraic balance without intentionally changing physics.", "Correct scaling for mesh-independent GMG."],
                ["FV1 Stokes", "UG4 Navier-Stokes, BJS wall operator, and an established lumen formulation.", "Interaction of equal-order stabilization with coupling."],
                ["P1 FE Darcy", "Simple elliptic operator and nodal trace compatible with pressure elimination; this version reached numRefs=2.", "Gradient velocity is nonconservative; FV is preferred long term."],
                ["Trace elimination", "Avoids a zero-diagonal multiplier and equates pressure values at shared interface vertices.", "Whether merged FE and FV rows form a consistent flux law."],
                ["Normal factor 1e3", "Added to damp an interface-normal mode.", "Numerical stabilization, not yet a justified membrane law."],
            ], [36*mm, 84*mm, 53*mm], font=6.35),
            Spacer(1, 5),
            Paragraph("Solver used in this file", h2),
            Paragraph("The saddle-point matrix is solved by <b>BiCGStab + geometric multigrid</b>. The GMG smoother is Uzawa: ILU approximates the velocity/Darcy main block, and an ILU-applied pressure-mass operator approximates the Stokes Schur complement. V-cycles use two pre- and post-smoothing steps. RAP/Galerkin coarse operators are enabled to retain algebraically modified interface rows. SuperLU is restricted to the level-0 coarse problem. We recompute <i>b-Ax</i> for every component and require absolute defect <= 1e-8.", body),
            Paragraph("Configurations tried and what motivated each change", h2),
            table([
                ["Configuration", "Observed behavior", "Reason for next change"],
                ["Full LU/SuperLU", "Reliable on small systems.", "Cannot support intended mesh sizes or MPI scaling."],
                ["BiCGStab + ILU", "Solved some coarse cases but deteriorated with refinement.", "A multilevel preconditioner was required."],
                ["GMG + Uzawa/ILU", "This FE version solved numRefs=2; pressure patches followed refinement generations and higher levels could stagnate.", "Investigate constraint transfer and smoother robustness."],
                ["GMG without RAP", "Constrained interface rows were not adequately preserved; convergence worsened.", "Use RAP/Galerkin coarse operators."],
                ["GMG + Element Vanka", "A later FV diagnostic on Long_tube_V3 converged through ref2 more robustly.", "Promising, but not proof that coupling is correct."],
            ], [42*mm, 79*mm, 52*mm], font=6.35),
            Spacer(1, 5),
            Spacer(1, 4),
            Paragraph("<b>Why GMG remains necessary.</b> The final nephron meshes will be too large for a global direct solve. We therefore need a formulation whose iteration count, pressure field, and interface flux remain stable under refinement and MPI partitioning. Convergence of the algebraic residual alone is not sufficient if the coupled physical flux is incorrect.", callout),
            PageBreak(),
            Paragraph("3. Observed problem, tests performed, and guidance requested", h1),
            RLImage(PRESSURE_FIG, width=173*mm, height=63*mm),
            Paragraph(f"Figure 2. Darcy pressure rendered directly from the membrane VTU at numRefs=2 (nodal range {pmin:.3e} to {pmax:.3e} mm2/s2). The concern is mesh-aligned banding: level 0 appeared smooth, while vertices created at levels 1 and 2 developed systematic differences.", caption),
            Paragraph("Tests already performed", h2),
            table([
                ["Test/change", "Result and present interpretation"],
                ["Constrained TripleJunction", "Removed endpoint trace mismatch but not all pressure banding."],
                ["Simultaneous slave elimination and full column scan", "Removed vertex-order dependence that could imprint coarse-grid patterns; artifacts remain."],
                ["Equal Membrane/Inter permeability", "Pattern persisted; a Basolateral permeability jump is not the sole cause."],
                ["Projector vs ProMesh control meshes", "Projector was suspected, but recent Long_tube_V3 statistics localized the largest deviations near end rings."],
                ["OuterWall pressure and traction-sign tests", "Changed mean pressure and leakage direction but did not cure refinement behavior; the exact sign convention still needs confirmation for this file."],
                ["Stokes vs reconstructed Darcy flux", "DarcyVelocityLinker did not close the balance; it may not represent the assembled conservative flux."],
            ], [63*mm, 110*mm], font=6.35),
            Spacer(1, 4),
            Paragraph("Questions on which we need UG4 and modeling guidance", h2),
            Paragraph("1. <b>Continuum interface law:</b> after neglecting viscous normal stress, should we impose both <i>p<sub>S</sub>=p<sub>D</sub></i> and Darcy-pressure traction, or only the normal-stress law together with normal velocity continuity and BJS slip?", bullet),
            Paragraph("2. <b>Discrete conservation:</b> is adding an FE Darcy residual row directly to an FV1 Stokes continuity row at a shared vertex mathematically valid? Are control-volume weights, normals, signs, or row scaling missing?", bullet),
            Paragraph("3. <b>Function spaces:</b> should Darcy use FV for local conservation, or should FE Darcy be coupled by a mortar, Lagrange multiplier, or Nitsche operator?", bullet),
            Paragraph("4. <b>Pressure variables:</b> should Stokes and Darcy retain separate scaled pressures, or is there a stable UG4 formulation with a shared interface pressure?", bullet),
            Paragraph("5. <b>Multigrid:</b> which smoother and transfer strategy preserves the mixed interface law under refinement? Is RAP sufficient, or must coupling be assembled explicitly at every level?", bullet),
            Paragraph("6. <b>Verification:</b> how can we obtain the actual assembled Darcy boundary flux in UG4 rather than an elementwise gradient reconstruction?", bullet),
            Paragraph("Our proposed next step is a minimal manufactured-solution or straight-tube benchmark with known interface flux. We would compare FE and FV Darcy formulations, verify flux balance from the assembled operators, and only then calibrate OuterWall pressure or a Robin tissue-resistance law. Guidance on the mathematically appropriate UG4 formulation would be extremely helpful.", callout),
        ]
        return content

    # PAGE 1
    story += [
        Paragraph("Monolithic Stokes-Darcy Coupling in a Tubular Geometry", title),
        Paragraph("Assembly, solver configuration, observed pressure artifact, and open questions for independent review", subtitle),
        Table([[Paragraph("FILE UNDER REVIEW", small),
                Paragraph("code/flow/working_flow_coupling_simple_long_tube_pressure_trace_darcy_fe_current_mesh.lua", small)],
               [Paragraph("REFERENCE OUTPUT", small),
                Paragraph("pressure_trace_darcy_fe_ref2_{lumen,membrane,inter}_t0000.vtu", small)]],
              colWidths=[34*mm, 139*mm], style=[("BACKGROUND",(0,0),(0,-1),PALE_BLUE),
              ("BOX",(0,0),(-1,-1),0.5,colors.HexColor("#B7C8D8")),
              ("INNERGRID",(0,0),(-1,-1),0.3,colors.HexColor("#D6E0E8")),
              ("VALIGN",(0,0),(-1,-1),"TOP"),("LEFTPADDING",(0,0),(-1,-1),5),
              ("RIGHTPADDING",(0,0),(-1,-1),5),("TOPPADDING",(0,0),(-1,-1),4),
              ("BOTTOMPADDING",(0,0),(-1,-1),4)]),
        Spacer(1, 5),
        RLImage(SCHEMATIC_FIG, width=173*mm, height=51*mm),
        Paragraph("Figure 1. Domain decomposition and interface conditions used by the Lua model.", caption),
        Paragraph("1. What the code solves", h1),
        Paragraph("The script assembles one steady linear system for five P1 unknowns: lumen velocity <i>(u,v,w)</i>, lumen kinematic pressure <i>p<sub>S</sub></i>, and a scaled porous pressure <i>p<sub>D,s</sub></i>. Stokes equations are assembled only in <b>Lumen</b>. Darcy equations are assembled in <b>Membrane</b> and <b>Inter</b>. All unknowns share one UG4 ApproximationSpace and one monolithic matrix.", body),
        Paragraph("Darcy scaling is algebraic: <i>p<sub>D</sub> = S p<sub>D,s</sub></i>, with <i>S = 10<sup>4</sup></i>. The Darcy diffusion coefficient is multiplied by the same factor, so scaling should change conditioning but not the physical pressure or Darcy velocity.", body),
        Paragraph("q<sub>D</sub> = -(K/nu) grad(p<sub>D</sub>), &nbsp;&nbsp; p<sub>S</sub> = p<sub>D</sub> on Apical, &nbsp;&nbsp; u<sub>S</sub>.n = q<sub>D</sub>.n", eq),
        Paragraph("At <b>Apical</b>, the implementation combines: (i) equality of the physical pressure traces using PressureTraceConstraint, (ii) normal conservation by merging the eliminated Darcy conservation row into the Stokes continuity row, (iii) Darcy-pressure traction in the three Stokes momentum rows, and (iv) Beavers-Joseph-Saffman tangential slip. The viscous normal-stress correction is neglected.", body),
        Paragraph("Boundary conditions", h2),
        table([
            ["Boundary/subset", "Condition in this file"],
            ["Inlet", "Prescribed parabolic velocity, vmax = 3 mm/s; radius = 0.01 mm"],
            ["Outlet", "Stokes zero normal-stress outflow"],
            ["Apical", "Open Stokes flux closure + BJS + pressure/normal coupling"],
            ["TripleJunction", "No-slip wall and included in pressure trace constraint"],
            ["OuterWall", "Darcy Dirichlet pressure = 28.7 Pa"],
            ["Basolateral", "Internal Membrane/Inter interface; no external pressure BC"],
        ], [38*mm, 135*mm], font=7.25),
    ]

    story.append(PageBreak())

    # PAGE 2
    story += [
        Paragraph("2. Discretization and linear solver", h1),
        table([
            ["Block", "Discretization / role"],
            ["Stokes momentum + continuity", "NavierStokesFV1, equal-order P1/P1, fields/cor stabilization, symmetric-gradient stress"],
            ["Darcy pressure", "ConvectionDiffusionFE, P1, separately assembled in Membrane and Inter"],
            ["Pressure trace", "Custom algebraic elimination: pD,s = pS/S on Apical and TripleJunction"],
            ["BJS", "NavierStokesWSBCFV1; tangential friction alpha*nu/sqrt(Kinterface); normal factor = 1e3"],
            ["Pressure reference", "OuterWall Dirichlet value; no Darcy condition on Basolateral"],
        ], [48*mm, 125*mm], font=7.25),
        Spacer(1, 5),
        Paragraph("Current solver path", h2),
        Paragraph("The assembled saddle-point system is solved by <b>BiCGStab</b> preconditioned with <b>geometric multigrid (GMG)</b>. Each V-cycle uses two pre- and two post-smoothing passes. The smoother is an <b>Uzawa</b> iteration: ILU approximates the velocity/Darcy main block and an ILU-applied pressure-mass operator approximates the Stokes Schur complement. Galerkin/RAP coarse operators are enabled so that algebraically modified interface rows are inherited. <b>SuperLU</b> is used only for the small level-0 coarse problem, not for the full refined system.", body),
        Paragraph("Stopping rule: maximum 400 BiCGStab iterations per cycle; global absolute defect <= 1e-8; global reduction <= 1e-8; up to eight restarted cycles; and an independently recomputed <i>b-Ax</i> defect <= 1e-8 for every component.", callout),
        Paragraph("Solver/configuration history", h2),
        table([
            ["Configuration tested", "Observed behavior", "Assessment"],
            ["Full-system LU / SuperLU", "Reliable for small meshes; memory and factorization cost grow rapidly; unsuitable for intended MPI-scale problem.", "Debug reference only"],
            ["BiCGStab + single-level ILU", "Worked on smaller cases but deteriorated with refinement and did not provide scalable coarse correction.", "Not mesh independent"],
            ["BiCGStab + GMG + Uzawa/ILU (this file)", "Historical configuration solved numRefs=2 and produced the supplied VTUs, but pressure patches track refinement generations; higher levels can stagnate.", "Current FE baseline"],
            ["GMG without RAP", "Coarse operators did not faithfully retain algebraically merged pressure-trace rows; convergence was worse.", "Rejected"],
            ["GMG with RAP", "Improved preservation of interface constraint and reduced defects, but did not by itself remove pressure patterning.", "Retained"],
            ["Element Vanka GMG (later FV diagnostic)", "More robust for thin anisotropic layers; recent Long_tube_V3 tests converged through ref2, but this is not yet the formulation in the reviewed file.", "Promising comparison"],
        ], [43*mm, 91*mm, 39*mm], font=6.65),
        Spacer(1, 5),
        Paragraph("Important implementation changes already tried", h2),
        Paragraph("- Included TripleJunction in PressureTraceConstraint. This removed the previously measurable pressure-trace mismatch at interface endpoint vertices but did not remove all refined-grid patterning.", bullet),
        Paragraph("- Reworked PressureTraceConstraint to eliminate all Darcy trace slaves simultaneously and scan complete matrix columns. The former sequential elimination could make the matrix depend on interface vertex ordering and imprint coarse-grid patches.", bullet),
        Paragraph("- Tested equal Membrane/Inter permeability to isolate the Basolateral interface. The visual pattern persisted, so a permeability jump alone is not the explanation.", bullet),
        Paragraph("- Compared projector-generated and ProMesh control meshes. Recent quantitative tests on Long_tube_V3 found newly refined interior vertices no rougher than inherited vertices; the largest deviations were concentrated near inlet/outlet rings.", bullet),
        Paragraph("- Tested OuterWall pressure values and traction signs. Outer pressure changes the mean pressure and leakage direction, but changing this right-hand-side value did not cure refinement-dependent solver behavior. The UG4 Neumann sign convention was also identified as a separate concern in the later FV version.", bullet),
    ]

    story.append(PageBreak())

    # PAGE 3
    story += [
        Paragraph("3. Observed pressure artifact and present diagnosis", h1),
        RLImage(PRESSURE_FIG, width=173*mm, height=68*mm),
        Paragraph(f"Figure 2. Re-rendered directly from the supplied membrane VTU at numRefs=2 (full nodal range {pmin:.3e} to {pmax:.3e} mm2/s2). The user-observed symptom is a banded or patch-like Darcy pressure aligned with coarse/refined mesh generations, especially near the ends and on the Basolateral side.", caption),
        Paragraph("What is established", h2),
        table([
            ["Evidence", "Interpretation"],
            ["Level 0 pressure appeared smooth; level 1 showed differences at new vertices; level 2 showed the pattern on both refinement generations.", "The artifact is correlated with hierarchy/refinement structure, although projector geometry is not proven to be the sole cause."],
            ["PressureTraceConstraint can reduce the Apical L2 pressure jump to near machine precision.", "Literal trace equality can be satisfied while the volume pressure field or flux reconstruction remains questionable."],
            ["Changing OuterWall pressure altered direction/magnitude but did not remove refined-level convergence problems.", "OuterWall pressure is physically important, but it is not sufficient as an explanation for the mesh-aligned pressure patches."],
            ["Raw DarcyVelocityLinker surface flux does not match the Stokes Apical integral or close global Darcy balance in later FV tests.", "The reconstructed gradient velocity is not necessarily the conservative discrete interface flux; it should not be used alone to validate coupling."],
            ["Recent ref0-ref2 tests showed lumen surface-integral imbalance decreasing with refinement.", "Some apparent end behavior is a boundary quadrature/discrete-flux diagnostic issue rather than a converged physical flux split."],
        ], [64*mm, 109*mm], font=6.75),
        Spacer(1, 4),
        Paragraph("Most likely unresolved causes (ranked)", h2),
        Paragraph("1. <b>Interface conservation implementation:</b> verify that the merged FE Darcy row and FV1 Stokes continuity row have compatible signs, units, control-volume weights, and transfer behavior on every GMG level.", bullet),
        Paragraph("2. <b>FE/FV mismatch and postprocessing:</b> this baseline couples FV1 Stokes to FE Darcy. The FE gradient-derived velocity is elementwise and nonconservative. A consistent FV Darcy formulation or a true assembled boundary-flux diagnostic is preferred for the long-term model.", bullet),
        Paragraph("3. <b>End-ring subset treatment:</b> the first ring after Inlet and last ring before Outlet showed velocity fans and the largest pressure deviations. TripleJunction, Inlet/Outlet, and Apical intersection vertices require a single unambiguous boundary classification.", bullet),
        Paragraph("4. <b>Normal interface impedance:</b> the current value 1e3 is a numerical stabilization, far below a simple physical membrane resistance estimate delta*nu/K (previously of order 5e6). A too-small value can permit excessive Apical flow; a large value can worsen conditioning unless represented consistently in GMG.", bullet),
        Paragraph("5. <b>OuterWall model:</b> a fixed pressure is mathematically refinement independent, but the chosen value controls leakage. After discrete normal conservation is verified, a physical Robin/resistance law q.n=(pD-pext)/Rext may be preferable to calibrating a Dirichlet value to a target split.", bullet),
        Paragraph("Questions for the reviewer", h2),
        Paragraph("- Is adding the unscaled FE Darcy residual row directly to the FV1 Stokes continuity row a consistent discrete enforcement of normal-flux continuity?", bullet),
        Paragraph("- Should the Apical normal condition be imposed with a mortar/Lagrange-multiplier or Nitsche interface formulation instead of algebraic trace elimination?", bullet),
        Paragraph("- Which block smoother and transfer operator best preserve this mixed Stokes-Darcy interface under geometric refinement?", bullet),
        Paragraph("- Can UG4 expose the assembled FV/FE boundary flux directly so conservation is tested without DarcyVelocityLinker reconstruction?", bullet),
    ]

    # The consultation version supersedes the earlier code-review-oriented
    # draft retained above while the document was being developed.
    story = consultation_story()
    doc.build(story)


if __name__ == "__main__":
    os.makedirs(TMP, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    pmin, pmax = create_pressure_figure()
    create_schematic()
    make_pdf(pmin, pmax)
    print(PDF_PATH)
