from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.platypus import (
    BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer, PageBreak,
    Table, TableStyle, KeepTogether, Preformatted, NextPageTemplate
)
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.graphics.shapes import Drawing, Rect, String, Line, Polygon, Circle


ROOT = Path("/Users/alhilaba/Documents/Promesh_swc/Model")
OUT = ROOT / "output/pdf/Urea_flow_projector_transport_technical_summary.pdf"

NAVY = colors.HexColor("#17324D")
BLUE = colors.HexColor("#2878B5")
CYAN = colors.HexColor("#55A8C7")
PALE = colors.HexColor("#EAF3F8")
INK = colors.HexColor("#25313C")
MUTED = colors.HexColor("#5E6B75")
GREEN = colors.HexColor("#2D7A56")
GREEN_BG = colors.HexColor("#E7F4ED")
AMBER = colors.HexColor("#A86708")
AMBER_BG = colors.HexColor("#FFF3D9")
RED = colors.HexColor("#A63A3A")
RED_BG = colors.HexColor("#FBE9E8")
LIGHT = colors.HexColor("#F5F7F9")
GRID = colors.HexColor("#C9D3DB")


def register_fonts():
    candidates = [
        ("/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Supplemental/Arial Bold.ttf"),
        ("/System/Library/Fonts/Supplemental/Helvetica.ttf", "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf"),
    ]
    for regular, bold in candidates:
        if Path(regular).exists() and Path(bold).exists():
            pdfmetrics.registerFont(TTFont("BodyFont", regular))
            pdfmetrics.registerFont(TTFont("BodyFont-Bold", bold))
            return "BodyFont", "BodyFont-Bold"
    return "Helvetica", "Helvetica-Bold"


BODY_FONT, BOLD_FONT = register_fonts()


class ReportDoc(BaseDocTemplate):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        body = Frame(18 * mm, 18 * mm, A4[0] - 36 * mm, A4[1] - 35 * mm, id="body")
        cover = Frame(0, 0, A4[0], A4[1], leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0, id="cover")
        self.addPageTemplates([
            PageTemplate(id="Cover", frames=[cover], onPage=self.cover_page),
            PageTemplate(id="Body", frames=[body], onPage=self.body_page),
        ])

    def cover_page(self, canvas, doc):
        canvas.saveState()
        canvas.setFillColor(NAVY)
        canvas.rect(0, 0, A4[0], A4[1], fill=1, stroke=0)
        canvas.setFillColor(CYAN)
        canvas.rect(0, 0, 16 * mm, A4[1], fill=1, stroke=0)
        canvas.setFillColor(colors.HexColor("#224663"))
        canvas.circle(A4[0] - 20 * mm, A4[1] - 22 * mm, 38 * mm, fill=1, stroke=0)
        canvas.circle(A4[0] - 8 * mm, 18 * mm, 28 * mm, fill=1, stroke=0)
        canvas.restoreState()

    def body_page(self, canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(GRID)
        canvas.line(18 * mm, A4[1] - 13 * mm, A4[0] - 18 * mm, A4[1] - 13 * mm)
        canvas.setFont(BODY_FONT, 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(18 * mm, A4[1] - 10 * mm, "UG4 nephron flow and urea transport - technical implementation record")
        canvas.drawRightString(A4[0] - 18 * mm, 10 * mm, f"Page {doc.page}")
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph):
            name = flowable.style.name
            if name in ("H1", "H2"):
                level = 0 if name == "H1" else 1
                text = flowable.getPlainText()
                key = f"section-{self.seq.nextf('section')}"
                self.canv.bookmarkPage(key)
                self.canv.addOutlineEntry(text, key, level=level, closed=False)
                self.notify("TOCEntry", (level, text, self.page, key))


styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="CoverTitle", fontName=BOLD_FONT, fontSize=25, leading=30, textColor=colors.white, alignment=TA_LEFT, leftIndent=28 * mm, rightIndent=24 * mm, spaceAfter=12))
styles.add(ParagraphStyle(name="CoverSub", fontName=BODY_FONT, fontSize=12, leading=18, textColor=colors.HexColor("#D8EAF2"), leftIndent=28 * mm, rightIndent=24 * mm, spaceAfter=10))
styles.add(ParagraphStyle(name="CoverMeta", fontName=BODY_FONT, fontSize=9, leading=14, textColor=colors.HexColor("#B9D4E2"), leftIndent=28 * mm, rightIndent=24 * mm))
styles.add(ParagraphStyle(name="H1", fontName=BOLD_FONT, fontSize=17, leading=21, textColor=NAVY, spaceBefore=7, spaceAfter=9, keepWithNext=True))
styles.add(ParagraphStyle(name="H2", fontName=BOLD_FONT, fontSize=12.5, leading=16, textColor=BLUE, spaceBefore=10, spaceAfter=5, keepWithNext=True))
styles.add(ParagraphStyle(name="H3", fontName=BOLD_FONT, fontSize=10.2, leading=13, textColor=INK, spaceBefore=7, spaceAfter=4, keepWithNext=True))
styles.add(ParagraphStyle(name="BodyTextX", fontName=BODY_FONT, fontSize=9.15, leading=13.2, textColor=INK, spaceAfter=6))
styles.add(ParagraphStyle(name="Small", fontName=BODY_FONT, fontSize=7.6, leading=10.5, textColor=MUTED, spaceAfter=4))
styles.add(ParagraphStyle(name="BulletX", fontName=BODY_FONT, fontSize=8.9, leading=12.7, textColor=INK, leftIndent=13, firstLineIndent=-7, bulletIndent=4, spaceAfter=4))
styles.add(ParagraphStyle(name="Callout", fontName=BODY_FONT, fontSize=9, leading=13, textColor=INK, leftIndent=7, rightIndent=7, borderColor=CYAN, borderWidth=0.8, borderPadding=7, backColor=PALE, spaceBefore=5, spaceAfter=8))
styles.add(ParagraphStyle(name="CodeX", fontName="Courier", fontSize=7.2, leading=9.4, textColor=INK, backColor=LIGHT, borderColor=GRID, borderWidth=0.5, borderPadding=6, spaceBefore=4, spaceAfter=8))
styles.add(ParagraphStyle(name="TableHead", fontName=BOLD_FONT, fontSize=7.4, leading=9, textColor=colors.white))
styles.add(ParagraphStyle(name="TableCell", fontName=BODY_FONT, fontSize=7.15, leading=9.4, textColor=INK))
styles.add(ParagraphStyle(name="TableCellBold", fontName=BOLD_FONT, fontSize=7.15, leading=9.4, textColor=INK))
styles.add(ParagraphStyle(name="Caption", fontName=BODY_FONT, fontSize=7.2, leading=9.4, textColor=MUTED, alignment=TA_CENTER, spaceBefore=3, spaceAfter=8))


def P(text, style="BodyTextX"):
    return Paragraph(text, styles[style])


def bullets(items):
    return [Paragraph(item, styles["BulletX"], bulletText="-") for item in items]


def section(title):
    return Paragraph(title, styles["H1"])


def subsection(title):
    return Paragraph(title, styles["H2"])


def table(data, widths, header=True, font_size=None):
    cooked = []
    for r, row in enumerate(data):
        cooked.append([P(str(cell), "TableHead" if header and r == 0 else "TableCell") for cell in row])
    t = Table(cooked, colWidths=widths, repeatRows=1 if header else 0, hAlign="LEFT")
    commands = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, GRID),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]
    if header:
        commands += [("BACKGROUND", (0, 0), (-1, 0), NAVY)]
        if len(data) > 1:
            commands += [("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT])]
    t.setStyle(TableStyle(commands))
    return t


def status_table():
    data = [
        ["Subsystem", "Current state", "Evidence"],
        ["Projector hierarchy", "WORKING", "Refinement is completed before distribution; level-1 MPI hierarchy is built."],
        ["Stationary lumen flow", "WORKING", "Serial ref0 and MPI-4 ref1 both converged."],
        ["Flow-to-transport velocity", "WORKING", "ExplicitGridFunctionVector supplies u,v,w as valid UG4 user data."],
        ["Short adaptive transport", "WORKING", "Serial and MPI short-horizon LIMEX checks accepted steps."],
        ["Full 0-60 s, refined MPI transport", "NOT YET PROVEN", "BiCGStab failed; GMRES is now configured but has not been validated for the full run."],
    ]
    t = table(data, [37 * mm, 31 * mm, 104 * mm])
    t.setStyle(TableStyle([
        ("TEXTCOLOR", (1, 1), (1, 4), GREEN),
        ("TEXTCOLOR", (1, 5), (1, 5), AMBER),
        ("FONTNAME", (1, 1), (1, -1), BOLD_FONT),
    ]))
    return t


def architecture_diagram():
    d = Drawing(485, 220)
    boxes = [
        (12, 163, 92, 38, "Original UGX\n+ projector data", PALE),
        (126, 163, 92, 38, "Projector refine\non one rank", PALE),
        (240, 163, 92, 38, "Distribute full\nGMG hierarchy", PALE),
        (375, 163, 98, 38, "Stokes FV1\nin lumen", GREEN_BG),
        (375, 92, 98, 38, "Fixed velocity\nuser data", GREEN_BG),
        (240, 21, 110, 43, "Adaptive LIMEX\ntransport stages", AMBER_BG),
        (92, 21, 112, 43, "Lumen / membrane /\ninterstitium fields", PALE),
        (12, 92, 112, 38, "Apical + basolateral\nmembrane fluxes", PALE),
    ]
    for x, y, w, h, label, fill in boxes:
        d.add(Rect(x, y, w, h, rx=6, ry=6, fillColor=fill, strokeColor=BLUE, strokeWidth=1))
        lines = label.split("\n")
        for i, line in enumerate(lines):
            d.add(String(x + w / 2, y + h / 2 + 5 - i * 11, line, textAnchor="middle", fontName=BOLD_FONT if i == 0 else BODY_FONT, fontSize=8, fillColor=INK))
    arrows = [(104, 182, 126, 182), (218, 182, 240, 182), (332, 182, 375, 182), (424, 163, 424, 130), (375, 111, 335, 64), (240, 42, 204, 42), (92, 42, 66, 92)]
    for x1, y1, x2, y2 in arrows:
        d.add(Line(x1, y1, x2, y2, strokeColor=NAVY, strokeWidth=1.4))
        angle = 5
        if abs(x2 - x1) >= abs(y2 - y1):
            s = 1 if x2 > x1 else -1
            d.add(Polygon([x2, y2, x2 - s * 6, y2 + angle, x2 - s * 6, y2 - angle], fillColor=NAVY, strokeColor=NAVY))
        else:
            s = 1 if y2 > y1 else -1
            d.add(Polygon([x2, y2, x2 + angle, y2 - s * 6, x2 - angle, y2 - s * 6], fillColor=NAVY, strokeColor=NAVY))
    return d


def subset_diagram():
    d = Drawing(485, 190)
    cx, cy = 235, 95
    d.add(Circle(cx, cy, 86, fillColor=colors.HexColor("#DDE9D9"), strokeColor=GREEN, strokeWidth=1.5))
    d.add(Circle(cx, cy, 61, fillColor=colors.white, strokeColor=CYAN, strokeWidth=8))
    d.add(Circle(cx, cy, 38, fillColor=colors.HexColor("#DDEFFC"), strokeColor=BLUE, strokeWidth=2))
    d.add(String(cx, cy + 2, "Lumen", textAnchor="middle", fontName=BOLD_FONT, fontSize=10, fillColor=NAVY))
    d.add(String(cx, cy + 52, "Membrane", textAnchor="middle", fontName=BOLD_FONT, fontSize=8, fillColor=NAVY))
    d.add(String(cx, cy + 74, "Interstitium", textAnchor="middle", fontName=BOLD_FONT, fontSize=8, fillColor=NAVY))
    d.add(Line(274, 62, 370, 31, strokeColor=BLUE, strokeWidth=1))
    d.add(String(374, 28, "Apical: lumen <-> membrane", fontName=BODY_FONT, fontSize=8, fillColor=INK))
    d.add(Line(294, 125, 370, 151, strokeColor=GREEN, strokeWidth=1))
    d.add(String(374, 148, "Basolateral: membrane <-> inter", fontName=BODY_FONT, fontSize=8, fillColor=INK))
    d.add(Line(190, 65, 87, 31, strokeColor=RED, strokeWidth=1))
    d.add(String(3, 18, "No interstitium/Apical coupling", fontName=BOLD_FONT, fontSize=8, fillColor=RED))
    return d


story = []

# Cover
story += [
    Spacer(1, 38 * mm),
    Paragraph("Flow, Projector Refinement,<br/>and Urea Transport in UG4", styles["CoverTitle"]),
    Paragraph("Implementation history, numerical design, failures, convergence evidence, and present validation status", styles["CoverSub"]),
    Spacer(1, 10 * mm),
    Paragraph("Primary implementation:<br/><b>Urea_trans_with_flow_coupled_gmg.lua</b>", styles["CoverMeta"]),
    Spacer(1, 90 * mm),
    Paragraph("Prepared from the current Lua source and recorded serial/MPI solver traces<br/>30 September 2026", styles["CoverMeta"]),
    NextPageTemplate("Body"), PageBreak(),
]

story += [section("Executive summary"),
    P("The final architecture is a <b>one-way sequential flow/transport model</b>. A stationary Stokes field is solved once in the lumen. That solved velocity is then held fixed while urea is advanced adaptively through lumen convection-diffusion, membrane diffusion, interstitial diffusion, and two membrane-permeability interfaces. This is not a monolithic flow-transport solve and it no longer recomputes an unchanged flow field at every transport step."),
    P("The stationary flow portion is working on the tested meshes. The projector-aware refinement and MPI distribution sequence is also working at one refinement level. Short adaptive transport tests have passed. The complete 0-60 s refined MPI transport run is <b>not yet demonstrated</b>: BiCGStab failed inside a LIMEX stage, and the replacement GMRES configuration has only passed a very short structural test."),
    status_table(), Spacer(4 * mm, 1),
    P("The current code was also reorganized to follow the visual structure of <b>2D_dif_mig_no_elec_pot.lua</b>: direct top-level inputs, literal solver tolerances, and one-line commands where the line remains no longer than the reference file's 181-character maximum.", "Callout"),
    PageBreak(),
]

toc = TableOfContents()
toc.levelStyles = [
    ParagraphStyle(name="TOC0", fontName=BOLD_FONT, fontSize=10, leading=15, leftIndent=0, textColor=NAVY, spaceBefore=3),
    ParagraphStyle(name="TOC1", fontName=BODY_FONT, fontSize=8.5, leading=12, leftIndent=13, textColor=INK),
]
story += [section("Contents"), toc, PageBreak()]

story += [section("1. Model scope and governing equations"),
    subsection("1.1 Stationary lumen flow"),
    P("The code uses UG4's FV1 Stokes formulation in the three-dimensional lumen. With inertia disabled, the continuous model is the incompressible stationary system"),
    Preformatted("-nu Laplacian(u) + grad(p) = 0\ndiv(u) = 0", styles["CodeX"]),
    P("where u = (u,v,w), p is pressure, and nu = 1.0 mm^2/s is the configured kinematic viscosity. The call <b>set_stokes(true)</b> removes the nonlinear advective momentum term. Consequently, with fixed geometry and boundary data, the velocity does not evolve in time and only needs to be solved once."),
    subsection("1.2 Urea transport"),
    P("The normalized concentration C* obeys convection-diffusion in the lumen and diffusion in the membrane and interstitium. In compact form:"),
    Preformatted("Lumen:       dC_L*/dt + div(u C_L*) = div(D grad(C_L*))\nMembrane:    dC_M*/dt              = div(D grad(C_M*)) + interface fluxes\nInterstitium:dC_I*/dt              = div(D grad(C_I*)) + interface fluxes", styles["CodeX"]),
    P("The diffusion coefficient is D = 1.67e-3 mm^2/s. The apical flux exchanges lumen and membrane concentration with permeability 3.78e-4 mm/s. The basolateral flux exchanges membrane and interstitium with permeability 3.70e-4 mm/s."),
    subsection("1.3 Concentration scaling"),
    P("The physical inlet concentration is 5.0e-6 mol/mm^3. The numerical field is scaled as C* = C / (5.0e-6 mol/mm^3), so the inlet value is exactly 1. This keeps the solution near order one and avoids asking iterative solvers and temporal error estimators to distinguish values clustered near 1e-6."),
    P("To recover physical units in ParaView, use <b>C = Urea_normalized x 5.0e-6 mol/mm^3</b>. Scaling improves numerical magnitude, but it does not by itself fix a wrong velocity direction, an invalid subset coupling, or a deficient preconditioner."),
    PageBreak(),
]

story += [section("2. End-to-end implementation"),
    architecture_diagram(),
    P("Figure 1. Data flow from the serialized projector grid to fixed-flow adaptive transport.", "Caption"),
    subsection("2.1 Direct run inputs"),
    P("The main inputs are now direct assignments at the beginning of the Lua file, matching the style of the 2D reference. The current defaults are grid runs/third_trace_mm.ugx, one refinement, end time 60 s, output every 1 s, initial LIMEX step 1e-3 s, minimum step 1e-12 s, maximum step 1 s, and temporal tolerance 1e-3."),
    P("This has an operational consequence: former command-line flags such as -numRefs or -endTime no longer change these values. To change a run, edit the corresponding top-of-file assignment. Solver absolute and relative tolerances are likewise literal values in the solver descriptions."),
    subsection("2.2 Two approximation spaces"),
    P("Flow and transport use separate approximation spaces. The flow space contains u, v, w, and p only on the lumen and its boundary subsets. The transport space contains three scalar fields with different supports. This separation prevents pressure or concentration components from being evaluated where they have no degrees of freedom."),
    table([
        ["Field", "Supported subsets", "Purpose"],
        ["u,v,w,p", "Lumen, Inlet, Outlet, Apical, TripleJunction", "Stationary Stokes field"],
        ["urea_lumen", "Lumen, Apical, Inlet, Outlet, TripleJunction", "Advective-diffusive lumen concentration"],
        ["urea_membrane", "Membrane, Apical, Basolateral, TripleJunction", "Membrane diffusion and both exchange interfaces"],
        ["urea_inter", "Inter, Basolateral, OuterWall, Inlet, Outlet, TripleJunction", "Interstitial diffusion; deliberately excluded from Apical"],
    ], [31 * mm, 84 * mm, 57 * mm]),
    PageBreak(),
]

story += [section("3. Mesh subsets and the TripleJunction correction"),
    subset_diagram(),
    P("Figure 2. Compartment topology and the only two permitted exchange paths.", "Caption"),
    subsection("3.1 Why TripleJunction was needed"),
    P("At the curve where lumen, membrane, and interstitial topology meet, assigning every shared entity to either Apical or Basolateral creates ambiguity. A dedicated <b>TripleJunction</b> subset lets the mesh represent the geometric junction without falsely declaring a membrane transport interface."),
    P("TripleJunction participates in the support of all required fields, and it is included in the flow wall condition. It is not itself configured as an apical or basolateral leak. This distinction prevents unsupported degrees of freedom and preserves the intended physical connectivity."),
    subsection("3.2 Why urea_inter must not be on Apical"),
    P("The apical interface separates lumen from membrane. Interstitial urea has no direct apical interaction. Adding urea_inter to Apical created a large defect because it introduced degrees of freedom or constraints on an interface where the interstitial compartment is not physically present. The correct exchanges are lumen <-> membrane on Apical and membrane <-> interstitium on Basolateral."),
    subsection("3.3 Boundary-subset repair"),
    P("Older coarse grids labeled one-sided terminal annuli as Basolateral. The current pipeline repairs these labels on rank 0 before refinement and again after every refinement level, before distribution. This keeps Basolateral restricted to a true two-sided membrane-interstitium interface and avoids unsupported interface DoFs."),
    PageBreak(),
]

story += [section("4. Projector-aware refinement and parallel consistency"),
    subsection("4.1 The original parallel failure"),
    P("The flow direction was correct with one process but some velocities pointed incorrectly in parallel. The main cause was not the Stokes equation: the curved mesh hierarchy was being refined after or during distribution in a way that made NeuriteProjector refinement partition-dependent. Different ranks could create or classify refined geometry inconsistently."),
    subsection("4.2 Current sequence"),
    Preformatted("1. Load once into a temporary domain to register projSH and serialized projector data.\n2. Load the real domain.\n3. Repair legacy terminal subset labels on rank 0.\n4. Set balancer.firstDistLvl = numRefs.\n5. Refine with GlobalDomainRefiner while the relevant geometry is still coherent.\n6. Repair newly inherited midpoint labels after every level.\n7. Rebalance/distribute the completed hierarchy.", styles["CodeX"]),
    P("The key is that projector-generated child geometry is created consistently before the hierarchy is partitioned. The parent-child hierarchy remains available to geometric multigrid. The serialized NeuriteProjector and its npSurfParams attachment must both be present; opening and re-saving a grid through a tool that drops the attachment can leave an apparently valid projector that cannot refine correctly."),
    subsection("4.3 What the projector does and does not do"),
    P("The projector places new refinement vertices onto the curved nephron geometry. It does not solve the flow and it does not determine the inlet velocity sign. Flow direction still comes from the inlet tangent vector and boundary conditions. Projector correctness is a geometric prerequisite for a stable multilevel discretization."),
    P("Source anchors in the reorganized file: projector registration and domain load around lines 60-71; first distribution level and explicit refiner around lines 86-103."),
    PageBreak(),
]

story += [section("5. How lumen flow is implemented"),
    subsection("5.1 Discretization and support"),
    P("The flow approximation uses first-order Lagrange functions for u, v, w, and p. NavierStokesFV1 is assembled only on Lumen. Stokes mode, Laplace form, viscosity, LPS upwinding, and field-based stabilization are enabled. The source anchors are lines 115-127."),
    subsection("5.2 Parabolic inlet profile"),
    P("The inlet function uses the inlet center and a unit centerline tangent pointing into the lumen. For a point x, the displacement d = x - x_center is split into axial and radial parts. With a = d dot t, the squared radial distance is |d|^2 - a^2. The profile is"),
    Preformatted("u(x) = Umax * max(0, 1 - r^2/R^2) * t", styles["CodeX"]),
    P("where R = 1e-3 mm, mean speed is 1.5 mm/s, and peak speed is 3.0 mm/s. Points at the shared inlet/wall rim are forced exactly to zero within a 1e-6 relative-radius tolerance to eliminate round-off leakage."),
    subsection("5.3 Outlet and wall conditions"),
    P("Inlet uses NavierStokesInflow. Outlet uses NavierStokesNoNormalStressOutflow. Apical and TripleJunction use NavierStokesWall. Including TripleJunction in the wall closes the no-slip boundary at the geometric seam."),
    subsection("5.4 Pressure reference"),
    P("A pressure reference point was considered as a way to control pressure scaling. It was not adopted in the working configuration. The no-normal-stress outlet provides the natural pressure reference for this open-flow problem. A pointwise pressure Dirichlet condition would address a pure pressure nullspace, but it would not repair a reversed tangent, inconsistent projector refinement, or a bad multigrid hierarchy."),
    PageBreak(),
]

story += [section("6. Flow geometric multigrid"),
    subsection("6.1 Hierarchy configuration"),
    table([
        ["Element", "Current value", "Reason"],
        ["Outer linear solver", "BiCGStab", "Handles the nonsymmetric/preconditioned saddle-point system."],
        ["Preconditioner", "GeometricMultiGrid", "Scales to the intended large refined mesh."],
        ["Level smoother", "ILU, damping 0.7", "Coupled smoothing is more robust than scalar Gauss-Seidel."],
        ["Parallel interfaces", "Overlap + consistent interfaces", "Required for coherent distributed smoothing."],
        ["Cycle", "V, 4 pre + 4 post", "Strong smoothing around each coarse correction."],
        ["Coarse operator", "RAP false", "Rediscretizes Stokes on coarse levels; Galerkin RAP was unstable."],
        ["Level-0 solver", "AgglomeratingSolver(SuperLU)", "Direct pivoting solve for the saddle-point coarse problem."],
        ["Transfer", "StdTransfer", "Uses the preserved geometric hierarchy."],
    ], [39 * mm, 48 * mm, 85 * mm]),
    subsection("6.2 Why RAP is false for flow"),
    P("The successful flow configuration rediscretizes the Stokes operator on coarse levels. Earlier scalar smoothing with a Galerkin RAP hierarchy diverged at higher refinement. For this FV1 saddle-point problem, the combination of rediscretization, coupled ILU, and a direct level-0 solver was materially more robust."),
    subsection("6.3 Why pressure averaging was not used"),
    P("AverageComponent(\"p\") and similar all-subset pressure operations were unsuitable because pressure is intentionally undefined outside the lumen support. Applying such a helper over the mixed domain can request pressure DoFs where none exist. The open outlet condition is therefore used without a global all-subset pressure operation."),
    PageBreak(),
]

story += [section("7. Flow convergence evidence"),
    P("The values below are taken from recorded terminal output, not estimated from plots."),
    table([
        ["Run", "Initial defect", "Linear result", "Newton result", "Status"],
        ["Serial, ref0, stationary", "6.615282e-4", "2 BiCGStab steps; 3.826123e-18", "1 step; 9.147669e-15", "PASS"],
        ["MPI-4, ref1, stationary", "4.441660e-4", "21 steps; 2.129542e-12; avg rate 0.4016473", "1 step; 2.129588e-12; relative 4.794575e-9", "PASS"],
        ["MPI-4, older repeated-flow step", "2.220830e-4", "21 steps; 1.147174e-12", "1 step; relative 5.165513e-9", "PASS for flow"],
        ["Older ref2 flow hierarchy", "about 1.19e-5", "Residual increased and oscillated", "No acceptable solve", "FAIL"],
    ], [31 * mm, 28 * mm, 52 * mm, 48 * mm, 18 * mm]),
    P("The MPI ref1 stationary solve required more iterations than serial ref0, but it satisfied the configured relative reduction of 1e-8. Its recorded runtime was 0.022417 minutes. The earlier ref2 failure motivated the projector/distribution and flow-GMG corrections; it is historical failure evidence, not evidence that the current ref2 configuration has passed."),
    P("The short post-refactor validation reproduced the serial ref0 result exactly: BiCGStab reached 3.826123e-18 in two iterations and Newton ended at 9.147669e-15. This confirms that the structural reformat did not change the stationary flow path."),
    PageBreak(),
]

story += [section("8. Flow defects encountered and their corrections"),
    table([
        ["Observed problem", "Cause", "Correction", "Outcome"],
        ["Wrong directions only in MPI", "Projector refinement was partition-sensitive after early distribution", "Refine before distribution using firstDistLvl=numRefs and explicit hierarchy handling", "MPI ref1 flow passed"],
        ["Negative urea concentration at ref0", "Velocity direction/field usage was wrong", "Correct tangent and supply the solved vector field", "Negative concentration issue disappeared"],
        ["Wall leak at seam", "Triple junction was omitted from wall BC", "Wall on Apical,TripleJunction", "Closed no-slip boundary"],
        ["Ref2 GMG divergence", "Scalar/Galerkin coarse treatment was weak for Stokes", "Damped coupled ILU, RAP false, SuperLU coarse solve", "Lower tested levels passed; current ref2 still needs a fresh proof"],
        ["Pressure-reference concern", "Pressure scale was suspected", "Retain no-normal-stress outlet; do not mask geometry/solver defects with a point constraint", "Working open-flow gauge"],
    ], [35 * mm, 47 * mm, 55 * mm, 35 * mm]),
    subsection("8.1 Important distinction"),
    P("A solution can converge algebraically and still be physically wrong. The MPI direction problem was discovered by inspecting velocity orientation, not by solver failure alone. The final validation therefore needs both convergence traces and ParaView checks of axial direction, wall no-slip behavior, and absence of disconnected or reversed patches."),
    PageBreak(),
]

story += [section("9. Passing the solved velocity into transport"),
    subsection("9.1 The API error"),
    P("The first attempt called <b>set_velocity({\"u\", \"v\", \"w\"})</b>. UG4 rejected this because ConvectionDiffusionFV1 does not accept a Lua table of component names. It accepts a vector user-data object, a numerical vector, or a registered string form."),
    subsection("9.2 The valid implementation"),
    Preformatted("local solvedFlowVelocity = ExplicitGridFunctionVector(flowSolution, \"u,v,w\")\nlumenTransport:set_velocity(solvedFlowVelocity)", styles["CodeX"]),
    P("ExplicitGridFunctionVector creates the coupled vector-valued user data expected by ConvectionDiffusionFV1 while retaining access to flowSolution. In the current design flowSolution is solved once and never altered afterward, so every LIMEX stage sees the same stationary velocity."),
    subsection("9.3 Spatial stabilization"),
    P("FullUpwind is applied to lumen transport. This adds numerical robustness for convection-dominated cells and reduces nonphysical undershoot. It cannot compensate for an incorrect velocity orientation, which is why the direction fix had to precede transport tuning."),
    PageBreak(),
]

story += [section("10. Transport coupling and ParaView output"),
    subsection("10.1 Interface coupling"),
    P("The apical Leak couples only urea_lumen and urea_membrane. The basolateral Leak couples only urea_membrane and urea_inter. The code explicitly documents that no urea_inter/Apical coupling exists. This matches the compartment topology and prevents the large defect observed when Apical was added to interstitial support incorrectly."),
    subsection("10.2 VTK failure and resolution"),
    P("UG4 raised 'function component 6 has 0 DoFs in a vertex of an element of subset 7' when a mixed grid function was written over subsets where a selected component did not exist. VTK nodal output requires exactly one DoF per selected component at every written vertex."),
    P("The resolution is three independent VTKOutput selections and three subset-specific writes: urea_lumen on Lumen, urea_membrane on Membrane, and urea_inter on Inter. The output array is Urea_normalized. Separate PVD collections provide the complete time series for each physical compartment without fabricating unsupported values."),
    subsection("10.3 Exact output times"),
    P("The outer loop advances LIMEX separately from one requested output time to the next. Adaptive substeps occur inside each interval, but the returned solution lands exactly at 0, 1, 2, ..., 60 s. Resetting the initial step at each interval is predictable for output but may be less efficient than carrying the last accepted step size across intervals."),
    PageBreak(),
]

story += [section("11. Attempts at coupled flow and transport"),
    subsection("11.1 Re-solving flow every transport step"),
    P("An earlier implementation advanced or solved Stokes at every transport step. The flow portion could converge, but transport at a 0.5 s step became unstable. A representative flow step converged in 21 BiCGStab iterations to about 1.147e-12, while the following transport solve produced large oscillatory residuals and failed."),
    P("Because the geometry, viscosity, and flow boundary data are stationary, re-solving Stokes returns the same field and adds cost without adding physical coupling. It also entangles two solver failures in one time loop and makes diagnosis harder."),
    subsection("11.2 Monolithic versus one-way coupling"),
    P("The current implementation is <b>not monolithic coupling</b>. Urea does not modify density, viscosity, pressure, or momentum, so there is no feedback from transport to flow. It is mathematically consistent to solve Stokes once and use that field as a fixed coefficient in transport. If future physics adds concentration-dependent viscosity or osmotic flow, a true coupled or staggered nonlinear solve would again be necessary."),
    subsection("11.3 What worked"),
    *bullets([
        "Stationary Stokes solve followed by ExplicitGridFunctionVector.",
        "Fixed-flow backward-Euler transport at sufficiently small steps.",
        "Short-horizon adaptive LIMEX with a subset-aware global error estimator.",
    ]),
    subsection("11.4 What did not work"),
    *bullets([
        "Repeated flow solves combined with a large 0.5 s transport step.",
        "Treating the problem as fully coupled when there was no transport-to-flow feedback.",
        "Using a generic Newton loop inside linearly implicit LIMEX stages; it repeated a frozen linearization instead of one intended correction.",
    ]),
    PageBreak(),
]

story += [section("12. Adaptive LIMEX design"),
    subsection("12.1 Two-stage estimate"),
    P("LimexTimeIntegrator(2) compares one backward-Euler step with two half steps. Their difference estimates local temporal error. The configured tolerance is 1e-3 in the normalized concentration scale, with initial dt 1e-3 s, dt_min 1e-12 s, dt_max 1 s, safety factor 0.25, reduction factor 0.5, increase factor 2.0, and at most ten reductions."),
    subsection("12.2 Error-estimator problems"),
    table([
        ["Estimator attempt", "Failure", "Resolution"],
        ["Component L2 over mixed space", "Requested components on subsets with zero DoFs", "Use subset-specific evaluators"],
        ["Composite relative error", "Initially zero membrane/inter fields produced NaN or infinity", "Use absolute normalized scale"],
        ["Algebraic NormInf", "Appeared correct in serial but was rank-local in MPI", "Use GridFunctionEstimator + SupErrorEvaluator"],
        ["Current estimator", "Global compartment maximum", "One SupErrorEvaluator per volume subset"],
    ], [39 * mm, 65 * mm, 68 * mm]),
    subsection("12.3 Matrix caching"),
    P("Matrix caching is enabled because diffusion, permeability, and the fixed velocity field do not change with time. This avoids unnecessary reassembly. It is safe for the current linear, fixed-coefficient transport model; it would need reassessment if coefficients became concentration- or time-dependent."),
    PageBreak(),
]

story += [section("13. Transport convergence evidence"),
    table([
        ["Test", "Key values", "Status", "Interpretation"],
        ["Fixed BE, MPI-4 ref1, dt=1e-3", "Initial 2.708517e-9; about 12 BiCGStab steps to 1.10e-11", "PASS under relaxed absolute check", "Small fixed step can be solved"],
        ["LIMEX serial ref0, horizon 4e-6", "Temporal error 3.836532e-4 < 1e-3", "PASS", "Short adaptive path accepted"],
        ["LIMEX MPI-4 ref1, proposed dt=4e-6", "Global error 4.948570e-3", "REJECTED correctly", "Estimator detected excess temporal error"],
        ["LIMEX MPI-4 ref1, dt=8.99063e-7", "Error 3.327e-4", "PASS", "Reduced step accepted; reached 4e-6 horizon"],
        ["Full 0-60 s start, MPI-4 ref1, BiCGStab", "Stage initial 2.708517e-9; iteration 200 ended 1.840953e-12", "FAIL", "Did not satisfy absolute 1e-14 or relative 1e-8"],
        ["Relaxed BiCGStab trial", "First full step passed in 12; first half passed in 6; second half failed by 100", "FAIL", "Tolerance relaxation alone is insufficient"],
        ["Current GMRES, tiny ref0 validation", "Three stage solves each reported convergence in one outer iteration", "PASS only as smoke test", "Does not prove refined 60 s run"],
    ], [36 * mm, 61 * mm, 29 * mm, 46 * mm]),
    P("The full-run BiCGStab trace is especially diagnostic: it briefly reached about 1.10e-11, then stagnated and showed large spikes. A residual floor several orders above the requested relative target means step reduction alone may not cure the Krylov/preconditioner interaction."),
    PageBreak(),
]

story += [section("14. Solver pass/fail matrix"),
    table([
        ["Solver or configuration", "Problem block", "Result", "Conclusion"],
        ["BiCGStab + flow GMG, RAP false", "Stationary Stokes", "PASS: serial ref0 and MPI-4 ref1", "Current flow baseline"],
        ["Scalar smoother + Galerkin/RAP hierarchy", "Refined Stokes", "FAIL", "Unsuitable coarse treatment"],
        ["SuperLU as level-0 flow solver", "Stokes coarse grid", "PASS", "Needed for robust saddle-point base solve"],
        ["Generic Newton in LIMEX", "Linearly implicit transport stage", "FAIL", "Wrong stage algorithm"],
        ["LimexNewtonSolver", "Transport stage", "PASS on short tests", "Correct one-correction stage solver"],
        ["BiCGStab + transport GMG", "Long/refined LIMEX stages", "FAIL", "Residual spikes and stagnation"],
        ["GMRES + transport GMG", "Current transport configuration", "UNPROVEN", "Tiny ref0 smoke test passed; full refined run pending"],
        ["Subset-aware SupErrorEvaluator", "MPI temporal error", "PASS", "Correct global compartment error"],
    ], [48 * mm, 43 * mm, 37 * mm, 44 * mm]),
    subsection("14.1 Flow and transport GMG differ deliberately"),
    P("Flow uses RAP false because rediscretized coarse Stokes operators were more robust. Transport uses RAP true so coarse corrections inherit the assembled mixed-compartment fine operator. Applying one RAP choice blindly to both blocks was not justified; their algebraic structures are different."),
    PageBreak(),
]

story += [section("15. Why the full transport run is not fully solving"),
    P("The remaining failure is not evidence that diffusion physics is impossible to solve. It is a linear-algebra robustness problem that appears inside some adaptive stage systems on the refined distributed mesh."),
    *bullets([
        "The transport operator is nonsymmetric because lumen convection uses upwinding.",
        "The mixed support of three concentration fields and two membrane interfaces produces a coupled block structure.",
        "Adaptive LIMEX generates several related stage matrices; a preconditioner adequate for one stage may be weak for another.",
        "BiCGStab is vulnerable to loss of bi-orthogonality and breakdown-like residual spikes on such systems.",
        "An absolute requirement below the attainable distributed residual floor can label a physically adequate solve as failure, but simply relaxing it did not cure later-stage breakdown.",
        "The 60-second request requires many successful stage solves, so a solver that passes only the first stage is not sufficient.",
    ]),
    P("The current switch to restarted GMRES is a reasonable next candidate because GMRES minimizes the residual over a Krylov subspace and is often less erratic than BiCGStab. However, the implementation's GMRES restart length and transport-GMG quality still need testing on MPI ref1/ref2 over meaningful time horizons."),
    P("It would be inaccurate to report the 0-60 s simulation as converged until the code reaches 60 s, writes all requested checkpoints, and the concentration fields pass physical checks."),
    PageBreak(),
]

story += [section("16. Recommended validation sequence"),
    P("The next experiments should isolate solver scaling before spending time on a full production run."),
    table([
        ["Step", "Run", "Acceptance evidence"],
        ["1", "Serial ref0 to 1 s with current GMRES", "Every LIMEX stage converges; t=1 output exists; C* remains bounded"],
        ["2", "MPI-4 ref1 to 0.01 s", "No rank-specific failure; accepted/rejected steps are plausible"],
        ["3", "MPI-4 ref1 to 1 s", "Checkpoint 1 exists; residuals do not stagnate repeatedly"],
        ["4", "MPI-4 ref1 to 60 s", "All 61 frames exist; no negative/overshoot concentrations beyond justified tolerance"],
        ["5", "MPI-4 ref2 short horizon", "Flow and transport both converge before extending runtime"],
        ["6", "Refinement comparison", "Fluxes and concentration profiles approach a mesh-independent result"],
    ], [13 * mm, 68 * mm, 91 * mm]),
    subsection("16.1 Diagnostics to retain"),
    *bullets([
        "Verbose linear convergence traces for flow and transport.",
        "Accepted/rejected LIMEX step sizes and error estimates.",
        "Minimum and maximum of each normalized concentration after every output interval.",
        "Integrated inlet/outlet advective flux and apical/basolateral exchange flux for mass balance.",
        "ParaView arrows or streamlines confirming velocity always follows the nephron lumen.",
    ]),
    subsection("16.2 If GMRES still fails"),
    P("The next investigation should focus on the transport preconditioner, not on reintroducing repeated flow solves. Candidate work includes a stronger block-aware smoother, flexible GMRES if the preconditioner varies, adjusted restart size, and separate examination of the three compartment/interface blocks. Any tolerance change should be justified against normalized concentration error and mass balance."),
    PageBreak(),
]

story += [section("17. Current source map and exact configuration"),
    table([
        ["Feature", "Current source location"],
        ["Direct run inputs", "Urea_trans_with_flow_coupled_gmg.lua:13-43"],
        ["Projector refinement pipeline", "lines 60-103"],
        ["Flow approximation/discretization", "lines 113-127"],
        ["Inlet profile and flow boundaries", "lines 129-160"],
        ["Flow GMG and solver", "lines 168-203"],
        ["Single stationary flow solve", "lines 205-221"],
        ["Transport field supports", "lines 226-230"],
        ["Transport PDE and membrane fluxes", "lines 233-266"],
        ["Transport GMG and GMRES", "lines 270-303"],
        ["LIMEX and global MPI error estimator", "lines 313-336"],
        ["Subset-specific VTK output", "lines 342-364"],
        ["Adaptive 1-second checkpoint loop", "lines 368-394"],
    ], [68 * mm, 104 * mm]),
    subsection("17.1 Current literal solver tolerances"),
    table([
        ["Block", "Iterations", "Absolute", "Relative"],
        ["Flow Newton", "50", "1e-12", "1e-8"],
        ["Flow BiCGStab", "150", "1e-14", "1e-8"],
        ["Transport GMRES", "100", "1e-30", "1e-8"],
        ["LIMEX temporal estimator", "adaptive", "1e-3 in C* sup norm", "not relative"],
    ], [48 * mm, 34 * mm, 45 * mm, 45 * mm]),
    P("The transport absolute linear threshold of 1e-30 intentionally leaves relative reduction as the effective stopping criterion. It should not be interpreted as an expectation that distributed arithmetic reaches a 1e-30 residual."),
    PageBreak(),
]

story += [section("18. Final assessment"),
    P("The most important achievement is the separation of concerns. Geometry is refined consistently with the NeuriteProjector before partitioning. Stokes flow is solved once with a flow-specific multigrid hierarchy. The velocity is passed through the correct UG4 vector-user-data interface. Transport uses physically correct compartment supports and membrane interfaces, normalized concentrations, subset-aware ParaView output, and an MPI-global adaptive error estimator."),
    P("The main unresolved item is the robustness of the refined long-time transport linear solver. Stationary flow is supported by clear serial and MPI convergence evidence. Short adaptive transport is supported by accepted steps and quantified error estimates. The full 60-second refined MPI transport trajectory is still pending and should remain labeled unvalidated until it completes with physically consistent outputs."),
    P("This distinction is central: the code architecture and flow implementation are now sound at the tested level, while long-horizon transport scalability remains an active numerical-solver task.", "Callout"),
    Spacer(8 * mm, 1),
    subsection("Evidence base"),
    P("Primary source: Urea_trans_with_flow_coupled_gmg.lua in the Model directory. Comparative source: Dif_itterative_flow_projector.lua. Recorded run evidence: the supplied serial/MPI terminal logs, including the MPI-4 ref1 stationary-flow pass and the 0-60 s LIMEX/BiCGStab failure trace. No claim in this report treats an unexecuted configuration as proven."),
]


OUT.parent.mkdir(parents=True, exist_ok=True)
doc = ReportDoc(
    str(OUT),
    pagesize=A4,
    title="Flow, Projector Refinement, and Urea Transport in UG4",
    author="Technical implementation summary",
    subject="UG4 nephron Stokes flow, NeuriteProjector refinement, urea transport, LIMEX, and GMG convergence",
    leftMargin=18 * mm,
    rightMargin=18 * mm,
    topMargin=18 * mm,
    bottomMargin=18 * mm,
)
doc.multiBuild(story)
print(OUT)
