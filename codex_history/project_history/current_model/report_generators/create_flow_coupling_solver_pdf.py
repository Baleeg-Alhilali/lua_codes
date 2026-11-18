from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


OUTPUT = "output/pdf/stokes_darcy_coupling_solver_explanation.pdf"

PAGE_W, PAGE_H = A4
NAVY = colors.HexColor("#17324D")
BLUE = colors.HexColor("#2B6F9F")
PALE_BLUE = colors.HexColor("#EAF3F8")
PALE_GRAY = colors.HexColor("#F4F6F7")
MID_GRAY = colors.HexColor("#687784")
RULE = colors.HexColor("#CBD6DD")


styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name="DocTitle",
    parent=styles["Title"],
    fontName="Helvetica-Bold",
    fontSize=23,
    leading=28,
    textColor=NAVY,
    alignment=TA_LEFT,
    spaceAfter=8,
))
styles.add(ParagraphStyle(
    name="Subtitle",
    parent=styles["Normal"],
    fontName="Helvetica",
    fontSize=11,
    leading=16,
    textColor=MID_GRAY,
    spaceAfter=20,
))
styles.add(ParagraphStyle(
    name="Section",
    parent=styles["Heading1"],
    fontName="Helvetica-Bold",
    fontSize=15,
    leading=19,
    textColor=NAVY,
    spaceBefore=13,
    spaceAfter=7,
    keepWithNext=True,
))
styles.add(ParagraphStyle(
    name="Subsection",
    parent=styles["Heading2"],
    fontName="Helvetica-Bold",
    fontSize=11.5,
    leading=15,
    textColor=BLUE,
    spaceBefore=9,
    spaceAfter=4,
    keepWithNext=True,
))
styles.add(ParagraphStyle(
    name="BodyFlow",
    parent=styles["BodyText"],
    fontName="Helvetica",
    fontSize=9.5,
    leading=14,
    textColor=colors.HexColor("#263640"),
    spaceAfter=7,
))
styles.add(ParagraphStyle(
    name="BulletFlow",
    parent=styles["BodyText"],
    fontName="Helvetica",
    fontSize=9.3,
    leading=13.4,
    leftIndent=13,
    firstLineIndent=-7,
    bulletIndent=4,
    textColor=colors.HexColor("#263640"),
    spaceAfter=3,
))
styles.add(ParagraphStyle(
    name="Equation",
    parent=styles["Code"],
    fontName="Courier",
    fontSize=8.6,
    leading=12,
    alignment=TA_CENTER,
    textColor=NAVY,
    borderColor=RULE,
    borderWidth=0.6,
    borderPadding=7,
    backColor=PALE_GRAY,
    spaceBefore=4,
    spaceAfter=9,
    keepWithNext=True,
))
styles.add(ParagraphStyle(
    name="Callout",
    parent=styles["BodyText"],
    fontName="Helvetica-Bold",
    fontSize=9.5,
    leading=14,
    textColor=NAVY,
    borderColor=BLUE,
    borderWidth=0.8,
    borderPadding=8,
    backColor=PALE_BLUE,
    spaceBefore=5,
    spaceAfter=10,
))
styles.add(ParagraphStyle(
    name="Small",
    parent=styles["BodyText"],
    fontName="Helvetica",
    fontSize=8,
    leading=11,
    textColor=MID_GRAY,
))


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(RULE)
    canvas.setLineWidth(0.5)
    canvas.line(18 * mm, 14 * mm, PAGE_W - 18 * mm, 14 * mm)
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(MID_GRAY)
    canvas.drawString(18 * mm, 9.5 * mm, "Stokes-Darcy flow coupling and solver setup")
    canvas.drawRightString(PAGE_W - 18 * mm, 9.5 * mm, f"Page {doc.page}")
    canvas.restoreState()


doc = BaseDocTemplate(
    OUTPUT,
    pagesize=A4,
    leftMargin=19 * mm,
    rightMargin=19 * mm,
    topMargin=18 * mm,
    bottomMargin=20 * mm,
    title="Stokes-Darcy Flow Coupling and Solver Setup",
    author="Promesh SWC Model",
    subject="Monolithic Stokes-Darcy-BJS formulation and GMG solver",
)
frame = Frame(
    doc.leftMargin,
    doc.bottomMargin,
    doc.width,
    doc.height,
    leftPadding=0,
    rightPadding=0,
    topPadding=0,
    bottomPadding=0,
)
doc.addPageTemplates(PageTemplate(id="main", frames=[frame], onPage=footer))


story = []


def p(text, style="BodyFlow"):
    story.append(Paragraph(text, styles[style]))


def bullet(text):
    story.append(Paragraph("- " + text, styles["BulletFlow"]))


def section(number, title):
    story.append(Paragraph(f"{number}. {title}", styles["Section"]))


def subsection(title):
    story.append(Paragraph(title, styles["Subsection"]))


def equation(text):
    story.append(Paragraph(text.replace("\n", "<br/>"), styles["Equation"]))


story.append(Spacer(1, 6 * mm))
p("Stokes-Darcy Flow Coupling", "DocTitle")
p("How the monolithic interface formulation and the BiCGStab/GMG solver work", "Subtitle")
p(
    "The current flow model solves the lumen Stokes velocity, lumen pressure, "
    "and porous Darcy pressure in one assembled algebraic system. It does not "
    "solve the lumen first and pass a frozen result to Darcy. Both interface "
    "directions remain active inside the same matrix.",
    "Callout",
)

summary_data = [
    [Paragraph("Component", styles["Small"]), Paragraph("Current choice", styles["Small"])],
    ["Flow model", "Steady Stokes in Lumen; Darcy in Membrane and Inter"],
    ["Interface", "Normal traction + normal flux continuity + BJS tangential slip"],
    ["Unknowns", "u, v, w, p_S, scaled p_D"],
    ["Outer solver", "BiCGStab"],
    ["Preconditioner", "Geometric multigrid with an Uzawa block smoother"],
    ["Coarse solve", "SuperLU on level 0 only"],
    ["Parallel use", "Distributed fine-grid solve; gathered coarse solve when needed"],
]
summary = Table(summary_data, colWidths=[40 * mm, 126 * mm], repeatRows=1)
summary.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), NAVY),
    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
    ("FONTNAME", (0, 1), (0, -1), "Helvetica-Bold"),
    ("FONTNAME", (1, 1), (1, -1), "Helvetica"),
    ("FONTSIZE", (0, 1), (-1, -1), 8.5),
    ("LEADING", (0, 0), (-1, -1), 11),
    ("BACKGROUND", (0, 1), (-1, -1), colors.white),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE_GRAY]),
    ("GRID", (0, 0), (-1, -1), 0.35, RULE),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("LEFTPADDING", (0, 0), (-1, -1), 6),
    ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ("TOPPADDING", (0, 0), (-1, -1), 5),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
]))
story.append(summary)
story.append(Spacer(1, 5 * mm))
p(
    "Implementation reference: code/flow/working_flow_coupling_simple_long_tube.lua. "
    "The curved-tube production script uses the same coupled formulation.",
    "Small",
)
story.append(PageBreak())


section(1, "Unknowns in the coupled system")
p("The global solution vector is")
equation("x = [ u, v, w, p_S, p_D_hat ]^T")
bullet("u, v, and w are the lumen Stokes velocity components.")
bullet("p_S is the lumen Stokes pressure.")
bullet("p_D_hat is the scaled Darcy-pressure unknown.")
p("The physical Darcy pressure is recovered from")
equation("p_D = S_p p_D_hat,     S_p = 10^4.")
p(
    "All five fields belong to one ApproximationSpace. Therefore, interface "
    "terms can contribute derivatives directly to the same global matrix."
)
p(
    "The pressure scaling reduces the large algebraic magnitude difference "
    "between the free-flow and porous-flow variables. To preserve the physical "
    "Darcy law, the diffusion coefficient is scaled by the same factor."
)
equation("M_D_hat = (K / nu) S_p\nu_D = -(K / nu) grad(p_D) = -(K / nu) S_p grad(p_D_hat)")
p(
    "Thus scaling changes the numerical unknown but does not change the "
    "physical Darcy velocity."
)

section(2, "Equations in each region")
subsection("Lumen: steady Stokes")
equation("-nu div(grad(u) + grad(u)^T) + grad(p_S) = 0\ndiv(u) = 0")
p(
    "The current parameter bstokes is true, so the convective Navier-Stokes "
    "term is disabled. The complete coupled problem is therefore linear."
)
subsection("Membrane and Inter: Darcy pressure")
equation("-div((K / nu) grad(p_D)) = 0")
p("The current permeabilities are:")
bullet("Membrane: K_M = 10^-15 m^2.")
bullet("Interstitium: K_I = 10^-14 m^2.")
p(
    "The same Darcy-pressure field spans Membrane and Inter. FV1 assembly "
    "therefore provides the pressure and conservative-flux transmission across "
    "the Basolateral material interface."
)


section(3, "Coupling at the Apical interface")
subsection("Normal stress and Darcy-pressure coupling")
p("Darcy pressure supplies the normal traction on the lumen:")
equation("sigma_S n_S = -p_D n_S")
p(
    "The code inserts p_D into the x, y, and z Stokes momentum rows as Neumann "
    "traction contributions. This is not literal equality p_S = p_D. It is a "
    "relationship between Darcy pressure and the full normal Stokes stress."
)

subsection("Normal velocity continuity")
equation("u_S . n_S = u_D . n_S     (with the mesh-orientation sign)")
p(
    "The Stokes velocity exported by stokesEq:velocity() becomes the Darcy "
    "boundary flux on the Membrane side of Apical. This exported value keeps "
    "its Jacobian derivatives: it is not a frozen velocity from a previous "
    "solve. Darcy pressure likewise enters the Stokes traction equation as "
    "linked data. The coupling therefore acts in both directions during the "
    "same linear solve."
)
equation("[ A_S   B_S^T  C_SD ] [ u_S     ]   [ f_S ]\n[ B_S     0       0   ] [ p_S     ] = [  0  ]\n[ C_DS    0      A_D  ] [ p_D_hat ]   [ f_D ]")
p(
    "C_SD transfers Darcy pressure into the Stokes traction, while C_DS "
    "transfers Stokes normal velocity into the Darcy conservation equation."
)

subsection("Tangential Beavers-Joseph-Saffman condition")
equation("(sigma_S n)_t = -[alpha nu / sqrt(K_BJS)] u_S,t")
p(
    "NavierStokesWSBCFV1 assembles this tangential friction law. Darcy velocity "
    "is a porous-volume average, so the model does not equate the complete "
    "Stokes and Darcy velocity vectors."
)


section(4, "Normal membrane resistance")
p("The interface also includes the normal Robin resistance")
equation("R_n = delta nu / K_M")
p(
    "Physically, this represents the effective resistance of a membrane of "
    "thickness delta. Numerically, it removes unresolved vertex-scale normal "
    "interface modes that previously appeared as isolated leakage jets."
)


section(5, "Boundary conditions that close the system")
bullet("A parabolic velocity profile is imposed at Inlet.")
bullet("A no-normal-stress condition is applied at Outlet.")
bullet("The open-boundary discretization is also present on Apical to close the Stokes continuity control volumes with the actual boundary-normal velocity.")
bullet("The BJS condition is applied on Apical.")
bullet("A Stokes wall condition is applied at TripleJunction.")
bullet("Darcy pressure is prescribed at OuterWall.")
bullet("Porous boundaries without an explicit condition retain the natural no-flux condition.")
p(
    "Unused field extensions are pinned only outside their physical domains. "
    "The Apical trace is never pinned because both Stokes and Darcy variables "
    "must remain active there for monolithic coupling."
)


section(6, "Why the outer solver is BiCGStab")
p("The assembled system is:")
bullet("indefinite because of incompressible Stokes;")
bullet("nonsymmetric because of the directional interface terms;")
bullet("strongly heterogeneous because the Darcy permeabilities differ from the lumen scales; and")
bullet("too large for a full direct solve at the intended refinements.")
p(
    "BiCGStab is therefore applied to the entire five-field matrix. It is not "
    "run independently for Stokes and Darcy."
)


section(7, "Geometric multigrid preconditioner")
p("The outer BiCGStab iteration is preconditioned by a geometric V-cycle:")
equation("2 pre-smooths -> restrict residual -> coarse correction\n-> prolong correction -> 2 post-smooths")
p(
    "GMG removes short-wavelength error by smoothing on each level and "
    "long-wavelength error through coarse-grid correction. RAP is enabled, so "
    "the coarse operator is constructed as"
)
equation("A_H = R A_h P")
p(
    "This is important for the coupled problem because the Stokes-Darcy "
    "off-diagonal blocks are inherited by coarse levels. The coarse system "
    "therefore continues to represent the interface coupling."
)


section(8, "Uzawa block smoother")
p(
    "A generic ILU smoother is unreliable for incompressible Stokes because "
    "the pressure-pressure block is essentially zero. The GMG smoother is "
    "therefore an Uzawa block iteration."
)
p("Only p_S is classified as the saddle-point pressure:")
equation("Main block:  [u, v, w, p_D_hat]\nSchur block: [p_S]")
p(
    "Keeping p_D_hat in the main block is essential. Darcy pressure has an "
    "elliptic diffusion operator and carries the interface coupling. Treating "
    "both pressures as Schur variables destroys this useful structure."
)
p("One smoothing operation approximately performs:")
bullet("A damped ILU correction for the main velocity-Darcy block.")
bullet("A Stokes-pressure correction using an approximate Schur complement.")
bullet("A second damped ILU correction for the main block.")
p(
    "The Stokes Schur complement is approximated by a lumen pressure-mass "
    "operator scaled by 1/nu."
)


section(9, "Coarse-grid solve and MPI")
p(
    "At level zero, the small coupled saddle-point system is solved with "
    "SuperLU. SuperLU is not applied to the complete refined problem. It is "
    "restricted to the coarse level, where pivoting gives a reliable global "
    "correction."
)
equation("BiCGStab + GMG + Uzawa/ILU smoothing + coarse SuperLU")
p(
    "The fine-grid matrix, vectors, smoothing, and Krylov iteration remain "
    "distributed under mpirun. When required, only the small coarse system is "
    "gathered for the direct solve. The projector-backed hierarchy is created "
    "consistently before distribution so that refinement preserves the curved "
    "tube and interface geometry."
)


story.append(PageBreak())
section(10, "Why the system is solvable now")
bullet("The free-flow and porous-flow equations are assembled monolithically.")
bullet("Both directions of the Apical coupling are present in the matrix.")
bullet("The Apical Stokes continuity control volumes are properly closed.")
bullet("BJS supplies the tangential interface law.")
bullet("Normal membrane resistance suppresses nonphysical leakage modes.")
bullet("TripleJunction has an explicit wall condition.")
bullet("OuterWall supplies a unique Darcy-pressure reference.")
bullet("Darcy-pressure scaling reduces the algebraic scale imbalance.")
bullet("Only Stokes pressure is treated as the saddle-point variable.")
bullet("RAP preserves the coupled blocks on coarse GMG levels.")
bullet("The small coarse saddle-point problem receives a pivoted direct solve.")
p(
    "These pieces address different failure modes. No single solver switch "
    "would have corrected missing interface equations, pressure nullspaces, "
    "poor scaling, unresolved normal modes, and an unsuitable block smoother "
    "at the same time.",
    "Callout",
)


section(11, "Convergence verification")
p(
    "Each BiCGStab cycle permits up to 400 iterations, with an internal relative "
    "reduction of 10^-4 and an absolute threshold of 10^-14. The code does not "
    "accept that combined convergence report by itself."
)
p("After each cycle, it explicitly recomputes the true residual")
equation("r = b - A x")
p("and prints separate norms for:")
bullet("u, v, and w;")
bullet("p_S; and")
bullet("p_D_hat.")
p(
    "BiCGStab is restarted up to eight times until every component meets the "
    "verified componentwise tolerance. This prevents one large or differently "
    "scaled field from hiding a poorly converged component."
)
p("The final physical checks report:")
bullet("inlet, outlet, and Apical flow;")
bullet("lumen mass-balance error;")
bullet("outlet/leakage split;")
bullet("Darcy interface-flow diagnostics; and")
bullet("Stokes and Darcy pressure diagnostics.")
p(
    "Algebraic convergence means that the discrete equations have been solved. "
    "Physical credibility still requires reviewing conservation diagnostics, "
    "the flow split, pressure levels, and the VTK fields.",
    "Callout",
)


doc.build(story)
print(OUTPUT)
