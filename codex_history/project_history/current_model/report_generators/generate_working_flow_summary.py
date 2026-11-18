#!/usr/bin/env python3
"""Create a concise 3-page boundary-and-solver summary for the long-tube model."""

from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer, PageBreak,
    Table, TableStyle, Flowable,
)

ROOT = Path('/Users/alhilaba/Documents/Promesh_swc/Model')
OUT = ROOT / 'output/pdf/working_flow_boundary_solver_summary.pdf'
W, H = A4
ML = MR = 17 * mm
MT = 17 * mm
MB = 16 * mm
CW = W - ML - MR

NAVY = colors.HexColor('#17324D')
BLUE = colors.HexColor('#2563A6')
CYAN = colors.HexColor('#1499A8')
GOLD = colors.HexColor('#C98718')
RED = colors.HexColor('#B7473E')
GREEN = colors.HexColor('#2C7A5A')
INK = colors.HexColor('#17212B')
MUTED = colors.HexColor('#5C6974')
LINE = colors.HexColor('#CCD7E0')
LIGHT = colors.HexColor('#F5F8FA')
PBLUE = colors.HexColor('#EAF3FA')
PCYAN = colors.HexColor('#E9F8F8')
PGOLD = colors.HexColor('#FFF5DC')
PRED = colors.HexColor('#FCECE9')
PGREEN = colors.HexColor('#EAF6EF')

SANS = 'Helvetica'
BOLD = 'Helvetica-Bold'
MONO = 'Courier'

ss = getSampleStyleSheet()
ss.add(ParagraphStyle(name='K', fontName=BOLD, fontSize=7.6, leading=9,
                      textColor=CYAN, tracking=1.1, spaceAfter=3))
ss.add(ParagraphStyle(name='T', fontName=BOLD, fontSize=20, leading=23,
                      textColor=NAVY, spaceAfter=6))
ss.add(ParagraphStyle(name='S', fontName=BOLD, fontSize=11.5, leading=14,
                      textColor=BLUE, spaceBefore=5, spaceAfter=4))
ss.add(ParagraphStyle(name='B', fontName=SANS, fontSize=8.3, leading=11.4,
                      textColor=INK, spaceAfter=4))
ss.add(ParagraphStyle(name='Sm', fontName=SANS, fontSize=6.8, leading=8.5,
                      textColor=INK))
ss.add(ParagraphStyle(name='TH', fontName=BOLD, fontSize=6.7, leading=8.2,
                      textColor=colors.white))
ss.add(ParagraphStyle(name='Eq', fontName=MONO, fontSize=7.5, leading=10,
                      textColor=NAVY, alignment=1, spaceBefore=3, spaceAfter=4))
ss.add(ParagraphStyle(name='CallT', fontName=BOLD, fontSize=8.2, leading=10,
                      textColor=NAVY, spaceAfter=2))
ss.add(ParagraphStyle(name='CallB', fontName=SANS, fontSize=7.3, leading=9.8,
                      textColor=INK))


def p(text, style='B'):
    return Paragraph(text, ss[style])


def heading(kicker, title, intro):
    return [p(kicker.upper(), 'K'), p(title, 'T'), p(intro, 'B'), Spacer(1, 1.5*mm)]


def table(rows, widths, small=True):
    body_style = 'Sm' if small else 'B'
    data = [[p(str(x), 'TH' if r == 0 else body_style) for x in row]
            for r, row in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1, hAlign='LEFT')
    cmds = [
        ('BACKGROUND', (0,0), (-1,0), NAVY),
        ('GRID', (0,0), (-1,-1), .35, LINE),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('LEFTPADDING', (0,0), (-1,-1), 4),
        ('RIGHTPADDING', (0,0), (-1,-1), 4),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ]
    for r in range(2, len(rows), 2):
        cmds.append(('BACKGROUND', (0,r), (-1,r), LIGHT))
    t.setStyle(TableStyle(cmds))
    return t


def callout(title, text, fill, stroke):
    t = Table([[[p(title, 'CallT'), p(text, 'CallB')]]], colWidths=[CW])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), fill),
        ('BOX', (0,0), (-1,-1), .55, stroke),
        ('LINEBEFORE', (0,0), (0,-1), 4, stroke),
        ('LEFTPADDING', (0,0), (-1,-1), 9),
        ('RIGHTPADDING', (0,0), (-1,-1), 8),
        ('TOPPADDING', (0,0), (-1,-1), 6),
        ('BOTTOMPADDING', (0,0), (-1,-1), 6),
    ]))
    return t


class BoundarySketch(Flowable):
    def __init__(self):
        super().__init__()
        self.width, self.height = CW, 47*mm

    def draw(self):
        c = self.canv
        x, y = 13*mm, 8*mm
        ww, hh = CW-26*mm, 25*mm
        wall = 7*mm
        c.setFillColor(PGOLD); c.setStrokeColor(GOLD)
        c.roundRect(x, y, ww, hh, 4, fill=1, stroke=1)
        c.setFillColor(PBLUE); c.setStrokeColor(BLUE)
        c.roundRect(x, y+wall/2, ww, hh-wall, 3, fill=1, stroke=1)
        c.setStrokeColor(CYAN); c.setLineWidth(2)
        c.line(x, y+wall/2, x+ww, y+wall/2)
        c.line(x, y+hh-wall/2, x+ww, y+hh-wall/2)
        c.setFillColor(NAVY); c.setFont(BOLD, 8)
        c.drawCentredString(x+ww/2, y+hh/2-3, 'Lumen: Stokes FV1')
        c.setFillColor(GOLD); c.setFont(BOLD, 6.8)
        c.drawCentredString(x+ww/2, y+hh-7, 'Inter: Darcy FV1')
        c.drawCentredString(x+ww/2, y+1, 'Inter: Darcy FV1')
        c.setFillColor(CYAN); c.drawCentredString(x+ww/2, y+hh+10, 'Apical: BJS + normal Robin; transferred mean flux')
        c.setFillColor(MUTED); c.setFont(SANS, 6.7)
        c.drawString(x-1, y+hh+17, 'Inlet: prescribed parabola')
        c.drawRightString(x+ww+1, y+hh+17, 'Outlet: open stress')
        c.drawCentredString(x+ww/2, y-6, 'OuterWall: p_D = p_outer')


class CouplingSketch(Flowable):
    def __init__(self):
        super().__init__()
        self.width, self.height = CW, 39*mm

    def box(self, c, x, y, w, h, title, detail, fill, stroke):
        c.setFillColor(fill); c.setStrokeColor(stroke)
        c.roundRect(x,y,w,h,4,fill=1,stroke=1)
        c.setFillColor(stroke); c.setFont(BOLD,7.5); c.drawString(x+6,y+h-11,title)
        c.setFillColor(INK); c.setFont(SANS,6.3)
        for i, s in enumerate(detail): c.drawString(x+6,y+h-20-i*8,s)

    def arrow(self,c,x1,y,x2,label):
        c.setStrokeColor(NAVY); c.setFillColor(NAVY)
        c.line(x1,y,x2,y); c.line(x2,y,x2-4,y+2); c.line(x2,y,x2-4,y-2)
        c.setFont(BOLD,6); c.drawCentredString((x1+x2)/2,y+5,label)

    def draw(self):
        c=self.canv; y=8*mm; bw=48*mm; bh=22*mm; gap=10*mm
        x1=1*mm; x2=x1+bw+gap; x3=x2+bw+gap
        self.box(c,x1,y,bw,bh,'1. STOKES',['Robin and BJS implicit','solve u and p'],PBLUE,BLUE)
        self.box(c,x2,y,bw,bh,'2. REDUCE',['integrate Q_S','q_D = -Q_S / A'],PCYAN,CYAN)
        self.box(c,x3,y,bw,bh,'3. DARCY',['uniform Neumann flux','solve p_D; derive u_D'],PGOLD,GOLD)
        self.arrow(c,x1+bw,y+bh/2,x2,'Q_S')
        self.arrow(c,x2+bw,y+bh/2,x3,'q_D')
        c.setFillColor(RED); c.setFont(BOLD,6.8)
        c.drawCentredString(CW/2,1.5*mm,'No computed Darcy pressure is iterated back to Stokes.')


class SolverSketch(Flowable):
    def __init__(self):
        super().__init__()
        self.width, self.height = CW, 58*mm

    def box(self,c,x,y,w,h,title,detail,fill,stroke):
        c.setFillColor(fill); c.setStrokeColor(stroke)
        c.roundRect(x,y,w,h,4,fill=1,stroke=1)
        c.setFillColor(stroke); c.setFont(BOLD,7.8); c.drawString(x+7,y+h-11,title)
        c.setFillColor(INK); c.setFont(SANS,6.5)
        for i,s in enumerate(detail): c.drawString(x+7,y+h-20-i*8,s)

    def down(self,c,x,y1,y2):
        c.setStrokeColor(NAVY); c.line(x,y1,x,y2); c.line(x,y2,x-3,y2+4); c.line(x,y2,x+3,y2+4)

    def draw(self):
        c=self.canv; cx=CW/2
        self.box(c,cx-50*mm,43*mm,100*mm,13*mm,'OUTER: BiCGStab',['1000 iters; abs 1e-14; reduction 1e-8'],PBLUE,BLUE)
        self.down(c,cx,43*mm,39*mm)
        self.box(c,cx-66*mm,23*mm,132*mm,16*mm,'PRECONDITIONER: GMG V-cycle',['4 ILU pre + 4 ILU post; damping 0.7; overlap','rap=false; base level 0; standard transfer'],PCYAN,CYAN)
        self.down(c,cx,23*mm,19*mm)
        self.box(c,cx-55*mm,2*mm,110*mm,17*mm,'COARSE: AgglomeratingSolver',['BiCGStab: 500; abs 1e-14; reduction 1e-10','ILUT drop 1e-8; Cuthill-McKee'],PGOLD,GOLD)


def hf(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(LINE); canvas.setLineWidth(.4)
    canvas.line(ML, H-11*mm, W-MR, H-11*mm)
    canvas.setFont(BOLD, 6.8); canvas.setFillColor(NAVY)
    canvas.drawString(ML, H-8.3*mm, 'WORKING FLOW COUPLING - BOUNDARY AND SOLVER SUMMARY')
    canvas.setFont(SANS, 6.7); canvas.setFillColor(MUTED)
    canvas.drawRightString(W-MR, H-8.3*mm, 'Current source snapshot: 12 Oct 2026')
    canvas.line(ML, 10*mm, W-MR, 10*mm)
    canvas.drawString(ML, 6.7*mm, 'working_flow_coupling_simple_long_tube.lua')
    canvas.drawRightString(W-MR, 6.7*mm, f'Page {doc.page} of 3')
    canvas.restoreState()


doc = BaseDocTemplate(str(OUT), pagesize=A4, leftMargin=ML, rightMargin=MR,
                      topMargin=MT, bottomMargin=MB,
                      title='Boundary and solver summary for working_flow_coupling_simple_long_tube.lua',
                      author='Generated for Baleeg Alhilali')
frame = Frame(ML, MB, CW, H-MT-MB, leftPadding=0, rightPadding=0,
              topPadding=0, bottomPadding=0)
doc.addPageTemplates(PageTemplate(id='all', frames=[frame], onPage=hf))

s=[]

# Page 1
s += heading('01 | Boundary map', 'Boundary conditions in the active model',
             'The model uses separate FV1 approximation spaces: P1/P1 velocity-pressure in the lumen and P1 Darcy pressure in the porous wall. Apical is the interface shared by their boundary descriptions.')
s += [BoundarySketch(), Spacer(1,1*mm)]
rows=[
 ['Boundary','Stokes side','Darcy side','Purpose'],
 ['Inlet','u = (0,0,vmax[1-(x^2+y^2)/R^2]); vmax=3 mm/s, R=1 mm.','No active Darcy boundary operator.','Prescribed parabolic axial inflow.'],
 ['Outlet','No-normal-stress/open outflow.','None.','Allows axial lumen discharge.'],
 ['Apical','Open-stress base plus WSBC: tangential BJS friction and normal Robin resistance.','Uniform Neumann flux q_D=-Q_S/A after Stokes is solved.','Reduced lumen/porous interface.'],
 ['OuterWall','None.','p_D=p_outer; default 0 Pa converted to kinematic pressure.','Porous pressure reference/sink.'],
 ['TripleJunction','u=v=w=0.','No explicit condition.','Closes velocity at the boundary meeting edge.'],
]
s += [table(rows,[25*mm,60*mm,47*mm,44*mm]), Spacer(1,3*mm)]
s += [p('The interface equations', 'S'),
      p('(sigma_S n)_t = - beta u_t,      beta = alpha nu / sqrt(K) = 200 mm/s', 'Eq'),
      p('(sigma_S n).n + R_D u_n = 0,      R_D = L_D/(K/nu) = 4 x 10^6 mm/s', 'Eq'),
      callout('Meaning', 'BJS friction resists tangential slip. The Robin term models the pressure resistance of a 4 mm porous layer in the normal direction. Both terms are assembled into the Stokes defect and Jacobian, so they are implicit in the Stokes matrix.', PGREEN, GREEN),
      PageBreak()]

# Page 2
s += heading('02 | Interface coupling', 'From Stokes leakage to Darcy flux',
             'The central idea is to replace simultaneous Stokes-Darcy coupling by a Stokes-side porous impedance, then reconstruct the Darcy field from the final lumen leakage.')
s += [CouplingSketch(), Spacer(1,1*mm), p('Implemented sequence', 'S')]
seq=[
 ['Step','Operation','Consequence'],
 ['1','Solve Stokes with BJS tangential friction and normal Robin resistance already in the Apical wall matrix.','Porous feedback is approximated before Darcy pressure is computed.'],
 ['2','Integrate Q_S = integral_Apical u_S . n_S dA.','Obtains the total signed leakage from the lumen.'],
 ['3','Use A = 2 pi R L = 40 pi = 125.6637 mm^2.','Converts total leakage to one mean normal velocity.'],
 ['4','Set q_D = -Q_S/A.','Minus sign accounts for opposite outward normals of the adjacent domains.'],
 ['5','Apply q_D uniformly as the Darcy Neumann condition; solve p_D and derive u_D=-(K/nu)grad(p_D).','Creates a porous pressure/velocity field consistent with the mean transferred flux.'],
]
s += [table(seq,[13*mm,94*mm,69*mm]), Spacer(1,3*mm)]
s += [callout('Correct name for this method', 'Partitioned Robin-reduced Stokes-Darcy coupling with mean-flux Darcy reconstruction. It is not a monolithic two-way solve: no computed Darcy pressure returns to the Stokes traction, and local Apical flux variation is replaced by its area average.', PCYAN, CYAN), Spacer(1,2*mm)]
s += [callout('Two accuracy cautions', '<b>Area:</b> the code divides by the analytic cylinder area, not the discrete mesh area, so integrated conservation depends on those areas matching. <b>Outer pressure:</b> a nonzero p_outer is applied to Darcy, but the corresponding forcing is not added to the homogeneous Stokes Robin equation.', PRED, RED), Spacer(1,2*mm)]
s += [p('Discretization at the interface', 'S'),
      p('FV1 uses vertex-centered control-volume balances. On Apical, WSBC adds traction terms to Stokes momentum rows, while the Darcy Neumann callback adds boundary flux. Equal-order Stokes P1/P1 uses the active <font name="Courier">fields</font> and <font name="Courier">cor</font> stabilization settings.', 'B'),
      PageBreak()]

# Page 3
s += heading('03 | Solver', 'Iterative hierarchy, stopping thresholds, and failed setups',
             'The same solver construction is instantiated separately for the Stokes matrix and the Darcy matrix. Each block is assembled, solved, and copied back to its grid function before the next stage.')
s += [SolverSketch(), Spacer(1,1*mm)]
rows=[
 ['Parameter','Active value','Interpretation'],
 ['Outer BiCGStab','1000 iterations; absolute 1e-14; reduction 1e-8','Main Krylov solve for each block.'],
 ['GMG cycle','V; 4 pre + 4 post; rap=false','Rediscretized coarse operators, not Galerkin R A P.'],
 ['ILU smoother','damping 0.7; overlap and consistent interfaces; inversionEps 1e-12','inversionEps is a pivot safeguard, not a residual tolerance.'],
 ['Coarse BiCGStab','500 iterations; absolute 1e-14; reduction 1e-10','Stricter relative coarse solve.'],
 ['Coarse ILUT','drop threshold 1e-8; Native Cuthill-McKee','Drop/fill control and ordering; 1e-8 is not convergence tolerance.'],
]
s += [table(rows,[37*mm,57*mm,82*mm]), Spacer(1,2*mm)]
s += [callout('Threshold interpretation', 'The absolute residual threshold is 10^-14, above double-precision machine epsilon (about 2.22 x 10^-16) but still strict. A converged algebraic residual does not prove mesh convergence or interface conservation; verify Q_S against the integrated Darcy boundary flux.', PGREEN, GREEN), Spacer(1,2*mm)]
fails=[
 ['Documented failed setup','Observed behavior','Fix now used'],
 ['Coarse solve left ambiguous instead of explicitly wrapped in AgglomeratingSolver.','Outer BiCGStab stagnated at numRefs=2 because level 0 was incomplete in parallel.','Explicit agglomerating wrapper around iterative coarse BiCGStab.'],
 ['Distribution delayed until the finest requested refinement.','Refined grid was duplicated as ghosts; numRefs=3 appeared stuck in assembly.','firstDistLvl=0: distribute base grid before refinement.'],
]
s += [p('Failures documented in this exact file', 'S'), table(fails,[55*mm,61*mm,60*mm]), Spacer(1,2*mm)]
s += [callout('Nonlinear caveat', 'The validated/default path is <font name="Courier">bstokes=true</font>. If it is set false, this file still performs only one linear assembly/solve and does not run a Newton loop, so it should not be treated as a complete nonlinear Navier-Stokes workflow without further implementation and testing.', PGOLD, GOLD)]

doc.build(s)
print(OUT)
