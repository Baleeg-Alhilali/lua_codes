#!/usr/bin/env python3
"""Generate a detailed, source-traceable explanation of the long-tube flow script."""

from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.platypus import (
    BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer, PageBreak,
    Table, TableStyle, Flowable, KeepTogether, Preformatted,
)


ROOT = Path('/Users/alhilaba/Documents/Promesh_swc/Model')
OUT = ROOT / 'output/pdf/working_flow_coupling_simple_long_tube_explained.pdf'
SOURCE = ROOT / 'code/flow/working_flow_coupling_simple_long_tube.lua'
SOURCE_HASH = 'f81a8d3cc01c3d4f829bac96cdae3546c88b05db8c92b03d13408165918bafd7'

PAGE_W, PAGE_H = A4
M_LEFT = 18 * mm
M_RIGHT = 18 * mm
M_TOP = 19 * mm
M_BOTTOM = 18 * mm
CONTENT_W = PAGE_W - M_LEFT - M_RIGHT

NAVY = colors.HexColor('#17324D')
BLUE = colors.HexColor('#2563A6')
CYAN = colors.HexColor('#18A6B8')
PALE_BLUE = colors.HexColor('#EAF3FA')
PALE_CYAN = colors.HexColor('#E9F8F8')
PALE_GOLD = colors.HexColor('#FFF5DC')
GOLD = colors.HexColor('#D69524')
RED = colors.HexColor('#B7463D')
PALE_RED = colors.HexColor('#FCECE9')
GREEN = colors.HexColor('#2C7A5A')
PALE_GREEN = colors.HexColor('#EAF6EF')
INK = colors.HexColor('#17212B')
MUTED = colors.HexColor('#596875')
LINE = colors.HexColor('#CCD7E0')
LIGHT = colors.HexColor('#F5F8FA')


def register_fonts():
    candidates = [
        ('Aptos', '/System/Library/Fonts/Supplemental/Arial.ttf'),
        ('AptosBold', '/System/Library/Fonts/Supplemental/Arial Bold.ttf'),
        ('Mono', '/System/Library/Fonts/Menlo.ttc'),
    ]
    for name, path in candidates:
        if Path(path).exists():
            try:
                pdfmetrics.registerFont(TTFont(name, path))
            except Exception:
                pass


register_fonts()
SANS = 'Aptos' if 'Aptos' in pdfmetrics.getRegisteredFontNames() else 'Helvetica'
SANS_B = 'AptosBold' if 'AptosBold' in pdfmetrics.getRegisteredFontNames() else 'Helvetica-Bold'
MONO = 'Courier'

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name='CoverTitle', fontName=SANS_B, fontSize=27, leading=31,
    textColor=colors.white, alignment=TA_LEFT, spaceAfter=8,
))
styles.add(ParagraphStyle(
    name='CoverSub', fontName=SANS, fontSize=12.5, leading=18,
    textColor=colors.HexColor('#DCEAF4'), alignment=TA_LEFT,
))
styles.add(ParagraphStyle(
    name='Kicker', fontName=SANS_B, fontSize=8.5, leading=10,
    textColor=CYAN, spaceAfter=4, tracking=1.2,
))
styles.add(ParagraphStyle(
    name='H1x', fontName=SANS_B, fontSize=22, leading=27,
    textColor=NAVY, spaceAfter=8,
))
styles.add(ParagraphStyle(
    name='H2x', fontName=SANS_B, fontSize=13.5, leading=17,
    textColor=BLUE, spaceBefore=7, spaceAfter=5,
))
styles.add(ParagraphStyle(
    name='H3x', fontName=SANS_B, fontSize=10.5, leading=13,
    textColor=NAVY, spaceBefore=5, spaceAfter=3,
))
styles.add(ParagraphStyle(
    name='Bodyx', fontName=SANS, fontSize=9.1, leading=13.1,
    textColor=INK, spaceAfter=5,
))
styles.add(ParagraphStyle(
    name='Smallx', fontName=SANS, fontSize=7.8, leading=10.6,
    textColor=MUTED, spaceAfter=3,
))
styles.add(ParagraphStyle(
    name='Bulletx', fontName=SANS, fontSize=8.8, leading=12.4,
    leftIndent=12, firstLineIndent=-7, bulletIndent=3, textColor=INK,
    spaceAfter=3,
))
styles.add(ParagraphStyle(
    name='TableHead', fontName=SANS_B, fontSize=7.5, leading=9.2,
    textColor=colors.white,
))
styles.add(ParagraphStyle(
    name='TableCell', fontName=SANS, fontSize=7.35, leading=9.4,
    textColor=INK,
))
styles.add(ParagraphStyle(
    name='TableCellSmall', fontName=SANS, fontSize=6.8, leading=8.3,
    textColor=INK,
))
styles.add(ParagraphStyle(
    name='Codex', fontName=MONO, fontSize=7.1, leading=9.2,
    textColor=INK, leftIndent=6, rightIndent=6, spaceBefore=3, spaceAfter=5,
    backColor=colors.HexColor('#F2F5F7'), borderColor=LINE,
    borderWidth=0.5, borderPadding=6,
))
styles.add(ParagraphStyle(
    name='CalloutTitle', fontName=SANS_B, fontSize=9.5, leading=12,
    textColor=NAVY, spaceAfter=3,
))
styles.add(ParagraphStyle(
    name='CalloutBody', fontName=SANS, fontSize=8.4, leading=11.5,
    textColor=INK,
))
styles.add(ParagraphStyle(
    name='Equation', fontName=MONO, fontSize=8.2, leading=11.2,
    textColor=NAVY, alignment=TA_CENTER, spaceBefore=5, spaceAfter=6,
))


def P(text, style='Bodyx'):
    return Paragraph(text, styles[style])


def bullet(text):
    return Paragraph('- ' + text, styles['Bulletx'])


def section(kicker, title, intro=None):
    parts = [P(kicker.upper(), 'Kicker'), P(title, 'H1x')]
    if intro:
        parts.append(P(intro, 'Bodyx'))
    parts.append(Spacer(1, 2 * mm))
    return parts


def callout(title, text, fill=PALE_BLUE, stripe=BLUE):
    inner = [P(title, 'CalloutTitle'), P(text, 'CalloutBody')]
    t = Table([[inner]], colWidths=[CONTENT_W - 1 * mm])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), fill),
        ('BOX', (0, 0), (-1, -1), 0.5, stripe),
        ('LINEBEFORE', (0, 0), (0, -1), 4, stripe),
        ('LEFTPADDING', (0, 0), (-1, -1), 11),
        ('RIGHTPADDING', (0, 0), (-1, -1), 9),
        ('TOPPADDING', (0, 0), (-1, -1), 8),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
    ]))
    return t


def data_table(rows, widths, header=True, small=False, aligns=None):
    cell_style = 'TableCellSmall' if small else 'TableCell'
    converted = []
    for r, row in enumerate(rows):
        converted.append([
            P(str(value), 'TableHead' if header and r == 0 else cell_style)
            for value in row
        ])
    t = Table(converted, colWidths=widths, repeatRows=1 if header else 0,
              hAlign='LEFT')
    cmds = [
        ('BACKGROUND', (0, 0), (-1, 0), NAVY if header else colors.white),
        ('GRID', (0, 0), (-1, -1), 0.35, LINE),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 5),
        ('RIGHTPADDING', (0, 0), (-1, -1), 5),
        ('TOPPADDING', (0, 0), (-1, -1), 5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
    ]
    for r in range(1 if header else 0, len(rows)):
        if r % 2 == (0 if header else 1):
            cmds.append(('BACKGROUND', (0, r), (-1, r), LIGHT))
    if aligns:
        for col, align in aligns.items():
            cmds.append(('ALIGN', (col, 0), (col, -1), align))
    t.setStyle(TableStyle(cmds))
    return t


class TubeDiagram(Flowable):
    def __init__(self, width=CONTENT_W, height=64 * mm):
        super().__init__()
        self.width, self.height = width, height

    def draw(self):
        c = self.canv
        x0, y0 = 12 * mm, 14 * mm
        w, h = self.width - 24 * mm, 28 * mm
        wall = 10 * mm
        # Porous annulus and lumen.
        c.setFillColor(PALE_GOLD)
        c.setStrokeColor(GOLD)
        c.roundRect(x0, y0, w, h, 5, fill=1, stroke=1)
        c.setFillColor(PALE_BLUE)
        c.setStrokeColor(BLUE)
        c.roundRect(x0, y0 + wall / 2, w, h - wall, 4, fill=1, stroke=1)
        c.setFont(SANS_B, 9)
        c.setFillColor(NAVY)
        c.drawCentredString(x0 + w / 2, y0 + h / 2 - 3, 'Lumen: Stokes / Navier-Stokes')
        c.setFillColor(colors.HexColor('#815B19'))
        c.setFont(SANS_B, 7.6)
        c.drawCentredString(x0 + w / 2, y0 + h - 8, 'Inter: porous Darcy layer')
        c.drawCentredString(x0 + w / 2, y0 + 2, 'Inter: porous Darcy layer')
        # Boundary labels.
        c.setFillColor(MUTED)
        c.setFont(SANS, 7)
        c.drawString(x0 - 1, y0 + h + 9, 'Inlet')
        c.drawRightString(x0 + w + 1, y0 + h + 9, 'Outlet')
        c.setStrokeColor(CYAN)
        c.setLineWidth(2)
        c.line(x0, y0 + wall / 2, x0 + w, y0 + wall / 2)
        c.line(x0, y0 + h - wall / 2, x0 + w, y0 + h - wall / 2)
        c.setFillColor(CYAN)
        c.setFont(SANS_B, 7.4)
        c.drawString(x0 + 2, y0 + h + 20, 'Apical interface Gamma')
        c.setStrokeColor(MUTED)
        c.setLineWidth(0.8)
        c.line(x0 + 37 * mm, y0 + h + 17, x0 + 37 * mm, y0 + h - wall / 2 + 1)
        # Flow arrows.
        c.setStrokeColor(BLUE)
        c.setFillColor(BLUE)
        for frac in (0.27, 0.47, 0.67):
            yy = y0 + h * frac
            xa, xb = x0 + 18 * mm, x0 + w - 18 * mm
            c.line(xa, yy, xb, yy)
            c.line(xb, yy, xb - 4, yy + 2)
            c.line(xb, yy, xb - 4, yy - 2)
        # Radial transfer arrows.
        c.setStrokeColor(RED)
        c.setFillColor(RED)
        for frac in (0.36, 0.62):
            xx = x0 + w * frac
            ya = y0 + h - wall / 2 - 1
            yb = y0 + h - 2
            c.line(xx, ya, xx, yb)
            c.line(xx, yb, xx - 2, yb - 4)
            c.line(xx, yb, xx + 2, yb - 4)
        # Legend.
        ly = 3 * mm
        c.setFillColor(BLUE); c.rect(x0, ly, 4, 4, fill=1, stroke=0)
        c.setFillColor(INK); c.setFont(SANS, 7)
        c.drawString(x0 + 7, ly, 'axial lumen flow')
        c.setFillColor(RED); c.rect(x0 + 43 * mm, ly, 4, 4, fill=1, stroke=0)
        c.setFillColor(INK); c.drawString(x0 + 43 * mm + 7, ly, 'normal interface flux')
        c.setFillColor(CYAN); c.rect(x0 + 98 * mm, ly, 4, 4, fill=1, stroke=0)
        c.setFillColor(INK); c.drawString(x0 + 98 * mm + 7, ly, 'BJS + Robin law')


class CouplingDiagram(Flowable):
    def __init__(self, width=CONTENT_W, height=76 * mm):
        super().__init__()
        self.width, self.height = width, height

    def box(self, c, x, y, w, h, title, lines, fill, stroke):
        c.setFillColor(fill); c.setStrokeColor(stroke); c.setLineWidth(1)
        c.roundRect(x, y, w, h, 5, fill=1, stroke=1)
        c.setFillColor(stroke); c.setFont(SANS_B, 8.5)
        c.drawString(x + 7, y + h - 13, title)
        c.setFillColor(INK); c.setFont(SANS, 6.9)
        yy = y + h - 25
        for line in lines:
            c.drawString(x + 7, yy, line)
            yy -= 9

    def arrow(self, c, x1, y1, x2, y2, label):
        c.setStrokeColor(NAVY); c.setFillColor(NAVY); c.setLineWidth(1.1)
        c.line(x1, y1, x2, y2)
        c.line(x2, y2, x2 - 5, y2 + 3)
        c.line(x2, y2, x2 - 5, y2 - 3)
        c.setFont(SANS_B, 6.8); c.drawCentredString((x1+x2)/2, y1 + 6, label)

    def draw(self):
        c = self.canv
        y = 23 * mm
        bw = 49 * mm
        bh = 35 * mm
        gap = 10 * mm
        x1 = 1 * mm
        x2 = x1 + bw + gap
        x3 = x2 + bw + gap
        self.box(c, x1, y, bw, bh, '1. STOKES BLOCK', [
            'FV1 on Lumen', 'BJS tangential friction',
            'implicit normal Robin resistance', 'solve with BiCGStab / GMG'], PALE_BLUE, BLUE)
        self.box(c, x2, y, bw, bh, '2. FLUX REDUCTION', [
            'integrate Q_S on Apical', 'divide by analytic cylinder area',
            'flip sign for opposite normals', 'produce one uniform flux value'], PALE_CYAN, CYAN)
        self.box(c, x3, y, bw, bh, '3. DARCY BLOCK', [
            'FV1 diffusion on Inter', 'Neumann flux on Apical',
            'p_D = p_outer on OuterWall', 'reconstruct Darcy velocity'], PALE_GOLD, GOLD)
        self.arrow(c, x1+bw, y+bh/2, x2, y+bh/2, 'Q_S')
        self.arrow(c, x2+bw, y+bh/2, x3, y+bh/2, '-Q_S / A')
        c.setFillColor(PALE_RED); c.setStrokeColor(RED)
        c.roundRect(33 * mm, 2 * mm, self.width - 66 * mm, 13 * mm, 4, fill=1, stroke=1)
        c.setFillColor(RED); c.setFont(SANS_B, 7.4)
        c.drawCentredString(self.width/2, 10 * mm,
            'No Darcy-pressure iterate returns to the Stokes solve: this is reduced one-way coupling.')
        c.setFillColor(INK); c.setFont(SANS, 6.8)
        c.drawCentredString(self.width/2, 6 * mm,
            'Porous feedback is represented by the Stokes-side Robin resistance, not by an outer fixed-point loop.')


class SolverDiagram(Flowable):
    def __init__(self, width=CONTENT_W, height=84 * mm):
        super().__init__()
        self.width, self.height = width, height

    def draw_box(self, c, x, y, w, h, fill, stroke, title, detail):
        c.setFillColor(fill); c.setStrokeColor(stroke)
        c.roundRect(x, y, w, h, 4, fill=1, stroke=1)
        c.setFillColor(stroke); c.setFont(SANS_B, 8)
        c.drawString(x+7, y+h-12, title)
        c.setFillColor(INK); c.setFont(SANS, 6.8)
        for i, line in enumerate(detail): c.drawString(x+7, y+h-22-i*8, line)

    def down(self, c, x, y1, y2):
        c.setStrokeColor(NAVY); c.line(x, y1, x, y2)
        c.line(x, y2, x-3, y2+5); c.line(x, y2, x+3, y2+5)

    def draw(self):
        c = self.canv
        cx = self.width/2
        w1, w2 = 104 * mm, 134 * mm
        self.draw_box(c, cx-w1/2, 63*mm, w1, 17*mm, PALE_BLUE, BLUE,
                      'OUTER KRYLOV: BiCGStab', ['max 1000; abs 1e-14; reduction 1e-8'])
        self.down(c, cx, 63*mm, 57*mm)
        self.draw_box(c, cx-w2/2, 37*mm, w2, 20*mm, PALE_CYAN, CYAN,
                      'PRECONDITIONER: geometric multigrid, V-cycle',
                      ['4 ILU pre-smooth + 4 ILU post-smooth; damping 0.7',
                       'rap=false; standard transfer; base level 0'])
        self.down(c, cx, 37*mm, 31*mm)
        self.draw_box(c, cx-w1/2, 8*mm, w1, 23*mm, PALE_GOLD, GOLD,
                      'COARSE LEVEL: AgglomeratingSolver',
                      ['BiCGStab: max 500; abs 1e-14; reduction 1e-10',
                       'ILUT drop 1e-8; Native Cuthill-McKee ordering'])
        c.setFillColor(MUTED); c.setFont(SANS, 6.8)
        c.drawCentredString(cx, 2*mm,
            'The same solver construction is instantiated separately for Stokes and Darcy.')


def header_footer(canvas, doc):
    canvas.saveState()
    if doc.page > 1:
        canvas.setStrokeColor(LINE)
        canvas.setLineWidth(0.45)
        canvas.line(M_LEFT, PAGE_H - 12 * mm, PAGE_W - M_RIGHT, PAGE_H - 12 * mm)
        canvas.setFont(SANS_B, 7.2)
        canvas.setFillColor(NAVY)
        canvas.drawString(M_LEFT, PAGE_H - 9 * mm, 'WORKING FLOW COUPLING - TECHNICAL EXPLANATION')
        canvas.setFont(SANS, 7)
        canvas.setFillColor(MUTED)
        canvas.drawRightString(PAGE_W - M_RIGHT, PAGE_H - 9 * mm,
                               'Source snapshot: 12 Oct 2026')
        canvas.line(M_LEFT, 11 * mm, PAGE_W - M_RIGHT, 11 * mm)
        canvas.setFont(SANS, 7)
        canvas.drawString(M_LEFT, 7 * mm, 'Baleeg Alhilali | UG4 long-tube model')
        canvas.drawRightString(PAGE_W - M_RIGHT, 7 * mm, f'Page {doc.page}')
    canvas.restoreState()


class ReportDoc(BaseDocTemplate):
    pass


doc = ReportDoc(
    str(OUT), pagesize=A4,
    leftMargin=M_LEFT, rightMargin=M_RIGHT,
    topMargin=M_TOP, bottomMargin=M_BOTTOM,
    title='Detailed implementation of working_flow_coupling_simple_long_tube.lua',
    author='Generated for Baleeg Alhilali',
    subject='UG4 FV1 Stokes-Darcy-BJS coupling and iterative solver explanation',
)
frame = Frame(M_LEFT, M_BOTTOM, CONTENT_W, PAGE_H - M_TOP - M_BOTTOM,
              leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
doc.addPageTemplates(PageTemplate(id='main', frames=[frame], onPage=header_footer))

story = []

# Cover
cover = Table([[
    [Spacer(1, 18 * mm),
     P('SOURCE-TRACED TECHNICAL REPORT', 'Kicker'),
     P('How the long-tube flow coupling is implemented', 'CoverTitle'),
     P('FV1 discretization, Stokes-Darcy interface treatment, BJS slip, solver hierarchy, thresholds, and failed configurations', 'CoverSub'),
     Spacer(1, 11 * mm),
     P('<b>File</b><br/>code/flow/working_flow_coupling_simple_long_tube.lua', 'CoverSub'),
     Spacer(1, 5 * mm),
     P('<b>Model</b><br/>A partitioned Robin-reduced Stokes-Darcy reconstruction in UG4', 'CoverSub'),
     Spacer(1, 5 * mm),
     P('<b>Snapshot</b><br/>Current working copy, 12 October 2026', 'CoverSub'),
     Spacer(1, 30 * mm),
     P('Prepared for Baleeg Alhilali', 'CoverSub')]
]], colWidths=[CONTENT_W], rowHeights=[PAGE_H - M_TOP - M_BOTTOM])
cover.setStyle(TableStyle([
    ('BACKGROUND', (0, 0), (-1, -1), NAVY),
    ('LEFTPADDING', (0, 0), (-1, -1), 16 * mm),
    ('RIGHTPADDING', (0, 0), (-1, -1), 16 * mm),
    ('TOPPADDING', (0, 0), (-1, -1), 0),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 0),
]))
story += [cover, PageBreak()]

# Page 2 - reading guide / executive summary
story += section('01 | Executive reading', 'What the program actually solves',
    'The most important point is architectural: the file does not assemble one monolithic Stokes-Darcy matrix. It solves a lumen problem, compresses its interface flux to one mean value, and then solves a porous-wall problem. The Stokes-side Robin term stands in for the pressure resistance of a finite porous layer.')
story += [callout('One-sentence interpretation',
    'A pressure-driven lumen flow loses fluid radially through an Apical interface whose tangential motion is damped by a Beavers-Joseph-Saffman law and whose normal motion is resisted by an equivalent Darcy-layer impedance; the resulting total radial flux is then imposed on a separate Darcy pressure solve.', PALE_CYAN, CYAN), Spacer(1, 5*mm)]
story += [P('Implemented pipeline', 'H2x')]
for txt in [
    '<b>Geometry:</b> load <font name="Courier">runs/long_tube.ugx</font>, verify required subsets and a cylindrical projector, repair axis labels, distribute, then refine.',
    '<b>Lumen:</b> solve an equal-order P1/P1 FV1 Stokes problem on <font name="Courier">Lumen</font> with prescribed parabolic inflow, open outlet, and interface stress terms on <font name="Courier">Apical</font>.',
    '<b>Interface reduction:</b> integrate the Stokes normal velocity over Apical, divide by the analytic cylinder area, and reverse the sign to respect the opposite porous-domain normal.',
    '<b>Porous wall:</b> solve a scalar P1 FV1 Darcy-pressure diffusion problem on <font name="Courier">Inter</font>, using the transferred uniform normal flux and an outer pressure.',
    '<b>Output:</b> write lumen velocity/pressure and porous pressure/Darcy velocity to separate VTK series.',
]: story.append(bullet(txt))
story += [Spacer(1, 4*mm), P('What this report means by "coupling"', 'H2x'),
          P('The word coupling is used in a reduced sense. The normal Robin coefficient in the Stokes wall condition encodes a Darcy-layer pressure drop model, and the final Stokes flux drives the Darcy reconstruction. There is no iteration in which the computed Darcy pressure is fed back into Stokes. Therefore, this is best described as <b>partitioned, Robin-reduced, and one-way after the Stokes solve</b>.', 'Bodyx'),
          callout('Scope and confidence',
            'Statements about implementation are traced to the current Lua source. Two failed configurations are explicitly documented by comments in this file. A third observation, about Gauss-Seidel plus Galerkin RAP, comes from a related flow script and is labeled separately later. No numerical run was required to prepare this code-reading report.', PALE_GOLD, GOLD),
          PageBreak()]

# Page 3 - source map
story += section('02 | Program anatomy', 'Where each responsibility lives in the Lua file')
rows = [
    ['Source lines', 'Responsibility', 'Why it matters'],
    ['1-20', 'Resolve Model root, change working directory, load UG4 utilities.', 'Makes grid and output paths independent of the shell launch directory.'],
    ['22-59', 'Read run parameters; define geometry, material values, result names, and load balancing.', 'Defines the physical scale and the discretization hierarchy.'],
    ['61-133', 'Load and validate the domain; repair axis labels; distribute and refine.', 'Produces a usable distributed mesh and avoids projector/refinement failures.'],
    ['135-159', 'Build separate Stokes and Darcy approximation spaces.', 'Prevents one global mixed cross-domain matrix.'],
    ['161-217', 'Convert permeability/pressure units; define inlet and auxiliary callbacks.', 'Connects SI inputs to the millimeter-based UG4 model.'],
    ['219-275', 'Assemble Stokes FV1, inlet/outlet terms, and WSBC interface terms.', 'Implements lumen equations and BJS/Robin interface physics.'],
    ['277-300', 'Assemble Darcy FV1, interface Neumann condition, and outer Dirichlet pressure.', 'Implements the porous pressure reconstruction.'],
    ['302-371', 'Define iterative BiCGStab/GMG/ILU/ILUT solver hierarchy.', 'Controls convergence and parallel coarse-level robustness.'],
    ['373-396', 'Assemble and solve Stokes; integrate interface flux; solve Darcy.', 'This is the actual sequential coupling algorithm.'],
    ['398-436', 'Construct Darcy velocity and write VTK output.', 'Converts pressure gradient to velocity and exports both subproblems.'],
]
story += [data_table(rows, [22*mm, 63*mm, 88*mm]), Spacer(1, 5*mm)]
story += [P('A compact execution trace', 'H2x'),
          Preformatted("""parse parameters -> load/repair/refine mesh
        -> create two approximation spaces
        -> assemble Stokes operators and interface WSBC
        -> solve Stokes
        -> integrate Q_S over Apical
        -> set q_D = -Q_S / A_Gamma
        -> assemble and solve Darcy pressure
        -> derive Darcy velocity from grad(p_D)
        -> export two VTK series""", styles['Codex']),
          callout('A subtle but important default',
            'The command-line option <font name="Courier">-bstokes</font> defaults to true. With that default, the momentum equation is linear Stokes. Although the file accepts false, its solve path still performs one linear assembly/solve and does not run a Newton loop. The false setting should not be treated as a validated nonlinear Navier-Stokes workflow without additional verification.', PALE_RED, RED),
          PageBreak()]

# Page 4 geometry
story += section('03 | Geometry and units', 'The two physical regions and their named boundaries')
story += [TubeDiagram(), Spacer(1, 2*mm)]
story += [P('Subset roles', 'H2x')]
rows = [
    ['Subset', 'Used by', 'Mathematical role'],
    ['Lumen', 'Stokes volume', 'Viscous incompressible flow region.'],
    ['Inter', 'Darcy volume', 'Porous-wall pressure region.'],
    ['Apical', 'Both sides', 'Lumen/porous interface Gamma; receives open stress plus BJS/Robin on Stokes and flux on Darcy.'],
    ['Inlet', 'Stokes boundary; space support', 'Parabolic axial velocity.'],
    ['Outlet', 'Stokes boundary; space support', 'No-normal-stress outflow.'],
    ['OuterWall', 'Darcy boundary', 'Prescribed kinematic outer pressure.'],
    ['TripleJunction', 'Stokes boundary; both supports', 'No-slip velocity; closes edge where boundary subsets meet.'],
]
story += [data_table(rows, [27*mm, 47*mm, 99*mm], small=True), Spacer(1, 4*mm)]
story += [P('Unit conversion performed by the code', 'H2x')]
for txt in [
    '<b>Permeability:</b> K = 10^-12 m^2 x 10^6 = 10^-6 mm^2.',
    '<b>Darcy mobility:</b> M = K / nu = 10^-6 s, because nu = 1 mm^2/s.',
    '<b>BJS friction:</b> beta = alpha nu / sqrt(K) = 0.2 x 1 / 0.001 = 200 mm/s.',
    '<b>Kinematic pressure:</b> p_kin = p_Pa / (rho x 10^-6), in mm^2/s^2. The default outer pressure is zero.',
    '<b>Normal resistance:</b> R_D = L_D / M = 4 / 10^-6 = 4 x 10^6 mm/s.',
]: story.append(bullet(txt))
story += [PageBreak()]

# Page 5 PDE
story += section('04 | Continuous model', 'Equations represented before spatial discretization')
story += [P('Lumen momentum and incompressibility', 'H2x'),
          P('With the default <font name="Courier">bstokes=true</font>, the code represents a steady Stokes system in the lumen. Pressure is kinematic pressure, consistent with the use of kinematic viscosity.', 'Bodyx'),
          P('-div(sigma_S(u,p)) = 0,       div(u) = 0       in Omega_S', 'Equation'),
          P('sigma_S(u,p) = 2 nu D(u) - p I,       D(u) = (grad u + grad u^T)/2', 'Equation'),
          P('The call <font name="Courier">set_laplace(false)</font> selects the symmetric-gradient stress form. This matters at a slip interface because the tangential traction should be derived from the physical viscous stress rather than a componentwise Laplacian shortcut.', 'Bodyx'),
          P('Porous pressure and Darcy velocity', 'H2x'),
          P('The scalar convection-diffusion discretization has zero mass and diffusion coefficient M = K/nu. In steady form it represents:', 'Bodyx'),
          P('-div(M grad p_D) = 0       in Omega_D', 'Equation'),
          P('u_D = -M grad p_D', 'Equation'),
          P('The output linker is configured with permeability K, viscosity nu, density 1, no gravity, and the pressure gradient; this produces the Darcy velocity used in the VTK output.', 'Bodyx'),
          P('External boundary data', 'H2x')]
for txt in [
    '<b>Inlet:</b> u = (0, 0, v_max [1 - (x^2+y^2)/R^2]) inside the radius; it is exactly zero on or outside the rim.',
    '<b>Outlet:</b> no-normal-stress/open outflow on <font name="Courier">Outlet</font>.',
    '<b>Outer porous wall:</b> p_D = p_outer on <font name="Courier">OuterWall</font>.',
    '<b>Triple junction:</b> u = 0.',
]: story.append(bullet(txt))
story += [callout('Pressure convention',
    'The code uses kinematic pressure. If a physical pressure difference is supplied in pascals, division by density converts it to m^2/s^2 and the extra 10^6 converts to mm^2/s^2. Mixing dynamic pressure in pascals directly with this system would give the wrong scale.', PALE_GOLD, GOLD), PageBreak()]

# Page 6 discretization
story += section('05 | FV1 discretization', 'How the continuous equations become algebraic systems')
story += [P('Approximation spaces', 'H2x'),
          P('Both subproblems use first-order Lagrange basis functions, but they live in separate approximation spaces. The Stokes space contains four nodal unknowns <font name="Courier">u, v, w, pressure</font>; the Darcy space contains one nodal unknown <font name="Courier">darcyPressure</font>. The separation means the two algebraic matrices are built and solved independently.', 'Bodyx')]
rows = [
    ['Block', 'Unknowns', 'Basis', 'Element scheme', 'Support'],
    ['Stokes', 'u, v, w, p', 'P1 nodal', 'NavierStokesFV1', 'Lumen and named boundary subsets'],
    ['Darcy', 'p_D', 'P1 nodal', 'ConvectionDiffusionFV1', 'Inter and named boundary subsets'],
]
story += [data_table(rows, [28*mm, 31*mm, 26*mm, 43*mm, 45*mm]), Spacer(1, 5*mm)]
story += [P('Vertex-centered control-volume idea', 'H2x'),
          P('FV1 uses the P1 finite-element geometry and basis to form a dual control volume around each mesh vertex. Each row is a balance over that dual volume. Volume terms become sums of fluxes over dual faces, and boundary terms become flux or traction contributions over the part of the dual boundary lying on a named subset.', 'Bodyx'),
          P('For a Darcy control volume V_i:', 'Bodyx'),
          P('sum over dual faces f:  [-M (grad p_D)|_f . n_f] |f|  = boundary-source contribution', 'Equation'),
          P('For the Stokes momentum rows, viscous traction and pressure traction cross the dual faces; the continuity row balances the velocity flux. The interface WSBC contributes directly to the momentum residual and Jacobian on Apical.', 'Bodyx'),
          P('Equal-order Stokes stabilization', 'H2x'),
          P('Velocity and pressure are both P1, so the code enables the stabilization controls named <font name="Courier">fields</font> and <font name="Courier">cor</font>. These are part of the selected UG4 FV1 formulation and are important because a naive equal-order velocity-pressure discretization may admit unstable pressure modes.', 'Bodyx'),
          callout('What is implicit at the interface',
            'The BJS tangential coefficient and the normal Robin resistance are assembled in <font name="Courier">NavierStokesWSBCFV1</font>. Because that operator contributes to both defect and Jacobian, these two Stokes-side interface terms are not post-processing corrections; they are part of the matrix solved by BiCGStab.', PALE_GREEN, GREEN),
          PageBreak()]

# Page 7 interface theory
story += section('06 | Interface law', 'Tangential BJS slip plus a reduced normal Darcy resistance')
story += [P('Tangential direction: Beavers-Joseph-Saffman', 'H2x'),
          P('The wall operator receives <font name="Courier">sliding_factor = beta</font> and <font name="Courier">sliding_limit = 0</font>. The code comments interpret this as:', 'Bodyx'),
          P('(sigma_S n)_t = - beta u_t,       beta = alpha nu / sqrt(K)', 'Equation'),
          P('At the default values, beta = 200 mm/s. Tangential traction therefore opposes tangential lumen velocity. Larger alpha, larger viscosity, or smaller permeability increases tangential resistance and moves the interface behavior toward no-slip.', 'Bodyx'),
          P('Normal direction: finite porous-layer reduction', 'H2x'),
          P('Assume normal Darcy flow crosses a porous layer of effective thickness L_D and ends at a prescribed outer pressure. In a one-dimensional normal approximation:', 'Bodyx'),
          P('u_n = -M (p_outer - p_interface)/L_D', 'Equation'),
          P('Rearranging gives a pressure-drop relation proportional to normal speed:', 'Bodyx'),
          P('p_interface - p_outer = (L_D/M) u_n = R_D u_n', 'Equation'),
          P('The Stokes wall operator is configured with <font name="Courier">normal_factor = R_D</font>, giving the homogeneous default form stated in the source:', 'Bodyx'),
          P('(sigma_S n).n + R_D u_n = 0', 'Equation'),
          P('This is the key idea: instead of solving the Darcy pressure simultaneously with Stokes, the porous wall is condensed to a local hydraulic impedance. It lets normal leakage occur while penalizing it according to permeability and layer thickness.', 'Bodyx'),
          callout('Important nonzero-pressure limitation',
            'The default outer pressure is zero, which is consistent with the homogeneous Robin equation above. If <font name="Courier">outerPressurePa</font> is changed, that value is applied to the Darcy Dirichlet boundary, but no corresponding nonzero forcing is added to the Stokes Robin condition. A physically consistent nonzero outer pressure would require the appropriate traction/pressure term on the Stokes side.', PALE_RED, RED),
          PageBreak()]

# Page 8 coupling algorithm
story += section('07 | Coupling algorithm', 'How information crosses Apical in the implemented solve')
story += [CouplingDiagram(), Spacer(1, 2*mm), P('Step-by-step transfer', 'H2x')]
for txt in [
    '<b>1.</b> Solve the Stokes block with the Robin normal impedance and BJS tangential friction already inside its matrix.',
    '<b>2.</b> Evaluate <font name="Courier">IntegralNormalComponentOnManifold</font> over Apical, using Lumen as the adjacent volume subset. This returns the total signed normal Stokes flux Q_S.',
    '<b>3.</b> Compute the analytic lateral cylinder area A_Gamma = 2 pi R L. With R = 1 mm and L = 20 mm, A_Gamma = 40 pi = approximately 125.6637 mm^2.',
    '<b>4.</b> Set the scalar callback value q_D = -Q_S/A_Gamma. The minus sign accounts for the fact that outward normals of adjacent lumen and porous regions point in opposite directions.',
    '<b>5.</b> The Darcy Neumann boundary calls that scalar callback, so every Apical integration point receives the same mean normal flux.',
    '<b>6.</b> Solve Darcy pressure, then reconstruct u_D from its gradient.',
]: story.append(bullet(txt))
story += [Spacer(1, 3*mm), callout('Conservation level',
    'The transfer preserves the integrated flux only if the analytic cylinder area used for division matches the discrete Apical area used by the boundary quadrature. It does not preserve the local spatial flux distribution: all interface variations are compressed to one mean value. A mesh-based area integral would remove the analytic-area mismatch risk.', PALE_GOLD, GOLD),
          PageBreak()]

# Page 9 boundary assembly
story += section('08 | Boundary assembly', 'Which operator contributes what on each boundary')
rows = [
    ['Boundary', 'Stokes contribution', 'Darcy contribution', 'Interpretation'],
    ['Inlet', 'Dirichlet parabolic velocity.', 'None in the assembled Darcy equation.', 'Drives axial lumen flow.'],
    ['Outlet', 'No-normal-stress outflow.', 'None.', 'Lets lumen flow leave axially.'],
    ['Apical', 'No-normal-stress/open base plus WSBC tangential BJS and normal Robin terms.', 'Neumann normal flux q_D = -Q_S/A.', 'Reduced interface between lumen and porous wall.'],
    ['OuterWall', 'None.', 'Dirichlet p_D = p_outer.', 'Pressure sink/reference for porous flow.'],
    ['TripleJunction', 'Velocity Dirichlet zero.', 'No explicit boundary equation.', 'Closes velocity at the meeting edge.'],
]
story += [data_table(rows, [25*mm, 60*mm, 46*mm, 42*mm], small=True), Spacer(1, 5*mm)]
story += [P('How the two Apical Stokes operators coexist', 'H2x'),
          P('The source adds <font name="Courier">NavierStokesNoNormalStressOutflow</font> on <font name="Courier">Outlet,Apical</font> and then adds <font name="Courier">NavierStokesWSBCFV1</font> on Apical. The intended composition is an open/no-prescribed-traction base plus the WSBC correction that supplies tangential friction and normal resistance. Exact sign and overlap behavior remain governed by UG4 operator conventions, so a manufactured boundary test is the safest verification.', 'Bodyx'),
          P('Callbacks that are present but unused', 'H2x'),
          P('The file defines radial normal functions for the lumen and Darcy sides, but it does not pass them to the boundary discretizations. The actual interface normals are therefore obtained internally from the mesh geometry by the UG4 boundary operators. The functions document the intended orientation but do not control the current assembly.', 'Bodyx'),
          P('Algebraic form of the sequential method', 'H2x'),
          P('[ A_S + B_BJS + B_Robin    G ] [ u ] = [ f_in ]', 'Equation'),
          P('[          D               0 ] [ p ]   [  0   ]', 'Equation'),
          P('Q_S = integral_Gamma u . n_S dA,        q_D = -Q_S/A_Gamma', 'Equation'),
          P('A_D p_D = b_outer + b_Neumann(q_D)', 'Equation'),
          callout('Not a monolithic interface matrix',
            'There are no off-diagonal matrix blocks mapping p_D into the Stokes traction or u directly into Darcy during one global solve. The only cross-block data passed after assembly is the scalar mean flux.', PALE_RED, RED),
          PageBreak()]

# Page 10 solver
story += section('09 | Linear solver', 'BiCGStab preconditioned by geometric multigrid')
story += [SolverDiagram(), Spacer(1, 2*mm)]
story += [P('What happens inside <font name="Courier">solveLinearSubproblem</font>', 'H2x')]
for txt in [
    'Create an assembled linear operator A, a right-hand side b, and a zero initial algebraic solution x.',
    'Let the domain discretization adjust the initial solution for constraints and Dirichlet values.',
    'Call <font name="Courier">assemble_linear</font> to build A and b.',
    'Build a solver descriptor tied to the corresponding approximation space and discretization.',
    'Initialize the linear solver with A, then call <font name="Courier">apply(x,b)</font>. Failure raises an assertion.',
    'Copy the algebraic solution back into the grid function.',
]: story.append(bullet(txt))
story += [P('Why multigrid helps', 'H2x'),
          P('The outer BiCGStab iteration handles the nonsymmetric or indefinite algebraic behavior. GMG supplies a hierarchy-aware preconditioner: ILU smoothing reduces high-frequency error on each mesh level, restriction transfers residual to coarser levels, the coarse solver addresses global low-frequency error, and prolongation corrects the fine level.', 'Bodyx'),
          P('The setting <font name="Courier">rap=false</font> requests rediscretized level operators instead of forming coarse operators by the Galerkin product R A P. Four pre-smoothing and four post-smoothing steps are used. Transfer damping is 1.0 in both directions, and P1 optimization is enabled.', 'Bodyx'),
          PageBreak()]

# Page 11 threshold glossary
story += section('10 | Thresholds', 'Every numerical tolerance in the active solver block')
rows = [
    ['Location / parameter', 'Value', 'Meaning'],
    ['Outer BiCGStab: iterations', '1000', 'Maximum Krylov iterations for each Stokes or Darcy solve.'],
    ['Outer BiCGStab: absolute', '1e-14', 'Absolute residual stopping scale requested by the convergence check.'],
    ['Outer BiCGStab: reduction', '1e-8', 'Required residual reduction relative to the starting residual.'],
    ['Coarse BiCGStab: iterations', '500', 'Maximum iterations inside the agglomerated level-0 solve.'],
    ['Coarse BiCGStab: absolute', '1e-14', 'Absolute coarse-solve residual threshold.'],
    ['Coarse BiCGStab: reduction', '1e-10', 'Stricter relative reduction for the coarse problem.'],
    ['ILU smoother: inversionEps', '1e-12', 'Small-pivot/inversion safeguard; not a convergence tolerance.'],
    ['ILUT coarse preconditioner: threshold', '1e-8', 'Drop threshold controlling fill/approximation; not a residual tolerance.'],
    ['ILU smoother: damping', '0.7', 'Scales the smoothing correction to improve robustness.'],
    ['GMG: pre / post smooth', '4 / 4', 'Number of ILU smoothing applications around coarse correction.'],
]
story += [data_table(rows, [67*mm, 25*mm, 81*mm], small=True), Spacer(1, 5*mm)]
story += [callout('How to read 1e-14',
    'The requested absolute threshold is 10^-14, safely above double-precision machine epsilon (about 2.22 x 10^-16). It is still very strict. Whether the absolute criterion is reached depends on equation scaling; the relative reduction criterion may terminate first. Always interpret the printed residual history together with physical balance errors.', PALE_GREEN, GREEN), Spacer(1, 4*mm)]
story += [P('Residual thresholds are not physical error bounds', 'H2x'),
          P('A small algebraic residual says that the discrete equations have been solved accurately relative to their scaling. It does not, by itself, prove mesh convergence, local interface conservation, or correctness of the reduced coupling model. Those require refinement studies and conservation checks.', 'Bodyx'),
          P('Verbose output', 'H2x'),
          P('Both convergence checks have <font name="Courier">verbose=true</font>, so their iteration histories should appear in the terminal. The Stokes and Darcy blocks are solved separately, which makes their convergence histories distinguishable in sequence.', 'Bodyx'),
          PageBreak()]

# Page 12 failure history
story += section('11 | What was tried and did not work', 'Evidence preserved in comments and related scripts')
story += [callout('Evidence rule used here',
    'Only configurations documented by source comments are called failures. The report does not infer that every alternate file or solver experiment failed.', PALE_BLUE, BLUE), Spacer(1, 4*mm)]
rows = [
    ['Configuration', 'Observed problem', 'Current remedy', 'Evidence'],
    ['Implicit/ambiguous gathered coarse solve without an explicit agglomerating wrapper.',
     'Outer BiCGStab worked on an incomplete distributed level-0 problem and stagnated at numRefs=2.',
     'Wrap the iterative coarse BiCGStab explicitly in AgglomeratingSolver; keep gatheredBaseSolverIfAmbiguous=true as an additional safeguard.',
     'Comments at source lines 329-331.'],
    ['Delay distribution until the finest requested refinement.',
     'Refined grids were duplicated as ghosts and numRefs=3 appeared stuck during assembly.',
     'Set firstDistLvl=0 so the base grid is distributed before refinement.',
     'Comments at source lines 50-52.'],
    ['Scalar Gauss-Seidel smoother with Galerkin RAP in a related projector-flow experiment.',
     'Documented divergence at refinement 2 in that related script.',
     'That experiment switched to damped ILU and rediscretized coarse operators. The present file likewise uses ILU and rap=false, but uses an iterative agglomerated coarse solver rather than SuperLU.',
     'Related file: code/coupled/Dif_itterative_flow_projector.lua, around line 240.'],
]
story += [data_table(rows, [45*mm, 46*mm, 58*mm, 24*mm], small=True), Spacer(1, 5*mm)]
story += [P('Why the explicit agglomeration is important', 'H2x'),
          P('At level 0 in parallel, no single process necessarily owns a complete coarse problem. A local coarse solve can therefore fail to represent the global low-frequency error. The agglomerating wrapper gathers or consolidates the coarse system so the nested iterative solver operates on a complete coarse problem.', 'Bodyx'),
          callout('What is not established by the source',
            'The repository contains monolithic and alternative coupling scripts, but the requested file does not document that those methods failed. It is therefore inaccurate to list them as failed attempts without run logs or explicit comments.', PALE_GOLD, GOLD),
          PageBreak()]

# Page 13 mesh processing
story += section('12 | Mesh preparation', 'Why the code repairs labels and distributes early')
story += [P('Required mesh metadata', 'H2x'),
          P('The program refuses to continue unless the UGX grid contains a <font name="Courier">ProjectionHandler</font> and a cylindrical projector. Curved boundary refinement relies on this metadata; silently refining without it could destroy the intended tube geometry.', 'Bodyx'),
          P('Axis-label repair', 'H2x'),
          P('Vertices near x = 0 and y = 0 can be mislabeled as <font name="Courier">Inter</font> even though they lie on the open tube axis. The function <font name="Courier">assign_axis_vertices_to_lumen</font> selects a narrow vertical cuboid around the axis and reassigns selected vertices to <font name="Courier">Lumen</font>. The repair runs once before and once after refinement.', 'Bodyx'),
          P('The half-width is max(10^-8, 10^-6 R). With R = 1 mm, this is 10^-6 mm. It is deliberately tiny so only axis vertices are targeted.', 'Bodyx'),
          P('Distribution order', 'H2x'),
          P('The active load-balancing configuration uses <font name="Courier">firstDistLvl=0</font>. The base mesh is distributed before the serial refinement loop. This avoids the documented ghost duplication/assembly stall that occurred when distribution was deferred to the final refinement level.', 'Bodyx'),
          P('Refinement structure', 'H2x')]
for txt in [
    'Create the domain with zero initial refinements in <font name="Courier">util.CreateDomain</font>.',
    'Repair axis vertices and distribute the base grid.',
    'Perform <font name="Courier">numRefs</font> global refinement steps through the global refiner.',
    'Repair axis vertices again because refinement may create new vertices on the axis.',
    'Create algebraic IDs, print layout statistics, and then initialize both approximation spaces.',
]: story.append(bullet(txt))
story += [callout('Mesh-dependent coupling check',
    'Because the transfer divides by an analytic area, compute the discrete Apical area independently when validating a new grid. A discrepancy changes the total Darcy flux even when Q_S is integrated correctly.', PALE_GREEN, GREEN),
          PageBreak()]

# Page 14 outputs / validation
story += section('13 | Outputs and verification', 'What is written and what should be checked')
rows = [
    ['VTK series', 'Fields', 'Where'],
    ['Stokes output', 'u, v, w as velocity; pressure as stokes_pressure', 'Nodal data on Lumen'],
    ['Darcy output', 'darcyPressure; DarcyVelocityLinker components as darcy_velocity', 'Pressure nodal on Inter; velocity element data on Inter'],
]
story += [data_table(rows, [38*mm, 73*mm, 62*mm]), Spacer(1, 5*mm)]
story += [P('Recommended verification ladder', 'H2x')]
checks = [
    ('1. Algebraic convergence', 'Confirm both verbose BiCGStab histories meet the declared stopping criteria before their iteration limits.'),
    ('2. Interface flux sign', 'Check that outward lumen leakage produces the expected sign of Q_S and that the transferred Darcy flux points into the porous region.'),
    ('3. Integrated conservation', 'Integrate the resulting Darcy normal flux over discrete Apical and compare it with -Q_S. Report relative mismatch.'),
    ('4. Pressure-drop consistency', 'Compare the mean interface pressure drop with R_D times mean normal speed for the reduced one-dimensional model.'),
    ('5. Mesh refinement', 'Repeat at increasing numRefs and track lumen flow rate, Q_S, porous pressure, and conservation mismatch.'),
    ('6. Limiting behavior', 'Decrease K: leakage should fall and BJS friction should rise. Increase L_D: normal leakage should fall. Set alpha larger: tangential slip should decrease.'),
    ('7. Boundary composition', 'Use a manufactured or simplified test to confirm the sign and combined effect of open-stress plus WSBC terms on Apical.'),
]
for title, detail in checks:
    story += [KeepTogether([P(title, 'H3x'), P(detail, 'Bodyx')])]
story += [callout('Best single improvement for coupling fidelity',
    'Transfer the spatial Stokes normal-flux field to the Darcy interface instead of reducing it to one scalar mean, then iterate Darcy pressure feedback into the Stokes normal traction. That would move the method from a reduced one-way reconstruction toward a genuine two-way partitioned Stokes-Darcy coupling.', PALE_CYAN, CYAN),
          PageBreak()]

# Page 15 limitations
story += section('14 | Assumptions and limitations', 'What the current formulation intentionally leaves out')
limitations = [
    ('One-way pressure feedback', 'The Darcy solution is computed after Stokes and never changes the Stokes solution. The Robin impedance is a surrogate for porous feedback.'),
    ('Uniform transferred flux', 'A single scalar q_D is imposed over the entire interface; local leakage structure is discarded.'),
    ('Analytic interface area', 'The area is hard-coded as 2 pi R L rather than integrated from the actual mesh.'),
    ('Homogeneous Stokes Robin forcing', 'A nonzero outer pressure is not included in the Stokes-side normal boundary forcing.'),
    ('Linear solve path', 'The active/default Stokes case is linear. Setting bstokes=false does not add the Newton loop normally needed for a consistent nonlinear solve.'),
    ('No gravity in Darcy velocity', 'The linker explicitly disables gravity; only pressure-gradient-driven Darcy velocity is reconstructed.'),
    ('Parameter constants', 'Material and geometry values are embedded in the file, except for numRefs, bstokes, and vmax.'),
    ('Operator convention dependency', 'The precise traction signs and the combination of open-boundary and WSBC terms depend on UG4 implementation conventions.'),
]
rows = [['Limitation', 'Consequence']] + [[a,b] for a,b in limitations]
story += [data_table(rows, [55*mm, 118*mm], small=True), Spacer(1, 5*mm)]
story += [P('Why the reduced idea is still useful', 'H2x'),
          P('The model is inexpensive and structurally simple: it avoids a global mixed saddle-point system across two domains, retains BJS tangential physics, and represents normal porous resistance implicitly in the Stokes matrix. It is a practical screening model when the porous layer is thin, its effective resistance is known, and detailed local pressure feedback is not the primary quantity of interest.', 'Bodyx'),
          callout('Use the right name',
            'For publications or presentations, describe this implementation as a "partitioned Robin-reduced Stokes-Darcy model with mean-flux Darcy reconstruction." Calling it a fully coupled monolithic Stokes-Darcy discretization would overstate what the code assembles.', PALE_GREEN, GREEN),
          PageBreak()]

# Page 16 parameter audit and traceability
story += section('15 | Parameter audit', 'Active defaults and source traceability')
rows = [
    ['Quantity', 'Symbol / setting', 'Default'],
    ['Tube inlet radius', 'R', '1 mm'],
    ['Tube length', 'L', '20 mm'],
    ['Kinematic viscosity', 'nu', '1 mm^2/s'],
    ['Permeability', 'K', '10^-12 m^2 = 10^-6 mm^2'],
    ['BJS coefficient', 'alpha', '0.2'],
    ['Density', 'rho', '1000 kg/m^3'],
    ['Outer pressure', 'p_outer', '0 Pa'],
    ['Equivalent porous resistance length', 'L_D', '4 mm'],
    ['Maximum inlet speed', 'v_max', '3 mm/s'],
    ['Global refinements', 'numRefs', '2'],
    ['Stokes mode', 'bstokes', 'true'],
    ['Distribution level', 'firstDistLvl', '0'],
]
story += [data_table(rows, [62*mm, 48*mm, 63*mm], small=True), Spacer(1, 5*mm)]
story += [P('Source identity', 'H2x'),
          P('<b>Absolute path:</b><br/><font name="Courier">/Users/alhilaba/Documents/Promesh_swc/Model/code/flow/working_flow_coupling_simple_long_tube.lua</font>', 'Bodyx'),
          P('<b>SHA-256:</b><br/><font name="Courier">' + SOURCE_HASH + '</font>', 'Bodyx'),
          P('<b>Snapshot date:</b> 12 October 2026. This report describes the current working-copy contents at that hash; later edits may make line references or behavior differ.', 'Bodyx'),
          P('Interpretation notes', 'H2x')]
for txt in [
    'Mathematical equations in this report are reconstructed from the operators and comments in the Lua source; they are not copied from external documentation.',
    'The exact discrete stencil is implemented inside UG4 classes such as NavierStokesFV1, NavierStokesWSBCFV1, and ConvectionDiffusionFV1, not expanded in the Lua file.',
    'The failed-solver statements are limited to comments preserved in this file plus one explicitly labeled related-projector script observation.',
]: story.append(bullet(txt))
story += [Spacer(1, 4*mm), callout('Bottom line',
    'The implementation first solves a stabilized P1/P1 FV1 lumen Stokes problem with BJS tangential friction and a Darcy-derived normal Robin resistance. It then conserves the integrated Stokes leakage approximately by imposing its area average, with reversed normal sign, on a separate P1 FV1 Darcy pressure problem. Both blocks use BiCGStab with a rediscretized GMG V-cycle and an explicitly agglomerated iterative coarse solve.', PALE_CYAN, CYAN)]

doc.build(story)
print(OUT)
