"""Reproduce the review flow diagram (counts unchanged; PRISMA counting units explicit).

Run from mixed_paternity_final with Python and reportlab installed:
    python scripts/build_prisma.py
The full-text totals count reports, not distinct studies or analytical data rows.
"""
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, white

OUT = Path(__file__).resolve().parents[1] / "Figures/prisma_reporting.pdf"
OUT.parent.mkdir(parents=True, exist_ok=True)
identified, duplicates, screened, excluded_records = 2288, 623, 1665, 1476
assessed, excluded_reports, included = 189, 46, 143
reasons = [
    ("No relevant paternity data", 10),
    ("No genetic parentage data", 8),
    ("Overlapping data", 14),
    ("Non-eligible study design", 6),
    ("Not a full article", 5),
    ("Case study", 2),
    ("Other", 1),
]
assert identified - duplicates == screened
assert screened - excluded_records == assessed
assert assessed - excluded_reports == included
assert sum(n for _, n in reasons) == excluded_reports

c = canvas.Canvas(str(OUT), pagesize=(842, 595))
c.setTitle("Literature search and selection: records and reports")
c.setAuthor("Avian mixed-paternity project")
blue = HexColor("#1c6580")
ink = HexColor("#17252d")

def box(x, y, w, h, lines, size=16):
    c.setStrokeColor(ink)
    c.setFillColor(white)
    c.setLineWidth(1.2)
    c.rect(x, y, w, h, stroke=1, fill=1)
    c.setFillColor(ink)
    c.setFont("Helvetica", size)
    leading = size * 1.3
    start = y + h / 2 + (len(lines) - 1) * leading / 2 - size * .34
    for i, line in enumerate(lines):
        assert c.stringWidth(line, "Helvetica", size) < w - 16
        c.drawCentredString(x + w / 2, start - i * leading, line)

def arrow(x1, y1, x2, y2):
    c.setStrokeColor(blue)
    c.setFillColor(blue)
    c.setLineWidth(2.3)
    c.line(x1, y1, x2, y2)
    p = c.beginPath()
    if x1 == x2:
        p.moveTo(x2, y2); p.lineTo(x2 - 4, y2 + 8); p.lineTo(x2 + 4, y2 + 8)
    else:
        p.moveTo(x2, y2); p.lineTo(x2 - 8, y2 - 4); p.lineTo(x2 - 8, y2 + 4)
    p.close(); c.drawPath(p, fill=1, stroke=0)

def stage(label, y, h):
    c.setFillColor(blue); c.rect(30, y, 38, h, stroke=0, fill=1)
    c.saveState(); c.translate(55, y + h / 2); c.rotate(90)
    c.setFillColor(white); c.setFont("Helvetica", 15)
    c.drawCentredString(0, 0, label); c.restoreState()

stage("Identification", 487, 93)
stage("Screening", 160, 310)
stage("Included", 49, 100)
box(100, 500, 310, 70, ["Records identified through", "database searching (n = 2,288)"])
box(490, 500, 320, 70, ["Duplicate records removed", "(n = 623)"])
box(100, 390, 310, 70, ["Records screened by title", "and abstract (n = 1,665)"])
box(490, 390, 320, 70, ["Records excluded", "(n = 1,476)"])
box(100, 275, 310, 70, ["Reports assessed for eligibility", "(n = 189)"])
box(100, 65, 310, 70, ["Reports included in the", "updated review (n = 143)"])
box(490, 170, 320, 175, [])
c.setFont("Helvetica-Bold", 15); c.drawString(504, 322, "Reports excluded (n = 46)")
c.setFont("Helvetica", 13)
for i, (reason, n) in enumerate(reasons):
    c.drawString(504, 298 - i * 17, f"{reason} (n = {n})")
for top, bottom in [(500, 460), (390, 345), (275, 135)]:
    arrow(255, top, 255, bottom)
arrow(410, 535, 490, 535)
arrow(410, 425, 490, 425)
arrow(410, 310, 490, 310)
c.setFillColor(ink); c.setFont("Helvetica", 10.5)
c.drawString(100, 28, "Full-text counts refer to reports, not unique studies or analytical database records.")
c.showPage(); c.save()
print(OUT)
