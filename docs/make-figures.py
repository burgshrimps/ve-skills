import json, math
from pathlib import Path

DOCS = Path("/home/vinzenzmay/dev/nebius_clawbio/ve-skills/docs")
FONT = "-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif"
INK, INK2, MUT, GRID, RULE = "#1F1D1A", "#5A544C", "#8A8278", "#DCD8D0", "#BBB5AA"
LOSS, GAIN, NEUT = "#d03b3b", "#2b6ca3", "#8A8278"
LOSS_SOFT = "#f7e4e2"
def esc(s): return s.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

# ---------------------------------------------------------------- chart 1
W, H = 760, 292
L, R = 74, 742
AX = 196                     # axis y
lo, hi = 1e-4, 1.0
def x_of(v): return L + (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * (R - L)

s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" role="img" '
     f'aria-label="Log-scale allele-frequency axis. The three rarity thresholds in common use (0.1%, 1%, 5%) '
     f'all sit below 12.5%, which is the smallest non-zero allele frequency a four-sample VCF can contain, '
     f'so all eight possible values are classed as common.">',
     f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>',
     f'<g font-family="{FONT}">']
s.append(f'<text x="{L}" y="26" font-size="14" font-weight="600" fill="{INK}">INFO/AF cannot be rare</text>')
s.append(f'<text x="{L}" y="44" font-size="10.5" fill="{MUT}">4 diploid samples → AN = 8 → the smallest non-zero AF the file can contain is 1/8</text>')
s.append(f'<text x="{L}" y="60" font-size="10" fill="{MUT}">dashed: rarity thresholds in use — 0.1% strict dominant · 1% standard rare-disease cut · 5% ACMG BA1</text>')

# discard band
bx = x_of(0.125)
s.append(f'<rect x="{bx:.1f}" y="68" width="{R-bx:.1f}" height="{AX-68:.1f}" fill="{LOSS_SOFT}"/>')
s.append(f'<text x="{bx-14:.1f}" y="86" text-anchor="end" font-size="11" font-weight="600" fill="{LOSS}">the only values INFO/AF can take  →</text>')
s.append(f'<text x="{bx-14:.1f}" y="101" text-anchor="end" font-size="10" fill="{INK2}">every one classed COMMON → discarded</text>')

# axis + ticks
s.append(f'<line x1="{L}" y1="{AX}" x2="{R}" y2="{AX}" stroke="{RULE}" stroke-width="1"/>')
for v, lab in [(1e-4,"0.01%"),(1e-3,"0.1%"),(1e-2,"1%"),(1e-1,"10%"),(1.0,"100%")]:
    x = x_of(v)
    s.append(f'<line x1="{x:.1f}" y1="{AX}" x2="{x:.1f}" y2="{AX+5}" stroke="{RULE}" stroke-width="1"/>')
    s.append(f'<text x="{x:.1f}" y="{AX+19}" text-anchor="middle" font-size="10" fill="{MUT}">{lab}</text>')
s.append(f'<text x="{(L+R)/2:.1f}" y="{AX+38}" text-anchor="middle" font-size="10.5" fill="{INK2}">population allele frequency (log scale)</text>')

# thresholds
for v, lab, sub in [(0.001,"0.1%","strict dominant"),(0.01,"1%","standard rare-disease cut"),(0.05,"5%","ACMG BA1")]:
    x = x_of(v)
    s.append(f'<line x1="{x:.1f}" y1="100" x2="{x:.1f}" y2="{AX}" stroke="{RULE}" stroke-width="1" stroke-dasharray="3 3"/>')
    s.append(f'<text x="{x:.1f}" y="{AX-8}" text-anchor="middle" font-size="10" font-weight="600" fill="{INK2}">{lab}</text>')
    _ = sub  # descriptors moved to the caption line; rotated text collided with the axis

# the 8 possible values
for k in range(1, 9):
    v = k/8; x = x_of(v)
    s.append(f'<circle cx="{x:.1f}" cy="{AX}" r="5" fill="{LOSS}" stroke="#FFFFFF" stroke-width="2"/>')
s.append(f'<text x="{x_of(0.125):.1f}" y="{AX-30}" text-anchor="middle" font-size="12" font-weight="700" fill="{LOSS}">12.5%</text>')
s.append(f'<text x="{x_of(0.125):.1f}" y="{AX-16}" text-anchor="middle" font-size="9.5" fill="{LOSS}">floor = 1/8</text>')
s.append(f'<text x="{x_of(1.0):.1f}" y="{AX-16}" text-anchor="middle" font-size="9.5" fill="{MUT}">8/8</text>')

s.append(f'<line x1="{L}" y1="{H-40}" x2="{R}" y2="{H-40}" stroke="{GRID}" stroke-width="1"/>')
s.append(f'<text x="{L}" y="{H-22}" font-size="11" fill="{INK2}">Read INFO/AF as population frequency and all 68 records score as common — including all 68 that segregate.</text>')
s.append(f'<text x="{L}" y="{H-8}" font-size="11" fill="{INK2}">ve-frequency returns <tspan font-weight="600" fill="{INK}">NO_DATA</tspan> instead: cohort AF, not population AF.</text>')
s.append("</g></svg>")
(DOCS/"frequency-floor.svg").write_text("\n".join(s), encoding="utf-8")
print("wrote frequency-floor.svg")

# ---------------------------------------------------------------- chart 2
res = json.load(open("/tmp/vs_demo/result.json"))["results"]
rows = []
for r in res:
    e = r.get("evidence") or {}
    rows.append(dict(gene=e.get("gene"), key=r["variant_key"], d=e.get("delta_bits"),
                     ref=e.get("mes_ref_bits"), alt=e.get("mes_alt_bits"),
                     dirn=r.get("direction"), ab=r.get("abstain_reason")))
scored = sorted([r for r in rows if r["d"] is not None], key=lambda r: r["d"])
absta  = [r for r in rows if r["d"] is None]

ROW, TOP = 23, 104
W2 = 760
H2 = TOP + ROW*(len(scored)+len(absta)) + 74
XL, XR = 250, 660
dmin, dmax = -10.0, 6.0
def x2(v): return XL + (v - dmin)/(dmax - dmin)*(XR - XL)
COL = {"damaging": LOSS, "strengthening": GAIN, "neutral": NEUT}

s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W2} {H2}" width="{W2}" height="{H2}" role="img" '
     f'aria-label="SnpEff labels all thirteen variants HIGH. MaxEntScan measures a change in splice-site strength '
     f'ranging from minus 8.75 bits to plus 5.23 bits, plus two abstentions, turning one label into a ranking.">',
     f'<rect width="{W2}" height="{H2}" fill="#FFFFFF"/>',
     f'<g font-family="{FONT}">']
s.append(f'<text x="40" y="26" font-size="14" font-weight="600" fill="{INK}">One HIGH label, actually measured</text>')
s.append(f'<text x="40" y="44" font-size="10.5" fill="{MUT}">MaxEntScan Δ in splice-site strength, ref → alt. Negative = the motif is weaker.</text>')
s.append(f'<text x="40" y="{TOP-26}" font-size="9.5" font-weight="600" fill="{MUT}" letter-spacing="0.08em">SNPEFF</text>')
s.append(f'<text x="104" y="{TOP-26}" font-size="9.5" font-weight="600" fill="{MUT}" letter-spacing="0.08em">GENE</text>')
s.append(f'<text x="{XL}" y="{TOP-38}" font-size="9.5" font-weight="600" fill="{MUT}" letter-spacing="0.08em">Δ BITS (MaxEntScan)</text>')

for v in range(-10, 7, 2):
    x = x2(v)
    s.append(f'<line x1="{x:.1f}" y1="{TOP-16}" x2="{x:.1f}" y2="{TOP+ROW*len(scored)-6}" stroke="{GRID}" stroke-width="1"/>')
    s.append(f'<text x="{x:.1f}" y="{TOP-20}" text-anchor="middle" font-size="9" fill="{MUT}">{v}</text>')
x0 = x2(0)
s.append(f'<line x1="{x0:.1f}" y1="{TOP-16}" x2="{x0:.1f}" y2="{TOP+ROW*len(scored)-6}" stroke="{RULE}" stroke-width="1.5"/>')

y = TOP
for r in scored:
    c = COL.get(r["dirn"], NEUT); xd = x2(r["d"])
    s.append(f'<text x="40" y="{y+4}" font-size="9.5" font-weight="600" fill="{MUT}">HIGH</text>')
    s.append(f'<text x="104" y="{y+4}" font-size="11" fill="{INK}">{esc(str(r["gene"]))}</text>')
    s.append(f'<line x1="{x0:.1f}" y1="{y}" x2="{xd:.1f}" y2="{y}" stroke="{c}" stroke-width="2"/>')
    s.append(f'<circle cx="{xd:.1f}" cy="{y}" r="4.5" fill="{c}" stroke="#FFFFFF" stroke-width="2"/>')
    anc, lx = ("end", xd-10) if r["d"] < 0 else ("start", xd+10)
    s.append(f'<text x="{lx:.1f}" y="{y+4}" text-anchor="{anc}" font-size="10" fill="{INK2}">{r["d"]:+.2f}</text>')
    if r["gene"] == "POLR3C":
        s.append(f'<text x="{x0+12:.1f}" y="{y+4}" font-size="9.5" fill="{LOSS}">canonical GT → AT — 94% of the motif lost</text>')
    y += ROW

s.append(f'<line x1="40" y1="{y-8}" x2="{XR}" y2="{y-8}" stroke="{GRID}" stroke-width="1" stroke-dasharray="3 3"/>')
for r in absta:
    s.append(f'<text x="40" y="{y+7}" font-size="9.5" font-weight="600" fill="{MUT}">HIGH</text>')
    s.append(f'<text x="104" y="{y+7}" font-size="11" fill="{MUT}">{esc(r["key"])}</text>')
    s.append(f'<circle cx="{XL-16:.1f}" cy="{y+3}" r="4.5" fill="#FFFFFF" stroke="{MUT}" stroke-width="1.5" stroke-dasharray="2 2"/>')
    short = ("no reference GT donor within 20 nt — nothing to measure a delta against"
             if "no canonical GT donor" in (r["ab"] or "")
             else "routed protein_truncating — outside this skill's domain")
    s.append(f'<text x="{XL}" y="{y+7}" font-size="10" fill="{MUT}">abstained — {esc(short)}</text>')
    y += ROW

s.append(f'<line x1="40" y1="{y+8}" x2="{XR}" y2="{y+8}" stroke="{GRID}" stroke-width="1"/>')
s.append(f'<text x="40" y="{y+28}" font-size="11" fill="{INK2}">'
         f'<tspan font-weight="600" fill="{LOSS}">7 weakened</tspan>  ·  '
         f'<tspan font-weight="600" fill="{GAIN}">1 strengthened</tspan>  ·  '
         f'<tspan font-weight="600" fill="{NEUT}">3 negligible</tspan>  ·  2 abstained '
         f'— a ranking where the label gave none.</text>')
s.append("</g></svg>")
(DOCS/"splice-delta.svg").write_text("\n".join(s), encoding="utf-8")
print(f"wrote splice-delta.svg  ({len(scored)} scored, {len(absta)} abstained, {H2}px tall)")
