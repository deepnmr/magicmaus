"""Figure 1 for the magicmaus Application Note: pipeline schematic + benchmark."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8})

fig = plt.figure(figsize=(7.2, 3.1), dpi=300)
gsA = fig.add_axes([0.015, 0.02, 0.55, 0.96]); gsA.axis('off'); gsA.set_xlim(0, 10); gsA.set_ylim(0, 10)
axB = fig.add_axes([0.655, 0.22, 0.335, 0.66])

BLUE, GREEN, ORANGE, GREY, RED = '#2c6fb5', '#2e9e6b', '#e08a2b', '#8a8f98', '#c0392b'


def rect(ax, x, y, w, h, fc, ec):
  ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.06,rounding_size=0.18',
                              linewidth=1.1, edgecolor=ec, facecolor=fc, zorder=2))


def tierbox(ax, x, y, w, h, title, sub, fc, ec, tc):
  rect(ax, x, y, w, h, fc, ec)
  ax.text(x + w / 2, y + h * 0.63, title, ha='center', va='center', fontsize=7.0,
          weight='bold', color=tc, zorder=3)
  ax.text(x + w / 2, y + h * 0.28, sub, ha='center', va='center', fontsize=6.2,
          color='#333', zorder=3)


def arrow(ax, x1, y1, x2, y2):
  ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=11,
                               lw=1.3, color='#444', zorder=1))


# --- Panel A: schematic ---
gsA.text(0.1, 9.6, 'A', fontsize=13, weight='bold')
# inputs
rect(gsA, 0.25, 6.55, 2.5, 2.5, '#eef4fb', BLUE)
gsA.text(1.5, 8.75, 'Inputs', ha='center', fontsize=7.4, weight='bold', color=BLUE)
gsA.text(1.5, 7.55, 'HMQC (2D peaks)\nNOESY (3D (H)CCH,\n1/r$^6$ intensity)\nPDB structure',
         ha='center', va='center', fontsize=6.6, color='#1a1a1a')
arrow(gsA, 2.85, 7.8, 3.5, 7.8)
# MAUS layer
rect(gsA, 3.55, 6.35, 3.15, 2.9, '#e8f6ee', GREEN)
gsA.text(5.12, 8.85, 'MAUS  (SAT)', ha='center', fontsize=8.2, weight='bold', color=GREEN)
gsA.text(5.12, 8.42, 'hard constraints', ha='center', fontsize=6.8, style='italic', color=GREEN)
gsA.text(5.12, 7.35, 'per-peak option set\ntruth never excluded\n(100% envelope)',
         ha='center', va='center', fontsize=6.6, color='#1a1a1a')
arrow(gsA, 5.12, 6.3, 5.12, 5.2)
gsA.text(5.35, 5.72, 'pruned domains\n(mostly 1–3 candidates)', ha='left', va='center', fontsize=6.0, color=GREY)
# MAGIC layer
rect(gsA, 3.55, 2.35, 3.15, 2.85, '#fdf0e2', ORANGE)
gsA.text(5.12, 4.82, 'MAGIC-style scoring', ha='center', fontsize=7.8, weight='bold', color=ORANGE)
gsA.text(5.12, 4.42, 'intensity · (1/r$^6$)', ha='center', fontsize=6.8, style='italic', color=ORANGE)
gsA.text(5.12, 3.35, 'SAT-feasible seed +\ncoordinate ascent\n→ single coherent map',
         ha='center', va='center', fontsize=6.6, color='#1a1a1a')
# output tiers
gsA.text(8.68, 6.35, 'single call + confidence', ha='center', fontsize=6.8, style='italic', color='#333')
tierbox(gsA, 7.55, 4.95, 2.25, 1.05, 'unique', 'forced, certain', '#eef4fb', BLUE, BLUE)
tierbox(gsA, 7.55, 3.30, 2.25, 1.05, 'scored', 'NOE-ranked', '#fdf0e2', ORANGE, ORANGE)
tierbox(gsA, 7.55, 1.65, 2.25, 1.05, 'ambiguous', 'true symmetry', '#f3f4f6', GREY, '#555')
arrow(gsA, 6.75, 4.3, 7.5, 5.45); arrow(gsA, 6.75, 3.8, 7.5, 3.82); arrow(gsA, 6.75, 3.3, 7.5, 2.2)

# --- Panel B: seven-target benchmark (same intensity NOESY) ---
# (target, methyls, MAGIC%, magicmaus+soft%, MAGIC-converged)
DATA = [
  ('Ubq', 43, 9.3, 90.7, True),
  ('HNH', 57, 12.3, 57.9, True),
  ('IL-2', 59, 8.5, 89.8, True),
  ('REC2', 63, None, 76.2, False),
  ('REC3', 85, None, 28.2, False),
  ('MBP', 192, 5.7, 79.7, True),
  ('MSG', 257, None, 33.5, False),
]
x = range(len(DATA))
w = 0.38
for i, (_lab, _n, mg, mm, conv) in enumerate(DATA):
  axB.bar(i - w / 2, mm, width=w, color=BLUE, zorder=3)
  if conv:
    axB.bar(i + w / 2, mg, width=w, color=RED, zorder=3)
  else:
    axB.bar(i + w / 2, 4.0, width=w, color='none', edgecolor=RED, hatch='///', lw=0.6, zorder=3)
    axB.text(i + w / 2, 5.5, 'n.c.', ha='center', fontsize=5.0, color=RED, rotation=90, va='bottom')
  axB.plot(i - w / 2, 100, marker='v', ms=4, color=GREEN, zorder=4)
axB.axhline(100, ls=':', lw=0.9, color='#999', zorder=1)
axB.text(len(DATA) - 0.5, 103, 'envelope = 100%', ha='right', fontsize=5.6, color='#666')
# legend
axB.bar(-9, 0, color=BLUE, label='magicmaus +soft')
axB.bar(-9, 0, color=RED, label='MAGIC')
axB.legend(loc='upper center', fontsize=5.6, frameon=False, ncol=2,
           bbox_to_anchor=(0.5, 1.18), handlelength=1.0, columnspacing=1.0)
axB.set_ylim(0, 108); axB.set_xlim(-0.6, len(DATA) - 0.4)
axB.set_xticks(list(x))
axB.set_xticklabels([f'{lab}\n{n}' for (lab, n, *_ ) in DATA], fontsize=5.6)
axB.set_ylabel('methyl-level correct (%)', fontsize=7.0)
axB.set_title('B   seven targets (methyls below)', fontsize=7.4, weight='bold', loc='left')
axB.spines[['top', 'right']].set_visible(False)
axB.tick_params(labelsize=6.2)

fig.savefig('figure1.png', dpi=300, bbox_inches='tight', facecolor='white')
print('wrote figure1.png')
