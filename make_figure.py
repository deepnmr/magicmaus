"""Figure 1 for the magicmaus Application Note: pipeline schematic + benchmark."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8})

fig = plt.figure(figsize=(7.2, 3.1), dpi=300)
gsA = fig.add_axes([0.015, 0.02, 0.60, 0.96]); gsA.axis('off'); gsA.set_xlim(0, 10); gsA.set_ylim(0, 10)
axB = fig.add_axes([0.715, 0.20, 0.275, 0.70])

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

# --- Panel B: MBP benchmark (same intensity NOESY) ---
labels = ['MAGIC', 'MAUS', 'magic\nmaus', 'magicmaus\n+soft']
acc = [5.7, 26.6, 72.9, 79.7]
cols = [RED, GREEN, BLUE, BLUE]
env = [None, 100, 100, 100]
x = range(len(labels))
axB.bar(x, acc, color=cols, width=0.66, zorder=3)
for i, v in enumerate(acc):
  axB.text(i, v + 2.5, f'{v:.1f}', ha='center', fontsize=7, weight='bold', color=cols[i])
# envelope markers
axB.axhline(100, ls=':', lw=1.0, color='#999', zorder=1)
axB.text(len(labels) - 0.5, 101.5, 'truth-in-envelope = 100%', ha='right', fontsize=6.0, color='#666')
for i, e in enumerate(env):
  if e:
    axB.plot(i, 100, marker='v', ms=5, color='#2e9e6b', zorder=4)
axB.text(1, 22, 'decisive\nonly', ha='center', va='top', fontsize=5.8, color='#555')
axB.set_ylim(0, 112); axB.set_xticks(list(x)); axB.set_xticklabels(labels, fontsize=6.6)
axB.set_ylabel('methyl-level correct (%)', fontsize=7.2)
axB.set_title('B   MBP, 192 methyls (same intensity NOESY)', fontsize=7.6, weight='bold', loc='left')
axB.spines[['top', 'right']].set_visible(False)
axB.tick_params(labelsize=6.6)

fig.savefig('figure1.png', dpi=300, bbox_inches='tight', facecolor='white')
print('wrote figure1.png')
