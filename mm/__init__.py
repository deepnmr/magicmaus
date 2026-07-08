"""magicmaus — methyl-NMR assignment fusing MAUS (hard SAT bounds) and
MAGIC (distance-weighted NOE scoring).

Clean-room fresh build (does not reuse the repo's older maus.py / magicmaus.py).

  * MAUS idea  — cast assignment as subgraph isomorphism, solve as SAT, return
    for every HMQC peak the *set* of methyls consistent with all hard
    constraints.  Never excludes the truth when firm edges are correct; abstains
    under degeneracy.
  * MAGIC idea — score a global peak->methyl map with a distance-weighted NOE
    objective (intensity ~ 1/r^6) and commit to the best scorer.

magicmaus uses MAUS to bound the space with certainty, then MAGIC-style scoring
to commit a single coherent call within the residual degeneracy.

Inputs (per experiment):
  * HMQC          2D methyl peak list : label  H_ppm  C_ppm  res_type
  * 3D HMBC-HMQC  (C2, C1, H1)        : detected methyl = (C1,H1); C2 = geminal
                                        partner carbon (same-residue Leu/Val link)
  * 3D NOESY      (C2, C1, H1)        : detected methyl = (C1,H1); C2 = through-
                                        space partner carbon (NOE edge)

Structure is an mmCIF (or PDB); a homo-oligomer contributes one chain image per
subunit and inter-subunit NOEs are taken as the min distance over images.
"""

from .structure import Methyl, parse_structure, build_structure_graph, min_dist
from .peaks import (
    Peak, load_hmqc, load_triple, load_truth, match_noesy, match_hmbc,
    THREE_TO_ONE, METHYL_ATOMS, parse_labeling,
)
from .engine import MagicMaus, assign

__all__ = [
    'Methyl', 'parse_structure', 'build_structure_graph', 'min_dist',
    'Peak', 'load_hmqc', 'load_triple', 'load_truth', 'match_noesy', 'match_hmbc',
    'THREE_TO_ONE', 'METHYL_ATOMS', 'parse_labeling',
    'MagicMaus', 'assign',
]
