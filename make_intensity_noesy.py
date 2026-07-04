"""Add physically-realistic intensities to a simulated NOESY peak list.

The MBP example NOESY is a *boolean* network (a cross peak exists iff two
methyls are within a cutoff), which carries no which-contact information beyond
the hard edge -- so MAGIC-style scoring has nothing to grade.  A real NOESY
cross peak, though, has an intensity ~ 1/r^6 in the underlying methyl-methyl
distance.  This script reconstructs that intensity for every existing row using
the truth key and the structure, writing a 5-column list
``label C1 C2 H2 intensity`` that magicmaus can exploit.

For each row (C1, C2, H2): the observed methyl is the truth of the HMQC peak at
(H2, C2); the partner is the truth methyl -- among the peaks whose carbon shift
matches C1 -- closest to the observed methyl (the contact that actually gave
rise to the peak).  intensity = (d0 / d)^6, normalised so the closest contact in
the list is 1.0.  Usage:

    python make_intensity_noesy.py PDB HMQC_true.tsv NOESY.tsv OUT.tsv [tol_c]
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'maus'))
import maus  # noqa: E402


def main(argv):
  pdb, truth_path, noesy_path, out_path = argv[1:5]
  tol_c = float(argv[5]) if len(argv) > 5 else 0.05

  lab = maus.parse_labeling('A;I;L;M;T;V')
  methyls = maus.parse_structure(Path(pdb).read_text().splitlines(), lab)
  coord = {m.label: m.coord for m in methyls}
  peaks = maus.load_hmqc(str(Path(noesy_path).with_name('hmqc.tsv')))
  truth = maus.load_truth(truth_path)
  peak_coord = {p.index: coord.get(truth.get(p.peak_id)) for p in peaks}

  rows = []
  for line in Path(noesy_path).read_text().splitlines():
    if not line.strip() or line.startswith('#') or line.startswith('label'):
      continue
    f = line.split('\t')
    label, c1, c2, h2 = f[0], float(f[1]), float(f[2]), float(f[3])
    obs = [p.index for p in peaks
           if abs(p.h_ppm - h2) <= 0.01 and abs(p.c_ppm - c2) <= tol_c]
    if not obs or peak_coord[obs[0]] is None:
      rows.append((label, c1, c2, h2, 1.0)); continue
    m_obs = peak_coord[obs[0]]
    partners = [peak_coord[p.index] for p in peaks
                if abs(p.c_ppm - c1) <= tol_c and peak_coord[p.index] is not None]
    dists = [math.dist(m_obs, pc) for pc in partners if pc is not m_obs]
    d = min(dists) if dists else None
    rows.append((label, c1, c2, h2, 1.0 / (d ** 6) if d else 1.0))

  dmax = max(r[4] for r in rows) or 1.0            # normalise closest contact -> 1
  lines = ['label\tC1\tC2\tH2\tintensity']
  for (label, c1, c2, h2, inten) in rows:
    lines.append(f'{label}\t{c1:.3f}\t{c2:.3f}\t{h2:.3f}\t{inten / dmax:.4e}')
  Path(out_path).write_text('\n'.join(lines) + '\n')
  print(f'wrote {out_path}: {len(rows)} rows with reconstructed 1/r^6 intensities')


if __name__ == '__main__':
  raise SystemExit(main(sys.argv))
