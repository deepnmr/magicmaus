"""TNF-alpha assignment on the given (assigned) peak lists, using NOESY symmetry.

A 3D (H)CCH NOESY row (C1_partner, C2_obs, H2_obs) gives the partner only by
carbon; its reciprocal (C2_partner, C1_obs, H1_obs) supplies the partner's
proton.  Pairing the two resolves BOTH endpoints of the NOE by full (C,H), so
every edge is a firm, correct methyl-methyl contact -> the MAUS envelope reaches
100% (never excludes the truth) with no wrong hard constraint.

The discarded carbon-only ambiguous NOEs are folded back as soft evidence for the
MAGIC-style commitment.  The given HMBC-HMQC list adds geminal hard links.

Usage:  python run_tnfa_symmetric.py
"""
import re
from pathlib import Path

import maus
import magicmaus as mm

D = Path('examples/TNFa')
CIF = D / 'fold_tnfa_trimer_model_0.cif'
TOL_H, TOL_C = 0.01, 0.05


def residue(label):
  m = re.match(r'^([A-Z])(\d+)', label or '')
  return (m.group(1), int(m.group(2))) if m else None


def symmetric_edges(peaks, rows, tol_c, tol_h):
  """Firm NOE edges resolved by full (C,H) on both ends via reciprocal pairing."""
  def uniq(c, h):
    cand = [p.index for p in peaks
            if abs(p.c_ppm - c) <= tol_c and abs(p.h_ppm - h) <= tol_h]
    return cand[0] if len(cand) == 1 else None

  ed = {}
  for (c1, c2, h2, it) in rows:
    a = uniq(c2, h2)                       # observed methyl of this row
    if a is None:
      continue
    for (cc1, cc2, hh, ii) in rows:        # its reciprocal supplies the partner proton
      if abs(cc2 - c1) <= tol_c and abs(cc1 - c2) <= tol_c:
        q = uniq(cc2, hh)                  # partner methyl, now by (C,H)
        if q is not None and q != a:
          key = tuple(sorted((a, q)))
          ed[key] = max(ed.get(key, 0.0), min(it, ii))
  return ed


def main():
  peaks = maus.load_hmqc(str(D / 'hmqc.tsv'))
  truth = maus.load_truth(str(D / 'hmqc_true.tsv'))
  lab = maus.parse_labeling('A;I;L;M;T;V')
  methyls = maus.parse_structure(CIF.read_text().splitlines(), lab)
  lbi = {m.index: m.label for m in methyls}
  gem, sg, lg = maus.build_structure_graph(methyls, 6.0, 10.0)
  rows = mm.load_noesy_rows(str(D / 'noesy_intensity.tsv'))

  ei = symmetric_edges(peaks, rows, TOL_C, TOL_H)
  noe = set()
  for (a, b) in ei:
    noe.add((a, b))
    noe.add((b, a))
  # carbon-only ambiguous rows -> soft evidence for the commitment
  _, _, amb, _ = mm.match_noe_intensity(peaks, rows, TOL_H, TOL_C)
  gl, gstat = maus.match_hmbc(peaks, maus.load_hmbc(str(D / 'hmbc.tsv')), TOL_H, TOL_C)
  glset = set()
  for a, b in gl:
    glset.add((a, b))
    glset.add((b, a))

  eng = mm.MagicMaus(methyls, peaks, gem, sg, lg, noe,
                     gem_links=glset or None, edge_intensity=ei)
  eng.set_soft_evidence(amb)
  chosen, options = eng.solve()

  n = len(peaks)
  inset = sum(1 for p in peaks if truth.get(p.peak_id) in [lbi[g] for g in options[p.index]])
  meth = sum(1 for p in peaks if lbi.get(chosen[p.index]) == truth.get(p.peak_id))
  resd = sum(1 for p in peaks if residue(lbi.get(chosen[p.index])) == residue(truth.get(p.peak_id)))
  print(f'symmetric NOE edges = {len(ei)}   HMBC gem-links = {gstat["firm"]}')
  print(f'MAUS envelope (assignment) = {inset}/{n} = {100 * inset / n:.1f}%')
  print(f'committed methyl = {meth}/{n} = {100 * meth / n:.1f}%   '
        f'residue = {resd}/{n} = {100 * resd / n:.1f}%')


if __name__ == '__main__':
  main()
