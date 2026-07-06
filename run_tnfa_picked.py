"""Run magicmaus on the peak-picked TNF-alpha spectra and score against the
known assignment by ppm (picked peak labels do not align with the truth key)."""
import re
from pathlib import Path

import maus
import magicmaus as mm

D = Path('examples/TNFa')
CIF = D / 'fold_tnfa_trimer_model_0.cif'
TOL_H, TOL_C = 0.02, 0.1


def residue(label):
  m = re.match(r'^([A-Z])(\d+)', label or '')
  return (m.group(1), int(m.group(2))) if m else None


def load_truth_peaks():
  out = []
  for l in (D / 'hmqc_true.tsv').read_text().splitlines()[1:]:
    f = l.split('\t')
    out.append((float(f[1]), float(f[2]), f[3], f[4]))   # h, c, type, label
  return out


def main(hmqc='hmqc_picked.tsv', noesy='noesy_picked.tsv'):
  peaks = maus.load_hmqc(str(D / hmqc))
  lab = maus.parse_labeling('A;I;L;M;T;V')
  methyls = maus.parse_structure(CIF.read_text().splitlines(), lab)
  lbi = {m.index: m.label for m in methyls}
  crosses = mm.load_noesy_rows(str(D / noesy))
  gem, sg, lg = maus.build_structure_graph(methyls, 6.0, 10.0)
  noe, ei, amb, stats = mm.match_noe_intensity(peaks, crosses, TOL_H, TOL_C)
  print(f'peaks={len(peaks)}  NOE match: {stats}')

  eng = mm.MagicMaus(methyls, peaks, gem, sg, lg, noe, edge_intensity=ei)
  eng.set_soft_evidence(amb)
  chosen, options = eng.solve()

  # map each picked peak -> nearest truth peak by ppm
  truth = load_truth_peaks()
  used = [False] * len(truth)
  meth = resd = inset = mapped = 0
  for p in peaks:
    best = None
    for ti, (th, tc, tt, tl) in enumerate(truth):
      if used[ti] or abs(p.c_ppm - tc) > TOL_C or abs(p.h_ppm - th) > TOL_H:
        continue
      d = abs(p.c_ppm - tc) + abs(p.h_ppm - th)
      if best is None or d < best[0]:
        best = (d, ti, tl)
    if not best:
      continue
    used[best[1]] = True
    mapped += 1
    true_label = best[2]
    call = lbi.get(chosen[p.index])
    opt = [lbi[g] for g in options[p.index]]
    if call == true_label:
      meth += 1
    if residue(call) == residue(true_label):
      resd += 1
    if true_label in opt:
      inset += 1
  nt = len(truth)
  print(f'mapped {mapped}/{nt} picked->true')
  print(f'envelope (truth in option set) = {inset}/{mapped} = {100*inset/mapped:.1f}%')
  print(f'magicmaus methyl = {meth}/{mapped} = {100*meth/mapped:.1f}%  '
        f'residue = {resd}/{mapped} = {100*resd/mapped:.1f}%')


if __name__ == '__main__':
  main()
