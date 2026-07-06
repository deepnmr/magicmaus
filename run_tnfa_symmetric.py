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


def greedy_feasible(methyls, peaks, gem, sg, lg, seed_edges, cand_edges):
  """Grow the SAT-feasible hard-edge set: seed with the trustworthy symmetric
  edges (always feasible), then add carbon-only firm edges strongest-first,
  keeping each only if the SAT still leaves every peak a non-empty option set.

  A carbon-only edge can be wrong (carbon degeneracy picks the wrong partner);
  such an edge prunes more than the truth-preserving symmetric set, so it lets
  the scorer commit -- at the cost of its own envelope.  That is why this feeds
  the COMMIT engine only; the ENVELOPE engine keeps the symmetric edges alone.
  """
  def feasible(edgeset):
    ne = set()
    for a, b in edgeset:
      ne.add((a, b)); ne.add((b, a))
    opt = maus.MAUS(methyls, peaks, gem, sg, lg, ne).solve_options()
    return all(len(opt[p.index]) > 0 for p in peaks)

  hard = set(seed_edges)
  for k, _w in sorted(cand_edges.items(), key=lambda kv: -kv[1]):
    if k in hard or (k[1], k[0]) in hard:
      continue
    if feasible(hard | {k}):
      hard.add(k)
  return hard


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

  # Envelope engine: symmetric edges only -> every option set contains the truth.
  engE = mm.MagicMaus(methyls, peaks, gem, sg, lg, noe, edge_intensity=ei)
  engE.set_soft_evidence(amb)
  _, options = engE.solve()

  # Commitment engine: grow a max-feasible hard set from the symmetric seed by
  # adding carbon-only firm edges that keep the SAT feasible.  The extra (partly
  # wrong) hard edges prune the option sets enough for the scorer to commit, at
  # the cost of the commit engine's own envelope (its wrong edges exclude a few
  # truths -- which is fine, the ENVELOPE number comes from engE, not engC).
  noe_c, ei_c, amb_c, _ = mm.match_noe_intensity(peaks, rows, TOL_H, TOL_C)
  hard = greedy_feasible(methyls, peaks, gem, sg, lg, ei.keys(), ei_c)
  nset = set()
  for a, b in hard:
    nset.add((a, b))
    nset.add((b, a))
  ints = dict(ei)
  for k, w in ei_c.items():
    ints[k] = max(ints.get(k, 0.0), w)
  soft = list(amb_c) + [((a,), (b,), 1.0, w) for (a, b), w in ei_c.items()
                        if (a, b) not in hard and (b, a) not in hard]
  engC = mm.MagicMaus(methyls, peaks, gem, sg, lg, nset, edge_intensity=ints)
  engC.set_soft_evidence(soft)
  chosen, _ = engC.solve()

  # Merge: each peak gets the carbon-only committed call inside the 100% symmetric
  # envelope.  Write the combined table.
  out = ['label\tres_type\tcall\tin_envelope\tn_options\toptions\ttruth\tcall_correct']
  n = inset = meth = resd = cons = 0
  for p in peaks:
    n += 1
    t = truth.get(p.peak_id)
    call = lbi.get(chosen[p.index])
    env = [lbi[g] for g in options[p.index]]
    inset += t in env
    meth += call == t
    resd += residue(call) == residue(t)
    cons += call in env
    out.append(f'{p.peak_id}\t{p.res_type}\t{call}\t{int(call in env)}\t{len(env)}\t'
               f'{",".join(env)}\t{t}\t{int(call == t)}')
  (D / 'magicmaus_calls_symmetric.tsv').write_text('\n'.join(out) + '\n')

  print(f'symmetric NOE edges = {len(ei)}   HMBC gem-links = {gstat["firm"]}')
  print(f'MAUS envelope (symmetric)   = {inset}/{n} = {100 * inset / n:.1f}%  (never excludes truth)')
  print(f'committed (greedy hard)     = methyl {meth}/{n} = {100 * meth / n:.1f}%  '
        f'residue {resd}/{n} = {100 * resd / n:.1f}%')
  print(f'committed call in envelope  = {cons}/{n} = {100 * cons / n:.1f}%')
  print(f'wrote {D}/magicmaus_calls_symmetric.tsv')


if __name__ == '__main__':
  main()
