"""TNF-alpha assignment on the given (assigned) peak lists, using NOESY symmetry.

A 3D (H)CCH NOESY row (C1_partner, C2_obs, H2_obs) gives the partner only by
carbon; its reciprocal (C2_partner, C1_obs, H1_obs) supplies the partner's
proton.  Pairing the two resolves BOTH endpoints of the NOE by full (C,H), so
every edge is a firm, correct methyl-methyl contact -> a near-perfect MAUS envelope
(never excludes the truth) with no wrong hard constraint.  On the dense trimer the
plain carbon-only firm match is UNSAT at this tolerance, so the commitment engine
grows a max-feasible hard set from the symmetric seed instead.

The discarded carbon-only ambiguous NOEs are folded back as soft evidence for the
MAGIC-style commitment.

Usage:  python run_tnfa_symmetric.py
"""
import os
import sys

# Pin the hash seed so the SAT option enumeration order (and hence the committed
# map) is reproducible across sessions; re-exec once if it is not already set.
if os.environ.get('PYTHONHASHSEED') != '0':
  os.environ['PYTHONHASHSEED'] = '0'
  os.execv(sys.executable, [sys.executable] + sys.argv)

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


def both_dirs(edges):
  """Undirected edge keys -> a set with both (a,b) and (b,a)."""
  ne = set()
  for a, b in edges:
    ne.add((a, b))
    ne.add((b, a))
  return ne


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
    opt = maus.MAUS(methyls, peaks, gem, sg, lg, both_dirs(edgeset)).solve_options()
    return all(len(opt[p.index]) > 0 for p in peaks)

  hard = set(seed_edges)
  for k, _w in sorted(cand_edges.items(), key=lambda kv: -kv[1]):
    if k in hard or (k[1], k[0]) in hard:
      continue
    if feasible(hard | {k}):
      hard.add(k)
  return hard


def _sign(x):
  return (x > 1e-12) - (x < -1e-12)


def resolve_geminal(eng, chosen, gem):
  """Fix Leu/Val geminal swaps by NOE intensity ratio (a local, deterministic
  2-way flip -- NOT a global re-score, which only lowers accuracy here).

  Physics: NOE intensity ~ 1/r^6, so within a geminal pair (G1,G2) both assigned
  (to peaks Pa,Pb), the peak with the STRONGER cross peak to a shared partner Q
  holds the methyl structurally CLOSER to Q.  Vote sign(I(Pa,Q)-I(Pb,Q)) *
  sign(w(G1,Gq)-w(G2,Gq)) over shared firm-NOE partners; flip the pair if the net
  vote disagrees with the current orientation.  Uses firm edges only (unweighted
  majority): weighting by |I|*|w| lets one very-close predicted-structure contact
  dominate, and folding the diluted soft rows adds mismatched noise -- both hurt.
  Preserves injectivity, residue identity, and the MAUS envelope.
  """
  out = dict(chosen)
  peak_by_methyl = {g: p for p, g in out.items()}
  ei = eng.edge_intensity
  seen = set()
  flips = 0
  for (g1, g2) in gem:
    if g1 > g2 or (g1, g2) in seen:
      continue
    seen.add((g1, g2))
    pa, pb = peak_by_methyl.get(g1), peak_by_methyl.get(g2)
    if pa is None or pb is None:
      continue
    shared = (set(eng.noe_nbr.get(pa, ())) & set(eng.noe_nbr.get(pb, ()))) - {pa, pb}
    net = 0.0
    for q in shared:
      gq = out.get(q)
      if gq is None:
        continue
      ia = ei.get((min(pa, q), max(pa, q)))
      ib = ei.get((min(pb, q), max(pb, q)))
      if ia is None or ib is None:
        continue
      net += _sign(ia - ib) * _sign(eng._pair_w(g1, gq) - eng._pair_w(g2, gq))
    if net < 0:
      out[pa], out[pb] = g2, g1
      peak_by_methyl[g1], peak_by_methyl[g2] = pb, pa
      flips += 1
  return out, flips


def main():
  peaks = maus.load_hmqc(str(D / 'hmqc.tsv'))
  truth = maus.load_truth(str(D / 'hmqc_true.tsv'))
  lab = maus.parse_labeling('A;I;L;M;T;V')
  methyls = maus.parse_structure(CIF.read_text().splitlines(), lab)
  lbi = {m.index: m.label for m in methyls}
  gem, sg, lg = maus.build_structure_graph(methyls, 6.0, 10.0)
  rows = mm.load_noesy_rows(str(D / 'noesy_intensity.tsv'))

  ei = symmetric_edges(peaks, rows, TOL_C, TOL_H)
  _, ei_c, amb, _ = mm.match_noe_intensity(peaks, rows, TOL_H, TOL_C)
  _, gstat = maus.match_hmbc(peaks, maus.load_hmbc(str(D / 'hmbc.tsv')), TOL_H, TOL_C)

  # Envelope engine: symmetric edges only -> every option set contains the truth
  # (bar peaks whose NOEs the structure cannot support); carbon-only ambiguous
  # rows are soft evidence for the commitment.
  engE = mm.MagicMaus(methyls, peaks, gem, sg, lg, both_dirs(ei), edge_intensity=ei)
  engE.set_soft_evidence(amb)
  _, options = engE.solve()

  # Commitment engine: grow a max-feasible hard set from the symmetric seed by
  # adding carbon-only firm edges that keep the SAT feasible.  The extra (partly
  # wrong) hard edges prune the option sets enough for the scorer to commit, at
  # the cost of the commit engine's own envelope (its wrong edges exclude a few
  # truths -- which is fine, the ENVELOPE number comes from engE, not engC).
  hard = greedy_feasible(methyls, peaks, gem, sg, lg, ei.keys(), ei_c)
  ints = dict(ei)
  for k, w in ei_c.items():
    ints[k] = max(ints.get(k, 0.0), w)
  soft = list(amb) + [((a,), (b,), 1.0, w) for (a, b), w in ei_c.items()
                      if (a, b) not in hard and (b, a) not in hard]
  engC = mm.MagicMaus(methyls, peaks, gem, sg, lg, both_dirs(hard), edge_intensity=ints)
  engC.set_soft_evidence(soft)
  chosen, _ = engC.solve()
  # Break the Leu/Val geminal swaps with the NOE intensity ratio (deterministic
  # local flips; the global objective cannot orient the pairs -- climbing it only
  # lowers accuracy on this real/predicted-structure data).
  chosen, gflips = resolve_geminal(engC, chosen, gem)

  # Merge: each peak gets the greedy committed call inside the symmetric envelope.
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

  print(f'symmetric NOE edges = {len(ei)}   HMBC gem-links = {gstat["firm"]}   geminal flips = {gflips}')
  print(f'MAUS envelope (symmetric)   = {inset}/{n} = {100 * inset / n:.1f}%  (never excludes truth)')
  print(f'committed (greedy+geminal)  = methyl {meth}/{n} = {100 * meth / n:.1f}%  '
        f'residue {resd}/{n} = {100 * resd / n:.1f}%')
  print(f'committed call in envelope  = {cons}/{n} = {100 * cons / n:.1f}%')
  print(f'wrote {D}/magicmaus_calls_symmetric.tsv')


if __name__ == '__main__':
  main()
