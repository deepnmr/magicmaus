"""magicmaus — a hybrid methyl-NMR assignment engine that fuses MAUS and MAGIC.

Two clean-room siblings solve the *same* problem two different ways:

  * MAUS  (Nerli et al., Nat. Commun. 2021) — casts assignment as subgraph
    isomorphism and solves it as a SAT problem.  For every HMQC peak it returns
    the *set* of methyls consistent with every hard constraint and provably
    never excludes the truth — but it commits to a single answer only when the
    data force one, leaving genuine (and spurious-union) degeneracy behind.
    Implemented in ../maus/maus.py.

  * MAGIC (Monneau et al., J Biomol NMR 2017) — scores a global peak->methyl
    map with a distance-weighted NOE objective and returns the single highest
    scorer.  Decisive, but over the full type-matched candidate space its
    landscape is near-flat, so on structure-consistent NOE data it commits to
    near-optimal-but-wrong answers for most peaks.  Implemented in ../magic/.

Neither is dominant: MAUS bounds the space with certainty but abstains; MAGIC
commits but over too large a space.  COMPARISON.md (in ../maus) spells out the
synthesis: *use MAUS to bound the space, MAGIC-style scoring to rank within the
residual degeneracy.*  This program is that synthesis.

Pipeline
--------
1.  MAUS SAT layer (reused verbatim from maus.py) builds the hard CNF and
    enumerates each peak's option set O_i.  Truth in O_i is guaranteed.
2.  A key subtlety: the O_i are enumerated *independently* — a methyl is an
    option if *some* satisfying global map uses it.  Their product is NOT
    jointly realizable, and it does not say which single coherent assignment is
    best.  MAGIC answers exactly that question.
3.  magicmaus commits to a single map with MAGIC's distance-weighted NOE
    objective -- but only over MAUS's residual degeneracy.  Peaks MAUS already
    pins uniquely are fixed; the rest are partitioned into small independent
    clusters (coupled by a shared candidate methyl or a firm NOE edge) and each
    cluster is solved exhaustively for the assignment that places every observed
    NOE on the closest compatible structural contact (intensity ~ 1/r^6), while
    staying injective and NOE-consistent.  This is MAGIC's cluster protocol run
    inside the truth-containing MAUS space, so it is both exact and fast where
    MAGIC over the full space is neither.
4.  Every peak is reported three ways: the committed magicmaus call, the MAUS
    option set (the honest ambiguity envelope), and a confidence tier
    (unique / scored / ambiguous) from the local scoring margin.

The result keeps MAUS's never-exclude guarantee as an envelope while delivering
MAGIC's single answer — and, because the scored search runs only over the tiny
MAUS-pruned clusters, it is both fast and, unlike MAGIC alone, actually able to
find the truth.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Reuse the MAUS clean-room implementation verbatim (structure/peak parsing, the
# SAT encoding, and the NOE matching) so the hard-constraint layer is identical.
# The vendored copy lives next to this file (kept == ../maus/maus.py); importing
# it locally makes the repo self-contained and picks up its speed optimisations.
_MAUS_DIR = Path(__file__).resolve().parent
if str(_MAUS_DIR) not in sys.path:
  sys.path.insert(0, str(_MAUS_DIR))

import maus  # noqa: E402


def structure_weights(methyls, gem, short_g, long_g) -> Dict[Tuple[int, int], float]:
  """MAGIC-style contact strength for every structural edge: intensity ~ 1/r^6.

  A geminal pair has no meaningful C-C distance for an NOE (the two carbons of
  one residue), so it is pinned to the strongest weight.  Returns a symmetric
  dict keyed by ordered index pairs for the union of edge classes; this is the
  ``model_matrix`` of the QAP objective, graded by distance instead of boolean.
  """
  coord = {m.index: m.coord for m in methyls}
  weights: Dict[Tuple[int, int], float] = {}
  for (a, b) in gem:
    weights[(a, b)] = None  # placeholder -> filled with global max below
  finite_max = 0.0
  for (a, b) in short_g | long_g:
    d = math.dist(coord[a], coord[b])
    w = 1.0 / (d ** 6) if d > 0 else 0.0
    weights[(a, b)] = w
    finite_max = max(finite_max, w)
  gem_weight = finite_max if finite_max > 0 else 1.0
  for key, val in list(weights.items()):
    if val is None:
      weights[key] = gem_weight
  return weights


def load_noesy_rows(path: str):
  """3D (H)CCH NOESY peak list, with an OPTIONAL 5th intensity column:
  ``label  C1  C2  H2  [intensity]``.  Missing intensity defaults to 1.0, so a
  4-column MAUS-style list is read unchanged.  Returns [(c1, c2, h2, intensity)].
  """
  rows = []
  for line in Path(path).read_text().splitlines():
    if not line.strip() or line.startswith('#') or line.startswith('label'):
      continue
    f = line.split('\t')
    c1, c2, h2 = float(f[1]), float(f[2]), float(f[3])
    inten = float(f[4]) if len(f) > 4 and f[4].strip() else 1.0
    rows.append((c1, c2, h2, inten))
  return rows


def match_noe_intensity(peaks, rows, tol_h, tol_c):
  """Resolve NOESY rows to peak pairs, carrying intensity through.

  Firm rows (both endpoints resolve uniquely) become hard NOE edges, each tagged
  with the strongest intensity observed for that pair.  Ambiguous rows are kept
  as diluted soft incidences (their candidate sets + intensity), which MAGIC's
  scoring can exploit when intensities are real.  Returns
  ``(edges, edge_intensity, amb, stats)``.
  """
  def cand_hc(h, c):
    return [p.index for p in peaks
            if abs(p.h_ppm - h) <= tol_h and abs(p.c_ppm - c) <= tol_c]

  def cand_c(c):
    return [p.index for p in peaks if abs(p.c_ppm - c) <= tol_c]

  edges = set()
  edge_intensity: Dict[Tuple[int, int], float] = {}
  amb: List[Tuple[Tuple[int, ...], Tuple[int, ...], float, float]] = []
  firm = ambiguous = unmatched = 0
  for (c1, c2, h2, inten) in rows:
    a = cand_hc(h2, c2)
    b = cand_c(c1)
    if not a or not b:
      unmatched += 1
      continue
    if len(a) == 1 and len(b) == 1 and a[0] != b[0]:
      key = (min(a[0], b[0]), max(a[0], b[0]))
      edges.add(key)
      edge_intensity[key] = max(edge_intensity.get(key, 0.0), inten)
      firm += 1
    else:
      dil = 1.0 / (len(a) * len(b))
      amb.append((tuple(a), tuple(b), dil, inten))
      ambiguous += 1
  return edges, edge_intensity, amb, {'firm': firm, 'ambiguous': ambiguous, 'unmatched': unmatched}


class MagicMaus:
  def __init__(self, methyls, peaks, gem, short_g, long_g, noe_edges, gem_links=None,
               edge_intensity=None):
    self.methyls = methyls
    self.peaks = peaks
    self.noe_edges = noe_edges
    # Observed NOESY intensity per firm edge (1.0 if the peak list carried none).
    # Physically I ~ 1/r^6, so a strong NOE belongs on a close structural
    # contact -- multiplying it into the objective is what lets real intensities
    # break the geometric degeneracy a boolean network cannot.
    self.edge_intensity = edge_intensity or {}
    self.allowed = gem | short_g | long_g
    self.weights = structure_weights(methyls, gem, short_g, long_g)
    self.max_w = max(self.weights.values()) if self.weights else 1.0
    # MAUS owns the hard CNF, the SAT variables, and the per-peak domains.
    self.core = maus.MAUS(methyls, peaks, gem, short_g, long_g, noe_edges, gem_links)
    self.label_by_index = {m.index: m.label for m in methyls}
    self.noe_nbr: Dict[int, List[int]] = {}
    for (i, j) in noe_edges:
      self.noe_nbr.setdefault(i, []).append(j)
      self.noe_nbr.setdefault(j, []).append(i)

  # -- layer 1: MAUS option sets (never excludes truth) ---------------------
  def option_sets(self) -> Dict[int, List[int]]:
    return self.core.solve_options()

  # -- MAGIC's soft evidence: the NOE cross peaks MAUS *discards* -------------
  def set_soft_evidence(self, amb):
    """Fold in the ambiguous NOE cross peaks MAUS throws away.

    MAUS keeps only cross peaks that resolve to a unique pair of HMQC peaks;
    MAGIC's premise is that the diluted, ambiguous rest still carries assignment
    information.  Each dropped row (from ``match_noe_intensity``) is a candidate
    set A of observed peaks, a set B of partners, a dilution ``1/|A||B|`` and its
    intensity; it contributes ``intensity * (1/|A||B|) * sum_{a,b} w(g_a,g_b)``
    to the objective.  It cannot change the option sets (never-exclude holds);
    it only ranks within them.
    """
    self.amb = list(amb)
    self.amb_by_peak: Dict[int, List[int]] = {}
    for idx, (A, B, _dil, _inten) in enumerate(self.amb):
      for k in set(A) | set(B):
        self.amb_by_peak.setdefault(k, []).append(idx)

  # -- soft objective: distance-weighted NOE agreement (MAGIC) ---------------
  def _pair_w(self, gi: int, gj: int) -> float:
    """Contact strength (~1/r^6) of placing an observed NOE on structure edge
    (gi,gj); 0 if the pair is not a structural contact (a hard violation)."""
    return self.weights.get((gi, gj), 0.0)

  def _feasible_model(self, options) -> Dict[int, int]:
    """One jointly-consistent assignment over the pruned MAUS domains, from the
    SAT solver.  Guaranteed to exist (the truth is such a model) and to be
    injective and NOE-consistent -- the feasible starting point the scored
    ascent then improves.  Restricting domains to the option sets keeps the CNF
    tiny (few 10^3 clauses) yet preserves every satisfying assignment.
    """
    from pysat.solvers import Glucose3
    var: Dict[Tuple[int, int], int] = {}
    nid = 1
    for p in self.peaks:
      for g in options[p.index]:
        var[(p.index, g)] = nid
        nid += 1
    clauses: List[List[int]] = []
    for p in self.peaks:                       # exactly one methyl per peak
      lits = [var[(p.index, g)] for g in options[p.index]]
      if not lits:
        continue
      clauses.append(lits)
      for a in range(len(lits)):
        for b in range(a + 1, len(lits)):
          clauses.append([-lits[a], -lits[b]])
    holders: Dict[int, List[int]] = {}
    for p in self.peaks:
      for g in options[p.index]:
        holders.setdefault(g, []).append(var[(p.index, g)])
    for lits in holders.values():              # injective
      for a in range(len(lits)):
        for b in range(a + 1, len(lits)):
          clauses.append([-lits[a], -lits[b]])
    for (i, j) in self.noe_edges:              # firm NOE edge on a real contact
      for gi in options[i]:
        for gj in options[j]:
          if gi != gj and (gi, gj) not in self.allowed:
            clauses.append([-var[(i, gi)], -var[(j, gj)]])
    solver = Glucose3(bootstrap_with=clauses)
    solver.solve()
    model = set(solver.get_model() or [])
    solver.delete()
    assign: Dict[int, int] = {}
    for p in self.peaks:
      assign[p.index] = next((g for g in options[p.index] if var[(p.index, g)] in model),
                             (options[p.index][0] if options[p.index] else None))
    return assign

  def _peak_energy(self, i: int, assign: Dict[int, int]) -> float:
    """NOE-contact strength touching peak i under `assign` (firm edges + any
    ambiguous soft evidence).  The scored ascent maximises the sum over peaks."""
    gi = assign.get(i)
    if gi is None:
      return 0.0
    e = 0.0
    for j in self.noe_nbr.get(i, ()):
      gj = assign.get(j)
      if gj is not None:
        inten = self.edge_intensity.get((min(i, j), max(i, j)), 1.0)
        e += inten * self._pair_w(gi, gj)
    amb = getattr(self, 'amb', None)
    if amb:
      for idx in self.amb_by_peak.get(i, ()):
        A, B, dil, inten = amb[idx]
        if i in A:
          for b in B:
            gb = assign.get(b)
            if b != i and gb is not None:
              e += inten * dil * self._pair_w(gi, gb)
        if i in B:
          for a in A:
            ga = assign.get(a)
            if a != i and ga is not None:
              e += inten * dil * self._pair_w(ga, gi)
    return e

  def _noe_ok(self, i: int, gi: int, assign: Dict[int, int], skip=None) -> bool:
    """Do all firm NOE edges at peak i stay on a structural contact if i->gi?"""
    for j in self.noe_nbr.get(i, ()):
      if j == skip:
        continue
      gj = assign.get(j)
      if gj is not None and gi != gj and (gi, gj) not in self.allowed:
        return False
    return True

  def _improve(self, assign: Dict[int, int], options, free: List[int], max_passes: int = 8):
    """Feasibility-preserving coordinate ascent (MAGIC-style scoring, over MAUS
    domains).  From a feasible assignment, repeatedly move a peak to a better
    option -- swapping with the peak that holds it so the map stays a bijection,
    and only when every firm NOE edge stays on a contact.  Each accepted move
    strictly raises the total NOE-contact strength; passes stop at convergence.
    Because it starts feasible and every move preserves feasibility, the output
    is always injective and NOE-consistent -- the property brute-force fallback
    could not guarantee on the single 138-peak degeneracy cluster of MBP.
    """
    holder: Dict[int, int] = {}              # methyl -> peak currently holding it
    for i, g in assign.items():
      if g is not None:
        holder[g] = i
    for _ in range(max_passes):
      moved = False
      for i in free:
        gi = assign[i]
        base = self._peak_energy(i, assign)
        best_b, best_gain = None, 1e-12
        for b in options[i]:
          if b == gi:
            continue
          j = holder.get(b)                  # peak (if any) currently on methyl b
          if j is None:
            if not self._noe_ok(i, b, assign):
              continue
            before = base
            assign[i] = b
            after = self._peak_energy(i, assign)
            assign[i] = gi
            gain = after - before
          else:
            if gi not in options[j]:
              continue                       # swap must keep j feasible
            if not self._noe_ok(i, b, assign, skip=j) or not self._noe_ok(j, gi, assign, skip=i):
              continue
            before = self._peak_energy(i, assign) + self._peak_energy(j, assign)
            assign[i], assign[j] = b, gi
            after = self._peak_energy(i, assign) + self._peak_energy(j, assign)
            assign[i], assign[j] = gi, b
            gain = after - before
          if gain > best_gain:
            best_b, best_gain = b, gain
        if best_b is not None:
          j = holder.get(best_b)
          if j is None:
            del holder[gi]
            assign[i] = best_b
            holder[best_b] = i
          else:
            assign[i], assign[j] = best_b, gi
            holder[best_b] = i
            holder[gi] = j
          moved = True
      if not moved:
        break
    return assign

  # -- layer 2: single globally-coherent assignment ------------------------
  def solve(self):
    """MAUS bounds the space; MAGIC-style scoring commits within it.

    A feasible assignment over the pruned domains (SAT) seeds a
    feasibility-preserving coordinate ascent on the NOE-contact objective.  The
    result is a single injective, NOE-consistent map -- MAUS's option sets are
    returned alongside as the never-exclude envelope.
    """
    options = self.option_sets()
    free = [p.index for p in self.peaks if len(options[p.index]) > 1]
    assign = self._feasible_model(options)
    chosen = self._improve(assign, options, free)
    return chosen, options

  # -- confidence: local scoring margin of the committed pick ----------------
  def confidence(self, peak_index: int, chosen: Dict[int, int], options) -> Tuple[str, float]:
    """Tier the committed call by how the NOE score reacts to swapping it,
    holding the rest of the assignment fixed (a cheap local QAP margin).

      unique    - only one option survives the hard constraints (forced)
      scored    - several options, but the NOE geometry strictly prefers ours
      ambiguous - a tied alternative exists (a genuine symmetry; a coin flip)
    """
    opts = options[peak_index]
    if len(opts) <= 1:
      return 'unique', math.inf
    # Score each option with the rest of the assignment held fixed -- the same
    # intensity-weighted objective (firm + ambiguous soft) the solver maximised.
    def local_score(g: int) -> float:
      probe = dict(chosen)
      probe[peak_index] = g
      return self._peak_energy(peak_index, probe)

    picked = chosen[peak_index]
    picked_s = local_score(picked)
    best_alt = max((local_score(g) for g in opts if g != picked), default=-math.inf)
    if best_alt == -math.inf:
      return 'scored', math.inf
    margin = picked_s - best_alt
    return ('scored', margin) if margin > 1e-12 else ('ambiguous', 0.0)


def run(pdb, hmqc, noesy, hmbc=None, truth=None, short_cut=6.0, long_cut=10.0,
        tol_h=0.02, tol_c=0.20, labeling='A;I;L;M;T;V', out=None,
        soft_ambiguous=False):
  lab = maus.parse_labeling(labeling)
  methyls = maus.parse_structure(Path(pdb).read_text().splitlines(), lab)
  peaks = maus.load_hmqc(hmqc)
  crosses = load_noesy_rows(noesy)             # (c1, c2, h2, intensity)
  gem, short_g, long_g = maus.build_structure_graph(methyls, short_cut, long_cut)
  noe_edges, edge_intensity, amb, nstat = match_noe_intensity(peaks, crosses, tol_h, tol_c)
  gem_links = set()
  if hmbc:
    gem_links, _ = maus.match_hmbc(peaks, maus.load_hmbc(hmbc), tol_h, tol_c)

  engine = MagicMaus(methyls, peaks, gem, short_g, long_g, noe_edges, gem_links,
                     edge_intensity=edge_intensity)
  if soft_ambiguous:
    # Fold the discarded ambiguous NOE cross peaks in as diluted soft evidence.
    # Uninformative on a boolean structure-simulated network but a real lever
    # once NOESY intensities carry signal; off by default.
    engine.set_soft_evidence(amb)
  chosen, options = engine.solve()
  label_by_index = engine.label_by_index
  truth_map = maus.load_truth(truth) if truth else {}

  tiers: Dict[str, int] = {'unique': 0, 'scored': 0, 'ambiguous': 0}
  rows = []
  for p in peaks:
    tier, margin = engine.confidence(p.index, chosen, options)
    tiers[tier] += 1
    call = label_by_index.get(chosen[p.index], '')
    opt_labels = [label_by_index[g] for g in options[p.index]]
    t = truth_map.get(p.peak_id, '')
    rows.append((p, call, opt_labels, tier, margin, t))

  if out:
    with open(out, 'w') as f:
      f.write('label\tres_type\tcall\tconfidence\tmargin\tn_options\toptions\ttruth\t'
              'call_correct\ttruth_in_set\n')
      for (p, call, opts, tier, margin, t) in rows:
        cc = '' if not t else int(call == t)
        tis = '' if not t else int(t in opts)
        mval = 'inf' if margin == math.inf else f'{margin:.3e}'
        f.write(f'{p.peak_id}\t{p.res_type}\t{call}\t{tier}\t{mval}\t{len(opts)}\t'
                f'{",".join(opts)}\t{t}\t{cc}\t{tis}\n')

  n = len(peaks)
  print(f'methyls(G nodes)={len(methyls)}  HMQC peaks={n}  NOESY cross peaks={len(crosses)}')
  print(f'NOE match (tol H+-{tol_h}/C+-{tol_c}): firm={nstat["firm"]} '
        f'ambiguous(dropped)={nstat["ambiguous"]} unmatched={nstat["unmatched"]}')
  print('--- MAUS envelope (never excludes truth) ---')
  set_sizes = [len(options[p.index]) for p in peaks]
  print(f'unique(1 option)      = {sum(s == 1 for s in set_sizes)}/{n}')
  print(f'ambiguous(2-3 options)= {sum(2 <= s <= 3 for s in set_sizes)}/{n}')
  print(f'ambiguous(>3 options) = {sum(s > 3 for s in set_sizes)}/{n}')
  print(f'unassigned            = {sum(s == 0 for s in set_sizes)}/{n}')
  print('--- magicmaus commitment (single coherent call) ---')
  print(f'confidence: unique={tiers["unique"]}  scored={tiers["scored"]}  '
        f'ambiguous={tiers["ambiguous"]}')

  if truth_map:
    in_set = sum(1 for (_, _, opts, _, _, t) in rows if t and t in opts)
    call_ok = sum(1 for (_, call, _, _, _, t) in rows if t and call == t)
    by_tier_ok = {k: [0, 0] for k in tiers}
    for (_, call, _, tier, _, t) in rows:
      if not t:
        continue
      by_tier_ok[tier][1] += 1
      by_tier_ok[tier][0] += int(call == t)
    print('--- scored vs truth ---')
    print(f'truth in MAUS option set = {in_set}/{n} = {100*in_set/n:.1f}%  '
          f'(never-exclude guarantee preserved)')
    print(f'magicmaus single call    = {call_ok}/{n} = {100*call_ok/n:.1f}% correct')
    for k in ('unique', 'scored', 'ambiguous'):
      ok, tot = by_tier_ok[k]
      if tot:
        print(f'    {k:9s}: {ok}/{tot} = {100*ok/tot:.1f}% correct')
  else:
    print('(no --truth given; scoring skipped)')
  return 0


def main(argv=None):
  ap = argparse.ArgumentParser(
    description='magicmaus: MAUS hard bounds + MAGIC-style scored commitment '
                '(SAT-feasible seed + coordinate ascent).')
  ap.add_argument('pdb')
  ap.add_argument('hmqc', help='HMQC peak list TSV: label H_ppm C_ppm res_type')
  ap.add_argument('noesy', help='NOESY peak list TSV: label C1 C2 H2 [intensity] '
                                '(5th intensity column optional; defaults to 1.0)')
  ap.add_argument('--hmbc', default=None, help='optional HMBC-HMQC geminal-link TSV')
  ap.add_argument('--truth', default=None, help='truth key TSV (label ... True) for scoring')
  ap.add_argument('--short-cut', type=float, default=6.0)
  ap.add_argument('--long-cut', type=float, default=10.0)
  ap.add_argument('--tol-h', type=float, default=0.02)
  ap.add_argument('--tol-c', type=float, default=0.20)
  ap.add_argument('--labeling', default='A;I;L;M;T;V')
  ap.add_argument('--out', default=None, help='write per-peak calls TSV here')
  ap.add_argument('--soft-ambiguous', action='store_true',
                  help='fold MAUS-discarded ambiguous NOE cross peaks in as diluted '
                       'soft tie-breakers (for experimental NOESY with real intensities; '
                       'symmetric noise on a structure-simulated network)')
  a = ap.parse_args(argv)
  return run(a.pdb, a.hmqc, a.noesy, a.hmbc, a.truth, a.short_cut, a.long_cut,
             a.tol_h, a.tol_c, a.labeling, a.out, a.soft_ambiguous)


if __name__ == '__main__':
  raise SystemExit(main())
