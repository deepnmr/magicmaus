"""Engine: MAUS hard-constraint SAT (per-peak option sets, never excludes the
truth) fused with MAGIC distance-weighted NOE scoring (one coherent call within
the residual degeneracy).

Layer 1 (MAUS) — a CNF over x(peak, methyl) variables enforces: exactly one
methyl per peak (of the peak's type), injectivity, and every firm NOE / HMBC
edge lands on a structural contact within long_cut.  For each peak the set of
methyls appearing in >=1 satisfying model is enumerated.  The truth is retained
whenever every firm edge is correct; under measurement noise a mis-resolved
"firm" edge can prune it, which the CLI flags rather than hides.

Layer 2 (MAGIC) — over MAUS's pruned option sets, a feasibility-preserving
search maximises the distance-weighted NOE objective (intensity ~ 1/r^6) and
commits a single injective, NOE-consistent map.  Because the search runs only
inside the truth-containing MAUS space it is both fast and able to find the
truth, unlike MAGIC over the full candidate space.
"""
from __future__ import annotations

import math
import random
from collections import deque
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pysat.card import CardEnc, EncType
from pysat.solvers import Glucose3

from . import peaks as pk
from . import structure as st


# ---------------------------------------------------------------------------
# Layer 1: MAUS hard SAT option sets
# ---------------------------------------------------------------------------
class _Maus:
    def __init__(self, methyls, peaks, gem, short_g, long_g, noe_edges, gem_links=None):
        self.methyls = methyls
        self.peaks = peaks
        self.gem, self.short_g, self.long_g = gem, short_g, long_g
        self.noe_edges = noe_edges
        self.gem_links = gem_links or set()
        self.sites_by_type: Dict[str, List[int]] = {}
        for m in methyls:
            self.sites_by_type.setdefault(m.res_type, []).append(m.index)
        idx_by_label = {m.label: m.index for m in methyls}

        # initial domains: methyls of the peak's type (or a tentative pin)
        self.domain: Dict[int, List[int]] = {}
        self.n_tentative = 0
        for p in peaks:
            dom = self.sites_by_type.get(p.res_type, [])
            if p.tentative and idx_by_label.get(p.tentative) in dom:
                dom = [idx_by_label[p.tentative]]
                self.n_tentative += 1
            self.domain[p.index] = list(dom)

        self._prune_domains()

        self.var: Dict[Tuple[int, int], int] = {}
        nid = 1
        for p in peaks:
            for g in self.domain[p.index]:
                self.var[(p.index, g)] = nid
                nid += 1
        self.n_vars = nid - 1

    def _prune_domains(self) -> None:
        """Arc-consistency prune on the firm NOE / HMBC constraints before the
        expensive per-candidate SAT: a methyl with no allowed structural partner
        in a linked peak's domain cannot be in any satisfying assignment, so
        dropping it is exact and preserves the never-exclude guarantee."""
        allowed = self.gem | self.short_g | self.long_g
        nbr: Dict[int, List[Tuple[int, frozenset]]] = {}
        for (i, j) in self.noe_edges:
            nbr.setdefault(i, []).append((j, allowed))
            nbr.setdefault(j, []).append((i, allowed))
        for (i, j) in self.gem_links:
            nbr.setdefault(i, []).append((j, self.gem))
            nbr.setdefault(j, []).append((i, self.gem))
        dom = {pi: set(d) for pi, d in self.domain.items()}
        queue = deque((x, y, al) for x, arcs in nbr.items() for (y, al) in arcs)
        while queue:
            x, y, al = queue.popleft()
            dy = dom[y]
            dropped = [gx for gx in dom[x]
                       if not any(gy != gx and (gx, gy) in al for gy in dy)]
            if dropped:
                dom[x].difference_update(dropped)
                for (z, al2) in nbr.get(x, ()):
                    if z != y:
                        queue.append((z, x, al2))
        for pi in self.domain:
            self.domain[pi] = [g for g in self.domain[pi] if g in dom[pi]]

    def _amo(self, clauses, lits):
        """At-most-one: pairwise for small domains, sequential-counter for large
        ones (O(d) clauses + aux vars, so a 130-methyl Leu domain is cheap)."""
        if len(lits) <= 1:
            return
        if len(lits) <= 8:
            for a, b in combinations(lits, 2):
                clauses.append([-a, -b])
            return
        enc = CardEnc.atmost(lits=lits, bound=1, top_id=self._top, encoding=EncType.seqcounter)
        clauses.extend(enc.clauses)
        if enc.nv > self._top:
            self._top = enc.nv

    def _support(self, clauses, i, j, allowed):
        """i->gi implies j maps to an allowed partner (both directions)."""
        di, dj = self.domain[i], self.domain[j]
        for gi in di:
            clauses.append([-self.var[(i, gi)]] +
                           [self.var[(j, gj)] for gj in dj if gj != gi and (gi, gj) in allowed])
        for gj in dj:
            clauses.append([-self.var[(j, gj)]] +
                           [self.var[(i, gi)] for gi in di if gi != gj and (gi, gj) in allowed])

    def _base_clauses(self):
        clauses: List[List[int]] = []
        self._top = self.n_vars
        for p in self.peaks:                          # exactly one methyl per peak
            lits = [self.var[(p.index, g)] for g in self.domain[p.index]]
            if not lits:
                continue
            clauses.append(lits)
            self._amo(clauses, lits)
        peaks_of: Dict[int, List[int]] = {}           # injective: each methyl <=1 peak
        for p in self.peaks:
            for g in self.domain[p.index]:
                peaks_of.setdefault(g, []).append(self.var[(p.index, g)])
        for lits in peaks_of.values():
            self._amo(clauses, lits)
        allowed = self.gem | self.short_g | self.long_g
        for (i, j) in self.noe_edges:                 # firm NOE on a structural contact
            self._support(clauses, i, j, allowed)
        for (i, j) in self.gem_links:                 # HMBC geminal link -> geminal edge
            self._support(clauses, i, j, self.gem)
        return clauses

    def solve_options(self) -> Dict[int, List[int]]:
        base = self._base_clauses()
        solver = Glucose3(bootstrap_with=base)
        owner = {v: k for k, v in self.var.items()}
        confirmed = set()
        for p in self.peaks:
            for g in self.domain[p.index]:
                if (p.index, g) in confirmed:
                    continue
                lit = self.var[(p.index, g)]
                # Cheap unit-propagation filter: if asserting this candidate
                # conflicts by propagation alone, it is invalid — skip the ~1000x
                # more expensive full solve.  This is what makes wide-domain
                # targets (MSG: 138-wide Leu/Val) tractable; most of the many
                # invalid candidates die here.
                try:
                    ok, _ = solver.propagate(assumptions=[lit])
                    if not ok:
                        continue
                except Exception:
                    pass
                if solver.solve(assumptions=[lit]):
                    for m in solver.get_model() or ():
                        if m > 0 and m in owner:
                            confirmed.add(owner[m])
        options = {p.index: [g for g in self.domain[p.index] if (p.index, g) in confirmed]
                   for p in self.peaks}
        solver.delete()
        return options


def structure_weights(methyls, gem, short_g, long_g) -> Dict[Tuple[int, int], float]:
    """MAGIC contact strength for every structural edge: intensity ~ 1/r^6, min
    over chain images.  A geminal pair (no meaningful C-C NOE distance) is pinned
    to the strongest finite weight."""
    by_idx = {m.index: m for m in methyls}
    weights: Dict[Tuple[int, int], float] = {}
    finite_max = 0.0
    for (a, b) in short_g | long_g:
        d = st.min_dist(by_idx[a], by_idx[b])
        w = 1.0 / (d ** 6) if d > 0 else 0.0
        weights[(a, b)] = w
        finite_max = max(finite_max, w)
    gem_weight = finite_max if finite_max > 0 else 1.0
    for (a, b) in gem:
        weights[(a, b)] = gem_weight
    return weights


# ---------------------------------------------------------------------------
# Layer 2: MAGIC-style scored commitment over the MAUS option sets
# ---------------------------------------------------------------------------
class MagicMaus:
    def __init__(self, methyls, peaks, gem, short_g, long_g, noe_edges,
                 gem_links=None, edge_intensity=None, enforce_contacts=True):
        self.methyls = methyls
        self.peaks = peaks
        self.noe_edges = noe_edges
        self.edge_intensity = edge_intensity or {}
        # Whether the SCORED commitment search hard-requires every firm NOE edge
        # to sit on a structural contact.  On for magicmaus: the hard NOE
        # constraint is the dominant disambiguating signal.  A pure-MAGIC baseline
        # sets it off (soft 1/r^6 reward only) to show that MAGIC without MAUS's
        # hard constraints scatters.  Geminal HMBC links stay hard regardless
        # (residue-local, never trap the search, and matter residue-wise).
        self.enforce_contacts = enforce_contacts
        self.gem = gem
        self.allowed = gem | short_g | long_g
        self.gem_links = gem_links or set()   # HMBC geminal same-residue links
        self.weights = structure_weights(methyls, gem, short_g, long_g)
        self.core = _Maus(methyls, peaks, gem, short_g, long_g, noe_edges, gem_links)
        self.n_tentative = self.core.n_tentative
        self.label_by_index = {m.index: m.label for m in methyls}
        self.noe_nbr: Dict[int, List[int]] = {}
        for (i, j) in noe_edges:
            self.noe_nbr.setdefault(i, []).append(j)
            self.noe_nbr.setdefault(j, []).append(i)
        # HMBC geminal links are hard constraints on the *committed* map too, not
        # just on the option sets: the two linked peaks must land on the two
        # geminal methyls of one residue.  Kept separate from noe_nbr because
        # their allowed set is geminal-only and they carry no scoring weight.
        self.gem_nbr: Dict[int, List[int]] = {}
        for (i, j) in self.gem_links:
            self.gem_nbr.setdefault(i, []).append(j)
            self.gem_nbr.setdefault(j, []).append(i)
        self.amb: List[Tuple[tuple, tuple, float, float]] = []
        self.amb_by_peak: Dict[int, List[int]] = {}

    def option_sets(self) -> Dict[int, List[int]]:
        return self.core.solve_options()

    def set_soft_evidence(self, ambiguous):
        """Fold the ambiguous NOE rows MAUS discards in as diluted soft evidence
        (each contributes intensity * dilution * sum of pair weights).  It cannot
        change the option sets — it only ranks within them."""
        self.amb = list(ambiguous)
        self.amb_by_peak = {}
        for idx, (A, B, _dil, _inten) in enumerate(self.amb):
            for k in set(A) | set(B):
                self.amb_by_peak.setdefault(k, []).append(idx)

    def _pair_w(self, gi, gj) -> float:
        return self.weights.get((gi, gj), 0.0)

    def _seed_cnf(self, options):
        """CNF whose models are the feasible assignments over `options`: exactly
        one methyl per peak, injective, firm NOE edges on contacts (when
        enforce_contacts), and HMBC geminal links on geminal edges.  Returns
        (clauses, var, nv)."""
        var: Dict[Tuple[int, int], int] = {}
        nid = 1
        for p in self.peaks:
            for g in options[p.index]:
                var[(p.index, g)] = nid
                nid += 1
        clauses: List[List[int]] = []
        for p in self.peaks:
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
        for lits in holders.values():
            for a in range(len(lits)):
                for b in range(a + 1, len(lits)):
                    clauses.append([-lits[a], -lits[b]])
        if self.enforce_contacts:
            for (i, j) in self.noe_edges:
                for gi in options[i]:
                    for gj in options[j]:
                        if gi != gj and (gi, gj) not in self.allowed:
                            clauses.append([-var[(i, gi)], -var[(j, gj)]])
        for (i, j) in self.gem_links:            # HMBC geminal link -> geminal edge only
            for gi in options[i]:
                for gj in options[j]:
                    if (gi, gj) not in self.gem:
                        clauses.append([-var[(i, gi)], -var[(j, gj)]])
        return clauses, var, nid - 1

    def _model_to_assign(self, options, var, model):
        return {p.index: next((g for g in options[p.index] if var[(p.index, g)] in model),
                              options[p.index][0] if options[p.index] else None)
                for p in self.peaks}

    def _feasible_model(self, options) -> Dict[int, int]:
        """One jointly-consistent assignment over the option sets, from SAT
        (guaranteed to exist — the truth is one)."""
        clauses, var, _ = self._seed_cnf(options)
        solver = Glucose3(bootstrap_with=clauses)
        solver.solve()
        model = set(solver.get_model() or [])
        solver.delete()
        return self._model_to_assign(options, var, model)

    def _diverse_seeds(self, options, restarts, seed):
        """Yield up to `restarts` DIVERSE feasible seeds by re-solving the seed
        CNF with randomised variable phases.  Greedy ascent from each lands in a
        different basin, so the best over many starts reaches the objective's
        optimum inside the pruned space — where a single seed's local search
        traps well below it."""
        clauses, var, nv = self._seed_cnf(options)
        rng = random.Random(seed)
        for r in range(restarts):
            s = Glucose3(bootstrap_with=clauses)
            if r > 0:                              # r==0 keeps the solver's default model
                try:
                    s.set_phases([(v if rng.random() < 0.5 else -v) for v in range(1, nv + 1)])
                except Exception:
                    pass
            sat = s.solve()
            model = set(s.get_model() or [])
            s.delete()
            if sat:
                yield self._model_to_assign(options, var, model)

    def _peak_energy(self, i, assign) -> float:
        gi = assign.get(i)
        if gi is None:
            return 0.0
        e = 0.0
        for j in self.noe_nbr.get(i, ()):
            gj = assign.get(j)
            if gj is not None:
                inten = self.edge_intensity.get((min(i, j), max(i, j)), 1.0)
                e += inten * self._pair_w(gi, gj)
        if self.amb:
            for idx in self.amb_by_peak.get(i, ()):
                A, B, dil, inten = self.amb[idx]
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

    def _noe_ok(self, i, gi, assign, skip=None) -> bool:
        """HMBC geminal links (always) and, when enforce_contacts is on, firm NOE
        edges at peak i stay satisfied if i->gi."""
        if self.enforce_contacts:
            for j in self.noe_nbr.get(i, ()):
                if j == skip:
                    continue
                gj = assign.get(j)
                if gj is not None and gi != gj and (gi, gj) not in self.allowed:
                    return False
        for j in self.gem_nbr.get(i, ()):
            if j == skip:
                continue
            gj = assign.get(j)
            if gj is not None and (gi, gj) not in self.gem:
                return False
        return True

    def _improve(self, assign, options, free, max_passes=8):
        """Feasibility-preserving coordinate ascent: move a peak to a better
        option, swapping with the peak that holds it so the map stays a bijection
        and every firm NOE edge stays on a contact."""
        holder: Dict[int, int] = {g: i for i, g in assign.items() if g is not None}
        for _ in range(max_passes):
            moved = False
            for i in free:
                gi = assign[i]
                best_b, best_gain = None, 1e-12
                for b in options[i]:
                    if b == gi:
                        continue
                    j = holder.get(b)
                    if j is None:
                        if not self._noe_ok(i, b, assign):
                            continue
                        before = self._peak_energy(i, assign)
                        assign[i] = b
                        after = self._peak_energy(i, assign)
                        assign[i] = gi
                        gain = after - before
                    else:
                        if gi not in options[j]:
                            continue
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

    def _total_obj(self, assign) -> float:
        return 0.5 * sum(self._peak_energy(p.index, assign) for p in self.peaks)

    def solve(self, restarts=4, iters=None, seed=0):
        options = self.option_sets()
        return self.optimize(options, restarts=restarts, iters=iters, seed=seed), options

    def optimize(self, options, restarts=24, iters=None, seed=0):
        """Commit a single map maximising the NOE objective over `options`.

        The objective is a QAP whose truth is a strong local optimum, but a single
        feasible seed's greedy ascent traps well below it (the pruned Leu/Val
        option graphs are rugged).  So run a DIVERSE-SEED multistart: many
        independent feasible seeds (randomised SAT phases), each greedy-ascended,
        keeping the highest-objective result.  This reaches the objective's
        optimum inside MAUS's small pruned space — which is exactly what MAGIC
        cannot do over its un-pruned full space, where the same multistart lands
        in one of astronomically many near-flat basins.  `iters` is accepted for
        API compatibility and ignored."""
        free = [p.index for p in self.peaks if len(options[p.index]) > 1]
        # only NOE/geminal-connected peaks can change the objective; a peak with
        # neither is objective-neutral, so leave it at the seed (its residue-wise
        # call is a coin flip regardless).  This is what keeps the ascent fast on
        # MAGIC's un-pruned full domains.
        movable = [i for i in free if self.noe_nbr.get(i) or self.gem_nbr.get(i)]
        best = self._improve(self._feasible_model(options), options, movable)
        if not movable:
            return best
        best_obj = self._total_obj(best)
        n_restarts = max(restarts, min(32, 2 * len(movable)))   # scale to problem size
        for cand in self._diverse_seeds(options, n_restarts, seed):
            c = self._improve(cand, options, movable)
            o = self._total_obj(c)
            if o > best_obj + 1e-18:
                best, best_obj = c, o
        return best

    def confidence(self, peak_index, chosen, options) -> Tuple[str, float]:
        """Tier the committed call by its local NOE-score margin over the best
        alternative: unique (forced) / scored (geometry prefers it) / ambiguous
        (a tied alternative — a genuine symmetry)."""
        opts = options[peak_index]
        if len(opts) == 0:
            return 'unassigned', -math.inf   # no methyl survives the hard constraints
        if len(opts) == 1:
            return 'unique', math.inf

        def local_score(g):
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


# ---------------------------------------------------------------------------
# Convenience driver
# ---------------------------------------------------------------------------
@dataclass
class Result:
    engine: MagicMaus
    chosen: Dict[int, int]
    options: Dict[int, List[int]]
    stats: dict


def assign(structure, hmqc, noesy, hmbc=None, *, short_cut=6.0, long_cut=10.0,
           tol_h=0.02, tol_c=0.10, labeling='A;I;L;M;T;V',
           soft_ambiguous=False, restarts=4, seed=0) -> Result:
    """Run magicmaus end to end on file paths and return a Result.

    structure : mmCIF or PDB path (homo-oligomer -> multimer contacts)
    hmqc      : HMQC TSV        (label H_ppm C_ppm res_type)
    noesy     : 3D NOESY TSV    (label C2 C1 H1 [intensity])
    hmbc      : 3D HMBC-HMQC TSV (label C2 C1 H1), geminal links; optional
    """
    lab = pk.parse_labeling(labeling)
    methyls = st.parse_structure(Path(structure).read_text().splitlines(), lab)
    peaks = pk.load_hmqc(hmqc)
    gem, short_g, long_g = st.build_structure_graph(methyls, short_cut, long_cut)

    noe_edges, edge_intensity, amb, nstat = pk.match_noesy(
        peaks, pk.load_triple(noesy), tol_h, tol_c, keep_amb=soft_ambiguous)
    gem_links, hstat = set(), None
    if hmbc:
        gem_links, hstat = pk.match_hmbc(peaks, pk.load_triple(hmbc), tol_h, tol_c)

    eng = MagicMaus(methyls, peaks, gem, short_g, long_g, noe_edges, gem_links,
                    edge_intensity=edge_intensity)
    if soft_ambiguous:
        eng.set_soft_evidence(amb)
    chosen, options = eng.solve(restarts=restarts, seed=seed)
    return Result(eng, chosen, options,
                  {'noe': nstat, 'hmbc': hstat,
                   'G': {'gem': len(gem) // 2, 'short': len(short_g) // 2, 'long': len(long_g) // 2},
                   'n_methyls': len(methyls), 'n_peaks': len(peaks), 'n_noe_edges': len(noe_edges)})
