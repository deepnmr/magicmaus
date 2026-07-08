"""Head-to-head of the three ideas on the SAME noisy peak lists:

  * MAUS   — hard SAT constraints only, no scoring.  Commits a peak only when the
             constraints force a unique methyl; on an ambiguous set it cannot
             rank, so a forced single call must tiebreak arbitrarily.  Its real
             strength is the option-set envelope (truth-in-set), not a single call.
  * MAGIC  — the 1/r^6 distance-weighted NOE score, optimised over the FULL
             type-matched candidate space (no SAT bounding).  Decisive but the
             landscape over that large space is near-flat, so it commits
             near-optimal-but-wrong.
  * magicmaus — the SAME MAGIC score, optimised only inside the MAUS-pruned
             option sets.  The bounding is what lets the scoring actually find
             the truth.

Metric: RESIDUE-WISE accuracy = fraction of ALL peaks committed to the correct
*residue* (res_type + res_num), ignoring the geminal CD1/CD2 (Leu) or CG1/CG2
(Val) swap, which is a near-symmetric coin flip and not the assignment problem's
hard part.  magicmaus should beat both baselines; MAUS also reports its honest
envelope (% of peaks whose true residue is retained in the option set).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from . import peaks as pk
from . import structure as st
from .engine import MagicMaus


@dataclass
class Compare:
    name: str
    n_peaks: int
    multimer: bool
    envelope: float          # % truth retained in MAUS option set
    maus_unique_cov: float   # % peaks MAUS pins uniquely
    maus_unique_acc: float   # % of those unique pins that are correct
    maus_single: float       # % all peaks correct, MAUS forced to commit (arbitrary tiebreak)
    magic_single: float      # % all peaks correct, MAGIC over full space
    magicmaus_single: float  # % all peaks correct, magicmaus (scored within MAUS bounds)


def three_way(structure, hmqc, noesy, hmbc, truth, *, short_cut=6.0, long_cut=10.0,
              tol_h=0.02, tol_c=0.10, labeling='A;I;L;M;T;V', restarts=4, seed=0,
              name='', enumerate_options=True) -> Compare:
    """enumerate_options=True runs MAUS's exact per-candidate option enumeration.
    Set False to use the arc-consistency-pruned domains instead (a looser but
    still never-exclude envelope) — needed for large targets like MSG (257 peaks,
    138-wide Leu/Val domains) where exact enumeration is intractable, and itself
    a point: magicmaus stays tractable there."""
    lab = pk.parse_labeling(labeling)
    methyls = st.parse_structure(Path(structure).read_text().splitlines(), lab)
    multimer = any(len(m.images) > 1 for m in methyls)
    peaks = pk.load_hmqc(hmqc)
    gem, short_g, long_g = st.build_structure_graph(methyls, short_cut, long_cut)
    noe_edges, edge_intensity, _amb, _ = pk.match_noesy(peaks, pk.load_triple(noesy), tol_h, tol_c)
    gem_links = set()
    if hmbc:
        gem_links, _ = pk.match_hmbc(peaks, pk.load_triple(hmbc), tol_h, tol_c)

    eng = MagicMaus(methyls, peaks, gem, short_g, long_g, noe_edges, gem_links,
                    edge_intensity=edge_intensity)
    label = eng.label_by_index
    truth_map = pk.load_truth(truth)
    n = len(peaks)

    # residue key (res_type, res_num) per methyl index, and per truth peak.
    res_of = {m.index: (m.res_type, m.res_num) for m in methyls}
    by_label = {m.label: m.index for m in methyls}
    truth_res = {p.index: res_of.get(by_label.get(truth_map.get(p.peak_id, '')))
                 for p in peaks}

    def correct(call: Dict[int, Optional[int]]) -> int:
        """Residue-wise: called methyl's residue == truth residue (atom ignored)."""
        return sum(1 for p in peaks
                   if call.get(p.index) is not None
                   and truth_res[p.index] is not None
                   and res_of[call[p.index]] == truth_res[p.index])

    # ---- MAUS layer: option sets (exact enumeration, or the pruned domains) ----
    options = eng.option_sets() if enumerate_options \
        else {p.index: list(eng.core.domain[p.index]) for p in peaks}
    in_set = sum(1 for p in peaks
                 if truth_res[p.index] is not None
                 and truth_res[p.index] in {res_of[g] for g in options[p.index]})
    uniq = [p for p in peaks if len({res_of[g] for g in options[p.index]}) == 1]
    uniq_ok = sum(1 for p in uniq
                  if truth_res[p.index] == res_of[options[p.index][0]])
    # MAUS forced single call: unique -> its methyl; ambiguous -> first (arbitrary,
    # since MAUS has no way to rank within the set).
    maus_call = {p.index: (options[p.index][0] if options[p.index] else None) for p in peaks}
    maus_ok = correct(maus_call)

    # ---- MAGIC alone: soft 1/r^6 scoring over the FULL type-matched space,
    #      WITHOUT MAUS's hard NOE-on-contact constraint (the faithful MAGIC:
    #      a global soft optimiser with no hard combinatorial bounding). ----
    full = {p.index: list(eng.core.sites_by_type.get(p.res_type, [])) for p in peaks}
    eng.enforce_contacts = False
    magic_call = eng.optimize(full, restarts=restarts, seed=seed)
    magic_ok = correct(magic_call)

    # ---- magicmaus: MAGIC scoring WITH MAUS's hard constraints — the firm NOE
    #      edges are enforced and the search runs only inside the pruned options. ----
    eng.enforce_contacts = True
    mm_call = eng.optimize(options, restarts=restarts, seed=seed)
    mm_ok = correct(mm_call)

    pct = lambda a, b: 100.0 * a / b if b else 0.0
    return Compare(
        name=name or Path(structure).stem, n_peaks=n, multimer=multimer,
        envelope=pct(in_set, n),
        maus_unique_cov=pct(len(uniq), n), maus_unique_acc=pct(uniq_ok, len(uniq)),
        maus_single=pct(maus_ok, n), magic_single=pct(magic_ok, n),
        magicmaus_single=pct(mm_ok, n))
