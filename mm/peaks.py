"""Peak-list layer: parse the experimental inputs and match 3D cross peaks back
to pairs of HMQC peaks by frequency.

The 3D lists use the column order  label  C2  C1  H1  [intensity] :
  * the *detected* methyl resonates at (C1, H1) — its own carbon and proton;
  * the *partner* methyl contributes carbon C2 only (no proton in the plane),
    which is the irreducible source of ambiguity (many carbons overlap).

A cross peak becomes a firm edge only when BOTH endpoints resolve to a unique
HMQC peak; otherwise it is ambiguous.  This is conservative: a valid assignment
is never excluded by dropping an ambiguous edge from the hard layer.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .structure import THREE_TO_ONE, METHYL_ATOMS  # re-exported for convenience


def parse_labeling(spec: str) -> Dict[str, List[Tuple[str, Optional[str]]]]:
    """'A;I;L;M;T;V' -> {one_letter: METHYL_ATOMS[one_letter]}."""
    labeling: Dict[str, List[Tuple[str, Optional[str]]]] = {}
    for chunk in spec.split(';'):
        one = chunk.split(',')[0].strip()
        if one:
            labeling[one] = METHYL_ATOMS.get(one, [])
    return labeling


@dataclass(frozen=True)
class Peak:
    index: int
    peak_id: str
    res_type: str
    h_ppm: float
    c_ppm: float
    tentative: str = ''   # normalized structure label (e.g. L45CD2) or ''


def parse_res_field(cell: str) -> Tuple[str, str]:
    """The HMQC res_type cell is a bare one-letter type ('L') or a tentative
    assignment ('L45D2' / 'L45CD2').  Returns (one_letter_type, tentative_label);
    tentative is '' for a bare type."""
    cell = cell.strip()
    m = re.match(r'^([A-Z])(\d+)([A-Za-z]\w*)$', cell)
    if not m:
        return cell, ''
    one, num, atom = m.group(1), m.group(2), m.group(3).upper()
    if not atom.startswith('C'):
        atom = 'C' + atom
    return one, f'{one}{num}{atom}'


def _rows(path: str):
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#') or line.startswith('label'):
            continue
        yield line.split('\t')


def load_hmqc(path: str) -> List[Peak]:
    """HMQC peak list: label  H_ppm  C_ppm  res_type."""
    peaks: List[Peak] = []
    for f in _rows(path):
        label, h, c, rtype = f[0], f[1], f[2], f[3]
        one, tent = parse_res_field(rtype)
        peaks.append(Peak(index=len(peaks), peak_id=label, res_type=one,
                          h_ppm=float(h), c_ppm=float(c), tentative=tent))
    return peaks


def load_triple(path: str) -> List[Tuple[float, float, float, float]]:
    """3D list  label  C2  C1  H1  [intensity]  -> [(c2, c1, h1, intensity)].
    Missing intensity defaults to 1.0 (both NOESY and HMBC share this layout)."""
    out = []
    for f in _rows(path):
        c2, c1, h1 = float(f[1]), float(f[2]), float(f[3])
        inten = float(f[4]) if len(f) > 4 and f[4].strip() else 1.0
        out.append((c2, c1, h1, inten))
    return out


def load_truth(path: str) -> Dict[str, str]:
    """Truth key: label  H_ppm  C_ppm  res_type  True  -> {label: True}."""
    truth: Dict[str, str] = {}
    for f in _rows(path):
        if len(f) >= 5:
            truth[f[0]] = f[4].strip()
    return truth


def _match(peaks: List[Peak], rows, tol_h: float, tol_c: float, keep_amb: bool):
    """Resolve each 3D row to a pair of HMQC peaks.  Detected methyl by (H1,C1);
    partner by C2 (carbon only).  Returns (firm_edges, edge_intensity, ambiguous,
    stats).  Firm edges are undirected index pairs; edge_intensity keeps the max
    observed intensity per firm pair.  `ambiguous` (only if keep_amb) holds
    (detected_cands, partner_cands, dilution, intensity) for soft scoring."""
    def cand_hc(h, c):
        return [p.index for p in peaks
                if abs(p.h_ppm - h) <= tol_h and abs(p.c_ppm - c) <= tol_c]

    def cand_c(c):
        return [p.index for p in peaks if abs(p.c_ppm - c) <= tol_c]

    edges = set()
    edge_intensity: Dict[Tuple[int, int], float] = {}
    ambiguous: List[Tuple[tuple, tuple, float, float]] = []
    firm = amb = unmatched = 0
    for (c2, c1, h1, inten) in rows:
        det = cand_hc(h1, c1)     # detected methyl (proton + carbon)
        par = cand_c(c2)          # partner methyl (carbon only)
        if not det or not par:
            unmatched += 1
            continue
        if len(det) == 1 and len(par) == 1 and det[0] != par[0]:
            i, j = det[0], par[0]
            key = (min(i, j), max(i, j))
            edges.add(key)
            edge_intensity[key] = max(edge_intensity.get(key, 0.0), inten)
            firm += 1
        else:
            amb += 1
            if keep_amb:
                dil = 1.0 / (len(det) * len(par))
                ambiguous.append((tuple(det), tuple(par), dil, inten))
    return edges, edge_intensity, ambiguous, {'firm': firm, 'ambiguous': amb, 'unmatched': unmatched}


def match_noesy(peaks, rows, tol_h, tol_c, keep_amb=False):
    return _match(peaks, rows, tol_h, tol_c, keep_amb)


def match_hmbc(peaks, rows, tol_h, tol_c):
    edges, _inten, _amb, stats = _match(peaks, rows, tol_h, tol_c, keep_amb=False)
    return edges, stats
