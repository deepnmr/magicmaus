"""Simulate the three experimental peak lists from a structure + a chemical-shift
source, with realistic measurement scatter.

Each observed methyl is drawn ONE measured shift = true_shift + N(0, sigma), with
sigma_H = 0.02 ppm and sigma_C = 0.1 ppm by default (the goal's error model).  A
single draw per methyl is used consistently across HMQC / NOESY / HMBC — each
methyl has one resonance in the sample; the assignment ambiguity comes from
distinct methyls whose measured shifts fall within the matching tolerance (the
carbon-only NOE partner is the dominant source, exactly as in a real dataset).

Shift source is a BMRB NMR-STAR `_Atom_chem_shift` loop, or a paired-shift TSV
(`resnum res_type atom H_ppm C_ppm`, or a maus-style hmqc_true.tsv).

A homo-oligomer collapses to one methyl per residue/atom (its subunits are
equivalent, so they share a shift); NOE distances and 1/r^6 intensities are the
min over all chain-image pairs, so inter-subunit contacts appear.

Outputs (into out_dir):  hmqc.tsv  hmqc_true.tsv  noesy.tsv  hmbc.tsv
The 3D lists use the  label  C2  C1  H1  [intensity]  column order.
"""
from __future__ import annotations

import re
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import peaks as pk
from . import structure as st

SIGMA_H = 0.02   # 1H measurement error (ppm), normal distribution
SIGMA_C = 0.10   # 13C measurement error (ppm), normal distribution


# ---------------------------------------------------------------------------
# chemical-shift sources
# ---------------------------------------------------------------------------
def parse_bmrb_shifts(star_path: str) -> Dict[Tuple[int, str], float]:
    """{(res_num, atom_id): val} from the BMRB `_Atom_chem_shift` loop.  Columns
    resolved by tag name (order varies); res_num prefers Auth_seq_ID."""
    lines = Path(star_path).read_text().splitlines()
    cols: Dict[str, int] = {}
    tag_end = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith('_Atom_chem_shift.'):
            cols[s.split('.', 1)[1]] = len(cols)
            tag_end = i
    if tag_end is None:
        raise ValueError('no _Atom_chem_shift loop found')
    ci_comp = cols['Comp_index_ID']
    ci_auth = cols.get('Auth_seq_ID')
    ci_atom = cols['Atom_ID']
    ci_val = cols['Val']
    ncol = len(cols)
    shifts: Dict[Tuple[int, str], float] = {}
    for ln in lines[tag_end + 1:]:
        s = ln.strip()
        if not s or s.startswith(('_', 'loop_')):
            continue
        if s in ('stop_', 'save_'):
            break
        f = s.split()
        if len(f) != ncol:
            continue
        resnum = None
        if ci_auth is not None:
            try:
                resnum = int(f[ci_auth])
            except ValueError:
                resnum = None
        if resnum is None:
            try:
                resnum = int(f[ci_comp])
            except ValueError:
                continue
        try:
            val = float(f[ci_val])
        except ValueError:
            continue
        shifts[(resnum, f[ci_atom])] = val
    return shifts


def _methyl_proton_shift(shifts, resnum, carbon) -> Optional[float]:
    """Mean of the methyl protons on `carbon` (CD1 -> HD1{1,2,3}, etc.)."""
    prefix = 'H' + carbon[1:]
    vals = [v for (rn, a), v in shifts.items() if rn == resnum and a.startswith(prefix)]
    return sum(vals) / len(vals) if vals else None


def bmrb_hc(shifts):
    def get(m: st.Methyl):
        c = shifts.get((m.res_num, m.atom))
        h = _methyl_proton_shift(shifts, m.res_num, m.atom)
        return (h, c) if (c is not None and h is not None) else None
    return get


def parse_shift_tsv(path: str) -> Dict[Tuple[int, str], Tuple[float, float]]:
    """Paired methyl shift table -> {(res_num, atom): (h, c)}.  Accepts either
    `resnum res_type atom H_ppm C_ppm` or a maus hmqc_true.tsv
    (`label H_ppm C_ppm res_type True`, atom taken from the True label)."""
    hc: Dict[Tuple[int, str], Tuple[float, float]] = {}
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#') or s.startswith(('resnum', 'label')):
            continue
        f = s.split('\t')
        m = re.match(r'^([A-Z])(\d+)(C[A-Z0-9]+)$', f[-1].strip())
        if m and len(f) >= 5:                       # hmqc_true.tsv row
            hc[(int(m.group(2)), m.group(3))] = (float(f[1]), float(f[2]))
        elif len(f) >= 5:                           # resnum res_type atom H C
            hc[(int(f[0]), f[2])] = (float(f[3]), float(f[4]))
    return hc


def tsv_hc(hc):
    def get(m: st.Methyl):
        return hc.get((m.res_num, m.atom))
    return get


# ---------------------------------------------------------------------------
# simulation
# ---------------------------------------------------------------------------
def observed_methyls(methyls, get_hc):
    """[(methyl, true_h, true_c)] for every methyl with both shifts available."""
    obs = []
    for m in methyls:
        v = get_hc(m)
        if v and v[0] is not None and v[1] is not None:
            obs.append((m, float(v[0]), float(v[1])))
    obs.sort(key=lambda t: (t[0].res_num, t[0].atom))
    return obs


def add_noise(obs, sigma_h, sigma_c, seed):
    """One measured shift per methyl: true + N(0, sigma).  Returns
    {methyl_index: (h_meas, c_meas)}."""
    rng = np.random.default_rng(seed)
    meas: Dict[int, Tuple[float, float]] = {}
    for m, h, c in obs:
        meas[m.index] = (round(h + rng.normal(0.0, sigma_h), 4),
                         round(c + rng.normal(0.0, sigma_c), 4))
    return meas


def write_lists(obs, meas, out_dir, noe_cut, with_intensity=True):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    peak_of = {m.index: f'P{i + 1}' for i, (m, _h, _c) in enumerate(obs)}

    with open(out / 'hmqc.tsv', 'w') as fh, open(out / 'hmqc_true.tsv', 'w') as ft:
        fh.write('label\tH_ppm\tC_ppm\tres_type\n')
        ft.write('label\tH_ppm\tC_ppm\tres_type\tTrue\n')
        for m, _h, _c in obs:
            h, c = meas[m.index]
            fh.write(f'{peak_of[m.index]}\t{h}\t{c}\t{m.res_type}\n')
            ft.write(f'{peak_of[m.index]}\t{h}\t{c}\t{m.res_type}\t{m.label}\n')

    # NOESY: each observed pair within noe_cut (min over chain images) emits two
    # cross peaks (both directions).  label  C2  C1  H1  [intensity]:
    # detected methyl = (C1,H1); partner carbon = C2.
    n = 0
    with open(out / 'noesy.tsv', 'w') as f:
        f.write('label\tC2\tC1\tH1\tintensity\n' if with_intensity else 'label\tC2\tC1\tH1\n')
        for (a, _ha, _ca), (b, _hb, _cb) in combinations(obs, 2):
            d = st.min_dist(a, b)
            if d >= noe_cut:
                continue
            n += 1
            inten = 1.0 / (d ** 6) if (with_intensity and d > 0) else None
            ha, ca = meas[a.index]
            hb, cb = meas[b.index]
            row = f'X{n}\t{cb}\t{ca}\t{ha}'          # detect a, partner b
            f.write(row + (f'\t{inten:.6e}\n' if inten is not None else '\n'))
            row = f'X{n}b\t{ca}\t{cb}\t{hb}'         # detect b, partner a
            f.write(row + (f'\t{inten:.6e}\n' if inten is not None else '\n'))

    # HMBC-HMQC: geminal Leu/Val links.  label  C2  C1  H1: detected methyl =
    # (C1,H1); C2 = the same-residue geminal partner carbon.
    by_res: Dict[Tuple[str, int], list] = {}
    for m, _h, _c in obs:
        by_res.setdefault((m.res_type, m.res_num), []).append(m)
    g = 0
    with open(out / 'hmbc.tsv', 'w') as f:
        f.write('label\tC2\tC1\tH1\n')
        for pair in by_res.values():
            if len(pair) != 2 or pair[0].res_type not in ('L', 'V'):
                continue
            a, b = pair
            ha, ca = meas[a.index]
            hb, cb = meas[b.index]
            g += 1
            f.write(f'B{g}\t{cb}\t{ca}\t{ha}\n')     # detect a, geminal partner b
            f.write(f'B{g}b\t{ca}\t{cb}\t{hb}\n')     # detect b, geminal partner a
    return len(obs), n, g


def simulate(structure: str, *, bmrb: str = None, shifts_tsv: str = None,
             out_dir: str, noe_cut: float = 8.0, labeling: str = 'A;I;L;M;T;V',
             sigma_h: float = SIGMA_H, sigma_c: float = SIGMA_C, seed: int = 0,
             with_intensity: bool = True):
    """Simulate hmqc/noesy/hmbc peak lists + truth key from a structure and a
    shift source (BMRB or paired TSV).  Returns (n_hmqc, n_noe_pairs, n_gem)."""
    lab = pk.parse_labeling(labeling)
    methyls = st.parse_structure(Path(structure).read_text().splitlines(), lab)
    if bmrb:
        get_hc = bmrb_hc(parse_bmrb_shifts(bmrb))
    elif shifts_tsv:
        get_hc = tsv_hc(parse_shift_tsv(shifts_tsv))
    else:
        raise ValueError('provide bmrb= or shifts_tsv=')
    obs = observed_methyls(methyls, get_hc)
    meas = add_noise(obs, sigma_h, sigma_c, seed)
    return write_lists(obs, meas, out_dir, noe_cut, with_intensity)
