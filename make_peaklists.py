"""Build maus/magicmaus peak lists from real data: a PDB structure + a BMRB
NMR-STAR chemical-shift file.  Nothing about the answer is baked into the
indexing -- HMQC peaks are anonymised (P1..Pn); the truth key is written
separately for scoring only.

Pipeline
--------
* Structure methyls come from the PDB (reusing maus.parse_structure): one node
  per methyl carbon of the labelled types, with coordinates and a canonical
  label (e.g. L7CD1).
* Chemical shifts come from the BMRB _Atom_chem_shift loop: for each methyl the
  carbon shift is the carbon Val; the proton shift is the mean of the methyl
  protons (HD11/HD12/HD13 for CD1, HG21.. for CG2, HB1.. for CB, etc.).
* A methyl is *observed* iff the structure has its carbon AND the BMRB has both
  its carbon and proton shift.  Observed methyls become HMQC peaks; the truth
  key maps each anonymous peak to its structural methyl.
* The NOESY is a simulated 3D (H)CCH network: every pair of observed methyls
  within ``--noe-cut`` in the structure emits two cross peaks (both directions),
  ``label C1 C2 H2`` = (partner carbon, observed carbon, observed proton) in
  chemical-shift space, so degeneracy is real and re-matched by frequency at run
  time -- exactly what maus.py consumes.
* An optional HMBC-HMQC list links geminal Leu/Val methyl pairs.

Outputs (into --out-dir): hmqc.tsv, hmqc_true.tsv, noesy.tsv, hmbc.tsv.
Run make_intensity_noesy.py afterwards for noesy_intensity.tsv.

Usage:
    python make_peaklists.py PDB BMRB.str --out-dir examples/xxx [--noe-cut 8.0]
"""
from __future__ import annotations

import argparse
import math
from itertools import combinations
from pathlib import Path

import maus  # parse_structure, parse_labeling, METHYL_ATOMS, THREE_TO_ONE


def parse_bmrb_shifts(star_path: str) -> dict:
  """Return {(res_num, atom_id): val} from the _Atom_chem_shift loop.

  NMR-STAR column *order* varies between entries (some carry extra tags such as
  Entity_assembly_asym_ID / Ambiguity_set_ID), so columns are resolved by tag
  name, not position.  Residue number prefers Auth_seq_ID (matches PDB author
  numbering) and falls back to Comp_index_ID."""
  lines = Path(star_path).read_text().splitlines()
  cols: dict = {}          # tag suffix -> column index
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
  ci_resn = cols['Comp_ID']
  ci_atom = cols['Atom_ID']
  ci_val = cols['Val']
  ncol = len(cols)

  shifts: dict = {}
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


def methyl_proton_shift(shifts: dict, resnum: int, carbon: str):
  """Mean shift of the three methyl protons attached to `carbon`.
  CD1 -> HD1{1,2,3}; CG2 -> HG2{1,2,3}; CB -> HB{,1,2,3}; CE -> HE{,1,2,3}.
  Proton Atom_IDs share the prefix 'H'+carbon[1:]."""
  prefix = 'H' + carbon[1:]
  vals = [v for (rn, a), v in shifts.items()
          if rn == resnum and a.startswith(prefix)]
  return sum(vals) / len(vals) if vals else None


def build(pdb: str, star: str, labeling: dict, noe_cut: float):
  methyls = maus.parse_structure(Path(pdb).read_text().splitlines(), labeling)
  shifts = parse_bmrb_shifts(star)

  observed = []          # (methyl, h_ppm, c_ppm)
  for m in methyls:
    c = shifts.get((m.res_num, m.atom))
    h = methyl_proton_shift(shifts, m.res_num, m.atom)
    if c is None or h is None:
      continue
    observed.append((m, round(h, 3), round(c, 3)))

  observed.sort(key=lambda t: (t[0].res_num, t[0].atom))
  peak_of = {id(m): f'P{i+1}' for i, (m, _h, _c) in enumerate(observed)}
  return observed, peak_of


def write_lists(observed, peak_of, out_dir: str, noe_cut: float):
  out = Path(out_dir)
  out.mkdir(parents=True, exist_ok=True)

  with open(out / 'hmqc.tsv', 'w') as fh, open(out / 'hmqc_true.tsv', 'w') as ft:
    fh.write('label\tH_ppm\tC_ppm\tres_type\n')
    ft.write('label\tH_ppm\tC_ppm\tres_type\tTrue\n')
    for m, h, c in observed:
      pid = peak_of[id(m)]
      fh.write(f'{pid}\t{h}\t{c}\t{m.res_type}\n')
      ft.write(f'{pid}\t{h}\t{c}\t{m.res_type}\t{m.label}\n')

  # NOESY: every observed pair within noe_cut -> two (H)CCH cross peaks.
  n = 0
  with open(out / 'noesy.tsv', 'w') as f:
    f.write('label\tC1\tC2\tH2\n')
    for (a, ha, ca), (b, hb, cb) in combinations(observed, 2):
      if math.dist(a.coord, b.coord) >= noe_cut:
        continue
      n += 1
      f.write(f'X{n}\t{cb}\t{ca}\t{ha}\n')   # observe a, partner b
      f.write(f'X{n}b\t{ca}\t{cb}\t{hb}\n')  # observe b, partner a

  # HMBC: geminal Leu/Val links (same residue, the two methyls of the pair).
  by_res = {}
  for m, h, c in observed:
    by_res.setdefault((m.res_type, m.res_num), []).append((m, h, c))
  g = 0
  with open(out / 'hmbc.tsv', 'w') as f:
    f.write('label\tC1\tC2\tH2\n')
    for pair in by_res.values():
      if len(pair) != 2:
        continue
      (a, ha, ca), (b, hb, cb) = pair
      if a.res_type not in ('L', 'V'):
        continue
      g += 1
      f.write(f'B{g}\t{cb}\t{ca}\t{ha}\n')
      f.write(f'B{g}b\t{ca}\t{cb}\t{hb}\n')
  return len(observed), n, g


def main(argv=None):
  ap = argparse.ArgumentParser(description='Build peak lists from PDB + BMRB.')
  ap.add_argument('pdb')
  ap.add_argument('bmrb', help='BMRB NMR-STAR .str file')
  ap.add_argument('--out-dir', required=True)
  ap.add_argument('--noe-cut', type=float, default=8.0,
                  help='methyl-methyl distance (A) that produces a NOESY cross peak')
  ap.add_argument('--labeling', default='A;I;L;M;T;V')
  args = ap.parse_args(argv)

  labeling = maus.parse_labeling(args.labeling)
  observed, peak_of = build(args.pdb, args.bmrb, labeling, args.noe_cut)
  n_hmqc, n_noe_pairs, n_gem = write_lists(observed, peak_of, args.out_dir, args.noe_cut)

  by_type = {}
  for m, _h, _c in observed:
    by_type[m.res_type] = by_type.get(m.res_type, 0) + 1
  comp = ' '.join(f'{by_type[t]} {t}' for t in sorted(by_type))
  print(f'observed methyls (HMQC peaks) = {n_hmqc}  ({comp})')
  print(f'NOESY pairs < {args.noe_cut} A = {n_noe_pairs}  -> {2*n_noe_pairs} cross peaks')
  print(f'HMBC geminal Leu/Val links   = {n_gem}')
  print(f'wrote hmqc.tsv hmqc_true.tsv noesy.tsv hmbc.tsv to {args.out_dir}/')
  return 0


if __name__ == '__main__':
  raise SystemExit(main())
