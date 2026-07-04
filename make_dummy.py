"""Generate a tiny, hand-designed DUMMY dataset for the magicmaus tutorial.

Eight methyls on a toy 'protein'.  Coordinates and shifts are chosen so the run
exhibits every outcome the tutorial explains:

  * peaks the hard constraints pin uniquely,
  * a geminal Leu pair that stays a genuine 2-way ambiguity (an achiral NOE
    network cannot tell CD1 from CD2), and
  * shift-degenerate peaks a boolean network cannot separate but intensities can.

Writes examples/dummy/{model.pdb, hmqc.tsv, hmqc_true.tsv, noesy.tsv,
noesy_intensity.tsv}.  Deterministic; no randomness.
"""
import math
from pathlib import Path

# label, res_type, res_num, atom, (x, y, z), H_ppm, C_ppm
METHYLS = [
  ('A10CB',  'A', 10, 'CB',  (0.0, 0.0, 0.0), 1.30, 19.0),
  ('A20CB',  'A', 20, 'CB',  (0.0, 0.0, 3.0), 1.42, 20.5),
  ('I30CD1', 'I', 30, 'CD1', (3.0, 0.0, 1.5), 0.80, 13.0),
  ('L40CD1', 'L', 40, 'CD1', (3.0, 2.0, 3.0), 0.75, 24.0),
  ('L40CD2', 'L', 40, 'CD2', (3.0, 2.0, 5.0), 0.75, 24.0),   # shift-degenerate with CD1 -> geminal ambiguity
  ('V50CG1', 'V', 50, 'CG1', (6.0, 1.0, 2.0), 0.90, 21.0),
  ('V50CG2', 'V', 50, 'CG2', (6.0, 1.0, 4.0), 0.95, 22.5),
  ('T60CG2', 'T', 60, 'CG2', (2.0, 4.0, 0.0), 1.10, 21.8),
]

GEMINAL = {'L40CD1': 'CD2', 'L40CD2': 'CD1', 'V50CG1': 'CG2', 'V50CG2': 'CG1'}
SHORT_CUT, LONG_CUT = 6.0, 10.0

out = Path('examples/dummy')
out.mkdir(parents=True, exist_ok=True)

# --- PDB (only the methyl carbons are needed by the parser) ---
pdb = []
serial = 1
for (label, rt, rn, atom, (x, y, z), h, c) in METHYLS:
  three = {'A': 'ALA', 'I': 'ILE', 'L': 'LEU', 'M': 'MET', 'T': 'THR', 'V': 'VAL'}[rt]
  pdb.append(
    f'ATOM  {serial:>5} {atom:<4} {three} A{rn:>4}    '
    f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C')
  serial += 1
pdb.append('END')
(out / 'model.pdb').write_text('\n'.join(pdb) + '\n')

# --- HMQC peak lists (anonymous P-labels; truth key keeps the real methyl) ---
hmqc, true = ['label\tH_ppm\tC_ppm\tres_type'], ['label\tH_ppm\tC_ppm\tres_type\tTrue']
peakid = {}
for i, (label, rt, rn, atom, xyz, h, c) in enumerate(METHYLS, start=1):
  pid = f'P{i}'
  peakid[label] = pid
  hmqc.append(f'{pid}\t{h:.3f}\t{c:.3f}\t{rt}')
  true.append(f'{pid}\t{h:.3f}\t{c:.3f}\t{rt}\t{label}')
(out / 'hmqc.tsv').write_text('\n'.join(hmqc) + '\n')
(out / 'hmqc_true.tsv').write_text('\n'.join(true) + '\n')

# --- NOESY 3D (H)CCH cross peaks from structural contacts (both directions) ---
coord = {m[0]: m[4] for m in METHYLS}
shift = {m[0]: (m[5], m[6]) for m in METHYLS}
contacts = []
for a in range(len(METHYLS)):
  for b in range(len(METHYLS)):
    if a == b:
      continue
    la, lb = METHYLS[a][0], METHYLS[b][0]
    same_res = METHYLS[a][2] == METHYLS[b][2] and METHYLS[a][1] == METHYLS[b][1]
    geminal = same_res and GEMINAL.get(la) == METHYLS[b][3]
    d = math.dist(coord[la], coord[lb])
    if geminal or d <= LONG_CUT:
      contacts.append((la, lb, d))

boolean, inten = ['label\tC1\tC2\tH2'], ['label\tC1\tC2\tH2\tintensity']
k = 1
dmax6 = max(1.0 / (d ** 6) for (_, _, d) in contacts if d > 0)
for (obs, par, d) in contacts:
  h2, c2 = shift[obs]           # observed methyl: proton + carbon
  _, c1 = shift[par]           # partner methyl: carbon only
  boolean.append(f'X{k}\t{c1:.3f}\t{c2:.3f}\t{h2:.3f}')
  I = (1.0 / (d ** 6)) / dmax6 if d > 0 else 1.0
  inten.append(f'X{k}\t{c1:.3f}\t{c2:.3f}\t{h2:.3f}\t{I:.4e}')
  k += 1
(out / 'noesy.tsv').write_text('\n'.join(boolean) + '\n')
(out / 'noesy_intensity.tsv').write_text('\n'.join(inten) + '\n')

print(f'wrote {out}/  ({len(METHYLS)} methyls, {len(contacts)} directed NOE cross peaks)')
