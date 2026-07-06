"""Convert the TNFa Sparky peak lists (examples/TNFa/*.list) into magicmaus
input TSVs, plus a separate assignment answer key.

Sparky lists (whitespace-separated, first line = header):
  13C_HMQC : Assignment  w1(C)  w2(H)  Height
  NOESY    : Assignment  w1(C1) w2(C2) w3(H)  Height  Note
  HMBC     : Assignment  w1(C2) w2(C1) w3(H)  Height  Note   (assign = C2-C1-H1)

magicmaus TSVs (see maus.load_hmqc / load_ccch):
  hmqc.tsv         : label  H_ppm  C_ppm  res_type            (res_type = bare letter)
  hmqc_true.tsv    : label  H_ppm  C_ppm  res_type  True      (answer, e.g. V1CG1)
  noesy.tsv        : label  C1  C2  H2                         (observed=(H2,C2), partner=C1)
  noesy_intensity.tsv : label  C1  C2  H2  intensity
  hmbc.tsv         : label  C1  C2  H2

Usage:  python make_tnfa_input.py [examples/TNFa]
"""

import re
import sys
from pathlib import Path

ILVATM = set('AILMTV')  # methyl-bearing types magicmaus knows


def read_list(path):
  """Yield (assignment, [float cols], note) for each data row."""
  for line in Path(path).read_text().splitlines():
    s = line.strip()
    if not s or s.startswith('Assignment'):
      continue
    tok = s.split()
    assign = tok[0]
    # numeric cols run until the first non-float (Data Height is an int -> float ok).
    nums, i = [], 1
    while i < len(tok):
      try:
        nums.append(float(tok[i]))
      except ValueError:
        break
      i += 1
    note = ' '.join(tok[i:])
    yield assign, nums, note


def w(x):
  return f'{x:.3f}'


def main(argv=None):
  d = Path(argv[0]) if argv else Path('examples/TNFa')

  # --- HMQC -> input + answer key ---
  hmqc, skipped = [], []
  for assign, n, _note in read_list(d / 'TNFa_ILVAT_13C_HMQC.list'):
    c_ppm, h_ppm = n[0], n[1]
    m = re.match(r'^([A-Z])(\d+)(C[A-Z0-9]*)-', assign)  # V1CG1-HG1
    if not m:
      skipped.append(assign)                              # ?-? , vat85C-H (ambiguous)
      continue
    one = m.group(1)
    if one not in ILVATM:
      skipped.append(assign)
      continue
    true = f'{m.group(1)}{m.group(2)}{m.group(3)}'        # V1CG1
    hmqc.append((h_ppm, c_ppm, one, true))

  lines_in = ['label\tH_ppm\tC_ppm\tres_type']
  lines_tr = ['label\tH_ppm\tC_ppm\tres_type\tTrue']
  for i, (h, c, one, true) in enumerate(hmqc, 1):
    lab = f'P{i}'
    lines_in.append(f'{lab}\t{w(h)}\t{w(c)}\t{one}')
    lines_tr.append(f'{lab}\t{w(h)}\t{w(c)}\t{one}\t{true}')
  (d / 'hmqc.tsv').write_text('\n'.join(lines_in) + '\n')
  (d / 'hmqc_true.tsv').write_text('\n'.join(lines_tr) + '\n')

  # --- NOESY -> C1=w1, C2=w2, H2=w3 ; drop diagonals ---
  noe, noe_i = ['label\tC1\tC2\tH2'], ['label\tC1\tC2\tH2\tintensity']
  k = 0
  for assign, n, note in read_list(d / 'TNFa_ILVAT_NOESY_HMQC.list'):
    c1, c2, h2, height = n[0], n[1], n[2], n[3]
    if 'Diagonal' in note or abs(c1 - c2) < 1e-3:          # loader drops these anyway
      continue
    k += 1
    lab = f'X{k}'
    noe.append(f'{lab}\t{w(c1)}\t{w(c2)}\t{w(h2)}')
    noe_i.append(f'{lab}\t{w(c1)}\t{w(c2)}\t{w(h2)}\t{height:.6g}')
  (d / 'noesy.tsv').write_text('\n'.join(noe) + '\n')
  (d / 'noesy_intensity.tsv').write_text('\n'.join(noe_i) + '\n')

  # --- HMBC -> C1(partner)=w2, C2(observed)=w1, H2=w3 ; drop diagonals ---
  hmbc = ['label\tC1\tC2\tH2']
  b = 0
  for assign, n, note in read_list(d / 'TNFa_ILVAT_HMBC_HMQC.list'):
    c2, c1, h2 = n[0], n[1], n[2]                          # assign order C2-C1-H1
    if 'Diagonal' in note or abs(c1 - c2) < 1e-3:
      continue
    b += 1
    hmbc.append(f'B{b}\t{w(c1)}\t{w(c2)}\t{w(h2)}')
  (d / 'hmbc.tsv').write_text('\n'.join(hmbc) + '\n')

  print(f'HMQC : {len(hmqc)} peaks  (skipped {len(skipped)}: {skipped})')
  print(f'NOESY: {k} cross-peaks')
  print(f'HMBC : {b} rows')
  print(f'wrote hmqc.tsv, hmqc_true.tsv, noesy.tsv, noesy_intensity.tsv, hmbc.tsv -> {d}/')
  print('NOTE: magicmaus also needs a TNFa PDB (not in the peak lists) for the run.')


if __name__ == '__main__':
  raise SystemExit(main(sys.argv[1:]))
