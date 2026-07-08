"""Convert a maus/magicmaus-format dataset into MAGIC's control-file bundle,
so all three engines can be scored on the *same* NOE network with intensities.

maus HMQC  : label  H_ppm  C_ppm  res_type
MAGIC HMQC : label  C_ppm  H_ppm  res_type            (carbon first)

maus NOESY : label  C1  C2  H2  [intensity]           (C1=partner, C2/H2=observed)
MAGIC NOESY: <nuclei header> / <tol header> / label donorC refC refH intensity
             (donorC=C1, refC=C2, refH=H2 -> identical column order)

Usage:
    python convert_to_magic.py HMQC.tsv NOESY_intensity.tsv PDB OUTDIR \
        [--tol-h 0.01] [--tol-c 0.05]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import maus  # noqa: E402

ONE_TO_THREE = {o: t for t, o in maus.THREE_TO_ONE.items()}


def render_pdb(pdb_lines) -> str:
  """Emit a multi-chain PDB whose per-atom records carry one image per subunit,
  so MAGIC (which keys atoms by residue only) sees every chain as a symmetry
  image. Reuses maus's parsers, so a .cif trimer becomes a MAGIC-readable PDB
  and a monomer round-trips to one chain. Space-separated to stay
  split()-parseable by MAGIC's structure reader."""
  lines = list(pdb_lines)
  is_cif = any(l.startswith('_atom_site.') for l in lines)
  coords = maus._coords_from_cif(lines) if is_cif else maus._coords_from_pdb(lines)
  out, serial = [], 1
  for (one, resi) in sorted(coords, key=lambda k: (k[1], k[0])):
    for atom, images in coords[(one, resi)].items():
      for k, (x, y, z) in enumerate(images):
        chain = chr(ord('A') + k)
        out.append(f'ATOM {serial:>6} {atom:<4} {ONE_TO_THREE[one]} {chain} '
                   f'{resi:>4} {x:8.3f} {y:8.3f} {z:8.3f} 1.00 0.00')
        serial += 1
  return '\n'.join(out) + '\n'


def main(argv=None):
  ap = argparse.ArgumentParser()
  ap.add_argument('hmqc'); ap.add_argument('noesy'); ap.add_argument('pdb')
  ap.add_argument('outdir')
  ap.add_argument('--tol-h', type=float, default=0.01)
  ap.add_argument('--tol-c', type=float, default=0.05)
  a = ap.parse_args(argv)
  out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

  # HMQC: swap to carbon-first.
  peaks = maus.load_hmqc(a.hmqc)
  hm = [f'{p.peak_id} {p.c_ppm:.3f} {p.h_ppm:.3f} {p.res_type}' for p in peaks]
  (out / 'hmqc.list').write_text('\n'.join(hm) + '\n')

  # NOESY: same column order, prepend MAGIC's two header lines.
  noe_body = []
  for line in Path(a.noesy).read_text().splitlines():
    if not line.strip() or line.startswith('#') or line.startswith('label'):
      continue
    f = line.split('\t')
    label, c1, c2, h2 = f[0], f[1], f[2], f[3]
    inten = f[4] if len(f) > 4 and f[4].strip() else '1.0'
    noe_body.append(f'{label} {c1} {c2} {h2} {inten}')
  header = f'13C;13C;1H\n{a.tol_c};{a.tol_c};{a.tol_h}'
  (out / 'noesy.list').write_text(header + '\n' + '\n'.join(noe_body) + '\n')

  # Sequence: one entry per methyl-bearing residue.
  lab = maus.parse_labeling('A;I;L;M;T;V')
  methyls = maus.parse_structure(Path(a.pdb).read_text().splitlines(), lab)
  residues = sorted({(m.res_type, m.res_num) for m in methyls}, key=lambda r: r[1])
  (out / 'seq.txt').write_text('\n'.join(f'{t}{n}' for t, n in residues) + '\n')

  # PDB: render a multi-chain PDB (all subunits as symmetry images) so MAGIC's
  # min-over-images contacts pick up inter-subunit NOEs, and .cif input works.
  (out / 'model.pdb').write_text(render_pdb(Path(a.pdb).read_text().splitlines()))
  control = (
    'HMQC = hmqc.list\nNOESY = noesy.list\nPDB = model.pdb\nSEQ = seq.txt\n'
    'LABELING = A;I,CD1;L,CD1,CD2;M;T;V,CG1,CG2\nGEMINAL = 1\n'
    'CUTOFF_FACTOR = 1.0\nDISTANCE_LIMITS = 6 10\nSCORE_TOL_END = true\n'
  )
  (out / 'control.txt').write_text(control)
  print(f'wrote MAGIC bundle to {out}/  ({len(peaks)} HMQC peaks, {len(noe_body)} NOE rows)')


if __name__ == '__main__':
  raise SystemExit(main())
