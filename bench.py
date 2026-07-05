"""Benchmark MAUS and magicmaus on one example directory, on the SAME intensity
NOESY.  MAGIC (sibling engine) is scored too if an assignments file is passed.

Usage:
    python bench.py EXDIR PDB [--magic assignments.tsv]
      EXDIR/{hmqc.tsv, hmqc_true.tsv, noesy_intensity.tsv, hmbc.tsv?}
"""
import argparse
import re
from pathlib import Path

import maus
import magicmaus as mm


def residue(label):
  m = re.match(r'^([A-Z])(\d+)', label or '')
  return (m.group(1), int(m.group(2))) if m else None


def run(exdir, pdb, magic_out=None, tol_h=0.01, tol_c=0.05):
  ex = Path(exdir)
  truth = maus.load_truth(ex / 'hmqc_true.tsv')
  peaks = maus.load_hmqc(str(ex / 'hmqc.tsv'))
  n = len(peaks)

  lab = maus.parse_labeling('A;I;L;M;T;V')
  methyls = maus.parse_structure(Path(pdb).read_text().splitlines(), lab)
  lbi = {m.index: m.label for m in methyls}
  crosses = mm.load_noesy_rows(str(ex / 'noesy_intensity.tsv'))
  gem, short_g, long_g = maus.build_structure_graph(methyls, 6.0, 10.0)
  noe, ei, amb, _ = mm.match_noe_intensity(peaks, crosses, tol_h, tol_c)

  # MAUS: option sets (intensity ignored by construction)
  core = maus.MAUS(methyls, peaks, gem, short_g, long_g, noe)
  opts = core.solve_options()
  maus_unique = sum(1 for p in peaks if len(opts[p.index]) == 1)
  maus_unique_ok = sum(1 for p in peaks if len(opts[p.index]) == 1
                       and lbi[opts[p.index][0]] == truth.get(p.peak_id))
  maus_inset = sum(1 for p in peaks if truth.get(p.peak_id) in [lbi[g] for g in opts[p.index]])

  def score_mm(soft):
    eng = mm.MagicMaus(methyls, peaks, gem, short_g, long_g, noe, edge_intensity=ei)
    if soft:
      eng.set_soft_evidence(amb)
    chosen, options = eng.solve()
    meth = sum(1 for p in peaks if chosen[p.index] is not None and lbi[chosen[p.index]] == truth.get(p.peak_id))
    res = sum(1 for p in peaks if residue(lbi.get(chosen[p.index])) == residue(truth.get(p.peak_id)))
    inset = sum(1 for p in peaks if truth.get(p.peak_id) in [lbi[g] for g in options[p.index]])
    return meth, res, inset

  mm_meth, mm_res, mm_inset = score_mm(False)
  mms_meth, mms_res, mms_inset = score_mm(True)

  magic_meth = magic_res = None
  if magic_out and Path(magic_out).exists():
    magic_call = {}
    for line in Path(magic_out).read_text().splitlines()[1:]:
      f = line.split('\t')
      if len(f) > 3:
        magic_call[f[0]] = f[3]
    magic_meth = sum(1 for p in peaks if magic_call.get(p.peak_id) == truth.get(p.peak_id))
    magic_res = sum(1 for p in peaks if residue(magic_call.get(p.peak_id)) == residue(truth.get(p.peak_id)))

  return {
    'n': n, 'maus_unique': maus_unique, 'maus_unique_ok': maus_unique_ok,
    'maus_inset': maus_inset, 'mm_meth': mm_meth, 'mm_res': mm_res,
    'mm_inset': mm_inset, 'mms_meth': mms_meth, 'mms_res': mms_res,
    'mms_inset': mms_inset, 'magic_meth': magic_meth, 'magic_res': magic_res,
  }


def main(argv=None):
  ap = argparse.ArgumentParser()
  ap.add_argument('exdir')
  ap.add_argument('pdb')
  ap.add_argument('--magic', default=None)
  a = ap.parse_args(argv)
  r = run(a.exdir, a.pdb, a.magic)
  n = r['n']
  pct = lambda x: f'{x}/{n}={100*x/n:.1f}%'
  print(f"n={n}")
  if r['magic_meth'] is not None:
    print(f"MAGIC        methyl {pct(r['magic_meth'])}  residue {pct(r['magic_res'])}")
  print(f"MAUS unique  methyl {pct(r['maus_unique_ok'])} (of {r['maus_unique']} decisive)  envelope {pct(r['maus_inset'])}")
  print(f"magicmaus    methyl {pct(r['mm_meth'])}  residue {pct(r['mm_res'])}  envelope {pct(r['mm_inset'])}")
  print(f"  +soft-amb  methyl {pct(r['mms_meth'])}  residue {pct(r['mms_res'])}  envelope {pct(r['mms_inset'])}")
  return 0


if __name__ == '__main__':
  raise SystemExit(main())
