"""Score MSG once, reusing a single MAUS option-set enumeration for the plain
and +soft calls (the hard-constraint option sets are identical; only the scoring
differs), plus one more enumeration for the +HMBC column. Writes the committed
call files and prints the benchmark row. MSG is large (257 methyls, weak firm-NOE
constraints), so the enumeration is the slow part -- doing it twice instead of
four times roughly halves the wall-clock.
"""
import re
from pathlib import Path

import maus
import magicmaus as mm

PDB = 'examples/msg/1D8C.pdb'
EX = 'examples/msg'
LAB = 'I;L;V'


def residue(label):
  m = re.match(r'^([A-Z])(\d+)', label or '')
  return (m.group(1), int(m.group(2))) if m else None


def main():
  truth = maus.load_truth(f'{EX}/hmqc_true.tsv')
  peaks = maus.load_hmqc(f'{EX}/hmqc.tsv')
  n = len(peaks)
  lab = maus.parse_labeling(LAB)
  methyls = maus.parse_structure(Path(PDB).read_text().splitlines(), lab)
  lbi = {m.index: m.label for m in methyls}
  gem, sg, lg = maus.build_structure_graph(methyls, 6.0, 10.0)
  cross = mm.load_noesy_rows(f'{EX}/noesy_intensity.tsv')
  noe, ei, amb, _ = mm.match_noe_intensity(peaks, cross, 0.01, 0.05)

  def score(chosen, options):
    meth = sum(1 for p in peaks if chosen[p.index] is not None
               and lbi[chosen[p.index]] == truth.get(p.peak_id))
    inset = sum(1 for p in peaks if truth.get(p.peak_id) in [lbi[g] for g in options[p.index]])
    return meth, inset

  def write_calls(path, chosen, options):
    with open(path, 'w') as f:
      f.write('label\tres_type\tcall\tn_options\toptions\ttruth\tcall_correct\ttruth_in_set\n')
      for p in peaks:
        g = chosen[p.index]
        call = lbi.get(g, '')
        opts = [lbi[x] for x in options[p.index]]
        t = truth.get(p.peak_id, '')
        f.write(f'{p.peak_id}\t{p.res_type}\t{call}\t{len(opts)}\t{",".join(opts)}\t'
                f'{t}\t{int(call == t)}\t{int(t in opts)}\n')

  # --- no-HMBC engine: one enumeration, reused for plain + soft ---
  eng = mm.MagicMaus(methyls, peaks, gem, sg, lg, noe, edge_intensity=ei)
  options = eng.option_sets()
  unique = sum(1 for p in peaks if len(options[p.index]) == 1)
  unique_ok = sum(1 for p in peaks if len(options[p.index]) == 1
                  and lbi[options[p.index][0]] == truth.get(p.peak_id))

  plain = eng.optimize(options)
  p_meth, p_inset = score(plain, options)

  eng.set_soft_evidence(amb)
  soft = eng.optimize(options)
  s_meth, s_inset = score(soft, options)
  write_calls(f'{EX}/magicmaus_calls.tsv', soft, options)

  # --- +HMBC engine: separate enumeration (different hard CNF) ---
  gem_links, _ = maus.match_hmbc(peaks, maus.load_hmbc(f'{EX}/hmbc.tsv'), 0.01, 0.05)
  engh = mm.MagicMaus(methyls, peaks, gem, sg, lg, noe, gem_links=gem_links, edge_intensity=ei)
  oh = engh.option_sets()
  engh.set_soft_evidence(amb)
  hmbc = engh.optimize(oh)
  h_meth, h_inset = score(hmbc, oh)
  write_calls(f'{EX}/magicmaus_calls_hmbc.tsv', hmbc, oh)

  pc = lambda x: f'{x}/{n}={100*x/n:.1f}%'
  print(f'MSG n={n}')
  print(f'MAUS unique  {pc(unique_ok)} (of {unique} decisive)  envelope {pc(p_inset)}')
  print(f'magicmaus    {pc(p_meth)}  envelope {pc(p_inset)}')
  print(f'  +soft      {pc(s_meth)}  envelope {pc(s_inset)}')
  print(f'  +HMBC      {pc(h_meth)}  envelope {pc(h_inset)}')


if __name__ == '__main__':
  main()
