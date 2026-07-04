"""Score MAGIC, MAUS and magicmaus on the SAME intensity NOESY, head to head."""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'maus'))
import maus  # noqa: E402
import magicmaus as mm  # noqa: E402

PDB = 'examples/mbp/1ANF.pdb'
HMQC = 'examples/mbp/hmqc.tsv'
NOESY = 'mbp_noesy_intensity.tsv'
TRUTH = 'examples/mbp/hmqc_true.tsv'
MAGIC_OUT = 'magic_run_intensity/out/Output/assignments.tsv'


def residue(label):
  """L75CD1 / L75CD -> ('L', 75); residue-level identity, geminal-agnostic."""
  m = re.match(r'^([A-Z])(\d+)', label or '')
  return (m.group(1), int(m.group(2))) if m else None


truth = maus.load_truth(TRUTH)
peaks = maus.load_hmqc(HMQC)
n = len(peaks)

# --- MAGIC ---
magic_call = {}
for line in Path(MAGIC_OUT).read_text().splitlines()[1:]:
  f = line.split('\t')
  magic_call[f[0]] = f[3]
magic_res = sum(1 for p in peaks if residue(magic_call.get(p.peak_id)) == residue(truth.get(p.peak_id)))
magic_meth = sum(1 for p in peaks if magic_call.get(p.peak_id) == truth.get(p.peak_id))

# --- MAUS + magicmaus on the same intensity data ---
lab = maus.parse_labeling('A;I;L;M;T;V')
methyls = maus.parse_structure(Path(PDB).read_text().splitlines(), lab)
crosses = mm.load_noesy_rows(NOESY)
gem, short_g, long_g = maus.build_structure_graph(methyls, 6.0, 10.0)
noe, ei, amb, _ = mm.match_noe_intensity(peaks, crosses, 0.01, 0.05)
lbi = {m.index: m.label for m in methyls}

# MAUS: option sets (intensity ignored by construction)
core = maus.MAUS(methyls, peaks, gem, short_g, long_g, noe)
opts = core.solve_options()
maus_unique = sum(1 for p in peaks if len(opts[p.index]) == 1)
maus_unique_ok = sum(1 for p in peaks if len(opts[p.index]) == 1
                     and lbi[opts[p.index][0]] == truth.get(p.peak_id))
maus_inset = sum(1 for p in peaks if truth.get(p.peak_id) in [lbi[g] for g in opts[p.index]])

# magicmaus with intensity (+ soft-ambiguous)
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

pct = lambda x: f'{x}/{n} = {100*x/n:4.1f}%'
print(f'{"":22s}{"methyl-level":>18s}{"residue-level":>18s}{"truth-in-envelope":>20s}')
print(f'{"MAGIC (scoring)":22s}{pct(magic_meth):>18s}{pct(magic_res):>18s}{"— (no envelope)":>20s}')
print(f'{"MAUS (SAT)":22s}{"unique "+pct(maus_unique_ok):>18s}{"":>18s}{pct(maus_inset):>20s}')
print(f'{"magicmaus":22s}{pct(mm_meth):>18s}{pct(mm_res):>18s}{pct(mm_inset):>20s}')
print(f'{"magicmaus +soft-amb":22s}{pct(mms_meth):>18s}{pct(mms_res):>18s}{pct(mms_inset):>20s}')
print(f'\nMAUS decisive on {maus_unique}/{n} peaks (rest abstain); MAGIC & magicmaus commit on all {n}.')
