"""Generate every number the manuscript's Tables 1-3 need, for one target, in the
fresh mm/ residue-wise framework.  Emits JSON.

    python manuscript_stats.py <target>

Metrics per target (all on the same 1/r^6 intensity NOESY + HMBC + HMQC):
  residue-wise (correct residue, geminal swap ignored): MAGIC, MAUS-forced,
    magicmaus, magicmaus+soft; envelope (truth residue retained).
  methyl-level (atom-exact): magicmaus+soft.
  MAUS residue-decisive fraction (option set collapses to one residue).
  per residue type (residue-wise, magicmaus+soft) and per methyl carbon
  (methyl-level, magicmaus+soft).
"""
import argparse
import json
import re
import tempfile
from pathlib import Path

from mm import simulate as sim, structure as st, peaks as pk
from mm.engine import MagicMaus
from compare_all import pdb_to_cif, EX, EXAMPLES


def build(name):
    cfg = {e[0]: e for e in EXAMPLES}[name]
    _, struct, src = cfg
    tmp = Path(tempfile.mkdtemp())
    sp = EX / name / struct
    cif = sp if sp.suffix == '.cif' else tmp / 's.cif'
    if sp.suffix != '.cif':
        pdb_to_cif(sp, cif)
    if src.get('real'):
        f = lambda k: str(EX / name / src[k])
        return str(cif), f('hmqc'), f('noesy'), f('hmbc'), f('truth'), True
    out = tmp / name
    src_kw = {k: str(EX / name / v) for k, v in src.items()}
    sim.simulate(str(cif), out_dir=str(out), noe_cut=8.0, seed=0, **src_kw)
    return (str(cif), str(out / 'hmqc.tsv'), str(out / 'noesy.tsv'),
            str(out / 'hmbc.tsv'), str(out / 'hmqc_true.tsv'), False)


def analyze(name):
    cif, hmqc, noesy, hmbc, truth, real = build(name)
    lab = pk.parse_labeling('A;I;L;M;T;V')
    methyls = st.parse_structure(Path(cif).read_text().splitlines(), lab)
    multimer = any(len(m.images) > 1 for m in methyls)
    peaks = pk.load_hmqc(hmqc)
    gem, sh, lo = st.build_structure_graph(methyls, 6.0, 10.0)
    noe, inten, amb, nstat = pk.match_noesy(peaks, pk.load_triple(noesy), 0.02, 0.10, keep_amb=True)
    gl, hstat = pk.match_hmbc(peaks, pk.load_triple(hmbc), 0.02, 0.10)
    eng = MagicMaus(methyls, peaks, gem, sh, lo, noe, gl, edge_intensity=inten)

    res_of = {m.index: (m.res_type, m.res_num) for m in methyls}
    atom_of = {m.index: m.atom for m in methyls}
    by_label = {m.label: m.index for m in methyls}
    truth_map = pk.load_truth(truth)
    tmethyl = {p.index: by_label.get(truth_map.get(p.peak_id, '')) for p in peaks}
    tres = {p.index: res_of.get(tmethyl[p.index]) for p in peaks}
    n = len(peaks)

    enum = (name != 'msg')            # MSG exact enumeration is intractable
    options = eng.option_sets() if enum else {p.index: list(eng.core.domain[p.index]) for p in peaks}

    def res_ok(call):
        return sum(1 for p in peaks if call.get(p.index) is not None
                   and tres[p.index] is not None and res_of[call[p.index]] == tres[p.index])

    def methyl_ok(call):
        return sum(1 for p in peaks if call.get(p.index) is not None
                   and tmethyl[p.index] is not None and call[p.index] == tmethyl[p.index])

    env = sum(1 for p in peaks if tres[p.index] is not None
              and tres[p.index] in {res_of[g] for g in options[p.index]})
    maus_call = {p.index: (options[p.index][0] if options[p.index] else None) for p in peaks}
    resdec = [p for p in peaks if options[p.index] and len({res_of[g] for g in options[p.index]}) == 1]
    resdec_ok = sum(1 for p in resdec if tres[p.index] == res_of[options[p.index][0]])

    eng.enforce_contacts = True
    mm_call = eng.optimize(options, seed=0)
    eng.set_soft_evidence(amb)
    mm_soft = eng.optimize(options, seed=0)
    eng.amb, eng.amb_by_peak = [], {}
    full = {p.index: list(eng.core.sites_by_type.get(p.res_type, [])) for p in peaks}
    eng.enforce_contacts = False
    magic_call = eng.optimize(full, seed=0)
    eng.enforce_contacts = True

    # per-type / per-carbon are computed from the BASE magicmaus call (the featured
    # config: HMQC + 3D NOESY intensity + 3D HMBC, no soft-ambiguous), so every
    # table reconciles with the Table-1 magicmaus column.
    TYPES = ['I', 'L', 'V', 'A', 'T', 'M']
    per_type = {}
    for t in TYPES:
        idxs = [p.index for p in peaks if tres[p.index] and tres[p.index][0] == t]
        if idxs:
            per_type[t] = [sum(1 for i in idxs if mm_call.get(i) is not None
                               and res_of[mm_call[i]] == tres[i]), len(idxs)]
    CARBONS = [('I', 'CD1'), ('L', 'CD1'), ('L', 'CD2'), ('V', 'CG1'), ('V', 'CG2'),
               ('A', 'CB'), ('T', 'CG2'), ('M', 'CE')]
    per_methyl = {}
    for (t, a) in CARBONS:
        idxs = [p.index for p in peaks if tmethyl[p.index] is not None
                and res_of[tmethyl[p.index]][0] == t and atom_of[tmethyl[p.index]] == a]
        if idxs:
            per_methyl[f'{t}{a}'] = [sum(1 for i in idxs if mm_call.get(i) == tmethyl[i]), len(idxs)]

    pct = lambda a, b: round(100.0 * a / b, 1) if b else 0.0
    return {
        'name': name, 'n': n, 'multimer': multimer, 'real': real,
        'firm_noe': len(noe), 'ambiguous': nstat['ambiguous'], 'hmbc_firm': hstat['firm'],
        'envelope': pct(env, n),
        'magic': pct(res_ok(magic_call), n),
        'maus_forced': pct(res_ok(maus_call), n),
        'maus_resdec_cov': pct(len(resdec), n), 'maus_resdec_acc': pct(resdec_ok, len(resdec)),
        'magicmaus': pct(res_ok(mm_call), n),
        'magicmaus_soft': pct(res_ok(mm_soft), n),
        'magicmaus_methyl': pct(methyl_ok(mm_call), n),
        'magicmaus_methyl_soft': pct(methyl_ok(mm_soft), n),
        'per_type': per_type, 'per_methyl': per_methyl,
    }


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('target')
    a = ap.parse_args()
    print(json.dumps(analyze(a.target)))
