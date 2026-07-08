"""magicmaus CLI:  python -m mm assign ...   |   python -m mm simulate ...

  assign   — run the MAUS+MAGIC engine on HMQC / 3D NOESY / 3D HMBC-HMQC lists.
  simulate — build noisy peak lists from a structure + BMRB (or shift TSV).
"""
from __future__ import annotations

import argparse
import math

from . import engine, simulate as sim
from . import peaks as pk
from .compare import three_way


def _cmd_assign(a):
    r = engine.assign(a.structure, a.hmqc, a.noesy, a.hmbc,
                      short_cut=a.short_cut, long_cut=a.long_cut,
                      tol_h=a.tol_h, tol_c=a.tol_c, labeling=a.labeling,
                      soft_ambiguous=a.soft_ambiguous, restarts=a.restarts, seed=a.seed)
    eng, chosen, options, sstat = r.engine, r.chosen, r.options, r.stats
    label = eng.label_by_index
    truth = pk.load_truth(a.truth) if a.truth else {}
    peaks = eng.peaks
    n = len(peaks)

    tiers = {'unique': 0, 'scored': 0, 'ambiguous': 0, 'unassigned': 0}
    rows = []
    for p in peaks:
        tier, margin = eng.confidence(p.index, chosen, options)
        tiers[tier] += 1
        call = label.get(chosen[p.index], '')
        opt_labels = [label[g] for g in options[p.index]]
        rows.append((p, call, opt_labels, tier, margin, truth.get(p.peak_id, '')))

    if a.out:
        with open(a.out, 'w') as f:
            f.write('label\tres_type\tcall\tconfidence\tmargin\tn_options\toptions\ttruth\t'
                    'call_correct\ttruth_in_set\n')
            for (p, call, opts, tier, margin, t) in rows:
                cc = '' if not t else int(call == t)
                tis = '' if not t else int(t in opts)
                mval = 'inf' if margin == math.inf else f'{margin:.3e}'
                f.write(f'{p.peak_id}\t{p.res_type}\t{call}\t{tier}\t{mval}\t{len(opts)}\t'
                        f'{",".join(opts)}\t{t}\t{cc}\t{tis}\n')

    ns = sstat['noe']
    print(f"methyls(G nodes)={sstat['n_methyls']}  HMQC peaks={n}  "
          f"firm NOE edges={sstat['n_noe_edges']}")
    print(f"G edges: geminal={sstat['G']['gem']} short={sstat['G']['short']} long={sstat['G']['long']}")
    print(f"NOE match (tol H+-{a.tol_h}/C+-{a.tol_c}): firm={ns['firm']} "
          f"ambiguous={ns['ambiguous']} unmatched={ns['unmatched']}")
    if sstat['hmbc']:
        hs = sstat['hmbc']
        print(f"HMBC geminal links: firm={hs['firm']} ambiguous={hs['ambiguous']} unmatched={hs['unmatched']}")
    if eng.n_tentative:
        print(f"tentative anchors used = {eng.n_tentative}")
    sizes = [len(options[p.index]) for p in peaks]
    print('--- MAUS envelope (never excludes truth) ---')
    print(f"unique(1 option)       = {sum(s == 1 for s in sizes)}/{n}")
    print(f"ambiguous(2-3 options) = {sum(2 <= s <= 3 for s in sizes)}/{n}")
    print(f"ambiguous(>3 options)  = {sum(s > 3 for s in sizes)}/{n}")
    print(f"unassigned             = {sum(s == 0 for s in sizes)}/{n}")
    print('--- magicmaus commitment (single coherent call) ---')
    tline = f"confidence: unique={tiers['unique']}  scored={tiers['scored']}  ambiguous={tiers['ambiguous']}"
    if tiers['unassigned']:
        tline += f"  unassigned={tiers['unassigned']}"
    print(tline)

    if truth:
        in_set = sum(1 for (_, _, opts, _, _, t) in rows if t and t in opts)
        call_ok = sum(1 for (_, call, _, _, _, t) in rows if t and call == t)
        by_tier = {k: [0, 0] for k in tiers}
        for (_, call, _, tier, _, t) in rows:
            if not t:
                continue
            by_tier[tier][1] += 1
            by_tier[tier][0] += int(call == t)
        print('--- scored vs truth ---')
        # never-exclude holds only when every firm edge is correct; measurement
        # noise can mis-resolve a "firm" NOE and prune the truth, so only claim it
        # when the option sets actually retained every truth.
        tag = '(never-exclude preserved)' if in_set == n else \
            f'(WARNING: {n - in_set} truth pruned by a wrong firm edge under noise)'
        print(f"truth in MAUS option set = {in_set}/{n} = {100 * in_set / n:.1f}%  {tag}")
        print(f"magicmaus single call    = {call_ok}/{n} = {100 * call_ok / n:.1f}% correct")
        for k in ('unique', 'scored', 'ambiguous', 'unassigned'):
            ok, tot = by_tier[k]
            if tot:
                print(f"    {k:9s}: {ok}/{tot} = {100 * ok / tot:.1f}% correct")
    else:
        print('(no --truth given; scoring skipped)')
    return 0


def _cmd_simulate(a):
    n_hmqc, n_noe, n_gem = sim.simulate(
        a.structure, bmrb=a.bmrb, shifts_tsv=a.shifts_tsv, out_dir=a.out_dir,
        noe_cut=a.noe_cut, labeling=a.labeling, sigma_h=a.sigma_h, sigma_c=a.sigma_c,
        seed=a.seed, with_intensity=not a.no_intensity)
    print(f"observed methyls (HMQC peaks) = {n_hmqc}")
    print(f"NOESY pairs < {a.noe_cut} A     = {n_noe}  -> {2 * n_noe} cross peaks")
    print(f"HMBC geminal Leu/Val links    = {n_gem}")
    print(f"noise: sigma_H={a.sigma_h} sigma_C={a.sigma_c} (normal), seed={a.seed}")
    print(f"wrote hmqc.tsv hmqc_true.tsv noesy.tsv hmbc.tsv to {a.out_dir}/")
    return 0


def _cmd_compare(a):
    c = three_way(a.structure, a.hmqc, a.noesy, a.hmbc, a.truth,
                  short_cut=a.short_cut, long_cut=a.long_cut, tol_h=a.tol_h,
                  tol_c=a.tol_c, labeling=a.labeling, seed=a.seed, name='')
    print(f"residue-wise accuracy (n={c.n_peaks} peaks, multimer={c.multimer}):")
    print(f"  MAUS      (hard bounds, no ranking)        = {c.maus_single:5.1f}%")
    print(f"  MAGIC     (soft score, un-pruned full space)= {c.magic_single:5.1f}%")
    print(f"  magicmaus (MAGIC score within MAUS bounds)  = {c.magicmaus_single:5.1f}%")
    print(f"  MAUS option-set envelope (truth retained)   = {c.envelope:5.1f}%")
    print(f"  MAUS residue-unique coverage {c.maus_unique_cov:.0f}% @ {c.maus_unique_acc:.1f}% correct")
    best = max(c.maus_single, c.magic_single, c.magicmaus_single)
    print(f"  -> best: {'magicmaus' if c.magicmaus_single >= best - 1e-9 else 'MAUS' if c.maus_single >= best - 1e-9 else 'MAGIC'}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog='mm', description='magicmaus: MAUS+MAGIC methyl assignment')
    sub = ap.add_subparsers(dest='cmd', required=True)

    pa = sub.add_parser('assign', help='assign methyls from peak lists')
    pa.add_argument('structure', help='mmCIF or PDB (homo-oligomer -> multimer contacts)')
    pa.add_argument('hmqc', help='HMQC TSV: label H_ppm C_ppm res_type')
    pa.add_argument('noesy', help='3D NOESY TSV: label C2 C1 H1 [intensity]')
    pa.add_argument('--hmbc', default=None, help='3D HMBC-HMQC TSV: label C2 C1 H1')
    pa.add_argument('--truth', default=None, help='truth key TSV for scoring')
    pa.add_argument('--short-cut', dest='short_cut', type=float, default=6.0)
    pa.add_argument('--long-cut', dest='long_cut', type=float, default=10.0)
    pa.add_argument('--tol-h', dest='tol_h', type=float, default=0.02)
    pa.add_argument('--tol-c', dest='tol_c', type=float, default=0.10)
    pa.add_argument('--labeling', default='A;I;L;M;T;V')
    pa.add_argument('--soft-ambiguous', dest='soft_ambiguous', action='store_true',
                    help='fold discarded ambiguous NOE rows in as diluted soft tie-breakers')
    pa.add_argument('--restarts', type=int, default=4)
    pa.add_argument('--seed', type=int, default=0)
    pa.add_argument('--out', default=None, help='write per-peak calls TSV here')
    pa.set_defaults(func=_cmd_assign)

    ps = sub.add_parser('simulate', help='simulate noisy peak lists from structure + shifts')
    ps.add_argument('structure', help='mmCIF or PDB')
    ps.add_argument('--bmrb', default=None, help='BMRB NMR-STAR .str shift file')
    ps.add_argument('--shifts-tsv', dest='shifts_tsv', default=None,
                    help='paired shift TSV (resnum res_type atom H C) or hmqc_true.tsv')
    ps.add_argument('--out-dir', dest='out_dir', required=True)
    ps.add_argument('--noe-cut', dest='noe_cut', type=float, default=8.0)
    ps.add_argument('--labeling', default='A;I;L;M;T;V')
    ps.add_argument('--sigma-h', dest='sigma_h', type=float, default=sim.SIGMA_H)
    ps.add_argument('--sigma-c', dest='sigma_c', type=float, default=sim.SIGMA_C)
    ps.add_argument('--seed', type=int, default=0)
    ps.add_argument('--no-intensity', dest='no_intensity', action='store_true',
                    help='omit the 1/r^6 intensity column (boolean NOE network)')
    ps.set_defaults(func=_cmd_simulate)

    pc = sub.add_parser('compare', help='compare MAUS vs MAGIC vs magicmaus residue-wise on one dataset')
    pc.add_argument('structure')
    pc.add_argument('hmqc'); pc.add_argument('noesy')
    pc.add_argument('truth', help='truth key TSV (required for the comparison)')
    pc.add_argument('--hmbc', default=None)
    pc.add_argument('--short-cut', dest='short_cut', type=float, default=6.0)
    pc.add_argument('--long-cut', dest='long_cut', type=float, default=10.0)
    pc.add_argument('--tol-h', dest='tol_h', type=float, default=0.02)
    pc.add_argument('--tol-c', dest='tol_c', type=float, default=0.10)
    pc.add_argument('--labeling', default='A;I;L;M;T;V')
    pc.add_argument('--seed', type=int, default=0)
    pc.set_defaults(func=_cmd_compare)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == '__main__':
    raise SystemExit(main())
