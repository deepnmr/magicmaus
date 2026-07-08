"""Run the fresh magicmaus (mm/) across ALL examples and compare the three ideas
residue-wise: MAUS (hard bounds, no ranking) vs MAGIC (scoring over the full
space) vs magicmaus (scoring within the MAUS bounds).

Every example is driven through the mmCIF structure path (PDB inputs are
converted to a minimal mmCIF first), the BMRB-backed ones get N(0,sigma) peak
lists (sigma_H=0.02, sigma_C=0.1), and TNF-alpha runs as a real homotrimer
(multimer, inter-subunit NOEs).

    python compare_all.py [--seed N]
"""
from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from mm import structure as st
from mm import simulate as sim
from mm.compare import three_way

ONE_TO_THREE = {o: t for t, o in st.THREE_TO_ONE.items()}
EX = Path('examples')

# name, structure file, source kwargs.  BMRB examples simulate noisy peaks
# (bmrb=.str or shifts_tsv).  TNF-α uses its REAL experimental peak lists picked
# from .ucsf spectra (real=True), so its row is genuine data, not simulation.
EXAMPLES = [
    ('ubq',  '1UBQ.pdb', dict(bmrb='bmr6457.str')),
    ('hnh',  '6O56.pdb', dict(bmrb='bmr27949.str')),
    ('il2',  '1M47.pdb', dict(bmrb='bmr28104.str')),
    ('mbp',  '1ANF.pdb', dict(bmrb='bmr7114.str')),
    ('rec2', '4CMP.pdb', dict(bmrb='bmr28105.str')),
    ('rec3', '4ZT0.pdb', dict(bmrb='bmr28110.str')),
    ('msg',  '1D8C.pdb', dict(shifts_tsv='msg_methyl_shifts.tsv')),
    ('TNFa', 'fold_tnfa_trimer_model_0.cif',
     dict(real=True, hmqc='hmqc.tsv', noesy='noesy_intensity.tsv',
          hmbc='hmbc.tsv', truth='hmqc_true.tsv')),
]


def pdb_to_cif(pdb_path: Path, cif_path: Path) -> None:
    """Emit a minimal mmCIF with just the _atom_site loop the parser reads.
    Keeps every chain (so a multimer stays a multimer) and the first model only."""
    rows, serial = [], 1
    for line in pdb_path.read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        if line[16] not in (' ', 'A'):
            continue
        atom = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip() or '.'
        resi = line[22:26].strip()
        x, y, z = line[30:38].strip(), line[38:46].strip(), line[46:54].strip()
        elem = (line[76:78].strip() or atom[:1])
        rows.append(f'ATOM {serial} {elem} {atom} . {resn} {chain} 1 {resi} ? '
                    f'{x} {y} {z} 1.00 0.00 {resi} {chain} 1')
        serial += 1
    header = ('data_x\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n'
              '_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n'
              '_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n'
              '_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n'
              '_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n'
              '_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n'
              '_atom_site.auth_asym_id\n_atom_site.pdbx_PDB_model_num\n')
    cif_path.write_text(header + '\n'.join(rows) + '\n')


def run(seed: int, only=None, skip=()):
    results = []
    todo = [e for e in EXAMPLES if (not only or e[0] in only) and e[0] not in skip]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for name, struct, src in todo:
            sp = EX / name / struct
            # always drive the mmCIF path (convert PDB -> minimal cif)
            if sp.suffix == '.cif':
                cif = sp
            else:
                cif = tmp / f'{name}.cif'
                pdb_to_cif(sp, cif)
            out = tmp / name
            if src.get('real'):                      # real experimental peak lists
                c = three_way(str(cif), str(EX / name / src['hmqc']),
                              str(EX / name / src['noesy']), str(EX / name / src['hmbc']),
                              str(EX / name / src['truth']), seed=seed, name=name)
            else:                                    # simulate noisy peaks
                src_kw = {k: str(EX / name / v) for k, v in src.items()}
                sim.simulate(str(cif), out_dir=str(out), noe_cut=8.0, seed=seed, **src_kw)
                # MSG's exact enumeration is intractable (257 peaks, 138-wide
                # Leu/Val domains); use the arc-consistency-pruned domains there.
                c = three_way(str(cif), str(out / 'hmqc.tsv'), str(out / 'noesy.tsv'),
                              str(out / 'hmbc.tsv'), str(out / 'hmqc_true.tsv'),
                              seed=seed, name=name, enumerate_options=(name != 'msg'))
            results.append(c)
            tag = ' [real data]' if src.get('real') else ''
            print(f'  {name:5s} done (n={c.n_peaks}, multimer={c.multimer}){tag}')
    return results


def report(results):
    real_names = {e[0] for e in EXAMPLES if e[2].get('real')}
    print('\n=== residue-wise assignment accuracy (cif structures) ===')
    print('data: sim = simulated peaks (true shifts + N(0,sigma)); real = experimental peak list')
    hdr = f'{"example":6s} {"n":>4s} {"mm?":>4s} {"data":>5s} | {"MAUS":>6s} {"MAGIC":>6s} {"magicmaus":>9s} | {"envelope":>8s}'
    print(hdr)
    print('-' * len(hdr))
    agg = {'maus': 0.0, 'magic': 0.0, 'mm': 0.0, 'env': 0.0}
    wins = 0
    for c in results:
        mm_flag = 'yes' if c.multimer else ''
        data = 'real' if c.name in real_names else 'sim'
        print(f'{c.name:6s} {c.n_peaks:4d} {mm_flag:>4s} {data:>5s} | '
              f'{c.maus_single:5.1f}% {c.magic_single:5.1f}% {c.magicmaus_single:8.1f}% | '
              f'{c.envelope:7.1f}%')
        agg['maus'] += c.maus_single
        agg['magic'] += c.magic_single
        agg['mm'] += c.magicmaus_single
        agg['env'] += c.envelope
        if c.magicmaus_single >= max(c.maus_single, c.magic_single) - 1e-9:
            wins += 1
    k = len(results)
    print('-' * len(hdr))
    print(f'{"MEAN":6s} {"":4s} {"":>4s} | '
          f'{agg["maus"]/k:5.1f}% {agg["magic"]/k:5.1f}% {agg["mm"]/k:8.1f}% | '
          f'{agg["env"]/k:7.1f}%')
    print(f'\nmagicmaus >= both baselines on {wins}/{k} examples.')
    print(f'mean residue-wise: MAUS {agg["maus"]/k:.1f}%  MAGIC {agg["magic"]/k:.1f}%  '
          f'magicmaus {agg["mm"]/k:.1f}%')
    return agg['mm'] / k, agg['maus'] / k, agg['magic'] / k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--only', nargs='*', help='run only these example names')
    ap.add_argument('--skip', nargs='*', default=['msg'],
                    help="skip these (default: msg, ~15 min; pass '' to include it)")
    a = ap.parse_args()
    results = run(a.seed, only=a.only, skip=set(a.skip or ()))
    mm, maus, magic = report(results)
    assert mm >= maus and mm >= magic, 'magicmaus must lead on mean residue-wise accuracy'
    print('\nOK: magicmaus leads residue-wise on the mean across the examples run.')


if __name__ == '__main__':
    main()
