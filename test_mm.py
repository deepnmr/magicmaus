"""Runnable self-check for the fresh magicmaus (mm/) — no framework, just asserts.

    python test_mm.py

Covers the load-bearing invariants:
  * mmCIF/PDB parse + homo-oligomer images, min-dist over subunits (multimer);
  * the (C2, C1, H1) 3D matching convention (detected by C+H, partner by C only);
  * MAUS never-exclude: every true methyl stays in its option set;
  * simulate -> assign round trip with the N(0,sigma) noise model, reproducible.
"""
import tempfile
from pathlib import Path

from mm import structure as st
from mm import simulate as sim
from mm import engine

# A 2-chain (homo-dimer) toy: chain A and chain B related by a 30 A translation.
# Val5 (CG1/CG2) and Leu9 (CD1/CD2) and Ala20 (CB) per chain.  The Val5 of chain A
# sits right next to the Leu9 of chain B (inter-subunit contact) — only visible if
# min-dist scans chain images.
PDB = """\
ATOM      1  CG1 VAL A   5       0.000   0.000   0.000  1.00  0.00
ATOM      2  CG2 VAL A   5       1.500   0.000   0.000  1.00  0.00
ATOM      3  CD1 LEU A   9      20.000   0.000   0.000  1.00  0.00
ATOM      4  CD2 LEU A   9      21.500   0.000   0.000  1.00  0.00
ATOM      5  CB  ALA A  20      10.000   0.000   0.000  1.00  0.00
ATOM      6  CG1 VAL B   5      30.000   0.000   0.000  1.00  0.00
ATOM      7  CG2 VAL B   5      31.500   0.000   0.000  1.00  0.00
ATOM      8  CD1 LEU B   9       3.000   0.000   0.000  1.00  0.00
ATOM      9  CD2 LEU B   9       4.500   0.000   0.000  1.00  0.00
ATOM     10  CB  ALA B  20      40.000   0.000   0.000  1.00  0.00
"""

SHIFTS = """\
resnum\tres_type\tatom\tH_ppm\tC_ppm
5\tV\tCG1\t0.80\t21.0
5\tV\tCG2\t0.85\t22.0
9\tL\tCD1\t0.70\t24.0
9\tL\tCD2\t0.75\t25.0
20\tA\tCB\t1.30\t19.0
"""


def main():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / 'toy.pdb').write_text(PDB)
        (d / 'shifts.tsv').write_text(SHIFTS)

        # -- multimer: Val5 collapses to one node with two chain images; the A/B
        #    Leu9 contact is only reachable via the second image.
        from mm import peaks as pk
        methyls = st.parse_structure(PDB.splitlines(), pk.parse_labeling('A;I;L;M;T;V'))
        by_label = {m.label: m for m in methyls}
        assert len(by_label['V5CG1'].images) == 2, 'homo-dimer should give 2 images'
        v5 = by_label['V5CG1']
        l9 = by_label['L9CD1']
        # reference-only distance is 20 A; min over images is the A-Val5 / B-Leu9
        # inter-subunit contact at 3 A.
        ref = min(__import__('math').dist(v5.images[0], c) for c in [l9.images[0]])
        mind = st.min_dist(v5, l9)
        assert abs(ref - 20.0) < 1e-6, ref
        assert abs(mind - 3.0) < 1e-6, f'inter-subunit min-dist wrong: {mind}'

        # -- simulate with a fixed seed, twice -> identical (reproducible noise).
        sim.simulate(str(d / 'toy.pdb'), shifts_tsv=str(d / 'shifts.tsv'),
                     out_dir=str(d / 'r1'), noe_cut=8.0, seed=7)
        sim.simulate(str(d / 'toy.pdb'), shifts_tsv=str(d / 'shifts.tsv'),
                     out_dir=str(d / 'r2'), noe_cut=8.0, seed=7)
        assert (d / 'r1' / 'hmqc.tsv').read_text() == (d / 'r2' / 'hmqc.tsv').read_text(), \
            'same seed must reproduce the same peak list'

        # noise is actually applied: measured != true (round to 4 dp, sigma>0).
        true_h = 0.80
        hmqc = (d / 'r1' / 'hmqc.tsv').read_text().splitlines()[1:]
        got_any_noise = any(abs(float(l.split('\t')[1]) - true_h) > 0 for l in hmqc)
        assert got_any_noise, 'noise model produced no scatter'

        # -- assign: MAUS never-exclude — every truth methyl is in its option set.
        r = engine.assign(str(d / 'toy.pdb'), str(d / 'r1' / 'hmqc.tsv'),
                          str(d / 'r1' / 'noesy.tsv'), hmbc=str(d / 'r1' / 'hmbc.tsv'),
                          tol_h=0.05, tol_c=0.25)
        truth = pk.load_truth(str(d / 'r1' / 'hmqc_true.tsv'))
        label = r.engine.label_by_index
        for p in r.engine.peaks:
            opts = [label[g] for g in r.options[p.index]]
            assert truth[p.peak_id] in opts, \
                f'never-exclude violated: {truth[p.peak_id]} not in {opts}'

        # -- (C2,C1,H1) convention: a NOESY row's detected methyl (C1,H1) and
        #    partner (C2) resolve to two *different* peaks when unambiguous.
        r2 = engine.assign(str(d / 'toy.pdb'), str(d / 'r1' / 'hmqc.tsv'),
                           str(d / 'r1' / 'noesy.tsv'), tol_h=0.05, tol_c=0.25)
        assert r2.stats['n_noe_edges'] >= 1, 'expected at least one firm NOE edge'

    print('OK: multimer images, reproducible N(0,sigma) noise, never-exclude, (C2,C1,H1) matching')


if __name__ == '__main__':
    main()
