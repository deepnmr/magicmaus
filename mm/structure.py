"""Structure layer: parse an mmCIF (or legacy PDB) into methyl nodes and build
the distance-classified structure graph, with homo-oligomer support.

A homo-oligomer's subunits are chemically identical, so its methyls collapse to
one node per residue/atom in the asymmetric unit; every chain contributes a
coordinate *image* of that node.  An NOE contact distance is the closest
approach over all image pairs, so an inter-subunit contact is picked up
automatically (this is what makes the assignment "multimer-aware").
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Optional, Tuple

# three-letter -> one-letter for the methyl-bearing residue types we assign
THREE_TO_ONE = {'ALA': 'A', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'THR': 'T', 'VAL': 'V'}

# residue one-letter -> [(methyl carbon atom, geminal-partner atom or None)]
METHYL_ATOMS: Dict[str, List[Tuple[str, Optional[str]]]] = {
    'A': [('CB', None)],
    'I': [('CD1', None)],
    'L': [('CD1', 'CD2'), ('CD2', 'CD1')],
    'M': [('CE', None)],
    'T': [('CG2', None)],
    'V': [('CG1', 'CG2'), ('CG2', 'CG1')],
}

Coord = Tuple[float, float, float]


@dataclass(frozen=True)
class Methyl:
    index: int
    label: str          # canonical structure id, e.g. 'L135CD1'
    res_type: str       # one-letter
    res_num: int
    atom: str           # methyl carbon atom name
    geminal_atom: str   # partner methyl atom name, or '' if none
    images: Tuple[Coord, ...]  # one coordinate per subunit (>=1); [0] is reference

    @property
    def coord(self) -> Coord:
        return self.images[0]


def min_dist(a: Methyl, b: Methyl) -> float:
    """Closest approach between two methyls over all chain images (intra- or
    inter-subunit)."""
    return min(math.dist(ca, cb) for ca in a.images for cb in b.images)


def _coords_from_pdb(lines):
    """(one_letter, resnum) -> atom -> [xyz per chain], plus the set of chain
    ids.  All chains kept (one image per subunit); altloc other than blank/A
    skipped.  Only the FIRST model is read, so a multi-model NMR ensemble does
    not fabricate bogus per-conformer 'chain images'."""
    coords: Dict[Tuple[str, int], Dict[str, list]] = {}
    chains = set()
    for line in lines:
        if line.startswith('ENDMDL'):
            break                                   # keep only the first model
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        resn = line[17:20].strip()
        if resn not in THREE_TO_ONE:
            continue
        if line[16] not in (' ', 'A'):
            continue
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        chains.add(line[21])
        atom = line[12:16].strip()
        coords.setdefault((THREE_TO_ONE[resn], resi), {}).setdefault(atom, []).append(
            (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords, chains


def _coords_from_cif(lines):
    """Parse the mmCIF `_atom_site` loop by column name (order varies between
    depositions).  Keeps ATOM records of the first model, every chain (one image
    per subunit), non-alt (or altloc A).  Returns (coords, chain_ids)."""
    coords: Dict[Tuple[str, int], Dict[str, list]] = {}
    chains = set()
    cols: Dict[str, int] = {}
    in_loop, first_model = False, None
    it = iter(lines)
    for line in it:
        s = line.strip()
        if s == 'loop_':
            cols, in_loop = {}, True
            continue
        if in_loop and s.startswith('_atom_site.'):
            cols[s.split('.', 1)[1]] = len(cols)
            continue
        if in_loop and cols and s and not s.startswith('_'):
            while True:
                if not s or s.startswith(('_', '#', 'loop_')):
                    break
                f = s.split()
                if len(f) >= len(cols):
                    def g(name, default=None):
                        i = cols.get(name)
                        return f[i] if i is not None else default
                    if g('group_PDB', 'ATOM') == 'ATOM':
                        resn = g('auth_comp_id') or g('label_comp_id')
                        model = g('pdbx_PDB_model_num')
                        alt = g('label_alt_id', '.')
                        if first_model is None:
                            first_model = model
                        if (resn in THREE_TO_ONE and model == first_model
                                and alt in ('.', '?', 'A')):
                            try:
                                resi = int(g('auth_seq_id') or g('label_seq_id'))
                            except (TypeError, ValueError):
                                resi = None
                            if resi is not None:
                                atom = g('label_atom_id') or g('auth_atom_id')
                                chains.add(g('auth_asym_id') or g('label_asym_id'))
                                coords.setdefault((THREE_TO_ONE[resn], resi), {}) \
                                    .setdefault(atom, []).append(
                                        (float(g('Cartn_x')), float(g('Cartn_y')),
                                         float(g('Cartn_z'))))
                s = next(it, '').strip()
            in_loop = False
    return coords, chains


def parse_coords(lines):
    """Return (coords, chain_ids)."""
    lines = list(lines)
    is_cif = any(l.startswith('_atom_site.') for l in lines)
    return _coords_from_cif(lines) if is_cif else _coords_from_pdb(lines)


def parse_structure(lines, labeling: Dict[str, List[Tuple[str, Optional[str]]]]) -> List[Methyl]:
    """Build the methyl node list for the labelled residue types.

    A homo-oligomer collapses to one node per (res_type, res_num) with one
    coordinate image per subunit, which requires the symmetric chains to share
    author residue numbering (the standard homo-oligomer convention, e.g. the
    AlphaFold trimer here).  ponytail: numbering-based collapse; chains with
    per-subunit continuous numbering would need sequence alignment — warned
    below rather than handled, since it is non-standard for a homo-oligomer.
    """
    coords, chains = parse_coords(lines)
    methyls: List[Methyl] = []
    for (one, resi) in sorted(coords, key=lambda k: (k[1], k[0])):
        for atom, gem in labeling.get(one, []):
            if atom not in coords[(one, resi)]:
                continue
            imgs = tuple(coords[(one, resi)][atom])   # one per chain, reference first
            methyls.append(Methyl(
                index=len(methyls),
                label=f'{one}{resi}{atom}',
                res_type=one, res_num=resi, atom=atom,
                geminal_atom=gem or '', images=imgs))
    if len(chains) > 1 and methyls and max(len(m.images) for m in methyls) < 2:
        import sys
        print(f'warning: {len(chains)} chains present but no methyl collapsed to '
              f'multiple images — subunits use non-shared residue numbering, so '
              f'inter-subunit (multimer) NOE contacts are NOT included. Renumber '
              f'the chains to a shared scheme to enable them.', file=sys.stderr)
    return methyls


def build_structure_graph(methyls: List[Methyl], short_cut: float, long_cut: float):
    """Return three symmetric edge-index sets: geminal / short (<short_cut) /
    long (<=long_cut).  Distances are min-over-images so inter-subunit contacts
    of a homo-oligomer are included."""
    gem, short, long = set(), set(), set()
    for a, b in combinations(methyls, 2):
        same_res = a.res_num == b.res_num and a.res_type == b.res_type
        if same_res and a.atom == b.geminal_atom:
            gem.add((a.index, b.index)); gem.add((b.index, a.index))
            continue
        d = min_dist(a, b)
        if d < short_cut:
            short.add((a.index, b.index)); short.add((b.index, a.index))
        elif d <= long_cut:
            long.add((a.index, b.index)); long.add((b.index, a.index))
    return gem, short, long
