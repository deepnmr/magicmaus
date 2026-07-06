"""Peak-pick the TNF-alpha .ucsf methyl spectra, type the peaks with the
label-selective spectra (removing tag peaks), and emit magicmaus HMQC input.

Type logic (labelings):
  ILVAT  = I,L,V,A,T  (tag-free master list)
  ILV    = I,L,V      (tag-free)
  Val    = V          (has an N-terminal tag -> extra Val peaks)
  Thr    = T,I        (has tag)
  => A,T = ILVAT - ILV ;  T = (A,T) matched in Thr ;  A = (A,T) not in Thr
     V   = ILV matched in Val
     I   = ILV matched in Thr
     L   = ILV - V - I
Tag peaks are the Val/Thr peaks with no ILVAT partner: they never match the
tag-free master, so they are dropped automatically.
"""
import sys
from pathlib import Path

import numpy as np

import ucsf

D = Path('examples/TNFa')


def pick(path, k=25.0, min_h=None):
  """2D local-maxima peak pick. Returns [(c_ppm, h_ppm, height)] sorted by height."""
  data, axes = ucsf.read(path)
  c_ax, h_ax = axes[0], axes[1]
  noise = np.median(np.abs(data))
  thr = k * noise if min_h is None else min_h
  # local maxima vs 8-neighbourhood
  d = data
  m = (d > thr)
  m[1:-1, 1:-1] &= (d[1:-1, 1:-1] >= d[:-2, 1:-1]) & (d[1:-1, 1:-1] >= d[2:, 1:-1]) \
      & (d[1:-1, 1:-1] >= d[1:-1, :-2]) & (d[1:-1, 1:-1] >= d[1:-1, 2:]) \
      & (d[1:-1, 1:-1] >= d[:-2, :-2]) & (d[1:-1, 1:-1] >= d[2:, 2:]) \
      & (d[1:-1, 1:-1] >= d[:-2, 2:]) & (d[1:-1, 1:-1] >= d[2:, :-2])
  m[0, :] = m[-1, :] = m[:, 0] = m[:, -1] = False
  peaks = []
  for i, j in zip(*np.where(m)):
    # parabolic sub-pixel refine in each axis
    di = _refine(d[i - 1, j], d[i, j], d[i + 1, j])
    dj = _refine(d[i, j - 1], d[i, j], d[i, j + 1])
    peaks.append((float(c_ax.ppm(i + di)), float(h_ax.ppm(j + dj)), float(d[i, j])))
  peaks.sort(key=lambda p: -p[2])
  return _dedup(peaks)


def _dedup(peaks, tol_c=0.06, tol_h=0.012):
  """Drop weaker peaks within (tol_c, tol_h) of a kept stronger one (same methyl
  picked twice as shoulders)."""
  kept = []
  for c, h, ht in peaks:                      # already height-desc
    if any(abs(c - kc) <= tol_c and abs(h - kh) <= tol_h for kc, kh, _ in kept):
      continue
    kept.append((c, h, ht))
  return kept


def _refine(ym, y0, yp):
  denom = ym - 2 * y0 + yp
  return 0.0 if denom == 0 else np.clip(0.5 * (ym - yp) / denom, -0.5, 0.5)


def pick_3d(path, k=40.0):
  """3D local-maxima pick. Returns [(ppm_ax0, ppm_ax1, ppm_ax2, height)]."""
  data, axes = ucsf.read(path)
  noise = np.median(np.abs(data))
  thr = k * noise
  c = data[1:-1, 1:-1, 1:-1]
  m = c > thr
  for di in (-1, 0, 1):
    for dj in (-1, 0, 1):
      for dk in (-1, 0, 1):
        if di == dj == dk == 0:
          continue
        sl = data[1 + di:c.shape[0] + 1 + di, 1 + dj:c.shape[1] + 1 + dj,
                  1 + dk:c.shape[2] + 1 + dk]
        m &= c >= sl
  idx = np.argwhere(m)
  peaks = []
  for (i, j, kk) in idx:
    i, j, kk = i + 1, j + 1, kk + 1
    peaks.append((float(axes[0].ppm(i)), float(axes[1].ppm(j)),
                  float(axes[2].ppm(kk)), float(data[i, j, kk])))
  peaks.sort(key=lambda p: -p[3])
  return peaks


def write_noesy(path, k=10.0, top=80, out='noesy_picked.tsv'):
  """Pick a 3D NOESY-HMQC (ax0=C1 partner, ax1=C2 observed, ax2=H); drop diagonal,
  keep the `top` strongest cross peaks.  A low `top` limits false firm hard-
  constraints: auto-picked noise makes the SAT UNSAT above ~15 firm edges."""
  pk = [p for p in pick_3d(path, k) if abs(p[0] - p[1]) >= 0.05][:top]
  lines = ['label\tC1\tC2\tH2\tintensity']
  for i, (c1, c2, h2, ht) in enumerate(pk, 1):
    lines.append(f'X{i}\t{c1:.3f}\t{c2:.3f}\t{h2:.3f}\t{ht:.6g}')
  (D / out).write_text('\n'.join(lines) + '\n')
  print(f'NOESY: kept top {len(pk)} cross peaks -> {D}/{out}')
  return len(pk)


def matches(peak, plist, tol_c=0.10, tol_h=0.02):
  c, h, _ = peak
  return any(abs(c - pc) <= tol_c and abs(h - ph) <= tol_h for pc, ph, _ in plist)


def geminal_partners(ilvat, hmbc, tol_c=0.30, tol_h=0.02):
  """For each master peak, find its geminal partner via the 2D HMBC (each HMBC
  peak = (partner_carbon, observed_proton)).  A pair (A,B) is geminal only if
  BOTH reciprocal correlations exist — HMBC near (Cb,Ha) and near (Ca,Hb) — which
  rejects the spurious one-directional links.  Only Leu/Val pair; Ile/Ala/Thr
  (single methyl) return None.  Returns {index -> partner index}."""
  def has(cq, hq):
    return any(abs(cp - cq) <= tol_c and abs(hp - hq) <= tol_h for cp, hp, _ in hmbc)
  part = {}
  for ai, (ca, ha, _) in enumerate(ilvat):
    found, bestd = None, 1e9
    for bi, (cb, hb, _) in enumerate(ilvat):
      if bi == ai or abs(cb - ca) < 0.3:
        continue
      if has(cb, ha) and has(ca, hb) and abs(cb - ca) < bestd:
        found, bestd = bi, abs(cb - ca)
    part[ai] = found
  return part


def validate(rows):
  """Compare picked+typed peaks to the known truth key by ppm match."""
  truth = []
  for l in (D / 'hmqc_true.tsv').read_text().splitlines()[1:]:
    f = l.split('\t')
    truth.append((float(f[1]), float(f[2]), f[3], f[4]))   # h, c, type, label
  picked = [(h, c, t) for h, c, t, _ in rows]
  used = [False] * len(picked)
  rec = typeok = 0
  for th, tc, tt, _ in truth:
    best = None
    for i, (h, c, t) in enumerate(picked):
      if used[i] or abs(c - tc) > 0.10 or abs(h - th) > 0.02:
        continue
      d = abs(c - tc) + abs(h - th)
      if best is None or d < best[0]:
        best = (d, i, t)
    if best:
      used[best[1]] = True
      rec += 1
      if best[2] == tt:
        typeok += 1
  print(f'validate: recovered {rec}/{len(truth)} true peaks; '
        f'type correct {typeok}/{rec}; spurious {sum(1 for u in used if not u)} picked-unmatched')


def main():
  kk = {'ILVAT': 65, 'ILV': 110, 'Val': 7, 'Thr': 10}
  # cap master list below the structural methyl count (89 in the trimer protomer
  # set) so the injective SAT stays feasible; extra picks are noise/tag.
  CAP = 88
  ilvat = pick(D / 'TNFa_ILVAT_13C_HMQC.ucsf', kk['ILVAT'])[:CAP]
  ilv = pick(D / 'TNFa_ILV_13C_HMQC.ucsf', kk['ILV'])
  val = pick(D / 'TNFa_Val_Methyl_HMQC.ucsf', kk['Val'])
  thr = pick(D / 'TNFa_Thr_Methyl_HMQC.ucsf', kk['Thr'])
  hmbc = pick(D / 'TNFa_HMBC_HMQC_13C_2D.ucsf', 8)
  part = geminal_partners(ilvat, hmbc)
  print(f'picked: ILVAT={len(ilvat)} ILV={len(ilv)} Val={len(val)} Thr={len(thr)} '
        f'HMBC={len(hmbc)} geminal-linked={sum(1 for v in part.values() if v is not None)}')

  # Ile delta1 is the only methyl with 13C below ~17 ppm (Ile 12.8-15.8;
  # every other type >= 18.3), so type it by chemical shift, not the Thr sample
  # (whose 13C window clips the low-delta1 Ile).
  # HMBC geminal link separates the paired types (Leu/Val) from the single-methyl
  # types (Ile/Ala/Thr): a peak with a geminal partner is L or V; without one it
  # is I/A/T.  Val is confirmed by the Val sample, and propagated across the
  # geminal link (if a peak's partner is Val, so is it).
  ILE_C_MAX = 17.0
  rows = []
  counts = {}
  for i, pk in enumerate(ilvat):
    c13 = pk[0]
    partner = part[i]
    is_val = matches(pk, val) or (partner is not None and matches(ilvat[partner], val))
    if c13 < ILE_C_MAX:
      t = 'I'
    elif is_val:
      t = 'V'
    elif partner is not None:  # geminal pair, not Val -> Leu
      t = 'L'
    elif matches(pk, thr):     # single methyl in Thr sample -> Thr
      t = 'T'
    elif matches(pk, ilv):     # in ILV, no partner picked -> Leu/Val, default Leu
      t = 'L'
    else:                      # single, not in ILV -> Ala
      t = 'A'
    counts[t] = counts.get(t, 0) + 1
    rows.append((pk[1], pk[0], t, pk[2]))   # h_ppm, c_ppm, res_type, height

  print('typed(raw):', counts, 'total', len(rows))

  # cap each type to its structural methyl capacity (injective SAT needs
  # peaks_of_type <= methyls_of_type); keep the strongest peaks per type.
  import maus
  lab = maus.parse_labeling('A;I;L;M;T;V')
  methyls = maus.parse_structure((D / 'fold_tnfa_trimer_model_0.cif').read_text().splitlines(), lab)
  cap = {}
  for m in methyls:
    cap[m.res_type] = cap.get(m.res_type, 0) + 1
  kept = []
  for t in set(r[2] for r in rows):
    same = sorted([r for r in rows if r[2] == t], key=lambda r: -r[3])   # height desc
    kept += same[:cap.get(t, 0)]
  print('capacity:', cap)
  print('typed(capped):', {t: sum(1 for r in kept if r[2] == t) for t in sorted(set(r[2] for r in kept))}, 'total', len(kept))
  validate(kept)
  rows = kept
  out = ['label\tH_ppm\tC_ppm\tres_type']
  for i, (h, c, t, _) in enumerate(sorted(rows, key=lambda r: (r[2], r[1])), 1):
    out.append(f'P{i}\t{h:.3f}\t{c:.3f}\t{t}')
  (D / 'hmqc_picked.tsv').write_text('\n'.join(out) + '\n')
  print(f'wrote {D}/hmqc_picked.tsv')
  write_noesy(D / 'TNFa_ILVAT_NOESY_HMQC.ucsf', top=80)


if __name__ == '__main__':
  main()
