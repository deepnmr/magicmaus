"""Minimal reader for UCSF/Sparky NMR spectrum files (.ucsf), numpy only.

Format: 180-byte file header, then 128 bytes per axis, then float32 (big-endian)
data stored tile by tile. ppm axis reconstructed from sf/sw/center.
"""
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Axis:
  nucleus: str
  npoints: int
  tile: int
  sf: float      # spectrometer freq MHz
  sw: float      # spectral width Hz
  center: float  # center ppm

  def ppm(self, idx):
    """Point index (may be fractional) -> ppm."""
    sw_ppm = self.sw / self.sf
    return self.center + sw_ppm * (0.5 - (np.asarray(idx) + 0.5) / self.npoints)


def read(path):
  """Return (data ndarray [axis0, axis1, ...], [Axis,...]).  data indexed in the
  file's axis order (axis 0 = first / w1)."""
  raw = Path(path).read_bytes()
  assert raw[:8] == b'UCSF NMR', 'not a UCSF file'
  ndim = raw[10]
  axes = []
  for a in range(ndim):
    off = 180 + 128 * a
    b = raw[off:off + 128]
    nucleus = b[:8].split(b'\x00')[0].decode('ascii', 'replace')
    npoints, = struct.unpack('>I', b[8:12])
    tile, = struct.unpack('>I', b[16:20])
    sf, sw, center = struct.unpack('>fff', b[20:32])
    axes.append(Axis(nucleus, npoints, tile, sf, sw, center))

  data_off = 180 + 128 * ndim
  buf = np.frombuffer(raw, dtype='>f4', offset=data_off)
  ntiles = [-(-ax.npoints // ax.tile) for ax in axes]          # ceil
  tiles = [ax.tile for ax in axes]
  padded = [nt * t for nt, t in zip(ntiles, tiles)]
  # file layout: tile-grid in C order; within a tile, C order over tile axes.
  # reshape to (nt0,nt1,..., t0,t1,...) then transpose to interleave.
  shp = ntiles + tiles
  arr = buf[:int(np.prod(padded))].reshape(shp)
  # move each tile axis next to its grid axis: (nt0,t0,nt1,t1,...)
  perm = []
  for i in range(ndim):
    perm += [i, ndim + i]
  arr = arr.transpose(perm).reshape(padded)
  # crop padding
  sl = tuple(slice(0, ax.npoints) for ax in axes)
  return np.ascontiguousarray(arr[sl]), axes


if __name__ == '__main__':
  import sys
  for p in sys.argv[1:]:
    data, axes = read(p)
    print(f'{Path(p).name}: shape={data.shape} '
          f'max={data.max():.3g} noise~{np.median(np.abs(data)):.3g}')
    for i, ax in enumerate(axes):
      lo, hi = ax.ppm(ax.npoints - 1), ax.ppm(0)
      print(f'  ax{i} {ax.nucleus:4s} n={ax.npoints} tile={ax.tile} '
            f'sf={ax.sf:.2f} sw={ax.sw:.1f} range=[{lo:.2f},{hi:.2f}] ppm')
