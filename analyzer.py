"""Welch spectrum, peaks and 3-D vibration direction (port of vibspec core/Analyzer.kt to numpy).

Input: evenly sampled x, y, z acceleration in m/s². Amplitudes are sine peak amplitudes.
The WebDAQ IEPE inputs are AC coupled, so there is no gravity in the data: the vertical direction is
taken from configuration (sensor axis, default z) instead of being measured.
"""
import math
from dataclasses import dataclass, field, replace
from typing import List, Optional

import numpy as np

G0 = 9.80665

UNITS = {  # id: (label, order, scale)
    'ACC': ('m/s²', 0, 1.0),
    'MG': ('mg', 0, 1000 / G0),
    'VEL': ('mm/s', 1, 1000.0),
    'DISP': ('µm', 2, 1e6),
}


def from_acc(a, f, unit):
    """Convert acceleration amplitude(s) a (m/s²) at frequency f to the given unit."""
    _, order, scale = UNITS[unit]
    a = np.asarray(a, dtype=float)
    if order == 0:
        return a * scale
    f = np.asarray(f, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = a * scale / (2 * np.pi * f) ** order
    return np.where(f > 0, out, 0.0)


@dataclass
class Params:
    avg_s: float = 2.0       # analysed window = averaging time (0.1-10 s)
    res_hz: float = 0.0      # requested resolution, 0 = auto (1 Hz)
    fmin: float = 1.0
    fmax: float = 100.0
    unit: str = 'ACC'
    log_y: bool = True
    n_peaks: int = 3

    def clamped(self):
        fmin = min(max(self.fmin, 0.0), 10000.0)
        return replace(
            self,
            avg_s=min(max(self.avg_s, 0.1), 10.0),
            res_hz=0.0 if self.res_hz <= 0 else min(max(self.res_hz, 0.1), 20.0),
            fmin=fmin,
            fmax=max(min(max(self.fmax, 1.0), 10000.0), fmin + 1),
            unit=self.unit if self.unit in UNITS else 'ACC',
            n_peaks=int(min(max(self.n_peaks, 1), 5)),
        )

    def segment_s(self):
        return min(1 / self.res_hz if self.res_hz > 0 else 1.0, self.avg_s)


@dataclass
class Peak:
    f: float
    amp: np.ndarray          # per-axis sine peak amplitude, m/s²
    amp_total: float
    dir: np.ndarray          # unit vector of the principal axis, sensor frame
    amp_dir: float
    linearity: float
    angle_from_vertical: float
    amp_extra: float = float('nan')   # single-axis channel, same band, m/s²


@dataclass
class SpectrumResult:
    fs: float
    n_samples: int
    duration_s: float
    seg_n: int
    nfft: int
    n_segments: int
    res_hz: float
    freq: np.ndarray
    ax: np.ndarray
    ay: np.ndarray
    az: np.ndarray
    atot: np.ndarray
    psd: np.ndarray          # (3, nb) (m/s²)²/Hz
    vertical: np.ndarray     # unit vector "up" in the sensor frame (configured, not measured)
    rms_band: np.ndarray     # (3,) m/s² in [fmin, fmax]
    peaks: List[Peak]
    params: Params
    t0: float = 0.0          # time of the first analysed sample (s since acquisition start)
    raw: Optional[np.ndarray] = field(default=None, repr=False)  # (n, 3 or 4) m/s², kept for snapshots
    ae: Optional[np.ndarray] = None    # single-axis channel amplitude spectrum, m/s²
    psd_e: Optional[np.ndarray] = None
    rms_e: float = float('nan')


def _detrend(seg):
    """Remove mean and linear trend along the last axis."""
    n = seg.shape[-1]
    t = np.arange(n) - (n - 1) / 2.0
    stt = float(np.sum(t * t))
    mean = seg.mean(axis=-1, keepdims=True)
    slope = (np.sum(seg * t, axis=-1, keepdims=True) / stt) if stt > 0 else 0.0
    return seg - mean - slope * t


def analyze(data, fs, p0: Params, vertical=(0.0, 0.0, 1.0), t0=0.0) -> Optional[SpectrumResult]:
    """data: array (n, 3) = x, y, z in m/s², or (n, 4) with a single-axis channel in the last column."""
    p = p0.clamped()
    data = np.asarray(data, dtype=float)
    xyz = data[:, :3]
    n = xyz.shape[0]
    if n < 8:
        return None
    seg_n = min(n, max(8, int(round(p.segment_s() * fs))))
    nfft = (1 << (seg_n - 1).bit_length()) * 2   # zero-pad x2 for smoother curve / better interpolation
    hop = max(1, seg_n // 2)
    n_seg = 1 + (n - seg_n) // hop
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(seg_n) / seg_n)   # periodic Hann
    sum_w, sum_w2 = w.sum(), (w * w).sum()
    idx = np.arange(n_seg)[:, None] * hop + np.arange(seg_n)[None, :]
    segs = xyz.T[:, idx]                                  # (3, n_seg, seg_n)
    X = np.fft.rfft(_detrend(segs) * w, n=nfft, axis=-1)   # (3, n_seg, nb)
    nb = nfft // 2 + 1
    # Averaged cross-spectra (real parts): 0=xx 1=yy 2=zz 3=xy 4=xz 5=yz
    c = np.empty((6, nb))
    for i, (a, b) in enumerate([(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]):
        c[i] = np.mean((X[a] * np.conj(X[b])).real, axis=0)

    df = fs / nfft
    freq = np.arange(nb) * df
    one = np.full(nb, 2.0)
    one[0] = 1.0
    one[-1] = 1.0
    amp = np.sqrt(np.maximum(c[:3], 0)) / sum_w * one
    atot = np.sqrt(np.sum(amp ** 2, axis=0))
    psd = c[:3] / (fs * sum_w2) * one
    var_k = one / (nfft * sum_w2)

    k_min = max(1, int(math.ceil(p.fmin / df)))
    k_max = min(nb - 1, int(p.fmax / df))
    rms = np.sqrt(np.sum(c[:3, k_min:k_max + 1] * var_k[k_min:k_max + 1], axis=1)) if k_max >= k_min else np.zeros(3)
    up = np.asarray(vertical, dtype=float)
    up = up / np.linalg.norm(up)
    res = fs / seg_n
    peaks = _find_peaks(c, freq, var_k, k_min, k_max, res, df, up, p.n_peaks)
    r = SpectrumResult(fs, n, n / fs, seg_n, nfft, n_seg, res, freq, amp[0], amp[1], amp[2], atot,
                       psd, up, rms, peaks, p, t0, data)
    if data.shape[1] > 3:   # single-axis channel: same segments and window, not part of the 3-D analysis
        E = np.fft.rfft(_detrend(data[:, 3][idx]) * w, n=nfft, axis=-1)
        ce = np.mean(np.abs(E) ** 2, axis=0)
        r.ae = np.sqrt(ce) / sum_w * one
        r.psd_e = ce / (fs * sum_w2) * one
        r.rms_e = float(np.sqrt(np.sum(ce[k_min:k_max + 1] * var_k[k_min:k_max + 1]))) if k_max >= k_min else 0.0
        for pk in peaks:
            lo = max(1, int(round((pk.f - 2 * res) / df)))
            hi = min(nb - 1, int(round((pk.f + 2 * res) / df)))
            pk.amp_extra = math.sqrt(2 * float(np.sum(ce[lo:hi + 1] * var_k[lo:hi + 1])))
    return r


def _find_peaks(c, freq, var_k, k_min, k_max, res, df, up, n_peaks):
    if k_max - k_min < 2:
        return []
    tot = c[0] + c[1] + c[2]
    floor = np.sort(tot[k_min:k_max + 1])[(k_max - k_min + 1) // 2]
    lo_k, hi_k = max(k_min, 1), min(k_max, len(tot) - 2)
    ks = np.arange(lo_k, hi_k + 1)
    m = (tot[ks] > 10 * floor) & (tot[ks] >= tot[ks - 1]) & (tot[ks] > tot[ks + 1])
    cand = ks[m][np.argsort(-tot[ks[m]], kind='stable')]
    half_band = 2 * res   # Hann main lobe: +-2 resolution bins
    chosen = []
    for k in cand:
        if len(chosen) >= n_peaks:
            break
        if all(abs(freq[j] - freq[k]) > 2 * half_band for j in chosen):
            chosen.append(int(k))
    peaks = []
    for k in chosen:
        l, mm, r = (math.log(max(tot[k + d], 1e-300)) for d in (-1, 0, 1))
        den = l - 2 * mm + r
        delta = min(max(0.5 * (l - r) / den, -0.5), 0.5) if den < 0 else 0.0
        f = (k + delta) * df
        lo = max(1, int(round((f - half_band) / df)))
        hi = min(len(freq) - 1, int(round((f + half_band) / df)))
        s = slice(lo, hi + 1)
        v = var_k[s]
        cxx, cyy, czz, cxy, cxz, cyz = (float(np.sum(c[i, s] * v)) for i in range(6))
        cov = np.array([[cxx, cxy, cxz], [cxy, cyy, cyz], [cxz, cyz, czz]])
        vals, vecs = np.linalg.eigh(cov)   # ascending
        lam, d = vals[-1], vecs[:, -1].copy()
        trace = cxx + cyy + czz
        dot = float(d @ up)
        flip = dot < 0 if abs(dot) > 0.05 else d[np.argmax(np.abs(d))] < 0
        if flip:
            d = -d
        peaks.append(Peak(
            f=f,
            amp=np.sqrt(2 * np.diag(cov).clip(0)),
            amp_total=math.sqrt(2 * max(trace, 0)),
            dir=d,
            amp_dir=math.sqrt(2 * max(0.0, lam)),
            linearity=lam / trace if trace > 0 else 0.0,
            angle_from_vertical=math.degrees(math.acos(min(abs(dot), 1.0))),
        ))
    return peaks


def az_el(d):
    """Azimuth (deg, in the x-y plane from +x towards +y) and elevation (deg, towards +z)."""
    return math.degrees(math.atan2(d[1], d[0])), math.degrees(math.atan2(d[2], math.hypot(d[0], d[1])))
