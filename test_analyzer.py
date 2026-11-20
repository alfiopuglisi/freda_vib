"""Synthetic-data checks of analyzer.py. Run: python test_analyzer.py"""
import math

import numpy as np

from analyzer import Params, analyze, from_acc

FS = 2560.0
rng = np.random.default_rng(1)


def signal(seconds, comps, noise=1e-5):
    t = np.arange(int(seconds * FS)) / FS
    xyz = rng.normal(0, noise, (t.size, 3))
    for f, a, d, ph in comps:
        d = np.asarray(d, float) / np.linalg.norm(d)
        xyz += np.outer(a * np.sin(2 * np.pi * f * t + ph), d)
    return xyz


def test_amplitude_frequency_direction():
    xyz = signal(10, [(50.3, 0.01, (1, 1, 0), 0.3), (12.0, 0.003, (0, 0, 1), 1.0)])
    r = analyze(xyz, FS, Params(avg_s=10, n_peaks=3))
    p1, p2 = r.peaks[0], r.peaks[1]
    assert abs(p1.f - 50.3) < 0.02, p1.f
    assert abs(p1.amp_dir - 0.01) / 0.01 < 0.02, p1.amp_dir
    assert abs(p1.amp_total - 0.01) / 0.01 < 0.02
    assert np.allclose(np.abs(p1.dir), [1 / math.sqrt(2), 1 / math.sqrt(2), 0], atol=0.01), p1.dir
    assert p1.linearity > 0.99 and abs(p1.angle_from_vertical - 90) < 1
    assert abs(p2.f - 12.0) < 0.02 and abs(p2.amp[2] - 0.003) / 0.003 < 0.02
    assert p2.dir[2] > 0.99 and p2.angle_from_vertical < 1   # sign points "up" (+z)
    # Spectrum value at the line ~ sine amplitude (Hann scalloping at most -1.4 dB)
    k = np.argmin(np.abs(r.freq - 50.3))
    assert 0.0085 < r.atot[k] < 0.0102, r.atot[k]


def test_circular_motion_linearity():
    xyz = signal(5, [(30, 0.01, (1, 0, 0), 0), (30, 0.01, (0, 1, 0), math.pi / 2)])
    p = analyze(xyz, FS, Params(avg_s=5)).peaks[0]
    assert abs(p.linearity - 0.5) < 0.02, p.linearity


def test_rms_and_units():
    xyz = signal(4, [(20, 0.02, (0, 0, 1), 0)], noise=0)
    r = analyze(xyz, FS, Params(avg_s=4, fmin=1, fmax=100))
    assert abs(r.rms_band[2] - 0.02 / math.sqrt(2)) / 0.0141 < 0.02, r.rms_band
    assert abs(from_acc(1.0, 0, 'MG') - 101.97) < 0.01
    assert abs(from_acc(2 * math.pi * 10, 10, 'VEL') - 1000) < 1e-9


def test_extra_channel():
    xyz = signal(5, [(30, 0.01, (1, 0, 0), 0)])
    t = np.arange(len(xyz)) / FS
    e = 0.004 * np.sin(2 * np.pi * 30 * t) + 0.002 * np.sin(2 * np.pi * 70 * t)
    r = analyze(np.column_stack([xyz, e]), FS, Params(avg_s=5))
    p = r.peaks[0]
    assert abs(p.f - 30) < 0.02 and abs(p.amp_dir - 0.01) / 0.01 < 0.02   # 3-D analysis unchanged
    assert abs(p.amp_extra - 0.004) / 0.004 < 0.02, p.amp_extra
    assert abs(r.rms_e - math.sqrt(0.004 ** 2 + 0.002 ** 2) / math.sqrt(2)) < 1e-4, r.rms_e
    k = np.argmin(np.abs(r.freq - 70))
    assert 0.0017 < r.ae[k] < 0.00205, r.ae[k]


def test_short_window_resolution():
    r = analyze(signal(0.5, [(40, 0.01, (1, 0, 0), 0)]), FS, Params(avg_s=0.5))
    assert abs(r.res_hz - 2.0) < 1e-9 and r.n_segments == 1
    assert abs(r.peaks[0].f - 40) < 0.2


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('ok', name)
