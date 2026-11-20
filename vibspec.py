"""VibSpec for a PC with a WebDAQ-504 + triaxial IEPE accelerometer.

Run:  python vibspec.py            (then the browser opens http://localhost:8080)
      python vibspec.py --help     for options (device address, job name, port, ...)
"""
import argparse
import dataclasses
import io
import json
import math
import os
import re
import signal
import threading
import time
import traceback
import webbrowser
import zipfile
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, unquote, urlparse

import numpy as np
import requests

import analyzer as A
from webdaq import RUNNING_CODES, WebDAQ

HERE = os.path.dirname(os.path.abspath(__file__))
RATES = [1280, 2560, 5120, 10240]
AXES = {'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1), '-x': (-1, 0, 0), '-y': (0, -1, 0), '-z': (0, 0, -1)}


def r4(v):
    v = float(v)
    return None if math.isnan(v) or math.isinf(v) else float(f'{v:.4g}')


class Engine:
    def __init__(self, args):
        self.args = args
        self.daq = WebDAQ(args.host, args.job)
        self.channels = tuple(int(c) for c in args.channels.split(','))          # triaxial x, y, z
        self.extra = None if args.extra in ('', 'none') else int(args.extra)       # single-axis channel
        self.extra_name = args.extra_name
        self.extra_on = True          # user setting: acquire and show the single-axis channel
        self.vertical = np.array(AXES[args.vertical], float)
        self.snap_root = os.path.abspath(args.snapshots)
        self.settings_file = os.path.join(HERE, 'settings.json')
        self.lock = threading.RLock()
        # user settings (persisted)
        self.params = A.Params()
        self.rate = 2560
        self.snapshot_delay = 5.0
        self.autostart = True
        self._load_settings()
        self._set_channels()
        # state
        self.want_run = self.autostart
        self.acquiring = False
        self.hold = False
        self.message = ''
        self.seq = 0
        self.latest = None
        self.last_snapshot = None
        self.device = f'WebDAQ @ {args.host}'
        self.free_bytes = None
        self.measured_rate = 0.0
        self.restart_req = False
        # sample buffer: rows = samples (m/s²); buf_start = absolute index of buf[0] in the current run
        self.buf = np.zeros((0, self.ncols))
        self.buf_start = 0
        self.run_started = None      # wall time of the current device run
        self.snap_target = None      # absolute sample index where the snapshot window starts
        self.snap_due_wall = None
        self.stop_evt = threading.Event()
        self.takeover = False
        self.closing = False

    def _set_channels(self):
        """Channels acquired in the next run: triaxial, plus the single-axis one if available and enabled."""
        names = self.args.names.split(',')
        spec = [(c, names[i], self.args.sensitivity) for i, c in enumerate(self.channels)]
        if self.extra is not None and self.extra_on:
            spec.append((self.extra, self.extra_name, self.args.extra_sensitivity))
        # the WebDAQ returns the channels sorted by number; column of each wanted channel in that data
        nums = sorted(c for c, _, _ in spec)
        self.chan_spec = spec
        self.cols = [nums.index(c) for c, _, _ in spec]
        self.ncols = len(spec)

    # ---------- settings ----------
    def _load_settings(self):
        try:
            with open(self.settings_file, encoding='utf-8') as f:
                s = json.load(f)
            self.params = A.Params(**s.get('params', {})).clamped()
            self.rate = int(s.get('rate', self.rate))
            self.snapshot_delay = float(s.get('delay', self.snapshot_delay))
            self.autostart = bool(s.get('autostart', self.autostart))
            self.extra_on = bool(s.get('extra', self.extra_on))
        except FileNotFoundError:
            pass
        except Exception as e:
            print('settings.json ignored:', e)

    def _save_settings(self):
        try:
            with open(self.settings_file, 'w', encoding='utf-8') as f:
                json.dump({'params': dataclasses.asdict(self.params), 'rate': self.rate,
                           'delay': self.snapshot_delay, 'autostart': self.autostart, 'extra': self.extra_on},
                          f, indent=1)
        except Exception as e:
            print('cannot save settings:', e)

    def apply_settings(self, q):
        with self.lock:
            p = self.params
            num = lambda k: float(q[k]) if k in q and q[k] != '' else None
            if num('avg') is not None: p = dataclasses.replace(p, avg_s=num('avg'))
            if num('res') is not None: p = dataclasses.replace(p, res_hz=num('res'))
            if num('fmin') is not None: p = dataclasses.replace(p, fmin=num('fmin'))
            if num('fmax') is not None: p = dataclasses.replace(p, fmax=num('fmax'))
            if q.get('unit', '').upper() in A.UNITS: p = dataclasses.replace(p, unit=q['unit'].upper())
            if 'log' in q: p = dataclasses.replace(p, log_y=q['log'] in ('1', 'true'))
            if num('peaks') is not None: p = dataclasses.replace(p, n_peaks=int(num('peaks')))
            self.params = p.clamped()
            if num('delay') is not None: self.snapshot_delay = min(max(num('delay'), 0.0), 3600.0)
            if 'autostart' in q: self.autostart = q['autostart'] in ('1', 'true')
            if num('rate') is not None and int(num('rate')) != self.rate:
                self.rate = int(min(max(num('rate'), 100), 51200))
                self.restart_req = True
            if 'extra' in q and (q['extra'] in ('1', 'true')) != self.extra_on:
                self.extra_on = q['extra'] in ('1', 'true')
                if self.acquiring:
                    self.restart_req = True   # the channel list of the WebDAQ job changes
                else:
                    self._set_channels()
                    self.buf = np.zeros((0, self.ncols))
            self._save_settings()
            self.seq += 1  # force the page to re-read even when held

    # ---------- commands ----------
    def start(self):
        with self.lock:
            self.want_run = True
            self.hold = False
            self.takeover = True   # an explicit Start may take over a running job of the same name

    def stop(self):
        with self.lock:
            self.want_run = False
            self.snap_target = None

    def set_hold(self, on):
        with self.lock:
            self.hold = on

    def request_snapshot(self, delay):
        with self.lock:
            self.want_run = True
            self.hold = False
            delay = max(0.0, float(delay))
            self.snap_target = self.buf_start + len(self.buf) + int(round(delay * self.rate))
            self.snap_due_wall = time.time() + delay
            self.message = ''

    def cancel_snapshot(self):
        with self.lock:
            self.snap_target = None

    # ---------- acquisition thread ----------
    def daq_loop(self):
        last_free = 0
        while not self.stop_evt.is_set():
            try:
                if self.want_run and not self.acquiring:
                    self._device_start()
                elif not self.want_run and self.acquiring:
                    self._device_stop()
                elif self.acquiring:
                    if self.restart_req or self._rotation_due():
                        self._device_stop(keep_message=True)
                        continue
                    self._poll()
                    if time.time() - last_free > 30:
                        self.free_bytes = self.daq.free_bytes()
                        last_free = time.time()
                time.sleep(0.25 if self.acquiring else 0.5)
            except Exception as e:
                self.message = f'WebDAQ error: {e}'
                print(traceback.format_exc())
                self.acquiring = False
                time.sleep(3)
        if self.acquiring:
            try:
                self._device_stop()
            except Exception as e:
                print('stop on exit failed:', e)

    def _rotation_due(self):
        """The WebDAQ keeps a log file while it streams; restart every N minutes so it can be deleted, or earlier if
        the device is running out of storage. Not while a snapshot window is pending, unless storage is critical."""
        critical = self.free_bytes is not None and self.free_bytes < self.args.min_free_mb * 1e6
        old = self.run_started and time.time() - self.run_started > self.args.rotate_min * 60
        return critical or (old and self.snap_target is None)

    def _device_start(self):
        st = self.daq.schedule_status()
        cur = st.get('currentJobName') or ''
        if cur and cur != self.daq.job:
            self.message = f'WebDAQ busy: job "{cur}" is running; stop it from the WebDAQ page first.'
            self.want_run = False
            return
        if cur == self.daq.job:
            # Either left over from a session that was killed, or another VibSpec (on this or another PC) is using
            # it. Never take it over silently: only after an explicit Start from this page.
            if not self.takeover:
                self.message = (f'Job "{cur}" is already running on the WebDAQ (another VibSpec, here or on another '
                                f'PC?). Press Start to take it over.')
                self.want_run = False
                return
            self.daq.stop()
            time.sleep(1)
        self.takeover = True   # from now on this instance owns the job (rotations, error recovery)
        info = self.daq.info()
        with self.lock:
            self._set_channels()
        chs = ', '.join(f'ch{c} {n} {s * 1000:g} mV/g' for c, n, s in self.chan_spec)
        self.device = f"{info.get('model', 'WebDAQ')} {info.get('name', '')} @ {self.args.host}, job {self.daq.job} ({chs})"
        self.daq.ensure_job(self.rate, self.chan_spec)
        self.daq.start()
        with self.lock:
            self.acquiring = True
            self.restart_req = False
            self.run_started = time.time()
            self.run_iter = None
            pending = None if self.snap_target is None else max(0, self.snap_target - (self.buf_start + len(self.buf)))
            self.buf = np.zeros((0, self.ncols))
            self.buf_start = 0
            self.fetched = 0
            if pending is not None:
                self.snap_target = pending
            self.message = ''
        self.free_bytes = self.daq.free_bytes()
        print(time.strftime('%H:%M:%S'), 'acquisition started at', self.rate, 'Hz')

    def _device_stop(self, keep_message=False):
        try:
            self.daq.stop()
            time.sleep(1.0)
        finally:
            self.acquiring = False
            try:
                self.daq.wipe_folder()
            except Exception as e:
                self.message = f'could not delete the WebDAQ log folder: {e}'
            self.free_bytes = self.daq.free_bytes()
        if not keep_message and self.message.startswith('WebDAQ error'):
            self.message = ''
        print(time.strftime('%H:%M:%S'), 'acquisition stopped')

    def _poll(self):
        st = self.daq.status()
        code = st.get('statusCode')
        if code not in RUNNING_CODES:
            self.acquiring = False
            self.message = f"WebDAQ job {st.get('status', code)}"
            return
        if self.run_iter is None:
            self.run_iter = st.get('iterationIndex')
        n = int(st.get('samplesAcquired', 0))
        if st.get('sampleRate'):
            self.measured_rate = float(st['sampleRate'])
        nch = int(st.get('channelCount', self.ncols))
        max_keep = int(self.rate * self.args.buffer_s)
        if n - self.fetched > max_keep:          # fell behind: skip ahead (and drop the buffer: discontinuity)
            self.fetched = n - max_keep
            with self.lock:
                self.buf = np.zeros((0, self.ncols))
                self.buf_start = self.fetched
        while self.fetched < n:
            cnt = min(n - self.fetched, int(self.rate * 2))
            a = self.daq.samples(self.fetched, cnt, nch)
            if len(a) == 0:
                break
            xyz = a[:, self.cols] * A.G0  # g -> m/s²
            with self.lock:
                self.buf = np.vstack([self.buf, xyz])
                self.fetched += len(a)
                drop = len(self.buf) - max_keep
                if drop > 0:
                    self.buf = self.buf[drop:]
                    self.buf_start += drop

    # ---------- analysis thread ----------
    def analysis_loop(self):
        while not self.stop_evt.is_set():
            try:
                self._analyze_once()
            except Exception:
                print(traceback.format_exc())
            time.sleep(0.4)

    def _analyze_once(self):
        with self.lock:
            p, buf, start, fs = self.params, self.buf, self.buf_start, float(self.rate)
            target = self.snap_target
        if len(buf) < 16:
            return
        end = start + len(buf)
        if target is not None:
            n_win = int(round(p.avg_s * fs))
            if target < start:   # window start no longer in memory (e.g. restart): take the oldest data
                target = start
            if end >= target + n_win:
                i0 = target - start
                r = A.analyze(buf[i0:i0 + n_win], fs, p, self.vertical, t0=target / fs)
                with self.lock:
                    self.snap_target = None
                if r is not None:
                    name = self.save_snapshot(r)
                    with self.lock:
                        self.latest, self.hold = r, True   # show what was saved
                        self.seq += 1
                        self.last_snapshot = name
                        self.message = f'snapshot {name} saved (held; press Live to resume)'
                return
        if self.hold or not self.acquiring:
            return
        n_win = min(len(buf), int(round(p.avg_s * fs)))
        r = A.analyze(buf[len(buf) - n_win:], fs, p, self.vertical, t0=(end - n_win) / fs)
        if r is not None:
            with self.lock:
                self.latest = r
                self.seq += 1

    # ---------- snapshots ----------
    def summary(self, r, title):
        p, u = r.params, r.params.unit
        lab = A.UNITS[u][0]
        g = lambda v: f'{v:.4g}'
        L = [title, f'Device: {self.device}',
             f'Sampling: {r.fs:.1f} Hz, {r.n_samples} samples, window {r.duration_s:.2f} s',
             f'Averaging: {r.n_segments} Hann segments of {r.seg_n / r.fs:.3f} s, 50% overlap, resolution {r.res_hz:.3f} Hz',
             f'Band {p.fmin:.1f}-{p.fmax:.1f} Hz, RMS x/y/z: {g(r.rms_band[0])} / {g(r.rms_band[1])} / {g(r.rms_band[2])} m/s²',
             *([f'Single-axis channel {self.extra_name} (ch{self.extra}): RMS {g(r.rms_e)} m/s² in the band'] if r.ae is not None else []),
             f'Axes: sensor frame (WebDAQ channels {",".join(map(str, self.channels))} = x, y, z); '
             f'vertical assumed along {self.args.vertical} (IEPE inputs are AC coupled: gravity is not measured).',
             f'Amplitudes: sine peak amplitude in {lab}.']
        for i, pk in enumerate(r.peaks):
            az, el = A.az_el(pk.dir)
            c = lambda v: g(float(A.from_acc(v, pk.f, u)))
            L += ['', f'Peak {i + 1}: {pk.f:.2f} Hz',
                  f'  amplitude x/y/z: {c(pk.amp[0])} / {c(pk.amp[1])} / {c(pk.amp[2])} {lab}',
                  f'  total {c(pk.amp_total)} {lab}, along main axis {c(pk.amp_dir)} {lab}',
                  f'  direction (x, y, z) = ({pk.dir[0]:.3f}, {pk.dir[1]:.3f}, {pk.dir[2]:.3f}), azimuth {az:.1f}°, elevation {el:.1f}°',
                  f'  angle from vertical {pk.angle_from_vertical:.1f}°, linearity {pk.linearity:.2f}']
            if r.ae is not None:
                L.append(f'  {self.extra_name}: {c(pk.amp_extra)} {lab}')
        if not r.peaks:
            L += ['', 'No peak above the noise floor.']
        return '\n'.join(L) + '\n'

    def raw_csv(self, r):
        """The samples r was computed from: time from the first sample (s), acceleration (m/s²) per channel."""
        bio = io.StringIO()
        bio.write('t_s,ax_ms2,ay_ms2,az_ms2' + (f',{self.extra_name}_ms2' if r.raw.shape[1] > 3 else '') + '\n')
        t = np.arange(len(r.raw)) / r.fs
        np.savetxt(bio, np.column_stack([t, r.raw]), fmt=['%.6f'] + ['%.6g'] * r.raw.shape[1], delimiter=',')
        return bio.getvalue()

    def latest_raw(self):
        """(filename, csv) of the samples behind the spectrum on screen, or None."""
        with self.lock:
            r, started = self.latest, self.run_started
        if r is None or r.raw is None:
            return None
        t_start = datetime.fromtimestamp((started or time.time()) + r.t0)   # approximate wall time of sample 0
        return f'vibspec_raw_{t_start:%Y%m%d_%H%M%S}.csv', self.raw_csv(r)

    def save_snapshot(self, r):
        name = 'vib_' + datetime.now().strftime('%Y%m%d_%H%M%S')
        d = os.path.join(self.snap_root, name)
        os.makedirs(d, exist_ok=True)
        u = r.params.unit
        with open(os.path.join(d, name + '_summary.txt'), 'w', encoding='utf-8') as f:
            f.write(self.summary(r, f'VibSpec snapshot {name}'))
        with open(os.path.join(d, name + '_spectrum.csv'), 'w', encoding='utf-8') as f:
            en = self.extra_name
            ex = f',{en}_{u},{en}_ms2,psd{en}_ms2^2/Hz' if r.ae is not None else ''
            f.write(f'f_Hz,x_{u},y_{u},z_{u},total_{u},x_ms2,y_ms2,z_ms2,total_ms2,psdx_ms2^2/Hz,psdy_ms2^2/Hz,psdz_ms2^2/Hz{ex}\n')
            f.write(f'# amplitudes are sine peak amplitudes ({A.UNITS[u][0]}); unit={A.UNITS[u][0]}\n')
            conv = [A.from_acc(a, r.freq, u) for a in (r.ax, r.ay, r.az, r.atot)]
            extra = [A.from_acc(r.ae, r.freq, u), r.ae, r.psd_e] if r.ae is not None else []
            cols = np.column_stack([r.freq, *conv, r.ax, r.ay, r.az, r.atot, *r.psd, *extra])
            np.savetxt(f, cols, fmt='%.6g', delimiter=',')
        with open(os.path.join(d, name + '_raw.csv'), 'w', encoding='utf-8') as f:
            f.write(self.raw_csv(r))
        try:
            self._plot_png(r, os.path.join(d, name + '_plot.png'), name)
        except Exception as e:
            print('PNG failed:', e)
        print('snapshot saved in', d)
        return name

    def _plot_png(self, r, path, title):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        p, u = r.params, r.params.unit
        fig = plt.figure(figsize=(14, 8), dpi=100)
        ax = fig.add_axes([0.06, 0.38, 0.6, 0.55])
        m = (r.freq >= max(p.fmin, r.freq[1])) & (r.freq <= min(p.fmax, r.fs / 2))
        lines = [(r.ax, '#e53935', 'x', '-'), (r.ay, '#43a047', 'y', '-'), (r.az, '#1e88e5', 'z', '-'), (r.atot, 'k', '|xyz|', '--')]
        if r.ae is not None:
            lines.append((r.ae, '#00acc1', self.extra_name, '-'))
        for a, col, lab, ls in lines:
            ax.plot(r.freq[m], A.from_acc(a[m], r.freq[m], u), ls, color=col, lw=1.0, label=lab)
        if p.log_y:
            ax.set_yscale('log')
        for i, pk in enumerate(r.peaks):
            ax.annotate(f'{pk.f:.1f}', (pk.f, float(A.from_acc(pk.amp_total, pk.f, u))), textcoords='offset points',
                        xytext=(0, 8), ha='center', fontsize=8, color=f'C{i + 1}')
        ax.set_xlim(p.fmin, min(p.fmax, r.fs / 2))
        ax.set_xlabel('Hz'); ax.set_ylabel(A.UNITS[u][0] + ' (sine peak)'); ax.grid(alpha=0.3, which='both')
        ax.legend(loc='upper right'); ax.set_title(title)
        a3 = fig.add_axes([0.68, 0.38, 0.3, 0.55], projection='3d')
        for i, (col, lab) in enumerate(zip(['#e53935', '#43a047', '#1e88e5'], 'xyz')):
            v = np.eye(3)[i]
            a3.quiver(0, 0, 0, *v, color=col, lw=1)
            a3.text(*(v * 1.1), lab, color=col)
        mx = max([pk.amp_dir for pk in r.peaks] + [1e-30])
        for i, pk in enumerate(r.peaks):
            l = 0.25 + 0.75 * pk.amp_dir / mx
            d = pk.dir * l
            a3.plot([-d[0], d[0]], [-d[1], d[1]], [-d[2], d[2]], lw=3, color=f'C{i + 1}', label=f'{pk.f:.1f} Hz')
        a3.set_xlim(-1, 1); a3.set_ylim(-1, 1); a3.set_zlim(-1, 1)
        a3.set_title('vibration directions (sensor frame)')
        if r.peaks:
            a3.legend(loc='lower left', fontsize=8)
        fig.text(0.02, 0.01, self.summary(r, '').strip(), family='monospace', fontsize=7.5, va='bottom')
        fig.savefig(path)
        plt.close(fig)

    def snapshots(self):
        if not os.path.isdir(self.snap_root):
            return []
        out = []
        for n in sorted(os.listdir(self.snap_root), reverse=True)[:100]:
            d = os.path.join(self.snap_root, n)
            if os.path.isdir(d) and n.startswith('vib_'):
                out.append({'name': n, 'files': [{'name': f, 'size': os.path.getsize(os.path.join(d, f))}
                                                 for f in sorted(os.listdir(d))]})
        return out

    # ---------- JSON for the page ----------
    def status(self):
        with self.lock:
            p = self.params
            snap_in = None if self.snap_target is None else (
                (self.snap_target - (self.buf_start + len(self.buf))) / self.rate + p.avg_s)
            return {
                'running': self.want_run and self.acquiring, 'starting': self.want_run and not self.acquiring,
                'closing': self.closing,
                'extraAvailable': self.extra is not None, 'extraOn': self.extra_on, 'extraName': self.extra_name,
                'hold': self.hold, 'rate': round(self.measured_rate or self.rate, 1), 'rateSet': self.rate,
                'rates': RATES, 'sensor': self.device, 'message': self.message,
                'freeMB': None if self.free_bytes is None else round(self.free_bytes / 1e6),
                'snapshotIn': snap_in, 'lastSnapshot': self.last_snapshot, 'snapDir': self.snap_root,
                'delay': self.snapshot_delay, 'autostart': self.autostart, 'vertical': self.args.vertical,
                'params': {'avg': p.avg_s, 'res': p.res_hz, 'fmin': p.fmin, 'fmax': p.fmax, 'unit': p.unit,
                           'unitLabel': A.UNITS[p.unit][0], 'log': p.log_y, 'peaks': p.n_peaks},
                'units': [{'id': k, 'label': v[0]} for k, v in A.UNITS.items()],
            }

    def data(self, since):
        with self.lock:
            o = {'status': self.status(), 'seq': self.seq}
            r = self.latest
            if r is None or since == self.seq:
                return o
        p, u = r.params, r.params.unit
        fmax = min(p.fmax, r.fs / 2)
        ks = np.nonzero((r.freq > 0) & (r.freq >= p.fmin) & (r.freq <= fmax))[0]
        f = r.freq[ks]
        arr = lambda a: [r4(v) for v in A.from_acc(a[ks], f, u)]
        o.update({'fs': r.fs, 'n': r.n_samples, 'segments': r.n_segments, 'resHz': r.res_hz, 'duration': r.duration_s,
                  'unit': A.UNITS[u][0], 'log': p.log_y, 'fmin': p.fmin, 'fmax': fmax,
                  'f': [r4(v) for v in f], 'x': arr(r.ax), 'y': arr(r.ay), 'z': arr(r.az), 't': arr(r.atot),
                  'vertical': [float(v) for v in r.vertical], 'rms': [r4(v) for v in r.rms_band]})
        if r.ae is not None:
            o.update({'e': arr(r.ae), 'extraName': self.extra_name, 'rmsE': r4(r.rms_e)})
        peaks = []
        for pk in r.peaks:
            az, el = A.az_el(pk.dir)
            c = lambda v: r4(A.from_acc(v, pk.f, u))
            peaks.append({'f': r4(pk.f), 'amp': [c(v) for v in pk.amp], 'total': c(pk.amp_total), 'main': c(pk.amp_dir),
                          'dir': [r4(v) for v in pk.dir], 'az': r4(az), 'el': r4(el),
                          'vert': r4(pk.angle_from_vertical), 'lin': r4(pk.linearity), 'extra': c(pk.amp_extra)})
        o['peaks'] = peaks
        return o


def make_handler(eng: Engine):
    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, ctype, body, extra=None):
            self.send_response(code)
            self.send_header('Content-Type', ctype)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            for k, v in (extra or {}).items():
                self.send_header(k, v)
            self.end_headers()
            self.wfile.write(body)

        def _json(self, o):
            self._send(200, 'application/json', json.dumps(o).encode())

        def do_GET(self):
            try:
                url = urlparse(self.path)
                path = unquote(url.path)
                q = {k: v[-1] for k, v in parse_qs(url.query).items()}
                if path in ('/', '/index.html'):
                    with open(os.path.join(HERE, 'index.html'), 'rb') as f:
                        return self._send(200, 'text/html; charset=utf-8', f.read())
                if path == '/api/status': return self._json(eng.status())
                if path == '/api/data': return self._json(eng.data(int(q.get('since', -1))))
                if path == '/api/start': eng.start(); return self._json(eng.status())
                if path == '/api/stop': eng.stop(); return self._json(eng.status())
                if path == '/api/hold':
                    eng.set_hold(q['on'] in ('1', 'true') if 'on' in q else not eng.hold)
                    return self._json(eng.status())
                if path == '/api/snapshot':
                    eng.request_snapshot(float(q.get('delay', eng.snapshot_delay))); return self._json(eng.status())
                if path == '/api/cancel': eng.cancel_snapshot(); return self._json(eng.status())
                if path == '/api/set': eng.apply_settings(q); return self._json(eng.status())
                if path == '/api/snapshots': return self._json(eng.snapshots())
                if path == '/api/raw.csv':
                    raw = eng.latest_raw()
                    if raw is None:
                        return self._send(404, 'text/plain', b'no data yet')
                    return self._send(200, 'text/csv', raw[1].encode(),
                                      {'Content-Disposition': f'attachment; filename="{raw[0]}"'})
                if path.startswith('/files/'):
                    return self._file(path[len('/files/'):])
                self._send(404, 'text/plain', b'not found')
            except Exception as e:
                traceback.print_exc()
                self._send(500, 'text/plain', f'Error: {e}'.encode())

        def _file(self, rel):
            ok = re.compile(r'^vib_\d{8}_\d{6}$')
            if rel.endswith('.zip') and ok.match(rel[:-4]):
                d = os.path.join(eng.snap_root, rel[:-4])
                if not os.path.isdir(d):
                    return self._send(404, 'text/plain', b'no snapshot')
                bio = io.BytesIO()
                with zipfile.ZipFile(bio, 'w', zipfile.ZIP_DEFLATED) as z:
                    for f in sorted(os.listdir(d)):
                        z.write(os.path.join(d, f), f'{rel[:-4]}/{f}')
                return self._send(200, 'application/zip', bio.getvalue(),
                                  {'Content-Disposition': f'attachment; filename="{rel}"'})
            parts = rel.split('/')
            if len(parts) == 2 and ok.match(parts[0]) and re.match(r'^[\w.\-]+$', parts[1]):
                fp = os.path.join(eng.snap_root, parts[0], parts[1])
                if os.path.isfile(fp):
                    with open(fp, 'rb') as f:
                        body = f.read()
                    png = fp.endswith('.png')
                    ctype = 'image/png' if png else 'text/plain; charset=utf-8' if fp.endswith('.txt') else 'text/csv'
                    disp = 'inline' if png or fp.endswith('.txt') else 'attachment'
                    return self._send(200, ctype, body, {'Content-Disposition': f'{disp}; filename="{parts[1]}"'})
            self._send(404, 'text/plain', b'not found')
    return H


def main():
    ap = argparse.ArgumentParser(description='Vibration spectrum analyser for a WebDAQ-504 + triaxial accelerometer')
    ap.add_argument('--host', default='192.168.34.61', help='WebDAQ address')
    ap.add_argument('--job', default='freda', help='WebDAQ job used for streaming (created if missing)')
    ap.add_argument('--channels', default='0,1,2', help='WebDAQ channels of the triaxial x,y,z')
    ap.add_argument('--names', default='board_x,board_y,board_z', help='names of the triaxial channels')
    ap.add_argument('--sensitivity', type=float, default=0.1, help='triaxial sensitivity, V/g (0.1 = 100 mV/g)')
    ap.add_argument('--extra', default='3', help='WebDAQ channel of the single-axis accelerometer ("none" to disable)')
    ap.add_argument('--extra-name', default='bench_z', help='name of the single-axis channel')
    ap.add_argument('--extra-sensitivity', type=float, default=0.1, help='single-axis sensitivity, V/g')
    ap.add_argument('--vertical', default='z', choices=list(AXES), help='sensor axis pointing up')
    ap.add_argument('--port', type=int, default=8080)
    ap.add_argument('--bind', default='127.0.0.1', help='0.0.0.0 to allow other computers to use the page')
    ap.add_argument('--snapshots', default=os.path.join(HERE, 'snapshots'), help='snapshot folder')
    ap.add_argument('--rotate-min', type=float, default=10, help='restart the WebDAQ job every N minutes to delete its log')
    ap.add_argument('--min-free-mb', type=float, default=40, help='restart earlier if the WebDAQ storage drops below this')
    ap.add_argument('--buffer-s', type=float, default=15, help='seconds of samples kept in memory')
    ap.add_argument('--no-browser', action='store_true')
    args = ap.parse_args()

    url = f'http://{"localhost" if args.bind in ("127.0.0.1", "0.0.0.0") else args.bind}:{args.port}'
    eng = Engine(args)
    # On Windows SO_REUSEADDR lets a second server bind a port that is in use: disable it, so that a second
    # copy of the app notices the first one instead of fighting over the WebDAQ job.
    ThreadingHTTPServer.allow_reuse_address = os.name != 'nt'
    srv, t_end = None, time.time() + 20
    while srv is None:
        try:
            srv = ThreadingHTTPServer((args.bind, args.port), make_handler(eng))
        except OSError:
            # Port in use: either another copy is running (open its page and quit), or a copy is still shutting
            # down after its window was closed (it says "closing", or no longer answers): wait for the port.
            other = None
            try:
                other = requests.get(url.replace('localhost', '127.0.0.1') + '/api/status', timeout=2).json()
            except Exception:
                pass
            if other is not None and not other.get('closing'):
                print(f'VibSpec is already running at {url}: opening it instead.')
                if not args.no_browser:
                    webbrowser.open(url)
                    time.sleep(3)
                return
            if time.time() > t_end:
                print(f'Port {args.port} is in use by another program: use --port to choose another one.')
                return
            print(f'waiting for port {args.port} (previous VibSpec still closing)...')
            time.sleep(1)
    srv.daemon_threads = True
    t_daq = threading.Thread(target=eng.daq_loop, name='daq', daemon=True)
    t_daq.start()
    threading.Thread(target=eng.analysis_loop, name='analysis', daemon=True).start()
    print(f'VibSpec: {url}  (WebDAQ {args.host}, job {args.job}). Ctrl+C to quit (stops the WebDAQ job).')
    if not args.no_browser:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()

    def shutdown(*_):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, shutdown)

    def stop_device_now():
        """Stop the WebDAQ job and delete its log (called when the console window is closed)."""
        eng.closing = True
        eng.want_run = False
        eng.stop_evt.set()
        t_daq.join(timeout=4.5)   # Windows allows ~5 s after a console close event

    if os.name == 'nt':
        import ctypes
        HANDLER = ctypes.WINFUNCTYPE(ctypes.c_bool, ctypes.c_uint)

        def on_console_event(ev):
            if ev in (2, 5, 6):   # CTRL_CLOSE_EVENT, CTRL_LOGOFF_EVENT, CTRL_SHUTDOWN_EVENT
                print('window closed: stopping the WebDAQ job...')
                stop_device_now()
                return True
            return False          # Ctrl+C / Ctrl+Break: default handling (KeyboardInterrupt)
        main.console_handler = HANDLER(on_console_event)   # keep a reference
        ctypes.windll.kernel32.SetConsoleCtrlHandler(main.console_handler, True)
    try:
        srv.serve_forever(poll_interval=0.5)
    except KeyboardInterrupt:
        pass
    finally:
        print('stopping...')
        eng.closing = True
        eng.want_run = False
        eng.stop_evt.set()
        t_daq.join(timeout=15)
        srv.server_close()


if __name__ == '__main__':
    main()
