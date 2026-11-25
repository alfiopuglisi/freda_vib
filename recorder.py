"""Continuous recording of raw samples and spectra to hourly CSV files."""
import io
import os
import shutil
import threading
import time
from datetime import datetime

import numpy as np


class Recorder:
    """Writes <root>/<YYYYMMDD>/raw_<YYYYMMDD_HHMMSS>.csv and spectrum_<...>.csv, a new pair every hour (and when the
    columns change, e.g. the single-axis channel is switched on/off). Stops itself when the disk gets full."""

    def __init__(self, root, min_free_gb):
        self.root = os.path.abspath(root)
        self.min_free = min_free_gb * 1e9
        self.on = False
        self.error = ''
        self.bytes = 0                  # written since the recording was switched on
        self.files = {}                 # kind -> [file, header, hour key, path]
        self.free = None
        self._last_check = 0.0
        self.lock = threading.Lock()

    def set_on(self, on):
        with self.lock:
            if on and not self.on:
                self.error, self.bytes = '', 0
                self._last_check = 0.0
            if not on:
                self._close_all()
            self.on = on
            return self._space_ok() if on else True

    def close(self):
        with self.lock:
            self._close_all()

    def current_files(self):
        with self.lock:
            return [v[3] for v in self.files.values()]

    def free_gb(self):
        try:
            return shutil.disk_usage(self.root if os.path.isdir(self.root) else os.path.dirname(self.root)).free / 1e9
        except OSError:
            return None

    # ---------- internals (call with the lock held) ----------
    def _close_all(self):
        for f, *_ in self.files.values():
            try:
                f.close()
            except OSError:
                pass
        self.files = {}

    def _space_ok(self):
        if time.time() - self._last_check > 30:
            self._last_check = time.time()
            try:
                os.makedirs(self.root, exist_ok=True)
                self.free = shutil.disk_usage(self.root).free
            except OSError as e:
                self._fail(f'cannot use {self.root}: {e}')
                return False
            if self.free < self.min_free:
                self._fail(f'recording stopped: less than {self.min_free / 1e9:g} GB free on the data disk')
                return False
        return True

    def _fail(self, msg):
        self.error = msg
        self.on = False
        self._close_all()
        print(time.strftime('%H:%M:%S'), msg)

    def _file(self, kind, header, t):
        """Open file for this kind, rotated every hour or when the header changes."""
        hour = time.strftime('%Y%m%d%H', time.localtime(t))
        cur = self.files.get(kind)
        if cur and cur[1] == header and cur[2] == hour:
            return cur[0]
        if cur:
            cur[0].close()
        dt = datetime.fromtimestamp(t)
        d = os.path.join(self.root, dt.strftime('%Y%m%d'))
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, f'{kind}_{dt:%Y%m%d_%H%M%S}.csv')
        f = open(path, 'a', encoding='ascii', newline='\n')
        if f.tell() == 0:
            f.write(header + '\n')
        self.files[kind] = [f, header, hour, path]
        return f

    def _write(self, kind, header, t, text):
        try:
            f = self._file(kind, header, t)
            f.write(text)
            f.flush()
            self.bytes += len(text)
        except OSError as e:
            self._fail(f'recording stopped: {e}')

    # ---------- public writers ----------
    def write_raw(self, t0, fs, data, names):
        """data: (n, k) samples in m/s², first sample at unix time t0."""
        with self.lock:
            if not self.on or len(data) == 0 or not self._space_ok():
                return
            t = t0 + np.arange(len(data)) / fs
            bio = io.StringIO()
            np.savetxt(bio, np.column_stack([t, data]), fmt=['%.5f'] + ['%.6g'] * data.shape[1], delimiter=',')
            header = 't_unix_s,' + ','.join(f'{n}_ms2' for n in names)
            self._write('raw', header, t0, bio.getvalue())

    def write_spectrum(self, t0, duration, freq, amps, names):
        """One spectrum (window starting at unix time t0): freq (m,), amps list of (m,) sine peak amplitudes, m/s²."""
        with self.lock:
            if not self.on or len(freq) == 0 or not self._space_ok():
                return
            m = len(freq)
            bio = io.StringIO()
            np.savetxt(bio, np.column_stack([np.full(m, t0), np.full(m, duration), freq, *amps]),
                       fmt=['%.3f', '%.3f', '%.4f'] + ['%.5g'] * len(amps), delimiter=',')
            header = 't_unix_s,window_s,f_Hz,' + ','.join(f'{n}_ms2' for n in names)
            self._write('spectrum', header, t0, bio.getvalue())
