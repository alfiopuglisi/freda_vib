# freda_vib: VibSpec for a WebDAQ-504

Vibration spectrum analyser for a PC connected to an MCC **WebDAQ-504** with a triaxial IEPE accelerometer
(default `192.168.34.61`, channels 0/1/2 = x/y/z) and a single-axis IEPE accelerometer (channel 3,
`bench_z`), all 100 mV/g. It is a port of
[vibspec](https://github.com/alfiopuglisi/vibspec) (an Android app that does the same with the phone's
accelerometer) and has the same web page, analysis and HTTP API.

## Run

Double-click `vibspec.bat` (or the desktop shortcut). It starts the app and opens `http://localhost:8080`.
Close the window (or Ctrl+C) to quit: the WebDAQ job is stopped and its temporary log deleted.

`python vibspec.py --help` lists the options: `--host`, `--job` (default `freda`), `--channels` (triaxial,
default `0,1,2`), `--sensitivity` (triaxial, V/g, default 0.1 = 100 mV/g), `--extra` (single-axis channel,
default `3`, `none` to disable), `--extra-name`, `--extra-sensitivity` (default 0.1 V/g), `--vertical` (sensor
axis pointing up, default `z`), `--port`, `--bind 0.0.0.0` (let other computers open the page; no password),
`--rotate-min`, `--min-free-mb`, `--snapshots`. The sensitivities are written into the `freda` job on every
start: the WebDAQ converts volts to g with them, so a wrong value scales every amplitude.

Only one copy runs at a time: starting it again (e.g. double-clicking the shortcut twice) just opens the page
of the running one. If `freda` is already running on the WebDAQ (another copy on another PC, or a copy that
was killed), the app does not touch it until you press **Start** on its page.

Needs Python 3.8+ with numpy, requests and matplotlib (Anaconda has them).

## What it shows

- **Spectrum** of x, y, z and |xyz|, log or linear, in m/s², mg, mm/s or µm, as *sine peak amplitudes*.
- **Averaging** 0.1–10 s (Welch, Hann, 50% overlap), resolution 1 Hz by default (0.1–20 Hz settable),
  range default 1–100 Hz, sample rate 1280–10240 Hz.
- **Peaks** (up to 5, above 10× the median noise power): interpolated frequency, amplitude per axis, total,
  and the **3-D vibration direction** (principal axis of the band-limited covariance) with azimuth/elevation,
  angle from vertical and linearity (1 = straight line, 0.5 = planar, 0.33 = isotropic).
- **Single-axis channel** (`bench_z`): its spectrum is drawn as an extra line, with its band RMS and its
  amplitude at each peak; it is not used for the 3-D direction.
- **Snapshots** (with a start delay, so you can walk away): saved in `snapshots/vib_YYYYMMDD_HHMMSS/` as
  `_plot.png`, `_summary.txt`, `_spectrum.csv`, `_raw.csv`, downloadable from the page (also as zip).

Differences from the phone app: the IEPE inputs are AC coupled, so gravity is not measured; the vertical is
assumed along the sensor axis given by `--vertical`.

HTTP API (all GET): `/api/status`, `/api/data`, `/api/start`, `/api/stop`, `/api/hold?on=1`,
`/api/snapshot?delay=S`, `/api/cancel`, `/api/set?avg=&res=&fmin=&fmax=&unit=ACC|MG|VEL|DISP&log=1&peaks=&rate=&delay=&autostart=`,
`/api/snapshots`, `/files/<snapshot>/<file>`, `/files/<snapshot>.zip`.

## How it talks to the WebDAQ

- On start it creates/updates the job `freda` on the WebDAQ (ch 0–3, IEPE, chosen rate, logging into
  `/storage/data/freda/`) and runs it alone, like "Start job" in the WebDAQ page. Other saved jobs are not
  touched; if another job is running, the app refuses to start.
- Live samples (`/api/v1.0/schedule/jobs/freda/samples/<start>/<count>/bin`, float64 in g) are only served
  while the job logs to disk, and the WebDAQ storage is small (~290 MB free; 2560 Hz × 4 ch ≈ 300 MB/h).
  So the app restarts the job every 10 minutes (`--rotate-min`, or earlier below `--min-free-mb`) and
  deletes `/storage/data/freda/` on every start/stop. A restart drops a few seconds of data.

## Files

- `analyzer.py`: Welch spectrum, peaks, direction (numpy port of vibspec `Analyzer.kt`); `test_analyzer.py`
  checks it on synthetic data (`python test_analyzer.py`).
- `webdaq.py`: WebDAQ REST client.
- `vibspec.py`: acquisition + analysis threads, snapshots, HTTP server.
- `index.html`: the page (no external resources).
