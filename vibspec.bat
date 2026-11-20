@echo off
rem Starts VibSpec for the WebDAQ and opens it in the browser. Close this window (or Ctrl+C) to stop:
rem the WebDAQ job is stopped and its temporary log folder deleted.
cd /d "%~dp0"
set PY=python
if exist "%USERPROFILE%\anaconda3\python.exe" set PY="%USERPROFILE%\anaconda3\python.exe"
%PY% vibspec.py %*
pause
