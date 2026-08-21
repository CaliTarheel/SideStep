@echo off
rem Start the tectonic globe viewer. Usage: run-viewer.cmd [run_dir] [port]
setlocal
set RUN=%~1
if "%RUN%"=="" set RUN=out\hr640
set PORT=%~2
if "%PORT%"=="" set PORT=8077
cd /d "%~dp0"
if not exist target\release\viewer.exe (
  echo viewer.exe not built. Run: cargo build --release
  pause
  exit /b 1
)
echo Serving %RUN% at http://127.0.0.1:%PORT%/  (close this window to stop)
start "" http://127.0.0.1:%PORT%/
target\release\viewer.exe "%RUN%" %PORT%
