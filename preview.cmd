@echo off
cd /d "%~dp0"
python -X utf8 build.py
if errorlevel 1 (
    echo Build failed. See the message above.
    pause
    exit /b 1
)
start "" "%~dp0docs\index.html"
