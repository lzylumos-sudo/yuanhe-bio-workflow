@echo off
rem Double-click on Windows. UNTESTED on a real Windows machine in this release.
chcp 65001 >nul
set "HERE=%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%HERE%run.ps1" web
if errorlevel 1 pause
