@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0START_WINDOWS.ps1"
if errorlevel 1 pause
