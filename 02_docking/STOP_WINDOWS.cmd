@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0PREPARE_MIGRATION.ps1"
if errorlevel 1 pause
