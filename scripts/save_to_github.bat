@echo off
echo ==============================================
echo Saving your Capstone Project to GitHub...
echo ==============================================

cd /d "%~dp0"
git add .
git commit -m "Auto-commit: Update capstone project files"
git push origin main

echo ==============================================
echo Done! Please close this window.
echo ==============================================
pause
