#!/bin/bash
set -e

echo "[phenosuite] Creating path symlinks for hardcoded source() calls..."

# Satisfies: source('/srv/shiny-server/phenomenalist/...')
ln -sfn /srv/shiny-server/phenosuite /srv/shiny-server/phenomenalist

# Satisfies: source('/srv/shiny-server/Phenoptics-Menu/...')
ln -sfn /srv/shiny-server/phenosuite /srv/shiny-server/Phenoptics-Menu

# Configure reticulate Python
export RETICULATE_PYTHON="${RETICULATE_PYTHON:-/opt/venv/bin/python}"
echo "[phenosuite] RETICULATE_PYTHON=${RETICULATE_PYTHON}"

# Ensure shiny user can write to log and tmp directories
chown -R shiny:shiny /var/log/shiny-server 2>/dev/null || true
mkdir -p /tmp/shiny-server
chown -R shiny:shiny /tmp/shiny-server

echo "[phenosuite] Starting Shiny Server on port 3838..."
exec shiny-server
