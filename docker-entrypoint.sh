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

# Some legacy apps (pcf-v2, masquerade) call use_virtualenv('r-reticulate'),
# which reticulate resolves under ~/.virtualenvs/r-reticulate. Symlink that
# path to the real venv at /opt/venv for every user that might run an app.
for home in /home/shiny /root; do
  mkdir -p "${home}/.virtualenvs"
  ln -sfn /opt/venv "${home}/.virtualenvs/r-reticulate"
done
chown -R shiny:shiny /home/shiny/.virtualenvs 2>/dev/null || true

# Ensure shiny user can write to log and tmp directories
chown -R shiny:shiny /var/log/shiny-server 2>/dev/null || true
mkdir -p /tmp/shiny-server
chown -R shiny:shiny /tmp/shiny-server

echo "[phenosuite] Starting Shiny Server on port 3838..."
exec shiny-server
