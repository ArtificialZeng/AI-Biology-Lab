# Source me: exports Antigravity agentapi env using the newest live daemon hub file (port changes on restart).
_H=$(ls -t ~/.gemini/antigravity/daemon/hub_*.json 2>/dev/null | head -1)
export ANTIGRAVITY_CSRF_TOKEN=$(python3 -c "import json;print(json.load(open('$_H'))['csrfToken'])")
export ANTIGRAVITY_LS_ADDRESS=127.0.0.1:$(python3 -c "import json;print(json.load(open('$_H'))['httpPort'])")
AGY=~/.gemini/antigravity/bin/agentapi
