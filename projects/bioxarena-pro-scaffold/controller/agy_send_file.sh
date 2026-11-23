#!/bin/bash
# Usage: bash controller/agy_send_file.sh <conversation_id> <message_file>
# Sends the file's text as one message (lets non-ASCII messages go through an ASCII-only terminal command line).
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
source controller/agy_env.sh
$AGY send-message "$1" "$(cat "$2")" >/dev/null 2>&1 && echo "SENT to $1 ($(wc -c < "$2") bytes)" || echo "SEND FAILED to $1"
