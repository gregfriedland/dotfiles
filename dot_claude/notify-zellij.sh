#!/bin/bash
# Custom Claude Code notification hook with Zellij tab name.
# Usage: notify-zellij.sh <event_type>
#   event_type: stop | idle | permission

EVENT_TYPE="${1:-unknown}"

# Get Zellij tab name (empty string if not in Zellij).
TAB_NAME=""
if [ -n "$ZELLIJ" ]; then
    TAB_NAME=$(
        zellij action query-pane-info 2>/dev/null \
            | jq -r '.tab_name // empty' 2>/dev/null
    )
fi

# Build notification title and message based on event type.
TITLE="${TAB_NAME:-Claude Code}"

case "$EVENT_TYPE" in
    stop)
        MESSAGE="Task completed"
        ;;
    idle)
        MESSAGE="Waiting for input"
        ;;
    permission)
        MESSAGE="Needs permission"
        ;;
    *)
        MESSAGE="Notification"
        ;;
esac

# Send notification via terminal-notifier (falls back to
# osascript if unavailable).
if command -v terminal-notifier &>/dev/null; then
    terminal-notifier \
        -title "$TITLE" \
        -message "$MESSAGE" \
        -sound default \
        -group "claude-code-$$"
else
    osascript -e \
        "display notification \"$MESSAGE\" with title \"$TITLE\" sound name \"default\""
fi
