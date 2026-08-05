#!/usr/bin/env bash
# Shared launcher registration -- the stop handle every campaign must write.
#
# Detached launchers (`setsid nohup bash launcher.sh ... &`) cannot be stopped
# via the `$!` captured at launch: that is the short-lived setsid PARENT, not
# the session-leader child actually running the script.  The only reliable
# handle is the launcher recording its own `$$` once it is running.
#
# Usage (one line near the top of every campaign launcher, after the runs dir
# is known):
#
#   source "${SCRIPT_DIR}/../lib/launcher_common.sh"
#   campaign_register_launcher "${RUNS_DIR}"
#
# `scripts/campaigns/stop_campaign.sh` reads this pid file first; its
# process-discovery fallbacks cover launchers that have not registered, but
# the pid file is the fast, unambiguous path.
#
# NOTE for wrappers that `exec` into a generic runner (run_fgeo.sh -> run.sh):
# `exec` preserves the PID, so registering BEFORE the exec records the correct
# final launcher pid.

campaign_register_launcher() {
    local runs_dir="$1"
    mkdir -p "${runs_dir}"
    echo $$ > "${runs_dir}/launcher.pid"
}
