#!/bin/sh
# Point this clone's Git hooks at the tracked .githooks/ directory.
#
#   scripts/install_hooks.sh
#
# Git does not share hooks through a clone, so every working copy has to opt in
# once. This sets core.hooksPath, which is a local config value — it affects
# only this clone and is never pushed.
#
# Installed hooks:
#   pre-commit   rebuilds docs/doxygen/html when a documentation input is
#                staged, fails the commit if Doxygen emits any warning, and
#                adds the refreshed reference to the commit.
#
# To turn the documentation hook off without uninstalling:
#   git config chromo.docs-hook false
# To uninstall entirely:
#   git config --unset core.hooksPath

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$repo_root"

git rev-parse --git-dir >/dev/null 2>&1 || {
    echo "install_hooks: not a Git repository: $repo_root" >&2
    exit 1
}

chmod +x .githooks/* 2>/dev/null || true
git config core.hooksPath .githooks

echo "install_hooks: core.hooksPath -> .githooks"
for hook in .githooks/*; do
    [ -f "$hook" ] && echo "install_hooks:   $(basename "$hook")"
done

command -v doxygen >/dev/null 2>&1 || {
    echo "install_hooks: NOTE - doxygen is not installed, so the pre-commit hook" >&2
    echo "install_hooks:        will block commits that touch documented sources." >&2
    echo "install_hooks:        Install it with 'brew install doxygen', or run" >&2
    echo "install_hooks:        'git config chromo.docs-hook false' to opt out." >&2
}
