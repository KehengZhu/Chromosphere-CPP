#!/bin/sh
# Build the Doxygen code reference.
#
#   scripts/build_docs.sh            build, then report the warning count
#   scripts/build_docs.sh --open     also open the result in the default browser
#   scripts/build_docs.sh --strict   exit non-zero if Doxygen emitted any warning
#
# Output:  docs/doxygen/html/index.html   (gitignored build artifact; publish it
#          with scripts/publish_docs.sh)
# Warnings: docs/doxygen/doxygen-warnings.log
#
# Doxygen must run from the repository root — every path in the Doxyfile is
# relative to it.

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$repo_root"

open_after=0
strict=0
for arg in "$@"; do
    case "$arg" in
        --open)   open_after=1 ;;
        --strict) strict=1 ;;
        *) echo "build_docs: unknown option '$arg'" >&2; exit 2 ;;
    esac
done

if ! command -v doxygen >/dev/null 2>&1; then
    echo "build_docs: doxygen not found. Install it with:" >&2
    echo "    brew install doxygen        # macOS" >&2
    echo "    sudo apt-get install doxygen # Debian/Ubuntu" >&2
    exit 1
fi

echo "build_docs: doxygen $(doxygen --version)"

# Wipe the output tree first. Doxygen overwrites files but never deletes them,
# so anything it no longer generates lingers -- which is how disabling
# SOURCE_BROWSER once left 34 stale *_source.html listings in place, ready to be
# published. The tree is tracked in git, so a failed build is recoverable with
# `git checkout -- docs/doxygen/html`.
rm -rf docs/doxygen/html
doxygen docs/doxygen/Doxyfile

# Harmless belt-and-braces for any static host that runs Jekyll, which would
# drop the files and directories Doxygen names with a leading underscore.
touch docs/doxygen/html/.nojekyll

log="docs/doxygen/doxygen-warnings.log"
warnings=0
if [ -s "$log" ]; then
    warnings=$(grep -c . "$log" || true)
fi

echo "build_docs: wrote docs/doxygen/html/index.html"
if [ "$warnings" -gt 0 ]; then
    echo "build_docs: $warnings warning line(s) in $log"
    [ "$strict" -eq 1 ] && exit 1
else
    echo "build_docs: no warnings"
fi

if [ "$open_after" -eq 1 ]; then
    if command -v open >/dev/null 2>&1; then
        open docs/doxygen/html/index.html
    elif command -v xdg-open >/dev/null 2>&1; then
        xdg-open docs/doxygen/html/index.html
    fi
fi
