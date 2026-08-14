#!/bin/sh
# Publish the Doxygen code reference to https://kehengphysics.site/docs/chromosphere/
#
#   scripts/publish_docs.sh              rsync the committed HTML as-is
#   scripts/publish_docs.sh --build      rebuild with Doxygen first (--strict)
#   scripts/publish_docs.sh --dry-run    list what would change, transfer nothing
#
# The remote tree is documented in /srv/www/kehengphysics.site/README.md on the
# server. It is served by /etc/nginx/snippets/kehengphysics-docs.conf and sits
# outside /var/www/wordpress, so WordPress never sees it.
#
# The reference is behind HTTP basic auth because Doxygen emits full source
# listings of this private repository.
#
# Override any of these in the environment:
#   CHROMO_DOCS_HOST    ssh destination            (default ubuntu@163.192.207.107)
#   CHROMO_DOCS_KEY     ssh identity file          (default ~/.ssh/id_ed25519)
#   CHROMO_DOCS_REMOTE  remote directory           (default /srv/www/.../docs/chromosphere)
#   CHROMO_DOCS_URL     URL printed on success

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$repo_root"

host=${CHROMO_DOCS_HOST:-ubuntu@163.192.207.107}
key=${CHROMO_DOCS_KEY:-$HOME/.ssh/id_ed25519}
remote=${CHROMO_DOCS_REMOTE:-/srv/www/kehengphysics.site/docs/chromosphere}
url=${CHROMO_DOCS_URL:-https://kehengphysics.site/docs/chromosphere/}

build=0
dry_run=0
for arg in "$@"; do
    case "$arg" in
        --build)   build=1 ;;
        --dry-run) dry_run=1 ;;
        -h|--help) sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "publish_docs: unknown option '$arg'" >&2; exit 2 ;;
    esac
done

[ -f "$key" ] || { echo "publish_docs: ssh key not found: $key" >&2; exit 1; }

src=docs/doxygen/html

# docs/doxygen/html/ is not tracked, so a fresh clone has nothing to publish.
# Build it automatically rather than failing on something we can just do.
if [ "$build" -eq 1 ]; then
    echo "publish_docs: rebuilding the reference"
    scripts/build_docs.sh --strict
elif [ ! -f "$src/index.html" ]; then
    echo "publish_docs: $src/ is absent (it is a build artifact, not tracked) -- building it"
    scripts/build_docs.sh --strict
fi

[ -f "$src/index.html" ] || {
    echo "publish_docs: $src/index.html still missing after build -- aborting" >&2
    exit 1
}

commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)
dirty=$(git status --porcelain -- ':!docs/doxygen/html' 2>/dev/null | head -1)

# DEPLOY-STAMP.txt is written on the server after the transfer, never into the
# tracked docs/doxygen/html/ tree -- a file that changes every deploy must not
# land in git. --exclude keeps --delete from removing it.
rsync_flags="-a --delete --exclude DEPLOY-STAMP.txt --human-readable --chmod=D755,F644"
[ "$dry_run" -eq 1 ] && rsync_flags="$rsync_flags --dry-run --itemize-changes"

echo "publish_docs: $src/  ->  $host:$remote/"
# shellcheck disable=SC2086
rsync $rsync_flags \
      -e "ssh -i $key -o StrictHostKeyChecking=accept-new" \
      "$src/" "$host:$remote/"

if [ "$dry_run" -eq 1 ]; then
    echo "publish_docs: dry run only, nothing transferred"
    exit 0
fi

ssh -i "$key" "$host" "cat > $remote/DEPLOY-STAMP.txt" <<EOF
Chromosphere2026 Doxygen reference -- GENERATED, do not hand-edit this directory.
deployed : $(date -u '+%Y-%m-%dT%H:%M:%SZ') UTC
commit   : $commit${dirty:+ (working tree had uncommitted changes)}
command  : scripts/publish_docs.sh
EOF

# Confirm the live site actually serves the new copy.
code=$(curl -sS -o /dev/null -w '%{http_code}' "$url" || echo "??")
echo "publish_docs: published commit $commit"
echo "publish_docs: $url -> HTTP $code"
[ "$code" = "200" ] || echo "publish_docs: WARNING expected HTTP 200 from the public docs site" >&2
