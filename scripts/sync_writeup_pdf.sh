#!/bin/sh
# Mirror a compiled writeup PDF into the main repository's docs/ directory.
#
#   main.pdf   ->  docs/chromosphere_writeup.pdf   (deep engineering record)
#   paper.pdf  ->  docs/chromosphere_paper.pdf     (concise current release)
#
# The build itself lives in docs/writeup-overleaf/latex-build/, which is
# gitignored in BOTH repositories, so the copies made here are the only tracked
# PDFs of the writeup.
#
# Usage:
#   scripts/sync_writeup_pdf.sh                 sync whichever PDFs already exist
#   scripts/sync_writeup_pdf.sh <path/to.pdf>   sync one specific build product
#
# It is invoked automatically by the $success_cmd hook in
# docs/writeup-overleaf/.latexmkrc, so every successful latexmk run — from the
# VS Code LaTeX Workshop recipe or from the command line — refreshes the copy.
# That hook is a no-op anywhere the main repository is not present (e.g. when
# Overleaf compiles the same sources on its own servers).
#
# Copies are content-compared first, so an unchanged rebuild does not touch the
# destination mtime and does not create a spurious git diff.

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir="$repo_root/docs/writeup-overleaf/latex-build"
docs_dir="$repo_root/docs"

target_for() {
    case "$1" in
        main)  echo "chromosphere_writeup.pdf" ;;
        paper) echo "chromosphere_paper.pdf" ;;
        *)     echo "" ;;
    esac
}

sync_one() {
    src=$1
    base=$(basename -- "$src" .pdf)
    name=$(target_for "$base")

    if [ -z "$name" ]; then
        # Not one of the two documents we mirror (e.g. a figure-only subbuild).
        return 0
    fi
    if [ ! -f "$src" ]; then
        echo "sync_writeup_pdf: no such PDF: $src" >&2
        return 1
    fi
    if [ ! -d "$docs_dir" ]; then
        echo "sync_writeup_pdf: missing destination directory: $docs_dir" >&2
        return 1
    fi

    dst="$docs_dir/$name"
    if [ -f "$dst" ] && cmp -s -- "$src" "$dst"; then
        echo "sync_writeup_pdf: docs/$name already current"
        return 0
    fi
    cp -- "$src" "$dst"
    echo "sync_writeup_pdf: updated docs/$name from $base.pdf"
}

status=0
if [ "$#" -gt 0 ]; then
    for pdf in "$@"; do
        sync_one "$pdf" || status=1
    done
else
    found=0
    for base in main paper; do
        if [ -f "$build_dir/$base.pdf" ]; then
            found=1
            sync_one "$build_dir/$base.pdf" || status=1
        fi
    done
    if [ "$found" -eq 0 ]; then
        echo "sync_writeup_pdf: nothing to sync; no main.pdf or paper.pdf in $build_dir" >&2
        echo "sync_writeup_pdf: build the writeup first (latexmk in docs/writeup-overleaf)" >&2
        status=1
    fi
fi
exit "$status"
