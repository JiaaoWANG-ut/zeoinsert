#!/usr/bin/env bash
# Build manuscript.pdf from manuscript.md using Pandoc + Tectonic.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FINAL="$(cd "$(dirname "$0")" && pwd)"
VENV_PY="$ROOT/.venv/bin/python"
TECTONIC="$ROOT/.bin/tectonic"

if [[ ! -x "$VENV_PY" ]]; then
  echo "Creating venv and installing pypandoc-binary..."
  cd "$ROOT"
  uv venv .venv
  uv pip install pypandoc-binary
fi

if [[ ! -x "$TECTONIC" ]]; then
  echo "Downloading Tectonic..."
  mkdir -p "$ROOT/.bin"
  curl -fsSL "https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.15.0/tectonic-0.15.0-x86_64-unknown-linux-gnu.tar.gz" \
    | tar -xz -C "$ROOT/.bin"
fi

cd "$FINAL"
"$VENV_PY" - <<'PY'
import pypandoc

extra_args = [
    "--pdf-engine=../.bin/tectonic",
    "--include-in-header=pandoc-header.tex",
    "-V", "geometry:margin=1in",
    "-V", "documentclass=article",
    "-V", "fontsize=11pt",
    "--resource-path=.:figures",
    "--standalone",
    "--toc",
    "--number-sections",
]
pypandoc.convert_file(
    "manuscript.md", "pdf", outputfile="manuscript.pdf", extra_args=extra_args
)
print("Wrote final/manuscript.pdf")
PY

ls -lh "$FINAL/manuscript.pdf"
