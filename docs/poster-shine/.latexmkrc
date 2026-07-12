# latexmk recipe for the SHINE poster.
# Always keep build artifacts (.aux/.log/.fls/.pdf/...) out of the source tree:
# every latexmk run in this directory writes into latex-build/, even a bare
# `latexmk poster.tex` with no -output-directory flag.
$out_dir = 'latex-build';

# Gemini/UM beamerposter theme needs fontspec fonts from the TeX-Live tree, so
# the engine must be LuaLaTeX (XeLaTeX/fontconfig does not see them by name).
$pdf_mode = 4;   # 4 = lualatex
