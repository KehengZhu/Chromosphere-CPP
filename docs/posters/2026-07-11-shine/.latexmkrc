# latexmk recipe for the SHINE poster.
# Output layout comes from ~/.latexmkrc: poster.pdf next to poster.tex
# (gitignored), every other build file in latex-build/.

# Gemini/UM beamerposter theme needs fontspec fonts from the TeX-Live tree, so
# the engine must be LuaLaTeX (XeLaTeX/fontconfig does not see them by name).
$pdf_mode = 4;   # 4 = lualatex
