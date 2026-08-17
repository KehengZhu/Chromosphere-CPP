# Assets for the 2026-07-11 talk

Moved here out of the throwaway `docs/presentation-build/` directory when the talk was
filed under `docs/presentations/`, so its inputs are not sitting inside a build artifact.

`poster_c05_closed.jpg` and `poster_o00_open.jpg` are tracked: `\includegraphics` needs
them, so the deck would not build without them.

The two `event_20240801_*_Te.mp4` movies are **deliberately not tracked** — this repository
keeps movies out of Git (`visualization/` is gitignored for the same reason). They are only
`\href{run:...}` launch targets, so the PDF builds fine without them; the links simply do
nothing on a machine that lacks them. Regenerate them from the commands catalogued in
`visualize_commands.md`.
