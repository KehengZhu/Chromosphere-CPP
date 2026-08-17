# From observations to model inputs — the 2024‑08‑01 flare

This note explains the four observables you asked about — the **STIX hard‑X‑ray spectrum**, the **active‑region (AR) location**, the **AIA 1600/1700 Å ribbon mask**, and **ribbon arrival timing** — what each one physically *is*, and exactly how we turn it into a number our 1.5‑D field‑aligned model consumes. It is written around the concrete event we are running: the **M8.2 flare of 2024‑08‑01 ~07:09 UT in NOAA AR 13768**.

Supporting papers are in [`docs/supporting-papers/`](supporting-papers/): Krucker 2020 (STIX instrument), Kontar 2011 & 2018 (HXR → electron spectrum), Brown 1971 / Emslie 1978 (thick‑target), Kazachenko 2017 (ribbon reconnection flux), Liu & Qiu 2013 and Qiu 2025 (ribbon photometry → per‑loop heating), Fisher, Canfield & McClymont 1985 (explosive evaporation).

---

## 0. The mapping at a glance

Our code does **not** read the atmosphere off the magnetogram. The event supplies *geometry* (which loop) and *forcing* (the beam); the pre‑flare stratification is the standard C7 model atmosphere that [`util/make_loop_dat.py`](../../util/make_loop_dat.py) lays down. So every observable below feeds one of two things: **where/what loop** or **how it is hit**.

| Code input (knob) | Physical meaning | Observable it comes from | Instrument |
|---|---|---|---|
| loop length `L`, topology, footpoint `\|B\|`, 3‑D path | which flux tube | **AR location** → PFSS trace | ADAPT‑GONG magnetogram + HEK/NOAA AR position |
| `FLARE_DELTA` (δ) | electron spectral index | **HXR photon index** γ, via δ = γ+1 | STIX |
| `FLARE_E_CUT` (E_c) | low‑energy cutoff [keV] | **HXR spectral fit** (warm‑target) | STIX |
| `FLARE_BEAM_FLUX` (F_e) | beam energy flux [W m⁻²] | **HXR nonthermal power** ÷ footpoint area | STIX (+ AIA ribbon area) |
| which lines get a beam | is this footpoint flaring? | **AIA 1600 Å ribbon mask** | SDO/AIA |
| `FLARE_T_ON` per line | ignition time | **AIA 1600 Å ribbon arrival** at that footpoint | SDO/AIA |
| `FLARE_DUR` per line | injection duration | **AIA 1600 Å brightening duration** | SDO/AIA |

The code knobs are read in [`scenarios/pfss_field_line.cpp`](../../scenarios/pfss_field_line.cpp) (the `PFSS_FLARE` block) and applied by the thick‑target beam in [`physics.hpp`](../../physics.hpp) (`beam_heating_rate`).

---

## 1. STIX hard‑X‑ray spectrum → the electron beam `(F_e, δ, E_c)`

### What STIX measures

The Spectrometer/Telescope for Imaging X‑rays (STIX) on Solar Orbiter is a hard‑X‑ray imaging spectrometer covering **4–150 keV** (Krucker 2020). Flare‑accelerated electrons stream down the loop and emit **bremsstrahlung** (braking radiation) when they collide with the dense chromosphere; STIX records that X‑ray spectrum. A flare HXR spectrum has two parts:

- a **thermal** component at low energies (a steep, roughly exponential fall‑off) — the >10 MK plasma that is *already heated*;
- a **nonthermal power‑law tail** at higher energies, `I(ε) ∝ ε^(−γ)` — the signature of the accelerated electron *beam* we want.

The beam is exactly the energy source our `model_flare` represents, so STIX is the instrument that tells us how hard the beam is and how much power it carries.

### Turning the photon spectrum into electron‑beam parameters

This is a standard inversion (Brown 1971; reviewed in Kontar 2011). For a **collisional thick target** — electrons fully stopped in the chromosphere, which is our case — the injected electron flux spectrum is a power law `F(E) ∝ E^(−δ)`, and it produces a photon spectrum of index

> **δ = γ + 1**   (electron index = photon index + 1)

confirmed explicitly in Kontar 2018. The fit returns three things we need:

1. **Photon index γ** → electron index **δ = γ + 1**. This sets how the beam attenuates as it goes down the loop.
2. **Low‑energy cutoff E_c** — the energy below which the power law turns over. It dominates the *energetics* (most electrons are near E_c), but the standard "cold‑target" fit cannot pin it; the **warm‑target** method (Kontar 2018) constrains it to ~7%. Typical values are **15–40 keV** (Holman 2003 found 20–40 keV for an X‑class event).
3. **Total nonthermal power P** above E_c [erg s⁻¹] — the integral of the beam.

### From the fit to our three knobs

| STIX product | Code knob | Relation |
|---|---|---|
| photon index γ | `FLARE_DELTA` = δ | δ = γ + 1 |
| low‑energy cutoff E_c | `FLARE_E_CUT` | direct; sets the stopping column `N_c = 2×10²¹·E_c²` m⁻² in [`physics.hpp`](../../physics.hpp) |
| nonthermal power P, footpoint area A | `FLARE_BEAM_FLUX` = F_e | **F_e = P / A** [W m⁻²]; in the code `∫Q ds = F_e·g(t)` (thick target = all absorbed) |

E_c enters through the **stopping column** `N_c`: a higher cutoff means more energetic electrons, which punch deeper before stopping, so the deposition layer sits lower. In the code the beam weight is `Q ∝ n_tot·(1 + N/N_c)^(−δ/2)` where `N(s)` is the column from the apex — so δ and E_c together set *where* on the loop the energy lands, and F_e sets *how much*.

### This event (2024‑08‑01 M8.2)

Solar Orbiter was at Stonyhurst longitude **+159°** on this date, so the AR (Earth‑side W74°) sat **84.9° from SolO's central meridian — just on STIX's east limb**, i.e. STIX did observe it. A production run would fetch the STIX science data (`stixpy` → STIX data center) and forward‑fit a thermal + thick‑target model (the OSPEX‑equivalent step) to get γ, E_c, P directly.

For the present run we use **literature‑/GOES‑anchored values for an M8.2** (the geometry is the genuinely event‑specific part; the bespoke STIX fit is the next upgrade):

- **δ = 4** — a hard M8 spectrum (photon index γ ≈ 3, so δ = γ+1 = 4); AR 13768 is rooted in large sunspots, which tend to produce hard spectra (Saqri 2024).
- **E_c = 20 keV** — mid‑range of the Holman/Kontar values.
- **F_e = 3×10¹⁰ erg cm⁻² s⁻¹ = 3×10⁷ W m⁻²** ("3F10"): solidly above Fisher, Canfield & McClymont's (1985) explosive‑evaporation threshold (~10¹⁰ erg cm⁻² s⁻¹) and in the 10¹⁰–10¹¹ band typical of M‑class footpoints. `F_e = P/A`: e.g. a nonthermal power `P ≈ 3×10²⁸ erg s⁻¹` spread over a footpoint area `A ≈ 10¹⁸ cm²` gives ≈3×10¹⁰ erg cm⁻² s⁻¹.

Run command (what produced the current outputs):

```bash
PFSS_FLARE=1 FLARE_BEAM_FLUX=3e7 FLARE_DELTA=4 FLARE_E_CUT=20 \
  build/chromo_main outputs/event_20240801_AR13768_flare.txt \
  full ionization pfss_field_line scenarios/data/event_20240801_AR13768.dat 1.0 cooling
```

The result is physical and clamp‑free: peak loop temperature **7.4 MK**, explosive upflow **655 km s⁻¹**, apex density filling **22×**.

---

## 2. AR location → where to seed the loop

### What it is

The "AR location" is simply *where on the Sun the flare happened*: the active region's heliographic position. Flare catalogs (the Heliophysics Event Knowledgebase, HEK; or NOAA's daily region summary) report it in **Stonyhurst** coordinates — longitude measured from the Earth‑facing central meridian (+ = west), latitude from the equator.

For 2024‑08‑01 we queried HEK (via `sunpy`'s `Fido`) for flares ≥ M1 and found the one at our target time:

> **M8.2, 2024‑08‑01 07:09 UT, AR NOAA 13768, Stonyhurst (lon = +74° W, lat = −16° S).**

PFSS works in the rotating **Carrington** frame, so we convert with `astropy` (no hand arithmetic): at 07:09 UT the sub‑Earth point is Carrington L0 = 291.4°, giving the AR a Carrington position of **(≈5.3°, −16°)**.

### Why we need it and where it goes

PFSS reconstructs the potential coronal field from the photospheric magnetogram. Seeded *at the AR*, a field‑line trace returns exactly the geometry our 1‑D loop model needs — and nothing more (PFSS is current‑free: it gives connectivity, not the eruption or the non‑potential core):

- **loop length `L`** (arc length of the traced line),
- **topology** (closed loop vs. open),
- **footpoint `|B|`** (field strength at the photospheric end),
- the **3‑D path** (for the eventual 3‑D movie).

That is what [`util/extract_event_loop.py`](../../util/extract_event_loop.py) does: it loads the ADAPT‑GONG map, reprojects it to the sine‑latitude (CEA) grid PFSS expects, converts the AR Stonyhurst position to Carrington, solves PFSS, traces a small cluster of seeds around the AR, and picks a representative **closed** loop. For this event it found **128.7 Mm total length (64.3 Mm half), footpoint |B| = 643 G**. Those two numbers feed straight into the loop builder:

```bash
python util/make_loop_dat.py --topology full \
  --half-length-mm 64.3 --fp-B-G 643 --apex-T-MK 2.0 \
  --output scenarios/data/event_20240801_AR13768.dat
```

| PFSS output | Code knob |
|---|---|
| half‑length L [Mm] | `make_loop_dat.py --half-length-mm` |
| footpoint \|B\| [G] | `make_loop_dat.py --fp-B-G` |
| topology | `make_loop_dat.py --topology` (→ `[META] topology=` → reflecting BC in `pfss_ic`) |

---

## 3. AIA 1600/1700 Å ribbon mask → *which* field lines flare

### What flare ribbons are

When the beam (and conductive flux) slam into the chromosphere, the footpoints light up in the ultraviolet. Imaged by SDO/AIA, these bright footpoints form **flare ribbons** — the chromospheric "scorch marks" of the loops that have just reconnected (Fletcher 2004; Kazachenko 2017). They are the single best map of *where* flare energy is being deposited, and therefore of *which* loops are flaring at all. This is the answer to "not every field line has evaporation": **a line flares only if its footpoint lies under a ribbon.**

### 1600 vs 1700 Å — why both

- **1600 Å** contains the C IV lines (formed ~10⁵ K, transition region) plus UV continuum — strongly enhanced during flares.
- **1700 Å** is essentially pure photospheric continuum — *not* flare‑sensitive.

So 1700 Å is the quiescent baseline and 1600 Å is the flare signal; the **1600/1700 ratio** (or 1600 relative to its own pre‑flare level) cleanly isolates the genuine flare brightening from ordinary photospheric/network emission (Qiu 2025 shows 1600 Å relative brightness is a quantitative chromospheric‑heating measure).

### Building the mask and using it

1. For each AIA frame, threshold the 1600 Å brightness (e.g. > 10× the quiescent level, after the 1700 Å baseline) → a **binary ribbon mask** of flaring pixels at that time.
2. Trace the PFSS ensemble; project each line's footpoint onto the mask. **Footpoint inside the mask → beam on; outside → no beam** (quiet line, conductive response only).
3. **Weighting:** the ribbon also sets *how much* energy each line gets. Kazachenko 2017 defines the **reconnection flux** Φ = ∫ B_n dA over the area the ribbon has swept, and finds the GOES peak X‑ray flux scales as `I_X ∝ Φ^1.5`. The local 1600 Å brightness (or local reconnection rate) gives the **relative** `FLARE_BEAM_FLUX` per line; the ensemble is then normalized so the total nonthermal power matches the STIX fit from §1.

| Ribbon product | Code use |
|---|---|
| binary 1600 Å mask at time t | per‑line beam on/off |
| local 1600 Å brightness / reconnection rate | relative `FLARE_BEAM_FLUX` weight (normalized to STIX total) |

---

## 4. AIA 1600 Å ribbon arrival → *when* each line ignites (`FLARE_T_ON`)

A flare is not lit all at once: the ribbons **sweep across the arcade** as reconnection progressively involves new field lines. A given footpoint stays dark until the ribbon *arrives* there, then brightens sharply. In the 1600 Å light curve of a single pixel this shows as an impulsive rise — half‑rise time < 2 min (Qiu 2025) — and that rise marks the instant energy deposition begins on the loop rooted at that pixel.

So the **time the 1600 Å ribbon reaches a footpoint is that line's `FLARE_T_ON`**, and the **duration of the brightening sets `FLARE_DUR`**. Staggering ignition this way across the ensemble reproduces the observed behaviour where each newly‑lit footpoint runs through the *same* "elementary flare kernel" evolution — ~300 km s⁻¹ evaporation upflow, ~40 km s⁻¹ condensation downflow — just offset in time (Graham & Cauzzi 2015). That is exactly the per‑footpoint signature our `model_flare` produces.

| Ribbon product | Code knob |
|---|---|
| 1600 Å arrival time at footpoint | per‑line `FLARE_T_ON` |
| 1600 Å brightening duration | per‑line `FLARE_DUR` |

---

## 5. Putting it together — the per‑line driver

The full event pipeline (Phase 2 of the campaign) is then:

```text
for each PFSS field line:
    footpoint ← line's photospheric root           # §2
    if footpoint ∈ AIA 1600 ribbon mask(t):         # §3  (else: quiet, no beam)
        FLARE_T_ON   ← ribbon arrival time          # §4
        FLARE_DUR    ← ribbon brightening duration  # §4
        FLARE_BEAM_FLUX ← local weight × F_total     # §3, normalized to STIX P  (§1)
        FLARE_DELTA  ← γ + 1     (STIX)              # §1
        FLARE_E_CUT  ← E_c       (STIX)              # §1
        run model_flare on this loop (geometry from §2)
```

**What is done now (the starter, single line):** the geometry of one representative AR‑13768 loop is real (§2); the beam values are M8.2‑appropriate STIX‑style numbers (§1); timing is uniform (a single `T_ON`/`DUR`, not yet ribbon‑staggered). **What remains:** the AIA ribbon mask and per‑footpoint timing (§3–§4) to turn the single line into an ensemble where the flare spreads loop‑by‑loop, then map the 1‑D solutions onto the 3‑D PFSS paths for the movie.

---

## 6. Data products & how to fetch — quick reference

| Product | Instrument / access | Gives | Status |
|---|---|---|---|
| ADAPT‑GONG global map | `gong.nso.edu/adapt/maps/gong/` (HTTP) | PFSS lower boundary | **fetched** (06:00 UT map) |
| AR position / flare list | HEK via `sunpy` `Fido` (`a.hek.FL`) | AR location, GOES class, peak time | **done** (M8.2, AR 13768) |
| GOES XRS | `Fido` (`a.Instrument('XRS')`) | flare class, impulsive‑phase timing → energy budget | to fetch |
| STIX HXR | `stixpy` → STIX data center | γ, E_c, P (the beam) | **method documented (§1)**; values literature‑anchored pending fit |
| AIA 1600 / 1700 Å | `Fido` (`a.Instrument('AIA')`, `a.Wavelength`) | ribbon mask, arrival timing | Phase 2 |
| AIA 131 / 304 Å | `Fido` | flare‑loop emission → validation of the synthetic movie | Phase 2 |

> Practical note: solar‑data servers are reachable from this environment, but their Apache autoindex pages can be large — write `curl` output to a file rather than piping to a terminal. `sunpy`'s networking needs `beautifulsoup4`, `lxml`, `zeep`, and `drms` (now installed in `.venv`).

---

### References (in `docs/supporting-papers/`)

- **Krucker et al. 2020**, *A&A* 642, A15 — STIX instrument (4–150 keV).
- **Brown 1971**, *Solar Phys.* 18, 489 — thick‑target bremsstrahlung.
- **Emslie 1978**, *ApJ* 224, 241 — collisional beam stopping (our `(1+N/N_c)^(−δ/2)`).
- **Kontar et al. 2011**, *Space Sci. Rev.* 159, 301 — deducing electron spectra from HXR.
- **Kontar et al. 2018**, *ApJ* 868, 109 — warm‑target low‑energy cutoff (δ = γ+1, E_c to ~7%).
- **Kazachenko et al. 2017**, *ApJ* 845, 49 — RibbonDB: reconnection flux from AIA 1600 Å, `I_X ∝ Φ^1.5`.
- **Liu, Qiu et al. 2013**, *ApJ* 770, 111 — per‑loop heating functions from UV 1600 + thick‑target HXR.
- **Fisher, Canfield & McClymont 1985**, *ApJ* 289, 414/425/434 — explosive chromospheric evaporation (the regime `model_flare` reproduces).
