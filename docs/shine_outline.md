problem:
awsom-r lacks dynamic chrompshere (not sure. please search literature)

model:
2 fluid 1.5D model + PFSS field lines
account for ionization&recombination
boundary conditions
heating and cooling

result:
ionization & recombination rates in chromosphere
chromospheric evaporation dynamics

Toward a Dynamic Chromosphere for AWSoM-R: A Field-Aligned Two-Fluid Model of Ionization and Chromospheric Evaporation

The Alfvén Wave Solar atmosphere Model (AWSoM) and its realtime variant (AWSoM-R), the SWMF's global model of the corona and solar wind, places its inner boundary in the upper chromosphere and injects Alfvén-wave Poynting flux through a thin, quasi-static chromospheric layer. As a result the partially ionized, radiatively cooled, dynamically evolving chromosphere — the reservoir that sets the mass and energy supplied to the corona — is prescribed rather than solved. We present a 1.5D field-aligned two-fluid (ion + neutral) chromosphere model, run along PFSS field lines, aimed at supplying this missing physics. The model evolves separate ion and neutral fluids coupled by collisional drag and two-temperature exchange, with field-aligned electron and neutral heat conduction, a physically grounded coronal upper boundary, and radiative cooling (optically-thick H I / Ca II / Mg II line losses plus optically-thin transition-region emission). Non-equilibrium hydrogen ionization and recombination — collisional, photoionization, and radiative/three-body channels — are integrated self-consistently. We report the resulting chromospheric ionization and recombination rates and their height structure, and demonstrate the model's dynamics through a flare-driven explosive chromospheric evaporation case, recovering the characteristic upflow and condensation-downflow signatures. The framework is designed as a drop-in dynamic lower boundary for AWSoM-R.

1. What your code needs vs. what the event gives you
Code input	What it is	Event source	How to convert
s, L, topology	field-aligned coord, loop length, open/closed	PFSS trace	pfsspy → trace field lines → arc length along 3-D path; closed if both ends hit r=R☉
B(s)	field strength along line	PFSS	sample `
footpoint (lat,lon)	where the line roots	PFSS + AR location	used to look up the ribbon mask and footpoint `
IC: T,n,V(s)	pre-flare atmosphere	model, NOT event	C7/FAL chromosphere anchored at footpoint + TR ramp + hydrostatic corona (make_loop_dat.py)
lower BC	deep-chromosphere wall	model	V=0 + numerical diffusion (already implemented)
upper/apex BC	apex or open top	PFSS topology	closed → reflecting apex (full-loop mode); open → coronal outflow BC
F_e, δ, E_c	nonthermal beam	STIX HXR spectrum	thick-target / warm-target fit → FLARE_FLUX / FLARE_DELTA / FLARE_E_CUT
t_on(s) per line	ignition timing	AIA 1600 ribbon arrival	ribbon reaches footpoint at t → stagger that line's beam on-time
which lines flare	spatial energy map	AIA 1600/1700 ribbon mask	footpoint-under-ribbon → beam on, weighted by brightness; else no beam
The key conceptual correction to "infer IC from event data": you don't. The event tells you the loop (length, shape, field, where/when it's hit). The pre-flare thermodynamic profile is a standard model — the chromosphere you care about is the C7/FAL stratification you already build, not something read off the magnetogram.