# AGU26 Abstract Draft

## Recommended title

**Toward a Dynamic Chromospheric Boundary for AWSoM: Ionization-Aware Thermodynamics and Conduction-Driven Evaporation**

## Abstract

The chromosphere mediates the exchange of mass and energy between the lower solar atmosphere and the corona, yet global coronal models such as the Alfvén Wave Solar atmosphere Model (AWSoM) represent this dynamic layer through a simplified lower boundary. We present a 1.5D field-aligned model that evolves the chromosphere’s time-dependent thermodynamic and hydrodynamic response and is designed for eventual coupling to AWSoM within the Space Weather Modeling Framework. The model evolves a gravitationally stratified hydrogen atmosphere with finite-volume hydrodynamics, field-aligned electron and neutral thermal conduction, radiative energy exchange, and physically consistent boundary conditions.

The thermodynamic closure consistently incorporates equilibrium hydrogen ionization into the pressure, conserved internal energy, and density- and temperature-dependent first adiabatic index Γ1(ρ,T). The ionization fraction is calculated from the Saha relation, and hydrogen ionization energy is included explicitly in the conserved internal energy. The resulting nonlinear relation between internal energy and temperature is used directly in the implicit treatment of thermal conduction. This closure is also applied consistently in hydrostatic initialization, state reconstruction, flux calculation, and boundary treatment.

In the hydrogen ionization region, ionization acts as an internal-energy buffer: part of the conductive energy changes the ionization state rather than immediately increasing the temperature, thereby modifying the pressure response and subsequent mass motion. The calculations develop a coherent evaporation flow whose upward mass flux is approximately uniform with height across much of the upper chromosphere. The model establishes a thermodynamically consistent platform for studying conduction-driven chromospheric evaporation and, in future coupled simulations, for supplying time-dependent mass and energy fluxes to the global corona and solar wind in AWSoM.


## Positioning

This version avoids numerical-development and debugging details. It presents the release as an integrated physical capability: an ionization-aware thermodynamic closure, consistently coupled to hydrodynamics, conduction, stratification, and boundaries. It does not claim a quantitatively converged evaporation solution or completed AWSoM coupling.
