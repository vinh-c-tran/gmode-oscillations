# g-mode Oscillations Hyperons

This repository consists of the code used to ultimately calculate the $g$-mode oscillations for compositions including nucleons and hyperons in the context of relativistic mean field (RMF) models. This is split up into three separate code segments that do three different things
1. Calculate the Equilibrium Parameters
2. Calculate the Adiabatic Sound Speed/Sound Speed Difference
3. Calculate Structure, g-mode oscillations and Tidal Deformability 


## Equilibrium Parameters
The first thing that we'd like to find is the equation of state: $p(\epsilon)$. In RMF models this reduces to determining the particle momentum $k_{F_i}$ and meson field values (we calculate these values as a function of baryon number density $n_B$) 
