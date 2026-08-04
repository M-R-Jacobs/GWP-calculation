Global warming potential (GWP) is a way of expressing how much a given mass of a gas contributes to warming relative to the same mass of CO2 over a chosen time horizon, 
so you can compare wildly different molecules on a single, reliable scale. It uses the radiative efficiency (RE), how strongly a molecule's IR absorption spectrum interacts 
with Earth's outgoing thermal radiation, in W m⁻² per ppb, combined with how long the molecule survives in the atmosphere before breaking down. A molecule can be a strong 
absorber but still have low GWP if it's destroyed quickly, and vice versa. A relevant distinguishing threshold that corrects a lifetime adjustment factor in the GWP has been 
identified in Hodnebrog et al. as ~10 years. 

CO2 is the reference gas specifically because it's the dominant anthropogenic driver and because its unique atmospheric decay behavior (the impulse response function, IRF)
is well characterized from carbon cycle models. It doesn't decay as one exponential as other chemicals might, but as a mix of processes on different timescales due to its
presence throughout earth's systems (ocean mixed layer, deep ocean uptake, terrestrial biosphere). Integrating its forcing over a single horizon isn't trivial, but the
function for the various processes (and cumulative IRF) exist and are not recalculated here - only utilized.

The global warming potential over a given time horizon, H is given by GWP(H) and is a ratio: the total radiative forcing produced by a 1 kg pulse of gas X, integrated from 
emission out to time horizon H, divided by the total radiative forcing a 1 kg pulse of CO2 would produce over that same H. The numerator of this term is the AGWPX(H) - the 
absolute global warming potential (also included is the ΔAGWPX(H) - the AR6 climate-carbon feedback term (from Eq. 7.SM.5.5), which accounts for additional CO2 released from 
natural sinks in response to the warming X causes.) It's "absolute" in the sense that it's not yet a ratio, just the raw integrated forcing in W m⁻²·yr per kg emitted. GWP(H) 
is then AGWPX(H) divided by the equivalent AGWPCO2(H). Both numerator and denominator are time-integrals of instantaneous forcing (W m⁻²·yr), with the standard horizon being 
H=100 years, though figures for H=20 and H=500 are also reported. Importantly, each side of this fraction is not the forcing itself - because X and CO2 decay at different 
rates and in different ways, you have to integrate their whole forcing trajectories over the time horizon, not just compare a snapshot like is done in the calculation of the RE.
