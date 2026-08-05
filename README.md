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

Here is the collection of formulas referenced in the calculation file:

**Radiative efficiency:** cross-section convolved with the Pinnock weighting function:

$$
RE_{inst} = \int \sigma(\nu) \cdot F(\nu)\, d\nu \times 10^{15}
$$

**Lifetime correction:** accounts for non-uniform mixing when lifetime, τ, is short:

$$
f(\tau) = \frac{a\tau^b}{1 + c\tau^d}, \qquad RE_{final} = RE_{inst} \times f(\tau)
$$

**AGWP of X:** integrate the decaying forcing of a 1 kg pulse of X out to time horizon, H:

$$
AGWP_X(H) = RE_{final} \times \frac{M_{CO_2}}{M_X} \times \tau \left(1 - e^{-H/\tau}\right)
$$

**AGWP of CO2:** same operation for the reference gas, using its multi-exponential impulse response C(t):

$$
AGWP_{CO_2}(H) = RE_{CO_2} \int_0^H C(t)\, dt
$$
$$
C(t) = a_0 + a_1e^{-t/\tau_1} + a_2e^{-t/\tau_2} + a_3e^{-t/\tau_3}
$$

**Climate-carbon feedback term:** the extra forcing from CO2 that natural sinks release in response to X's warming:

$$
\Delta AGWP_X(H) = \gamma \int_0^H AGTP_X(t)\, r_F(H-t)\, dt
$$

**GWP:** combine numerator and denominator:

$$
GWP(H) = \frac{AGWP_X(H) + \Delta AGWP_X(H)}{AGWP_{CO_2}(H)}
$$

*Note for the CO2 impulse response function and climate-carbon feedback terms:*

- `a0 = 0.2173` — the fraction of a CO2 pulse that stays airborne essentially permanently (no decay term attached to it; represents ocean/land carbon cycle saturation on human timescales)
- `a1, a2, a3` — fractions removed by three different processes, each with its own timescale:
  - τ1 = 394.4 yr — slow deep-ocean equilibration
  - τ2 = 36.54 yr — intermediate ocean mixing
  - τ3 = 4.304 yr — fast surface ocean/biosphere uptake

All four `a` coefficients sum to 1 (they're fractional weights of total pulse mass). These are from Joos et al. (2013), adopted by AR6 as the standard CO2 IRF

- `γ` = 11.06×10¹² kgCO2 yr⁻¹ K⁻¹ — the climate-carbon feedback coefficient (Gasser et al. 2017b). Converts a temperature increase into an equivalent mass of CO2 released from natural sinks per year per degree of warming.
- `AGTPX(t)` — absolute global temperature potential: the temperature response at time t caused by a 1 kg pulse of X emitted at t=0 (computed via convolution of X's forcing decay with the two-layer climate response, this is the `AGTP_X` function using `lambda`, `a_f/a_s`, `tau_f/tau_s` in the script.)
- `rF(t)` — the impulse response of *carbon release* per unit warming, itself a 3-exponential sum (`alpha_C`/`tau_C` in the script.) Weights 0.6368/0.3322/0.0031 at timescales 2.376/30.14/490.1 yr) which represent multiple processes on different timescales (like the CO2 IRF, just for the feedback pathway instead of the original pulse.)
