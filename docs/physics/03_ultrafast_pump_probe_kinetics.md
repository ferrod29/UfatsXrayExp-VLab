# Ultrafast Pump-Probe Kinetics: IRF Convolution and the Spin-Crossover Cascade

Code: [`src/vlab_xray/kinetics.py`](../../src/vlab_xray/kinetics.py) ·
Notebooks: [`notebooks/06_static_and_transient_xas_modeling.ipynb`](../../notebooks/06_static_and_transient_xas_modeling.ipynb),
[`notebooks/05_xes_modeling.ipynb`](../../notebooks/05_xes_modeling.ipynb)

## 1. The pump-probe method

A short optical ("pump") laser pulse launches a photoreaction at time `t0`.
A short X-ray ("probe") pulse arrives at a variable delay `t` and records a
snapshot spectrum. Repeating this over many delays `t` builds up a movie of
the reaction. Because the pump pulse itself has a finite duration, the
measured "reaction start time" is smeared out: the population actually
excited at delay `t` is not a sharp step at `t = t0`, but that step
convolved with the pump pulse's temporal profile -- the **instrument
response function (IRF)**, modeled here as a Gaussian of width `sigma`
centered at `t0`.

## 2. From a bare rate equation to an IRF-convolved population

Ignoring the pulse duration for a moment, a state populated instantaneously
at `t0` and decaying with lifetime `tau` would simply follow
`N(t) = N0 * exp(-(t-t0)/tau)` for `t > t0` and `N(t) = 0` before. To account
for the finite pump duration, replace the instantaneous excitation at `t0`
with a continuous distribution of "micro-excitation" events spread out over
time according to the Gaussian pulse shape, and sum (integrate) the decay
that started at each of those times:

```
N_excited(t) = integral over x of
    Gaussian(x; 0, sigma) * Heaviside(t - t0 - x) * exp(-(t - t0 - x)/tau) dx
```

Here `x` is a dummy variable: "how much earlier than the pulse center did
this particular micro-excitation happen". This is
`kinetics._excited_state_kernel`, integrated by
`kinetics.excited_state_population` (`scipy.integrate.quad`, over
`x in [-10*sigma, +10*sigma]` -- wide enough that the Gaussian factor makes
anything further out numerically negligible, regardless of `t`, `t0`, or
`tau`). This is a standard building block in ultrafast spectroscopy, often
called an "exponentially modified Gaussian" (EMG) response.

The complementary quantity, the population that has **transferred out** of
that first excited state into the next state down the cascade, replaces the
decaying exponential with its complement:

```
N_product(t) = integral over x of
    Gaussian(x; 0, sigma) * Heaviside(t - t0 - x) * (1 - exp(-(t - t0 - x)/tau)) dx
```

(`kinetics._product_state_kernel` / `kinetics.product_state_population`).
By construction, `N_excited + N_product -> I0` for `t` many `tau` after
`t0` (everything has transferred), and both `-> 0` well before `t0` (nothing
has happened yet) -- see `tests/test_kinetics.py` for these sanity checks.

## 3. Two-state case: the transient XAS map

For the [Fe(terpy)<sub>2</sub>]<sup>2+</sup> XAS measurement, a simple
two-state (ground -> excited) model suffices: at pump-probe delay `t`, the
sample is a population-weighted mix of the pure ground- and excited-state
spectra,

```
spectrum(E, t) = ground_state(E) + N_excited(t) * (excited_state(E) - ground_state(E))
```

`kinetics.build_transient_xas_map` returns the pump-induced *change* -- the
second term alone,

```
transient(E, t) = N_excited(t) * (excited_state(E) - ground_state(E))
```

which is what a pump-probe measurement actually reports; add
`ground_state(E)` back to recover the full spectrum at a delay. It uses the
two static spectra fit in `xas_model` (see
[`01_xas_xanes.md`](01_xas_xanes.md)) as `ground_state`/`excited_state`.

## 4. Three-state cascade: the transient XES map

The XES experiment resolves an intermediate step: singlet ground state ->
triplet intermediate -> quintet product (see
[`02_xes_spin_states.md`](02_xes_spin_states.md)). Treating "triplet" as the
directly-pumped state and "quintet" as the state it transfers into:

```
spectrum(E, t) = N_triplet(t) * triplet(E)
               + N_quintet(t) * quintet(E)
               - (N_triplet(t) + N_quintet(t)) * singlet(E)
```

i.e. the singlet ground-state population is depleted by exactly the amount
that has moved into the triplet and quintet populations (particle-number
conservation). This is `kinetics.build_transient_xes_map`. The parameters
`I0, t0, sigma, tau` (amplitude, IRF center and width, and the
singlet->triplet->quintet transfer time -- in that order, as
`notebooks/05_xes_modeling.ipynb` unpacks them) are obtained from an
independent kinetics fit to the time-resolved data (not reproduced in this
project; the resulting numbers are hardcoded where they're used).

## 5. The quintet's own, much slower decay

Within the cascade window above, the quintet is treated as effectively
infinitely long-lived -- valid for delays up to about a nanosecond. On
longer timescales the quintet state itself relaxes back to the singlet
ground state (full spin-crossover recovery), with a lifetime `tau2` far
longer than the singlet->triplet->quintet transfer time (`660 ps` in
`notebooks/05_xes_modeling.ipynb`, vs. the ~100 fs IRF and ~sub-ps transfer
time). `kinetics.build_transient_xes_map_with_quintet_decay` therefore
evaluates the cascade normally up to a `cascade_window` after `t0`, then
**freezes** the transient spectrum at its value at that cutoff and lets the
whole thing decay back towards zero with the long lifetime:

```
cutoff = t0 + cascade_window
spectrum(E, t > cutoff) = spectrum(E, cutoff) * exp(-(t - cutoff)/tau2)
```

The slow decay is clocked from the `cutoff`, not from `t0`, so the two
branches meet exactly at the window edge whatever the ratio of
`cascade_window` to `tau2`.

This is a simplification (it assumes the *shape* of the transient spectrum
doesn't change during the quintet decay, only its overall amplitude, i.e.
that essentially only the quintet population is decaying back to the
singlet at this stage) but a good one, since by the end of the fast cascade
essentially all of the transferred population has settled into the quintet
state.
