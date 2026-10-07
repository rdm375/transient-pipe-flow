# M7a transient-dynamics characterization

Run `make dynamics` from the repository root. The harness uses a
100-km pipe, N=40, dt=60 s, and a fine-time N=40 reference with
dt=1.5 s, theta=0.5. It compares all seven M6 theta values over
1800 s, with 100 -> 105 kg/s outlet demand changes.

Forcing IDs: 1 = instantaneous demand step; 2 = rectangular pulse
(300 <= t < 900 s); 3 = half-cosine ramp ending at 600 s.

`history.csv` stores sampled outlet pressure, inlet mass flow, and
linepack at 60-s intervals. `dynamics.csv` stores maxima of relative
errors over the whole trajectory, and extrema of outlet pressure.
The fine reference is **not** an exact solution. The 60-s output
sampling does not resolve sub-step oscillations or wavefronts.

The checks are deliberately loose regression envelopes, not
independent physical validation or evidence of unconditional stability.
The current reduced-momentum equations exclude convective acceleration.
This study does not yet establish propagation phase/amplitude accuracy.
