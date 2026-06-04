# EPG_GRE Gradient Spoiling — Discussion Notes
**Date:** 2026-06-04  
**Repository:** gabrielaBelsley/EPG-X  
**Branch:** `claude/fervent-carson-pygwc`

---

## Context

Script `ModelSpoilCorrAllT1_CorrData.m` computes a spoiling correction factor:

```
Correction = STheory / SEPG
```

where:
- **STheory** = ideal perfectly-spoiled SS-SPGR signal (`ssSPGR.m`)
- **SEPG** = signal simulated with EPG accounting for incomplete RF spoiling (`EPG_GRE.m`)

The EPG signal is computed with RF phase cycling (φ₀ = 150°):

```matlab
phi = RF_phase_cycle(NTR, phi0);
[s0,~,~] = EPG_GRE(d2r(alphaArray(1,ialpha))*ones(NTR,1), phi, TR, T1(iT1), T2);
```

The question was how to additionally incorporate **gradient spoiling** into this simulation.

---

## Q1 — Does EPG already model gradient spoiling?

**Yes.** The EPG shift matrix `S` (built by `EPG_shift_matrices.m`) already represents the action of a spoiler gradient each TR:

```
S: F+k → F+(k+1),  F-k* → F-(k-1)*,  Zk → Zk  (unchanged)
```

Each application of `S` per TR pushes transverse coherences to higher k-space orders, exactly modelling the dephasing caused by the spoiler gradient. This is already the default behaviour in `EPG_GRE.m` via `SE = S * E`.

---

## Q2 — Why does gradient amplitude NOT affect the EPG signal (without diffusion)?

**Key result:** Without diffusion, the EPG signal is **completely independent of spoiler gradient amplitude**. This is a fundamental property of the EPG formalism.

### Proof

The EPG state at order `k` evolves as:

```
FF(kidx) = SE(kidx,kidx) * F(kidx,jj) + b(kidx)
```

The relaxation matrix `E` is diagonal with:
- `E2 = exp(-TR/T2)` for all transverse states (F+k and F-k, every k)  
- `E1 = exp(-TR/T1)` for all longitudinal states (Zk, every k)

**E is k-independent** — the same relaxation applies at every EPG order. Similarly, the RF matrix `T` is block-diagonal with identical 3×3 blocks at every k-order.

Therefore, applying `S^ngrad` instead of `S` only **relabels** k-space positions (F+0→F+ngrad instead of F+0→F+1) but does not change any amplitude or relaxation history. The signal at F+0 follows the same path regardless of ngrad.

### Empirical verification

Running `EPG_GRE` with `ngrad=1` and `ngrad=10` (200 RF pulses) gives **identical signals** — the difference is identically zero. This confirms the theoretical result.

### Physical interpretation

In SPGR the spoiler gradient dephases transverse magnetisation within the voxel. As long as diffusion is negligible, each isochromat evolves independently and the **ensemble-averaged signal** depends only on the phase distribution, not on the gradient amplitude that created it. With a uniform isochromat distribution (what EPG assumes), the signal is the same for any non-zero gradient.

---

## Q3 — How does gradient amplitude affect the EPG signal?

The **only** mechanism by which gradient amplitude changes the EPG signal is through **diffusion** (T2* effects aside, which EPG does not model).

When `diff` is provided:

```matlab
d.G   = [G_pre, G_read, G_spoil];  % gradient amplitudes, mT/m
d.tau = [tau_pre, tau_read, tau_spoil];  % durations, ms
d.D   = 1e-9;  % diffusion coefficient, m^2/s
```

The `E_diff` function computes b-values for each EPG order k:

```
b(k) ∝ k^2 * (γ * G * τ)^2
```

Higher gradient amplitude → larger b-values → stronger attenuation of high-order states → effectively better spoiling. This is the correct way to model gradient amplitude effects.

---

## Q4 — How are readout and spoiler gradients handled separately?

In SPGR, the readout gradient is **balanced**: the pre-winding gradient cancels the readout gradient from echo to echo, so the net shift per TR = 0. Only the **spoiler gradient** contributes to the EPG shift `S`.

- **Without diffusion**: the readout gradient is irrelevant to EPG — only the existence of a spoiler matters (not its amplitude).
- **With diffusion** (`diff` parameter): all gradient events (pre-winder, readout, spoiler) contribute b-values. The `diff.G` and `diff.tau` vectors must list all segments including zero-gradient periods so that `sum(tau) = TR`.

---

## Q5 — What does this mean for the correction factor STheory/SEPG?

**Without diffusion:**

The EPG correction factor `STheory/SEPG` depends **only on the RF phase cycling pattern** (φ₀ = 150°), TR, T1, T2, and flip angle. It is **independent of the spoiler gradient amplitude**.

The standard call:

```matlab
[s0,~,~] = EPG_GRE(d2r(alphaArray(1,ialpha))*ones(NTR,1), phi, TR, T1(iT1), T2);
```

already correctly models the incomplete spoiling due to RF phase cycling. No modification is needed for gradient spoiling when diffusion is negligible.

**With diffusion:**

If diffusion effects are important (large gradients, long TR, high D), use the `diff` parameter. For the actual sequence parameters (spoiler ~1.3 rad/voxel/TR, TR = 4.1 ms, in-vivo tissue), diffusion attenuation of EPG states will be small and the correction factor will be very close to the no-diffusion result.

---

## Comparison with Hargreaves epg_rfspoil.m

Hargreaves' `epg_rfspoil.m` uses `epg_grelax.m` with `kg=1` (the `kg` parameter feeds only diffusion formulas, never changes the number of gradient shifts). `epg_grad.m` applies exactly one shift per TR via `circshift`. This is mathematically identical to `EPG_GRE.m` with the default `S` shift.

Both codes produce the same correction factor for the same RF phase cycling parameters.

---

## Files modified

| File | Change |
|---|---|
| `EPGX-src/EPG_GRE.m` | Reverted to original Malik code — removed incorrect `ngrad` parameter |
| `EPG_GRE_GradientSpoiling_Discussion.md` | Updated to reflect correct understanding |

---

## Summary

| Question | Answer |
|---|---|
| Does EPG model gradient spoiling? | Yes — shift matrix `S` already does this |
| Does gradient amplitude change SEPG (no diffusion)? | **No** — EPG is scale-invariant in gradient amplitude without diffusion |
| How to model gradient amplitude effects? | Use the `diff` parameter (diffusion) |
| Is the original `EPG_GRE.m` call correct for this sequence? | **Yes** — no modification needed |
| What drives the correction factor? | RF phase cycling (φ₀) only, when diffusion is negligible |

---

## References

- Zur Y, Wood ML, Neuringer LJ. *Spoiling of transverse magnetization in steady-state sequences.* Magn Reson Med. 1991;21(2):251-263.
- Weigel M, Schwenk S, Kiselev VG, Scheffler K, Hennig J. *Extended phase graphs with anisotropic diffusion.* J Magn Reson. 2010;205(2):276-285.
- Yarnykh VL. *Actual flip-angle imaging in the pulsed steady state: A method for rapid three-dimensional mapping of the transmitted radiofrequency field.* Magn Reson Med. 2007;57(1):192-200.
- Preibisch C, Deichmann R. *Influence of RF spoiling on the stability and accuracy of T1 mapping based on spoiled FLASH with varying flip angles.* Magn Reson Med. 2009;61(1):125-135.
- Malik SJ, Teixeira RPAG, Hajnal JV. *Extended phase graph formalism for systems with magnetization transfer and exchange.* Magn Reson Med. 2018;80(2):767-779.
