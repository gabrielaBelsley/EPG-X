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

The relaxation matrix `E` is diagonal with:
- `E2 = exp(-TR/T2)` for all transverse states (F+k and F-k*, every k)
- `E1 = exp(-TR/T1)` for all longitudinal states (Zk, every k)

**E is k-independent** — the same relaxation applies at every EPG order. The RF matrix `T` is also block-diagonal with **identical** 3×3 blocks at every k-order (the same flip angle and phase is applied to F+k, F-k*, Zk regardless of k).

Therefore, applying `S^ngrad` instead of `S` only **relabels** k-space positions (e.g. F+0 → F+2 instead of F+0 → F+1 per TR) but does not change any amplitude or relaxation history. The signal at F+0 is unchanged.

### Concrete pathway calculation

Consider the simplest refocusing pathway (the one-TR echo) for two cases.

**ngrad = 1 (one shift per TR):**

| Step | Operation | State |
|---|---|---|
| Start TR1 | Initial Mz | Z0 = 1 |
| TR1 RF (flip α, phase p₁) | Z0 → F+0 | F+0 = −½i·e^(−ip₁)·sin α |
| TR1 relax + shift (×1) | F+0 → F+1·E2 | F+1 = E2·(−½i·e^(−ip₁)·sin α) |
| TR2 RF (flip α, phase p₂) | F+1 → F−1\* | F−1\* = e^(2ip₂)·sin²(α/2)·F+1 |
| TR2 relax + shift (×1) | F−1\* → F+0·E2 | F+0 = E2²·½i·e^(ip₁−2ip₂)·sin²(α/2)·sin α |

**ngrad = 2 (two shifts per TR):**

| Step | Operation | State |
|---|---|---|
| Start TR1 | Initial Mz | Z0 = 1 |
| TR1 RF (flip α, phase p₁) | Z0 → F+0 | F+0 = −½i·e^(−ip₁)·sin α |
| TR1 relax + shift (×2) | F+0 → F+2·E2 | F+2 = E2·(−½i·e^(−ip₁)·sin α) |
| TR2 RF (flip α, phase p₂) | F+2 → F−2\* | F−2\* = e^(2ip₂)·sin²(α/2)·F+2 |
| TR2 relax + shift (×2) | F−2\* → F+0·E2 | F+0 = E2²·½i·e^(ip₁−2ip₂)·sin²(α/2)·sin α |

**The signal amplitude is identical in both cases:**

```
F+0 = ½i · E2² · exp(ip₁ − 2ip₂) · sin²(α/2) · sin α
```

The pathway visits k=2 (ngrad=2) instead of k=1 (ngrad=1) as an intermediate state, but it crosses **the same number of TR periods** and therefore accumulates **exactly the same T2 decay** (one factor of E2 per TR leg). The k-label is different; the amplitude is not.

The same argument extends to every pathway of every length: for each pathway contributing to F+0 with ngrad=1 (visiting orders k₀, k₁, k₂, …), there is a corresponding pathway with ngrad=2 visiting 2k₀, 2k₁, 2k₂, …, with identical amplitude. The total signal — the sum over all pathways — is therefore the same.

### Empirical verification

Running `EPG_GRE` with ngrad=1 and ngrad=10 (200 RF pulses) gives **identical signals** — the difference is identically zero. This confirms the theoretical result.

### Why the intuition "higher k-order = more decay" is wrong

It is tempting to reason: a stronger gradient pushes pathways to higher k-space orders, so they are further from F+0 and harder to refocus, giving a smaller signal. This reasoning is incorrect because **T2 decay in EPG is per TR, not per unit of k-space distance**. A state at k=2 after one TR has experienced the same T2 decay as a state at k=1 after one TR — the relaxation matrix E applies `E2` equally to both. The gradient amplitude only determines which k-label is attached to the magnetisation, not how much it has decayed.

---

## Q3 — Is it possible at all to model a stronger or weaker gradient spoiling in EPG?

**Without diffusion: no.** This is not a limitation of the code — it is a physical fact. EPG (like the isochromat model) assumes spins are **uniformly distributed** in phase across the voxel after any non-zero gradient. Whether the gradient creates 1 cycle/voxel or 100 cycles/voxel, the ensemble-averaged signal is identical, because only the phase *distribution* (uniform) matters, not the absolute gradient amplitude. Any non-zero gradient achieves the same uniform distribution.

**With diffusion: yes.** Diffusion introduces k-dependent attenuation:

```
attenuation(k) = exp(−b(k)·D),    b(k) ∝ k² · (γ·G·τ)²
```

Higher gradient amplitude → larger b-values → stronger attenuation of high-k states → those pathways genuinely cannot refocus. This is the only physical mechanism by which gradient amplitude affects the EPG signal.

To model this, use the `diff` parameter:

```matlab
d.G   = [G_prewinder, G_readout, G_spoiler];   % mT/m — all gradient events in TR
d.tau = [tau_pre,     tau_read,  tau_spoil];    % ms   — must sum to TR
d.D   = 0.8e-9;                                 % m²/s — tissue water diffusion

[s0,~,~] = EPG_GRE(theta, phi, TR, T1, T2, 'diff', d);
```

---

## Q4 — How are readout and spoiler gradients handled separately?

In SPGR, the readout gradient is **balanced**: the pre-winding gradient cancels the readout gradient echo-to-echo, so the net k-space shift per TR = 0. Only the **spoiler gradient** contributes to the EPG shift `S`.

- **Without diffusion**: the readout gradient is entirely irrelevant — only the *existence* of a spoiler matters (not its amplitude).
- **With diffusion** (`diff` parameter): all gradient events (pre-winder, readout, spoiler, and any zero-gradient periods) contribute to b-values. `diff.G` and `diff.tau` must list every segment so that `sum(diff.tau) = TR`.

---

## Q5 — What does this mean for the correction factor STheory/SEPG?

**Without diffusion:**

The correction factor depends **only on the RF phase cycling pattern** (φ₀ = 150°), TR, T1, T2, and flip angle. It is **independent of the spoiler gradient amplitude**.

The standard call:

```matlab
[s0,~,~] = EPG_GRE(d2r(alphaArray(1,ialpha))*ones(NTR,1), phi, TR, T1(iT1), T2);
```

already correctly models incomplete spoiling due to RF phase cycling. No modification is needed for gradient spoiling when diffusion is negligible.

**With diffusion:**

If diffusion effects are important (large gradients, long TR, high D), use the `diff` parameter. For the actual sequence parameters (spoiler ~1.3 rad/voxel/TR, TR = 4.1 ms, in-vivo tissue D ≈ 0.8×10⁻⁹ m²/s), diffusion attenuation will be small and the correction factor will be very close to the no-diffusion result.

---

## Q6 — Comparison with Hargreaves epg_rfspoil.m

Hargreaves' `epg_rfspoil.m` uses `epg_grelax.m` with `kg=1`. The `kg` parameter feeds only the diffusion b-value formula — it never changes the number of gradient shifts. `epg_grad.m` always applies exactly one shift per TR via `circshift`, identical to Malik's shift matrix `S`. The two implementations are mathematically equivalent and produce the same correction factor for the same RF phase cycling parameters.

---

## Summary

| Question | Answer |
|---|---|
| Does EPG model gradient spoiling? | Yes — shift matrix `S` already does this |
| Does gradient amplitude change SEPG (no diffusion)? | **No** — T2 decay is per TR, not per k-unit; all pathways are relabelled but unchanged in amplitude |
| Can a stronger/weaker gradient be modelled? | Only through diffusion (`diff` parameter) |
| How to model gradient amplitude effects? | Use the `diff` parameter with actual G and τ values |
| Is the original `EPG_GRE.m` call correct for this sequence? | **Yes** — no modification needed |
| What drives the correction factor? | RF phase cycling (φ₀) only, when diffusion is negligible |

---

## References

- Zur Y, Wood ML, Neuringer LJ. *Spoiling of transverse magnetization in steady-state sequences.* Magn Reson Med. 1991;21(2):251-263.
- Weigel M, Schwenk S, Kiselev VG, Scheffler K, Hennig J. *Extended phase graphs with anisotropic diffusion.* J Magn Reson. 2010;205(2):276-285.
- Yarnykh VL. *Actual flip-angle imaging in the pulsed steady state: A method for rapid three-dimensional mapping of the transmitted radiofrequency field.* Magn Reson Med. 2007;57(1):192-200.
- Preibisch C, Deichmann R. *Influence of RF spoiling on the stability and accuracy of T1 mapping based on spoiled FLASH with varying flip angles.* Magn Reson Med. 2009;61(1):125-135.
- Malik SJ, Teixeira RPAG, Hajnal JV. *Extended phase graph formalism for systems with magnetization transfer and exchange.* Magn Reson Med. 2018;80(2):767-779.
