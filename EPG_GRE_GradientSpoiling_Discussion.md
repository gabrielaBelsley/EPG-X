# EPG_GRE Gradient Spoiling — Discussion Notes
**Date:** 2026-05-27  
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

## Q1 — How to incorporate gradient spoiling into `EPG_GRE.m`?

### EPG theory background

In the EPG formalism, a gradient lobe is represented as a **shift** of the state vector:

```
F⁺ₙ  →  F⁺ₙ₊₁   (shift matrix S applied once per TR)
```

- **RF spoiling** is already modelled via the quadratic phase increment `phi` (φ₀ = 150°).
- **Gradient spoiling** means applying a spoiler gradient with `ngrad` times the unit gradient area per TR. In EPG terms this replaces the single shift `S` with `S^ngrad` (S applied ngrad times per TR).

With larger `ngrad`, transverse coherences are pushed to higher k-space orders per TR, making them harder to refocus. As ngrad → ∞, SEPG → STheory (perfect spoiling), so the correction factor → 1.

### Change made to `EPG_GRE.m`

A new optional parameter `'ngrad'` (positive integer, **default = 1**) was added. `ngrad = 1` reproduces the original behaviour exactly.

**Key code changes (summary):**

```matlab
% 1. Parse new parameter
ngrad = 1; % default
if strcmpi(varargin{ii},'ngrad')
    ngrad = round(varargin{ii+1});
end

% 2. Scale kmax and kmax_per_pulse by ngrad
kmax = (np - 1) * ngrad;
kmax_per_pulse = ngrad * [1:ceil(np/2) (floor(np/2)):-1:1];

% 3. Compute S^ngrad (once, before the TR loop)
Sn = speye(N);
for k = 1:ngrad
    Sn = S * Sn;
end

% 4. Use S^ngrad in the composite relax-shift operator
SE = Sn * E;   % was: SE = S * E
```

Because `S` is a sparse permutation-like matrix, `S^ngrad` is also sparse and cheap to compute.

**Usage in `ModelSpoilCorrAllT1_CorrData.m`:**

```matlab
% Original (RF spoiling only):
[s0,~,~] = EPG_GRE(d2r(alphaArray(1,ialpha))*ones(NTR,1), phi, TR, T1(iT1), T2);

% With combined RF + gradient spoiling:
[s0,~,~] = EPG_GRE(d2r(alphaArray(1,ialpha))*ones(NTR,1), phi, TR, T1(iT1), T2, 'ngrad', 3);
```

> **Computational cost** scales linearly with `ngrad` (state-vector size grows proportionally).

---

## Q2 — For a given physical dephasing, what `ngrad` to use?

### Convention: 1 EPG unit = 2π rad/voxel

If defining the EPG unit as one complete phase cycle across a voxel (2π rad/voxel):

$$\text{ngrad} = \frac{\text{dephasing [rad/voxel per TR]}}{2\pi}$$

| Dephasing per TR | Cycles/voxel | ngrad |
|---|---|---|
| 2π | 1 | 1 |
| 4π | 2 | 2 |
| 8π | 4 | 4 |
| 16π | 8 | 8 |

### Multi-axis spoiling (read + phase + slice)

`EPG_GRE.m` is a **1D EPG** — it tracks dephasing along a single gradient axis. For spoilers in all three directions, a coherence pathway (k_read, k_phase, k_slice) can only contribute to the signal if **all three** k-components are simultaneously zero, making 3D spoiling more effective than 1D alone.

**Practical 1D approximation** (Yarnykh 2007, Preibisch & Deichmann 2009):

$$\text{ngrad}_\text{eff} = \text{ngrad}_R + \text{ngrad}_P + \text{ngrad}_S$$

| Simulation | Implication |
|---|---|
| `ngrad` = dominant direction only | Conservative — slightly overestimates correction needed |
| `ngrad` = sum of all 3 directions | Best approximation of actual 3D spoiling |
| True 3D EPG (not in this codebase) | Exact — requires 3D state vector |

> **Note on phase encoding:** In 3D acquisitions the phase-encode gradient itself contributes varying dephasing each TR as the phase-encode table steps, adding natural spoiling beyond the explicit spoiler. The effective ngrad_P may therefore be larger than the spoiler gradient alone, and the actual in-vivo correction factor is likely even closer to 1.

---

## Q3 — Actual dephasing is 0.1256 rad per axis (not 4π)

### Key numbers

| Quantity | Per axis | All 3 axes |
|---|---|---|
| Dephasing per TR | 0.1256 rad | 0.3768 rad |
| As fraction of 1 cycle | 0.02 cycles | 0.06 cycles |
| TRs to accumulate 1 full cycle | **50 TRs** | ~17 TRs |
| Total over NTR = 2000 | 251 rad ≈ 40 cycles | — |

Note: 0.1256 ≈ 2π/50 = π/25.

### Why the "4π → ngrad=2" rule does NOT apply here

Under the 2π-per-voxel convention:

$$\frac{0.1256}{2\pi} \approx 0.02 \quad \Rightarrow \quad \text{less than 1 EPG unit}$$

Setting `ngrad = 2, 4, 8` etc. would model a gradient **50–400× stronger** than the actual physical gradient. This would be incorrect.

### Correct interpretation

The EPG unit is **not fixed** — it is whatever the physical gradient step is. In this sequence:

- 1 EPG unit = 0.1256 rad/voxel per TR = the actual gradient step per axis
- `ngrad = 1` represents the real single-axis gradient
- For all three axes combined: **`ngrad = 3`** (arithmetic sum, one unit per axis)

### Physical implication

With only 0.02 cycles/voxel/TR, the gradient spoiling is **very weak**. The sequence relies almost entirely on **RF spoiling** (φ₀ = 150°). The dominant source of incomplete spoiling is the RF phase cycling pattern, not the gradient.

The correction factor with `ngrad = 3` will be only marginally closer to 1 than with `ngrad = 1`, because even three axes combined accumulate only ~0.06 cycles per TR — far from a strong-gradient regime.

### Recommended simulation call

```matlab
% Best approximation for 3-axis spoiling, each 0.1256 rad/TR:
[s0,~,~] = EPG_GRE(d2r(alphaArray(1,ialpha))*ones(NTR,1), phi, TR, T1(iT1), T2, 'ngrad', 3);

% Conservative (single axis only):
[s0,~,~] = EPG_GRE(d2r(alphaArray(1,ialpha))*ones(NTR,1), phi, TR, T1(iT1), T2);  % ngrad=1 default
```

---

## Summary table

| Scenario | `ngrad` | Notes |
|---|---|---|
| RF spoiling only (original code) | 1 (default) | Baseline — φ₀=150° quadratic phase cycling |
| + gradient spoiling, 1 axis, 4π/TR | 2 | Unit = 2π/voxel convention |
| + gradient spoiling, 3 axes, 4π/TR each | 6 | Arithmetic sum |
| + gradient spoiling, 1 axis, 0.1256 rad/TR | 1 | Unit = physical step |
| **+ gradient spoiling, 3 axes, 0.1256 rad/TR each** | **3** | **This sequence — recommended** |
| Perfect gradient spoiling (theoretical limit) | large | Correction → 1 |

---

## Files modified

| File | Change |
|---|---|
| `EPGX-src/EPG_GRE.m` | Added `'ngrad'` optional parameter for gradient spoiling |

**Commit:** `ee18297` — *Add ngrad parameter to EPG_GRE for combined RF + gradient spoiling*

---

## References

- Zur Y, Wood ML, Neuringer LJ. *Spoiling of transverse magnetization in steady-state sequences.* Magn Reson Med. 1991;21(2):251-263.
- Weigel M, Schwenk S, Kiselev VG, Scheffler K, Hennig J. *Extended phase graphs with anisotropic diffusion.* J Magn Reson. 2010;205(2):276-285.
- Yarnykh VL. *Actual flip-angle imaging in the pulsed steady state: A method for rapid three-dimensional mapping of the transmitted radiofrequency field.* Magn Reson Med. 2007;57(1):192-200.
- Preibisch C, Deichmann R. *Influence of RF spoiling on the stability and accuracy of T1 mapping based on spoiled FLASH with varying flip angles.* Magn Reson Med. 2009;61(1):125-135.
- Malik SJ, Teixeira RPAG, Hajnal JV. *Extended phase graph formalism for systems with magnetization transfer and exchange.* Magn Reson Med. 2018;80(2):767-779.
