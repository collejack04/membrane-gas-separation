# Membrane Gas Separation — Theory

Single-stage, binary, complete-mixing model. This is the model in
`src/membrane.py` and in `assets/app.js` (same solver, same numbers).

Example system used in the demo: **CO₂ (A, fast) / N₂ (B, slow)**,
flue-gas-like feed, 15% CO₂. The math is identical for H₂/CH₄,
O₂/N₂, CO₂/CH₄ — only the numbers for α and permeance change.

## 1. Picture

```
                retentate R, x (CO2 fraction x)
  feed F, xF  ──► [  membrane  ] ──►
                │  area A     │
                └──────┬──────┘
                       ▼ permeate P, y (CO2 fraction y)
```

- `F` = feed molar flow, `P` = permeate flow, `R` = retentate flow
- `θ = P / F` = stage cut (fraction of feed that permeates)
- Feed side pressure `pF`, permeate side pressure `pP`
- Pressure ratio `r = pP / pF` (always < 1; r = 0 is vacuum permeate)

## 2. The one equation to remember

For each component i, flux through the membrane is permeance
times partial-pressure driving force:

```
J_i = Q_i · (pF·x_i − pP·y_i)
```

- `J_i` [mol/m²·s] — molar flux of i
- `Q_i` [mol/m²·s·Pa] — permeance of i
- `pF·x_i` is the partial pressure of i on the feed side,
  `pP·y_i` on the permeate side. No driving force, no flux.

This is the solution–diffusion picture in one line: gas dissolves
in the polymer, diffuses across, desorbs. All the material science
is inside `Q_i`.

## 3. Permeability vs permeance — do not mix these up

| Quantity | Symbol | Units | What it is |
|---|---|---|---|
| Permeability | `P_i` | Barrer | Material property only |
| Permeance | `Q_i` | GPU | Material + how thin you made it |
| Selectivity | `α = P_A/P_B = Q_A/Q_B` | — | How much faster A goes than B |

```
Q_i [GPU] = P_i [Barrer] / thickness [µm]
```

1 Barrer = 10⁻¹⁰ cm³(STP)·cm / (cm²·s·cmHg).
1 GPU = 10⁻⁶ cm³(STP) / (cm²·s·cmHg) = 3.35×10⁻¹⁰ mol/(m²·s·Pa).

So 100 Barrer in a 0.1 µm selective layer is 1000 GPU.
Same polymer, 1 µm layer, is only 100 GPU — ten times the area
for the same separation. **Permeability sets area, selectivity
sets separation.** That is why the demo has separate sliders.

Typical orders of magnitude: α(CO₂/N₂) ≈ 20–50 for commercial
polymers, CO₂ permeance ≈ 100–2000 GPU for thin-film composites.

## 4. Balances

Steady state, no reaction, no accumulation:

```
Overall:  F = P + R
On A:     F·xF = P·y + R·x
```

Divide by `F`, use `θ = P/F`:

```
x = (xF − θ·y) / (1 − θ)
```

Recovery of A (fraction of fed A that ends up in the permeate):

```
Recovery = P·y / (F·xF) = θ·y / xF
```

Note the hard limit: `y ≤ xF/θ`, otherwise `x < 0`. You cannot
recover more A than you fed.

## 5. Complete-mixing model (Weller–Steiner)

Assumptions:

1. Binary, isothermal, steady state.
2. Both sides well mixed: feed-side composition is `x` everywhere,
   permeate-side composition is `y` everywhere. This is the
   simplest textbook model (Seader et al., Geankoplis). A real
   cross-flow module does slightly better because composition
   changes along the membrane — that is a follow-up, not this model.
3. `Q_i` constant — no plasticisation, no concentration
   polarisation, no pressure dependence, no defects/pinholes.
4. No pressure drop along the membrane, ideal gas.

Take the ratio of the two fluxes. Permeate composition is flux ratio:

```
y / (1−y) = J_A / J_B
          = (Q_A / Q_B) · (pF·x − pP·y) / (pF·(1−x) − pP·(1−y))
          = α · (x − r·y) / ((1−x) − r·(1−y))
```

This is the Weller–Steiner equation. Substitute `x` from §4 and
solve for `y`. It is implicit — `y` is on both sides — so the code
solves it by bisection on `y` in `(xF, xF/θ)`. Bisection is slow
but bulletproof, which matters for a live slider demo: it never
diverges at extreme θ or r.

Why `r` matters: if `r` is close to 1 (permeate pressure almost
equal to feed pressure), the driving force `(x − r·y)` collapses
and even a high-α membrane cannot separate. Vacuum or sweep on
the permeate side (small `r`) helps a lot. In the demo, permeate
is fixed at 1 bar, so raising feed pressure lowers `r`.

## 6. Membrane area

Total permeate flow is `P = θ·F`. Total flux is `J_A + J_B`
evaluated at the (mixed) compositions `x, y`:

```
A = P / (J_A + J_B)
```

Area scales linearly with feed flow and inversely with permeance.
It does **not** depend on α except through `x, y`.

## 7. Worked example (default demo case)

Inputs: `xF = 0.15`, `α = 30`, `P_CO₂ = 100 Barrer`,
thickness `0.1 µm` → `Q_CO₂ = 1000 GPU`,
`pF = 5 bar`, `pP = 1 bar` → `r = 0.20`,
`θ = 0.50`, `F = 10 mol/s`.

Solve §5 + §4:

```
y = 0.243  (permeate is 24.3% CO2)
x = 0.057  (retentate is 5.7% CO2)
Recovery = 0.50 × 0.243 / 0.15 = 0.810  (81.0%)
```

Fluxes:

```
Q_A = 1000 × 3.35e-10 = 3.35e-7 mol/m²·s·Pa
J_A = Q_A·(5e5×0.057 − 1e5×0.243) = 1.42e-3 mol/m²·s
J_B = (Q_A/30)·(5e5×0.943 − 1e5×0.757) = 4.42e-3 mol/m²·s
J_tot = 5.84e-3 mol/m²·s
A = (0.50 × 10) / 5.84e-3 = 857 m²
```

Check with: `python3 src/membrane.py` — it prints exactly these
numbers. The JavaScript demo prints the same.

Tradeoff at the same α, pressures and feed (vary θ only):

| θ | y (purity) | Recovery |
|---|---:|---:|
| 0.02 | 0.552 | 0.074 |
| 0.14 | 0.430 | 0.391 |
| 0.25 | 0.347 | 0.584 |
| 0.37 | 0.289 | 0.711 |
| 0.48 | 0.248 | 0.800 |
| 0.60 | 0.216 | 0.866 |
| 0.72 | 0.192 | 0.916 |
| 0.83 | 0.172 | 0.956 |
| 0.95 | 0.156 | 0.988 |

Low θ skims only the fastest gas: high purity, low recovery.
High θ forces almost everything through: high recovery, purity
falls toward the feed composition. One stage cannot give both.

## 8. What this model does not do (good presentation slide)

- No cross-flow / counter-current flow pattern (real modules beat
  complete mixing by 10–20% in purity at the same recovery).
- No concentration polarisation, fouling, or plasticisation.
- No multicomponent feed, no water vapour, no non-ideal gas.
- No cost, no compressor work, no two-stage/recycle design —
  which is exactly what a real plant needs, because one stage
  at 15% feed and α = 30 only reaches ~55% purity even at
  very low recovery.

## References

- Seader, Henley & Roper, *Separation Process Principles*, Ch. on membrane separations; Weller–Steiner complete-mixing model.
- Geankoplis, *Transport Processes and Separation Process Principles*.
- Weller, S. & Steiner, W. A. (1950). Separation of gases by fractional permeation through membranes. *J. Appl. Phys.* 21, 279–283.
- Baker, R. W. *Membrane Technology and Applications* (for realistic α and GPU ranges).
