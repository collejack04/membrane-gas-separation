# Membrane Gas Separation — Interactive Demo

Single-stage CO₂ / N₂ membrane separation, live in the browser.
An AI project exploring AI capabilities, for Separations and Processes.

**Live demo:** https://collejack04.github.io/membrane-gas-separation/
**Repo:** https://github.com/collejack04/membrane-gas-separation
(open `index.html` locally — no build step, no server needed)

## What it does

Sliders for the things that actually matter:

- **Feed CO₂ fraction** (5–50%)
- **Selectivity α** (5–100) — moves the purity/recovery curve up
- **CO₂ permeability** (10–1000 Barrer) + **selective-layer thickness** (0.1–2 µm) — sets permeance, hence membrane area, not purity
- **Feed pressure** (2–20 bar, permeate fixed at 1 bar) — sets pressure ratio `r`
- **Stage cut θ** (5–90%) — the tradeoff slider
- **Feed flow** (1–20 mol/s) — scales area linearly

Outputs, live: permeate purity, CO₂ recovery, retentate composition,
membrane area, permeance, total flux — plus a **purity vs recovery
tradeoff curve** with your current operating point marked.

Default case: 15% CO₂, α = 30, 100 Barrer / 0.1 µm (= 1000 GPU),
5 bar → 1 bar, θ = 50%, 10 mol/s → **24.3% purity, 81.0% recovery,
857 m²**. Check it with `python3 src/membrane.py`.

## The math in one line

```
J_i = Q_i · (p_feed·x_i − p_perm·y_i)      flux = permeance × driving force
α = Q_A / Q_B                              selectivity
Q [GPU] = P [Barrer] / thickness [µm]
y/(1−y) = α·(x − r·y) / ((1−x) − r·(1−y))  Weller–Steiner, complete mixing
```

Full derivation, assumptions, worked example and limits:
[`docs/theory.md`](docs/theory.md)

## Repo layout

```
index.html              static demo (GitHub Pages, no build)
assets/app.js           browser solver — same math as the Python model
assets/style.css
src/membrane.py         reference Python model + tradeoff curve + CLI
docs/theory.md          full theory, worked example, references
```

The Python and JavaScript solvers are intentionally identical
(same bisection on the Weller–Steiner equation) and verified to
agree to 4 significant figures — see the numbers above.

## Run it

```bash
python3 src/membrane.py
python3 -m http.server 8000   # then open http://localhost:8000
```

## Publish to GitHub Pages

Repo is https://github.com/collejack04/membrane-gas-separation

```bash
git remote add origin https://github.com/collejack04/membrane-gas-separation.git
git push -u origin main
```

Then GitHub → Settings → Pages → Deploy from branch → `main`, folder `/ (root)`.
Or use the included workflow (`.github/workflows/pages.yml`) — it deploys
on every push to `main`. The site is static, so there is no build step;
a standalone WebAssembly build is not required for the assignment's
“publicly accessible website demo” requirement.

## Assignment mapping

1. **Website demo/documentation** — this repo / GitHub Pages site.
2. **YouTube video** — link to be added in `index.html` → “About this AI project”.
3. **Verbatim transcript** — viewable section to be added on the site after recording.
4. **Presentation/documentation** — theory above + AI process notes in `index.html`.
5. **Compare ≥2 AI models** — plan from one model, critique from another,
   implement with at least two, justify accept/reject. Document in the site.
6. **Attribution** — add group members / contributions here before submission.
7. **In-class presentations** — last day of class + final exam time/location.

## What I learned / want to learn next

- Selectivity alone does not set purity — pressure ratio and stage cut matter as much.
- Permeability moves area, not the separation curve.
- One stage cannot give high purity *and* high recovery at 15% feed.
- Next: two-stage design with recycle, and cross-flow (not complete-mixing) modelling.
