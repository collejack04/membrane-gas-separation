# Requirement 5 — Two-model comparison (draft, not yet run)

Project: Membrane gas separation, single-stage CO2 / N2.
Live site: https://collejack04.github.io/membrane-gas-separation/
Repo: https://github.com/collejack04/membrane-gas-separation
Status 2026-11-19: site is live and verified — default case 24.3% purity,
81.0% recovery, 857 m2. This file is a local draft only; it is NOT on
GitHub yet. Upload it (Add file > Upload files) when the critique below
has actually been run, and fill in the names/dates in the table.

## 1. Model 1 — the plan that was actually built

Model 1: Muse (the assistant in this chat), 2026-11-19.
Topic chosen after McCabe-Thiele was rejected — Jack had already done
it for two videos this semester and wanted something new to learn.

Plan Model 1 gave, and what was built from it:
- Single-stage, binary, complete-mixing (Weller-Steiner) model.
  J_i = Q_i * (pF*x_i - pP*y_i); alpha = Q_A/Q_B;
  Q [GPU] = P [Barrer] / thickness [um];
  y/(1-y) = alpha * (x - r*y) / ((1-x) - r*(1-y)), r = pP/pF.
  Solved by bisection on y.
- Static site (index.html + assets/app.js + assets/style.css), no build
  step, GitHub Pages. Python reference model in src/membrane.py with
  the identical solver, verified to agree.
- Sliders: feed CO2% (15% default), selectivity alpha (30), permeability
  (100 Barrer), thickness (0.1 um = 1000 GPU), feed pressure (5 bar,
  permeate 1 bar), stage cut theta (50%), feed flow (10 mol/s).
- Live outputs + purity-vs-recovery curve with the current point marked.
- Full derivation, assumptions and limits in docs/theory.md.
- Deliberate scope decisions in the plan: complete-mixing, not
  cross-flow; one stage, not two-stage/recycle; static JS, not a
  WebAssembly build (the brief recommends a standalone WASM build on
  GitHub Pages, it does not require one — the static demo satisfies
  "publicly accessible website demo").

Verified numbers (python3 src/membrane.py, same as the live site):
- Default (theta=0.50): y = 0.243, recovery = 0.810, area = 857 m2,
  retentate x = 0.057.
- Low cut (theta=0.02): y = 0.552, recovery = 0.074.
- High cut (theta=0.83): y = 0.172, recovery = 0.956.
- Live check 2026-11-19: moving theta 50% -> 62% gave y = 21.2%,
  recovery = 87.5%.

## 2. Model 2 — paste-ready critique prompt (run this in a DIFFERENT model)

Copy everything below the line into a second model (e.g. ChatGPT,
Claude, or Gemini — pick one Jack has access to, and record its name
in the table in section 3). Attach docs/theory.md and src/membrane.py,
or paste the plan in section 1 above.

---
You are critiquing a plan for an undergraduate Separations and
Processes AI project, not building it from scratch. The plan is in
section 1 above: a single-stage, binary, complete-mixing membrane
gas separation demo (CO2/N2, Weller-Steiner model, bisection solver),
as a static GitHub Pages site with sliders for feed CO2%, selectivity,
permeability, thickness, feed pressure, stage cut and feed flow, plus
a live purity-vs-recovery curve.

Critique it in detail, and be specific enough that I can accept or
reject each point:
1. Math: is the Weller-Steiner equation stated and solved correctly?
   Check the driving-force terms, the pressure ratio r, the balance
   x = (xF - theta*y)/(1-theta), recovery = theta*y/xF, and the hard
   limit y <= xF/theta. Flag any sign, unit, or definition error.
2. Units: permeability (Barrer) vs permeance (GPU) — is
   Q [GPU] = P [Barrer] / thickness [um] right, and is the 100 Barrer
   / 0.1 um = 1000 GPU example right? Is the area calculation
   (area = permeate flow / total flux) dimensionally consistent?
3. Model limits: what does complete-mixing miss vs a real cross-flow
   module? Does it over- or under-predict purity, and by roughly how
   much? What assumptions (constant Q, no plasticisation, no pressure
   drop, ideal gas, isothermal) are most likely to break first for
   CO2/N2, and should any be stated more strongly?
4. Demo design: are the slider ranges sensible (alpha 5-100, feed
   2-20 bar, theta 5-90%)? Is there a range where the solver breaks,
   gives y > 1, or a misleading result? Is the default case
   (15% CO2, alpha=30, 5 bar -> 1 bar, theta=50% -> 24.3% purity,
   81.0% recovery, 857 m2) realistic for a real polymer membrane?
5. Presentation story: is "one stage cannot give high purity AND
   high recovery" the right takeaway, and what is the single most
   important number or slider move to show live in class?
6. What is missing for an undergraduate to explain this cold, and
   what would you change, add, or remove — with a reason for each?
   Rank your changes: must-fix, should-fix, nice-to-have.
---

## 3. Accept / reject log (fill in after Model 2 replies)

| # | Model 2 point | Accept / reject | Reason / what changed |
|---|---------------|-----------------|-----------------------|
| 1 | (paste the point) | | |
| 2 | | | |
| 3 | | | |

Decisions already made by Jack / Model 1 (record these as the baseline,
so the comparison is honest):
- Complete-mixing over cross-flow: ACCEPTED (for now) — explainable in
  one equation, textbook model (Seader et al., Geankoplis); cross-flow
  is the stated "want to learn next".
- One stage over two-stage/recycle: ACCEPTED — the tradeoff is the
  point; two-stage is the follow-up, not this demo.
- Static JS over WebAssembly: ACCEPTED — brief recommends WASM, does
  not require it; static site deploys as-is and the Python model is
  the auditable reference.
- Bisection over a closed-form solve: ACCEPTED — never diverges at
  extreme theta or r, which matters for live sliders.

Model 2 name: __________________  Date run: ____________
Model 1 name: Muse  Date: 2026-11-19

## 4. Presentation outline (for the YouTube video, req. 2-3)

Jack's usual pattern: record, verbatim transcript, critique against
the source, re-record if needed. Same shape works here — record a
practice run first, then the final.
Suggested run-through, ~6-8 minutes, demoing the live site:
1. (30s) Why membranes, why CO2/N2 — carbon capture / flue gas, 15% feed.
2. (1 min) The one equation: flux = permeance x driving force; what
   lives inside Q (solution-diffusion, in one sentence).
3. (1 min) Permeability vs permeance — 100 Barrer / 0.1 um = 1000 GPU;
   permeability sets area, selectivity sets separation.
4. (2 min) Live demo: default 24.3% / 81.0% / 857 m2; move theta to
   ~5% (purity up, recovery collapses) and to ~85% (recovery up,
   purity falls to near feed); move alpha and feed pressure and show
   the curve move.
5. (1 min) The limit: one stage cannot do both; what a second stage
   or recycle would do.
6. (1 min) AI process: plan from Model 1, critique from Model 2,
   what was accepted/rejected and why (section 3 above); what I
   learned / want to learn next.

After recording: paste the practice transcript + AI critique in the
YouTube comments (req. 2), and put the verbatim auto-generated
transcript as a viewable section on the site (req. 3) — both are still
placeholders in index.html, "About this AI project".

## 5. Remaining checklist (req. 1-7)

- [x] 1. Website demo — live at collejack04.github.io/membrane-gas-separation,
      all 7 files in the repo (index.html, README.md, assets/app.js,
      assets/style.css, docs/theory.md, src/membrane.py, .gitignore).
      Note: .github/workflows/pages.yml exists locally but was NOT
      uploaded; not needed — Pages was turned on manually in Settings.
      Upload it only if you want auto-deploy on future pushes.
- [ ] 2. YouTube presentation video — not recorded yet; link to be
      added to index.html.
- [ ] 3. Verbatim transcript — viewable section to be added on the site
      after recording.
- [ ] 4. Documentation — theory.md is done; fill in the "About this
      AI project" section in index.html with the Model 1 / Model 2
      names, prompts, iteration notes and this file's accept/reject log.
- [ ] 5. Two-model comparison — this file is the draft; run the prompt
      in section 2 in a second model, then fill in section 3.
- [ ] 6. Attribution — solo or group? Not yet stated. If solo, record
      "individual submission, max 50 points"; if a group, add each
      member + contribution before submission.
- [ ] 7. In-class presentations — last day of class + final exam
      time/location; dates not yet in the brief Jack pasted.

Due: "in a couple of weeks" from 2026-11-19 — no exact date in the
brief Jack pasted. Confirm the exact deadline and the solo/group
question before recording, so the video and attribution match.
