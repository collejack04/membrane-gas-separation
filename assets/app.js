// Membrane gas separation — browser solver, mirrors src/membrane.py exactly.
// Complete-mixing (Weller-Steiner) model, binary CO2 (A, fast) / N2 (B, slow).
const GPU_TO_SI = 3.35e-10;
const BAR_TO_PA = 1e5;

function residual(y, x, alpha, r) {
  const num = x - r * y;
  const den = (1 - x) - r * (1 - y);
  if (den <= 1e-12 || num <= 0) return num <= 0 ? Infinity : -Infinity;
  return y / (1 - y) - alpha * num / den;
}

function solveStage(xF, alpha, r, theta) {
  const yHi = Math.min(0.999999, xF / theta * 0.999999);
  const yLo = xF * 1.000001 + 1e-9;
  const f = (y) => {
    const x = (xF - theta * y) / (1 - theta);
    if (x <= 1e-12 || x >= 1 - 1e-12) return y > yLo ? Infinity : -Infinity;
    return residual(y, x, alpha, r);
  };
  // Scan for a sign change, then bisect.
  let prevY = yLo, prevF = f(prevY), bracket = null;
  for (let i = 1; i <= 200; i++) {
    const yy = yLo + (yHi - yLo) * i / 200;
    const ff = f(yy);
    if (ff === 0) { bracket = [yy, yy]; break; }
    if ((prevF < 0 && ff > 0) || (prevF > 0 && ff < 0)) { bracket = [prevY, yy]; break; }
    prevY = yy; prevF = ff;
  }
  if (!bracket) {
    const y = yHi, x = Math.max(0, (xF - theta * y) / (1 - theta));
    return [y, x];
  }
  let [lo, hi] = bracket, flo = f(lo);
  for (let k = 0; k < 100; k++) {
    const mid = 0.5 * (lo + hi), fm = f(mid);
    if (fm === 0 || hi - lo < 1e-10) { lo = hi = mid; break; }
    if ((flo < 0 && fm > 0) || (flo > 0 && fm < 0)) hi = mid;
    else { lo = mid; flo = fm; }
  }
  const y = 0.5 * (lo + hi);
  const x = Math.max(0, Math.min(1, (xF - theta * y) / (1 - theta)));
  return [y, x];
}

function evaluate(p) {
  const r = p.pP / p.pF;
  const [y, x] = solveStage(p.xF, p.alpha, r, p.theta);
  const permeanceGPU = p.permBarrer / p.thickUm;
  const qA = permeanceGPU * GPU_TO_SI;
  const qB = qA / p.alpha;
  const jA = qA * (p.pF * BAR_TO_PA * x - p.pP * BAR_TO_PA * y);
  const jB = qB * (p.pF * BAR_TO_PA * (1 - x) - p.pP * BAR_TO_PA * (1 - y));
  const jTot = jA + jB;
  const permFlow = p.theta * p.feedFlow;
  return {
    y, x, r, permeanceGPU, jA, jB, jTot,
    recovery: p.theta * y / p.xF,
    area: jTot > 0 ? permFlow / jTot : Infinity,
    permFlow, retFlow: p.feedFlow - permFlow
  };
}

function tradeoff(p, n = 60) {
  const r = p.pP / p.pF, out = [];
  for (let i = 0; i < n; i++) {
    const theta = 0.02 + (0.95 - 0.02) * i / (n - 1);
    const [y, x] = solveStage(p.xF, p.alpha, r, theta);
    out.push({ theta, y, x, recovery: theta * y / p.xF });
  }
  return out;
}

const $ = (id) => document.getElementById(id);
const inputs = ["xF", "alpha", "permBarrer", "thickUm", "pF", "theta", "feedFlow"];

function readParams() {
  return {
    xF: +$("xF").value / 100,
    alpha: +$("alpha").value,
    permBarrer: +$("permBarrer").value,
    thickUm: +$("thickUm").value,
    pF: +$("pF").value,
    pP: 1.0,
    theta: +$("theta").value / 100,
    feedFlow: +$("feedFlow").value
  };
}

function fmtPct(v) { return (v * 100).toFixed(1) + "%"; }

function drawChart(curve, current) {
  const c = $("chart"), ctx = c.getContext("2d");
  const W = c.width, H = c.height, pad = 46;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = "#0f172a"; ctx.font = "12px system-ui";
  // axes: x = recovery 0-1, y = purity 0-1 (purity axis starts at 0)
  const X = (rec) => pad + rec * (W - pad - 16);
  const Y = (pur) => H - pad - pur * (H - pad - 16);
  ctx.strokeStyle = "#cbd5e1"; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(pad, 14); ctx.lineTo(pad, H - pad); ctx.lineTo(W - 16, H - pad); ctx.stroke();
  ctx.fillStyle = "#64748b";
  for (let t = 0; t <= 1.0001; t += 0.25) {
    ctx.fillText((t * 100).toFixed(0) + "%", X(t) - 10, H - pad + 18);
    ctx.fillText((t * 100).toFixed(0) + "%", 10, Y(t) + 4);
  }
  ctx.fillStyle = "#0f172a"; ctx.font = "bold 12px system-ui";
  ctx.fillText("Recovery of CO2 →", W / 2 - 50, H - 8);
  ctx.save(); ctx.translate(12, H / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText("Permeate purity (CO2) →", -70, 0); ctx.restore();
  // feed composition reference line
  const p = readParams();
  ctx.strokeStyle = "#e2e8f0"; ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(pad, Y(p.xF)); ctx.lineTo(W - 16, Y(p.xF)); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#64748b"; ctx.font = "11px system-ui";
  ctx.fillText("feed " + fmtPct(p.xF), W - 78, Y(p.xF) - 5);
  // curve
  ctx.strokeStyle = "#2563eb"; ctx.lineWidth = 2.5; ctx.beginPath();
  curve.forEach((pt, i) => { const px = X(pt.recovery), py = Y(pt.y); i ? ctx.lineTo(px, py) : ctx.moveTo(px, py); });
  ctx.stroke();
  // current point
  ctx.fillStyle = "#dc2626"; ctx.beginPath();
  ctx.arc(X(current.recovery), Y(current.y), 6, 0, Math.PI * 2); ctx.fill();
  ctx.strokeStyle = "#fff"; ctx.lineWidth = 2; ctx.stroke();
}

function update() {
  const p = readParams();
  $("xFVal").textContent = fmtPct(p.xF);
  $("alphaVal").textContent = "α = " + p.alpha.toFixed(0);
  $("permVal").textContent = p.permBarrer.toFixed(0) + " Barrer";
  $("thickVal").textContent = p.thickUm.toFixed(1) + " µm";
  $("pFVal").textContent = p.pF.toFixed(1) + " bar  (r = " + (1 / p.pF).toFixed(2) + ")";
  $("thetaVal").textContent = fmtPct(p.theta);
  $("feedVal").textContent = p.feedFlow.toFixed(0) + " mol/s";
  const res = evaluate(p);
  $("outPurity").textContent = fmtPct(res.y);
  $("outRecovery").textContent = fmtPct(res.recovery);
  $("outRetentate").textContent = fmtPct(res.x);
  $("outArea").textContent = res.area.toLocaleString(undefined, { maximumFractionDigits: 0 }) + " m²";
  $("outPermeance").textContent = res.permeanceGPU.toLocaleString(undefined, { maximumFractionDigits: 0 }) + " GPU";
  $("outFlux").textContent = (res.jTot * 1000).toFixed(2) + " mmol/m²·s";
  const set = (id, v) => { const el = $(id); if (el) el.textContent = v; };
  set("cardPurity", fmtPct(res.y));
  set("cardRecovery", fmtPct(res.recovery));
  set("cardArea", res.area.toLocaleString(undefined, { maximumFractionDigits: 0 }) + " m²");
  set("cardPermeance", res.permeanceGPU.toLocaleString(undefined, { maximumFractionDigits: 0 }) + " GPU");
  set("dimline", "α = " + p.alpha.toFixed(0) + "  ·  r = " + (1 / p.pF).toFixed(2) + "  ·  θ = " + fmtPct(p.theta) + "  ·  xF = " + fmtPct(p.xF) + "  ·  Q = " + res.permeanceGPU.toLocaleString(undefined, { maximumFractionDigits: 0 }) + " GPU  ·  pF = " + p.pF.toFixed(1) + " bar  ·  pP = 1 bar");
  $("explain").textContent =
    `At stage cut ${(p.theta * 100).toFixed(0)}%, ${fmtPct(res.recovery)} of the CO2 in the feed ends up in the permeate, ` +
    `at ${fmtPct(res.y)} purity. The retentate still holds ${fmtPct(res.x)} CO2. ` +
    `Push the stage cut up and recovery rises but purity falls — that is the tradeoff. ` +
    `Raise selectivity or feed pressure and the whole curve moves up. ` +
    `Raise permeability (or thin the layer) and purity barely moves, but area falls: you need ${res.area.toLocaleString(undefined, { maximumFractionDigits: 0 })} m² here.`;
  drawChart(tradeoff(p), res);
}

inputs.forEach(id => $(id).addEventListener("input", update));
$("resetBtn").addEventListener("click", () => {
  $("xF").value = 15; $("alpha").value = 30; $("permBarrer").value = 100;
  $("thickUm").value = 0.1; $("pF").value = 5; $("theta").value = 50; $("feedFlow").value = 10;
  update();
});
// Tabs, chetools style: App / Theory / Python / AI Project
function showTab(name) {
  document.querySelectorAll(".tab").forEach(b => b.classList.toggle("active", b.dataset.tab === name));
  document.querySelectorAll(".tabpanel").forEach(p => p.classList.toggle("active", p.id === "panel-" + name));
  const sidebar = $("sidebar");
  if (sidebar) sidebar.style.display = name === "app" ? "" : "none";
  if (name === "app") { const c = $("chart"); if (c) drawChart(tradeoff(readParams()), evaluate(readParams())); }
}
document.querySelectorAll(".tab").forEach(b => b.addEventListener("click", () => showTab(b.dataset.tab)));
// Python sub-tabs
function showPySub(name) {
  document.querySelectorAll(".subtab").forEach(b => b.classList.toggle("active", b.dataset.sub === name));
  document.querySelectorAll(".pysub").forEach(p => p.classList.toggle("active", p.id === "py-" + name));
}
document.querySelectorAll(".subtab").forEach(b => b.addEventListener("click", () => showPySub(b.dataset.sub)));
update();
