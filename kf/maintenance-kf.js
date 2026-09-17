/* maintenance-kf.js — Kalman-filter maintenance (TDEE) estimator for HDH Log.
 *
 * Model (one step per calendar day; a morning weigh-in observes the state at the start of that day):
 *   state x = [ W  tissue weight (lb, glycogen/water excluded)
 *               M  maintenance (kcal/day)
 *               U  AR(1) water/gut noise (lb) ]
 *   W[t+1] = W[t] + (I[t] − M[t]) / rho              I[t] = intake eaten during day t, rho = kcal per lb of tissue
 *   M[t+1] = M[t] + beta·(W[t+1] − W[t]) + drift      beta ≈ 10 kcal/day per lb (maintenance falls as you get lighter)
 *   U[t+1] = phi·U[t] + eta                           water carryover, same AR(1) idea as the OLS band in the app
 *   scale[t] = W[t] + U[t] + G[t] + v                 G[t] = carb-driven glycogen+water offset computed from logged carbs
 *
 * What this buys over windowed OLS:
 *   • a prior on M (formula or last phase's estimate) so day-12 output is ±200, not ±800, and data takes over smoothly
 *   • every weigh-in weighted by how much it can tell you; robust gate stops one bad reading moving M by 250 kcal
 *   • phase-aware energy density (3,200 in a cut; 3,500 is pure fat)
 *   • unlogged days widen the band instead of silently averaging; no hard 10-day discard at phase changes
 *   • a daily time series of M with an honest band, and a projection of how many days until the band is ±100
 *
 * Usage:  const R = estimateMaintKF(rows, opts)
 *   rows: [{ d:int dayIndex | k:'YYYY-MM-DD', weight:lb|null, cals:kcal|null, carbs:g|null }]
 *         (from the app: build from allData() + dayTotals(), see adapterFromLog below)
 *   opts: { end, prior:{tdee, sd}, phase:'cut'|'maintain'|'bulk', phaseStarts:[{d|k}], kcalPerLb, ... see DEFAULTS }
 * Returns: { status, tdee, se, ci, trendW, tissueW, lbWk, pctWk, balance, avgIn, nW, nI, coverage, series, daysTo100, warnings, params }
 */
(function(root, factory){
  if(typeof module==='object' && module.exports) module.exports=factory();
  else root.MaintKF=factory();
})(typeof self!=='undefined'?self:this, function(){
'use strict';

const DEFAULTS = {
  // energy density of tissue change, kcal/lb. Pure fat is 3,500; a high-protein cut loses ~85–90% fat → ~3,200.
  kcalPerLbByPhase: { cut:3200, maintain:3500, bulk:3000 },
  beta: 10,            // kcal/day of maintenance per lb of tissue change
  phi: 0.5,            // AR(1) water carryover
  sigmaU: 0.7,         // stationary SD of AR(1) water, lb
  sigmaV: 0.35,        // scale resolution + i.i.d. daily noise, lb
  sigmaM: 6,           // maintenance random-walk drift, kcal/day per day (activity/adaptation you can't see)
  sigmaW: 0.03,        // unmodeled tissue process noise, lb/day
  priorSd: 1000,       // used when no prior given (≈ flat)
  imputeSd: 350,       // SD of true intake on an unlogged day, kcal — widens the band across gaps
  kappa: 0.01,         // lb of glycogen+water per g/day of carbs above the running reference
  alpha: 0.5,          // daily relaxation of the glycogen compartment toward its carb target
  carbRefHalfLife: 21, // days; reference carb level is an EMA so G is mean-zero over long windows
  gateC: 2.5,          // robust gate: normalized innovation beyond this is downweighted
  phaseInflateW: 1.0,  // lb, one-time inflation of W uncertainty at a phase change (absorbs water re-equilibration)
  phaseInflateM: 75,   // kcal, one-time inflation of M uncertainty at a phase change (activity may change)
  settleDays: 10,      // after a phase start, tissue dynamics are distrusted for this many days (glycogen/water shift)
  settleSd: 1.0,       // lb/day tissue process noise during settle (W free to track the scale; decoupled from M)
  unloggedBump: 0,     // kcal added to the imputed intake on unlogged days (set ~200–300 if unlogged days are usually bigger)
  firstDayIsPhaseStart: false, // set true when the first row is the phase start (app: phase.start == first weigh-in)
  minWeighs: 4,        // below this, status 'building'
  minIntake: 5,        // logged intake days needed before a number is shown (without intake, M is not identifiable)
  calibrate: true,     // inflate the band if innovations run larger than the model expects (never shrinks it)
};

function keyToDay(k){ const [y,m,d]=k.split('-').map(Number); return Math.round(Date.UTC(y,m-1,d)/86400000); }
function median(a){ if(!a.length) return null; const s=a.slice().sort((x,y)=>x-y); const m=s.length>>1; return s.length%2?s[m]:(s[m-1]+s[m])/2; }

// Katch–McArdle prior: BMR from lean mass × activity. sd 300 is honest for a formula.
function katchPrior(weightLb, bfPct, activity){
  const lbmKg=weightLb*(1-bfPct/100)*0.45359;
  const bmr=370+21.6*lbmKg;
  return { tdee: Math.round(bmr*(activity||1.55)), sd: 300 };
}
// Prior from a previous estimate (e.g. last phase's KF or OLS result): tdee ± ci (95%).
function carriedPrior(tdee, ci95, extraSd){ return { tdee, sd: Math.sqrt((ci95/1.96)**2 + (extraSd||100)**2) }; }

// 3×3 helpers (tiny, explicit; no library)
function matMul(A,B){ const C=[[0,0,0],[0,0,0],[0,0,0]]; for(let i=0;i<3;i++)for(let j=0;j<3;j++){let s=0;for(let k=0;k<3;k++)s+=A[i][k]*B[k][j];C[i][j]=s;} return C; }
function matT(A){ return [[A[0][0],A[1][0],A[2][0]],[A[0][1],A[1][1],A[2][1]],[A[0][2],A[1][2],A[2][2]]]; }
function matAdd(A,B){ return A.map((r,i)=>r.map((v,j)=>v+B[i][j])); }

function estimateMaintKF(rowsIn, opts){
  const P=Object.assign({}, DEFAULTS, opts||{});
  const rows=rowsIn.map(r=>Object.assign({}, r, { d: r.d!=null ? r.d : keyToDay(r.k) })).sort((a,b)=>a.d-b.d);
  const byDay={}; rows.forEach(r=>{ byDay[r.d]=r; });
  const weighDays=rows.filter(r=>r.weight>0).map(r=>r.d);
  const R={ status:'building', nW:weighDays.length, warnings:[], params:P };
  if(weighDays.length<P.minWeighs) return R;

  const d0=weighDays[0];
  const end = P.end!=null ? (typeof P.end==='string'?keyToDay(P.end):P.end) : weighDays[weighDays.length-1];
  // phases: [{d|k, phase:'cut'|'maintain'|'bulk'}] or bare days/keys. Sets energy density per day and the settle windows.
  const phases=(P.phaseStarts||[]).map(p=>{
    const o = (typeof p==='object') ? p : { d:p };
    return { d: o.k!=null ? keyToDay(o.k) : (typeof o.d==='string' ? keyToDay(o.d) : o.d), phase:o.phase };
  }).filter(p=>p.d<=end).sort((a,b)=>a.d-b.d);
  if(P.firstDayIsPhaseStart && !phases.some(p=>p.d===d0)) phases.unshift({ d:d0, phase:P.phase });
  const phaseAt=(d)=>{ let ph=P.phase||'cut'; for(const p of phases){ if(p.d<=d && p.phase) ph=p.phase; } return ph; };
  const rhoAt=(d)=> P.kcalPerLb || P.kcalPerLbByPhase[phaseAt(d)] || 3500;
  const settlingAt=(d)=> phases.some(p=>d>=p.d && d<p.d+P.settleDays);
  const phaseStarts=phases.map(p=>p.d).filter(d=>d>d0);
  let rho=rhoAt(d0);
  const prior = P.prior || null;

  // ── state & covariance ──
  let x=[byDay[d0].weight, prior?prior.tdee:2500, 0];
  let Pm=[[9,0,0],[0,(prior?prior.sd:P.priorSd)**2,0],[0,0,P.sigmaU**2]];
  const qU=P.sigmaU**2*(1-P.phi**2);

  // ── carb reference (EMA) and glycogen compartment ──
  const lam=Math.pow(0.5,1/P.carbRefHalfLife);
  const carbRows=rows.filter(r=>r.carbs>0 && r.d<d0+7);
  let cRef = carbRows.length ? carbRows.reduce((s,r)=>s+r.carbs,0)/carbRows.length : null;
  let G=0;
  const anyCarbs=rows.some(r=>r.carbs>0);

  // ── intake imputation: trailing median of logged days ──
  const loggedIn=[];
  function intakeFor(d){
    const r=byDay[d];
    if(r && r.cals>0){ loggedIn.push(r.cals); if(loggedIn.length>14) loggedIn.shift(); return { I:r.cals, logged:true }; }
    const allLogged=rows.filter(r=>r.cals>0 && r.d<d).map(r=>r.cals);
    const m = loggedIn.length>=3 ? median(loggedIn) : (allLogged.length ? median(allLogged) : null);
    return { I:m, logged:false };
  }

  const series=[]; let nI=0, nDaysIn=0, nisList=[], gated=0, sumIn14=[], missingIntake=0;
  const H=[1,0,1];

  for(let d=d0; d<=end; d++){
    let gatedToday=false, nu=null, nis=null;
    if(d>d0){
      // ── predict d−1 → d using intake eaten on day d−1 ──
      const it=intakeFor(d-1); nDaysIn++;
      let I=it.I, qW=P.sigmaW**2;
      rho=rhoAt(d-1);
      const settling = settlingAt(d-1);
      if(I==null){ I=x[1]; qW+= (P.imputeSd/rho)**2; }          // nothing logged yet: assume at maintenance, wide
      else if(!it.logged){ I+=P.unloggedBump||0; qW+=(P.imputeSd/rho)**2; missingIntake++; }
      else { nI++; sumIn14.push(I); if(sumIn14.length>14) sumIn14.shift(); }
      let F, Bu, Q;
      if(settling){
        // Water re-equilibrating after a phase change: the scale's movement is not energy balance.
        // Decouple W from M (weigh-ins still anchor W and U, but teach nothing about M) — same idea as the
        // app's 10-day skip, except the prior and pre-change information are kept.
        const qS=P.settleSd**2+qW;
        F=[[1,0,0],[0,1,0],[0,0,P.phi]]; Bu=[0,0,0];
        Q=[[qS,0,0],[0,P.sigmaM**2,0],[0,0,qU]];
      } else {
        F=[[1,-1/rho,0],[0,1-P.beta/rho,0],[0,0,P.phi]];
        Bu=[I/rho, P.beta*I/rho, 0];
        Q=[[qW,P.beta*qW,0],[P.beta*qW,P.beta**2*qW+P.sigmaM**2,0],[0,0,qU]];
      }
      x=[F[0][0]*x[0]+F[0][1]*x[1]+Bu[0], F[1][1]*x[1]+Bu[1], F[2][2]*x[2]];
      Pm=matAdd(matMul(matMul(F,Pm),matT(F)),Q);
      if(phaseStarts.includes(d)){ Pm[0][0]+=P.phaseInflateW**2; Pm[1][1]+=P.phaseInflateM**2; }
      // ── glycogen/water offset from yesterday's carbs ──
      if(anyCarbs){
        const c=byDay[d-1] && byDay[d-1].carbs>0 ? byDay[d-1].carbs : null;
        if(c!=null){ if(cRef==null) cRef=c; G=P.alpha*G+(1-P.alpha)*P.kappa*(c-cRef); cRef=lam*cRef+(1-lam)*c; }
      }
    }
    // ── update with today's weigh-in ──
    const r=byDay[d];
    if(r && r.weight>0){
      const z=r.weight-G;
      const HP=[Pm[0][0]+Pm[2][0], Pm[0][1]+Pm[2][1], Pm[0][2]+Pm[2][2]];   // H·P
      let S=HP[0]+HP[2]+P.sigmaV**2;                                          // H·P·Hᵀ + R
      nu=z-(x[0]+x[2]);
      nis=nu*nu/S; nisList.push(nis);
      if(P.gateC>0 && nis>P.gateC**2){ S=nu*nu/(P.gateC**2); gatedToday=true; gated++; }   // inflate R so |ν|/√S = gateC
      const K=[HP[0]/S, HP[1]/S, HP[2]/S];
      x=[x[0]+K[0]*nu, x[1]+K[1]*nu, x[2]+K[2]*nu];
      // Joseph form for numerical safety: (I−KH)P(I−KH)ᵀ + K R Kᵀ, with R the (possibly inflated) S − HPHᵀ
      const Rk=S-(HP[0]+HP[2]);
      const IKH=[[1-K[0],0,-K[0]],[-K[1],1,-K[1]],[-K[2],0,1-K[2]]];
      Pm=matAdd(matMul(matMul(IKH,Pm),matT(IKH)), [[K[0]*K[0]*Rk,K[0]*K[1]*Rk,K[0]*K[2]*Rk],[K[1]*K[0]*Rk,K[1]*K[1]*Rk,K[1]*K[2]*Rk],[K[2]*K[0]*Rk,K[2]*K[1]*Rk,K[2]*K[2]*Rk]]);
    }
    series.push({ d, tdee:x[1], se:Math.sqrt(Pm[1][1]), W:x[0], G, U:x[2], weight:r&&r.weight>0?r.weight:null, nu, gated:gatedToday });
  }

  // ── calibration: if innovations run larger than the model expects, widen the band (never narrow it) ──
  let scale=1;
  if(P.calibrate && nisList.length>=6){ const m=median(nisList)/0.4549; scale=Math.max(1, Math.sqrt(m)); }
  const se=Math.sqrt(Pm[1][1])*scale, tdee=x[1];
  const avgIn = sumIn14.length ? sumIn14.reduce((a,b)=>a+b,0)/sumIn14.length : null;
  const trendW=x[0]+G, tissueW=x[0];
  const lbWk = avgIn!=null ? (avgIn-tdee)/rho*7 : null;

  // ── projection: with daily weigh-ins + logged intake from here, when does the 95% band reach ±100? ──
  let daysTo100=null;
  { let Pp=Pm.map(r=>r.slice()); const F=[[1,-1/rho,0],[0,1-P.beta/rho,0],[0,0,P.phi]];
    const Q=[[P.sigmaW**2,P.beta*P.sigmaW**2,0],[P.beta*P.sigmaW**2,P.beta**2*P.sigmaW**2+P.sigmaM**2,0],[0,0,qU]];
    for(let k=1;k<=180;k++){
      Pp=matAdd(matMul(matMul(F,Pp),matT(F)),Q);
      const HP=[Pp[0][0]+Pp[2][0],Pp[0][1]+Pp[2][1],Pp[0][2]+Pp[2][2]], S=HP[0]+HP[2]+P.sigmaV**2, K=[HP[0]/S,HP[1]/S,HP[2]/S];
      const IKH=[[1-K[0],0,-K[0]],[-K[1],1,-K[1]],[-K[2],0,1-K[2]]];
      Pp=matAdd(matMul(matMul(IKH,Pp),matT(IKH)), [[K[0]*K[0]*P.sigmaV**2,K[0]*K[1]*P.sigmaV**2,K[0]*K[2]*P.sigmaV**2],[K[1]*K[0]*P.sigmaV**2,K[1]*K[1]*P.sigmaV**2,K[1]*K[2]*P.sigmaV**2],[K[2]*K[0]*P.sigmaV**2,K[2]*K[1]*P.sigmaV**2,K[2]*K[2]*P.sigmaV**2]]);
      if(1.96*Math.sqrt(Pp[1][1])*scale<=100){ daysTo100=k; break; }
    }
    if(daysTo100==null && 1.96*se<=100) daysTo100=0;
  }

  if(gated) R.warnings.push(`${gated} weigh-in${gated>1?'s':''} looked like water/gut noise (>${P.gateC}σ) and were downweighted.`);
  if(missingIntake) R.warnings.push(`Intake missing on ${missingIntake} of ${nDaysIn} days; filled with your trailing median and the band widened for it.`);
  if(scale>1.15) R.warnings.push(`Day-to-day scatter is ${(scale*100-100).toFixed(0)}% larger than the model assumes; band widened to match.`);
  if(prior==null) R.warnings.push('No prior supplied: the first ~2 weeks will be wide. Pass katchPrior() or the previous phase estimate.');

  Object.assign(R,{ status: nI>=P.minIntake?'ok':'building', tdee, se, ci:1.96*se, trendW, tissueW, avgIn, lbWk, pctWk: lbWk!=null?lbWk/trendW*100:null,
    balance: avgIn!=null?avgIn-tdee:null, nI, coverage: nDaysIn?nI/nDaysIn:0, series, daysTo100, gated, nisScale:scale, rho, prior });
  if(nI<P.minIntake) R.warnings.unshift(`Only ${nI} day${nI===1?'':'s'} of intake logged; maintenance can't be separated from intake until ${P.minIntake}+ days are logged.`);
  return R;
}

// Build rows from the app's synced day-map: allData() → { 'YYYY-MM-DD': { body:{weight}, nutrition|food } }
function adapterFromLog(data, dayTotals){
  return Object.keys(data).filter(k=>/^\d{4}-\d{2}-\d{2}$/.test(k)).sort().map(k=>{
    const t=dayTotals(data[k]);
    return { k, weight: data[k].body&&data[k].body.weight>0 ? +data[k].body.weight : null, cals: t.cals>0?t.cals:null, carbs: t.carb>0?t.carb:null };
  });
}

return { estimateMaintKF, katchPrior, carriedPrior, adapterFromLog, DEFAULTS, keyToDay };
});
