// Monte Carlo benchmark: the app's OLS-AR1 estimator vs the Kalman filter, on simulated cuts where
// true maintenance is known. The simulator is deliberately richer than the filter's model
// (adaptive thermogenesis, AR(1) water + carb glycogen with different constants, outliers, missing
// weigh-ins, weekend under-logging, logging error) so neither estimator is fitting its own generator.
'use strict';
const { estimateOLS } = require('./ols_ar1.js');
const KF = require('./maintenance-kf.js');

// ── seeded RNG ──
function mulberry32(a){ return function(){ a|=0; a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; }
function makeRng(seed){ const u=mulberry32(seed); return { u, n(){ let a=u()||1e-12, b=u(); return Math.sqrt(-2*Math.log(a))*Math.cos(2*Math.PI*b); }, uni(lo,hi){ return lo+(hi-lo)*u(); } }; }

// ── simulator ──
function simulate(rng, S){
  const N=S.days;
  const tdee0=rng.uni(2300,2800), intakeMean=rng.uni(1700,2000), W0=rng.uni(145,175);
  const rhoTrue=rng.uni(3100,3350), betaTrue=10, atMax=rng.uni(40,120);
  const carbBase=rng.uni(80,140), rampCarb = rng.u()<0.5 ? rng.uni(40,80) : 0, rampDay=Math.floor(rng.uni(8,20));
  let W=W0, U=0, G=0, glyStep = S.phaseStart ? -rng.uni(2,3) : 0;   // scenario B: glycogen/water drops over first days
  const days=[], truth=[];
  for(let d=0; d<N; d++){
    const dow=(d+S.dow0)%7, weekend=(dow===0||dow===6);
    const at = Math.min(1, d/28)*atMax;
    const tdee = tdee0 + betaTrue*(W-W0) - at;
    truth.push({ d, tdee, W });
    // scale reading (morning of day d)
    const transient = glyStep*(1-Math.exp(-d/2.5));       // phase-start: glycogen/water leaves over ~a week
    let z = W + U + G + transient + 0.3*rng.n();
    if(rng.u()<0.04) z += (rng.u()<0.5?-1:1)*rng.uni(1.5,2.5);
    const weighed = rng.u()>=0.08;
    // intake eaten during day d
    let logged = weekend ? intakeMean+150+200*rng.n() : intakeMean+150*rng.n();
    const unlogged = rng.u() < (weekend?0.30:0.05);
    let Itrue = unlogged ? intakeMean+250+250*rng.n() : logged+75*rng.n();
    const carbs = carbBase + (d>=rampDay?rampCarb:0) + (weekend?40:0) + 30*rng.n();
    days.push({ d, weight: weighed?+z.toFixed(1):null, cals: unlogged?null:Math.round(logged), carbs: unlogged?null:Math.round(carbs) });
    // advance
    W += (Itrue - tdee)/rhoTrue;
    U = 0.5*U + 0.65*rng.n();
    G = 0.6*G + 0.4*0.012*(carbs-carbBase);
  }
  return { days, truth, meta:{ tdee0, intakeMean, W0, rhoTrue, atMax } };
}

function summarize(errs){
  const n=errs.length; if(!n) return null;
  const rmse=Math.sqrt(errs.reduce((s,e)=>s+e.err*e.err,0)/n), bias=errs.reduce((s,e)=>s+e.err,0)/n;
  const cov=errs.filter(e=>Math.abs(e.err)<=e.ci).length/n, ci=errs.reduce((s,e)=>s+e.ci,0)/n;
  return { rmse, bias, cov, ci, n };
}

function run(scenario, nSims, seed){
  const horizons=[12,21,28,42,60];
  const est={
    'OLS (app, ρ3500)':      (rows,h,s,pri)=>estimateOLS(rows,{end:h-1,win:42,kcalPerLb:3500,startFloor:s.phaseStart?10:null}),
    'OLS ρ3200':             (rows,h,s,pri)=>estimateOLS(rows,{end:h-1,win:42,kcalPerLb:3200,startFloor:s.phaseStart?10:null}),
    'KF no prior':           (rows,h,s,pri)=>KF.estimateMaintKF(rows,{end:h-1,phase:'cut',firstDayIsPhaseStart:!!s.phaseStart}),
    'KF prior ±300':         (rows,h,s,pri)=>KF.estimateMaintKF(rows,{end:h-1,phase:'cut',prior:pri,firstDayIsPhaseStart:!!s.phaseStart}),
    'KF prior +300 biased':  (rows,h,s,pri)=>KF.estimateMaintKF(rows,{end:h-1,phase:'cut',prior:{tdee:pri.tdee+300,sd:300},firstDayIsPhaseStart:!!s.phaseStart}),
    'KF prior, no gate/carb':(rows,h,s,pri)=>KF.estimateMaintKF(rows,{end:h-1,phase:'cut',prior:pri,gateC:0,kappa:0,firstDayIsPhaseStart:!!s.phaseStart}),
  };
  const res={}; Object.keys(est).forEach(k=>{ res[k]={}; horizons.forEach(h=>res[k][h]=[]); });
  const rng=makeRng(seed);
  for(let i=0;i<nSims;i++){
    const sim=simulate(rng, Object.assign({days:60, dow0:Math.floor(rng.uni(0,7))}, scenario));
    const pri={ tdee: sim.meta.tdee0 + 300*rng.n(), sd:300 };   // formula prior: right on average, ±300 off
    for(const h of horizons){
      const rows=sim.days.filter(r=>r.d<h);
      const truthT=sim.truth[h-1].tdee;
      for(const [name,f] of Object.entries(est)){
        const R=f(rows,h,scenario,pri);
        if(R.status!=='ok') continue;
        res[name][h].push({ err:R.tdee-truthT, ci:R.ci });
      }
    }
  }
  return { res, horizons, nSims };
}

function table(title,{res,horizons,nSims}){
  console.log(`\n${title}  (${nSims} simulated cuts; error = estimate − true current maintenance, kcal/day)`);
  console.log('estimator'.padEnd(26)+horizons.map(h=>`day ${h}`.padStart(22)).join(''));
  console.log(' '.repeat(26)+horizons.map(()=>'RMSE  bias  cov  ±CI'.padStart(22)).join(''));
  for(const [name,byH] of Object.entries(res)){
    let line=name.padEnd(26);
    for(const h of horizons){
      const s=summarize(byH[h]);
      line += s ? `${s.rmse.toFixed(0).padStart(5)} ${s.bias>=0?'+':''}${s.bias.toFixed(0).padStart(4)} ${(s.cov*100).toFixed(0).padStart(3)}% ${s.ci.toFixed(0).padStart(5)}`.padStart(22)
                : `(building ${byH[h].length}/${nSims})`.padStart(22);
    }
    console.log(line);
  }
}

const N=+process.argv[2]||300;
table('A. Steady cut (glycogen already settled, weigh-ins from day 0)', run({phaseStart:false}, N, 11));
table('B. Cut starts day 0 (2–3 lb water drop over first week; app skips 10 days, KF distrusts tissue slope for 10)', run({phaseStart:true}, N, 23));
