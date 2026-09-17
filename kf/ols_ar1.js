// Faithful port of the maintenance estimator currently in index.html (OLS slope on daily
// weigh-ins, AR(1)-inflated slope SE, intake weighted to match the slope's day weights).
// Pure function on integer day indices so it can be benchmarked head-to-head.
'use strict';

function t95(df){
  df=Math.max(df,3);
  const z=1.959964, z3=z*z*z, z5=z3*z*z, z7=z5*z*z;
  return z + (z3+z)/(4*df) + (5*z5+16*z3+3*z)/(96*df*df) + (3*z7+19*z5+17*z3-15*z)/(384*df*df*df);
}

const RHO_MARGIN_Z = 1.0;
function ar1Fit(xs, e, sxx){
  const n=xs.length, xb=xs.reduce((a,b)=>a+b,0)/n;
  const u0=1/Math.sqrt(n), u1=xs.map(x=>(x-xb)/Math.sqrt(sxx));
  const pairs=[]; for(let i=1;i<n;i++) if(xs[i]-xs[i-1]<=2) pairs.push(i);
  const sse=e.reduce((a,v)=>a+v*v,0);
  let num=0; pairs.forEach(i=>{ num+=e[i]*e[i-1]; });
  const obs = sse>0 ? num/sse : 0;
  function stats(rho){
    const SU0=new Array(n).fill(0), SU1=new Array(n).fill(0);
    for(let i=0;i<n;i++) for(let j=0;j<n;j++){ const s=rho===0?(i===j?1:0):Math.pow(rho,Math.abs(xs[i]-xs[j])); SU0[i]+=s*u0; SU1[i]+=s*u1[j]; }
    let G00=0,G01=0,G11=0; for(let i=0;i<n;i++){ G00+=u0*SU0[i]; G01+=u0*SU1[i]; G11+=u1[i]*SU1[i]; }
    const trMS=n-(G00+G11);
    const msm=(i,j)=>{ const s=rho===0?(i===j?1:0):Math.pow(rho,Math.abs(xs[i]-xs[j]));
      return s - (u0*SU0[j]+u1[i]*SU1[j]) - (SU0[i]*u0+SU1[i]*u1[j]) + (u0*G00*u0 + u0*G01*u1[j] + u1[i]*G01*u0 + u1[i]*G11*u1[j]); };
    let tAM=0; pairs.forEach(i=>{ tAM+=msm(i,i-1); });
    return { g:tAM/trMS, trMS, cSc:G11/sxx };
  }
  let rho;
  if(obs<=stats(0).g) rho=0;
  else if(obs>=stats(0.95).g) rho=0.95;
  else { let lo=0, hi=0.95; for(let it=0;it<30;it++){ const mid=(lo+hi)/2; if(stats(mid).g<obs) lo=mid; else hi=mid; } rho=(lo+hi)/2; }
  const rhoUsed=Math.min(0.95, rho + RHO_MARGIN_Z*Math.sqrt((1-rho*rho)/n));
  const f=stats(rhoUsed);
  return { rho, rhoUsed, trMS:f.trMS, cSc:f.cSc, obs };
}

// rows: [{d:int day, weight:number|null, cals:number|null}], any order.
// opts: { end:int (default max d), win:int (default 42; 0 = all), kcalPerLb (3500), startFloor:int|null }
function estimateOLS(rows, opts){
  opts=opts||{};
  const KCAL_PER_LB = opts.kcalPerLb||3500;
  const MT_MIN = { weighs: 7, span: 10, intake: 7 };
  const byDay={}; rows.forEach(r=>{ byDay[r.d]=r; });
  const end = opts.end!=null ? opts.end : Math.max(...rows.map(r=>r.d));
  const win = opts.win==null ? 42 : opts.win;
  let startKey = win>0 ? end-(win-1) : -Infinity;
  if(opts.startFloor!=null && opts.startFloor>startKey) startKey=opts.startFloor;
  const keys=Object.keys(byDay).map(Number).filter(d=>d>=startKey && d<=end).sort((a,b)=>a-b);
  const weighKeys=keys.filter(d=>byDay[d].weight>0);
  const R={ nW:weighKeys.length };
  R.first=weighKeys[0]; R.last=weighKeys[weighKeys.length-1];
  R.span = weighKeys.length>1 ? R.last-R.first : 0;
  const intake = weighKeys.length>1
    ? keys.filter(d=>d>=R.first && d<R.last).map(d=>({d, c:byDay[d].cals||0})).filter(p=>p.c>0)
    : [];
  R.nI=intake.length;
  R.coverage = R.span>0 ? R.nI/R.span : 0;
  if(R.nW<MT_MIN.weighs || R.span<MT_MIN.span || R.nI<MT_MIN.intake){ R.status='building'; return R; }

  const pts=weighKeys.map(d=>({ d, x:d-R.first, y:+byDay[d].weight }));
  const n=pts.length;
  const xb=pts.reduce((s,p)=>s+p.x,0)/n, yb=pts.reduce((s,p)=>s+p.y,0)/n;
  let sxx=0, sxy=0; pts.forEach(p=>{ sxx+=(p.x-xb)**2; sxy+=(p.x-xb)*(p.y-yb); });
  const slope=sxy/sxx, icpt=yb-slope*xb;
  pts.forEach(p=>{ p.fit=icpt+slope*p.x; p.e=p.y-p.fit; });
  const sse=pts.reduce((s,p)=>s+p.e*p.e,0);
  const ar=ar1Fit(pts.map(p=>p.x), pts.map(p=>p.e), sxx);
  const r1=ar.rho;
  const sigma2=Math.max(sse/ar.trMS, 0.04);
  const seSlope=Math.sqrt(sigma2*ar.cSc);
  const nEff=n*(1-r1)/(1+r1);
  const df=Math.max(nEff-2, 1);

  const inByDay={}; intake.forEach(p=>{ inByDay[p.d-R.first]=p.c; });
  let wSum=0, wcSum=0;
  for(let j=0;j<R.span;j++){
    if(inByDay[j]==null) continue;
    let w=0; pts.forEach(p=>{ if(p.x>j) w+=(p.x-xb)/sxx; });
    wSum+=w; wcSum+=w*inByDay[j];
  }
  const avgIn = wSum>0 ? wcSum/wSum : intake.reduce((a,p)=>a+p.c,0)/R.nI;
  const tdee=avgIn - slope*KCAL_PER_LB;
  const ci=t95(df)*seSlope*KCAL_PER_LB;
  const trendW=icpt+slope*R.span;
  Object.assign(R,{ status:'ok', slope, icpt, sigma:Math.sqrt(sigma2), r1, nEff, tdee, ci, se:ci/t95(df), avgIn, trendW,
    lbWk:slope*7, balance:slope*KCAL_PER_LB, residSD:Math.sqrt(sse/(n-2)) });
  return R;
}

module.exports = { estimateOLS, t95, ar1Fit };
