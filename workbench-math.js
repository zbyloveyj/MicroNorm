/* Shared, deterministic calculations for the MicroNorm interface. */
(function(g){
  'use strict';
  const mean=a=>a.reduce((s,x)=>s+x,0)/a.length;
  const median=a=>{const s=[...a].sort((x,y)=>x-y),n=s.length;return n%2?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2;};
  function referenceAt(m,age,sex,j){
    const s=sex==='male'?0:1;
    let lo=m.mids.findIndex(v=>v>age)-1;
    if(lo<0)lo=age>=m.mids.at(-1)?m.mids.length-1:0;
    const hi=Math.min(lo+1,m.mids.length-1),t=hi===lo?0:Math.max(0,Math.min(1,(age-m.mids[lo])/(m.mids[hi]-m.mids[lo])));
    function at(b){const n=m.pos[s][b][j],w=n/(n+m.positive_prior);return {p:(n+m.detection_prior*m.pg[j])/(m.n[s][b][0]+m.detection_prior),q:m.locq[s][b][j].map((v,k)=>w*v+(1-w)*m.qg[j][k])};}
    const a=at(lo),b=at(hi);return {age,p:a.p*(1-t)+b.p*t,q:a.q.map((v,k)=>Math.expm1(v*(1-t)+b.q[k]*t))};
  }
  function match(rows,groups,demo=false){
    if(demo||groups.length!==2||rows.some(r=>!r.study))return {reason:demo?'demo':'study',pairs:[],studies:[],effects:[]};
    const used=new Set(),pairs=[],controls=rows.filter(r=>r.group===groups[0]);
    for(const item of rows.filter(r=>r.group===groups[1]).sort((a,b)=>a.age-b.age||a.sample_id.localeCompare(b.sample_id))){
      let best=null,gap=Infinity;
      for(const c of controls){if(used.has(c.sample_id)||c.study!==item.study||c.sex!==item.sex||(c.age<18)!==(item.age<18))continue;const d=Math.abs(c.age-item.age);if(d<gap&&d<=(item.age<18?2:5)){best=c;gap=d;}}
      if(best){used.add(best.sample_id);pairs.push({case:item.sample_id,control:best.sample_id,study:item.study,age_gap:gap,d:item.result.modules.map((v,k)=>v-best.result.modules[k]).concat(item.result.CMD-best.result.CMD)});}
    }
    if(!pairs.length)return {reason:'none',pairs,studies:[],effects:[]};
    const studies=[...new Set(pairs.map(p=>p.study))],by=studies.map(s=>pairs.filter(p=>p.study===s));
    const effects=[0,1,2,3].map(k=>mean(by.map(ps=>mean(ps.map(p=>p.d[k])))));
    const out={reason:pairs.length<30?'small':'tested',pairs,studies,effects,age_gap:mean(pairs.map(p=>p.age_gap)),p:[],q:[]};
    if(pairs.length<30)return out;
    let state=20261026;const rand=()=>{state^=state<<13;state^=state>>>17;state^=state<<5;return(state>>>0)/4294967296;};
    const exceed=[0,0,0,0];
    for(let b=0;b<10000;b++){const draw=by.map(ps=>{const sums=[0,0,0,0];for(const p of ps){const sign=rand()<.5?-1:1;for(let k=0;k<4;k++)sums[k]+=sign*p.d[k];}return sums.map(v=>v/ps.length);});for(let k=0;k<4;k++)if(Math.abs(mean(draw.map(v=>v[k])))>=Math.abs(effects[k])-1e-12)exceed[k]++;}
    out.p=exceed.map(v=>(v+1)/10001);out.q=[null,null,null,null];const order=[0,1,2].sort((a,b)=>out.p[a]-out.p[b]);let lowest=1;for(let i=2;i>=0;i--){const k=order[i];lowest=Math.min(lowest,out.p[k]*3/(i+1));out.q[k]=lowest;}return out;
  }
  g.MicroNormMath={mean,median,referenceAt,match};
})(window);
