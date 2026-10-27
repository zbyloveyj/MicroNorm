(function () {
  'use strict';
  const $ = selector => document.querySelector(selector);
  const C = window.MicroNormCore;
  const BLUE = '#278cbd', PINK = '#d75482', GOLD = '#dda922';
  let model = null, singleResult = null, cohortResults = [];
  const fmt = (v, digits = 3) => Number(v).toFixed(digits);
  const esc = s => String(s).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  const mean = a => a.reduce((x, y) => x + y, 0) / a.length;
  const median = a => { const b = [...a].sort((x, y) => x - y); return b.length % 2 ? b[(b.length - 1) / 2] : (b[b.length / 2 - 1] + b[b.length / 2]) / 2; };
  function download(content, name, type) {
    const url = URL.createObjectURL(new Blob([content], {type}));
    const a = document.createElement('a'); a.href = url; a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1500);
  }
  function setMessage(which, message, error = false) {
    const el = $('#' + which + '-message'); el.textContent = message; el.classList.toggle('error', error);
  }
  function downloadSingle() {
    if (!singleResult) return;
    download(JSON.stringify(singleResult, null, 2), 'micronorm-sample-report.json', 'application/json');
  }
  function downloadSpecies() {
    if (!singleResult) return;
    const lines = ['species\tabundance_percent\tdetection_probability\tcentile'];
    for (const r of singleResult.rows) lines.push([r.species, r.abundance_percent, r.detection_probability, r.centile].join('\t'));
    download(lines.join('\n') + '\n', 'micronorm-species-centiles.tsv', 'text/tab-separated-values');
  }
  function renderSingle(r) {
    const notable = r.rows.filter(x => x.abundance_percent > 0 && x.detection_probability >= .1)
      .sort((a, b) => Math.abs(b.centile - .5) - Math.abs(a.centile - .5)).slice(0, 8);
    const status = r.age_support < 300 ? '<span class="badge caution">年龄样本数较少</span>' : '<span class="badge">该年龄段样本数 ≥300</span>';
    $('#single-result').innerHTML = `<div class="result-top"><div><span class="eyebrow">ECOLOGICAL DEVIATION REPORT</span><h3>相对于健康参照的位置</h3></div>${status}</div>
      <div class="metrics">${r.modules.map((v, i) => `<div class="metric"><b>${Math.round(v * 100)}%</b><span>M${i + 1} 平均物种百分位</span></div>`).join('')}<div class="metric"><b>${fmt(r.CMD, 2)}</b><span>CMD 联合距离</span></div></div>
      <div class="axis-panel"><div class="axis-row"><span>M1–M2 相对分化</span><div class="axis-track"><i style="left:${Math.max(0,Math.min(100,50+(r.modules[0]-r.modules[1])*100))}%"></i></div><b>${fmt(r.differentiation_axis, 3)}</b></div><div class="axis-row"><span>M1+M2 共同位置</span><div class="axis-track"><i style="left:${Math.max(0,Math.min(100,r.common_axis/Math.SQRT2*100))}%"></i></div><b>${fmt(r.common_axis, 3)}</b></div></div>
      <p class="taxa-title">检出物种中距中位位置较远的示例 / DESCRIPTIVE ONLY</p><table class="taxa-table"><thead><tr><th>物种</th><th>丰度</th><th>百分位</th></tr></thead><tbody>${notable.map(x => `<tr><td>${esc(x.species.replace(/^s__/, '').replaceAll('_', ' '))}</td><td>${fmt(x.abundance_percent, 2)}%</td><td>${Math.round(x.centile * 100)}%</td></tr>`).join('')}</tbody></table>
      <ul class="warn-list">${r.warnings.map(x => `<li>${esc(x)}</li>`).join('')}</ul>
      <div class="result-actions"><button id="download-report" type="button">下载完整报告 JSON ↓</button><button id="download-species" type="button">下载801物种表 TSV ↓</button></div>`;
    $('#download-report').onclick = downloadSingle;
    $('#download-species').onclick = downloadSpecies;
  }
  function normalizedAbundance(d) {
    const sum = Object.values(d).reduce((a, b) => a + b, 0);
    return Object.fromEntries(Object.entries(d).map(([k, v]) => [k, sum > 0 ? v * 95 / sum : 0]));
  }
  async function demoSingle() {
    try {
      const demos = window.MicroNormDemos;
      if(!Array.isArray(demos)) throw Error('离线演示数据不可用。');
      const d = demos.find(x => x.id === 'adult') || demos[0];
      $('#single-age').value = d.age; $('#single-sex').value = d.sex;
      singleResult = C.score(model, d.age, d.sex, d.abundance);
      renderSingle(singleResult);
      setMessage('single', '已载入合成演示谱；这是界面示例，不是真实受试者或疾病检测。');
    } catch (error) { singleResult=null; $('#single-result').innerHTML='<div class="empty-state"><h3>未生成结果</h3><p>请检查输入后重试。</p></div>'; setMessage('single', error.message, true); }
  }
  function svgScatter(rows, groups) {
    const w = 330, h = 210, x0 = 43, x1 = 310, y0 = 24, y1 = 176;
    const xp = x => x0 + Math.max(0, Math.min(1, x)) * (x1 - x0);
    const yp = y => y1 - Math.max(0, Math.min(1, y)) * (y1 - y0);
    const grid = [0,.25,.5,.75,1].map(v => `<line x1="${xp(v)}" x2="${xp(v)}" y1="${y0}" y2="${y1}" stroke="#e5ebeb"/><line x1="${x0}" x2="${x1}" y1="${yp(v)}" y2="${yp(v)}" stroke="#e5ebeb"/>`).join('');
    const points = rows.map((r, i) => `<circle cx="${xp(r.result.modules[0])}" cy="${yp(r.result.modules[1])}" r="${rows.length > 60 ? 3 : 5}" fill="${r.group === groups[0] ? BLUE : PINK}" fill-opacity=".72"><title>${esc(r.sample_id)} · ${esc(r.group)} · M1 ${fmt(r.result.modules[0],3)} · M2 ${fmt(r.result.modules[1],3)}</title></circle>`).join('');
    return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="M1与M2样本分布">${grid}<line x1="${x0}" x2="${x1}" y1="${y1}" y2="${y1}" stroke="#91a5a7"/><line x1="${x0}" x2="${x0}" y1="${y0}" y2="${y1}" stroke="#91a5a7"/>${points}<text x="${x1}" y="205" text-anchor="end" font-size="10" fill="#657984">M1 percentile</text><text x="14" y="${y0}" font-size="10" fill="#657984">M2</text></svg>`;
  }
  function svgCMD(rows, groups) {
    const w = 330, h = 210, x0 = 43, x1 = 310, max = Math.max(2, ...rows.map(r => r.result.CMD)) * 1.08;
    const xp = x => x0 + x / max * (x1 - x0);
    const baseline = [75, 145];
    const guides = [0,.25,.5,.75,1].map(f => `<line x1="${x0+(x1-x0)*f}" x2="${x0+(x1-x0)*f}" y1="25" y2="166" stroke="#e6ebeb"/><text x="${x0+(x1-x0)*f}" y="195" text-anchor="middle" fill="#6b7b81" font-size="9">${fmt(max*f,1)}</text>`).join('');
    const points = rows.map((r, i) => { const gi = r.group === groups[0] ? 0 : 1; const jitter = ((i * 17) % 21) - 10; return `<circle cx="${xp(r.result.CMD)}" cy="${baseline[gi]+jitter}" r="${rows.length > 60 ? 3 : 4.5}" fill="${gi ? PINK : BLUE}" fill-opacity=".68"><title>${esc(r.sample_id)} · CMD ${fmt(r.result.CMD,2)}</title></circle>`; }).join('');
    return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="各组CMD样本分布">${guides}<text x="4" y="78" font-size="9" fill="${BLUE}">${esc(groups[0].slice(0,10))}</text><text x="4" y="148" font-size="9" fill="${PINK}">${esc((groups[1]||'').slice(0,10))}</text>${points}<text x="${x1}" y="208" text-anchor="end" font-size="10" fill="#657984">CMD</text></svg>`;
  }
  function cohortTable(rows, groups) {
    const summary = groups.map(g => {
      const a = rows.filter(x => x.group === g), ages = a.map(x => x.age);
      return {group:g,n:a.length,age:mean(ages),female:a.filter(x => x.sex === 'female').length,
        m:[0,1,2].map(i => mean(a.map(x => x.result.modules[i]))),cmd:median(a.map(x => x.result.CMD))};
    });
    const body = summary.map(s => `<tr><td>${esc(s.group)}</td><td>${s.n}</td><td>${fmt(s.age,1)}</td><td>${s.female}/${s.n}</td>${s.m.map(v=>`<td>${fmt(v,3)}</td>`).join('')}<td>${fmt(s.cmd,2)}</td></tr>`).join('');
    const diff = summary.length === 2 ? `<tr><td>差值：${esc(groups[1])} − ${esc(groups[0])}</td><td>—</td><td>${fmt(summary[1].age-summary[0].age,1)}</td><td>—</td>${[0,1,2].map(i=>`<td>${fmt(summary[1].m[i]-summary[0].m[i],3)}</td>`).join('')}<td>${fmt(summary[1].cmd-summary[0].cmd,2)}</td></tr>` : '';
    return `<div class="table-scroll"><table class="compare-table"><thead><tr><th>组别</th><th>n</th><th>平均年龄</th><th>女性</th><th>M1均值</th><th>M2均值</th><th>M3均值</th><th>CMD中位数</th></tr></thead><tbody>${body}${diff}</tbody></table></div>`;
  }
  function matchedComparison(rows, groups, demo) {
    if (demo || groups.length !== 2 || rows.some(r=>!r.study)) return {html:'<p class="form-note">提供完整study列后，可进行同研究、同性别、近年龄配对。合成演示数据不运行显著性检验。</p>',count:0};
    const reference=rows.filter(r=>r.group===groups[0]);
    const comparison=rows.filter(r=>r.group===groups[1]).sort((a,b)=>a.age-b.age);
    const used=new Set(), pairs=[];
    for(const item of comparison){
      let best=null,bestGap=Infinity;
      for(const candidate of reference){
        if(used.has(candidate.sample_id)||candidate.study!==item.study||candidate.sex!==item.sex||(candidate.age<18)!==(item.age<18))continue;
        const gap=Math.abs(candidate.age-item.age),maxGap=item.age<18||candidate.age<18?2:5;
        if(gap<=maxGap&&gap<bestGap){best=candidate;bestGap=gap;}
      }
      if(best){used.add(best.sample_id);pairs.push({case:item,control:best,study:item.study,ageGap:bestGap,d:[0,1,2].map(i=>item.result.modules[i]-best.result.modules[i]).concat([item.result.CMD-best.result.CMD])});}
    }
    if(!pairs.length)return {html:'<p class="form-note">未找到同研究、同性别、近年龄且参照组不重复使用的配对。表格和图仅作描述。</p>',count:0};
    const studyNames=[...new Set(pairs.map(x=>x.study))];
    const means=[0,1,2,3].map(k=>mean(studyNames.map(s=>mean(pairs.filter(p=>p.study===s).map(p=>p.d[k])))));
    const n=pairs.length, studyN=studyNames.length,ageGap=mean(pairs.map(p=>p.ageGap));
    if(n<30) return {html:`<p class="form-note">已找到${n}个配对、${studyN}项研究，平均年龄差${fmt(ageGap,1)}岁。配对数不足30，暂只展示研究等权差值：M1 ${fmt(means[0])}，M2 ${fmt(means[1])}，M3 ${fmt(means[2])}，CMD ${fmt(means[3],2)}。</p>`,count:n};
    // Pair-label sign permutation; study-equal statistic. Fixed seed for reproducibility.
    let state=20261026;const rand=()=>{state^=state<<13;state^=state>>>17;state^=state<<5;return(state>>>0)/4294967296;};
    const byStudy=studyNames.map(s=>pairs.filter(p=>p.study===s));
    const exceed=[0,0,0,0],iterations=10000;
    for(let b=0;b<iterations;b++){
      const draw=byStudy.map(ps=>{
        const sums=[0,0,0,0];
        for(const p of ps){const sign=rand()<.5?-1:1;for(let k=0;k<4;k++)sums[k]+=sign*p.d[k];}
        return sums.map(v=>v/ps.length);
      });
      for(let k=0;k<4;k++)if(Math.abs(mean(draw.map(v=>v[k])))>=Math.abs(means[k])-1e-12)exceed[k]++;
    }
    const p=exceed.map(v=>(v+1)/(iterations+1));
    const order=[0,1,2].sort((a,b)=>p[a]-p[b]),q=[null,null,null];let smallest=1;
    for(let j=2;j>=0;j--){const i=order[j];smallest=Math.min(smallest,p[i]*3/(j+1));q[i]=smallest;}
    const labels=['M1','M2','M3','CMD'];
    const table=labels.map((label,k)=>`<tr><td>${label}</td><td>${fmt(means[k],k===3?2:3)}</td><td>${fmt(p[k],4)}</td><td>${k===3?'单一预设检验':fmt(q[k],4)}</td></tr>`).join('');
    const note=studyN<3?'仅有少于3项研究，检验只说明已上传配对内部的偏移，不支持跨研究方向推广。':'研究数达到3项；仍需逐项留一研究与独立队列验证。';
    return {html:`<div class="matched-block"><h4>同背景配对：${esc(groups[1])} − ${esc(groups[0])}</h4><p>${n}对 · ${studyN}项研究 · 平均年龄差${fmt(ageGap,1)}岁。每项研究先取配对均值，再研究等权汇总。</p><div class="table-scroll"><table class="compare-table"><thead><tr><th>维度</th><th>配对差</th><th>置换P</th><th>病种/组内q</th></tr></thead><tbody>${table}</tbody></table></div><p>${note} 10,000次配对符号置换；M1–M3在本次比较内部BH校正，CMD报告未校正P。</p></div>`,count:n};
  }
  function cohortTSV() {
    const head = ['sample_id','group','study','age','sex','M1','M2','M3','common_axis','differentiation_axis','CMD','age_support','recognized_species','warnings'];
    const textCell = v => { const t=String(v).replace(/[\t\r\n]/g,' '); return /^[\s]*[=+@-]/.test(t) ? "'"+t : t; };
    const rows = cohortResults.map(x => [textCell(x.sample_id),textCell(x.group),textCell(x.study||''),x.age,x.sex,...x.result.modules,x.result.common_axis,x.result.differentiation_axis,x.result.CMD,x.result.age_support,x.result.recognized_species,x.result.warnings.join(' | ')].map(v=>String(v).replace(/[\t\r\n]/g,' ')).join('\t'));
    return [head.join('\t'),...rows].join('\n')+'\n';
  }
  function renderCohort(rows, demo = false) {
    let groups = [...new Set(rows.map(x => x.group))];
    if (groups.length > 2) throw Error('本版两组比较最多支持两个组别；请筛选后重新上传。');
    const requested = demo ? '' : $('#reference-group').value.trim();
    if (requested && !groups.includes(requested)) throw Error('参照组名称与文件group列不一致：' + requested);
    if (requested) groups = [requested, ...groups.filter(g => g !== requested)];
    const ages = rows.map(x => x.age);
    const agesByGroup = groups.map(g=>rows.filter(x=>x.group===g).map(x=>x.age));
    const imbalance = groups.length === 2 && Math.abs(mean(agesByGroup[0])-mean(agesByGroup[1])) > 5;
    const matched = matchedComparison(rows,groups,demo);
    $('#cohort-result').innerHTML = `<div class="cohort-head"><div><span class="eyebrow">GROUP-LEVEL ECOLOGICAL COMPARISON</span><h3>${rows.length}份样本 · ${groups.length}个组</h3><p>模型 ${esc(model.version)} · 年龄范围 ${fmt(Math.min(...ages),1)}–${fmt(Math.max(...ages),1)} 岁</p></div><div class="legend">${groups.map((g,i)=>`<span><i style="background:${i?PINK:BLUE}"></i>${esc(g)}</span>`).join('')}</div></div>
      <div class="chart-grid"><div class="chart-box"><h4>M1 × M2 样本位置</h4>${svgScatter(rows,groups)}</div><div class="chart-box"><h4>联合距离 CMD · 原始点</h4>${svgCMD(rows,groups)}</div></div>
      ${cohortTable(rows,groups)}${matched.html}
      <ul class="warn-list"><li>“${esc(groups[0])}”作为参照，配对效应定义为另一组减参照组。</li><li>上表全组均值差是描述性比较；配对检验只使用找到的${matched.count}对，不能作为个人疾病诊断或发生风险。</li>${imbalance?'<li>两个组的平均年龄差超过5岁；请优先解释同背景配对结果。</li>':''}${demo?'<li>这是人工构造的演示队列；差异由演示数据生成，不是研究发现。</li>':''}<li>年龄段参考人数少于100的样本：${rows.filter(x=>x.result.age_support<100).length}；匹配物种记录少于100的样本：${rows.filter(x=>x.result.recognized_species<100).length}，该数量本身不代表谱不完整。</li></ul>
      <div class="result-actions"><button id="download-cohort-result" type="button">下载全部样本结果 TSV ↓</button><button id="download-cohort-json" type="button">下载结果 JSON ↓</button></div>`;
    $('#download-cohort-result').onclick = () => download(cohortTSV(),'micronorm-cohort-scores.tsv','text/tab-separated-values');
    $('#download-cohort-json').onclick = () => download(JSON.stringify(cohortResults.map(x=>({sample_id:x.sample_id,group:x.group,age:x.age,sex:x.sex,...Object.fromEntries(Object.entries(x.result).filter(([k])=>k!=='rows'))})),null,2),'micronorm-cohort-scores.json','application/json');
  }
  async function demoCohort() {
    try {
      const demos = window.MicroNormDemos;
      if(!Array.isArray(demos)) throw Error('离线演示数据不可用。');
      const base = demos.find(x => x.id === 'adult') || demos[0];
      const members = model.modules;
      const variant = (m1,m2,seed) => {
        const a = {...base.abundance};
        members[0].forEach(j=>{const key=model.species[j];if(a[key])a[key]*=m1*(1+seed*.01);});
        members[1].forEach(j=>{const key=model.species[j];if(a[key])a[key]*=m2*(1-seed*.01);});
        return normalizedAbundance(a);
      };
      const samples = [
        {sample_id:'DEMO-A1',age:45,sex:'female',group:'reference-like',abundance:variant(1,1,0)},
        {sample_id:'DEMO-A2',age:45,sex:'female',group:'reference-like',abundance:variant(1.04,.97,1)},
        {sample_id:'DEMO-B1',age:45,sex:'female',group:'shifted-demo',abundance:variant(1.9,.56,2)},
        {sample_id:'DEMO-B2',age:45,sex:'female',group:'shifted-demo',abundance:variant(1.7,.62,3)}
      ];
      cohortResults = samples.map(s=>({sample_id:s.sample_id,group:s.group,study:'DEMO',age:s.age,sex:s.sex,result:C.score(model,s.age,s.sex,s.abundance)}));
      renderCohort(cohortResults,true);
      setMessage('cohort','已载入同年龄、同性别的人工演示队列；不代表真实疾病结果。');
    } catch(error) { cohortResults=[]; $('#cohort-result').innerHTML='<div class="empty-state"><h3>未生成结果</h3><p>请检查输入后重试。</p></div>'; setMessage('cohort',error.message,true); }
  }
  async function init() {
    try {
      model = window.MicroNormModel;
      C.validateModel(model);
      $('#load-status').textContent = '模型就绪 · '+model.version;
    } catch(error) { $('#load-status').textContent = '模型不可用'; $('#load-status').classList.add('error'); setMessage('single',error.message,true); setMessage('cohort',error.message,true); }
  }
  document.addEventListener('DOMContentLoaded',()=>{
    init();
    $('#single-form').addEventListener('submit',async event=>{
      event.preventDefault(); if(!model)return;
      const file=$('#single-file').files[0];if(!file){setMessage('single','请先选择物种谱文件。',true);return;}
      try {singleResult=C.score(model,$('#single-age').value,$('#single-sex').value,C.parseSingle(await file.text()));renderSingle(singleResult);setMessage('single','计算已完成；样本文件未离开本地浏览器。');}
      catch(error){singleResult=null;$('#single-result').innerHTML='<div class="empty-state"><h3>未生成结果</h3><p>请检查输入后重试。</p></div>';setMessage('single',error.message,true);}
    });
    $('#cohort-form').addEventListener('submit',async event=>{
      event.preventDefault();if(!model)return;
      const file=$('#cohort-file').files[0];if(!file){setMessage('cohort','请先选择队列文件。',true);return;}
      try{const samples=C.parseCohort(await file.text());cohortResults=samples.map(s=>({sample_id:s.sample_id,group:s.group,study:s.study,age:s.age,sex:s.sex,result:C.score(model,s.age,s.sex,s.abundance)}));renderCohort(cohortResults);setMessage('cohort',`已完成${samples.length}份样本的本地评分。`);}
      catch(error){cohortResults=[];$('#cohort-result').innerHTML='<div class="empty-state"><h3>未生成结果</h3><p>请检查输入后重试。</p></div>';setMessage('cohort',error.message,true);}
    });
    $('#demo-single').onclick=demoSingle;$('#demo-cohort').onclick=demoCohort;
    $('#download-single-template').onclick=()=>download('species,abundance_percent\ns__Faecalibacterium_prausnitzii,2.0\ns__Akkermansia_muciniphila,0.5\n','single-sample-format.csv','text/csv');
    $('#download-cohort-template').onclick=()=>download('sample_id,age,sex,group,study,species,abundance_percent\nS01,45,female,group-A,study-01,s__Faecalibacterium_prausnitzii,2.0\nS01,45,female,group-A,study-01,s__Akkermansia_muciniphila,0.5\nS02,45,female,group-B,study-01,s__Faecalibacterium_prausnitzii,1.4\nS02,45,female,group-B,study-01,s__Akkermansia_muciniphila,0.8\n','cohort-long-format.csv','text/csv');
  });
})();

