/* MicroNorm research-preview scoring engine. Percentiles and CMD are not clinical probabilities. */
(function (global) {
  'use strict';
  function positiveRank(x, q, levels) {
    if (x <= q[0]) return clamp(levels[0] * x / Math.max(q[0], .001), .001, .999);
    for (let k = 0; k < levels.length - 1; k++) {
      if (x <= q[k + 1]) return clamp(levels[k] + (levels[k + 1] - levels[k]) * (x - q[k]) / Math.max(q[k + 1] - q[k], .001), .001, .999);
    }
    return clamp(levels[levels.length - 1] + (1 - levels[levels.length - 1]) *
      (1 - Math.exp(-(x - q[q.length - 1]) / Math.max(q[q.length - 1] - q[q.length - 2], .05))), .001, .999);
  }
  function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }
  function speciesName(raw) {
    const name = String(raw || '').trim();
    if (!name) throw Error('物种名称不能为空。');
    return name.startsWith('s__') ? name : 's__' + name;
  }
  function abundanceValue(raw, label) {
    if (raw === null || raw === undefined || String(raw).trim() === '') throw Error('丰度不能为空：' + label);
    const n = Number(raw);
    if (!Number.isFinite(n) || n < 0 || n > 100) throw Error('丰度必须在0–100之间，单位为百分比：' + label);
    return n;
  }
  function normalizeSex(value) {
    const s = String(value).trim().toLowerCase();
    if (['female', 'f', 'woman', '女', '1'].includes(s)) return 'female';
    if (['male', 'm', 'man', '男', '0'].includes(s)) return 'male';
    throw Error('登记性别需为 female/male（或 1/0）；当前预览模型未验证缺失性别评分。');
  }
  function validateModel(model) {
    if (!model || !Array.isArray(model.species) || model.species.length !== 801 ||
      !Array.isArray(model.modules) || model.modules.length !== 3 || !Array.isArray(model.age_counts)) {
      throw Error('模型文件缺少冻结的801物种、三模块或年龄支持数据。');
    }
    return true;
  }
  function score(model, ageInput, sexInput, abundance) {
    validateModel(model);
    if (ageInput === null || ageInput === undefined || String(ageInput).trim() === '') throw Error('年龄不能为空。');
    const age = Number(ageInput), sex = normalizeSex(sexInput);
    if (!Number.isFinite(age) || age < 0 || age > 100) throw Error('年龄必须为0–100岁；其他年龄超出本版模型范围。');
    if (!abundance || typeof abundance !== 'object') throw Error('缺少物种丰度谱。');
    const s = sex === 'male' ? 0 : 1;
    const index = new Map(model.species.map((name, i) => [name, i]));
    const x = new Float64Array(model.species.length), unknown = [];
    let recognized = 0, total = 0, allTotal = 0;
    const canonical = new Set();
    for (const [rawName, rawValue] of Object.entries(abundance)) {
      const name = speciesName(rawName);
      if (canonical.has(name)) throw Error('输入中物种名称重复：' + name);
      canonical.add(name);
      const v = abundanceValue(rawValue, rawName); allTotal += v;
      const j = index.get(name);
      if (j === undefined) { unknown.push(rawName); continue; }
      x[j] = v; recognized++; total += v;
    }
    if (!recognized) throw Error('输入中没有匹配冻结801物种清单的物种。');
    if (!(total > 0)) throw Error('匹配到的物种丰度全部为零，无法判断完整物种谱。');
    if (allTotal > 100.5) throw Error('输入的物种丰度总和超过100.5%；请检查是否混入属级/菌株级条目或单位不是百分比。');
    let lo = model.mids.findIndex(v => v > age) - 1;
    if (lo < 0) lo = age >= model.mids[model.mids.length - 1] ? model.mids.length - 1 : 0;
    const hi = Math.min(lo + 1, model.mids.length - 1);
    const t = hi === lo ? 0 : clamp((age - model.mids[lo]) / (model.mids[hi] - model.mids[lo]), 0, 1);
    const centiles = new Float64Array(x.length), rows = [];
    for (let j = 0; j < x.length; j++) {
      const at = b => {
        const positives = model.pos[s][b][j];
        const p = (positives + model.detection_prior * model.pg[j]) /
          (model.n[s][b][0] + model.detection_prior);
        const w = positives / (positives + model.positive_prior);
        const q = model.locq[s][b][j].map((v, k) => w * v + (1 - w) * model.qg[j][k]);
        return [p, q];
      };
      const [p0, q0] = at(lo), [p1, q1] = at(hi);
      const p = p0 * (1 - t) + p1 * t;
      const q = q0.map((v, k) => v * (1 - t) + q1[k] * t);
      const c = x[j] > 0 ? 1 - p + p * positiveRank(Math.log1p(x[j]), q, model.levels) : (1 - p) / 2;
      centiles[j] = c;
      rows.push({species: model.species[j], abundance_percent: x[j], detection_probability: p, centile: c});
    }
    const modules = model.modules.map(members => members.reduce((sum, j) => sum + centiles[j], 0) / members.length);
    const delta = modules.map((v, i) => v - model.reference_mean[i]);
    let squared = 0;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) squared += delta[i] * model.inverse_covariance[i][j] * delta[j];
    let bin = model.age_edges.findIndex(v => v > age) - 1;
    if (bin < 0) bin = model.age_counts.length - 1;
    bin = clamp(bin, 0, model.age_counts.length - 1);
    const support = model.age_counts[bin];
    const warnings = [];
    if (support < 100) warnings.push('此年龄段健康参考不足100人；结果对队列构成敏感。');
    else if (support < 300) warnings.push('此年龄段健康参考不足300人，建议结合不确定性阅读。');
    if (recognized < 100) warnings.push('匹配到的物种记录少于100个；该数量本身不代表谱不完整，请核对输入是否覆盖所有已检出的物种。');
    if (unknown.length) warnings.push(`${unknown.length}个输入物种不在冻结801物种清单中。`);
    warnings.push('国家、测序平台及独立临床域尚无针对本样本的外部支持判定。');
    return {
      model_version: model.version, model_status: 'empirical-bin research preview', age, sex,
      modules, common_axis: (modules[0] + modules[1]) / Math.SQRT2,
      differentiation_axis: (modules[0] - modules[1]) / Math.SQRT2,
      CMD: Math.sqrt(Math.max(0, squared)), age_support: support, recognized_species: recognized,
      unknown_species: unknown, recognized_abundance_sum: total,
      exploratory_tail_count: rows.filter(r => r.abundance_percent > 0 && r.detection_probability >= .1 && (r.centile < .05 || r.centile > .95)).length,
      exploratory_tail_evaluable: rows.filter(r => r.abundance_percent > 0 && r.detection_probability >= .1).length,
      warnings, rows
    };
  }
  function parseDelimited(text) {
    const lines = text.replace(/^\uFEFF/, '').split(/\r?\n/).filter(line => line.trim());
    if (!lines.length) throw Error('输入文件为空。');
    const delimiter = lines[0].includes('\t') ? '\t' : ',';
    const parseLine = line => {
      const cells = []; let value = '', quoted = false;
      for (let i = 0; i < line.length; i++) {
        const ch = line[i];
        if (ch === '"') { if (quoted && line[i + 1] === '"') { value += '"'; i++; } else quoted = !quoted; }
        else if (ch === delimiter && !quoted) { cells.push(value.trim()); value = ''; }
        else value += ch;
      }
      if (quoted) throw Error('CSV引号未闭合。');
      cells.push(value.trim()); return cells;
    };
    const table = lines.map(parseLine);
    const header = table.shift().map(x => x.trim().toLowerCase());
    if (new Set(header).size !== header.length) throw Error('输入文件存在重复列名。');
    return {header, records: table.map((cells, row) => {
      if (cells.length !== header.length) throw Error(`第${row + 2}行列数与表头不符。`);
      return Object.fromEntries(header.map((name, i) => [name, cells[i]]));
    })};
  }
  function parseSingle(text) {
    const lines = text.replace(/^\uFEFF/, '').split(/\r?\n/);
    const native = lines.findIndex(line => /^#?clade_name\t/i.test(line) && /relative_abundance/i.test(line));
    let parsed;
    if (native >= 0) {
      const header = lines[native].replace(/^#/, '');
      parsed = parseDelimited([header, ...lines.slice(native + 1).filter(line => line.trim() && !line.startsWith('#'))].join('\n'));
    } else parsed = parseDelimited(text);
    const {header, records} = parsed;
    const isNative = header.includes('clade_name') && header.includes('relative_abundance');
    if (!isNative && (!header.includes('species') || !header.includes('abundance_percent'))) throw Error('单样本文件需包含species和abundance_percent，或使用MetaPhlAn原生物种丰度表。');
    const abundance = {};
    for (const row of records) {
      let name = row.species;
      if (isNative) {
        const leaf = String(row.clade_name || '').split('|').pop().trim();
        if (!leaf.startsWith('s__')) continue;
        name = leaf;
      }
      if (!name) continue;
      const key = speciesName(name);
      if (Object.prototype.hasOwnProperty.call(abundance, key)) throw Error('文件有重复物种：' + key);
      abundance[key] = abundanceValue(isNative ? row.relative_abundance : row.abundance_percent, key);
    }
    return abundance;
  }
  function parseCohort(text) {
    const {header, records} = parseDelimited(text);
    const required = ['sample_id', 'age', 'sex', 'group', 'species', 'abundance_percent'];
    for (const col of required) if (!header.includes(col)) throw Error('队列文件缺少列：' + col);
    const samples = new Map();
    for (const row of records) {
      if (!row.sample_id || !row.species || !row.group) throw Error('sample_id、group和species不能为空。');
      if (String(row.age).trim() === '') throw Error('年龄不能为空：' + row.sample_id);
      const age = Number(row.age), sex = normalizeSex(row.sex);
      if (!Number.isFinite(age) || age < 0 || age > 100) throw Error('无效年龄：' + row.sample_id);
      const v = abundanceValue(row.abundance_percent, row.sample_id);
      const study = (row.study || '').trim();
      let sample = samples.get(row.sample_id);
      if (!sample) { sample = {sample_id: row.sample_id, age, sex, group: row.group, study, abundance: {}}; samples.set(row.sample_id, sample); }
      if (sample.age !== age || sample.sex !== sex || sample.group !== row.group || sample.study !== study) throw Error('同一样本的年龄、性别、组别或研究编号不一致：' + row.sample_id);
      const key = speciesName(row.species);
      if (Object.prototype.hasOwnProperty.call(sample.abundance, key)) throw Error('同一样本有重复物种：' + row.sample_id + ' / ' + key);
      sample.abundance[key] = v;
    }
    if (samples.size > 250) throw Error('本地预览一次最多处理250份样本，请分批上传。');
    if (!samples.size) throw Error('队列文件没有有效样本。');
    return Array.from(samples.values());
  }
  global.MicroNormCore = {score, parseSingle, parseCohort, normalizeSex, validateModel};
})(window);

