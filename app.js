/* ============================================================
 * A股股票价格周期分析 —— 傅里叶变换 (FFT) 套利周期
 * 纯前端实现，数据源：东方财富 K 线接口（支持跨域）
 * ============================================================ */

'use strict';

/* ---------- 工具函数 ---------- */

const $ = (id) => document.getElementById(id);

function nextPow2(n) {
  let p = 1;
  while (p < n) p <<= 1;
  return p;
}

function mean(arr) {
  let s = 0;
  for (let i = 0; i < arr.length; i++) s += arr[i];
  return s / arr.length;
}

function hannWindow(n) {
  const w = new Array(n);
  if (n === 1) { w[0] = 1; return w; }
  for (let i = 0; i < n; i++) w[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (n - 1)));
  return w;
}

/* 迭代式基-2 快速傅里叶变换（Cooley-Tukey） */
function fft(re, im) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * Math.PI) / len;
    const wr = Math.cos(ang);
    const wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wR = 1, wI = 0;
      for (let j = 0; j < (len >> 1); j++) {
        const a = i + j, b = i + j + (len >> 1);
        const uR = re[a], uI = im[a];
        const vR = re[b] * wR - im[b] * wI;
        const vI = re[b] * wI + im[b] * wR;
        re[a] = uR + vR; im[a] = uI + vI;
        re[b] = uR - vR; im[b] = uI - vI;
        const nR = wR * wr - wI * wi;
        wI = wR * wi + wI * wr;
        wR = nR;
      }
    }
  }
}

/* 线性去趋势：返回斜率/截距/趋势/残差 */
function detrendLinear(y) {
  const n = y.length;
  const xm = (n - 1) / 2;
  const ym = mean(y);
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    const x = i - xm;
    sxx += x * x;
    sxy += x * (y[i] - ym);
  }
  const slope = sxy / sxx;
  const intercept = ym - slope * xm;
  const trend = new Array(n), resid = new Array(n);
  for (let i = 0; i < n; i++) {
    trend[i] = slope * i + intercept;
    resid[i] = y[i] - trend[i];
  }
  return { slope, intercept, trend, resid };
}

/* 滞后1自相关 */
function lag1Autocorr(x) {
  const n = x.length;
  const m = mean(x);
  let num = 0, den = 0;
  for (let i = 0; i < n - 1; i++) num += (x[i] - m) * (x[i + 1] - m);
  for (let i = 0; i < n; i++) den += (x[i] - m) * (x[i] - m);
  if (den === 0) return 0;
  return num / den;
}

/* 解 3x3 线性方程组（高斯消元） */
function solve3(A, b) {
  const a = A.map((r) => r.slice());
  const v = b.slice();
  for (let i = 0; i < 3; i++) {
    let p = i;
    for (let r = i + 1; r < 3; r++) if (Math.abs(a[r][i]) > Math.abs(a[p][i])) p = r;
    [a[i], a[p]] = [a[p], a[i]];
    [v[i], v[p]] = [v[p], v[i]];
    const piv = a[i][i];
    if (Math.abs(piv) < 1e-15) continue;
    for (let r = i + 1; r < 3; r++) {
      const f = a[r][i] / piv;
      for (let c = i; c < 3; c++) a[r][c] -= f * a[i][c];
      v[r] -= f * v[i];
    }
  }
  const x = [0, 0, 0];
  for (let i = 2; i >= 0; i--) {
    let s = v[i];
    for (let c = i + 1; c < 3; c++) s -= a[i][c] * x[c];
    x[i] = (Math.abs(a[i][i]) < 1e-15) ? 0 : s / a[i][i];
  }
  return x;
}

/* 对给定周期 T（天）做最小二乘正弦拟合（作用于去趋势序列）
   返回 {A, B, C, amplitude, phase, r2} */
function fitSinusoid(t, y, T) {
  const n = y.length;
  const w = (2 * Math.PI) / T;
  let scc = 0, scs = 0, sss = 0, sc = 0, ss = 0, syc = 0, sys = 0, sy = 0;
  for (let i = 0; i < n; i++) {
    const c = Math.cos(w * t[i]);
    const s = Math.sin(w * t[i]);
    scc += c * c; scs += c * s; sss += s * s;
    sc += c; ss += s;
    syc += y[i] * c; sys += y[i] * s; sy += y[i];
  }
  const A = [
    [scc, scs, sc],
    [scs, sss, ss],
    [sc, ss, n],
  ];
  const sol = solve3(A, [syc, sys, sy]);
  const amplitude = Math.hypot(sol[0], sol[1]);
  const phase = Math.atan2(sol[1], sol[0]);
  // R^2
  const ym = mean(y);
  let sst = 0, ssr = 0;
  for (let i = 0; i < n; i++) {
    const fit = sol[0] * Math.cos(w * t[i]) + sol[1] * Math.sin(w * t[i]) + sol[2];
    ssr += (y[i] - fit) * (y[i] - fit);
    sst += (y[i] - ym) * (y[i] - ym);
  }
  const r2 = sst > 0 ? Math.max(0, 1 - ssr / sst) : 0;
  return { A: sol[0], B: sol[1], C: sol[2], amplitude, phase, r2 };
}

/* ---------- 数据获取 ---------- */

function detectMarket(code) {
  const c = code.trim()[0];
  if (c === '6' || c === '5' || c === '9') return '1'; // 沪市
  return '0'; // 深市 / 北交
}

/* 根据代码前缀识别板块 */
function detectBoard(code) {
  const c = String(code).trim();
  if (/^688|^689/.test(c)) return '科创板';
  if (/^300|^301|^302/.test(c)) return '创业板';
  if (/^600|^601|^603|^605/.test(c)) return '沪市主板';
  if (/^000|^001|^002|^003/.test(c)) return '深市主板';
  if (/^430|^43|^83|^87|^88|^920/.test(c)) return '北交所';
  if (/^5/.test(c)) return '沪市基金';
  if (/^1/.test(c)) return '深市基金';
  return 'A股';
}

/* 通用 JSONP 请求：script 标签注入，绕过跨域(CORS)限制。
   国内浏览器对东方财富接口的 CORS 头支持不稳定，fetch 常报 "Failed to fetch"；
   改用 JSONP（东方财富接口原生支持 cb 回调）后，任何站点/任何网络都能稳定取数。 */
function jsonp(url, timeoutMs) {
  return new Promise((resolve, reject) => {
    const cbName = '__emCb_' + Date.now() + '_' + Math.floor(Math.random() * 1e6);
    const fullUrl = url + (url.indexOf('?') >= 0 ? '&' : '?') + 'cb=' + cbName;
    const script = document.createElement('script');
    let done = false;
    const cleanup = () => {
      try { delete window[cbName]; } catch (e) {}
      if (script.parentNode) script.parentNode.removeChild(script);
    };
    window[cbName] = (data) => {
      if (done) return;
      done = true;
      cleanup();
      resolve(data);
    };
    script.onerror = () => {
      if (done) return;
      done = true;
      cleanup();
      reject(new Error('网络请求失败'));
    };
    script.src = fullUrl;
    document.head.appendChild(script);
    setTimeout(() => {
      if (!done) {
        done = true;
        cleanup();
        reject(new Error('请求超时'));
      }
    }, timeoutMs || 10000);
  });
}

async function fetchKline(code, market, lmt) {
  const secid = market + '.' + code;
  // 将“N 个交易日”换算为日历日（约 1.7 倍）并留 15 天冗余，覆盖节假日/停牌
  const calDays = Math.ceil(lmt * 1.7) + 15;
  const d = new Date();
  d.setDate(d.getDate() - calDays);
  const beg = d.toISOString().slice(0, 10).replace(/-/g, '');
  const url = 'https://push2his.eastmoney.com/api/qt/stock/kline/get' +
    '?secid=' + secid +
    '&fields1=f1,f2,f3,f4,f5,f6' +
    '&fields2=f51,f52,f53,f54,f55,f56,f57,f58,f59,f60,f61' +
    '&klt=101&fqt=1&beg=' + beg + '&end=20500101';
  const json = await jsonp(url, 10000);
  if (!json || json.rc !== 0 || !json.data) throw new Error('接口返回异常');
  return json.data;
}

async function loadStock(code, marketChoice, lmt) {
  code = code.trim();
  if (!/^\d{6}$/.test(code)) throw new Error('请输入6位股票代码');

  let markets = [marketChoice];
  if (marketChoice === 'auto') {
    markets = [detectMarket(code)];
  }
  // 尝试候选市场，失败则尝试另一市场
  const tried = [];
  let data = null, lastErr = null;
  for (const m of markets) {
    if (tried.includes(m)) continue;
    tried.push(m);
    try {
      data = await fetchKline(code, m, lmt);
      if (data && data.klines && data.klines.length) break;
    } catch (e) { lastErr = e; }
  }
  if (!data || !data.klines || !data.klines.length) {
    // 自动模式用另一市场再试一次
    if (marketChoice === 'auto') {
      const alt = markets[0] === '1' ? '0' : '1';
      try { data = await fetchKline(code, alt, lmt); } catch (e) { lastErr = e; }
    }
  }
  if (!data || !data.klines || !data.klines.length) {
    throw new Error('未获取到行情数据，请确认股票代码是否正确（' + (lastErr ? lastErr.message : '') + '）');
  }

  // 解析 kline：date,open,close,high,low,volume,amount,振幅,涨跌幅,涨跌额,换手率
  const allDates = [], allCloses = [], allOpens = [], allHighs = [], allLows = [];
  for (const line of data.klines) {
    const p = line.split(',');
    allDates.push(p[0]);
    allOpens.push(parseFloat(p[1]));
    allCloses.push(parseFloat(p[2]));
    allHighs.push(parseFloat(p[3]));
    allLows.push(parseFloat(p[4]));
  }
  // 截取最近 lmt 个交易日
  const n = allCloses.length;
  const start = Math.max(0, n - lmt);
  const slice = (a) => a.slice(start);
  return {
    name: data.name,
    code: data.code,
    market: data.market,
    preClose: data.preKPrice,
    dates: slice(allDates), closes: slice(allCloses),
    opens: slice(allOpens), highs: slice(allHighs), lows: slice(allLows),
  };
}

/* ---------- 股票搜索（拼音首字母 / 中文 / 代码） ---------- */

/* 搜索接口不支持跨域，改用 JSONP（script 标签注入）绕过 CORS */
function searchSuggestJSONP(keyword) {
  const url = 'https://searchapi.eastmoney.com/api/suggest/get?input=' +
    encodeURIComponent(keyword) + '&type=14&count=8';
  return jsonp(url, 8000);
}

async function searchSuggest(keyword) {
  const data = await searchSuggestJSONP(keyword);
  const list = (data && data.QuotationCodeTable && data.QuotationCodeTable.Data) || [];
  // 仅保留 A 股：市场前缀 0(深A/北交) 或 1(沪A)，代码为 6 位数字，排除港股/美股等
  return list
    .filter((x) => {
      const m = (x.QuoteID || '').split('.')[0];
      return /^\d{6}$/.test(x.Code) && (m === '0' || m === '1');
    })
    .map((x) => ({
      code: x.Code,
      name: x.Name,
      pinyin: x.PinYin,
      market: (x.QuoteID || '').split('.')[0],
      board: detectBoard(x.Code),
    }));
}

/* 输入联想下拉 */
let suggestTimer = null;
let suggestItems = [];

function hideSuggest() {
  const box = $('suggestList');
  if (box) { box.classList.remove('show'); box.innerHTML = ''; }
  suggestItems = [];
  activeIdx = -1;
}

function renderSuggest(items) {
  suggestItems = items;
  activeIdx = -1;
  const box = $('suggestList');
  if (!items.length) { hideSuggest(); return; }
  box.innerHTML = items.map((it, i) =>
    '<div class="suggest-item" data-i="' + i + '">' +
      '<span class="s-name">' + it.name + '</span>' +
      '<span class="s-code">' + it.code + '</span>' +
      '<span class="s-board">' + it.board + '</span>' +
    '</div>').join('');
  box.classList.add('show');
}

function onInput() {
  const v = $('codeInput').value.trim();
  if (/^\d{6}$/.test(v) || v.length === 0) { hideSuggest(); return; }
  clearTimeout(suggestTimer);
  suggestTimer = setTimeout(async () => {
    try {
      const items = await searchSuggest(v);
      // 若用户已继续输入，忽略过期结果
      if ($('codeInput').value.trim() !== v) return;
      renderSuggest(items);
    } catch (e) { hideSuggest(); }
  }, 220);
}

function selectSuggestion(it) {
  $('codeInput').value = it.code;
  hideSuggest();
  run({ code: it.code, market: it.market, name: it.name });
}

/* ---------- 核心分析：FFT 找周期 ---------- */

function analyze(closes) {
  const M = closes.length;
  // 对数价格
  const logP = closes.map((c) => Math.log(c));
  // 去趋势
  const { trend, resid } = detrendLinear(logP);
  const residMean = mean(resid);
  for (let i = 0; i < M; i++) resid[i] -= residMean; // 残差零均值

  // 汉宁窗
  const win = hannWindow(M);
  const winSq = win.reduce((a, b) => a + b * b, 0);
  const x = new Array(M);
  for (let i = 0; i < M; i++) x[i] = resid[i] * win[i];

  // 零填充到 N（>= 4M 的 2 的幂，提高频率插值分辨率）
  let N = nextPow2(Math.max(4 * M, 512));

  const re = new Array(N).fill(0);
  const im = new Array(N).fill(0);
  for (let i = 0; i < M; i++) re[i] = x[i];
  fft(re, im);

  // 单边功率谱（k = 1 .. N/2），归一化使之可与红噪声谱直接比较
  const half = N >> 1;
  const power = new Array(half + 1).fill(0);
  for (let k = 1; k <= half; k++) {
    // 归一化：|F|^2 / (N * sum(w^2))，其均值约等于方差
    power[k] = (re[k] * re[k] + im[k] * im[k]) / (N * winSq);
  }

  // AR(1) 红噪声谱 + 95% 显著性阈值
  const r1 = lag1Autocorr(resid);
  const r1c = Math.max(-0.99, Math.min(0.99, r1));
  const pRed = new Array(half + 1).fill(0);
  const thresh = new Array(half + 1).fill(0);
  let sumRed = 0, sumPow = 0;
  for (let k = 1; k <= half; k++) {
    const f = k / N; // cycles/day
    const c = Math.cos(2 * Math.PI * f);
    pRed[k] = (1 - r1c * r1c) / (1 - 2 * r1c * c + r1c * r1c);
    sumRed += pRed[k];
    sumPow += power[k];
  }
  // 将红噪声谱缩放到与观测功率谱同均值
  const scale = sumPow / sumRed;
  for (let k = 1; k <= half; k++) {
    pRed[k] *= scale;
    thresh[k] = pRed[k] * (-Math.log(0.05)); // 95% 置信水平
  }

  // 找出局部极大值（峰）
  const peaks = [];
  for (let k = 2; k < half; k++) {
    if (power[k] > power[k - 1] && power[k] >= power[k + 1]) {
      // 抛物线插值细化峰值（对数功率）
      const lp0 = Math.log(Math.max(power[k - 1], 1e-12));
      const lp1 = Math.log(Math.max(power[k], 1e-12));
      const lp2 = Math.log(Math.max(power[k + 1], 1e-12));
      const denom = lp0 - 2 * lp1 + lp2;
      let delta = 0;
      if (Math.abs(denom) > 1e-12) delta = 0.5 * (lp0 - lp2) / denom;
      const kPeak = k + delta;
      const period = N / kPeak;
      peaks.push({ k: kPeak, period, power: power[k] });
    }
  }
  peaks.sort((a, b) => b.power - a.power);

  // 筛选有效周期：period 在 [3, M/2] 之间
  const valid = peaks.filter((p) => p.period >= 3 && p.period <= M / 2);

  // 对每个有效峰做最小二乘拟合，得到振幅/相位/规律度
  const t = closes.map((_, i) => i);
  const result = valid.slice(0, 12).map((p) => {
    const fit = fitSinusoid(t, resid, p.period);
    const sig = p.power > thresh[Math.round(p.k)] ? true : false;
    const sigLevel = thresh[Math.round(p.k)] || 0;
    return {
      period: p.period,
      power: p.power,
      threshold: sigLevel,
      significant: sig,
      amplitude: fit.amplitude,      // 对数幅度 ≈ 百分比
      amplitudePct: fit.amplitude * 100,
      phase: fit.phase,             // 弧度
      r2: fit.r2,
      cycles: M / p.period,
      A: fit.A, B: fit.B,
    };
  });

  // 选出“最有规律的套利周期”：优先显著、且周期适中、振幅较大
  const candidates = result.filter((r) => r.period >= 4 && r.cycles >= 2);
  const top = (candidates.length ? candidates : result)[0];

  return { M, N, resid, trend, logP, power, thresh, result, top };
}

/* ---------- 周期重建 ---------- */

function buildRecon(logP, trend, result, topK) {
  const M = logP.length;
  const t = logP.map((_, i) => i);
  // 优先选显著的周期用于重建
  const ordered = result.slice().sort((a, b) => (b.significant - a.significant) || (b.power - a.power));
  const chosen = ordered.slice(0, topK);
  const detrendRecon = new Array(M).fill(0);
  const comps = chosen.map((r) => {
    const w = (2 * Math.PI) / r.period;
    const series = t.map((ti) => r.A * Math.cos(w * ti) + r.B * Math.sin(w * ti));
    for (let i = 0; i < M; i++) detrendRecon[i] += series[i];
    return { period: r.period, series };
  });
  const fittedLog = new Array(M);
  const fittedPrice = new Array(M);
  for (let i = 0; i < M; i++) {
    fittedLog[i] = trend[i] + detrendRecon[i];
    fittedPrice[i] = Math.exp(fittedLog[i]);
  }
  return { detrendRecon, fittedPrice, comps };
}

/* ---------- 时点测算 ---------- */

function nextTurningPoints(top, dates) {
  const M = dates.length;
  const T = top.period;
  const phi = top.phase; // 峰值相位
  const w = (2 * Math.PI) / T;
  const tLast = M - 1;
  // 下一个峰值：ωt = phi + 2πk，t > tLast
  let tPeak = ((phi - w * tLast) / w) % T;
  while (tPeak <= 0) tPeak += T;
  // 下一个谷值：ωt = phi + π + 2πk
  let tTrough = ((phi + Math.PI - w * tLast) / w) % T;
  while (tTrough <= 0) tTrough += T;
  const addDays = (iso, days) => {
    const d = new Date(iso + 'T00:00:00');
    d.setDate(d.getDate() + Math.round(days));
    return d.toISOString().slice(0, 10);
  };
  const lastDate = dates[dates.length - 1];
  return {
    peakDays: tPeak,
    troughDays: tTrough,
    peakDate: addDays(lastDate, tPeak),
    troughDate: addDays(lastDate, tTrough),
  };
}

function humanizePeriod(T) {
  const weeks = T / 5, months = T / 21;
  if (T < 8) return T.toFixed(1) + ' 天（约1周）';
  if (weeks < 4) return '约 ' + weeks.toFixed(1) + ' 周';
  if (months < 12) return '约 ' + months.toFixed(1) + ' 个月';
  return '约 ' + (T / 250).toFixed(2) + ' 年';
}

/* ---------- 图表渲染 ---------- */

let priceChart = null, spectrumChart = null, cyclesChart = null;
let chartsReady = false;

function initCharts() {
  if (chartsReady) return;
  if (typeof echarts === 'undefined') {
    throw new Error('图表库 echarts 加载失败，请检查网络后刷新');
  }
  priceChart = echarts.init($('chartPrice'));
  spectrumChart = echarts.init($('chartSpectrum'));
  cyclesChart = echarts.init($('chartCycles'));
  chartsReady = true;
  window.addEventListener('resize', () => {
    priceChart && priceChart.resize();
    spectrumChart && spectrumChart.resize();
    cyclesChart && cyclesChart.resize();
  });
}

const COLORS = ['#3aa0ff', '#ff8a3d', '#3fb950', '#e3b341', '#c084fc'];

function renderPriceChart(dates, closes, fittedPrice) {
  priceChart.setOption({
    backgroundColor: 'transparent',
    tooltip: { trigger: 'axis' },
    legend: { data: ['收盘价', '周期拟合'], textStyle: { color: '#8b98a9' }, top: 0 },
    grid: { left: 55, right: 20, top: 34, bottom: 40 },
    xAxis: {
      type: 'category', data: dates,
      axisLine: { lineStyle: { color: '#2b3542' } },
      axisLabel: { color: '#8b98a9' },
    },
    yAxis: {
      type: 'value', scale: true,
      splitLine: { lineStyle: { color: '#1c232e' } },
      axisLabel: { color: '#8b98a9' },
    },
    dataZoom: [{ type: 'inside' }, { type: 'slider', height: 16, bottom: 8 }],
    series: [
      {
        name: '收盘价', type: 'line', data: closes, smooth: true,
        showSymbol: false, lineStyle: { width: 1.6, color: '#e6edf3' },
        areaStyle: { opacity: 0.06 },
      },
      {
        name: '周期拟合', type: 'line', data: fittedPrice, smooth: true,
        showSymbol: false, lineStyle: { width: 2.4, color: '#ff8a3d', type: 'dashed' },
      },
    ],
  });
}

function renderSpectrumChart(ana) {
  const M = ana.M;
  const half = ana.power.length - 1;
  const data = [], threshData = [];
  for (let k = 1; k <= half; k++) {
    const period = ana.N / k;
    if (period > M || period < 1.9) continue; // 只显示周期 ≤ 数据长度的可信区间
    data.push([period, ana.power[k]]);
    threshData.push([period, ana.thresh[k]]);
  }
  // 仅标注靠前的峰，按显著性着色
  const peaks = ana.result.slice(0, 8).map((r) => ({
    value: [r.period, r.power],
    name: r.period.toFixed(1) + '天',
    itemStyle: { color: r.significant ? '#e3b341' : '#5b6470' },
  }));

  spectrumChart.setOption({
    backgroundColor: 'transparent',
    tooltip: {
      trigger: 'axis',
      formatter: (params) => {
        const p = params[0];
        if (!p) return '';
        return '周期 ≈ ' + Number(p.value[0]).toFixed(1) + ' 天<br/>功率 ' + Number(p.value[1]).toFixed(5);
      },
    },
    legend: { data: ['功率谱', '95%显著性阈值', '周期峰'], textStyle: { color: '#8b98a9' }, top: 0 },
    grid: { left: 55, right: 30, top: 34, bottom: 50 },
    xAxis: {
      type: 'log', logBase: 10, min: 2, max: Math.max(M, 5),
      name: '周期（交易日，对数刻度）', nameLocation: 'middle', nameGap: 32,
      nameTextStyle: { color: '#8b98a9' },
      axisLine: { lineStyle: { color: '#2b3542' } },
      axisLabel: { color: '#8b98a9', formatter: (v) => Math.round(v) + '天' },
      splitLine: { lineStyle: { color: '#1c232e' } },
    },
    yAxis: {
      type: 'log',
      name: '功率', nameTextStyle: { color: '#8b98a9' },
      splitLine: { lineStyle: { color: '#1c232e' } },
      axisLabel: { color: '#8b98a9' },
    },
    series: [
      {
        name: '功率谱', type: 'line', data, showSymbol: false,
        lineStyle: { color: '#3aa0ff', width: 1.2 },
      },
      {
        name: '95%显著性阈值', type: 'line', data: threshData, showSymbol: false,
        lineStyle: { color: '#f85149', width: 1.4 },
      },
      {
        name: '周期峰', type: 'scatter', data: peaks,
        symbolSize: 14, itemStyle: { borderColor: '#08111f', borderWidth: 1 },
        label: { show: true, position: 'top', color: '#e3b341', fontSize: 11, formatter: (p) => p.name },
      },
    ],
  });
}

function renderCyclesChart(dates, resid, comps, detrendRecon) {
  const series = [
    {
      name: '去趋势波动', type: 'line', data: resid, smooth: true,
      showSymbol: false, lineStyle: { width: 1.2, color: '#4b5563' },
    },
  ];
  comps.forEach((c, i) => {
    series.push({
      name: c.period.toFixed(1) + '天分量', type: 'line', data: c.series, smooth: true,
      showSymbol: false, lineStyle: { width: 1.4, color: COLORS[i % COLORS.length], type: 'dashed' },
    });
  });
  series.push({
    name: '分量叠加', type: 'line', data: detrendRecon, smooth: true,
    showSymbol: false, lineStyle: { width: 2.2, color: '#e6edf3' },
  });

  cyclesChart.setOption({
    backgroundColor: 'transparent',
    tooltip: { trigger: 'axis' },
    legend: { type: 'scroll', textStyle: { color: '#8b98a9' }, top: 0 },
    grid: { left: 55, right: 20, top: 34, bottom: 40 },
    xAxis: {
      type: 'category', data: dates,
      axisLine: { lineStyle: { color: '#2b3542' } },
      axisLabel: { color: '#8b98a9' },
    },
    yAxis: {
      type: 'value', scale: true,
      splitLine: { lineStyle: { color: '#1c232e' } },
      axisLabel: { color: '#8b98a9' },
    },
    dataZoom: [{ type: 'inside' }, { type: 'slider', height: 16, bottom: 8 }],
    series,
  });
}

/* ---------- UI 渲染 ---------- */

function fmtSign(n, withPct = true) {
  const sign = n > 0 ? '+' : '';
  return sign + n.toFixed(2) + (withPct ? '%' : '');
}

function renderResults(stock, ana) {
  const { M, resid, trend, logP, result, top } = ana;

  // 股票信息
  const lastClose = stock.closes[stock.closes.length - 1];
  const firstClose = stock.closes[0];
  const chg = (lastClose / firstClose - 1) * 100;
  $('siName').textContent = stock.name;
  $('siBoard').textContent = detectBoard(stock.code);
  $('siCode').textContent = stock.code;
  $('siClose').textContent = lastClose.toFixed(2);
  $('siChg').textContent = fmtSign(chg);
  $('siChg').className = chg >= 0 ? 'up' : 'down';
  $('siCount').textContent = M + ' 天';
  $('siRange').textContent = stock.dates[0] + ' ~ ' + stock.dates[M - 1];
  $('stockInfo').classList.remove('hidden');

  // hero 最优周期
  const tp = nextTurningPoints(top, stock.dates);
  $('hcPeriod').textContent = top.period.toFixed(1);
  $('hcDesc').textContent = humanizePeriod(top.period) + '的周期性波动';
  $('hcAmp').textContent = '±' + top.amplitudePct.toFixed(1) + '%';
  $('hcSig').textContent = top.significant ? '显著（95%）' : '未达显著';
  $('hcSig').style.color = top.significant ? 'var(--green)' : 'var(--muted)';
  $('hcR2').textContent = (top.r2 * 100).toFixed(1) + '%';
  $('hcCycles').textContent = '约 ' + Math.floor(top.cycles) + ' 次';
  $('hcTiming').innerHTML =
    '按该周期模型推算（统计外推，仅供参考）：下一个潜在<b>低点</b>约在 <b>' + tp.troughDate +
    '</b>（约 ' + Math.round(tp.troughDays) + ' 个交易日后），下一个潜在<b>高点</b>约在 <b>' + tp.peakDate +
    '</b>（约 ' + Math.round(tp.peakDays) + ' 个交易日后）。';
  $('heroCycle').classList.remove('hidden');

  // 图表（先显示容器再初始化，避免尺寸为 0）
  $('chartsArea').classList.remove('hidden');
  initCharts();
  const reconObj = buildRecon(logP, trend, result, 3);
  renderPriceChart(stock.dates, stock.closes, reconObj.fittedPrice);
  renderSpectrumChart(ana);
  renderCyclesChart(stock.dates, resid, reconObj.comps, reconObj.detrendRecon);

  // 表格
  const tbody = $('cycleTable').querySelector('tbody');
  tbody.innerHTML = '';
  result.slice(0, 10).forEach((r, i) => {
    const tr = document.createElement('tr');
    const badge = i < 3 ? 'rank-badge rank-' + (i + 1) : 'rank-badge';
    const sig = r.significant ? '<span class="sig-yes">显著</span>' : '<span class="sig-no">未达</span>';
    tr.innerHTML =
      '<td><span class="' + badge + '">' + (i + 1) + '</span></td>' +
      '<td>' + r.period.toFixed(1) + ' 天</td>' +
      '<td>' + humanizePeriod(r.period) + '</td>' +
      '<td>±' + r.amplitudePct.toFixed(1) + '%</td>' +
      '<td>' + (r.r2 * 100).toFixed(1) + '%</td>' +
      '<td>约 ' + Math.floor(r.cycles) + ' 次</td>' +
      '<td>' + sig + '</td>';
    tbody.appendChild(tr);
  });
  $('tableArea').classList.remove('hidden');
}

/* ---------- 主流程 ---------- */

async function run(target) {
  const btn = $('analyzeBtn');
  const status = $('statusBar');
  const lmt = parseInt($('rangeSelect').value, 10);

  btn.disabled = true;
  status.className = 'status-bar loading';
  status.textContent = '⏳ 正在获取行情并计算傅里叶频谱…';
  if (typeof echarts === 'undefined') {
    status.className = 'status-bar err';
    status.textContent = '❌ 图表库 echarts 加载失败，请检查网络后刷新页面重试';
    btn.disabled = false;
    return;
  }
  try {
    let code, market;
    if (target && target.code) {
      // 来自联想下拉选中的股票
      code = target.code;
      market = target.market || 'auto';
      $('codeInput').value = code;
    } else {
      const raw = $('codeInput').value.trim();
      if (/^\d{6}$/.test(raw)) {
        code = raw;
        market = $('marketSelect').value;
      } else {
        if (!raw) throw new Error('请输入股票代码、中文名或拼音首字母');
        status.textContent = '🔍 正在搜索「' + raw + '」…';
        const items = await searchSuggest(raw);
        if (!items.length) throw new Error('未找到「' + raw + '」对应的A股，请换关键词或直接输入6位代码');
        const first = items[0];
        code = first.code;
        market = first.market || 'auto';
        $('codeInput').value = first.code;
        hideSuggest();
      }
    }
    const stock = await loadStock(code, market, lmt);
    if (stock.closes.length < 20) {
      throw new Error('数据量过少（' + stock.closes.length + ' 天），无法进行可靠的周期分析');
    }
    const ana = analyze(stock.closes);
    if (!ana.result || !ana.result.length) {
      throw new Error('未检测到明显周期');
    }
    renderResults(stock, ana);
    status.className = 'status-bar ok';
    status.textContent = '✅ 分析完成：' + stock.name + '（' + stock.code + '）共 ' +
      stock.closes.length + ' 个交易日，检出 ' + ana.result.length + ' 个主要周期。';
  } catch (e) {
    status.className = 'status-bar err';
    status.textContent = '❌ ' + e.message;
  } finally {
    btn.disabled = false;
  }
}

function bindEvents() {
  $('analyzeBtn').addEventListener('click', () => run());
  $('codeInput').addEventListener('input', onInput);
  $('codeInput').addEventListener('keydown', (e) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      // 优先选中键盘高亮的项，其次选第一项，否则直接分析
      if ($('suggestList').classList.contains('show') && suggestItems.length) {
        const idx = activeIdx >= 0 ? activeIdx : 0;
        selectSuggestion(suggestItems[idx]);
      } else {
        run();
      }
    } else if (e.key === 'ArrowDown' && suggestItems.length) {
      e.preventDefault();
      moveActive(1);
    } else if (e.key === 'ArrowUp' && suggestItems.length) {
      e.preventDefault();
      moveActive(-1);
    } else if (e.key === 'Escape') {
      hideSuggest();
    }
  });
  $('suggestList').addEventListener('click', (e) => {
    const item = e.target.closest('.suggest-item');
    if (!item) return;
    selectSuggestion(suggestItems[parseInt(item.dataset.i, 10)]);
  });
  document.querySelectorAll('.tag').forEach((t) => {
    t.addEventListener('click', () => {
      $('codeInput').value = t.dataset.code;
      hideSuggest();
      run();
    });
  });
  // 点击输入区以外关闭下拉
  document.addEventListener('click', (e) => {
    if (!e.target.closest('.search-wrap')) hideSuggest();
  });
}

/* 键盘上下选择联想项 */
let activeIdx = -1;
function moveActive(step) {
  const items = document.querySelectorAll('#suggestList .suggest-item');
  if (!items.length) return;
  activeIdx = (activeIdx + step + items.length) % items.length;
  items.forEach((el, i) => el.classList.toggle('active', i === activeIdx));
  items[activeIdx].scrollIntoView({ block: 'nearest' });
}

document.addEventListener('DOMContentLoaded', () => {
  bindEvents();
});

/* ---------- Service Worker 注册（PWA 离线启动，仅 http/https 生效） ---------- */
if ('serviceWorker' in navigator && /^https?:$/.test(location.protocol)) {
  window.addEventListener('load', () => {
    navigator.serviceWorker.register('./sw.js').catch(() => {});
  });
}
