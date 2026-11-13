'use strict';
/*
 * Sayu Three-Body Problem
 * A graph view driven by a three-body gravity simulation.
 * The three most linked notes are suns; every other note is a planet.
 */
const obsidian = require('obsidian');

const VIEW_TYPE = 'sayu-three-body-problem';
const ICON_ID = 'sayu-3body';
const ICON_SVG =
  '<ellipse cx="50" cy="50" rx="43" ry="17" fill="none" stroke="currentColor" stroke-width="6" transform="rotate(-28 50 50)"/>' +
  '<circle cx="50" cy="50" r="10" fill="currentColor"/>' +
  '<circle cx="17" cy="62" r="7" fill="currentColor"/>' +
  '<circle cx="83" cy="38" r="7" fill="currentColor"/>';

const DEFAULTS = {
  mode: 'figure8',
  speed: 1,
  trail: 0.6,
  showLinks: false,
  showLabels: true,
  panelOpen: true,
  maxNodes: 3000,
  exclude: '',
};

/* ------------------------------------------------------------------ i18n */

const I18N = {
  en: {
    title: 'Three-Body Problem',
    open: 'Open Three-Body Problem graph',
    mode: 'Mode',
    figure8: 'Figure-8 (periodic)',
    triangle: 'Triangle (order → collapse)',
    chaos: 'Chaos',
    speed: 'Speed',
    trail: 'Trails',
    links: 'Show all links',
    labels: 'Show labels',
    pause: 'Pause',
    play: 'Play',
    reset: 'Reset',
    hint: 'Click: open note / Shift+click: follow / Double-click: reset view',
    notes: 'notes',
    linksCount: 'links',
    maxNodes: 'Maximum notes',
    maxNodesDesc: 'Notes with the most links are kept first.',
    exclude: 'Excluded folders',
    excludeDesc: 'One path prefix per line.',
  },
  ja: {
    title: '三体問題',
    open: '三体問題グラフを開く',
    mode: 'モード',
    figure8: '8の字（周期解）',
    triangle: '三角形（秩序→崩壊）',
    chaos: 'カオス',
    speed: '速度',
    trail: '軌跡',
    links: 'すべてのリンクを表示',
    labels: 'ラベルを表示',
    pause: '一時停止',
    play: '再生',
    reset: 'リセット',
    hint: 'クリック: ノートを開く ／ Shift+クリック: 追従 ／ ダブルクリック: 視点リセット',
    notes: 'ノート',
    linksCount: 'リンク',
    maxNodes: '最大ノート数',
    maxNodesDesc: 'リンクの多いノートから優先して表示します。',
    exclude: '除外フォルダ',
    excludeDesc: '1行に1つ、パスの先頭一致で指定します。',
  },
  ko: {
    title: '삼체 문제',
    open: '삼체 문제 그래프 열기',
    mode: '모드',
    figure8: '8자 궤도 (주기해)',
    triangle: '삼각형 (질서 → 붕괴)',
    chaos: '카오스',
    speed: '속도',
    trail: '궤적',
    links: '모든 링크 표시',
    labels: '라벨 표시',
    pause: '일시정지',
    play: '재생',
    reset: '리셋',
    hint: '클릭: 노트 열기 / Shift+클릭: 따라가기 / 더블클릭: 시점 초기화',
    notes: '노트',
    linksCount: '링크',
    maxNodes: '최대 노트 수',
    maxNodesDesc: '링크가 많은 노트부터 표시합니다.',
    exclude: '제외 폴더',
    excludeDesc: '한 줄에 하나씩, 경로 접두사로 지정합니다.',
  },
};

function currentLang() {
  let l = 'en';
  try {
    l = (window.localStorage.getItem('language') || 'en').slice(0, 2);
  } catch (e) {
    /* ignore */
  }
  return I18N[l] ? l : 'en';
}

function t(key) {
  return I18N[currentLang()][key] || I18N.en[key] || key;
}

/* --------------------------------------------------------------- physics */

const TAU = Math.PI * 2;
const P = {
  timeScale: 0.25, // simulation time units per real second at speed 1
  hMax: 0.0015, // largest integration sub-step
  maxSub: 12,
  epsP2: 0.04 * 0.04, // planet-sun softening (squared)
  wallR: 2.8, // soft outer wall radius
  wallK: 12,
  drag: 1.2, // "recapture" drag on stray planets
  rc2: 0.42 * 0.42, // drag is off inside this distance of a sun (squared)
};

/**
 * Bodies 0..2 are the suns (full three-body gravity on each other).
 * Bodies 3..n-1 are massless planets pulled by the three suns.
 * kind: 0 sun, 1 planet bound to a home sun, 2 outer-belt planet.
 */
class ThreeBodySim {
  constructor() {
    this.n = 0;
    this.M = new Float64Array(3);
    this.epsS2 = 0;
    this.kTrap = 0;
    this.X = this.Y = this.VX = this.VY = this.AX = this.AY = new Float64Array(0);
    this.kind = new Uint8Array(0);
    this.home = new Int8Array(0);
    this.depth = new Uint16Array(0);
  }

  /**
   * topo: { n, kind, home, depth }
   * weights: relative importance of the three suns (0..1)
   * prev: optional function(i) -> index in the previous state, or -1
   */
  build(topo, mode, weights, prev) {
    const old = { X: this.X, Y: this.Y, VX: this.VX, VY: this.VY };
    const n = topo.n;
    this.n = n;
    this.kind = topo.kind;
    this.home = topo.home;
    this.depth = topo.depth;
    this.X = new Float64Array(n);
    this.Y = new Float64Array(n);
    this.VX = new Float64Array(n);
    this.VY = new Float64Array(n);
    this.AX = new Float64Array(n);
    this.AY = new Float64Array(n);

    if (prev) {
      for (let i = 0; i < 3; i++) this.copyFrom(old, i, i);
    } else {
      this.initSuns(mode, weights);
    }
    for (let i = 3; i < n; i++) {
      const j = prev ? prev(i) : -1;
      if (j >= 0) this.copyFrom(old, i, j);
      else this.placePlanet(i);
    }
    this.computeAcc();
  }

  copyFrom(old, i, j) {
    this.X[i] = old.X[j];
    this.Y[i] = old.Y[j];
    this.VX[i] = old.VX[j];
    this.VY[i] = old.VY[j];
  }

  initSuns(mode, w) {
    const { X, Y, VX, VY, M } = this;
    if (mode === 'figure8') {
      // Chenciner–Montgomery figure-eight choreography (equal masses, G = 1).
      this.epsS2 = 0;
      this.kTrap = 0;
      const px = 0.97000436;
      const py = -0.24308753;
      const vx = -0.93240737;
      const vy = -0.86473146;
      M[0] = M[1] = M[2] = 1;
      X[0] = px; Y[0] = py;
      X[1] = -px; Y[1] = -py;
      X[2] = 0; Y[2] = 0;
      VX[0] = VX[1] = -vx / 2;
      VY[0] = VY[1] = -vy / 2;
      VX[2] = vx; VY[2] = vy;
      return;
    }
    this.epsS2 = 0.06 * 0.06;
    if (mode === 'triangle') {
      // Rotating equilateral triangle: looks stable, then falls apart.
      this.kTrap = 0.2;
      const s2 = 3 + this.epsS2;
      const v = Math.sqrt(3 / (s2 * Math.sqrt(s2)) + this.kTrap);
      for (let i = 0; i < 3; i++) {
        const th = Math.PI / 2 + (i * TAU) / 3;
        const k = v * (1 + (Math.random() - 0.5) * 2e-3);
        M[i] = 1;
        X[i] = Math.cos(th);
        Y[i] = Math.sin(th);
        VX[i] = -k * Math.sin(th);
        VY[i] = k * Math.cos(th);
      }
    } else {
      this.kTrap = 0.3;
      for (let i = 0; i < 3; i++) {
        const th = (i * TAU) / 3 + Math.random() * 1.2;
        const r = 0.55 + Math.random() * 0.55;
        const sp = 0.25 + Math.random() * 0.45;
        const ph = th + Math.PI / 2 + (Math.random() - 0.5) * 1.6;
        M[i] = 0.75 + 0.5 * (w ? w[i] : 0.5);
        X[i] = r * Math.cos(th);
        Y[i] = r * Math.sin(th);
        VX[i] = sp * Math.cos(ph);
        VY[i] = sp * Math.sin(ph);
      }
    }
    // Put the centre of mass at rest at the origin.
    let mt = 0, cx = 0, cy = 0, cvx = 0, cvy = 0;
    for (let i = 0; i < 3; i++) {
      mt += M[i];
      cx += M[i] * X[i]; cy += M[i] * Y[i];
      cvx += M[i] * VX[i]; cvy += M[i] * VY[i];
    }
    for (let i = 0; i < 3; i++) {
      X[i] -= cx / mt; Y[i] -= cy / mt;
      VX[i] -= cvx / mt; VY[i] -= cvy / mt;
    }
  }

  placePlanet(i) {
    const { X, Y, VX, VY, M } = this;
    const th = Math.random() * TAU;
    const c = Math.cos(th);
    const s = Math.sin(th);
    if (this.kind[i] === 1) {
      // Circular orbit around the home sun; directly linked notes sit closer.
      const h = this.home[i];
      const u = this.depth[i] <= 1 ? Math.random() * 0.5 : 0.3 + Math.random() * 0.7;
      const r = 0.1 + 0.25 * Math.sqrt(u);
      const v = Math.sqrt((M[h] * r * r) / Math.pow(r * r + P.epsP2, 1.5));
      const dir = h === 1 ? -1 : 1;
      X[i] = X[h] + r * c;
      Y[i] = Y[h] + r * s;
      VX[i] = VX[h] - dir * v * s;
      VY[i] = VY[h] + dir * v * c;
    } else {
      // Unlinked notes: a slow belt around the whole system.
      const r = 1.75 + 0.6 * Math.random();
      const v = Math.sqrt((M[0] + M[1] + M[2]) / r);
      X[i] = r * c;
      Y[i] = r * s;
      VX[i] = -v * s;
      VY[i] = v * c;
    }
  }

  computeAcc() {
    const { X, Y, VX, VY, AX, AY, M, n, kind } = this;
    const eS = this.epsS2;
    const k = this.kTrap;
    for (let i = 0; i < 3; i++) {
      let ax = -k * X[i];
      let ay = -k * Y[i];
      for (let j = 0; j < 3; j++) {
        if (j === i) continue;
        const dx = X[j] - X[i];
        const dy = Y[j] - Y[i];
        const r2 = dx * dx + dy * dy + eS;
        const inv = M[j] / (r2 * Math.sqrt(r2));
        ax += dx * inv;
        ay += dy * inv;
      }
      AX[i] = ax;
      AY[i] = ay;
    }
    const x0 = X[0], y0 = Y[0], x1 = X[1], y1 = Y[1], x2 = X[2], y2 = Y[2];
    const m0 = M[0], m1 = M[1], m2 = M[2];
    const eP = P.epsP2;
    const wallR = P.wallR;
    const wallR2 = wallR * wallR;
    for (let i = 3; i < n; i++) {
      const x = X[i];
      const y = Y[i];
      const dx0 = x0 - x, dy0 = y0 - y;
      const dx1 = x1 - x, dy1 = y1 - y;
      const dx2 = x2 - x, dy2 = y2 - y;
      const d0 = dx0 * dx0 + dy0 * dy0;
      const d1 = dx1 * dx1 + dy1 * dy1;
      const d2 = dx2 * dx2 + dy2 * dy2;
      let r2 = d0 + eP;
      let inv = m0 / (r2 * Math.sqrt(r2));
      let ax = dx0 * inv;
      let ay = dy0 * inv;
      r2 = d1 + eP;
      inv = m1 / (r2 * Math.sqrt(r2));
      ax += dx1 * inv;
      ay += dy1 * inv;
      r2 = d2 + eP;
      inv = m2 / (r2 * Math.sqrt(r2));
      ax += dx2 * inv;
      ay += dy2 * inv;

      const rr = x * x + y * y;
      if (rr > wallR2) {
        const r = Math.sqrt(rr);
        const f = (P.wallK * (r - wallR)) / r;
        ax -= x * f;
        ay -= y * f;
      }
      if (kind[i] === 1) {
        // Stray planets slowly match the nearest sun and get recaptured.
        let j = 0;
        let dm = d0;
        if (d1 < dm) { dm = d1; j = 1; }
        if (d2 < dm) { dm = d2; j = 2; }
        if (dm > P.rc2) {
          ax -= P.drag * (VX[i] - VX[j]);
          ay -= P.drag * (VY[i] - VY[j]);
        }
      }
      AX[i] = ax;
      AY[i] = ay;
    }
  }

  /** Leapfrog (kick-drift-kick). */
  advance(dt) {
    const { X, Y, VX, VY, AX, AY, n } = this;
    if (!n || dt <= 0) return;
    const nsub = Math.min(P.maxSub, Math.max(1, Math.ceil(dt / P.hMax)));
    const h = dt / nsub;
    const hh = h / 2;
    for (let s = 0; s < nsub; s++) {
      for (let i = 0; i < n; i++) {
        VX[i] += AX[i] * hh;
        VY[i] += AY[i] * hh;
        X[i] += VX[i] * h;
        Y[i] += VY[i] * h;
      }
      this.computeAcc();
      for (let i = 0; i < n; i++) {
        VX[i] += AX[i] * hh;
        VY[i] += AY[i] * hh;
      }
    }
  }
}

/* ----------------------------------------------------------------- scene */

const SUN_COLORS = ['#ffb14a', '#6cc8ff', '#ff6f91'];
const HEAD_COLORS = ['#ffd9a0', '#bfe6ff', '#ffc0cf', '#8f9bc4'];
const TRAIL_COLORS = [
  'rgba(255,177,74,0.5)',
  'rgba(108,200,255,0.5)',
  'rgba(255,111,145,0.5)',
  'rgba(143,155,196,0.32)',
];
const VIEW_HALF = 2.5; // world units visible from centre along the short axis
const TOP_LABELS = 6;

function el(doc, tag, cls, parent, text) {
  const e = doc.createElement(tag);
  if (cls) e.className = cls;
  if (text != null) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}

/**
 * Rendering + interaction. Independent of the Obsidian API so it can be
 * exercised in a plain browser page.
 * opts: { settings, onChange(), onOpen(path, evt), requestGraph() }
 */
class ThreeBodyScene {
  constructor(host, opts) {
    this.host = host;
    this.o = opts;
    this.s = opts.settings;
    this.doc = host.ownerDocument;
    this.win = this.doc.defaultView || window;
    this.sim = new ThreeBodySim();
    this.nodes = [];
    this.edges = new Int32Array(0);
    this.adj = [];
    this.deg = new Uint32Array(0);
    this.groups = [[], [], [], []];
    this.RAD = new Float32Array(0);
    this.SX = new Float32Array(0);
    this.SY = new Float32Array(0);
    this.PX = new Float32Array(0);
    this.PY = new Float32Array(0);
    this.sunPaths = null;
    this.pathIndex = new Map();

    this.w = 0;
    this.h = 0;
    this.dpr = 1;
    this.zoom = 100;
    this.zoom0 = 100;
    this.userZoom = false;
    this.panX = 0;
    this.panY = 0;
    this.follow = -1;
    this.hover = -1;
    this.mouseX = -1;
    this.mouseY = -1;
    this.paused = false;
    this.trailReset = true;
    this.fadeAcc = 0;
    this.drag = null;

    this.root = el(this.doc, 'div', 'sayu-3body-root', host);
    this.bg = el(this.doc, 'canvas', 'sayu-3body-layer', this.root);
    this.trail = el(this.doc, 'canvas', 'sayu-3body-layer', this.root);
    this.top = el(this.doc, 'canvas', 'sayu-3body-layer', this.root);
    this.bctx = this.bg.getContext('2d');
    this.tctx = this.trail.getContext('2d');
    this.ctx = this.top.getContext('2d');
    this.status = el(this.doc, 'div', 'sayu-3body-status', this.root);
    this.glow = SUN_COLORS.map((c) => this.makeGlow(c));
    this.font = this.win.getComputedStyle(host).fontFamily || 'sans-serif';

    this.buildPanel();
    this.bindEvents();

    this.ro = new this.win.ResizeObserver(() => this.resize());
    this.ro.observe(this.root);
    this.resize();

    this.running = true;
    this.last = 0;
    this.frame = this.frame.bind(this);
    this.raf = this.win.requestAnimationFrame(this.frame);
  }

  destroy() {
    this.running = false;
    this.win.cancelAnimationFrame(this.raf);
    this.ro.disconnect();
    this.root.remove();
  }

  /* ---------------------------------------------------------- graph data */

  /** graph: { nodes: [{path, name}], edges: [a0, b0, a1, b1, ...] } */
  setGraph(graph, fresh) {
    const src = graph.nodes;
    const e = graph.edges;
    const n0 = src.length;
    const deg0 = new Uint32Array(n0);
    for (let i = 0; i < e.length; i++) deg0[e[i]]++;
    const order = Array.from({ length: n0 }, (_, i) => i).sort((a, b) => deg0[b] - deg0[a]);

    // Keep the current suns on a refresh so the system does not jump.
    let suns = [];
    if (!fresh && this.sunPaths) {
      const byPath = new Map();
      src.forEach((nd, i) => byPath.set(nd.path, i));
      for (const p of this.sunPaths) {
        const i = p == null ? undefined : byPath.get(p);
        if (i !== undefined) suns.push(i);
      }
    }
    if (suns.length < 3) {
      suns = [];
      fresh = true;
    }
    for (const i of order) {
      if (suns.length >= 3) break;
      if (!suns.includes(i)) suns.push(i);
    }

    const max = Math.max(3, this.s.maxNodes | 0);
    const keep = suns.slice();
    for (const i of order) {
      if (keep.length >= max) break;
      if (i !== suns[0] && i !== suns[1] && i !== suns[2]) keep.push(i);
    }
    const remap = new Int32Array(n0).fill(-1);
    keep.forEach((o, i) => { remap[o] = i; });
    const nodes = keep.map((o) => ({ path: src[o].path, name: src[o].name }));
    while (nodes.length < 3) nodes.push({ path: null, name: '' });
    const n = nodes.length;

    const ed = [];
    const deg = new Uint32Array(n);
    const adj = Array.from({ length: n }, () => []);
    for (let i = 0; i + 1 < e.length; i += 2) {
      const a = remap[e[i]];
      const b = remap[e[i + 1]];
      if (a < 0 || b < 0) continue;
      ed.push(a, b);
      deg[a]++;
      deg[b]++;
      adj[a].push(b);
      adj[b].push(a);
    }

    // Each linked note belongs to the sun it is closest to in the link graph.
    const home = new Int8Array(n).fill(-1);
    const depth = new Uint16Array(n);
    const kind = new Uint8Array(n);
    const q = [0, 1, 2];
    home[0] = 0; home[1] = 1; home[2] = 2;
    for (let qi = 0; qi < q.length; qi++) {
      const a = q[qi];
      const na = adj[a];
      for (let k = 0; k < na.length; k++) {
        const b = na[k];
        if (home[b] < 0) {
          home[b] = home[a];
          depth[b] = depth[a] + 1;
          q.push(b);
        }
      }
    }
    for (let i = 3; i < n; i++) kind[i] = home[i] >= 0 ? 1 : 2;

    let prev = null;
    if (!fresh && this.sim.n) {
      const oldIdx = this.pathIndex;
      const oldKind = this.sim.kind;
      const oldHome = this.sim.home;
      prev = (i) => {
        const j = oldIdx.get(nodes[i].path);
        if (j === undefined || oldKind[j] !== kind[i] || oldHome[j] !== home[i]) return -1;
        return j;
      };
    }
    const dmax = Math.max(1, deg[0], deg[1], deg[2]);
    this.weights = [deg[0] / dmax, deg[1] / dmax, deg[2] / dmax];
    this.sim.build({ n, kind, home, depth }, this.s.mode, this.weights, prev);

    this.nodes = nodes;
    this.edges = Int32Array.from(ed);
    this.adj = adj;
    this.deg = deg;
    this.sunPaths = [nodes[0].path, nodes[1].path, nodes[2].path];
    this.pathIndex = new Map();
    nodes.forEach((nd, i) => { if (nd.path != null) this.pathIndex.set(nd.path, i); });
    this.groups = [[], [], [], []];
    this.RAD = new Float32Array(n);
    for (let i = 3; i < n; i++) {
      this.groups[kind[i] === 1 ? home[i] : 3].push(i);
      this.RAD[i] = Math.min(4.5, 0.9 + 0.45 * Math.sqrt(deg[i]));
    }
    this.SX = new Float32Array(n);
    this.SY = new Float32Array(n);
    this.PX = new Float32Array(n);
    this.PY = new Float32Array(n);
    if (this.follow >= 3 || this.follow >= n) this.follow = -1;
    this.hover = -1;
    this.trailReset = true;
    const real = nodes.filter((nd) => nd.path != null).length;
    this.status.textContent = `${real} ${t('notes')} · ${ed.length / 2} ${t('linksCount')}`;
  }

  /** Restart the simulation with the current graph (used on mode change). */
  restart() {
    const sim = this.sim;
    if (!sim.n) return;
    sim.build({ n: sim.n, kind: sim.kind, home: sim.home, depth: sim.depth }, this.s.mode, this.weights, null);
    this.trailReset = true;
  }

  /* ------------------------------------------------------------------ UI */

  buildPanel() {
    const d = this.doc;
    const s = this.s;
    const changed = () => this.o.onChange && this.o.onChange();
    const panel = el(d, 'div', 'sayu-3body-panel', this.root);
    panel.classList.toggle('is-collapsed', !s.panelOpen);
    const head = el(d, 'div', 'sayu-3body-head', panel);
    el(d, 'span', 'sayu-3body-title', head, t('title'));
    const chev = el(d, 'span', 'sayu-3body-chevron', head, s.panelOpen ? '▾' : '☰');
    head.addEventListener('click', () => {
      s.panelOpen = !s.panelOpen;
      panel.classList.toggle('is-collapsed', !s.panelOpen);
      chev.textContent = s.panelOpen ? '▾' : '☰';
      changed();
    });
    const body = el(d, 'div', 'sayu-3body-body', panel);

    let row = el(d, 'div', 'sayu-3body-row', body);
    el(d, 'label', '', row, t('mode'));
    const mode = el(d, 'select', 'dropdown', row);
    for (const m of ['figure8', 'triangle', 'chaos']) {
      const o = el(d, 'option', '', mode, t(m));
      o.value = m;
    }
    mode.value = s.mode;
    mode.addEventListener('change', () => {
      s.mode = mode.value;
      this.restart();
      changed();
    });

    const slider = (label, key, min, max, step) => {
      const r = el(d, 'div', 'sayu-3body-row', body);
      el(d, 'label', '', r, label);
      const inp = el(d, 'input', '', r);
      inp.type = 'range';
      inp.min = String(min);
      inp.max = String(max);
      inp.step = String(step);
      inp.value = String(s[key]);
      inp.addEventListener('input', () => {
        s[key] = parseFloat(inp.value);
        changed();
      });
    };
    slider(t('speed'), 'speed', 0.1, 4, 0.05);
    slider(t('trail'), 'trail', 0, 1, 0.02);

    const check = (label, key) => {
      const l = el(d, 'label', 'sayu-3body-check', body);
      const inp = el(d, 'input', '', l);
      inp.type = 'checkbox';
      inp.checked = !!s[key];
      el(d, 'span', '', l, label);
      inp.addEventListener('change', () => {
        s[key] = inp.checked;
        changed();
      });
    };
    check(t('links'), 'showLinks');
    check(t('labels'), 'showLabels');

    row = el(d, 'div', 'sayu-3body-buttons', body);
    const pause = el(d, 'button', '', row, t('pause'));
    pause.addEventListener('click', () => {
      this.paused = !this.paused;
      pause.textContent = this.paused ? t('play') : t('pause');
    });
    const reset = el(d, 'button', '', row, t('reset'));
    reset.addEventListener('click', () => {
      if (this.o.requestGraph) this.setGraph(this.o.requestGraph(), true);
      else this.restart();
      this.resetView();
    });
    el(d, 'div', 'sayu-3body-hint', body, t('hint'));
    // Keep canvas gestures from firing through the panel.
    for (const type of ['pointerdown', 'wheel', 'dblclick']) {
      panel.addEventListener(type, (ev) => ev.stopPropagation());
    }
  }

  bindEvents() {
    const c = this.root;
    const pos = (ev) => {
      const r = c.getBoundingClientRect();
      return [ev.clientX - r.left, ev.clientY - r.top];
    };
    c.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      const [x, y] = pos(ev);
      this.drag = { x, y, moved: false };
      c.setPointerCapture(ev.pointerId);
    });
    c.addEventListener('pointermove', (ev) => {
      const [x, y] = pos(ev);
      this.mouseX = x;
      this.mouseY = y;
      const d = this.drag;
      if (!d) return;
      if (!d.moved && Math.hypot(x - d.x, y - d.y) < 4) return;
      d.moved = true;
      this.panX -= (x - d.x) / this.zoom;
      this.panY -= (y - d.y) / this.zoom;
      d.x = x;
      d.y = y;
      this.trailReset = true;
    });
    const end = (ev) => {
      const d = this.drag;
      this.drag = null;
      if (!d || d.moved || ev.type !== 'pointerup') return;
      const [x, y] = pos(ev);
      const i = this.pick(x, y);
      if (i < 0) return;
      if (ev.shiftKey) {
        this.follow = this.follow === i ? -1 : i;
        this.panX = 0;
        this.panY = 0;
        this.trailReset = true;
      } else if (this.nodes[i].path != null && this.o.onOpen) {
        this.o.onOpen(this.nodes[i].path, ev);
      }
    };
    c.addEventListener('pointerup', end);
    c.addEventListener('pointercancel', end);
    c.addEventListener('pointerleave', () => {
      this.mouseX = this.mouseY = -1;
    });
    c.addEventListener(
      'wheel',
      (ev) => {
        ev.preventDefault();
        const [x, y] = pos(ev);
        const z0 = this.zoom;
        const z1 = Math.min(this.zoom0 * 40, Math.max(this.zoom0 * 0.3, z0 * Math.exp(-ev.deltaY * 0.0015)));
        this.panX += (x - this.w / 2) * (1 / z0 - 1 / z1);
        this.panY += (y - this.h / 2) * (1 / z0 - 1 / z1);
        this.zoom = z1;
        this.userZoom = true;
        this.trailReset = true;
      },
      { passive: false }
    );
    c.addEventListener('dblclick', (ev) => {
      const [x, y] = pos(ev);
      if (this.pick(x, y) < 0) this.resetView();
    });
  }

  resetView() {
    this.follow = -1;
    this.panX = 0;
    this.panY = 0;
    this.zoom = this.zoom0;
    this.userZoom = false;
    this.trailReset = true;
  }

  pick(mx, my) {
    const { SX, SY, RAD } = this;
    const n = this.sim.n;
    let best = -1;
    let bd = Infinity;
    for (let i = 0; i < n; i++) {
      const dx = SX[i] - mx;
      const dy = SY[i] - my;
      const d = dx * dx + dy * dy;
      const r = i < 3 ? 16 : Math.max(8, RAD[i] * this.rs + 4);
      // Suns win ties so they stay easy to grab inside their swarm.
      const score = i < 3 ? d * 0.25 : d;
      if (d < r * r && score < bd) {
        bd = score;
        best = i;
      }
    }
    return best;
  }

  makeGlow(color) {
    const c = this.doc.createElement('canvas');
    c.width = c.height = 128;
    const g = c.getContext('2d');
    const grad = g.createRadialGradient(64, 64, 0, 64, 64, 64);
    grad.addColorStop(0, 'rgba(255,255,255,1)');
    grad.addColorStop(0.1, color);
    grad.addColorStop(0.3, color + '55');
    grad.addColorStop(1, color + '00');
    g.fillStyle = grad;
    g.fillRect(0, 0, 128, 128);
    return c;
  }

  resize() {
    const w = this.root.clientWidth;
    const h = this.root.clientHeight;
    if (!w || !h) return;
    const dpr = this.win.devicePixelRatio || 1;
    this.w = w;
    this.h = h;
    this.dpr = dpr;
    for (const cv of [this.bg, this.trail, this.top]) {
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
      cv.getContext('2d').setTransform(dpr, 0, 0, dpr, 0, 0);
    }
    this.zoom0 = Math.min(w, h) / (2 * VIEW_HALF);
    if (!this.userZoom) this.zoom = this.zoom0;
    this.trailReset = true;
    // Static star field.
    const b = this.bctx;
    b.clearRect(0, 0, w, h);
    const stars = Math.round((w * h) / 2600);
    for (let i = 0; i < stars; i++) {
      const r = Math.random();
      b.fillStyle = `rgba(210,222,255,${0.08 + r * r * 0.5})`;
      const sz = r > 0.95 ? 1.6 : 1;
      b.fillRect(Math.random() * w, Math.random() * h, sz, sz);
    }
  }

  /* ---------------------------------------------------------------- loop */

  frame(ts) {
    if (!this.running) return;
    this.raf = this.win.requestAnimationFrame(this.frame);
    const dt = Math.min(0.034, Math.max(0, (ts - this.last) / 1000));
    this.last = ts;
    if (this.host.offsetParent === null || !this.w) return; // hidden tab
    const moving = !this.paused && this.sim.n > 0;
    if (moving) this.sim.advance(dt * P.timeScale * this.s.speed);
    this.draw(moving);
  }

  draw(moving) {
    const sim = this.sim;
    const n = sim.n;
    const { w, h, SX, SY, PX, PY, RAD } = this;
    const X = sim.X;
    const Y = sim.Y;
    const z = this.zoom;
    let cx = this.panX;
    let cy = this.panY;
    if (this.follow >= 0 && this.follow < n) {
      cx += X[this.follow];
      cy += Y[this.follow];
    }
    const ox = w / 2 - cx * z;
    const oy = h / 2 - cy * z;
    for (let i = 0; i < n; i++) {
      SX[i] = X[i] * z + ox;
      SY[i] = Y[i] * z + oy;
    }
    const rs = Math.min(2.6, Math.max(0.7, Math.sqrt(z / this.zoom0)));
    this.rs = rs;

    /* trails ---------------------------------------------------------- */
    const tc = this.tctx;
    const trail = this.s.trail;
    if (this.trailReset || trail <= 0) {
      tc.clearRect(0, 0, w, h);
      PX.set(SX);
      PY.set(SY);
      this.trailReset = false;
    } else if (moving) {
      // Fade in fixed 12% steps, less often for longer trails.
      this.fadeAcc += 1.5 * Math.pow(0.05, trail);
      if (this.fadeAcc >= 1) {
        tc.globalCompositeOperation = 'destination-out';
        tc.fillStyle = 'rgba(0,0,0,0.12)';
        while (this.fadeAcc >= 1) {
          tc.fillRect(0, 0, w, h);
          this.fadeAcc -= 1;
        }
        tc.globalCompositeOperation = 'source-over';
      }
      tc.lineCap = 'round';
      tc.lineWidth = 1;
      for (let g = 0; g < 4; g++) {
        const list = this.groups[g];
        tc.strokeStyle = TRAIL_COLORS[g];
        tc.beginPath();
        for (let k = 0; k < list.length; k++) {
          const i = list[k];
          const x = SX[i];
          const y = SY[i];
          if (x < -20 || y < -20 || x > w + 20 || y > h + 20) continue;
          const dx = x - PX[i];
          const dy = y - PY[i];
          if (dx * dx + dy * dy > 6400) continue;
          tc.moveTo(PX[i], PY[i]);
          tc.lineTo(x, y);
        }
        tc.stroke();
      }
      tc.lineWidth = 2.4 * rs;
      for (let i = 0; i < 3 && i < n; i++) {
        tc.strokeStyle = SUN_COLORS[i];
        tc.beginPath();
        tc.moveTo(PX[i], PY[i]);
        tc.lineTo(SX[i], SY[i]);
        tc.stroke();
      }
      PX.set(SX);
      PY.set(SY);
    }

    /* foreground ------------------------------------------------------ */
    const c = this.ctx;
    c.clearRect(0, 0, w, h);
    if (!n) return;

    if (this.s.showLinks) {
      const e = this.edges;
      c.strokeStyle = 'rgba(150,172,255,0.08)';
      c.lineWidth = 1;
      c.beginPath();
      for (let i = 0; i < e.length; i += 2) {
        const a = e[i];
        const b = e[i + 1];
        const ax = SX[a], ay = SY[a], bx = SX[b], by = SY[b];
        if ((ax < 0 && bx < 0) || (ay < 0 && by < 0) || (ax > w && bx > w) || (ay > h && by > h)) continue;
        c.moveTo(ax, ay);
        c.lineTo(bx, by);
      }
      c.stroke();
    }

    for (let g = 0; g < 4; g++) {
      const list = this.groups[g];
      c.fillStyle = HEAD_COLORS[g];
      c.beginPath();
      for (let k = 0; k < list.length; k++) {
        const i = list[k];
        const x = SX[i];
        const y = SY[i];
        if (x < -8 || y < -8 || x > w + 8 || y > h + 8) continue;
        const r = RAD[i] * rs;
        if (r < 1.7) {
          c.rect(x - r, y - r, r * 2, r * 2);
        } else {
          c.moveTo(x + r, y);
          c.arc(x, y, r, 0, TAU);
        }
      }
      c.fill();
    }

    c.globalCompositeOperation = 'lighter';
    for (let i = 0; i < 3; i++) {
      const R = (22 + 10 * sim.M[i]) * rs;
      c.drawImage(this.glow[i], SX[i] - R, SY[i] - R, R * 2, R * 2);
    }
    c.globalCompositeOperation = 'source-over';

    // Hover follows the pointer even while the bodies move underneath it.
    const hov = this.mouseX >= 0 && !(this.drag && this.drag.moved) ? this.pick(this.mouseX, this.mouseY) : -1;
    this.hover = hov;
    this.root.style.cursor = hov >= 0 ? 'pointer' : this.drag && this.drag.moved ? 'grabbing' : 'default';

    if (hov >= 0) {
      const nb = this.adj[hov];
      c.strokeStyle = 'rgba(255,255,255,0.45)';
      c.lineWidth = 1;
      c.beginPath();
      for (let k = 0; k < nb.length; k++) {
        c.moveTo(SX[hov], SY[hov]);
        c.lineTo(SX[nb[k]], SY[nb[k]]);
      }
      c.stroke();
      c.fillStyle = '#ffffff';
      c.beginPath();
      for (let k = 0; k < nb.length; k++) {
        const j = nb[k];
        if (j < 3) continue;
        const r = RAD[j] * rs + 0.8;
        c.moveTo(SX[j] + r, SY[j]);
        c.arc(SX[j], SY[j], r, 0, TAU);
      }
      c.fill();
    }
    const ring = (i, color) => {
      const r = (i < 3 ? 12 : RAD[i] + 4) * rs;
      c.strokeStyle = color;
      c.lineWidth = 1.5;
      c.beginPath();
      c.arc(SX[i], SY[i], r, 0, TAU);
      c.stroke();
    };
    if (this.follow >= 0 && this.follow < n) ring(this.follow, 'rgba(255,255,255,0.35)');
    if (hov >= 0) ring(hov, 'rgba(255,255,255,0.9)');

    /* labels ---------------------------------------------------------- */
    c.textBaseline = 'middle';
    c.lineJoin = 'round';
    const label = (i, size, alpha) => {
      const name = this.nodes[i].name;
      if (!name) return;
      const x = SX[i] + (i < 3 ? 14 : RAD[i] + 6) * rs;
      const y = SY[i];
      if (x < -200 || y < -10 || x > w || y > h + 10) return;
      c.font = `${size}px ${this.font}`;
      c.lineWidth = 3;
      c.strokeStyle = `rgba(4,5,11,${alpha * 0.85})`;
      c.strokeText(name, x, y);
      c.fillStyle = `rgba(232,238,255,${alpha})`;
      c.fillText(name, x, y);
    };
    if (this.s.showLabels) {
      const top = Math.min(n, 3 + TOP_LABELS);
      for (let i = 3; i < top; i++) if (this.deg[i] > 0 && i !== hov) label(i, 11, 0.55);
      for (let i = 0; i < 3; i++) if (i !== hov) label(i, 13, 0.92);
    }
    if (hov >= 0) {
      const nb = this.adj[hov];
      for (let k = 0; k < nb.length && k < 16; k++) label(nb[k], 11, 0.8);
      label(hov, 13, 1);
    }
  }
}

/* -------------------------------------------------------------- obsidian */

function buildGraphData(app, settings) {
  const ex = String(settings.exclude || '')
    .split('\n')
    .map((s) => s.trim())
    .filter(Boolean);
  const files = app.vault.getMarkdownFiles().filter((f) => !ex.some((p) => f.path.startsWith(p)));
  const index = new Map();
  files.forEach((f, i) => index.set(f.path, i));
  const count = files.length;
  const seen = new Set();
  const edges = [];
  const resolved = app.metadataCache.resolvedLinks || {};
  for (const src in resolved) {
    const a = index.get(src);
    if (a === undefined) continue;
    for (const dst in resolved[src]) {
      const b = index.get(dst);
      if (b === undefined || b === a) continue;
      const key = a < b ? a * count + b : b * count + a;
      if (seen.has(key)) continue;
      seen.add(key);
      edges.push(a, b);
    }
  }
  return { nodes: files.map((f) => ({ path: f.path, name: f.basename })), edges };
}

class ThreeBodyView extends obsidian.ItemView {
  constructor(leaf, plugin) {
    super(leaf);
    this.plugin = plugin;
    this.scene = null;
    this.sig = '';
  }

  getViewType() {
    return VIEW_TYPE;
  }

  getDisplayText() {
    return t('title');
  }

  getIcon() {
    return ICON_ID;
  }

  async onOpen() {
    const plugin = this.plugin;
    this.contentEl.empty();
    this.contentEl.addClass('sayu-3body-view');
    this.scene = new ThreeBodyScene(this.contentEl, {
      settings: plugin.settings,
      onChange: () => plugin.saveSoon(),
      requestGraph: () => buildGraphData(this.app, plugin.settings),
      onOpen: (path, ev) => {
        const file = this.app.vault.getAbstractFileByPath(path);
        if (!(file instanceof obsidian.TFile)) return;
        const split = ev.ctrlKey || ev.metaKey;
        this.app.workspace.getLeaf(split ? 'split' : 'tab').openFile(file);
      },
    });
    this.app.workspace.onLayoutReady(() => this.reload(true));
    this.registerEvent(
      this.app.metadataCache.on('resolved', obsidian.debounce(() => this.reload(false), 2000, true))
    );
  }

  reload(fresh) {
    if (!this.scene) return;
    const g = buildGraphData(this.app, this.plugin.settings);
    const sig = g.nodes.length + ':' + g.edges.length;
    if (!fresh && sig === this.sig) return;
    this.sig = sig;
    this.scene.setGraph(g, fresh);
  }

  async onClose() {
    if (this.scene) this.scene.destroy();
    this.scene = null;
  }
}

class ThreeBodySettingTab extends obsidian.PluginSettingTab {
  constructor(app, plugin) {
    super(app, plugin);
    this.plugin = plugin;
  }

  display() {
    const { containerEl } = this;
    const plugin = this.plugin;
    containerEl.empty();
    new obsidian.Setting(containerEl)
      .setName(t('maxNodes'))
      .setDesc(t('maxNodesDesc'))
      .addText((tx) =>
        tx.setValue(String(plugin.settings.maxNodes)).onChange((v) => {
          const num = parseInt(v, 10);
          if (!(num >= 3)) return;
          plugin.settings.maxNodes = Math.min(20000, num);
          plugin.saveSoon();
          plugin.reloadViews();
        })
      );
    new obsidian.Setting(containerEl)
      .setName(t('exclude'))
      .setDesc(t('excludeDesc'))
      .addTextArea((ta) =>
        ta.setValue(plugin.settings.exclude).onChange((v) => {
          plugin.settings.exclude = v;
          plugin.saveSoon();
          plugin.reloadViews();
        })
      );
  }
}

class ThreeBodyPlugin extends obsidian.Plugin {
  async onload() {
    this.settings = Object.assign({}, DEFAULTS, await this.loadData());
    this.saveSoon = obsidian.debounce(() => this.saveData(this.settings), 600, true);
    this.reloadViews = obsidian.debounce(
      () => {
        for (const leaf of this.app.workspace.getLeavesOfType(VIEW_TYPE)) {
          if (leaf.view instanceof ThreeBodyView) leaf.view.reload(true);
        }
      },
      800,
      true
    );
    obsidian.addIcon(ICON_ID, ICON_SVG);
    this.registerView(VIEW_TYPE, (leaf) => new ThreeBodyView(leaf, this));
    this.addRibbonIcon(ICON_ID, t('open'), () => this.activateView());
    this.addCommand({ id: 'open', name: t('open'), callback: () => this.activateView() });
    this.addSettingTab(new ThreeBodySettingTab(this.app, this));
  }

  async activateView() {
    const { workspace } = this.app;
    let leaf = workspace.getLeavesOfType(VIEW_TYPE)[0];
    if (!leaf) {
      leaf = workspace.getLeaf('tab');
      await leaf.setViewState({ type: VIEW_TYPE, active: true });
    }
    workspace.revealLeaf(leaf);
  }
}

module.exports = ThreeBodyPlugin;
module.exports.internals = { ThreeBodySim, ThreeBodyScene, P };
