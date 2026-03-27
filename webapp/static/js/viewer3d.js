/**
 * viewer3d.js - Three.js 3D Electrode Viewer for DEM Analysis
 * Uses ES module imports via importmap. Self-contained single file.
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

/* ── colour constants ──────────────────────────────────────── */
const COL = {
  AM_P: 0x222222, AM_S: 0x888888, SE: 0xf5d77a,
  SE_TOP_REACH: 0x34d399, SE_NON_REACH: 0xf87171,
  SE_BOTTOM: 0xfbbf24, SE_TOP: 0x22d3ee,
  PATH: 0xffd700, BG: 0xf5f5f5,
};
const OPA = { SE: 0.3, SE_REACH: 0.6, SE_NON: 0.15, SE_BOUND: 0.6 };

/* ── control-panel HTML ────────────────────────────────────── */
function buildControls(container) {
  const div = document.createElement('div');
  div.className = 'viewer-controls';
  div.innerHTML = `
    <label><input type="checkbox" data-layer="AM_P" checked> AM_P</label>
    <label><input type="checkbox" data-layer="AM_S" checked> AM_S</label>
    <label><input type="checkbox" data-layer="SE" checked> SE</label>
    <label><input type="checkbox" data-layer="percolation"> SE Percolation</label>
    <label><input type="checkbox" data-layer="forceChains"> Force Chains</label>
    <hr>
    <button data-action="resetView">Reset View</button>
    <button data-action="screenshot">Screenshot</button>
    <div class="viewer-info" id="viewer-info"></div>`;
  container.appendChild(div);
  return div;
}

/* ── inject CSS (once) ─────────────────────────────────────── */
function injectCSS() {
  if (document.getElementById('viewer3d-css')) return;
  const s = document.createElement('style');
  s.id = 'viewer3d-css';
  s.textContent = `
.viewer-container{position:relative;width:100%;height:600px;border-radius:10px;overflow:hidden;background:#0f1117}
.viewer-controls{position:absolute;top:10px;right:10px;background:rgba(22,25,46,.92);
  border:1px solid #2a2d3e;border-radius:8px;padding:10px 14px;display:flex;flex-direction:column;gap:4px;
  font:13px/1.5 'Inter',sans-serif;color:#e4e6f0;z-index:10;user-select:none}
.viewer-controls label{display:flex;align-items:center;gap:6px;cursor:pointer}
.viewer-controls hr{border:none;border-top:1px solid #2a2d3e;margin:4px 0}
.viewer-controls button{background:#6c8cff;color:#fff;border:none;border-radius:5px;padding:5px 10px;
  cursor:pointer;font-size:12px;margin-top:2px}
.viewer-controls button:hover{background:#8ba3ff}
.viewer-info{margin-top:6px;font-size:11px;color:#7c8194;max-width:200px;word-wrap:break-word}`;
  document.head.appendChild(s);
}

/* ── main init function ────────────────────────────────────── */
export function initElectrodeViewer(containerId, dataUrl) {
  injectCSS();
  const container = document.getElementById(containerId);
  if (!container) { console.error('viewer3d: container not found:', containerId); return; }
  container.classList.add('viewer-container');

  /* renderer */
  const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setClearColor(COL.BG);
  container.appendChild(renderer.domElement);

  /* scene, camera */
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(50, container.clientWidth / container.clientHeight, 0.1, 10000);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.12;
  controls.zoomSpeed = 0.5;  // 줌 속도 절반으로

  /* lights */
  scene.add(new THREE.AmbientLight(0xffffff, 0.4));
  const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
  dirLight.position.set(1, 1.5, 1);
  scene.add(dirLight);

  /* state */
  const state = {
    data: null, meshes: {}, percolationOn: false, pathGroup: null,
    selectedComponent: null, infoEl: null,
  };

  /* controls panel */
  const ctrlDiv = buildControls(container);
  state.infoEl = ctrlDiv.querySelector('#viewer-info');

  /* ── data fetch & build ──────────────────────────────────── */
  fetch(dataUrl).then(r => r.json()).then(data => {
    state.data = data;
    buildScene(scene, camera, controls, data, state);
    wireControls(ctrlDiv, renderer, camera, controls, scene, state);
    animate();
  }).catch(err => {
    console.error('viewer3d: failed to load data', err);
    state.infoEl.textContent = 'Error loading data';
  });

  /* ── animation loop ──────────────────────────────────────── */
  function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
  }

  /* resize */
  const ro = new ResizeObserver(() => {
    const w = container.clientWidth, h = container.clientHeight;
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    renderer.setSize(w, h);
  });
  ro.observe(container);

  /* raycaster for click */
  const raycaster = new THREE.Raycaster();
  const pointer = new THREE.Vector2();
  renderer.domElement.addEventListener('click', e => {
    if (!state.data || !state.meshes.SE) return;
    const rect = renderer.domElement.getBoundingClientRect();
    pointer.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    pointer.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(pointer, camera);
    const hit = raycaster.intersectObject(state.meshes.SE);
    if (hit.length > 0) onSEClick(hit[0].instanceId, scene, state);
  });
}

/* ── build scene from data ─────────────────────────────────── */
function buildScene(scene, camera, controls, data, state) {
  const box = data.box;
  // Z-up coordinate system: Three.js Y-up → swap Y↔Z for display
  // Data: x, y (horizontal), z (up=electrode height)
  // Three.js: x, y(=data z, up), z(=data y)
  const cx = (box.x_min + box.x_max) / 2;
  const cy = (box.z_min + box.z_max) / 2;  // data Z → Three.js Y (up)
  const cz = (box.y_min + box.y_max) / 2;  // data Y → Three.js Z
  const bw = box.x_max - box.x_min;
  const bh = box.z_max - box.z_min;  // height = data Z range
  const bd = box.y_max - box.y_min;
  const maxDim = Math.max(bw, bh, bd);

  /* camera position - isometric-ish view */
  camera.position.set(cx + maxDim * 1.2, cy + maxDim * 0.8, cz + maxDim * 1.2);
  controls.target.set(cx, cy, cz);
  controls.update();
  state.defaultCamPos = camera.position.clone();
  state.defaultTarget = controls.target.clone();

  /* bounding box wireframe */
  const bbGeo = new THREE.BoxGeometry(bw, bh, bd);
  const bbMat = new THREE.LineBasicMaterial({ color: 0x888888 });
  const bbEdges = new THREE.EdgesGeometry(bbGeo);
  const bbLine = new THREE.LineSegments(bbEdges, bbMat);
  bbLine.position.set(cx, cy, cz);
  scene.add(bbLine);

  /* grid at bottom (Y=0 in Three.js = Z=0 in data) */
  const gridSize = Math.max(bw, bd) * 1.2;
  const grid = new THREE.GridHelper(gridSize, 20, 0xcccccc, 0xe0e0e0);
  grid.position.set(cx, box.z_min, cz);
  scene.add(grid);

  /* axis labels (Z-up convention) */
  addAxisLabels(scene, box);

  /* group particles by type */
  const groups = { AM_P: [], AM_S: [], SE: [] };
  const idIndex = {};
  data.particles.forEach((p, i) => {
    if (groups[p.type]) groups[p.type].push(p);
    idIndex[p.id] = p;
  });
  state.idIndex = idIndex;

  /* instanced meshes */
  state.meshes.AM_P = createInstancedSpheres(groups.AM_P, 16, COL.AM_P, 1.0, false);
  state.meshes.AM_S = createInstancedSpheres(groups.AM_S, 16, COL.AM_S, 1.0, false);
  state.meshes.SE = createInstancedSpheres(groups.SE, 12, COL.SE, OPA.SE, true);
  state.seParticles = groups.SE;

  Object.values(state.meshes).forEach(m => { if (m) scene.add(m); });
}

/* ── instanced sphere builder ──────────────────────────────── */
function createInstancedSpheres(particles, segments, color, opacity, transparent) {
  if (!particles.length) return null;
  const geo = new THREE.SphereGeometry(1, segments, segments);
  const mat = new THREE.MeshPhongMaterial({
    color, transparent, opacity, depthWrite: !transparent, side: THREE.FrontSide,
  });
  const mesh = new THREE.InstancedMesh(geo, mat, particles.length);
  const dummy = new THREE.Object3D();
  const col = new THREE.Color();
  particles.forEach((p, i) => {
    dummy.position.set(p.x, p.z, p.y);  // Z-up: swap Y↔Z
    dummy.scale.setScalar(p.r);
    dummy.updateMatrix();
    mesh.setMatrixAt(i, dummy.matrix);
    mesh.setColorAt(i, col.setHex(color));
  });
  mesh.instanceMatrix.needsUpdate = true;
  if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
  mesh.userData.particles = particles;
  return mesh;
}

/* ── axis labels using sprite text ─────────────────────────── */
function addAxisLabels(scene, box) {
  // Z-up: X→right, Y→depth(Three.js Z), Z→up(Three.js Y)
  const labels = [
    { text: 'X (μm)', pos: [box.x_max + 5, box.z_min, (box.y_min+box.y_max)/2] },
    { text: 'Y (μm)', pos: [(box.x_min+box.x_max)/2, box.z_min, box.y_max + 5] },
    { text: 'Z (μm)', pos: [box.x_min - 5, box.z_max + 5, (box.y_min+box.y_max)/2] },
  ];
  labels.forEach(l => {
    const canvas = document.createElement('canvas');
    canvas.width = 160; canvas.height = 48;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#333333';
    ctx.font = 'bold 22px Inter, sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(l.text, 80, 32);
    const tex = new THREE.CanvasTexture(canvas);
    const mat = new THREE.SpriteMaterial({ map: tex, depthWrite: false });
    const sprite = new THREE.Sprite(mat);
    sprite.position.set(...l.pos);
    sprite.scale.set(20, 8, 1);
    scene.add(sprite);
  });
}

/* ── SE percolation colouring ──────────────────────────────── */
function applyPercolation(state, on) {
  state.percolationOn = on;
  const mesh = state.meshes.SE;
  if (!mesh) return;
  const particles = mesh.userData.particles;
  const perc = state.data.percolation || {};
  const topSet = new Set(perc.top_reachable || []);
  const botSet = new Set(perc.bottom_se || []);
  const topSESet = new Set(perc.top_se || []);
  const col = new THREE.Color();

  particles.forEach((p, i) => {
    if (on) {
      if (botSet.has(p.id)) { col.setHex(COL.SE_BOTTOM); }
      else if (topSESet.has(p.id)) { col.setHex(COL.SE_TOP); }
      else if (topSet.has(p.id)) { col.setHex(COL.SE_TOP_REACH); }
      else { col.setHex(COL.SE_NON_REACH); }
    } else {
      col.setHex(COL.SE);
    }
    mesh.setColorAt(i, col);
  });
  mesh.instanceColor.needsUpdate = true;

  /* adjust per-instance opacity via material - we use a single material so set average */
  if (on) {
    mesh.material.opacity = OPA.SE_REACH;
  } else {
    mesh.material.opacity = OPA.SE;
  }
}

/* ── SE click handler ──────────────────────────────────────── */
function onSEClick(instanceId, scene, state) {
  const particles = state.seParticles;
  if (!particles || instanceId >= particles.length) return;
  const clicked = particles[instanceId];
  const info = state.infoEl;

  /* find path containing this particle */
  const paths = state.data.paths || [];
  const path = paths.find(p => p.ids && p.ids.includes(clicked.id));

  /* clear old path */
  if (state.pathGroup) { scene.remove(state.pathGroup); state.pathGroup = null; }

  if (!path) {
    info.textContent = `SE #${clicked.id} - no tortuosity path`;
    return;
  }

  /* draw tube along path */
  const pts = path.ids.map(id => {
    const p = state.idIndex[id];
    return p ? new THREE.Vector3(p.x, p.z, p.y) : null;  // Z-up swap
  }).filter(Boolean);

  if (pts.length < 2) return;

  const curve = new THREE.CatmullRomCurve3(pts);
  const tubeGeo = new THREE.TubeGeometry(curve, pts.length * 4, 0.8, 8, false);
  const tubeMat = new THREE.MeshPhongMaterial({ color: COL.PATH, emissive: COL.PATH, emissiveIntensity: 0.4 });
  const tubeMesh = new THREE.Mesh(tubeGeo, tubeMat);

  const group = new THREE.Group();
  group.add(tubeMesh);
  scene.add(group);
  state.pathGroup = group;

  info.innerHTML = `<b>Tortuosity path</b><br>
    tau = ${path.tortuosity.toFixed(2)}<br>
    path = ${path.path_length.toFixed(1)} um<br>
    z dist = ${path.z_distance.toFixed(1)} um`;
}

/* ── wire up control panel ─────────────────────────────────── */
function wireControls(ctrlDiv, renderer, camera, controls, scene, state) {
  ctrlDiv.querySelectorAll('input[type=checkbox]').forEach(cb => {
    cb.addEventListener('change', () => {
      const layer = cb.dataset.layer;
      if (layer === 'percolation') {
        applyPercolation(state, cb.checked);
      } else if (layer === 'forceChains') {
        /* placeholder for force chains */
      } else if (state.meshes[layer]) {
        state.meshes[layer].visible = cb.checked;
      }
    });
  });

  ctrlDiv.querySelectorAll('button').forEach(btn => {
    const action = btn.dataset.action;
    btn.addEventListener('click', () => {
      if (action === 'resetView') {
        camera.position.copy(state.defaultCamPos);
        controls.target.copy(state.defaultTarget);
        controls.update();
      } else if (action === 'screenshot') {
        renderer.render(scene, camera);
        const link = document.createElement('a');
        link.download = 'electrode_3d.png';
        link.href = renderer.domElement.toDataURL('image/png');
        link.click();
      }
    });
  });
}
