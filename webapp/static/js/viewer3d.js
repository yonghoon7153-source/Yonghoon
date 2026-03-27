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
const OPA = { SE: 0.85, SE_CLUSTER: 0.5 };

/* ── control-panel HTML ────────────────────────────────────── */
function buildControls(container) {
  const div = document.createElement('div');
  div.className = 'viewer-controls';
  div.innerHTML = `
    <label><input type="checkbox" data-layer="AM_P" checked> AM_P</label>
    <label><input type="checkbox" data-layer="AM_S" checked> AM_S</label>
    <label><input type="checkbox" data-layer="SE" checked> SE</label>
    <hr>
    <label><input type="checkbox" id="path-toggle"> <span style="font-size:11px">Percolating Path</span></label>
    <div id="path-controls" style="display:none">
      <div style="display:flex;gap:4px;align-items:center;margin-top:3px">
        <button id="path-prev" style="background:#555;color:#fff;border:none;border-radius:3px;padding:1px 6px;cursor:pointer;font-size:12px">&lt;</button>
        <span id="path-current" style="font-size:11px;color:#e4e6f0;min-width:30px;text-align:center">-</span>
        <button id="path-next" style="background:#555;color:#fff;border:none;border-radius:3px;padding:1px 6px;cursor:pointer;font-size:12px">&gt;</button>
        <span id="path-total" style="font-size:10px;color:#7c8194">/ -</span>
      </div>
      <div id="cluster-info" style="font-size:10px;color:#e4e6f0;margin-top:3px;line-height:1.5"></div>
    </div>
    <hr>
    <button data-action="resetView">Reset</button>
    <button data-action="screenshot">Screenshot</button>`;
  container.appendChild(div);
  // Separate info panel
  const infoDiv = document.createElement('div');
  infoDiv.className = 'viewer-info';
  infoDiv.id = 'viewer-info';
  container.appendChild(infoDiv);
  div._infoEl = infoDiv;
  return div;
}

/* ── inject CSS (once) ─────────────────────────────────────── */
function injectCSS() {
  if (document.getElementById('viewer3d-css')) return;
  const s = document.createElement('style');
  s.id = 'viewer3d-css';
  s.textContent = `
.viewer-container canvas{display:block}
.viewer-controls{position:absolute;top:10px;right:10px;background:rgba(22,25,46,.9);
  border:1px solid #2a2d3e;border-radius:8px;padding:8px 12px;display:inline-flex;flex-direction:column;gap:3px;
  font:12px/1.4 'Inter',sans-serif;color:#e4e6f0;z-index:10;user-select:none;width:140px}
.viewer-controls label{display:flex;align-items:center;gap:5px;cursor:pointer;font-size:11px}
.viewer-controls hr{border:none;border-top:1px solid #2a2d3e;margin:3px 0}
.viewer-controls button{background:#555;color:#fff;border:none;border-radius:4px;padding:3px 8px;
  cursor:pointer;font-size:10px;margin-top:1px}
.viewer-controls button:hover{background:#777}
.viewer-info{position:absolute;bottom:12px;left:12px;background:rgba(22,25,46,.9);
  border:1px solid #2a2d3e;border-radius:8px;padding:8px 12px;
  font:11px/1.5 'JetBrains Mono',monospace;color:#e4e6f0;z-index:10;max-width:240px;display:none}`;
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
  controls.zoomSpeed = 0.5;
  controls.minDistance = 30;
  controls.maxDistance = 400;

  /* prevent page scroll when zooming inside viewer */
  renderer.domElement.addEventListener('wheel', e => e.preventDefault(), {passive: false});

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
  state.infoEl = ctrlDiv._infoEl || document.getElementById('viewer-info');

  /* ── data fetch & build ──────────────────────────────────── */
  fetch(dataUrl).then(r => r.json()).then(data => {
    state.data = data;
    buildScene(scene, camera, controls, data, state);
    wireControls(ctrlDiv, renderer, camera, controls, scene, state);

    // Setup percolating path navigation
    const clusters = (data.clusters || {}).clusters || [];
    const percClusters = clusters.filter(c => c.percolating && c.path);
    state.percClusters = percClusters;
    state.percIdx = 0;

    const totalEl = ctrlDiv.querySelector('#path-total');
    const currentEl = ctrlDiv.querySelector('#path-current');
    const infoEl = ctrlDiv.querySelector('#cluster-info');
    const prevBtn = ctrlDiv.querySelector('#path-prev');
    const nextBtn = ctrlDiv.querySelector('#path-next');

    if (totalEl) totalEl.textContent = `/ ${percClusters.length}`;

    const pathControls = ctrlDiv.querySelector('#path-controls');
    const pathToggle = ctrlDiv.querySelector('#path-toggle');

    state.currentPathIdx = 0;

    function showPercCluster(clusterIdx, pathIdx) {
      if (!percClusters.length) { if (infoEl) infoEl.innerHTML = 'No percolating'; return; }
      clusterIdx = ((clusterIdx % percClusters.length) + percClusters.length) % percClusters.length;
      state.percIdx = clusterIdx;
      const origIdx = clusters.indexOf(percClusters[clusterIdx]);
      const cluster = percClusters[clusterIdx];
      const allPaths = cluster.paths || (cluster.path ? [cluster.path] : []);
      const pi = pathIdx !== undefined ? pathIdx : 0;
      highlightCluster(origIdx, scene, state, infoEl, pi);
      if (currentEl) currentEl.textContent = `${clusterIdx+1}-${pi+1}`;
      if (totalEl) totalEl.textContent = `/ ${percClusters.length} (${allPaths.length}경로)`;
    }

    function clearPath() {
      if (state.pathGroup) { scene.remove(state.pathGroup); state.pathGroup = null; }
      resetSEColors(state);
      if (infoEl) infoEl.innerHTML = '';
      if (currentEl) currentEl.textContent = '-';
    }

    if (pathToggle) pathToggle.addEventListener('change', () => {
      if (pathToggle.checked) {
        pathControls.style.display = 'block';
        state.currentPathIdx = 0;
        if (percClusters.length > 0) showPercCluster(state.percIdx, 0);
      } else {
        pathControls.style.display = 'none';
        clearPath();
      }
    });

    if (prevBtn) prevBtn.addEventListener('click', () => {
      const pi = (state.currentPathIdx || 0) - 1;
      if (pi < 0) {
        // 이전 클러스터의 마지막 path
        const prevCluster = percClusters[((state.percIdx - 1) + percClusters.length) % percClusters.length];
        const prevPaths = prevCluster.paths || [prevCluster.path];
        showPercCluster(state.percIdx - 1, prevPaths.length - 1);
      } else {
        showPercCluster(state.percIdx, pi);
      }
    });
    if (nextBtn) nextBtn.addEventListener('click', () => {
      const curCluster = percClusters[state.percIdx];
      const curPaths = curCluster.paths || [curCluster.path];
      const pi = (state.currentPathIdx || 0) + 1;
      if (pi >= curPaths.length) {
        // 다음 클러스터의 첫 path
        showPercCluster(state.percIdx + 1, 0);
      } else {
        showPercCluster(state.percIdx, pi);
      }
    });

    animate();
  }).catch(err => {
    console.error('viewer3d: failed to load data', err);
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

  /* cluster input is handled via event listener in data fetch callback */
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

  /* set zoom limits based on data size */
  controls.minDistance = maxDim * 0.3;
  controls.maxDistance = maxDim * 4;

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
    canvas.width = 512; canvas.height = 128;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#000000';
    ctx.font = 'bold 72px Arial, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(l.text, 256, 64);
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

/* ── Cluster highlight by index ────────────────────────────── */
function highlightCluster(idx, scene, state, infoEl, pathIdx) {
  const particles = state.seParticles;
  if (!particles) return;

  /* clear previous */
  if (state.pathGroup) { scene.remove(state.pathGroup); state.pathGroup = null; }
  resetSEColors(state);

  const clusterList = ((state.data.clusters || {}).clusters) || [];
  if (idx < 0 || idx >= clusterList.length) {
    if (infoEl) infoEl.innerHTML = '';
    return;
  }

  const cluster = clusterList[idx];
  const allPaths = cluster.paths || (cluster.path ? [cluster.path] : []);
  state.currentClusterPaths = allPaths;
  state.currentPathIdx = pathIdx || 0;

  /* highlight cluster SE particles */
  const clusterSet = new Set(cluster.ids);
  const mesh = state.meshes.SE;
  const col = new THREE.Color();
  particles.forEach((p, i) => {
    if (clusterSet.has(p.id)) {
      col.setHex(0x1565C0);
    } else {
      col.setHex(COL.SE);
    }
    mesh.setColorAt(i, col);
  });
  mesh.instanceColor.needsUpdate = true;
  mesh.material.transparent = true;
  mesh.material.opacity = OPA.SE_CLUSTER;
  mesh.material.depthWrite = false;

  /* info text */
  const pi = state.currentPathIdx;
  let html = `<b>#${idx}</b> ${cluster.size}개`;
  html += cluster.percolating ? ` <span style="color:#34D399">✓</span>` : ` <span style="color:#F87171">✗</span>`;
  if (allPaths.length > 1) {
    html += `<br>Path ${pi+1}/${allPaths.length}`;
  }

  /* draw selected path */
  const path = allPaths[pi];
  if (path && path.ids) {
    const pts = path.ids.map(id => {
      const p = state.idIndex[id];
      return p ? new THREE.Vector3(p.x, p.z, p.y) : null;
    }).filter(Boolean);

    if (pts.length >= 2) {
      const group = new THREE.Group();

      /* straight tube segments between SE centers */
      for (let j = 0; j < pts.length - 1; j++) {
        const seg = new THREE.TubeGeometry(
          new THREE.LineCurve3(pts[j], pts[j+1]), 1, 0.5, 6, false
        );
        const mat = new THREE.MeshPhongMaterial({
          color: COL.PATH, emissive: COL.PATH, emissiveIntensity: 0.3,
        });
        group.add(new THREE.Mesh(seg, mat));
      }

      /* start(bottom)/end(top) markers */
      const mkSphere = (pos, color) => {
        const g = new THREE.SphereGeometry(1.8, 12, 12);
        const m = new THREE.MeshPhongMaterial({ color });
        const s = new THREE.Mesh(g, m);
        s.position.copy(pos);
        return s;
      };
      group.add(mkSphere(pts[0], 0x22D3EE));
      group.add(mkSphere(pts[pts.length - 1], 0xF87171));

      scene.add(group);
      state.pathGroup = group;

      html += `<br>τ=${path.tortuosity} L=${path.path_length}μm`;
    }
  }

  if (infoEl) infoEl.innerHTML = html;
}

function resetSEColors(state) {
  const mesh = state.meshes.SE;
  if (!mesh) return;
  const col = new THREE.Color(COL.SE);
  const particles = state.seParticles || [];
  particles.forEach((p, i) => mesh.setColorAt(i, col));
  mesh.instanceColor.needsUpdate = true;
  mesh.material.transparent = false;
  mesh.material.opacity = OPA.SE;
  mesh.material.depthWrite = true;
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
