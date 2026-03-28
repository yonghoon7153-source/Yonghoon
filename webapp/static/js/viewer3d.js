/**
 * viewer3d.js - Three.js 3D Electrode Viewer for DEM Analysis
 * Uses ES module imports via importmap. Self-contained single file.
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

/* ── colour constants ──────────────────────────────────────── */
const COL = {
  AM_P: 0x222222, AM_S: 0x888888, SE: 0xf5e6a3,
  SE_TOP_REACH: 0x34d399, SE_NON_REACH: 0xf87171,
  SE_BOTTOM: 0xfbbf24, SE_TOP: 0x22d3ee,
  PATH: 0xffd700, BG: 0xf5f5f5,
};
const OPA = { SE: 0.85 };

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
    <button data-action="pathOnly">Path Only View</button>
    <button data-action="resetView">Reset</button>
    <button data-action="screenshot">Screenshot</button>`;
  container.appendChild(div);
  // Zoom slider (bottom-right)
  const zoomDiv = document.createElement('div');
  zoomDiv.className = 'viewer-zoom';
  zoomDiv.innerHTML = `
    <button id="zoom-out">−</button>
    <input type="range" id="zoom-slider" min="30" max="350" value="200" step="5">
    <button id="zoom-in">+</button>`;
  container.appendChild(zoomDiv);
  div._zoomDiv = zoomDiv;
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
.viewer-info{position:absolute;bottom:50px;left:12px;background:rgba(22,25,46,.9);
  border:1px solid #2a2d3e;border-radius:8px;padding:8px 12px;
  font:11px/1.5 'JetBrains Mono',monospace;color:#e4e6f0;z-index:10;max-width:240px;display:none}
.viewer-zoom{position:absolute;bottom:12px;right:12px;background:rgba(22,25,46,.9);
  border:1px solid #2a2d3e;border-radius:8px;padding:6px 10px;z-index:10;
  display:flex;align-items:center;gap:6px}
.viewer-zoom button{background:#555;color:#fff;border:none;border-radius:4px;width:24px;height:24px;
  cursor:pointer;font-size:16px;line-height:1;display:flex;align-items:center;justify-content:center}
.viewer-zoom button:hover{background:#777}
.viewer-zoom input[type=range]{width:100px;accent-color:#6c8cff}
.path-modal-overlay{position:fixed;inset:0;background:rgba(0,0,0,.6);z-index:1000;display:flex;align-items:center;justify-content:center}
.path-modal{background:#fff;border-radius:12px;padding:20px;max-width:90vw;max-height:90vh;position:relative}
.path-modal img{max-width:100%;max-height:75vh;border-radius:8px;border:1px solid #ddd}
.path-modal-info{margin-top:10px;font:12px/1.5 'JetBrains Mono',monospace;color:#333}
.path-modal-actions{display:flex;gap:8px;margin-top:12px;justify-content:flex-end}
.path-modal-actions button{background:#6c8cff;color:#fff;border:none;border-radius:6px;padding:6px 14px;cursor:pointer;font-size:12px}
.path-modal-actions button:hover{background:#8ba3ff}
.path-modal-close{position:absolute;top:8px;right:12px;background:none;border:none;font-size:20px;cursor:pointer;color:#888}`;
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
  controls.enableZoom = true;
  controls.zoomSpeed = 1.0;

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

  /* Path ON: cluster = blue, rest = yellow, all opacity 0.1 */
  const clusterSet = new Set(cluster.ids);
  const mesh = state.meshes.SE;
  const col = new THREE.Color();
  particles.forEach((p, i) => {
    if (clusterSet.has(p.id)) {
      col.setHex(0x2196F3);  // blue
    } else {
      col.setHex(COL.SE);    // yellow
    }
    mesh.setColorAt(i, col);
  });
  mesh.instanceColor.needsUpdate = true;
  mesh.material.opacity = 0.05;

  /* info text */
  const pi = state.currentPathIdx;
  let html = `<b>#${idx}</b> ${cluster.size}개`;
  html += cluster.percolating ? ` <span style="color:#34D399">✓</span>` : ` <span style="color:#F87171">✗</span>`;
  if (allPaths.length > 1) {
    const cat = allPaths[pi]?.category || '';
    const catLabel = cat === 'best' ? '🟢best' : cat === 'worst' ? '🔴worst' : '🟡mean';
    html += `<br>Path ${pi+1}/${allPaths.length} ${catLabel}`;
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
      const box = state.data.box;
      const halfX = (box.x_max - box.x_min) / 2;
      const halfY = (box.y_max - box.y_min) / 2;

      /* draw segments, detect periodic jumps */
      const mkSphere = (pos, color, size) => {
        const g = new THREE.SphereGeometry(size || 1.8, 12, 12);
        const m = new THREE.MeshPhongMaterial({ color });
        const s = new THREE.Mesh(g, m);
        s.position.copy(pos);
        return s;
      };

      for (let j = 0; j < pts.length - 1; j++) {
        const a = pts[j], b = pts[j+1];
        // Check periodic jump: x or z (=data y) distance > half box
        const dx = Math.abs(a.x - b.x);
        const dz = Math.abs(a.z - b.z);  // Three.js z = data y
        const isPeriodic = dx > halfX || dz > halfY;

        if (isPeriodic) {
          // Mark both ends with red spheres, skip the tube
          group.add(mkSphere(a, 0xFF0000, 0.5));
          group.add(mkSphere(b, 0xFF0000, 0.5));
        } else {
          const seg = new THREE.TubeGeometry(
            new THREE.LineCurve3(a, b), 1, 0.5, 6, false
          );
          const mat = new THREE.MeshPhongMaterial({
            color: COL.PATH, emissive: COL.PATH, emissiveIntensity: 0.3,
          });
          group.add(new THREE.Mesh(seg, mat));
        }
      }

      /* start(bottom cyan) / end(top red) markers */
      group.add(mkSphere(pts[0], 0x22D3EE, 1.8));
      group.add(mkSphere(pts[pts.length - 1], 0xF87171, 1.8));

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
  mesh.material.opacity = OPA.SE;  // back to 0.85
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
      } else if (action === 'pathOnly') {
        showPathOnlyView(renderer, scene, camera, state);
      }
    });
  });

  /* zoom slider */
  const zoomDiv = ctrlDiv._zoomDiv;
  if (zoomDiv) {
    const slider = zoomDiv.querySelector('#zoom-slider');
    const zoomIn = zoomDiv.querySelector('#zoom-in');
    const zoomOut = zoomDiv.querySelector('#zoom-out');

    function setZoom(dist) {
      dist = Math.max(30, Math.min(350, dist));
      const dir = camera.position.clone().sub(controls.target).normalize();
      camera.position.copy(controls.target).addScaledVector(dir, dist);
      controls.update();
      slider.value = dist;
    }

    slider.addEventListener('input', () => setZoom(parseInt(slider.value)));
    zoomIn.addEventListener('click', () => setZoom(parseInt(slider.value) - 20));
    zoomOut.addEventListener('click', () => setZoom(parseInt(slider.value) + 20));

    // Sync slider when wheel zoom changes camera distance
    controls.addEventListener('change', () => {
      const dist = camera.position.distanceTo(controls.target);
      slider.value = Math.max(30, Math.min(350, Math.round(dist)));
    });
  }
}

/* ── Path Only View (interactive 3D popup) ───────────────── */
function showPathOnlyView(renderer, scene, camera, state) {
  if (!state.pathGroup) {
    alert('먼저 Percolating Path를 선택하세요.');
    return;
  }

  const clusters = ((state.data.clusters || {}).clusters) || [];
  const cidx = state.currentClusterIdx || 0;
  const cluster = clusters[cidx];
  const pathIdx = state.currentPathIdx || 0;
  const allPaths = cluster ? (cluster.paths || (cluster.path ? [cluster.path] : [])) : [];
  const pathData = allPaths[pathIdx];

  if (!pathData || !pathData.ids) {
    alert('경로 데이터가 없습니다.');
    return;
  }

  // Build unwrapped path
  const box = state.data.box;
  const bx = box.x_max - box.x_min, by = box.y_max - box.y_min;
  const rawPts = pathData.ids.map(id => state.idIndex[id]).filter(Boolean);
  const unwrapped = [];
  let offX = 0, offY = 0;
  for (let i = 0; i < rawPts.length; i++) {
    const p = rawPts[i];
    let x = p.x + offX, y = p.y + offY;
    if (i > 0) {
      const prev = unwrapped[i-1];
      const dx = (p.x + offX) - prev.x;
      const dy = (p.y + offY) - prev.z;
      if (Math.abs(dx) > bx * 0.5) { offX -= Math.sign(dx) * bx; x = p.x + offX; }
      if (Math.abs(dy) > by * 0.5) { offY -= Math.sign(dy) * by; y = p.y + offY; }
    }
    unwrapped.push(new THREE.Vector3(x, p.z, y));
  }
  // Center path
  if (unwrapped.length > 0) {
    let mnX=Infinity,mxX=-Infinity,mnZ=Infinity,mxZ=-Infinity;
    unwrapped.forEach(v=>{mnX=Math.min(mnX,v.x);mxX=Math.max(mxX,v.x);mnZ=Math.min(mnZ,v.z);mxZ=Math.max(mxZ,v.z);});
    const sx=(box.x_min+box.x_max)/2-(mnX+mxX)/2, sz=(box.y_min+box.y_max)/2-(mnZ+mxZ)/2;
    unwrapped.forEach(v=>{v.x+=sx;v.z+=sz;});
  }

  // Create modal with canvas
  const overlay = document.createElement('div');
  overlay.className = 'path-modal-overlay';
  const cat = pathData.category || '';
  const catLabel = cat === 'best' ? 'Best' : cat === 'worst' ? 'Worst' : 'Mean';
  overlay.innerHTML = `
    <div class="path-modal" style="width:700px;max-width:90vw">
      <button class="path-modal-close" onclick="this.closest('.path-modal-overlay').remove()">&times;</button>
      <div style="font-size:14px;font-weight:bold;margin-bottom:8px;text-align:center">Li⁺ Ion Path (${catLabel})</div>
      <div id="path-viewer-container" style="width:100%;height:500px;border-radius:8px;overflow:hidden;background:#f5f5f5"></div>
      <div class="path-modal-info" style="text-align:center;margin-top:10px">
        Cluster #${cidx} | ${cluster.size} SE | τ = ${pathData.tortuosity} | Path: ${pathData.path_length} μm | Z: ${pathData.z_distance} μm
      </div>
      <div class="path-modal-actions">
        <button id="path-screenshot-btn">PNG 다운로드</button>
        <button onclick="this.closest('.path-modal-overlay').remove()">닫기</button>
      </div>
    </div>`;
  document.body.appendChild(overlay);
  overlay.onclick = (e) => { if (e.target === overlay) { cancelAnimationFrame(pathAnimId); overlay.remove(); }};

  // Create separate Three.js scene for path
  const container = document.getElementById('path-viewer-container');
  const r2 = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  r2.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  r2.setSize(container.clientWidth, container.clientHeight);
  r2.setClearColor(0xf5f5f5);
  container.appendChild(r2.domElement);

  const s2 = new THREE.Scene();
  const c2 = new THREE.PerspectiveCamera(50, container.clientWidth / container.clientHeight, 0.1, 10000);
  const ctrl2 = new OrbitControls(c2, r2.domElement);
  ctrl2.enableDamping = true;
  ctrl2.dampingFactor = 0.12;
  ctrl2.enableZoom = true;
  ctrl2.zoomSpeed = 1.0;

  s2.add(new THREE.AmbientLight(0xffffff, 0.5));
  const dl = new THREE.DirectionalLight(0xffffff, 0.8);
  dl.position.set(1, 1.5, 1);
  s2.add(dl);

  // Bounding box
  const bw = box.x_max-box.x_min, bh = box.z_max-box.z_min, bd = box.y_max-box.y_min;
  const cx = (box.x_min+box.x_max)/2, cy = (box.z_min+box.z_max)/2, cz = (box.y_min+box.y_max)/2;
  const bbEdges = new THREE.EdgesGeometry(new THREE.BoxGeometry(bw, bh, bd));
  const bbLine = new THREE.LineSegments(bbEdges, new THREE.LineBasicMaterial({color: 0x999999}));
  bbLine.position.set(cx, cy, cz);
  s2.add(bbLine);

  // Grid
  const grid = new THREE.GridHelper(Math.max(bw,bd)*1.2, 20, 0xcccccc, 0xe0e0e0);
  grid.position.set(cx, box.z_min, cz);
  s2.add(grid);

  // Path tubes
  for (let j = 0; j < unwrapped.length - 1; j++) {
    const seg = new THREE.TubeGeometry(new THREE.LineCurve3(unwrapped[j], unwrapped[j+1]), 1, 0.6, 6, false);
    s2.add(new THREE.Mesh(seg, new THREE.MeshPhongMaterial({color: COL.PATH, emissive: COL.PATH, emissiveIntensity: 0.5})));
  }
  // Start/end
  const mkS = (pos, color) => {
    const s = new THREE.Mesh(new THREE.SphereGeometry(2, 12, 12), new THREE.MeshPhongMaterial({color}));
    s.position.copy(pos); return s;
  };
  s2.add(mkS(unwrapped[0], 0x22D3EE));
  s2.add(mkS(unwrapped[unwrapped.length-1], 0xF87171));

  // Axis labels
  addAxisLabels(s2, box);

  // Camera
  const maxDim = Math.max(bw, bh, bd);
  c2.position.set(cx + maxDim*1.2, cy + maxDim*0.8, cz + maxDim*1.2);
  ctrl2.target.set(cx, cy, cz);
  ctrl2.update();

  // Zoom slider
  const zoomDiv = document.createElement('div');
  zoomDiv.style.cssText = 'position:absolute;bottom:10px;right:10px;display:flex;gap:4px;align-items:center;background:rgba(22,25,46,.8);padding:4px 8px;border-radius:6px;z-index:10';
  zoomDiv.innerHTML = '<button id="pv-zo" style="background:#555;color:#fff;border:none;border-radius:3px;width:20px;height:20px;cursor:pointer">−</button><input id="pv-zs" type="range" min="30" max="350" value="200" style="width:80px"><button id="pv-zi" style="background:#555;color:#fff;border:none;border-radius:3px;width:20px;height:20px;cursor:pointer">+</button>';
  container.style.position = 'relative';
  container.appendChild(zoomDiv);

  function setZoom2(d) {
    d = Math.max(30, Math.min(350, d));
    const dir = c2.position.clone().sub(ctrl2.target).normalize();
    c2.position.copy(ctrl2.target).addScaledVector(dir, d);
    ctrl2.update();
    document.getElementById('pv-zs').value = d;
  }
  document.getElementById('pv-zs').addEventListener('input', e => setZoom2(parseInt(e.target.value)));
  document.getElementById('pv-zi').addEventListener('click', () => setZoom2(parseInt(document.getElementById('pv-zs').value) - 20));
  document.getElementById('pv-zo').addEventListener('click', () => setZoom2(parseInt(document.getElementById('pv-zs').value) + 20));

  // Sync slider when wheel zoom changes camera distance
  ctrl2.addEventListener('change', () => {
    const dist = c2.position.distanceTo(ctrl2.target);
    document.getElementById('pv-zs').value = Math.max(30, Math.min(350, Math.round(dist)));
  });

  // Screenshot
  document.getElementById('path-screenshot-btn').addEventListener('click', () => {
    r2.render(s2, c2);
    const a = document.createElement('a');
    a.download = `li_ion_path_${catLabel.toLowerCase()}_tau${pathData.tortuosity}.png`;
    a.href = r2.domElement.toDataURL('image/png');
    a.click();
  });

  // Animate
  let pathAnimId;
  function animPath() {
    pathAnimId = requestAnimationFrame(animPath);
    ctrl2.update();
    r2.render(s2, c2);
  }
  animPath();

  // Cleanup on close
  overlay.querySelector('.path-modal-close').addEventListener('click', () => {
    cancelAnimationFrame(pathAnimId);
    r2.dispose();
    overlay.remove();
  });
}
