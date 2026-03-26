import { useEffect, useRef, useState } from 'react';

// Element colors (CPK coloring)
const ELEMENT_COLORS = {
  H: 0xffffff, He: 0xd9ffff, Li: 0xcc80ff, Be: 0xc2ff00, B: 0xffb5b5,
  C: 0x909090, N: 0x3050f8, O: 0xff0d0d, F: 0x90e050, Ne: 0xb3e3f5,
  Na: 0xab5cf2, Mg: 0x8aff00, Al: 0xbfa6a6, Si: 0xf0c8a0, P: 0xff8000,
  S: 0xffff30, Cl: 0x1ff01f, Ar: 0x80d1e3, K: 0x8f40d4, Ca: 0x3dff00,
  Sc: 0xe6e6e6, Ti: 0xbfc2c7, V: 0xa6a6ab, Cr: 0x8a99c7, Mn: 0x9c7ac7,
  Fe: 0xe06633, Co: 0xf090a0, Ni: 0x50d050, Cu: 0xc88033, Zn: 0x7d80b0,
  Ga: 0xc28f8f, Ge: 0x668f8f, As: 0xbd80e3, Se: 0xffa100, Br: 0xa62929,
  Rb: 0x702eb0, Sr: 0x00ff00, Y: 0x94ffff, Zr: 0x94e0e0, Nb: 0x73c2c9,
  Mo: 0x54b5b5, Ru: 0x248f8f, Rh: 0x0a7d8c, Pd: 0x006985, Ag: 0xc0c0c0,
  Cd: 0xffd98f, In: 0xa67573, Sn: 0x668080, Sb: 0x9e63b5, Te: 0xd47a00,
  I: 0x940094, Cs: 0x57178f, Ba: 0x00c900, La: 0x70d4ff, Ce: 0xffffc7,
  Pt: 0xd0d0e0, Au: 0xffd123, Pb: 0x575961, Bi: 0x9e4fb5,
};

const ELEMENT_RADII = {
  H: 0.31, He: 0.28, Li: 1.28, Be: 0.96, B: 0.84, C: 0.76, N: 0.71, O: 0.66,
  F: 0.57, Na: 1.66, Mg: 1.41, Al: 1.21, Si: 1.11, P: 1.07, S: 1.05, Cl: 1.02,
  K: 2.03, Ca: 1.76, Ti: 1.60, V: 1.53, Cr: 1.39, Mn: 1.39, Fe: 1.32, Co: 1.26,
  Ni: 1.24, Cu: 1.32, Zn: 1.22, Ga: 1.22, Ge: 1.20, As: 1.19, Se: 1.20, Br: 1.20,
  Zr: 1.75, Nb: 1.64, Mo: 1.54, Ag: 1.45, Cd: 1.44, In: 1.42, Sn: 1.39, Sb: 1.39,
  Te: 1.38, I: 1.39, Cs: 2.44, Ba: 2.15, La: 2.07, Pt: 1.36, Au: 1.36, Pb: 1.46,
};

export default function StructureViewer({ cifText, parsed }) {
  const canvasRef = useRef(null);
  const [viewMode, setViewMode] = useState('ball-stick');
  const [showBonds, setShowBonds] = useState(true);
  const [rotation, setRotation] = useState({ x: -20, y: 30 });
  const [zoom, setZoom] = useState(1);
  const isDragging = useRef(false);
  const lastMouse = useRef({ x: 0, y: 0 });

  useEffect(() => {
    if (!parsed || !parsed.atoms.length || !canvasRef.current) return;
    drawStructure();
  }, [parsed, viewMode, showBonds, rotation, zoom]);

  const drawStructure = () => {
    const canvas = canvasRef.current;
    const ctx = canvas.getContext('2d');
    const w = canvas.width;
    const h = canvas.height;

    ctx.clearRect(0, 0, w, h);

    // Dark background
    ctx.fillStyle = '#1a1a2e';
    ctx.fillRect(0, 0, w, h);

    if (!parsed || !parsed.atoms.length) return;

    const { v1, v2, v3 } = parsed.cellParams;

    // Convert fractional to Cartesian
    const toCart = (fx, fy, fz) => [
      fx * v1[0] + fy * v2[0] + fz * v3[0],
      fx * v1[1] + fy * v2[1] + fz * v3[1],
      fx * v1[2] + fy * v2[2] + fz * v3[2],
    ];

    // Rotation matrix
    const rx = rotation.x * Math.PI / 180;
    const ry = rotation.y * Math.PI / 180;
    const cosX = Math.cos(rx), sinX = Math.sin(rx);
    const cosY = Math.cos(ry), sinY = Math.sin(ry);

    const rotate = ([x, y, z]) => {
      // Rotate Y
      let x1 = x * cosY + z * sinY;
      let z1 = -x * sinY + z * cosY;
      // Rotate X
      let y1 = y * cosX - z1 * sinX;
      let z2 = y * sinX + z1 * cosX;
      return [x1, y1, z2];
    };

    // Center of cell
    const center = toCart(0.5, 0.5, 0.5);

    // Project atoms
    const scale = zoom * Math.min(w, h) / (Math.max(parsed.a, parsed.b, parsed.c) * 2.5);

    const project = (cart) => {
      const centered = [cart[0] - center[0], cart[1] - center[1], cart[2] - center[2]];
      const rotated = rotate(centered);
      return {
        x: w / 2 + rotated[0] * scale,
        y: h / 2 - rotated[1] * scale,
        z: rotated[2],
      };
    };

    // Draw unit cell edges
    const corners = [
      [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
      [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1],
    ].map(([a, b, c]) => toCart(a, b, c));

    const edges = [
      [0, 1], [0, 2], [0, 3], [1, 4], [1, 5], [2, 4],
      [2, 6], [3, 5], [3, 6], [4, 7], [5, 7], [6, 7],
    ];

    ctx.strokeStyle = 'rgba(100, 180, 255, 0.4)';
    ctx.lineWidth = 1;
    for (const [i, j] of edges) {
      const p1 = project(corners[i]);
      const p2 = project(corners[j]);
      ctx.beginPath();
      ctx.moveTo(p1.x, p1.y);
      ctx.lineTo(p2.x, p2.y);
      ctx.stroke();
    }

    // Draw axis labels
    const axisLabels = [
      { pos: toCart(1.1, 0, 0), label: 'a', color: '#ff6b6b' },
      { pos: toCart(0, 1.1, 0), label: 'b', color: '#51cf66' },
      { pos: toCart(0, 0, 1.1), label: 'c', color: '#339af0' },
    ];
    ctx.font = 'bold 14px sans-serif';
    for (const ax of axisLabels) {
      const p = project(ax.pos);
      ctx.fillStyle = ax.color;
      ctx.fillText(ax.label, p.x, p.y);
    }

    // Prepare atoms with projected coords
    const atoms3d = parsed.atoms.map((atom) => {
      const cart = toCart(atom.x, atom.y, atom.z);
      const proj = project(cart);
      return { ...atom, cart, proj, radius: (ELEMENT_RADII[atom.symbol] || 1.0) };
    });

    // Sort by z for painter's algorithm
    atoms3d.sort((a, b) => a.proj.z - b.proj.z);

    // Draw bonds
    if (showBonds) {
      ctx.lineWidth = 2;
      for (let i = 0; i < atoms3d.length; i++) {
        for (let j = i + 1; j < atoms3d.length; j++) {
          const a = atoms3d[i], b = atoms3d[j];
          const dx = a.cart[0] - b.cart[0];
          const dy = a.cart[1] - b.cart[1];
          const dz = a.cart[2] - b.cart[2];
          const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
          const maxBond = (a.radius + b.radius) * 1.3;
          if (dist < maxBond && dist > 0.5) {
            ctx.strokeStyle = 'rgba(200, 200, 200, 0.3)';
            ctx.beginPath();
            ctx.moveTo(a.proj.x, a.proj.y);
            ctx.lineTo(b.proj.x, b.proj.y);
            ctx.stroke();
          }
        }
      }
    }

    // Draw atoms
    for (const atom of atoms3d) {
      const color = ELEMENT_COLORS[atom.symbol] || 0xcccccc;
      const r = viewMode === 'ball-stick'
        ? atom.radius * scale * 0.25
        : atom.radius * scale * 0.5;
      const clampedR = Math.max(3, Math.min(r, 25));

      const hex = '#' + color.toString(16).padStart(6, '0');

      // Sphere shading
      const gradient = ctx.createRadialGradient(
        atom.proj.x - clampedR * 0.3, atom.proj.y - clampedR * 0.3, clampedR * 0.1,
        atom.proj.x, atom.proj.y, clampedR
      );
      gradient.addColorStop(0, lightenColor(hex, 60));
      gradient.addColorStop(0.7, hex);
      gradient.addColorStop(1, darkenColor(hex, 40));

      ctx.beginPath();
      ctx.arc(atom.proj.x, atom.proj.y, clampedR, 0, Math.PI * 2);
      ctx.fillStyle = gradient;
      ctx.fill();

      // Label
      if (zoom > 0.8 && clampedR > 8) {
        ctx.fillStyle = '#fff';
        ctx.font = `${Math.max(9, clampedR * 0.7)}px sans-serif`;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(atom.symbol, atom.proj.x, atom.proj.y);
      }
    }

    // Legend
    const elements = [...new Set(parsed.atoms.map(a => a.symbol))];
    const legendX = 12;
    let legendY = h - elements.length * 22 - 10;
    ctx.font = '12px sans-serif';
    for (const el of elements) {
      const color = ELEMENT_COLORS[el] || 0xcccccc;
      const hex = '#' + color.toString(16).padStart(6, '0');
      ctx.fillStyle = hex;
      ctx.beginPath();
      ctx.arc(legendX + 7, legendY + 7, 6, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillStyle = '#ccc';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      ctx.fillText(el, legendX + 18, legendY + 7);
      legendY += 22;
    }
  };

  // Mouse interaction
  const handleMouseDown = (e) => {
    isDragging.current = true;
    lastMouse.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseMove = (e) => {
    if (!isDragging.current) return;
    const dx = e.clientX - lastMouse.current.x;
    const dy = e.clientY - lastMouse.current.y;
    setRotation(prev => ({
      x: prev.x - dy * 0.5,
      y: prev.y + dx * 0.5,
    }));
    lastMouse.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseUp = () => { isDragging.current = false; };

  const handleWheel = (e) => {
    e.preventDefault();
    setZoom(prev => Math.max(0.3, Math.min(3, prev - e.deltaY * 0.001)));
  };

  if (!parsed || !parsed.atoms.length) return null;

  return (
    <section className="card viewer-card">
      <h2>Structure Viewer</h2>
      <div className="viewer-controls">
        <select value={viewMode} onChange={e => setViewMode(e.target.value)}>
          <option value="ball-stick">Ball &amp; Stick</option>
          <option value="spacefill">Space Fill</option>
        </select>
        <label className="viewer-toggle">
          <input type="checkbox" checked={showBonds} onChange={e => setShowBonds(e.target.checked)} />
          Bonds
        </label>
        <button className="btn btn-sm" onClick={() => { setRotation({ x: -20, y: 30 }); setZoom(1); }}>
          Reset View
        </button>
        <span className="viewer-hint">Drag to rotate, scroll to zoom</span>
      </div>
      <canvas
        ref={canvasRef}
        width={800}
        height={500}
        className="viewer-canvas"
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onWheel={handleWheel}
      />
      <div className="viewer-info">
        {parsed.nat} atoms | {parsed.ntyp} species ({parsed.elements.join(', ')}) |
        a={parsed.a.toFixed(3)} b={parsed.b.toFixed(3)} c={parsed.c.toFixed(3)} Å |
        α={parsed.alpha}° β={parsed.beta}° γ={parsed.gamma}°
      </div>
    </section>
  );
}

function lightenColor(hex, percent) {
  const num = parseInt(hex.slice(1), 16);
  const r = Math.min(255, ((num >> 16) & 0xff) + percent);
  const g = Math.min(255, ((num >> 8) & 0xff) + percent);
  const b = Math.min(255, (num & 0xff) + percent);
  return `rgb(${r},${g},${b})`;
}

function darkenColor(hex, percent) {
  const num = parseInt(hex.slice(1), 16);
  const r = Math.max(0, ((num >> 16) & 0xff) - percent);
  const g = Math.max(0, ((num >> 8) & 0xff) - percent);
  const b = Math.max(0, (num & 0xff) - percent);
  return `rgb(${r},${g},${b})`;
}
