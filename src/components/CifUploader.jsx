import { useRef, useState } from 'react';

export default function CifUploader({ cifText, setCifText }) {
  const fileInput = useRef(null);
  const [dragging, setDragging] = useState(false);

  const readFile = (file) => {
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => setCifText(ev.target.result);
    reader.readAsText(file);
  };

  const handleFile = (e) => {
    readFile(e.target.files[0]);
  };

  const handleDragOver = (e) => {
    e.preventDefault();
    e.stopPropagation();
    setDragging(true);
  };

  const handleDragLeave = (e) => {
    e.preventDefault();
    e.stopPropagation();
    setDragging(false);
  };

  const handleDrop = (e) => {
    e.preventDefault();
    e.stopPropagation();
    setDragging(false);
    const file = e.dataTransfer.files[0];
    if (file) readFile(file);
  };

  return (
    <section className="card">
      <h2>1. Structure Input (CIF)</h2>
      <div className="cif-upload-area">
        <button className="btn btn-secondary" onClick={() => fileInput.current.click()}>
          Upload CIF File
        </button>
        <input
          ref={fileInput}
          type="file"
          accept=".cif,.txt"
          onChange={handleFile}
          style={{ display: 'none' }}
        />
        <span className="or-divider">or drag &amp; drop / paste below</span>
      </div>
      <textarea
        className={`cif-textarea${dragging ? ' dragging' : ''}`}
        placeholder="Paste your CIF file content here, or drag & drop a file..."
        value={cifText}
        onChange={(e) => setCifText(e.target.value)}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        style={{ height: cifText ? `${Math.min(Math.max(cifText.split('\n').length * 1.5, 12), 40)}em` : undefined }}
      />
      {cifText && (
        <div className="cif-preview">
          {cifText.split('\n').length} lines loaded
        </div>
      )}
    </section>
  );
}
