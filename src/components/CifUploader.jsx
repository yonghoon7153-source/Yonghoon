import { useRef } from 'react';

export default function CifUploader({ cifText, setCifText }) {
  const fileInput = useRef(null);

  const handleFile = (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => setCifText(ev.target.result);
    reader.readAsText(file);
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
          accept=".cif"
          onChange={handleFile}
          style={{ display: 'none' }}
        />
        <span className="or-divider">or paste below</span>
      </div>
      <textarea
        className="cif-textarea"
        rows={12}
        placeholder="Paste your CIF file content here..."
        value={cifText}
        onChange={(e) => setCifText(e.target.value)}
      />
      {cifText && (
        <div className="cif-preview">
          {cifText.length > 0 && `${cifText.split('\n').length} lines loaded`}
        </div>
      )}
    </section>
  );
}
