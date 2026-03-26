import { useState } from 'react';
import JSZip from 'jszip';
import { saveAs } from 'file-saver';
import { parseCIF } from '../utils/cifParser';
import { generateQEScripts } from '../generators/qeGenerator';
import { generateORCAScripts } from '../generators/orcaGenerator';
import { generateUMAScripts, generateMACEScripts } from '../generators/mlipGenerator';
import { DEPENDENCIES, CALC_TYPES } from '../config/calculations';

export default function GenerateButton({
  cifText, selectedCode, calcTypes, postProcessing, params
}) {
  const [preview, setPreview] = useState(null);
  const [previewFile, setPreviewFile] = useState(null);

  const generate = () => {
    if (!cifText.trim()) {
      alert('Please provide a CIF file first.');
      return;
    }

    const parsed = parseCIF(cifText);
    if (parsed.atoms.length === 0) {
      alert('Could not parse any atoms from the CIF. Please check the format.');
      return;
    }

    // Resolve dependencies
    const resolvedCalcTypes = new Set(calcTypes);
    for (const pp of postProcessing) {
      const deps = DEPENDENCIES[pp] || [];
      for (const dep of deps) {
        if (CALC_TYPES.find(ct => ct.id === dep)) {
          resolvedCalcTypes.add(dep);
        }
      }
    }

    let files;
    switch (selectedCode) {
      case 'qe':
        files = generateQEScripts(parsed, params, [...resolvedCalcTypes], postProcessing);
        break;
      case 'orca':
        files = generateORCAScripts(parsed, params, [...resolvedCalcTypes], postProcessing);
        break;
      case 'uma':
        files = generateUMAScripts(parsed, params, [...resolvedCalcTypes], postProcessing);
        break;
      case 'mace':
        files = generateMACEScripts(parsed, params, [...resolvedCalcTypes], postProcessing);
        break;
      default:
        alert('Please select a calculation code.');
        return;
    }

    setPreview(files);
    setPreviewFile(Object.keys(files)[0]);
  };

  const download = () => {
    if (!preview) return;

    const zip = new JSZip();
    const folder = zip.folder(`${selectedCode}_workflow`);

    for (const [name, content] of Object.entries(preview)) {
      folder.file(name, content);
    }

    // Also include the original CIF
    folder.file('structure.cif', cifText);

    zip.generateAsync({ type: 'blob' }).then((blob) => {
      saveAs(blob, `${selectedCode}_workflow.zip`);
    }).catch((err) => {
      alert('Error generating ZIP: ' + err.message);
    });
  };

  return (
    <section className="card">
      <h2>7. Generate Scripts</h2>
      <div className="generate-actions">
        <button className="btn btn-primary btn-large" onClick={generate}>
          Generate Scripts
        </button>
        {preview && (
          <button className="btn btn-success btn-large" onClick={download}>
            Download ZIP
          </button>
        )}
      </div>

      {preview && (
        <div className="preview-section">
          <h3>Generated Files Preview</h3>
          <div className="preview-layout">
            <div className="file-list">
              {Object.keys(preview).map((name) => (
                <button
                  key={name}
                  className={`file-item ${previewFile === name ? 'active' : ''}`}
                  onClick={() => setPreviewFile(name)}
                >
                  {name}
                </button>
              ))}
            </div>
            <div className="file-preview">
              <div className="file-preview-header">
                <span>{previewFile}</span>
                <button
                  className="btn btn-sm"
                  onClick={() => {
                    navigator.clipboard.writeText(preview[previewFile]);
                  }}
                >
                  Copy
                </button>
              </div>
              <pre className="file-content">
                <code>{preview[previewFile]}</code>
              </pre>
            </div>
          </div>
        </div>
      )}
    </section>
  );
}
