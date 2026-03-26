import { useState, useEffect, useRef } from 'react';
import CifUploader from './components/CifUploader';
import PseudopotentialSelector from './components/PseudopotentialSelector';
import CodeSelector from './components/CodeSelector';
import CalcTypeSelector from './components/CalcTypeSelector';
import PostProcessing from './components/PostProcessing';
import ParameterSettings from './components/ParameterSettings';
import GenerateButton from './components/GenerateButton';
import DependencyInfo from './components/DependencyInfo';
import StructureViewer from './components/StructureViewer';
import { parseCIF } from './utils/cifParser';
import { analyzeStructure } from './utils/structureAnalyzer';
import './App.css';

const DEFAULT_PARAMS = {
  ecutwfc: 60,
  ecutrho: 480,
  kx: 4,
  ky: 4,
  kz: 4,
  smearing: 'gaussian',
  degauss: 0.02,
  ppPath: './pseudo',
  ppLibrary: 'sssp_efficiency',
  hubbardU: '',
  functional: 'pbe',
  nprocs: 4,
  dftd3: false,
  spinPolarized: false,
  basisSet: 'def2-SVP',
  charge: 0,
  multiplicity: 1,
  mlipModel: '',
  mdTemp: 300,
  mdSteps: 5000,
  mdTimestep: 1.0,
};

function App() {
  const [cifText, setCifText] = useState('');
  const [selectedCode, setSelectedCode] = useState('qe');
  const [calcTypes, setCalcTypes] = useState(['scf']);
  const [postProcessing, setPostProcessing] = useState([]);
  const [params, setParams] = useState(DEFAULT_PARAMS);
  const [analysis, setAnalysis] = useState(null);
  const [parsedCIF, setParsedCIF] = useState(null);
  const prevCifRef = useRef('');

  // Analyze CIF when it changes
  useEffect(() => {
    if (!cifText.trim() || cifText === prevCifRef.current) return;
    prevCifRef.current = cifText;

    try {
      const parsed = parseCIF(cifText);
      if (parsed.atoms.length === 0) {
        setAnalysis(null);
        setParsedCIF(null);
        return;
      }

      setParsedCIF(parsed);
      const result = analyzeStructure(parsed);
      setAnalysis(result);

      // Auto-update params with recommendations
      setParams(prev => ({
        ...prev,
        ecutwfc: result.paramRecommendations.ecutwfc || prev.ecutwfc,
        ecutrho: result.paramRecommendations.ecutrho || prev.ecutrho,
        kx: result.paramRecommendations.kx || prev.kx,
        ky: result.paramRecommendations.ky || prev.ky,
        kz: result.paramRecommendations.kz || prev.kz,
        smearing: result.paramRecommendations.smearing || prev.smearing,
        degauss: result.paramRecommendations.degauss || prev.degauss,
        spinPolarized: result.paramRecommendations.spinPolarized ?? prev.spinPolarized,
        dftd3: result.paramRecommendations.dftd3 ?? prev.dftd3,
        mdTemp: result.paramRecommendations.mdTemp || prev.mdTemp,
        mdSteps: result.paramRecommendations.mdSteps || prev.mdSteps,
      }));
    } catch {
      setAnalysis(null);
    }
  }, [cifText]);

  return (
    <div className="app">
      <header className="app-header">
        <div className="header-content">
          <h1>DFT / MLIP Script Generator</h1>
          <p className="subtitle">
            Quantum ESPRESSO &middot; VASP &middot; ORCA &middot; UMA &middot; MACE
          </p>
        </div>
      </header>

      <main className="main-content">
        <CifUploader cifText={cifText} setCifText={setCifText} />

        <StructureViewer cifText={cifText} parsed={parsedCIF} />

        {analysis && (
          <div className="analysis-banner">
            <span className="analysis-icon">&#9881;</span>
            <div>
              <strong>{analysis.description}</strong>
              <span className="analysis-hint"> — Recommended settings auto-applied. Items marked with &#x25CF; are suggested.</span>
            </div>
          </div>
        )}

        <CodeSelector selectedCode={selectedCode} setSelectedCode={setSelectedCode} />
        <PseudopotentialSelector params={params} setParams={setParams} selectedCode={selectedCode} />
        <CalcTypeSelector calcTypes={calcTypes} setCalcTypes={setCalcTypes} analysis={analysis} />
        <PostProcessing selectedCode={selectedCode} postProcessing={postProcessing} setPostProcessing={setPostProcessing} analysis={analysis} />
        <DependencyInfo postProcessing={postProcessing} calcTypes={calcTypes} />
        <ParameterSettings params={params} setParams={setParams} selectedCode={selectedCode} />
        <GenerateButton cifText={cifText} selectedCode={selectedCode} calcTypes={calcTypes} postProcessing={postProcessing} params={params} />
      </main>

      <footer className="app-footer">
        <p>DFT/MLIP Script Generator &mdash; Auto-generate input scripts for computational materials science</p>
      </footer>
    </div>
  );
}

export default App;
