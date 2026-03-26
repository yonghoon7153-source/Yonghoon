import { useState } from 'react';
import CifUploader from './components/CifUploader';
import PseudopotentialSelector from './components/PseudopotentialSelector';
import CodeSelector from './components/CodeSelector';
import CalcTypeSelector from './components/CalcTypeSelector';
import PostProcessing from './components/PostProcessing';
import ParameterSettings from './components/ParameterSettings';
import GenerateButton from './components/GenerateButton';
import DependencyInfo from './components/DependencyInfo';
import NotionChatbot from './components/NotionChatbot';
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
  const [currentPage, setCurrentPage] = useState('generator');

  if (currentPage === 'chatbot') {
    return (
      <div className="app">
        <nav className="app-nav">
          <button className={currentPage === 'generator' ? 'active' : ''} onClick={() => setCurrentPage('generator')}>
            Script Generator
          </button>
          <button className={currentPage === 'chatbot' ? 'active' : ''} onClick={() => setCurrentPage('chatbot')}>
            Research Chatbot
          </button>
        </nav>
        <NotionChatbot />
      </div>
    );
  }

  return (
    <div className="app">
      <nav className="app-nav">
        <button className={currentPage === 'generator' ? 'active' : ''} onClick={() => setCurrentPage('generator')}>
          Script Generator
        </button>
        <button className={currentPage === 'chatbot' ? 'active' : ''} onClick={() => setCurrentPage('chatbot')}>
          Research Chatbot
        </button>
      </nav>

      <header className="app-header">
        <div className="header-content">
          <h1>DFT / MLIP Script Generator</h1>
          <p className="subtitle">
            Quantum ESPRESSO &middot; ORCA &middot; UMA &middot; MACE
          </p>
        </div>
      </header>

      <main className="main-content">
        <CifUploader cifText={cifText} setCifText={setCifText} />
        <CodeSelector selectedCode={selectedCode} setSelectedCode={setSelectedCode} />
        <PseudopotentialSelector params={params} setParams={setParams} selectedCode={selectedCode} />
        <CalcTypeSelector calcTypes={calcTypes} setCalcTypes={setCalcTypes} />
        <PostProcessing selectedCode={selectedCode} postProcessing={postProcessing} setPostProcessing={setPostProcessing} />
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
