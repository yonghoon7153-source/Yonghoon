// All supported calculation codes
export const CODES = [
  { id: 'vasp', label: 'VASP' },
  { id: 'qe', label: 'Quantum ESPRESSO' },
  { id: 'orca', label: 'ORCA' },
  { id: 'uma', label: 'UMA (MLIP)' },
  { id: 'mace', label: 'MACE (MLIP)' },
];

// Basic calculation types
export const CALC_TYPES = [
  { id: 'scf', label: 'SCF', description: 'Self-consistent field calculation' },
  { id: 'nscf', label: 'NSCF', description: 'Non-self-consistent field calculation' },
  { id: 'relax', label: 'Relaxation', description: 'Atomic position optimization' },
  { id: 'vc-relax', label: 'VC-Relax', description: 'Variable-cell relaxation (cell + atomic positions)' },
];

// Pseudopotential libraries
export const PP_LIBRARIES = [
  { id: 'sssp_efficiency', label: 'SSSP Efficiency' },
  { id: 'sssp_precision', label: 'SSSP Precision' },
  { id: 'pslibrary_paw', label: 'PSlibrary (PAW)' },
  { id: 'pslibrary_us', label: 'PSlibrary (US)' },
  { id: 'gbrv', label: 'GBRV' },
  { id: 'custom', label: 'Custom path' },
];

// Exchange-correlation functionals
export const FUNCTIONALS = [
  { id: 'pbe', label: 'PBE' },
  { id: 'pbesol', label: 'PBEsol' },
  { id: 'lda', label: 'LDA' },
  { id: 'hse06', label: 'HSE06 (hybrid)' },
  { id: 'b3lyp', label: 'B3LYP (hybrid)' },
  { id: 'scan', label: 'SCAN (meta-GGA)' },
];

// Post-processing categories and items
export const POST_PROCESSING = [
  {
    category: '전자구조 / 밴드 (Electronic Structure / Bands)',
    items: [
      { id: 'dos', label: 'DOS', description: 'Density of states — 전체 상태밀도 계산', codes: ['vasp', 'qe', 'orca'] },
      { id: 'pdos', label: 'PDOS', description: 'Projected density of states — 원자/오비탈별 분해 상태밀도', codes: ['vasp', 'qe', 'orca'] },
      { id: 'bands', label: 'Band structure', description: '에너지 밴드 구조 계산 및 플롯', codes: ['vasp', 'qe'] },
      { id: 'bandgap', label: 'Band gap extraction', description: 'VBM/CBM에서 밴드갭 추출', codes: ['vasp', 'qe'] },
      { id: 'effmass', label: 'Effective mass', description: '밴드 곡률로부터 유효질량 계산', codes: ['vasp', 'qe'] },
      { id: 'fermi_surface', label: 'Fermi surface', description: '페르미 면 계산 및 시각화', codes: ['vasp', 'qe'] },
      { id: 'wannier', label: 'Wannier interpolation (Wannier90)', description: 'Wannier 함수 기반 밴드 보간', codes: ['vasp', 'qe'] },
    ],
  },
  {
    category: '결합 분석 (Bonding Analysis)',
    items: [
      { id: 'bader', label: 'Bader charge analysis', description: 'Bader 전하 분석 (원자별 전하 할당)', codes: ['vasp', 'qe'] },
      { id: 'cohp', label: 'COHP', description: 'Crystal Orbital Hamilton Population — 결합/반결합 분석', codes: ['vasp', 'qe'] },
      { id: 'coop', label: 'COOP', description: 'Crystal Orbital Overlap Population', codes: ['vasp', 'qe'] },
      { id: 'cobi', label: 'COBI', description: 'Crystal Orbital Bond Index', codes: ['vasp', 'qe'] },
      { id: 'bond_length', label: 'Bond length analysis', description: '결합 길이 분석', codes: ['vasp', 'qe', 'orca', 'uma', 'mace'] },
      { id: 'cdd', label: 'Charge density difference (CDD)', description: '전하밀도 차이 시각화', codes: ['vasp', 'qe'] },
    ],
  },
  {
    category: '에너지 / 열역학 (Energy / Thermodynamics)',
    items: [
      { id: 'eos', label: 'EOS (Equation of State)', description: '상태방정식 피팅 → 체적탄성률(B₀) 추출', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'formation_energy', label: 'Formation energy', description: '생성 에너지 계산', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'adhesion_energy', label: 'Adhesion energy', description: '접착 에너지 계산', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'bonding_energy', label: 'Bonding energy', description: '결합 에너지 계산', codes: ['vasp', 'qe', 'orca', 'uma', 'mace'] },
      { id: 'surface_energy', label: 'Surface energy', description: '표면 에너지 계산', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'vacancy_energy', label: 'Vacancy / defect formation energy', description: '공공/결함 생성 에너지', codes: ['vasp', 'qe', 'uma', 'mace'] },
    ],
  },
  {
    category: '분광학 (Spectroscopy)',
    items: [
      { id: 'xps', label: 'XPS', description: 'X-ray photoelectron spectroscopy — 코어 레벨 시뮬레이션', codes: ['vasp', 'qe', 'orca'] },
      { id: 'raman', label: 'Raman spectrum', description: '라만 스펙트럼 계산', codes: ['vasp', 'qe', 'orca'] },
      { id: 'ir', label: 'IR spectrum', description: '적외선 스펙트럼 계산', codes: ['vasp', 'qe', 'orca'] },
    ],
  },
  {
    category: '전자밀도 / 전위 (Electron Density / Potential)',
    items: [
      { id: 'elf', label: 'ELF', description: 'Electron Localization Function — 전자 국소화 함수', codes: ['vasp', 'qe'] },
      { id: 'bvse', label: 'BVSE', description: 'Bond Valence Site Energy — 이온 이동 경로 분석', codes: ['vasp', 'qe'] },
      { id: 'work_function', label: 'Work function', description: '일함수 계산 (슬랩 모델)', codes: ['vasp', 'qe'] },
      { id: 'stm', label: 'STM image simulation', description: 'STM 이미지 시뮬레이션', codes: ['vasp', 'qe'] },
    ],
  },
  {
    category: '광학 물성 (Optical Properties)',
    items: [
      { id: 'optical_absorption', label: 'Optical absorption spectrum', description: '광흡수 스펙트럼', codes: ['vasp', 'qe', 'orca'] },
      { id: 'refractive_index', label: 'Refractive index', description: '굴절률 계산', codes: ['vasp', 'qe'] },
      { id: 'dielectric', label: 'Dielectric function', description: '유전 함수 (VASP: LOPTICS / QE: epsilon.x)', codes: ['vasp', 'qe'] },
    ],
  },
  {
    category: '역학 물성 (Mechanical Properties)',
    items: [
      { id: 'elastic', label: 'Elastic tensor (Cij)', description: '탄성 텐서 → Bulk/Shear/Young modulus, Poisson ratio, anisotropy', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'piezoelectric', label: 'Piezoelectric coefficients', description: '압전 계수 계산', codes: ['vasp', 'qe'] },
    ],
  },
  {
    category: '격자동역학 / 열물성 (Lattice Dynamics / Thermal)',
    items: [
      { id: 'phonon', label: 'Phonon dispersion + phonon DOS', description: '포논 분산 및 포논 상태밀도 (phonopy)', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'thermo', label: 'Thermodynamic properties', description: 'Cv, entropy, Helmholtz free energy', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'thermal_conductivity', label: 'Thermal conductivity', description: '열전도도 (phono3py)', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'born', label: 'Born effective charge / dielectric tensor', description: 'Born 유효 전하 및 유전 텐서 (VASP: LEPSILON)', codes: ['vasp', 'qe'] },
      { id: 'boltztrap', label: 'BoltzTraP', description: 'Seebeck coefficient, electrical conductivity, thermoelectric', codes: ['vasp', 'qe'] },
    ],
  },
  {
    category: '분자동역학 (Molecular Dynamics)',
    items: [
      { id: 'aimd', label: 'AIMD', description: 'Ab initio molecular dynamics', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'rdf', label: 'RDF', description: 'Radial distribution function — 방사 분포 함수', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'msd', label: 'MSD → diffusion → ionic conductivity', description: '평균 제곱 변위 → 확산 계수 → 이온 전도도', codes: ['vasp', 'qe', 'uma', 'mace'] },
      { id: 'stress_fluct', label: 'Stress fluctuation → elastic constants', description: '응력 요동법을 이용한 유한 온도 탄성 상수', codes: ['vasp', 'qe', 'uma', 'mace'] },
    ],
  },
  {
    category: '자성 (Magnetism)',
    items: [
      { id: 'magnetic_moments', label: 'Magnetic moments', description: '자기 모멘트 계산', codes: ['vasp', 'qe', 'orca'] },
      { id: 'mae', label: 'Magnetic anisotropy energy (MAE)', description: '자기 이방성 에너지', codes: ['vasp', 'qe'] },
    ],
  },
  {
    category: '반응경로 (Reaction Pathway)',
    items: [
      { id: 'neb', label: 'NEB', description: 'Nudged Elastic Band — 최소 에너지 경로 및 활성화 장벽', codes: ['vasp', 'qe', 'uma', 'mace'] },
    ],
  },
];

// Dependencies: selecting certain post-processing auto-includes prerequisites
export const DEPENDENCIES = {
  dos: ['scf', 'nscf'],
  pdos: ['scf', 'nscf'],
  bands: ['scf', 'nscf'],
  bandgap: ['scf', 'nscf'],
  effmass: ['scf', 'nscf', 'bands'],
  fermi_surface: ['scf', 'nscf'],
  wannier: ['scf', 'nscf'],
  bader: ['scf'],
  cohp: ['scf', 'nscf'],
  coop: ['scf', 'nscf'],
  cobi: ['scf', 'nscf'],
  cdd: ['scf'],
  eos: ['scf'],
  formation_energy: ['scf'],
  adhesion_energy: ['scf', 'relax'],
  bonding_energy: ['scf'],
  surface_energy: ['scf', 'relax'],
  vacancy_energy: ['scf', 'relax'],
  xps: ['scf'],
  raman: ['scf'],
  ir: ['scf'],
  elf: ['scf'],
  bvse: ['scf'],
  work_function: ['scf'],
  stm: ['scf', 'nscf'],
  optical_absorption: ['scf', 'nscf'],
  refractive_index: ['scf', 'nscf'],
  dielectric: ['scf', 'nscf'],
  elastic: ['scf'],
  piezoelectric: ['scf'],
  phonon: ['scf'],
  thermo: ['scf', 'phonon'],
  thermal_conductivity: ['scf', 'phonon'],
  born: ['scf'],
  boltztrap: ['scf', 'nscf'],
  aimd: ['scf'],
  rdf: ['scf', 'aimd'],
  msd: ['scf', 'aimd'],
  stress_fluct: ['scf', 'aimd'],
  magnetic_moments: ['scf'],
  mae: ['scf'],
  neb: ['scf', 'relax'],
  bond_length: ['scf'],
};
