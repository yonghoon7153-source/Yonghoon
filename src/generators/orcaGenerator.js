export function generateORCAScripts(parsed, params, calcTypes, postProcessing) {
  const files = {};

  const functional = params.functional || 'PBE';
  const basisSet = params.basisSet || 'def2-SVP';
  const charge = params.charge || 0;
  const multiplicity = params.multiplicity || 1;
  const nprocs = params.nprocs || 4;
  const dftd3 = params.dftd3 ? ' D3BJ' : '';

  // Convert fractional coords to cartesian (approximate for molecules)
  const { v1, v2, v3 } = parsed.cellParams;
  const cartesianAtoms = parsed.atoms.map(a => {
    const x = a.x * v1[0] + a.y * v2[0] + a.z * v3[0];
    const y = a.x * v1[1] + a.y * v2[1] + a.z * v3[1];
    const z = a.x * v1[2] + a.y * v2[2] + a.z * v3[2];
    return { symbol: a.symbol, x, y, z };
  });

  const coordBlock = cartesianAtoms.map(a =>
    `  ${a.symbol.padEnd(4)} ${a.x.toFixed(8)}  ${a.y.toFixed(8)}  ${a.z.toFixed(8)}`
  ).join('\n');

  const palBlock = `%pal
  nprocs ${nprocs}
end`;

  // --- SCF ---
  if (calcTypes.includes('scf')) {
    files['scf.inp'] = `# ORCA SCF calculation
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF PrintBasis

${palBlock}

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // --- Geometry Optimization ---
  if (calcTypes.includes('relax') || calcTypes.includes('vc-relax')) {
    files['opt.inp'] = `# ORCA Geometry Optimization
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF Opt

${palBlock}

%geom
  MaxIter 200
end

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // --- Frequency (IR/Raman) ---
  if (postProcessing.includes('ir') || postProcessing.includes('raman')) {
    files['freq.inp'] = `# ORCA Frequency calculation (IR${postProcessing.includes('raman') ? '/Raman' : ''})
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF Freq${postProcessing.includes('raman') ? '' : ' NumFreq'}

${palBlock}

%elprop
  Polar 1
end

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // --- DOS/PDOS ---
  if (postProcessing.includes('dos') || postProcessing.includes('pdos')) {
    files['dos.inp'] = `# ORCA DOS calculation
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF

${palBlock}

%output
  Print[P_MOs] 1
  Print[P_Overlap] 1
end

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // --- Optical / UV-Vis ---
  if (postProcessing.includes('optical_absorption')) {
    files['tddft.inp'] = `# ORCA TD-DFT calculation (optical absorption)
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF

${palBlock}

%tddft
  NRoots 30
  MaxDim 100
end

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // --- XPS ---
  if (postProcessing.includes('xps')) {
    files['xps.inp'] = `# ORCA XPS core-level calculation
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF

${palBlock}

%method
  CoreLevel true
end

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // --- Magnetic moments ---
  if (postProcessing.includes('magnetic_moments')) {
    files['mag.inp'] = `# ORCA Magnetic properties
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF UKS

${palBlock}

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // --- Bonding energy ---
  if (postProcessing.includes('bonding_energy')) {
    files['energy.inp'] = `# ORCA Single-point energy
! ${functional.toUpperCase()} ${basisSet}${dftd3} TightSCF

${palBlock}

* xyz ${charge} ${multiplicity}
${coordBlock}
*
`;
  }

  // Run script
  files['run_all.sh'] = `#!/bin/bash
# ORCA Workflow Script
set -e
ORCA_PATH=\${ORCA_PATH:-orca}

${Object.keys(files).filter(f => f.endsWith('.inp')).map(f => {
  const name = f.replace('.inp', '');
  return `echo "Running ${name}..."
\${ORCA_PATH} ${f} > ${name}.out 2>&1
echo "  ${name} completed."`;
}).join('\n\n')}

echo "All ORCA calculations completed!"
`;

  return files;
}
