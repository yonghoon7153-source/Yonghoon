"""
DEM Analysis Web Application
- Single mode: Upload one case → full analysis pipeline → figures + MD report
- Group mode: Upload multiple cases → comparison plots + summary report
"""
import os
import json
import uuid
import shutil
import subprocess
import glob as globmod
from datetime import datetime
from pathlib import Path

from flask import (
    Flask, render_template, request, jsonify, send_from_directory,
    redirect, url_for, send_file
)

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB max
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(__file__), 'uploads')
app.config['RESULTS_FOLDER'] = os.path.join(os.path.dirname(__file__), 'results')
app.config['SCRIPTS_FOLDER'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'scripts')
app.config['ARCHIVE_FOLDER'] = os.path.join(os.path.dirname(__file__), 'archive')

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)
os.makedirs(app.config['ARCHIVE_FOLDER'], exist_ok=True)

# ─── Helpers ────────────────────────────────────────────────────────────────

def get_case_dir(case_id):
    d = os.path.join(app.config['UPLOAD_FOLDER'], case_id)
    os.makedirs(d, exist_ok=True)
    return d

def get_results_dir(case_id):
    d = os.path.join(app.config['RESULTS_FOLDER'], case_id)
    os.makedirs(d, exist_ok=True)
    return d

def detect_mode(case_dir):
    """Detect if bimodal (3 types) or standard (2 types) from input script."""
    for f in os.listdir(case_dir):
        if f.endswith('.liggghts') and 'input' in f.lower():
            with open(os.path.join(case_dir, f)) as fh:
                content = fh.read()
                if 'AM_P' in content or 'type 3' in content.lower():
                    return 'bimodal'
    # Check atom file for 3 types
    for f in os.listdir(case_dir):
        if f.startswith('atom') and f.endswith('.liggghts'):
            with open(os.path.join(case_dir, f)) as fh:
                lines = fh.readlines()
                types = set()
                for line in lines:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            t = int(parts[1])
                            types.add(t)
                        except ValueError:
                            continue
                if len(types) >= 3:
                    return 'bimodal'
    return 'standard'

def list_cases():
    """List all uploaded cases with metadata."""
    cases = []
    upload_dir = app.config['UPLOAD_FOLDER']
    if not os.path.exists(upload_dir):
        return cases
    for case_id in sorted(os.listdir(upload_dir), reverse=True):
        case_dir = os.path.join(upload_dir, case_id)
        if not os.path.isdir(case_dir):
            continue
        meta_file = os.path.join(case_dir, 'meta.json')
        if os.path.exists(meta_file):
            with open(meta_file) as f:
                meta = json.load(f)
        else:
            meta = {'name': case_id, 'created': '', 'mode': 'unknown', 'status': 'uploaded'}
        meta['id'] = case_id
        # Check for results
        results_dir = os.path.join(app.config['RESULTS_FOLDER'], case_id)
        meta['has_results'] = os.path.isdir(results_dir) and len(os.listdir(results_dir)) > 0
        figures_dir = os.path.join(results_dir, 'figures')
        meta['has_figures'] = os.path.isdir(figures_dir) and len(globmod.glob(os.path.join(figures_dir, '*.png'))) > 0
        report_file = os.path.join(results_dir, 'report.md')
        meta['has_report'] = os.path.exists(report_file)
        cases.append(meta)
    return cases

def run_pipeline(case_id, mode, type_map, scale=1000):
    """Run the DEM analysis pipeline for a case."""
    case_dir = get_case_dir(case_id)
    results_dir = get_results_dir(case_id)
    figures_dir = os.path.join(results_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    scripts = app.config['SCRIPTS_FOLDER']

    # Find atom, contact, mesh, and input files
    atom_files = sorted(globmod.glob(os.path.join(case_dir, 'atom_*.liggghts')))
    contact_files = sorted(globmod.glob(os.path.join(case_dir, 'contact_*.liggghts')))
    mesh_files = sorted(globmod.glob(os.path.join(case_dir, '*.stl')))
    input_files = sorted(globmod.glob(os.path.join(case_dir, 'input_*.liggghts')))

    if not atom_files or not contact_files:
        return {'error': 'atom_*.liggghts 또는 contact_*.liggghts 파일을 찾을 수 없습니다.'}

    log = []

    # Step 1: Parse (atom + contact + mesh)
    cmd = ['python3', os.path.join(scripts, 'parse_liggghts.py')]
    cmd += atom_files + contact_files + mesh_files + input_files + ['-o', results_dir]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    log.append({'step': 'Parse', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})
    if result.returncode != 0:
        return {'error': f'Parse failed: {result.stderr}', 'log': log}

    atoms_csv = os.path.join(results_dir, 'atoms.csv')
    contacts_csv = os.path.join(results_dir, 'contacts.csv')

    if mode == 'bimodal':
        # Step 2: Bimodal contact analysis
        cmd = ['python3', os.path.join(scripts, 'analyze_contacts_bimodal.py'),
               atoms_csv, contacts_csv, '-o', results_dir,
               '-t', type_map, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Bimodal Contact Analysis', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

        # Step 3: Basic figures
        cmd = ['python3', os.path.join(scripts, 'generate_figures_bimodal.py'),
               results_dir, '-o', figures_dir, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Basic Figures', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

        # Step 4: Advanced
        atoms_analyzed = os.path.join(results_dir, 'atoms_analyzed.csv')
        contacts_analyzed = os.path.join(results_dir, 'contacts_analyzed.csv')
        cmd = ['python3', os.path.join(scripts, 'advanced_analysis_bimodal.py'),
               atoms_analyzed, contacts_analyzed, '-o', results_dir, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Advanced Analysis', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

        cmd = ['python3', os.path.join(scripts, 'generate_advanced_figures_bimodal.py'),
               results_dir, '-o', figures_dir, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Advanced Figures', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

        # Step 5: Bimodal specific
        cmd = ['python3', os.path.join(scripts, 'bimodal_specific_analysis.py'),
               atoms_analyzed, contacts_analyzed, '-o', results_dir, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Bimodal Specific', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

    else:
        # Standard mode
        cmd = ['python3', os.path.join(scripts, 'analyze_contacts.py'),
               atoms_csv, contacts_csv, '-o', results_dir,
               '-t', type_map, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Contact Analysis', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

        cmd = ['python3', os.path.join(scripts, 'generate_figures.py'),
               results_dir, '-o', figures_dir, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Basic Figures', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

        atoms_analyzed = os.path.join(results_dir, 'atoms_analyzed.csv')
        contacts_analyzed = os.path.join(results_dir, 'contacts_analyzed.csv')

        cmd = ['python3', os.path.join(scripts, 'advanced_analysis.py'),
               atoms_analyzed, contacts_analyzed, '-o', results_dir, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Advanced Analysis', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

        cmd = ['python3', os.path.join(scripts, 'generate_advanced_figures.py'),
               results_dir, '-o', figures_dir, '-s', str(scale)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        log.append({'step': 'Advanced Figures', 'stdout': result.stdout, 'stderr': result.stderr, 'rc': result.returncode})

    return {'success': True, 'log': log}

def generate_report(case_id, case_name='', notes=''):
    """Generate markdown report from analysis results."""
    results_dir = get_results_dir(case_id)
    figures_dir = os.path.join(results_dir, 'figures')

    now = datetime.now().strftime('%y%m%d')
    title = case_name or case_id

    lines = []
    lines.append(f'# {now}_{title}_DEM_Analysis\n')
    lines.append(f'> **날짜**: {datetime.now().strftime("%Y-%m-%d")}')
    lines.append(f'> **케이스**: {title}')
    if notes:
        lines.append(f'> **메모**: {notes}')
    lines.append('')
    lines.append('---\n')

    # Load summary CSV if exists
    summary_file = os.path.join(results_dir, 'contact_summary.csv')
    if os.path.exists(summary_file):
        import pandas as pd
        df = pd.read_csv(summary_file)
        lines.append('## 1. Contact Summary\n')
        lines.append(df.to_markdown(index=False))
        lines.append('')

    # Load atom stats
    atom_stats = os.path.join(results_dir, 'atom_statistics.csv')
    if os.path.exists(atom_stats):
        import pandas as pd
        df = pd.read_csv(atom_stats)
        lines.append('## 2. Particle Statistics\n')
        lines.append(df.to_markdown(index=False))
        lines.append('')

    # Figures
    if os.path.isdir(figures_dir):
        pngs = sorted(globmod.glob(os.path.join(figures_dir, '*.png')))
        if pngs:
            lines.append('## 3. Generated Figures\n')
            for png in pngs:
                fname = os.path.basename(png)
                lines.append(f'### {fname}\n')
                lines.append(f'![{fname}](figures/{fname})\n')

    # Tags
    lines.append('---\n')
    lines.append('#DEM #analysis #ASSB #composite-cathode\n')

    report = '\n'.join(lines)
    report_path = os.path.join(results_dir, 'report.md')
    with open(report_path, 'w') as f:
        f.write(report)
    return report

# ─── Routes ─────────────────────────────────────────────────────────────────

@app.route('/')
def index():
    cases = list_cases()
    return render_template('index.html', cases=cases)

@app.route('/upload', methods=['POST'])
def upload():
    """Upload files for a new case."""
    case_name = request.form.get('case_name', '').strip()
    mode = request.form.get('mode', 'auto')
    type_map = request.form.get('type_map', '')
    ps_ratio = request.form.get('ps_ratio', '').strip()
    scale = request.form.get('scale', '1000')

    case_id = datetime.now().strftime('%y%m%d_%H%M%S') + '_' + str(uuid.uuid4())[:6]
    case_dir = get_case_dir(case_id)

    files = request.files.getlist('files')
    if not files:
        return jsonify({'error': '파일을 선택해주세요.'}), 400

    filenames = []
    for f in files:
        if f.filename:
            safe_name = f.filename.replace('/', '_').replace('\\', '_')
            f.save(os.path.join(case_dir, safe_name))
            filenames.append(safe_name)

    # Detect mode
    if mode == 'auto':
        mode = detect_mode(case_dir)

    # Default type maps
    if not type_map:
        type_map = '1:AM_P,2:AM_S,3:SE' if mode == 'bimodal' else '1:AM,2:SE'

    meta = {
        'name': case_name or case_id,
        'created': datetime.now().isoformat(),
        'mode': mode,
        'type_map': type_map,
        'ps_ratio': ps_ratio,
        'scale': int(scale),
        'files': filenames,
        'status': 'uploaded'
    }
    with open(os.path.join(case_dir, 'meta.json'), 'w') as f:
        json.dump(meta, f, indent=2)

    return jsonify({'case_id': case_id, 'mode': mode, 'files': filenames})

@app.route('/analyze/<case_id>', methods=['POST'])
def analyze(case_id):
    """Run analysis pipeline for a case."""
    case_dir = get_case_dir(case_id)
    meta_file = os.path.join(case_dir, 'meta.json')
    if not os.path.exists(meta_file):
        return jsonify({'error': '케이스를 찾을 수 없습니다.'}), 404

    with open(meta_file) as f:
        meta = json.load(f)

    meta['status'] = 'running'
    with open(meta_file, 'w') as f:
        json.dump(meta, f, indent=2)

    result = run_pipeline(case_id, meta['mode'], meta['type_map'], meta.get('scale', 1000))

    meta['status'] = 'done' if result.get('success') else 'error'
    meta['analysis_log'] = result.get('log', [])
    with open(meta_file, 'w') as f:
        json.dump(meta, f, indent=2)

    # Generate report
    if result.get('success'):
        generate_report(case_id, meta.get('name', ''))

    return jsonify(result)

@app.route('/single/<case_id>')
def single(case_id):
    """View single case results."""
    case_dir = get_case_dir(case_id)
    results_dir = get_results_dir(case_id)
    meta_file = os.path.join(case_dir, 'meta.json')

    if not os.path.exists(meta_file):
        return redirect(url_for('index'))

    with open(meta_file) as f:
        meta = json.load(f)
    meta['id'] = case_id

    # Collect figures
    figures = []
    figures_dir = os.path.join(results_dir, 'figures')
    if os.path.isdir(figures_dir):
        for png in sorted(globmod.glob(os.path.join(figures_dir, '*.png'))):
            figures.append(os.path.basename(png))

    # Load report
    report = ''
    report_path = os.path.join(results_dir, 'report.md')
    if os.path.exists(report_path):
        with open(report_path) as f:
            report = f.read()

    # Load CSVs for tables (atom_statistics first, no force_summary)
    tables = {}
    for csv_name in ['atom_statistics', 'contact_summary', 'coordination_summary',
                     'network_summary']:
        csv_path = os.path.join(results_dir, f'{csv_name}.csv')
        if os.path.exists(csv_path):
            import pandas as pd
            df = pd.read_csv(csv_path)
            tables[csv_name] = {
                'columns': df.columns.tolist(),
                'data': df.values.tolist()
            }

    # Load full_metrics.json for header info
    metrics = {}
    metrics_path = os.path.join(results_dir, 'full_metrics.json')
    if os.path.exists(metrics_path):
        with open(metrics_path) as f:
            metrics = json.load(f)

    # Load input_params.json
    input_params = {}
    params_path = os.path.join(results_dir, 'input_params.json')
    if os.path.exists(params_path):
        with open(params_path) as f:
            input_params = json.load(f)

    return render_template('single.html', case=meta, figures=figures,
                         report=report, tables=tables, metrics=metrics,
                         input_params=input_params)

@app.route('/group', methods=['GET', 'POST'])
def group():
    """Group comparison page."""
    cases = list_cases()
    selected = request.args.getlist('cases')

    comparison_data = {}
    if selected:
        # Key metrics for comparison (from full_metrics.json)
        display_keys = [
            ('P:S', '', 'ps_ratio'),
            ('Porosity', '(%)', 'porosity'),
            ('두께', '(μm)', 'thickness_um'),
            ('AM-SE Total', '(μm²)', 'area_AM전체_SE_total'),
            ('SE-SE Total', '(μm²)', 'area_SE_SE_total'),
            ('SE-SE CN', '', 'se_se_cn'),
            ('SE Cluster', '', 'n_components'),
            ('Percolation', '(%)', 'percolation_pct'),
            ('Top Reachable', '(%)', 'top_reachable_pct'),
            ('Tortuosity', '', 'tortuosity_mean'),
            ('Ionic Active', '(%)', 'ionic_active_pct'),
            ('Coverage P', '(%)', 'coverage_AM_P_mean'),
            ('Coverage S', '(%)', 'coverage_AM_S_mean'),
        ]
        rows = []
        for cid in selected:
            results_dir = get_results_dir(cid)
            meta_file = os.path.join(get_case_dir(cid), 'meta.json')
            if not os.path.exists(meta_file):
                continue
            with open(meta_file) as f:
                meta = json.load(f)

            metrics_path = os.path.join(results_dir, 'full_metrics.json')
            metrics = {}
            if os.path.exists(metrics_path):
                with open(metrics_path) as f:
                    metrics = json.load(f)

            # Standard 모드: coverage_AM_mean → P:S에 따라 P 또는 S에 매핑
            if 'coverage_AM_mean' in metrics:
                ps = metrics.get('ps_ratio', '')
                if ps in ('P only', '10:0'):
                    metrics.setdefault('coverage_AM_P_mean', metrics['coverage_AM_mean'])
                else:
                    metrics.setdefault('coverage_AM_S_mean', metrics['coverage_AM_mean'])
            # Standard 모드: area_AM_SE_total → P:S에 따라 매핑
            if 'area_AM_SE_total' in metrics and 'area_AM_P_SE_total' not in metrics:
                ps = metrics.get('ps_ratio', '')
                if ps in ('P only', '10:0'):
                    metrics['area_AM_P_SE_total'] = metrics['area_AM_SE_total']
                else:
                    metrics['area_AM_S_SE_total'] = metrics['area_AM_SE_total']

            row = {'케이스': meta.get('name', cid)}
            for label, unit, key in display_keys:
                val = metrics.get(key, '')
                if isinstance(val, float):
                    val = round(val, 2)
                row[label] = val if val != '' else '-'
            rows.append(row)

        if rows:
            # Build columns with unit subtitles
            col_headers = [{'name': '케이스', 'unit': ''}]
            for label, unit, key in display_keys:
                col_headers.append({'name': label, 'unit': unit})

            comparison_data = {
                'columns': [c['name'] for c in col_headers],
                'units': [c['unit'] for c in col_headers],
                'rows': rows
            }

    return render_template('group.html', cases=cases, selected=selected,
                         comparison=comparison_data)

@app.route('/group/plots', methods=['POST'])
def group_plots():
    """Generate comparison plots for selected cases."""
    selected = request.form.getlist('cases')
    plots = request.form.getlist('plots')
    if not selected or not plots:
        return jsonify({'error': '케이스와 플롯을 선택하세요.'}), 400

    session_id = datetime.now().strftime('%y%m%d_%H%M%S') + '_' + str(uuid.uuid4())[:4]
    plot_dir = os.path.join(app.config['RESULTS_FOLDER'], 'group_plots', session_id)
    os.makedirs(plot_dir, exist_ok=True)

    # Collect metrics files and names
    input_files = []
    names = []
    for cid in selected:
        results_dir = get_results_dir(cid)
        metrics_path = os.path.join(results_dir, 'full_metrics.json')
        meta_file = os.path.join(get_case_dir(cid), 'meta.json')
        if os.path.exists(metrics_path) and os.path.exists(meta_file):
            input_files.append(metrics_path)
            with open(meta_file) as f:
                meta = json.load(f)
            names.append(meta.get('name', cid))

    if not input_files:
        return jsonify({'error': '메트릭 데이터가 없습니다.'}), 400

    scripts = app.config['SCRIPTS_FOLDER']
    cmd = ['python3', os.path.join(scripts, 'generate_comparison_plots.py'),
           '-i'] + input_files + ['-n'] + names + ['-o', plot_dir, '-p'] + plots
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    if result.returncode != 0:
        return jsonify({'error': f'Plot 생성 실패: {result.stderr}'}), 500

    # Load plot info
    info_path = os.path.join(plot_dir, 'plot_info.json')
    plot_list = []
    if os.path.exists(info_path):
        with open(info_path) as f:
            info = json.load(f)
        for key in plots:
            if key in info:
                plot_list.append(info[key])

    return jsonify({'session': session_id, 'plots': plot_list})

@app.route('/group/plot-image/<session>/<filename>')
def serve_group_plot(session, filename):
    plot_dir = os.path.join(app.config['RESULTS_FOLDER'], 'group_plots', session)
    return send_from_directory(plot_dir, filename)

@app.route('/group/report', methods=['POST'])
def group_report():
    """Generate group comparison markdown report."""
    selected = request.form.getlist('cases')
    title = request.form.get('title', 'Group Comparison')
    notes = request.form.get('notes', '')

    now = datetime.now().strftime('%y%m%d')
    lines = []
    lines.append(f'# {now}_{title}\n')
    lines.append(f'> **날짜**: {datetime.now().strftime("%Y-%m-%d")}')
    lines.append(f'> **비교 케이스**: {len(selected)}개')
    if notes:
        lines.append(f'> **메모**: {notes}')
    lines.append('\n---\n')

    import pandas as pd
    all_summaries = []
    for cid in selected:
        results_dir = get_results_dir(cid)
        meta_file = os.path.join(get_case_dir(cid), 'meta.json')
        if not os.path.exists(meta_file):
            continue
        with open(meta_file) as f:
            meta = json.load(f)

        summary = {'Case': meta.get('name', cid)}
        for csv_name in ['contact_summary', 'coordination_summary']:
            csv_path = os.path.join(results_dir, f'{csv_name}.csv')
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                for col in df.columns:
                    summary[col] = df[col].iloc[0] if len(df) == 1 else df[col].tolist()
        all_summaries.append(summary)

    if all_summaries:
        df_comp = pd.DataFrame(all_summaries)
        lines.append('## Comparison Table\n')
        lines.append(df_comp.to_markdown(index=False))
        lines.append('')

    lines.append('\n---\n')
    lines.append('#DEM #group-comparison #ASSB\n')

    report = '\n'.join(lines)

    # Save report
    group_dir = os.path.join(app.config['RESULTS_FOLDER'], 'group_reports')
    os.makedirs(group_dir, exist_ok=True)
    report_path = os.path.join(group_dir, f'{now}_{title.replace(" ", "_")}.md')
    with open(report_path, 'w') as f:
        f.write(report)

    return jsonify({'report': report, 'path': report_path})

@app.route('/results/<case_id>/figures/<filename>')
def serve_figure(case_id, filename):
    figures_dir = os.path.join(get_results_dir(case_id), 'figures')
    return send_from_directory(figures_dir, filename)

@app.route('/results/<case_id>/report')
def serve_report(case_id):
    report_path = os.path.join(get_results_dir(case_id), 'report.md')
    if os.path.exists(report_path):
        return send_file(report_path, mimetype='text/markdown',
                        as_attachment=True,
                        download_name=f'{case_id}_report.md')
    return 'Report not found', 404

@app.route('/rename/<case_id>', methods=['POST'])
def rename_case(case_id):
    """Rename a case."""
    new_name = request.form.get('name', '').strip()
    if not new_name:
        return jsonify({'error': '이름을 입력하세요.'}), 400
    case_dir = get_case_dir(case_id)
    meta_file = os.path.join(case_dir, 'meta.json')
    if not os.path.exists(meta_file):
        return jsonify({'error': '케이스를 찾을 수 없습니다.'}), 404
    with open(meta_file) as f:
        meta = json.load(f)
    meta['name'] = new_name
    with open(meta_file, 'w') as f:
        json.dump(meta, f, indent=2)
    return jsonify({'success': True})

@app.route('/delete/<case_id>', methods=['POST'])
def delete_case(case_id):
    case_dir = get_case_dir(case_id)
    results_dir = get_results_dir(case_id)
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    if os.path.exists(results_dir):
        shutil.rmtree(results_dir)
    return jsonify({'success': True})

# ─── Archive (보관함) ───────────────────────────────────────────────────────

def _archive_root():
    return app.config['ARCHIVE_FOLDER']

def _safe_path(rel):
    """Prevent path traversal."""
    base = os.path.realpath(_archive_root())
    target = os.path.realpath(os.path.join(base, rel))
    if not target.startswith(base):
        return None
    return target

def _scan_folder(abs_path, rel_prefix=''):
    """Recursively scan a folder and return tree structure."""
    items = []
    if not os.path.isdir(abs_path):
        return items
    for name in sorted(os.listdir(abs_path)):
        full = os.path.join(abs_path, name)
        rel = os.path.join(rel_prefix, name) if rel_prefix else name
        if os.path.isdir(full):
            children = _scan_folder(full, rel)
            file_count = sum(1 for c in children if c['type'] == 'file') + \
                         sum(c.get('file_count', 0) for c in children if c['type'] == 'folder')
            items.append({
                'type': 'folder', 'name': name, 'path': rel,
                'children': children, 'file_count': file_count
            })
        else:
            size = os.path.getsize(full)
            items.append({
                'type': 'file', 'name': name, 'path': rel,
                'size': size, 'ext': os.path.splitext(name)[1].lower()
            })
    return items

@app.route('/archive')
def archive():
    tree = _scan_folder(_archive_root())
    current = request.args.get('folder', '')
    # List files in current folder
    if current:
        target = _safe_path(current)
        if not target or not os.path.isdir(target):
            current = ''
    folder_path = _safe_path(current) if current else _archive_root()
    files = []
    folders = []
    for name in sorted(os.listdir(folder_path)):
        full = os.path.join(folder_path, name)
        rel = os.path.join(current, name) if current else name
        if os.path.isdir(full):
            cnt = len([f for f in os.listdir(full) if os.path.isfile(os.path.join(full, f))])
            sub = len([f for f in os.listdir(full) if os.path.isdir(os.path.join(full, f))])
            folders.append({'name': name, 'path': rel, 'file_count': cnt, 'subfolder_count': sub})
        else:
            size = os.path.getsize(full)
            files.append({'name': name, 'path': rel, 'size': size,
                         'ext': os.path.splitext(name)[1].lower()})

    # Breadcrumb
    breadcrumb = []
    if current:
        parts = current.split(os.sep)
        for i, p in enumerate(parts):
            breadcrumb.append({'name': p, 'path': os.sep.join(parts[:i+1])})

    return render_template('archive.html', tree=tree, folders=folders, files=files,
                         current=current, breadcrumb=breadcrumb)

@app.route('/archive/create-folder', methods=['POST'])
def archive_create_folder():
    parent = request.form.get('parent', '')
    name = request.form.get('name', '').strip()
    if not name:
        return jsonify({'error': '폴더 이름을 입력하세요.'}), 400
    # Sanitize
    name = name.replace('/', '_').replace('\\', '_').replace('..', '')
    base = _safe_path(parent) if parent else _archive_root()
    if not base:
        return jsonify({'error': '잘못된 경로입니다.'}), 400
    target = os.path.join(base, name)
    os.makedirs(target, exist_ok=True)
    return jsonify({'success': True, 'path': os.path.join(parent, name) if parent else name})

@app.route('/archive/delete', methods=['POST'])
def archive_delete():
    path = request.form.get('path', '')
    target = _safe_path(path)
    if not target or target == os.path.realpath(_archive_root()):
        return jsonify({'error': '삭제할 수 없습니다.'}), 400
    if os.path.isdir(target):
        shutil.rmtree(target)
    elif os.path.isfile(target):
        os.remove(target)
    return jsonify({'success': True})

@app.route('/archive/rename', methods=['POST'])
def archive_rename():
    old_path = request.form.get('path', '')
    new_name = request.form.get('new_name', '').strip()
    if not new_name:
        return jsonify({'error': '새 이름을 입력하세요.'}), 400
    new_name = new_name.replace('/', '_').replace('\\', '_').replace('..', '')
    target = _safe_path(old_path)
    if not target:
        return jsonify({'error': '잘못된 경로입니다.'}), 400
    parent = os.path.dirname(target)
    new_full = os.path.join(parent, new_name)
    os.rename(target, new_full)
    return jsonify({'success': True})

@app.route('/archive/move', methods=['POST'])
def archive_move():
    src = request.form.get('src', '')
    dst_folder = request.form.get('dst', '')
    src_full = _safe_path(src)
    dst_full = _safe_path(dst_folder) if dst_folder else _archive_root()
    if not src_full or not dst_full:
        return jsonify({'error': '잘못된 경로입니다.'}), 400
    name = os.path.basename(src_full)
    shutil.move(src_full, os.path.join(dst_full, name))
    return jsonify({'success': True})

@app.route('/archive/save-case', methods=['POST'])
def archive_save_case():
    """Save a case's results to archive folder."""
    case_id = request.form.get('case_id', '')
    folder = request.form.get('folder', '')
    results_dir = get_results_dir(case_id)
    case_dir = get_case_dir(case_id)
    meta_file = os.path.join(case_dir, 'meta.json')

    if not os.path.exists(meta_file):
        return jsonify({'error': '케이스를 찾을 수 없습니다.'}), 404

    with open(meta_file) as f:
        meta = json.load(f)

    case_name = meta.get('name', case_id)
    dst = _safe_path(folder) if folder else _archive_root()
    if not dst:
        return jsonify({'error': '잘못된 경로입니다.'}), 400

    save_dir = os.path.join(dst, case_name)
    if os.path.exists(save_dir):
        save_dir = save_dir + '_' + datetime.now().strftime('%H%M%S')
    shutil.copytree(results_dir, save_dir, dirs_exist_ok=True)
    # Also copy meta
    shutil.copy2(meta_file, os.path.join(save_dir, 'meta.json'))
    return jsonify({'success': True, 'saved_to': save_dir})

@app.route('/archive/download/<path:filepath>')
def archive_download(filepath):
    target = _safe_path(filepath)
    if not target or not os.path.isfile(target):
        return 'File not found', 404
    return send_file(target, as_attachment=True)

@app.route('/archive/preview/<path:filepath>')
def archive_preview(filepath):
    target = _safe_path(filepath)
    if not target or not os.path.isfile(target):
        return 'File not found', 404
    ext = os.path.splitext(target)[1].lower()
    if ext == '.md':
        with open(target) as f:
            content = f.read()
        return jsonify({'type': 'markdown', 'content': content})
    elif ext == '.csv':
        import pandas as pd
        df = pd.read_csv(target)
        return jsonify({'type': 'csv', 'columns': df.columns.tolist(),
                       'data': df.head(100).values.tolist()})
    elif ext == '.png':
        return send_file(target, mimetype='image/png')
    elif ext == '.json':
        with open(target) as f:
            content = f.read()
        return jsonify({'type': 'json', 'content': content})
    return jsonify({'type': 'unknown', 'message': '미리보기를 지원하지 않는 파일 형식입니다.'})

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
