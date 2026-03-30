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
import threading
from datetime import datetime
from pathlib import Path

# Load .env file if exists (for local development)
_env_path = os.path.join(os.path.dirname(__file__), '.env')
if os.path.exists(_env_path):
    with open(_env_path) as _f:
        for _line in _f:
            _line = _line.strip()
            if _line and not _line.startswith('#') and '=' in _line:
                _k, _v = _line.split('=', 1)
                os.environ.setdefault(_k.strip(), _v.strip())

from flask import (
    Flask, render_template, request, jsonify, send_from_directory,
    redirect, url_for, send_file
)
import storage_sync

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB max
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(__file__), 'uploads')
app.config['RESULTS_FOLDER'] = os.path.join(os.path.dirname(__file__), 'results')
app.config['SCRIPTS_FOLDER'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'scripts')
app.config['ARCHIVE_FOLDER'] = os.path.join(os.path.dirname(__file__), 'archive')

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)
os.makedirs(app.config['ARCHIVE_FOLDER'], exist_ok=True)

# ─── Supabase Storage: restore on startup ─────────────────────────────────
storage_sync.init()
storage_sync.restore_all(
    app.config['UPLOAD_FOLDER'],
    app.config['RESULTS_FOLDER'],
    app.config['ARCHIVE_FOLDER'],
)

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
    """Detect if bimodal (3 types) or standard (2 types) from atom file type count."""
    # Count unique atom types from atom dump file (most reliable)
    for f in sorted(os.listdir(case_dir)):
        if f.startswith('atom') and f.endswith('.liggghts'):
            with open(os.path.join(case_dir, f)) as fh:
                lines = fh.readlines()
                types = set()
                in_data = False
                for line in lines:
                    stripped = line.strip()
                    if stripped.startswith('ITEM: ATOMS'):
                        in_data = True
                        continue
                    if stripped.startswith('ITEM:'):
                        in_data = False
                        continue
                    if in_data:
                        parts = stripped.split()
                        if len(parts) >= 2:
                            try:
                                t = int(parts[1])
                                types.add(t)
                            except ValueError:
                                continue
                # 3 types = bimodal (AM_P + AM_S + SE)
                # 2 types = standard (AM + SE)
                if len(types) >= 3:
                    return 'bimodal'
                return 'standard'
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
        # Check for warnings
        metrics_file = os.path.join(results_dir, 'full_metrics.json')
        if os.path.exists(metrics_file):
            with open(metrics_file) as f:
                m = json.load(f)
            meta['warning_count'] = m.get('warning_count', 0)
            meta['warning_msgs'] = [w['msg'] for w in m.get('warnings', [])]
        cases.append(meta)
    return cases

def run_pipeline(case_id, mode, type_map, scale=1000):
    """Run the DEM analysis pipeline for a case."""
    # Clear pyc cache to ensure latest code runs
    import glob as globmod
    scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'scripts')
    for pyc in globmod.glob(os.path.join(scripts_dir, '__pycache__', '*.pyc')):
        os.remove(pyc)

    case_dir = get_case_dir(case_id)
    results_dir = get_results_dir(case_id)
    figures_dir = os.path.join(results_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    scripts = app.config['SCRIPTS_FOLDER']

    # Find atom, contact, mesh, and input files
    atom_files = sorted(globmod.glob(os.path.join(case_dir, 'atom_*.liggghts')))
    contact_files = sorted(globmod.glob(os.path.join(case_dir, 'contact_*.liggghts')))
    mesh_files = sorted(globmod.glob(os.path.join(case_dir, '*.stl')))
    input_files = sorted(globmod.glob(os.path.join(case_dir, 'input*.liggghts')))

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

# ─── Claude AI Analysis ────────────────────────────────────────────────────

def _generate_ai_analysis(all_metrics, case_names, title, notes):
    """Use Claude API to generate deep analysis of comparison data."""
    api_key = os.environ.get('ANTHROPIC_API_KEY', '')
    if not api_key:
        return None

    try:
        import anthropic
        client = anthropic.Anthropic(api_key=api_key)
    except ImportError:
        return None
    except Exception:
        return None

    # Build data summary for Claude
    import pandas as pd
    display_keys = [
        ('P:S', 'ps_ratio'), ('Porosity(%)', 'porosity'),
        ('Thickness(μm)', 'thickness_um'),
        ('AM-SE Total(μm²)', 'area_AM전체_SE_total'),
        ('SE-SE Total(μm²)', 'area_SE_SE_total'),
        ('SE-SE N', 'area_SE_SE_n'),
        ('SE-SE Mean Area(μm²)', 'area_SE_SE_mean'),
        ('SE-SE CN', 'se_se_cn'),
        ('SE Cluster', 'n_components'),
        ('Percolation(%)', 'percolation_pct'),
        ('Top Reachable(%)', 'top_reachable_pct'),
        ('Tortuosity', 'tortuosity_mean'),
        ('Ionic Active AM(%)', 'ionic_active_pct'),
        ('Coverage AM_P(%)', 'coverage_AM_P_mean'),
        ('Coverage AM_S(%)', 'coverage_AM_S_mean'),
    ]
    rows = []
    for i, name in enumerate(case_names):
        row = {'Case': name}
        for label, key in display_keys:
            val = all_metrics[i].get(key, '-')
            if isinstance(val, float):
                val = round(val, 2)
            row[label] = val
        rows.append(row)
    df = pd.DataFrame(rows)
    data_table = df.to_markdown(index=False)

    prompt = f"""당신은 고체전지 복합양극 DEM 시뮬레이션 전문가입니다.
아래는 bimodal AM (AM_P: 대립자 6μm, AM_S: 소립자 2μm) + SE (고체전해질) 복합 양극의 DEM 분석 비교 데이터입니다.

제목: {title}
{f'메모: {notes}' if notes else ''}

## 데이터

{data_table}

## 분석 원칙 (반드시 준수)

1. **자의적 가중치 사용 금지**: "종합 점수"를 만들 때 근거 없는 가중치(예: 0.4/0.6)를 부여하지 마세요.
2. **물리적 근거 기반 판단**: 각 주장에는 반드시 물리적 메커니즘 설명이 필요합니다.
3. **Percolation 포화 효과**: 95% 이상에서는 diminishing return. 93%→97% 차이보다 tortuosity 차이가 실제 성능에 더 큰 영향.
4. **병목 결정 원칙**: 고체전지에서는 일반적으로 ionic transport이 charge transfer보다 훨씬 큰 저항. AM-SE 계면이 "충분"하면 SE 네트워크 효율(tortuosity↓)이 핵심.
5. **데이터에 없는 값을 추정하지 마세요**: 주어진 수치만으로 분석.
6. **핵심 지표 우선순위**: Tortuosity > SE-SE Total Area > Percolation > Porosity > AM-SE Total (ionic transport limited 시스템 기준)

## 분석 내용 (한국어)

### 1. 핵심 발견
- AM_P 비율 변화에 따른 주요 경향 3-5개
- 각 경향의 물리적 원인 (입자 크기, 공간 배치 관점)

### 2. SE Contact Network Trade-off
- SE-SE 접촉 개수 vs 평균 면적의 trade-off 관계
- Total Area = N × Mean → 최적점 분석
- 직관적 비유 활용

### 3. AM-SE vs SE-SE 상위 Trade-off
- AM-SE 계면 (charge transfer) vs SE 네트워크 (ionic transport)
- 어느 쪽이 실제 병목인지 데이터 기반으로 판단
- 최적 P:S 비율 제안 (근거 명시)

### 4. 결론 및 제언
- 종합 최적 조건과 그 이유
- 주의사항, 추가 검증 필요 사항

마크다운 형식으로, 표/수치를 적극 활용하되 근거 없는 수치 생성은 금지."""

    try:
        message = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=3000,
            messages=[{"role": "user", "content": prompt}]
        )
        return message.content[0].text
    except Exception as e:
        return f"*AI 분석 생성 실패: {str(e)}*"


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
        if mode == 'bimodal':
            type_map = '1:AM_P,2:AM_S,3:SE'
        else:
            # Standard: detect AM_P vs AM_S from radius in atom file
            am_type_name = 'AM_S'  # default: small AM
            for f in sorted(os.listdir(case_dir)):
                if f.startswith('atom') and f.endswith('.liggghts'):
                    with open(os.path.join(case_dir, f)) as fh:
                        in_data = False
                        for line in fh:
                            stripped = line.strip()
                            if stripped.startswith('ITEM: ATOMS'):
                                in_data = True
                                continue
                            if stripped.startswith('ITEM:'):
                                in_data = False
                                continue
                            if in_data:
                                parts = stripped.split()
                                if len(parts) >= 6:
                                    try:
                                        t = int(parts[1])
                                        r = float(parts[5])
                                        if t == 1:
                                            # sim r > 0.004 → AM_P (6μm), else AM_S (2μm)
                                            am_type_name = 'AM_P' if r > 0.004 else 'AM_S'
                                            break
                                    except (ValueError, IndexError):
                                        continue
                    break
            type_map = f'1:{am_type_name},2:SE'

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

    # Sync to Supabase
    storage_sync.sync_dir_to_remote(case_dir, f'uploads/{case_id}')

    return jsonify({'case_id': case_id, 'mode': mode, 'files': filenames})

@app.route('/analyze/<case_id>', methods=['POST'])
def analyze(case_id):
    """Run analysis pipeline for a case (background thread)."""
    case_dir = get_case_dir(case_id)
    meta_file = os.path.join(case_dir, 'meta.json')
    if not os.path.exists(meta_file):
        return jsonify({'error': '케이스를 찾을 수 없습니다.'}), 404

    with open(meta_file) as f:
        meta = json.load(f)

    meta['status'] = 'running'
    with open(meta_file, 'w') as f:
        json.dump(meta, f, indent=2)

    def _run():
        # Clear previous results for clean re-analysis
        results_dir = get_results_dir(case_id)
        if os.path.exists(results_dir):
            shutil.rmtree(results_dir)
        os.makedirs(results_dir, exist_ok=True)

        result = run_pipeline(case_id, meta['mode'], meta['type_map'], meta.get('scale', 1000))

        meta['status'] = 'done' if result.get('success') else 'error'
        meta['analysis_log'] = result.get('log', [])
        with open(meta_file, 'w') as f:
            json.dump(meta, f, indent=2)

        if result.get('success'):
            generate_report(case_id, meta.get('name', ''))

        # Sync results + updated meta to Supabase
        storage_sync.sync_dir_to_remote(case_dir, f'uploads/{case_id}')
        storage_sync.sync_dir_to_remote(results_dir, f'results/{case_id}')

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    return jsonify({'success': True, 'status': 'running'})

@app.route('/analyze-status/<case_id>')
def analyze_status(case_id):
    """Check if analysis is still running."""
    meta_file = os.path.join(get_case_dir(case_id), 'meta.json')
    if os.path.exists(meta_file):
        with open(meta_file) as f:
            meta = json.load(f)
        return jsonify({'status': meta.get('status', 'unknown')})
    return jsonify({'status': 'unknown'})

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

    # Inject section headers into network_summary if missing
    if 'network_summary' in tables:
        data = tables['network_summary']['data']
        has_headers = any(str(row[0]).startswith('──') for row in data)
        if not has_headers:
            section_map = {
                'Porosity(%)': '── 구조 ──',
                'AM-SE Total(μm²)': '── 계면 ──',
                'SE-SE CN mean': '── 이온경로: 연결성 ──',
                'Tortuosity mean': '── 이온경로: 경로 효율 ──',
                'Path Hop Area mean(μm²)': '── 이온경로: 경로 품질 ──',
                'Ionic Active AM(%)': '── 활성도 ──',
                'Stress CV(%)': '── 응력 ──',
            }
            new_data = []
            for row in data:
                label = str(row[0])
                if label in section_map:
                    new_data.append([section_map[label], ''])
                new_data.append(row)
            tables['network_summary']['data'] = new_data

    # Load full_metrics.json for header info
    metrics = {}
    metrics_path = os.path.join(results_dir, 'full_metrics.json')
    if os.path.exists(metrics_path):
        with open(metrics_path) as f:
            metrics = json.load(f)

    # Patch network_summary with values from full_metrics.json
    if 'network_summary' in tables and metrics:
        # SE Cluster: plain number → large/total format
        n_large = metrics.get('n_large_components')
        if n_large is not None:
            for row in tables['network_summary']['data']:
                if str(row[0]) == 'SE Cluster 수' and '≥10' not in str(row[1]):
                    row[1] = f"{n_large}(≥10) / {row[1]}"
        # Fill placeholder '-' values from full_metrics
        placeholder_map = {
            'GB Density(hops/μm)': 'gb_density_mean',
            'Path Hop Area mean(μm²)': 'path_hop_area_mean',
            'Path Bottleneck(μm²)': 'path_hop_area_min_mean',
            'Path Conductance(μm²)': 'path_conductance_mean',
        }
        for row in tables['network_summary']['data']:
            label = str(row[0])
            if label in placeholder_map and str(row[1]).strip() in ('-', ''):
                val = metrics.get(placeholder_map[label])
                if val is not None:
                    row[1] = val

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
    case_groups_param = request.args.get('case_groups', '[]')

    comparison_data = {}
    case_groups_parsed = []
    try:
        case_groups_parsed = json.loads(case_groups_param)
    except:
        pass
    if selected:
        # Key metrics grouped by category
        # (label, unit, key, category)
        display_keys_raw = [
            # ── 구조/계면 ──
            ('P:S', '', 'ps_ratio', '구조/계면'),
            ('Porosity', '(%)', 'porosity', '구조/계면'),
            ('두께', '(μm)', 'thickness_um', '구조/계면'),
            ('AM-SE Total', '(μm²)', 'area_AM전체_SE_total', '구조/계면'),
            ('SE-SE N', '', 'area_SE_SE_n', '구조/계면'),
            ('SE-SE Mean', '(μm²)', 'area_SE_SE_mean', '구조/계면'),
            ('SE-SE Total', '(μm²)', 'area_SE_SE_total', '구조/계면'),
            ('Coverage P', '(%)', 'coverage_AM_P_mean', '구조/계면'),
            ('Coverage S', '(%)', 'coverage_AM_S_mean', '구조/계면'),
            # ── 이온경로 ──
            ('SE-SE CN', '', 'se_se_cn', '이온경로'),
            ('SE Cluster', '', 'n_components', '이온경로'),
            ('Large(≥10)', '', 'n_large_components', '이온경로'),
            ('Percolation', '(%)', 'percolation_pct', '이온경로'),
            ('Top Reachable', '(%)', 'top_reachable_pct', '이온경로'),
            ('Tortuosity', '', 'tortuosity_mean', '이온경로'),
            ('τ std', '', 'tortuosity_std', '이온경로'),
            ('GB Density', '(hops/μm)', 'gb_density_mean', '이온경로'),
            ('Hop Area', '(μm²)', 'path_hop_area_mean', '이온경로'),
            ('Bottleneck', '(μm²)', 'path_hop_area_min_mean', '이온경로'),
            ('Path Conductance', '(μm²)', 'path_conductance_mean', '이온경로'),
            # ── 활성도/전도 ──
            ('Ionic Active', '(%)', 'ionic_active_pct', '활성도'),
            ('AM-SE CN', '', 'am_se_cn_mean', '활성도'),
            ('Vulnerable', '(%)', 'am_vulnerable_pct', '활성도'),
            ('φ_SE', '', 'phi_se', '활성도'),
            ('σ_eff/σ_bulk', '', 'sigma_ratio', '활성도'),
            # ── 접촉력/응력 ──
            ('Fn AM-AM', '(μN)', 'fn_AM_P_AM_P_mean', '접촉력/응력'),
            ('Fn AM-SE', '(μN)', 'fn_AM_P_SE_mean', '접촉력/응력'),
            ('Fn SE-SE', '(μN)', 'fn_SE_SE_mean', '접촉력/응력'),
            ('CP mean', '(MPa)', 'contact_pressure_mean', '접촉력/응력'),
            ('CP max', '(MPa)', 'contact_pressure_max', '접촉력/응력'),
            ('Stress CV', '(%)', 'stress_cv', '접촉력/응력'),
            ('σ_AM_P/σ_mean', '', 'stress_ratio_AM_P', '접촉력/응력'),
            ('σ_AM_S/σ_mean', '', 'stress_ratio_AM_S', '접촉력/응력'),
            ('σ_SE/σ_mean', '', 'stress_ratio_SE', '접촉력/응력'),
        ]
        display_keys = [(l, u, k) for l, u, k, _ in display_keys_raw]
        # Track category boundaries for column separators
        col_categories = [''] + [cat for _, _, _, cat in display_keys_raw]
        rows = []
        for cid in selected:
            # Handle archive: prefix
            if cid.startswith('archive:'):
                archive_rel = cid[len('archive:'):]
                case_path = os.path.join(app.config['ARCHIVE_FOLDER'], archive_rel)
                meta_file = os.path.join(case_path, 'meta.json')
                metrics_path = os.path.join(case_path, 'full_metrics.json')
                case_name = os.path.basename(archive_rel)
            else:
                case_path = get_results_dir(cid)
                meta_file = os.path.join(get_case_dir(cid), 'meta.json')
                metrics_path = os.path.join(case_path, 'full_metrics.json')
                case_name = cid

            if os.path.exists(meta_file):
                with open(meta_file) as f:
                    meta = json.load(f)
                case_name = meta.get('name', case_name)

            metrics = {}
            if os.path.exists(metrics_path):
                with open(metrics_path) as f:
                    metrics = json.load(f)
            else:
                continue

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

            # Force metric fallbacks: AM_P↔AM_S
            if 'fn_AM_P_AM_P_mean' not in metrics and 'fn_AM_S_AM_S_mean' in metrics:
                metrics['fn_AM_P_AM_P_mean'] = metrics['fn_AM_S_AM_S_mean']
            if 'fn_AM_P_SE_mean' not in metrics and 'fn_AM_S_SE_mean' in metrics:
                metrics['fn_AM_P_SE_mean'] = metrics['fn_AM_S_SE_mean']

            row = {'케이스': case_name}
            for label, unit, key in display_keys:
                val = metrics.get(key, '')
                if isinstance(val, float):
                    if key in ('path_conductance_mean', 'path_hop_area_min_mean'):
                        val = f"{val:.2e}" if val > 0 else '-'
                    else:
                        val = round(val, 2)
                row[label] = val if val != '' else '-'
            rows.append(row)

        if rows:
            # Build columns with unit subtitles
            col_headers = [{'name': '케이스', 'unit': ''}]
            for label, unit, key in display_keys:
                col_headers.append({'name': label, 'unit': unit})

            # Split into case group tables
            if case_groups_parsed and len(case_groups_parsed) > 1:
                tables = []
                row_idx = 0
                for gi, g in enumerate(case_groups_parsed):
                    gname = g.get('name', '') or f"Case {chr(65+gi)}"
                    group_rows = []
                    for _ in g.get('cases', []):
                        if row_idx < len(rows):
                            group_rows.append(rows[row_idx])
                            row_idx += 1
                    # Mark best values per column
                    # lower_better: Porosity, 두께, Tortuosity, τ std, Stress CV, GB Density, Vulnerable, CP mean, CP max
                    # higher_better: everything else (except P:S, 케이스 which are labels)
                    lower_better = {'Porosity', '두께', 'Tortuosity', 'τ std', 'Stress CV',
                                    'GB Density', 'Vulnerable', 'CP mean', 'CP max', 'SE Cluster'}
                    skip_cols = {'P:S', '케이스'}
                    best_marks = {}  # col -> best row index
                    for label, unit, key in display_keys:
                        if label in skip_cols:
                            continue
                        vals = []
                        for ri, r in enumerate(group_rows):
                            v = r.get(label, '-')
                            try:
                                vals.append((ri, float(str(v).replace('e', 'E').strip())))
                            except (ValueError, TypeError):
                                pass
                        if vals:
                            if label in lower_better:
                                best_val = min(vals, key=lambda x: x[1])
                                worst_val = max(vals, key=lambda x: x[1])
                            else:
                                best_val = max(vals, key=lambda x: x[1])
                                worst_val = min(vals, key=lambda x: x[1])
                            # Skip if all same value
                            if best_val[1] != worst_val[1]:
                                best_marks[(label, best_val[0])] = True
                    for ri, r in enumerate(group_rows):
                        r['__best__'] = {label for (label, idx), _ in best_marks.items() if idx == ri}

                    tables.append({'name': gname, 'rows': group_rows, 'color': ['#6c8cff','#ff6b6b','#51cf66','#ffd43b'][gi % 4]})
                comparison_data = {
                    'columns': [c['name'] for c in col_headers],
                    'units': [c['unit'] for c in col_headers],
                    'categories': col_categories,
                    'tables': tables
                }
            else:
                comparison_data = {
                    'columns': [c['name'] for c in col_headers],
                    'units': [c['unit'] for c in col_headers],
                    'categories': col_categories,
                    'tables': [{'name': '', 'rows': rows, 'color': ''}]
                }

    # Scan archive for folders with full_metrics.json
    archive_folders = []
    archive_root = app.config['ARCHIVE_FOLDER']
    if os.path.isdir(archive_root):
        for dirpath, dirnames, filenames in os.walk(archive_root):
            # Count cases in this folder (subfolders with full_metrics.json)
            case_count = 0
            for d in dirnames:
                if os.path.exists(os.path.join(dirpath, d, 'full_metrics.json')):
                    case_count += 1
            if case_count > 0:
                rel = os.path.relpath(dirpath, archive_root)
                archive_folders.append({'path': rel if rel != '.' else '(최상위)', 'case_count': case_count})

    return render_template('group.html', cases=cases, selected=selected,
                         comparison=comparison_data, archive_folders=archive_folders,
                         case_groups_json=case_groups_param)

@app.route('/group/archive-cases')
def group_archive_cases():
    """Return cases in an archive folder that have full_metrics.json."""
    folder = request.args.get('folder', '')
    archive_root = app.config['ARCHIVE_FOLDER']
    if folder == '(최상위)':
        base = archive_root
    else:
        base = os.path.join(archive_root, folder)
    if not os.path.isdir(base):
        return jsonify({'cases': []})

    cases = []
    for name in sorted(os.listdir(base)):
        case_dir = os.path.join(base, name)
        metrics_path = os.path.join(case_dir, 'full_metrics.json')
        if os.path.isdir(case_dir) and os.path.exists(metrics_path):
            with open(metrics_path) as f:
                m = json.load(f)
            # Use archive: prefix to distinguish from dashboard cases
            case_id = f"archive:{os.path.relpath(case_dir, archive_root)}"
            cases.append({
                'id': case_id,
                'name': name,
                'ps_ratio': m.get('ps_ratio', ''),
                'warning_count': m.get('warning_count', 0),
                'warning_msgs': [w['msg'] for w in m.get('warnings', [])],
            })
    return jsonify({'cases': cases})

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
        if cid.startswith('archive:'):
            archive_rel = cid[len('archive:'):]
            case_path = os.path.join(app.config['ARCHIVE_FOLDER'], archive_rel)
            metrics_path = os.path.join(case_path, 'full_metrics.json')
            case_name = os.path.basename(archive_rel)
        else:
            case_path = get_results_dir(cid)
            metrics_path = os.path.join(case_path, 'full_metrics.json')
            meta_file = os.path.join(get_case_dir(cid), 'meta.json')
            case_name = cid
            if os.path.exists(meta_file):
                with open(meta_file) as f:
                    case_name = json.load(f).get('name', cid)
        if os.path.exists(metrics_path):
            input_files.append(metrics_path)
            names.append(case_name)

    if not input_files:
        return jsonify({'error': '메트릭 데이터가 없습니다.'}), 400

    # Pass case groups info
    case_groups_json = request.form.get('case_groups', '[]')
    try:
        case_groups_raw = json.loads(case_groups_json)
    except:
        case_groups_raw = []

    # Build group sizes string: "3,5" means first 3 cases = group A, next 5 = group B
    # Build group names string
    group_sizes = []
    group_names_list = []
    for g in case_groups_raw:
        group_sizes.append(str(len(g.get('cases', []))))
        group_names_list.append(g.get('name', '') or f"Case {chr(65 + len(group_names_list))}")

    global_rgb = request.form.get('global_rgb', '')

    scripts = app.config['SCRIPTS_FOLDER']
    cmd = ['python3', os.path.join(scripts, 'generate_comparison_plots.py'),
           '-i'] + input_files + ['-n'] + names + ['-o', plot_dir, '-p'] + plots
    if group_sizes and len(group_sizes) > 1:
        cmd += ['--group-sizes', ','.join(group_sizes),
                '--group-names', ','.join(group_names_list)]
    if global_rgb:
        cmd += ['--global-rgb', global_rgb]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    if result.returncode != 0:
        return jsonify({'error': f'Plot 생성 실패: {result.stderr}'}), 500

    # Load plot info
    info_path = os.path.join(plot_dir, 'plot_info.json')
    plot_list = []
    if os.path.exists(info_path):
        with open(info_path) as f:
            info = json.load(f)
        # particle_info always first
        if 'particle_info' in info:
            plot_list.append(info['particle_info'])
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
    """Generate comprehensive group comparison markdown report."""
    selected = request.form.getlist('cases')
    title = request.form.get('title', 'DEM_Bimodal_Comparison')
    notes = request.form.get('notes', '')

    now = datetime.now().strftime('%y%m%d')
    L = []  # lines

    # Header
    L.append(f'# {now}_{title}\n')
    L.append(f'> **날짜**: {datetime.now().strftime("%Y-%m-%d")}')
    L.append(f'> **비교 케이스**: {len(selected)}개')
    if notes:
        L.append(f'> **메모**: {notes}')
    L.append('')
    L.append('---\n')

    # Load all metrics
    all_metrics = []
    case_names = []
    for cid in selected:
        if cid.startswith('archive:'):
            archive_rel = cid[len('archive:'):]
            case_path = os.path.join(app.config['ARCHIVE_FOLDER'], archive_rel)
            meta_file = os.path.join(case_path, 'meta.json')
            metrics_path = os.path.join(case_path, 'full_metrics.json')
            name = os.path.basename(archive_rel)
        else:
            case_path = get_results_dir(cid)
            meta_file = os.path.join(get_case_dir(cid), 'meta.json')
            metrics_path = os.path.join(case_path, 'full_metrics.json')
            name = cid
        if os.path.exists(meta_file):
            with open(meta_file) as f:
                name = json.load(f).get('name', name)
        metrics = {}
        if os.path.exists(metrics_path):
            with open(metrics_path) as f:
                metrics = json.load(f)
        all_metrics.append(metrics)
        case_names.append(name)

    if not all_metrics:
        return jsonify({'report': '메트릭 데이터가 없습니다.', 'path': ''})

    # 1. System Overview
    L.append('## 1. System Overview\n')
    import pandas as pd
    overview_rows = []
    display_keys = [
        ('P:S', 'ps_ratio'), ('Porosity(%)', 'porosity'),
        ('Thickness(μm)', 'thickness_um'),
        ('AM-SE Total(μm²)', 'area_AM전체_SE_total'),
        ('SE-SE Total(μm²)', 'area_SE_SE_total'),
        ('SE-SE CN', 'se_se_cn'), ('SE Cluster', 'n_components'),
        ('Percolation(%)', 'percolation_pct'),
        ('Top Reachable(%)', 'top_reachable_pct'),
        ('Tortuosity', 'tortuosity_mean'),
        ('Ionic Active(%)', 'ionic_active_pct'),
    ]
    for i, name in enumerate(case_names):
        row = {'Case': name}
        for label, key in display_keys:
            val = all_metrics[i].get(key, '-')
            if isinstance(val, float):
                val = round(val, 2)
            row[label] = val
        overview_rows.append(row)
    df = pd.DataFrame(overview_rows)
    L.append(df.to_markdown(index=False))
    L.append('')

    # 2. Key Findings (자동 분석)
    L.append('\n## 2. Key Findings\n')

    porosities = [(case_names[i], m.get('porosity', 0)) for i, m in enumerate(all_metrics) if m.get('porosity')]
    if porosities:
        min_p = min(porosities, key=lambda x: x[1])
        L.append(f'- **Porosity 최저**: {min_p[0]} ({min_p[1]:.2f}%)')

    percs = [(case_names[i], m.get('percolation_pct', 0)) for i, m in enumerate(all_metrics) if m.get('percolation_pct')]
    if percs:
        max_perc = max(percs, key=lambda x: x[1])
        L.append(f'- **Percolation 최대**: {max_perc[0]} ({max_perc[1]:.1f}%)')

    torts = [(case_names[i], m.get('tortuosity_mean', 99)) for i, m in enumerate(all_metrics) if m.get('tortuosity_mean')]
    if torts:
        min_tort = min(torts, key=lambda x: x[1])
        L.append(f'- **Tortuosity 최저**: {min_tort[0]} ({min_tort[1]:.2f})')

    ionics = [(case_names[i], m.get('ionic_active_pct', 0)) for i, m in enumerate(all_metrics) if m.get('ionic_active_pct')]
    if ionics:
        max_ionic = max(ionics, key=lambda x: x[1])
        L.append(f'- **Ionic Active 최대**: {max_ionic[0]} ({max_ionic[1]:.1f}%)')

    se_totals = [(case_names[i], m.get('area_SE_SE_total', 0)) for i, m in enumerate(all_metrics) if m.get('area_SE_SE_total')]
    if se_totals:
        max_se = max(se_totals, key=lambda x: x[1])
        L.append(f'- **SE-SE Total Area 최대**: {max_se[0]} ({max_se[1]:,.1f} μm²)')

    L.append('')

    # 3. Trade-off Analysis
    L.append('\n## 3. Trade-off Analysis\n')
    L.append('### AM-SE vs SE-SE Trade-off\n')
    L.append('| AM_P ↑ | AM-SE 계면 | SE-SE 네트워크 |')
    L.append('|---|---|---|')
    L.append('| 변화 | **감소** | **개선** (percolation↑, tortuosity↓) |')
    L.append('| 의미 | 반응 면적 감소 | 이온 경로 확보 |')
    L.append('| 제한 요인 | charge transfer | ionic transport |')
    L.append('')
    L.append('> **전극 성능 = 병목(bottleneck)이 결정**')
    L.append('> - SE 네트워크 부족 (P 적음) → ionic transport limited')
    L.append('> - AM-SE 계면 부족 (P 많음) → charge transfer limited')
    L.append('')

    # 4. SE-SE Contact Network Trade-off
    L.append('### SE-SE Contact: Quality vs Quantity\n')
    L.append('AM_P↑ → SE-SE 접촉 개수↑ (넓은 공간에 분산) + 개별 접촉 면적↓ (느슨하게 배치)')
    L.append('')
    se_n = [(case_names[i], m.get('area_SE_SE_n', 0), m.get('area_SE_SE_mean', 0))
            for i, m in enumerate(all_metrics) if m.get('area_SE_SE_n')]
    if se_n:
        L.append('| Case | SE-SE N | SE-SE Mean Area(μm²) | SE-SE Total(μm²) |')
        L.append('|---|---|---|---|')
        for name, n, mean in se_n:
            total = n * mean
            L.append(f'| {name} | {int(n):,} | {mean:.4f} | {total:,.1f} |')
    L.append('')

    # 5. Conclusion
    L.append('\n## 4. Conclusion\n')
    if porosities and len(porosities) > 2:
        L.append(f'- Porosity: {porosities[0][0]} ({porosities[0][1]:.1f}%) → {porosities[-1][0]} ({porosities[-1][1]:.1f}%)')
    if percs and torts:
        L.append(f'- 이온경로: Percolation {percs[0][1]:.0f}%→{percs[-1][1]:.0f}%, Tortuosity {torts[0][1]:.1f}→{torts[-1][1]:.1f}')
    L.append('')

    # Tags
    L.append('---\n')

    # Claude AI 심층 분석
    ai_analysis = _generate_ai_analysis(all_metrics, case_names, title, notes)
    if ai_analysis:
        L.append('\n## 5. AI 심층 분석 (Claude)\n')
        L.append(ai_analysis)
        L.append('')

    L.append('---\n')
    L.append('#DEM #bimodal #comparison #percolation #tortuosity #coverage #ASSB\n')

    report = '\n'.join(L)

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

@app.route('/results/<case_id>/3d-data')
def serve_3d_data(case_id):
    """Serve particle + percolation data for 3D viewer."""
    import pandas as pd
    results_dir = get_results_dir(case_id)
    case_dir = get_case_dir(case_id)
    meta_file = os.path.join(case_dir, 'meta.json')

    if not os.path.exists(meta_file):
        return jsonify({'error': 'Case not found'}), 404

    with open(meta_file) as f:
        meta = json.load(f)
    scale = meta.get('scale', 1000)
    type_map_str = meta.get('type_map', '1:AM,2:SE')
    type_map = {}
    for item in type_map_str.split(','):
        k, v = item.split(':')
        type_map[int(k)] = v.strip()

    # Load atoms
    atoms_csv = os.path.join(results_dir, 'atoms.csv')
    if not os.path.exists(atoms_csv):
        return jsonify({'error': 'No atom data'}), 404

    df = pd.read_csv(atoms_csv)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    particles = []
    for _, row in df.iterrows():
        t = int(row['type'])
        particles.append({
            'id': int(row['id']),
            'type': type_map.get(t, f'T{t}'),
            'x': round(float(row['x']) * scale, 2),
            'y': round(float(row['y']) * scale, 2),
            'z': round(float(row['z']) * scale, 2),
            'r': round(float(row['radius']) * scale, 2),
        })

    # Box bounds
    box = {
        'x_min': 0, 'x_max': round(0.05 * scale, 1),
        'y_min': 0, 'y_max': round(0.05 * scale, 1),
        'z_min': 0,
    }
    # plate_z from mesh_info
    mesh_file = os.path.join(results_dir, 'mesh_info.json')
    if os.path.exists(mesh_file):
        with open(mesh_file) as f:
            box['z_max'] = round(json.load(f)['plate_z'] * scale, 1)
    else:
        box['z_max'] = round(max(p['z'] + p['r'] for p in particles), 1)

    # Percolation data from full_metrics
    percolation = {'top_reachable': [], 'bottom_se': [], 'top_se': []}
    metrics_path = os.path.join(results_dir, 'full_metrics.json')

    # Try to load percolation sets from a saved file
    perc_path = os.path.join(results_dir, 'percolation_sets.json')
    if os.path.exists(perc_path):
        with open(perc_path) as f:
            percolation = json.load(f)

    # Paths (tortuosity sample paths)
    paths = []
    paths_path = os.path.join(results_dir, 'tortuosity_paths.json')
    if os.path.exists(paths_path):
        with open(paths_path) as f:
            paths = json.load(f)

    # SE clusters for click interaction
    clusters = {}
    clusters_path = os.path.join(results_dir, 'se_clusters.json')
    if os.path.exists(clusters_path):
        with open(clusters_path) as f:
            clusters = json.load(f)

    return jsonify({
        'particles': particles,
        'box': box,
        'percolation': percolation,
        'paths': paths,
        'clusters': clusters,
    })

@app.route('/toggle-warning/<case_id>', methods=['POST'])
def toggle_warning(case_id):
    """Toggle a warning on/off in full_metrics.json."""
    data = request.get_json()
    warn_type = data.get('warning_type', '')
    if not warn_type:
        return jsonify({'error': 'No warning type'}), 400

    # Try dashboard results first, then archive
    metrics_path = os.path.join(get_results_dir(case_id), 'full_metrics.json')
    if case_id.startswith('archive:'):
        archive_rel = case_id[len('archive:'):]
        target = _safe_path(archive_rel)
        if target:
            metrics_path = os.path.join(target, 'full_metrics.json')

    if not os.path.exists(metrics_path):
        return jsonify({'error': 'Metrics not found'}), 404

    with open(metrics_path) as f:
        metrics = json.load(f)

    disabled = metrics.get('disabled_warnings', [])
    if warn_type in disabled:
        disabled.remove(warn_type)
        is_disabled = False
    else:
        disabled.append(warn_type)
        is_disabled = True
    metrics['disabled_warnings'] = disabled

    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2, default=str)

    return jsonify({'success': True, 'disabled': is_disabled})


@app.route('/results/<case_id>/save-screenshot', methods=['POST'])
def save_screenshot(case_id):
    """Save 3D viewer screenshot to figures folder."""
    import base64
    data = request.get_json()
    if not data or 'image' not in data:
        return jsonify({'error': 'No image data'}), 400
    figures_dir = os.path.join(get_results_dir(case_id), 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    filename = data.get('filename', 'screenshot.png')
    filename = filename.replace('/', '_').replace('\\', '_')
    img_data = data['image'].split(',')[1]  # strip data:image/png;base64,
    with open(os.path.join(figures_dir, filename), 'wb') as f:
        f.write(base64.b64decode(img_data))
    return jsonify({'success': True, 'filename': filename})


@app.route('/results/<case_id>/force-chains')
def serve_force_chains(case_id):
    """Serve force chain data for 3D viewer."""
    results_dir = get_results_dir(case_id)
    fc_path = os.path.join(results_dir, 'force_chains.json')
    if os.path.exists(fc_path):
        with open(fc_path) as f:
            return jsonify(json.load(f))
    return jsonify([])

@app.route('/results/<case_id>/report')
def serve_report(case_id):
    """Generate MD report on-the-fly from analysis CSVs."""
    import pandas as pd
    results_dir = get_results_dir(case_id)
    case_dir = get_case_dir(case_id)
    meta_file = os.path.join(case_dir, 'meta.json')

    meta = {}
    if os.path.exists(meta_file):
        with open(meta_file) as f:
            meta = json.load(f)

    metrics = {}
    metrics_path = os.path.join(results_dir, 'full_metrics.json')
    if os.path.exists(metrics_path):
        with open(metrics_path) as f:
            metrics = json.load(f)

    input_params = {}
    params_path = os.path.join(results_dir, 'input_params.json')
    if os.path.exists(params_path):
        with open(params_path) as f:
            input_params = json.load(f)

    name = meta.get('name', case_id)
    now = datetime.now().strftime('%y%m%d')
    L = []

    # Header
    L.append(f'# {now}_{name}_DEM_Analysis\n')
    L.append(f'> **날짜**: {datetime.now().strftime("%Y-%m-%d")}')
    L.append(f'> **케이스**: {name}')
    L.append(f'> **모드**: {meta.get("mode", "-")}')
    if metrics.get('ps_ratio'):
        L.append(f'> **P:S**: {metrics["ps_ratio"]}')
    if input_params.get('am_se_ratio'):
        L.append(f'> **AM:SE**: {input_params["am_se_ratio"]}')
    if input_params.get('target_press_sim'):
        L.append(f'> **Target Pressure**: {input_params["target_press_sim"] * 1000:.1f} MPa')
    L.append('')
    L.append('---\n')

    # Load and append each CSV as table
    csv_names = {
        'atom_statistics': '입자 정보',
        'contact_summary': '접촉 요약',
        'coordination_summary': '배위수',
        'network_summary': '네트워크 지표',
    }
    section_num = 1
    for csv_name, title in csv_names.items():
        csv_path = os.path.join(results_dir, f'{csv_name}.csv')
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            L.append(f'## {section_num}. {title}\n')
            L.append(df.to_markdown(index=False))
            L.append('')
            section_num += 1

    # Key metrics summary
    if metrics:
        L.append(f'\n## {section_num}. 핵심 지표 요약\n')
        summary_items = [
            ('Porosity', f"{metrics.get('porosity', '-')}%"),
            ('전극 두께', f"{metrics.get('thickness_um', '-')} μm"),
            ('SE-SE CN', f"{metrics.get('se_se_cn', '-')}"),
            ('SE Percolation', f"{metrics.get('percolation_pct', '-')}%"),
            ('Top Reachable', f"{metrics.get('top_reachable_pct', '-')}%"),
            ('Tortuosity', f"{metrics.get('tortuosity_mean', '-')}"),
            ('Ionic Active AM', f"{metrics.get('ionic_active_pct', '-')}%"),
        ]
        if metrics.get('stress_cv'):
            summary_items.append(('Stress CV', f"{metrics.get('stress_cv', '-')}%"))
        for label, val in summary_items:
            L.append(f'- **{label}**: {val}')
        L.append('')

    L.append('---\n')
    L.append('#DEM #analysis #ASSB #composite-cathode\n')

    report = '\n'.join(L)

    # Return as downloadable MD
    from io import BytesIO
    buf = BytesIO(report.encode('utf-8'))
    return send_file(buf, mimetype='text/markdown', as_attachment=True,
                    download_name=f'{name}_report.md')

@app.route('/download-file/<case_id>/<filename>')
def download_case_file(case_id, filename):
    """Download an uploaded file from a case."""
    case_dir = get_case_dir(case_id)
    return send_from_directory(case_dir, filename, as_attachment=True)

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
    storage_sync.upload_file(f'uploads/{case_id}/meta.json', meta_file)
    return jsonify({'success': True})

@app.route('/delete/<case_id>', methods=['POST'])
def delete_case(case_id):
    case_dir = get_case_dir(case_id)
    results_dir = get_results_dir(case_id)
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    if os.path.exists(results_dir):
        shutil.rmtree(results_dir)
    storage_sync.delete_prefix(f'uploads/{case_id}')
    storage_sync.delete_prefix(f'results/{case_id}')
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
            has_metrics = os.path.exists(os.path.join(full, 'full_metrics.json'))
            folders.append({'name': name, 'path': rel, 'file_count': cnt, 'subfolder_count': sub,
                           'has_metrics': has_metrics})
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
    rel = os.path.relpath(target, app.config['ARCHIVE_FOLDER'])
    if os.path.isdir(target):
        storage_sync.delete_prefix(f'archive/{rel}')
        shutil.rmtree(target)
    elif os.path.isfile(target):
        storage_sync.delete_path(f'archive/{rel}')
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
    # Also copy original uploaded files (atom, contact, mesh, input .liggghts)
    raw_dir = os.path.join(save_dir, 'raw_files')
    os.makedirs(raw_dir, exist_ok=True)
    for fname in os.listdir(case_dir):
        fpath = os.path.join(case_dir, fname)
        if os.path.isfile(fpath) and fname != 'meta.json':
            shutil.copy2(fpath, os.path.join(raw_dir, fname))
    # Sync archive to Supabase
    rel = os.path.relpath(save_dir, app.config['ARCHIVE_FOLDER'])
    storage_sync.sync_dir_to_remote(save_dir, f'archive/{rel}')
    return jsonify({'success': True, 'saved_to': save_dir})

@app.route('/archive/reanalyze/<path:folder>', methods=['POST'])
def archive_reanalyze(folder):
    """Re-run analysis on an archive case using raw_files."""
    target = _safe_path(folder)
    if not target or not os.path.isdir(target):
        return jsonify({'error': 'Not found'}), 404

    meta_file = os.path.join(target, 'meta.json')
    if not os.path.exists(meta_file):
        return jsonify({'error': 'No meta.json'}), 404

    with open(meta_file) as f:
        meta = json.load(f)

    # Find source files: raw_files/ or directly in folder
    raw_dir = os.path.join(target, 'raw_files')
    source_dir = raw_dir if os.path.isdir(raw_dir) else target

    atom_files = sorted(globmod.glob(os.path.join(source_dir, 'atom*.liggghts')))
    contact_files = sorted(globmod.glob(os.path.join(source_dir, 'contact*.liggghts')))
    mesh_files = sorted(globmod.glob(os.path.join(source_dir, '*.stl')))
    input_files = sorted(globmod.glob(os.path.join(source_dir, 'input*.liggghts')))

    if not atom_files or not contact_files:
        return jsonify({'error': 'No atom/contact files in raw_files/'}), 400

    # Write a status file for polling
    status_file = os.path.join(target, '.reanalyze_status')
    with open(status_file, 'w') as f:
        f.write('running')

    def _run():
        scripts = app.config['SCRIPTS_FOLDER']
        mode = meta.get('mode', 'standard')
        type_map = meta.get('type_map', '1:AM,2:SE')
        scale = meta.get('scale', 1000)

        for pyc in globmod.glob(os.path.join(scripts, '__pycache__', '*.pyc')):
            os.remove(pyc)

        cmd = ['python3', os.path.join(scripts, 'parse_liggghts.py')]
        cmd += atom_files + contact_files + mesh_files + input_files + ['-o', target]
        subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        atoms_csv = os.path.join(target, 'atoms.csv')
        contacts_csv = os.path.join(target, 'contacts.csv')
        script = 'analyze_contacts_bimodal.py' if mode == 'bimodal' else 'analyze_contacts.py'
        cmd = ['python3', os.path.join(scripts, script),
               atoms_csv, contacts_csv, '-o', target,
               '-t', type_map, '-s', str(scale)]
        subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        with open(status_file, 'w') as f:
            f.write('done')

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    return jsonify({'success': True, 'status': 'running'})


@app.route('/archive/reanalyze-status/<path:folder>')
def archive_reanalyze_status(folder):
    target = _safe_path(folder)
    if not target:
        return jsonify({'status': 'unknown'})
    status_file = os.path.join(target, '.reanalyze_status')
    if os.path.exists(status_file):
        with open(status_file) as f:
            return jsonify({'status': f.read().strip()})
    return jsonify({'status': 'done'})


@app.route('/archive/view/<path:folder>')
def archive_view(folder):
    """View archive case results like single page."""
    import pandas as pd
    target = _safe_path(folder)
    if not target or not os.path.isdir(target):
        return redirect(url_for('archive'))

    results_dir = target
    case_name = os.path.basename(folder)

    # Load meta
    meta = {'name': case_name, 'mode': 'unknown', 'status': 'done'}
    meta_file = os.path.join(results_dir, 'meta.json')
    if os.path.exists(meta_file):
        with open(meta_file) as f:
            meta = json.load(f)
    meta['id'] = f'archive:{folder}'
    meta['name'] = meta.get('name', case_name)

    # Figures
    figures = []
    figures_dir = os.path.join(results_dir, 'figures')
    if os.path.isdir(figures_dir):
        for png in sorted(globmod.glob(os.path.join(figures_dir, '*.png'))):
            figures.append(os.path.basename(png))

    # Report
    report = ''
    report_path = os.path.join(results_dir, 'report.md')
    if os.path.exists(report_path):
        with open(report_path) as f:
            report = f.read()

    # CSVs
    tables = {}
    for csv_name in ['atom_statistics', 'contact_summary', 'coordination_summary', 'network_summary']:
        csv_path = os.path.join(results_dir, f'{csv_name}.csv')
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            tables[csv_name] = {'columns': df.columns.tolist(), 'data': df.values.tolist()}

    # Metrics
    metrics = {}
    metrics_path = os.path.join(results_dir, 'full_metrics.json')
    if os.path.exists(metrics_path):
        with open(metrics_path) as f:
            metrics = json.load(f)

    # Section headers
    if 'network_summary' in tables:
        data = tables['network_summary']['data']
        has_headers = any(str(row[0]).startswith('──') for row in data)
        if not has_headers:
            section_map = {
                'Porosity(%)': '── 구조 ──',
                'AM-SE Total(μm²)': '── 계면 ──',
                'SE-SE CN mean': '── 이온경로: 연결성 ──',
                'Tortuosity mean': '── 이온경로: 경로 효율 ──',
                'Path Hop Area mean(μm²)': '── 이온경로: 경로 품질 ──',
                'Ionic Active AM(%)': '── 활성도 ──',
                'Stress CV(%)': '── 응력 ──',
            }
            new_data = []
            for row in data:
                label = str(row[0])
                if label in section_map:
                    new_data.append([section_map[label], ''])
                new_data.append(row)
            tables['network_summary']['data'] = new_data

    # Patch placeholder values
    if 'network_summary' in tables and metrics:
        n_large = metrics.get('n_large_components')
        if n_large is not None:
            for row in tables['network_summary']['data']:
                if str(row[0]) == 'SE Cluster 수' and '≥10' not in str(row[1]):
                    row[1] = f"{n_large}(≥10) / {row[1]}"
        placeholder_map = {
            'GB Density(hops/μm)': 'gb_density_mean',
            'Path Hop Area mean(μm²)': 'path_hop_area_mean',
            'Path Bottleneck(μm²)': 'path_hop_area_min_mean',
            'Path Conductance(μm²)': 'path_conductance_mean',
        }
        for row in tables['network_summary']['data']:
            label = str(row[0])
            if label in placeholder_map and str(row[1]).strip() in ('-', ''):
                val = metrics.get(placeholder_map[label])
                if val is not None:
                    row[1] = val

    input_params = {}
    params_path = os.path.join(results_dir, 'input_params.json')
    if os.path.exists(params_path):
        with open(params_path) as f:
            input_params = json.load(f)

    return render_template('single.html', case=meta, figures=figures,
                         report=report, tables=tables, metrics=metrics,
                         input_params=input_params, archive_path=folder)


@app.route('/archive/results/<path:folder>/figures/<filename>')
def serve_archive_figure(folder, filename):
    target = _safe_path(os.path.join(folder, 'figures'))
    if not target:
        return 'Not found', 404
    return send_from_directory(target, filename)


@app.route('/archive/results/<path:folder>/3d-data')
def serve_archive_3d_data(folder):
    """Serve 3D data for archive case."""
    import pandas as pd
    target = _safe_path(folder)
    if not target:
        return jsonify({'error': 'Not found'}), 404

    meta = {}
    meta_file = os.path.join(target, 'meta.json')
    if os.path.exists(meta_file):
        with open(meta_file) as f:
            meta = json.load(f)
    scale = meta.get('scale', 1000)
    type_map_str = meta.get('type_map', '1:AM,2:SE')
    type_map = {}
    for item in type_map_str.split(','):
        k, v = item.split(':')
        type_map[int(k)] = v.strip()

    atoms_csv = os.path.join(target, 'atoms.csv')
    if not os.path.exists(atoms_csv):
        return jsonify({'error': 'No atom data'}), 404

    df = pd.read_csv(atoms_csv)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    particles = []
    for _, row in df.iterrows():
        t = int(row['type'])
        particles.append({
            'id': int(row['id']),
            'type': type_map.get(t, f'T{t}'),
            'x': round(float(row['x']) * scale, 2),
            'y': round(float(row['y']) * scale, 2),
            'z': round(float(row['z']) * scale, 2),
            'r': round(float(row['radius']) * scale, 2),
        })

    box = {
        'x_min': 0, 'x_max': round(0.05 * scale, 1),
        'y_min': 0, 'y_max': round(0.05 * scale, 1),
        'z_min': 0,
    }
    mesh_file = os.path.join(target, 'mesh_info.json')
    if os.path.exists(mesh_file):
        with open(mesh_file) as f:
            box['z_max'] = round(json.load(f)['plate_z'] * scale, 1)
    else:
        box['z_max'] = round(max(p['z'] + p['r'] for p in particles), 1)

    percolation = {'top_reachable': [], 'bottom_se': [], 'top_se': []}
    perc_path = os.path.join(target, 'percolation_sets.json')
    if os.path.exists(perc_path):
        with open(perc_path) as f:
            percolation = json.load(f)

    paths = []
    paths_path = os.path.join(target, 'tortuosity_paths.json')
    if os.path.exists(paths_path):
        with open(paths_path) as f:
            paths = json.load(f)

    clusters = {}
    clusters_path = os.path.join(target, 'se_clusters.json')
    if os.path.exists(clusters_path):
        with open(clusters_path) as f:
            clusters = json.load(f)

    return jsonify({
        'particles': particles, 'box': box,
        'percolation': percolation, 'paths': paths, 'clusters': clusters,
    })


@app.route('/archive/results/<path:folder>/save-screenshot', methods=['POST'])
def archive_save_screenshot(folder):
    """Save 3D screenshot to archive figures folder."""
    import base64
    target = _safe_path(folder)
    if not target:
        return jsonify({'error': 'Not found'}), 404
    data = request.get_json()
    if not data or 'image' not in data:
        return jsonify({'error': 'No image data'}), 400
    figures_dir = os.path.join(target, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    filename = data.get('filename', 'screenshot.png')
    filename = filename.replace('/', '_').replace('\\', '_')
    img_data = data['image'].split(',')[1]
    with open(os.path.join(figures_dir, filename), 'wb') as f:
        f.write(base64.b64decode(img_data))
    return jsonify({'success': True, 'filename': filename})


@app.route('/archive/results/<path:folder>/force-chains')
def serve_archive_force_chains(folder):
    target = _safe_path(folder)
    if not target:
        return jsonify([])
    fc_path = os.path.join(target, 'force_chains.json')
    if os.path.exists(fc_path):
        with open(fc_path) as f:
            return jsonify(json.load(f))
    return jsonify([])


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

@app.route('/download-doc/<path:filename>')
def download_doc(filename):
    """Download documentation files from project root."""
    doc_dir = os.path.dirname(os.path.dirname(__file__))
    return send_from_directory(doc_dir, filename, as_attachment=True)


@app.route('/group/fitting-report', methods=['POST'])
def group_fitting_report():
    """Generate GB correction fitting analysis report (downloadable MD)."""
    selected = request.form.getlist('cases')
    if not selected:
        return jsonify({'error': '케이스를 선택하세요.'}), 400

    # Collect metrics
    input_files = []
    names = []
    for cid in selected:
        if cid.startswith('archive:'):
            archive_rel = cid[len('archive:'):]
            case_path = os.path.join(app.config['ARCHIVE_FOLDER'], archive_rel)
            metrics_path = os.path.join(case_path, 'full_metrics.json')
            case_name = os.path.basename(archive_rel)
        else:
            case_path = get_results_dir(cid)
            metrics_path = os.path.join(case_path, 'full_metrics.json')
            meta_file = os.path.join(get_case_dir(cid), 'meta.json')
            case_name = cid
            if os.path.exists(meta_file):
                with open(meta_file) as f:
                    case_name = json.load(f).get('name', cid)
        if os.path.exists(metrics_path):
            input_files.append(metrics_path)
            names.append(case_name)

    if len(input_files) < 3:
        return jsonify({'error': '최소 3개 케이스가 필요합니다.'}), 400

    session_id = datetime.now().strftime('%y%m%d_%H%M%S') + '_fit'
    report_dir = os.path.join(app.config['RESULTS_FOLDER'], 'fitting_reports', session_id)
    os.makedirs(report_dir, exist_ok=True)

    scripts = app.config['SCRIPTS_FOLDER']
    cmd = ['python3', os.path.join(scripts, 'generate_fitting_report.py'),
           '-i'] + input_files + ['-n'] + names + ['-o', report_dir]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

    if result.returncode != 0:
        return jsonify({'error': f'리포트 생성 실패: {result.stderr}'}), 500

    report_path = os.path.join(report_dir, 'fitting_report.md')
    if not os.path.exists(report_path):
        return jsonify({'error': '리포트 파일이 생성되지 않았습니다.'}), 500

    with open(report_path, encoding='utf-8') as f:
        content = f.read()

    return jsonify({
        'success': True,
        'content': content,
        'download_url': f'/group/fitting-report-download/{session_id}'
    })


@app.route('/group/fitting-report-download/<session_id>')
def download_fitting_report(session_id):
    """Download fitting report as markdown file."""
    report_dir = os.path.join(app.config['RESULTS_FOLDER'], 'fitting_reports', session_id)
    return send_from_directory(report_dir, 'fitting_report.md', as_attachment=True,
                              download_name='GB_correction_fitting_report.md')


if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(debug=True, host='0.0.0.0', port=port)
