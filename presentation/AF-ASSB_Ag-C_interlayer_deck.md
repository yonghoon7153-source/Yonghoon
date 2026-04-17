---
marp: true
theme: default
paginate: true
size: 16:9
style: |
  section { font-size: 22px; }
  h1 { color: #1a3a7a; }
  h2 { color: #1a3a7a; border-bottom: 2px solid #1a3a7a; padding-bottom: 4px; }
  table { font-size: 18px; }
  .small { font-size: 16px; }
  .hl { background: #fff3b0; padding: 2px 4px; border-radius: 3px; }
---

<!-- _class: lead -->

# 균일 Ag 분산 Ag–C Interlayer를 이용한
# 무음극 전고체 전지(AF-ASSB) 계면 안정화

**랩 실험 결과 → 글로벌 ASSB 산업 동향으로의 확장**

발표자: Yonghoon | 2026.04

---

## 발표 흐름

### Part A. 연구 (Lab work)
1. 문제 정의 — 무음극 ASSB 계면의 비균일 Li 핵형성
2. 제안: **AgNO₃–PVP + In-situ 전기화학 환원(EC)** interlayer
3. Ag 분산 / 입자 크기 (SEM·TEM)
4. 계면 화학 (XPS) + 초기 전기화학 / EIS
5. 풀셀 사이클링 (Ag 15 wt%, 250 cycles)

### Part B. 시장 (Zoom-out)
6. ASSB 상용화의 3대 장벽 — 본 연구가 겨냥한 지점
7. 글로벌 양산 로드맵 (한/일/미·유/중)
8. 비용·공정 트렌드와 본 연구의 접점
9. 결론 & Next Steps

---

<!-- _class: lead -->

# Part A
## 랩 실험: Ag–C Interlayer for AF-ASSB

---

## 1. 문제 정의 — 왜 무음극에서 "Ag 분산"인가

### AF-ASSB의 근본 난제
- 음극 Li foil 없음 → **매 사이클 in-situ Li plating / stripping**
- Li 핵형성이 **불균일**하면 → 국부 과증착 → Li/SE 계면 파괴 → void → 급격한 capacity fade

### Ag interlayer의 역할
- Ag가 Li와 합금(Li–Ag) 형성 → **Li 핵형성 사이트** 제공
- 단, Ag가 **응집**되면 역효과: 응집 지점에 Li 몰림 → 오히려 단락 유발

### 기존 공정(PM, Physical Mixing)의 한계
- Ag 입자 수 μm 단위 응집
- → 불균일 Li<sup>+</sup> flux / current density → 계면 열화

→ **핵심 질문**: *간단한 슬러리 공정으로 Ag를 수십 nm 단위로 균일 분산시킬 수 있는가?*

---

## 2. 제안 — AgNO₃–PVP 전구체 + In-situ EC 환원

### 설계 포인트
1. **AgNO₃ 용액** → carbon matrix 위에 원자 단위로 균질 도포
2. **PVP** = Ag 성장 억제제 (입자 크기 <30 nm)
3. 셀 내부 첫 방전 시 **전기화학적으로 Ag⁰ 환원** (별도 열처리·환원제 불필요)
4. 부산물 NO₃⁻ → 사이클 중 **nitrate-derived SEI** 자발 형성 (뒤에서 XPS로 입증)

### Fabrication flow
PVDF(NMP) + PVP + Carbon + AgNO₃ → **Mixing → Slurry casting → Drying**
→ 단일 step, 기존 산업 공정과 호환

| 구성 | Ag 분산 | 비고 |
|---|---|---|
| PM Ag-C | 응집(수 μm) | 기준 |
| EC Ag-C | 개선 | PVP 無 |
| **EC Ag-C-PVP** | **<30 nm 균일** | 본 연구 |

---

## 3. Ag 분산 결과 — SEM·TEM

### SEM/EDS (포스터 section 2)
- **PM Ag-C**: Ag가 수 μm 덩어리로 응집, 탄소 매트릭스 일부 영역 Ag 無
- **EC Ag-C**: Ag 입자 크게 분산, 하지만 국부 집중 존재
- **EC Ag-C-PVP**: Ag EDS 신호가 필드 전체에 **균일** 분포

### TEM (고배율)
- EC Ag-C-PVP에서 **Ag 입자 < 30 nm** 확인
- PVP의 steric hindrance 효과로 Ag nuclei 성장 억제

> **Takeaway**: 공정 스텝 수를 늘리지 않고 단일 슬러리 캐스팅으로 nano-dispersion 달성.

---

## 4. 계면 화학 — XPS (nitrate-derived SEI)

### Before cycling
- N 1s: **AgNO₃** 피크 + **C–N bonding** (PVP 유래)

### After cycling
- N 1s: **LiNO₂, Li₃N** 등 **lithium nitrate/nitride 계 SEI**로 전환
- 즉, 별도 전해질 첨가제 없이 **interlayer 자체가 in-situ SEI 전구체** 역할

### 의미
- LiNO₂/Li₃N은 **고속 Li⁺ 전도** + **계면 안정화**로 알려진 전형적 유익 SEI
- 기존 문헌 대비 차별점: *외부 nitrate 첨가제 없이 interlayer 설계만으로 확보*

---

## 5. 전기화학 거동 — 초기 / EIS / Li plating-stripping

### 초기 방전 overpotential
- **EC Ag-C-PVP**: 가장 낮은 nucleation overpotential
- PM > EC Ag-C > EC Ag-C-PVP 순으로 감소

### EIS (Nyquist, 45 °C)
- R<sub>int</sub>(계면 저항): PM ≫ EC Ag-C > **EC Ag-C-PVP (최저)**
- <span class="small">※ 포스터 to-do: EIS fitting 정량화 → equivalent circuit (R<sub>b</sub>/R<sub>int</sub>/CPE) 값 표 추가 예정</span>

### Li plating / stripping (0.7 → 1.5 → 2.0 mAh cm⁻²)
- 고전류에서도 EC Ag-C-PVP는 안정 전압 프로파일
- PM은 고전류에서 스파이크 / 단락

→ **균일 Ag 분산 = 낮은 과전압 + 낮은 계면 저항 + 고전류 내구성** 3박자.

---

## 6. 풀셀 — AF-ASSB 사이클링 (Ag 15 wt%, 0.2 C, 45 °C)

| Cell | 초기 용량 | 100 cyc | 250 cyc | 관찰 |
|---|---|---|---|---|
| PM Ag-C | 낮음 | — | **단락** | Ag 응집 → 국부 Li deposit |
| EC Ag-C | ~250 mAh g⁻¹ | 유지 | 200+ cyc 안정 | 개선 |
| **EC Ag-C-PVP** | **~250 mAh g⁻¹** | **안정** | **250+ cyc, ~80% retention** | **최고** |

### 핵심
- **Ag 15 wt%라는 저로딩**에서 250 cycles 달성 — 종래 문헌의 30–50 wt% 대비 절반 이하
- Coulombic efficiency도 200 cycle 이후까지 안정
- PVP 효과가 특히 **저 Ag 로딩 영역**에서 극대화

---

## 7. Part A 요약

<div class="small">

| 층위 | 관찰 | 메커니즘 |
|---|---|---|
| **Ag 분산** | <30 nm 균일 | PVP steric + AgNO₃ EC 환원 |
| **Li 핵형성** | 균질 | 분산된 Ag seed |
| **SEI** | LiNO₂/Li₃N | AgNO₃/PVP 유래 nitrate conversion |
| **계면 저항** | 최저 | 균일 접촉 + 유익 SEI |
| **풀셀** | 250+ cyc @ 15 wt% | 위 4층의 synergy |

</div>

> **학술적 novelty** — interlayer 하나가 *(i)* Li 핵형성 사이트, *(ii)* in-situ SEI 전구체, *(iii)* Ag 원가 절감 3역할을 동시에 수행.

---

<!-- _class: lead -->

# Part B
## Zoom-out: 이 결과가 시장에서 왜 중요한가

---

## 8. ASSB 상용화의 3대 기술 장벽 — 본 연구의 위치

### (1) 계면 불안정 / 비균일 Li 핵형성  ← **Part A가 정면으로 다룸**
- ACS Energy Letters (2025.4, 황화물-AF 리뷰): **#1 장벽으로 지목**
- 본 연구 = interlayer 레벨 해법

### (2) Stack pressure & 기계적 열화
- Discover Electrochem. (2026): 1→20 MPa에서 contact loss –94%
- Nat. Comm. (2025): 저스택압 운영 = cathode chemomechanics

### (3) 공정·대기·비용 안정성
- 황화물 SE 33% RH 2일 안정 (Nat. Comm. 2025)
- **건식 전극이 ASSB에 필수** (슬러리 ↔ 황화물 부조화)

> Part A의 슬러리 공정은 (3)과 절충점; (1)은 직접 해결, (2)는 균일 계면 덕에 간접 수혜.

---

## 9. 글로벌 양산 로드맵 — 본 연구가 올라탈 길

<div class="small">

| 권역 | 대표 | 파일럿 | 양산 | 키 포인트 |
|---|---|---|---|---|
| 🇰🇷 | **Samsung SDI** | 수원 S-Line ('23) | **2027** (세계 최초 목표) | 황화물, 900 Wh/L, BMW 검증 |
| 🇰🇷 | **LG ES** | 청주 건식 파일럿 (Q4) | 폴리머 '26 / 황화물 **'30** | 건식 전극 선행 |
| 🇰🇷 | **SK On** | 대전 ('25.9) | **'29** (1년 단축) | WIP-free 공정, Solid Power 협력 |
| 🇰🇷 | **Hyundai** | 의왕 ('25.3) | **'30** | "Dream Battery" |
| 🇯🇵 | **Toyota** | Idemitsu 공동 | **'27–28 BEV** | 1,000 km, 10분 충전 |
| 🇯🇵 | **Nissan** | 요코하마 ('25.1) | **'28 FY** | 800 Wh/L, $75→$65/kWh |
| 🇯🇵 | **Honda** | 사쿠라시 27,400 m² | 20년대 후반 | ¥430억, 전공정 |
| 🇺🇸 | **QuantumScape / Factorial** | Eagle Line / FEST | '26 샘플 / 데모 | 375 Wh/kg, EQS 1,205 km |
| 🇨🇳 | **CATL / BYD / WeLion** | 소량 시범 | '27 소량 / **'30 주류** | Semi → Full solid |

</div>

- **한국 3사 + 현대가 AF 방향으로 공통 진입** → **interlayer 전략의 수요처**
- 학술 돌파: KIST/UNIST (2026.1) — AF 1,500 cycle 75% retention (Korea가 선도)

---

## 10. 비용·공정 트렌드와 본 연구의 접점

| 산업 이슈 | 현재 상황 | 본 연구의 기여 |
|---|---|---|
| **팩 단가** | $108/kWh ('25, BNEF) → Nissan 목표 $65/kWh | **Ag 15 wt% 저로딩** (통상 30–50 wt%의 절반 이하) |
| **현 ASSB 비용** | LIB 대비 **3–5×** ($400–800/kWh, IDTechEx 2026) | 원가 민감도 완화 방향 |
| **건식/Simple 공정** | LG ES 17–30% 원가↓, 4680D, Tesla 4680 | **Slurry-casting 1-step** → 이식 가능 |
| **황화물 SE 공급** | POSCO 7,200톤 / Lotte 1,200톤 (2027) | 공급 ramp와 시점 정렬 |
| **Li 핵형성 균질화** | ACS Energy Lett. 2025가 꼽은 #1 장벽 | **Part A가 직접 해결** |
| **한국 AF 선도** | KIST/UNIST 1,500 cyc | interlayer 차원의 보완적 전략 |

---

## 11. 결론

### ✅ 랩 성과 (Part A)
- **AgNO₃–PVP + EC 환원**: 단일 슬러리 공정으로 <30 nm Ag 균일 분산
- **In-situ nitrate-derived SEI** (LiNO₂/Li₃N)로 계면 자가 안정화
- AF-ASSB 풀셀: **Ag 15 wt% / 250 cycles / ~80% retention**

### 🌐 시장 의의 (Part B)
- AF-ASSB의 **#1 장벽(핵형성 불균일)** 정면 공략
- **저 Ag 로딩 + 슬러리 단일 공정** = 한국 3사 2027–2030 로드맵에 이식 가능
- 황화물 SE 국내 공급 ramp(POSCO·Lotte)와 타이밍 정렬

### 🚀 Next Steps
- [ ] **EIS fitting 정량화** (포스터 to-do)
- [ ] **Scheme 수정** + 전기화학 / 특성 분석 데이터 분리
- [ ] 저스택압(<5 MPa) 운영 curve 확보 → 장벽 (2) 검증
- [ ] Ag 대체재(Cu/Sn/Bi) 탐색 → 원가·공급망 리스크 추가 완화

---

## References (요약)

- **학술 — AF-ASSB**: Nat. Mater. (2025), ACS Energy Lett. (2025), Nat. Commun. (2025), Energy Mater. Adv. (2025)
- **학술 — Sulfide SE**: Nat. Commun. (2024, 2025), ACS Energy Lett. (2024), ACS Nano (2025), Cell Rep. Phys. Sci. (2025)
- **학술 — Stack pressure / Interface**: Adv. Energy Mater. (2025), Nat. Commun. (2025), Discover Electrochem. (2026), Adv. Sci. (2025)
- **산업 로드맵**: Samsung SDI, LG ES, SK On, Hyundai, Toyota, Nissan, Honda, QuantumScape, Factorial, CATL/BYD 공식 PR
- **시장**: BloombergNEF 2025, IDTechEx 2026, Interact Analysis, EnergyTrend
- **Full DB**: `refs/SSB_reference_DB.md`

---

<!-- _class: lead -->

# Thank you

Q & A

<span class="small">Part A (lab) → Part B (market) 흐름.
근거 DB: `refs/SSB_reference_DB.md` (from `claude/ssb-market-research-ZiEJ4`)</span>
