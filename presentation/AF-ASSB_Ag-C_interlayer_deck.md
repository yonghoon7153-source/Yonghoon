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

**글로벌 ASSB 산업 동향과 본 연구의 포지셔닝**

발표자: Yonghoon | 2026.04

---

## 발표 목차

1. **왜 지금 전고체 배터리(ASSB)인가** — 시장의 요구
2. **글로벌 양산 로드맵** — 한/일/미·유/중 4각 경쟁
3. **ASSB 상용화의 기술 장벽** — 3대 한계
4. **무음극(Anode-Free) 전략** — 에너지밀도 × 비용
5. 🎯 **본 연구: Ag–C Interlayer 전략** (포스터 핵심)
6. **시장·공정 트렌드와의 연결점**
7. **결론 & Next Steps**

> 포스터의 실험 결과를 *현재 산업이 풀어야 할 문제*의 맥락에서 해석합니다.

---

## 1. 왜 지금 ASSB인가 — 시장의 압박

| 항목 | 현재 (LIB) | 2027–2030 목표 (ASSB) |
|---|---|---|
| 팩 단가 | $108/kWh ('25 BNEF) → $105 ('26E) | Nissan $75→$65/kWh (가솔린 동등) |
| 에너지밀도 | ~250 Wh/kg / 700 Wh/L | **375–500 Wh/kg / 900–1,000 Wh/L** |
| 안전성 | 유기 전해액 발화 위험 | 불연 고체전해질 |
| 충전 | 20–80% 30 min | 10–80% **10–18 min** (Toyota / Factorial) |

- **2027 = 첫 소량 SSB EV**, **2030 = 본격 양산**이 글로벌 컨센서스 <span class="small">(Interact Analysis, IDTechEx 2026)</span>
- 현 ASSB pack 비용은 LIB 대비 **3–5×** ($400–800/kWh) → **소재·공정 혁신으로 비용 구조를 깨야 함**

---

## 2. 글로벌 양산 로드맵 — 4각 경쟁 구도

<div class="small">

| 권역 | 대표 | 파일럿 | 양산 | 키 포인트 |
|---|---|---|---|---|
| 🇰🇷 | **Samsung SDI** | 수원 S-Line ('23) | **2027** (세계 최초 목표) | 황화물, 900 Wh/L, BMW 검증 |
| 🇰🇷 | **LG ES** | 청주 건식 파일럿 (Q4) | 폴리머 '26 / 황화물 **'30** | 건식 전극 선행 |
| 🇰🇷 | **SK On** | 대전 ('25.9) | **'29** (1년 단축) | WIP-free 공정, Solid Power 협력 |
| 🇰🇷 | **Hyundai** | 의왕 ('25.3) | **'30** | "Dream Battery", 자체 화학 |
| 🇯🇵 | **Toyota** | Idemitsu 공동 대형 | **'27–28 BEV** | 1,000 km, 10분 충전 |
| 🇯🇵 | **Nissan** | 요코하마 ('25.1) | **'28 FY** | 800 Wh/L, LiCAP 전극 |
| 🇯🇵 | **Honda** | 사쿠라시 27,400 m² | 20년대 후반 | ¥430억, 전공정 검증 |
| 🇺🇸 | **QuantumScape** | Eagle Line ('26.2) | 라이선스 (PowerCo) | QSE-5, Cobra 공정 |
| 🇺🇸 | **Factorial** | FEST 77 Ah | '26 데모 (Stellantis/MB) | **375 Wh/kg, 600 cyc, EQS 1,205 km** |
| 🇨🇳 | **CATL/BYD** | 소량 시범 | '27 소량 / **'30 주류** | Condensed + 황화물 |

</div>

→ **Samsung SDI 2027 양산**이 현재 선두. 한국 진영의 공통 병목 = **Li metal anode 계면 안정성**.

---

## 3. ASSB 상용화의 3대 기술 장벽

### (1) 계면 불안정 / 비균일 Li 핵형성
- ACS Energy Letters (2025.4): 황화물-AF 시스템에서 **불균일 Li 핵 → Li/SE 계면 파괴 → void 형성 → 급격한 용량 페이드**

### (2) Stack pressure & 기계적 열화
- Discover Electrochem. (2026): 1→20 MPa 운영 시 부피팽창 –41.8%, contact loss –94%
- Nat. Comm. (2025): 저스택압 운영의 핵심 변수 = **cathode chemomechanics**

### (3) 공정/대기 안정성
- 황화물 SE: 33% RH에서 2일 안정화 (표면 분자 엔지니어링, Nat. Comm. 2025)
- 건식 전극: 액체 슬러리와 황화물 부조화 → ASSB에 **dry process가 필수적**

> **본 연구가 타겟하는 문제 = (1) 계면 불안정 + 비균일 Li 핵형성**

---

## 4. 무음극(AF-ASSB) — 왜 고에너지의 마지막 퍼즐인가

| 이슈 | 기존 (Li-metal 음극) | 무음극 (AF) |
|---|---|---|
| 에너지밀도 | 300–350 Wh/kg | **400+ Wh/kg** (음극 무게·두께 제거) |
| 비용 | Li foil 가공 비용 | **Li foil 제거 → 원가↓** |
| 제조 난이도 | 불활성 환경 Li 취급 | 일반 공정 호환 |
| 문제 | Li 덴드라이트 | **초기 Li 핵형성 균질성이 전부** |

### 최신 무음극 학술 동향
- Nature Materials (2025.1): AF-SSB의 electro-chemo-mechanics 종합 리뷰
- **KIST/UNIST (2026.1): 1,500 cycle 후 75% 용량 유지** — 한국이 선도
- Nat. Comm. (2025.9): AI 기반 전해질 후보 7종 screening

→ **핵형성 균질화를 위한 interlayer 설계**가 AF-ASSB의 핵심 돌파구.

---

<!-- _class: lead -->

# 🎯 본 연구
## Uniformly Dispersed Ag in Ag–C Interlayer
## for Anode-Free All-Solid-State Batteries

<span class="small">*(포스터 기반 — 시장 관점에서 재해석)*</span>

---

## 5-1. 문제 정의 — "Ag 분산이 왜 중요한가"

### 기존 방식(PM, Physical Mixing)의 한계
- Ag 입자의 **심한 응집(agglomeration)**
- → 불균일 Li<sup>+</sup> flux / current density
- → **국부 Li 과핵형성 → 계면 파괴 → cycle life 급락**

### 본 연구의 해법: **AgNO₃–PVP 전구체 + In-situ 전기화학 환원(EC)**
- Slurry-casting 한 단계 공정 (산업 공정과 호환)
- PVP가 Ag 성장 억제 → **Ag 입자 30 nm 이하**
- 계면 전반에 **균질한 Li 핵형성 사이트** 제공

| 구성 | Ag 분산 | 결과 |
|---|---|---|
| PM Ag-C | 응집 (수 μm) | 비균일 Li 증착 |
| EC Ag-C | 양호 | 개선 |
| **EC Ag-C-PVP** | **<30 nm 균일** | **우수** |

---

## 5-2. 계면 안정성 — XPS로 본 "why it works"

### 사이클링 후 나이트레이트 유도 Interphase 형성
- **AgNO₃ → LiNO₂ / Li₃N 계 계면층** 자발 형성
- Li<sup>+</sup> 전도 향상 + 계면 저항 감소
- C–N bonding (PVP 유래) → 추가 안정화

### 초기 전기화학 거동
- **EC Ag-C-PVP**: 가장 낮은 nucleation overpotential
- 가장 낮은 R<sub>int</sub> (EIS)
- 고전류 Li plating/stripping 하에서도 안정

> **학술적 novelty**: interlayer 자체가 *in-situ*로 nitrate-derived SEI를 형성 — 외부 첨가 없이 계면 자가 치유.

---

## 5-3. AF-ASSB 풀셀 성능 — 저 Ag 로딩 (15 wt%)

| Cell | 초기 용량 | 수명 | 관찰 |
|---|---|---|---|
| PM Ag-C | 낮음 | **단락(short)** | Ag 응집 → 국부 Li deposit |
| EC Ag-C | ~250 mAh/g | 200+ cyc | 안정 |
| **EC Ag-C-PVP** | ~250 mAh/g | **250+ cyc, ~80% retention** | **최고** |

### 시장 관점에서의 임팩트
- **15 wt%라는 낮은 Ag 로딩** → **원가 민감도↓** (Ag 시세 변동 완화)
- Slurry-casting 단일 공정 → **LG ES 건식 전극 파일럿 / SK On WIP-free 라인**에 이식 가능
- 황화물 SE + AF 조합 → Samsung SDI·Toyota·Nissan의 2027–28 로드맵과 정렬

---

## 6. 본 연구 × 산업 트렌드 맵

<div class="small">

| 산업 이슈 | 본 연구의 기여 | 관련 외부 지표 |
|---|---|---|
| **비용 구조 파괴** ($400–800 → $75/kWh) | Ag 15 wt% 저로딩 + 단일 공정 | BNEF '25, Nissan target |
| **황화물 SE 대량 공급** | POSCO 7,200톤 / Lotte 1,200톤 ramp과 호환 | §2-3 DB |
| **건식/Simple 공정** | Slurry-casting 1-step, PVDF/PVP만 사용 | LG ES Ochang, Tesla 4680D |
| **Li 핵형성 균질화** | Ag 나노분산 + nitrate SEI | ACS Energy Lett. 2025 3대 한계 중 #1, #2 |
| **Stack pressure 민감도** | 균일 계면 → 저스택압 운영 잠재력 | Discover Electrochem. 2026 |
| **한국 AF 선도** | 1,500 cyc KIST/UNIST와 보완적 전략 | Korea AF 경쟁력 |

</div>

---

## 7. 결론 & Next Steps

### ✅ 달성
- **AgNO₃–PVP 전구체 + EC 전환**이라는 단일 공정 전략으로 <30 nm Ag 균일 분산 interlayer 제작
- Ag 15 wt% 저로딩에서도 **250+ cycles, 안정 사이클링** 입증
- In-situ **nitrate-derived SEI**로 계면 자가 안정화

### 🚀 다음 실험 (랩 내 후속)
- [ ] EIS fitting 정량화 (포스터 to-do 반영)
- [ ] 산화물/폴리머 SE 시스템으로 전이 가능성 테스트
- [ ] Scale-up: pouch 사이즈 + 저스택압 (<5 MPa) 운영 curve
- [ ] Ag 대체(Cu/Sn/Bi) 탐색 → 원가 추가 절감

### 🎯 포지셔닝 메시지
> *"한국 AF-ASSB 진영이 2027–2030 양산으로 가는 길목에서, **interlayer 차원의 nano-dispersion 전략**은 핵형성·계면·공정 3가지 이슈를 동시에 푸는 solution이다."*

---

## References (요약)

- **산업 로드맵**: Samsung SDI, LG ES, SK On, Hyundai, Toyota, Nissan, Honda, QuantumScape, Factorial, CATL/BYD 공식 PR (DB §1)
- **학술 — AF-ASSB**: Nat. Mater. (2025), ACS Energy Lett. (2025), Nat. Commun. (2025), Energy Mater. Adv. (2025)
- **학술 — Sulfide SE**: Nat. Commun. (2024, 2025), ACS Energy Lett. (2024), ACS Nano (2025), Cell Rep. Phys. Sci. (2025)
- **학술 — Stack pressure / Interface**: Adv. Energy Mater. (2025), Nat. Commun. (2025), Discover Electrochem. (2026), Adv. Sci. (2025)
- **시장**: BloombergNEF 2025, IDTechEx 2026, Interact Analysis, EnergyTrend
- **Full DB**: `refs/SSB_reference_DB.md` (본 repo)

---

<!-- _class: lead -->

# Thank you

Q & A

<span class="small">본 슬라이드는 `claude/ssb-market-research-ZiEJ4` 브랜치의 reference DB와
랩 내 Ag–C interlayer 포스터 데이터를 결합해 작성되었습니다.</span>
