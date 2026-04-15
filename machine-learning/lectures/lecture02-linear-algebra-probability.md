# Lecture 2: Review on Linear Algebra & Probability

## Part 1: Linear Algebra

### 기본 구조

| 이름 | 차원 | 예시 |
|------|------|------|
| **Scalar** | 0-D | a, b, c (이탤릭) |
| **Vector** | 1-D | 숫자의 1차원 배열 |
| **Matrix** | 2-D | 숫자의 2차원 배열 |
| **Tensor** | n-D | 0차원 이상의 일반적 배열 |

### Matrix Transpose
- (Aᵀ)ᵢⱼ = Aⱼᵢ

### Matrix Product
- C = AB → Cᵢⱼ = Σₖ Aᵢₖ Bₖⱼ

### Identity Matrix
- AI = IA = A

### 연립방정식 (Systems of Equations)
- Ax = b
- 해의 종류: 없음 / 무한 / 유일 (1개)
- 유일해가 존재 → A가 invertible

### Matrix Inverse
- A⁻¹A = I
- Ax = b → x = A⁻¹b

---

### Norms (벡터의 크기)

| Norm | 수식 | 특징 |
|------|------|------|
| **Lp norm** | \|\|x\|\|_p = (Σ\|xᵢ\|ᵖ)^(1/p) | 일반 형태 |
| **L2 norm** | \|\|x\|\|₂ = √(Σxᵢ²) | 가장 많이 사용 (유클리드 거리) |
| **L1 norm** | \|\|x\|\|₁ = Σ\|xᵢ\| | 절댓값의 합 |
| **Max norm** | \|\|x\|\|∞ = max\|xᵢ\| | 가장 큰 원소의 절댓값 |

---

### Special Matrices

| 이름 | 성질 |
|------|------|
| **Unit vector** | \|\|x\|\|₂ = 1 |
| **Symmetric matrix** | A = Aᵀ |
| **Orthogonal matrix** | AᵀA = AAᵀ = I → A⁻¹ = Aᵀ |

---

### Eigendecomposition (고유값 분해)

- **고유벡터/고유값**: Av = λv
- **분해**: A = VΛV⁻¹ (V: 고유벡터 행렬, Λ: 고유값 대각행렬)
- **실수 대칭행렬**: 항상 실수 직교 고유분해 가능 → A = QΛQᵀ

---

### SVD (Singular Value Decomposition)

- A = UΣVᵀ
- Eigendecomposition보다 **일반적** — 정방행렬이 아니어도 가능
- U: 왼쪽 특이벡터, Σ: 특이값, V: 오른쪽 특이벡터

---

### Moore-Penrose Pseudoinverse

- A⁺ = VΣ⁺Uᵀ (Σ의 0이 아닌 원소의 역수를 취함)
- Ax = b에서:
  - 유일해 → 역행렬과 동일
  - 해 없음 → 오차가 가장 작은 해
  - 무한 해 → norm이 가장 작은 해

> **Linear Regression의 (XᵀX)⁻¹Xᵀ가 바로 pseudoinverse!**

---

### Trace
- tr(A) = Σᵢ Aᵢᵢ (대각 원소의 합)

---

## Part 2: Probability

### PMF vs PDF

| | PMF (이산) | PDF (연속) |
|---|---|---|
| **조건** | 0 ≤ P(x) ≤ 1 | p(x) ≥ 0 (1 초과 가능) |
| **정규화** | Σ P(x) = 1 | ∫ p(x)dx = 1 |
| **균등분포** | P(x=xᵢ) = 1/k | u(x;a,b) = 1/(b-a) |

### 핵심 규칙들

| 규칙 | 수식 |
|------|------|
| **Marginal (주변확률)** | P(x) = Σ_y P(x,y) |
| **Conditional (조건부확률)** | P(y\|x) = P(x,y) / P(x) |
| **Chain Rule** | P(x,y) = P(y\|x)P(x) |
| **Bayes' Rule** | P(y\|x) = P(x\|y)P(y) / P(x) |

### Independence (독립)
- P(x,y) = P(x)P(y)
- Conditional independence: P(x,y|z) = P(x|z)P(y|z)

### 기대값 / 분산

| 개념 | 의미 |
|------|------|
| **Expectation** | E[x] = Σ x·P(x), 선형성 성립: E[ax+b] = aE[x]+b |
| **Variance** | Var(x) = E[(x-E[x])²] |
| **Covariance** | Cov(x,y) = E[(x-E[x])(y-E[y])] |

### 주요 분포

| 분포 | 용도 |
|------|------|
| **Bernoulli** | 이진 변수 (성공/실패) |
| **Gaussian (Normal)** | 연속 변수, 가장 흔히 사용 |
| **Multivariate Gaussian** | 다변량 정규분포 — N(μ, Σ) |

### Gaussian Distribution
- 1D: N(μ, σ²) — 평균 μ, 분산 σ²
- Multivariate: N(μ, Σ) — 평균 벡터 μ, 공분산 행렬 Σ
- Precision matrix: Σ⁻¹ (공분산의 역행렬)

### Information Theory
- **Information**: I(x) = -log P(x) — 확률이 낮을수록 정보량 높음
- **Entropy**: H(x) = -Σ P(x)log P(x) — 불확실성의 척도
