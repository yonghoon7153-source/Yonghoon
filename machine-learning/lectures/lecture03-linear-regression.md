# Lecture 3: Linear Regression

## 1. Machine Learning이란?

> "A computer program is said to **learn from experience E** with respect to some class of **tasks T** and **performance measure P**, if its performance at tasks in T as measured by P, improves with experience E."
> — Tom Mitchell, 1997

핵심 3요소: **Task (T)**, **Performance (P)**, **Experience (E)**

---

## 2. ML의 핵심 개념들

### 2.1 Data 구분
| 용어 | 역할 |
|------|------|
| **Training data** | 모델 학습에 사용 |
| **Validation data** | 모델 선택(하이퍼파라미터 튜닝)에 사용 |
| **Test data** | 최종 평가(성능 측정)에 사용, 학습에 절대 사용하지 않음 |

### 2.2 Classification 종류
- **Binary classification**: 2개 클래스 (예: dog vs cat)
- **Multiclass classification**: 3개 이상 클래스

### 2.3 Supervised vs Unsupervised

| | Supervised Learning | Unsupervised Learning |
|---|---|---|
| **데이터** | Labeled data | Unlabeled data |
| **예시** | SVM, Linear Regression, Logistic Regression, k-NN, Neural Networks, Naive Bayes, LDA, Decision Trees | PCA, k-means clustering |

---

## 3. 핵심 도전 과제: Generalization (일반화)

- **Generalization**: 학습 시 보지 못한 새로운 데이터에 대해 잘 수행하는 능력
- **Generalization error**: 새로운 데이터에 대한 기대 오차 (직접 계산 불가)
- **Training error**: 훈련 세트에서 측정 → 일반화 오차의 나쁜 대리지표
- **Test error**: 테스트 세트에서 측정 → 일반화 오차의 더 나은 대리지표

---

## 4. ML의 수학적 정형화 (Math Formulation)

### 4.1 기본 설정
- Feature vector: **x** = (x₁, ..., x_d)
- Label: **y**
- Training data: {(x⁽¹⁾, y⁽¹⁾), ..., (x⁽ⁿ⁾, y⁽ⁿ⁾)} — **i.i.d.** from distribution D

> **i.i.d.** = independently identically distributed (독립 동일 분포)

### 4.2 목표
Hypothesis class **H** 에서 함수 f ∈ H를 찾되, **expected loss**를 최소화:

```
min_{f ∈ H} E[L(f(x), y)]
```

### 4.3 Loss Functions

| Loss Function | 수식 | 용도 |
|---------------|------|------|
| **0-1 loss** | L(f(x), y) = 1 if f(x) ≠ y, 0 otherwise | Classification |
| **ℓ₂ loss (MSE)** | L(f(x), y) = (f(x) - y)² | Regression |

### 4.4 Empirical Loss (경험적 손실)
Expected loss를 직접 계산할 수 없으므로, training data로 근사:

```
L_emp(f) = (1/n) Σᵢ L(f(x⁽ⁱ⁾), y⁽ⁱ⁾)
```

---

## 5. Machine Learning 1-2-3

1. **Collect data** and extract features
2. **Build model**: hypothesis class H와 loss function L 선택
3. **Optimization**: empirical loss 최소화

---

## 6. Classification vs Regression

| | Classification | Regression |
|---|---|---|
| **출력** | 이산값 (클래스 레이블) | 연속값 |
| **예시** | 고양이/개 분류 | 집값 예측, 암 크기 예측 |

---

## 7. Linear Regression (선형 회귀)

### 7.1 설정
- Input: **x** = (x₁, ..., x_d) (예: 임상 측정값 8개)
- Output: y (예: 암의 크기/정도)
- Hypothesis class: **f(x) = wᵀx** (선형 함수)
- Loss: **ℓ₂ loss (MSE)**

### 7.2 목적 함수
```
min_w  (1/n) Σᵢ (wᵀx⁽ⁱ⁾ - y⁽ⁱ⁾)²
```

### 7.3 행렬 표기
- **X**: n × d 행렬 (i번째 행 = x⁽ⁱ⁾)
- **y**: n × 1 벡터

목적 함수를 행렬로 쓰면:
```
min_w  (1/n) ||Xw - y||²
```

### 7.4 최적해 (Closed-form Solution)

Gradient를 0으로 놓으면:

```
∇_w ||Xw - y||² = 0
2Xᵀ(Xw - y) = 0
XᵀXw = Xᵀy
```

**Normal Equation:**
```
w* = (XᵀX)⁻¹ Xᵀy
```

> - XᵀX는 **symmetric** 행렬
> - X가 tall matrix (n > d)이면 XᵀX가 invertible할 수 있음
> - (XᵀX)⁻¹Xᵀ 를 **pseudo-inverse**라고 부름

### 7.5 Bias term 포함

bias가 있는 경우:
```
f(x) = wᵀx + b
```

이를 bias 없는 경우로 변환:
- x̃ = [x, 1] (feature에 1을 추가)
- w̃ = [w, b]
- 그러면 f(x) = w̃ᵀx̃

→ 이전의 모든 유도가 그대로 적용됨

---

## 8. 핵심 정리 (시험 포인트)

1. **ML 정의**: Task, Performance, Experience 세 가지로 정의
2. **i.i.d. 가정**: training과 test data가 같은 분포에서 독립적으로 추출
3. **Generalization**: ML의 핵심 도전과제, test error로 근사 측정
4. **Linear Regression의 closed-form**: w* = (XᵀX)⁻¹Xᵀy
5. **Bias trick**: feature에 1을 추가하면 bias 항을 weight에 흡수 가능
6. **Empirical loss**: 실제 기대 손실을 training data로 근사한 것
