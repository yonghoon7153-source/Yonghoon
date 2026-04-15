# Lecture 5: SVM (continued), Features, Multiclass Classification

## 1. Recap: Separating Hyperplanes
- 여러 hyperplane이 training data를 완벽 분리 가능
- Empirical loss 같아도 test/expected loss는 다름
- → Large margin hyperplane 선택

## 2. Optimal Separating Hyperplane
- 두 클래스를 분리하면서 margin을 최대화
- 데이터에 너무 가깝지 않은 classifier → better generalization

## 3. Support Vector Machine
- 모든 점이 아닌 **boundary points**에만 집중
- **Support vectors** = margin boundary(±1) 위의 점들
- 이 점들만이 decision boundary를 결정

## 4. SVM: Hypothesis Class 원칙
- ||w||₂/2 ≤ R, s.t. y⁽ⁱ⁾(w₀+wᵀx⁽ⁱ⁾) ≥ 1
- R을 줄여가며 가장 작은 hypothesis class 찾기
- ŵ*는 가장 작은 R을 만족하는 최적 weight
- 타원 시각화: R↓ → 타원 축소 → ŵ*가 경계에 도달 = 최적

### 핵심 원칙
- **"올바른 해를 포함하는 가장 작은 hypothesis class를 선택하라"**
- Empirical loss ≠ Expected loss → large margin = high confidence
- SVM뿐 아니라 ML 전반에 적용
- 수학적 근거: VC-dimension theory
- Prior knowledge → constraint/regularizer로 추가

## 5. Soft-Margin SVM

### 문제: 데이터가 linearly separable하지 않을 때
- Hard-margin SVM 적용 불가

### 해결: Slack variables ξᵢ
- 일부 점이 margin 안에 들어오거나 잘못 분류되는 것 허용
- ξᵢ = 0: margin 밖 (올바름)
- 0 < ξᵢ < 1: margin 안이지만 올바른 쪽
- ξᵢ > 1: 잘못 분류됨

### Soft-Margin 목적함수
```
min  ||w||²/2 + γ·Σξᵢ
s.t. y⁽ⁱ⁾(w₀ + wᵀx⁽ⁱ⁾) ≥ 1 - ξᵢ  ∀i
     ξᵢ ≥ 0  ∀i
```

### γ (hyperparameter)
- γ = 0: constraint 무시, w=0, margin=∞ (의미없음)
- γ 작음: margin 최대화에 집중 (underfitting 위험)
- γ 큼: error 최소화에 집중 (overfitting 위험)
- γ → ∞: separable이면 hard-margin과 동일

## 6. Soft-Margin → Hinge Loss 변환

ξᵢ = max{0, 1 - y⁽ⁱ⁾(w₀ + wᵀx⁽ⁱ⁾)}

제약 없는 형태:
```
min Σᵢ max{0, 1 - y⁽ⁱ⁾(w₀+wᵀx⁽ⁱ⁾)} + ||w||²/(2γ)
    ────────────────────────────────   ──────────
         Hinge Loss                    Regularizer
```

**Hinge Loss: L_H(ŷ, y) = max{0, 1 - ŷy}**

## 7. Non-Linear Data → Kernel Trick

### Feature mapping φ
- 원래 공간에서 선형 분리 불가 → φ로 고차원 mapping → 선형 분리 가능

### XOR 예시
- 2D: (-1,1)●, (1,1)○, (-1,-1)○, (1,-1)● → linearly non-separable
- φ로 5차원 mapping → linearly separable!
- 핵심: √2·x₁x₂ 차원이 XOR 패턴 포착

### Feature mapping 예시
- 이미지 → Color Histogram (φ = feature extraction)
- 2D 원형 데이터 → 3D로 mapping → 평면으로 분리

## 8. Features as Part of the Model

```
전체: y = wᵀφ(x)
     = [non-linear φ] + [linear wᵀ·]
     = non-linear model
```

- 기존 ML: φ를 사람이 설계 (handcraft)
- **Deep Learning: φ도 데이터로부터 자동 학습!**
- Deep Learning ⊆ Representation Learning ⊆ ML ⊆ AI

## 9. Multiclass Classification

### One-hot encoding
- K개 클래스 → y = (0,...,0,1,0,...,0)
- k번째 원소만 1
- D개 input, K개 output → **K×D weight matrix W** + K-dim bias b

### 방법 1: 1-vs-all (1-vs-rest)
- K-1개 binary classifier
- 각각 "이 클래스 vs 나머지"
- **문제**: 애매한 영역 발생 (여러 답이 가능)

### 방법 2: 1-vs-1
- K(K-1)/2개 binary classifier (모든 쌍)
- 다수결 투표
- **문제**: transitivity 보장 안 됨 (가위바위보 현상)
  - C2 > C1, C1 > C3이라도 C2 > C3이 아닐 수 있음
