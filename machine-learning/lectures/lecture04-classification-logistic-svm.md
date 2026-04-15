# Lecture 4: Classification, Logistic Regression, k-NN, SVM

## 1. Recap: ML 1-2-3
- Collect data + extract features (Experience + Prior Knowledge)
- Build model: choose H and loss l (Prior Knowledge)
- Optimization: minimize empirical loss (Prior Knowledge)

## 2. Why l₂ loss?
- Empirical: easy to optimize (closed-form for linear case)
- Theoretical: Gaussian + MLE → l₂ loss

## 3. MLE (Maximum Likelihood Estimation)
- likelihood(θ) = Πᵢ P_θ(x⁽ⁱ⁾, y⁽ⁱ⁾)  (i.i.d.이므로 곱)
- log 씌움 → Σᵢ log P_θ(...)  (곱→합)
- 부호 뒤집기 → Negative Log-Likelihood (NLL) = loss function
- Conditional: P(y|x)만 모델링 (discriminative)
- Gaussian 가정 → MSE, Bernoulli 가정 → Cross-Entropy

| 분포 가정 | MLE 결과 | 용도 |
|----------|---------|------|
| Gaussian | l₂ loss (MSE) | Linear Regression |
| Bernoulli (sigmoid) | Cross-Entropy | Logistic Regression |

## 4. Classification vs Regression
- Classification: 카테고리 예측 (이산값)
- Regression: 연속값 예측

## 5. Computer Vision
- ML의 subfield
- 컴퓨터는 이미지를 숫자 배열로 봄
- CV Task 계층: Classification → Localization → Scene Understanding → Action Recognition

## 6. Image Classification Challenges (7가지)
1. Viewpoint variation (각도)
2. Background clutter (배경)
3. Illumination (조명)
4. Occlusion (가려짐)
5. Deformation (형태 변형)
6. Intraclass variation (같은 클래스 내 다양성)
7. Scale variation (크기)

## 7. Representation
- 이미지 → 벡터 x ∈ ℝᵈ
- 벡터로 변환해야 선형대수 적용 가능
- Deep Learning: meaningful feature를 자동 학습

## 8. Non-parametric vs Parametric

| | Non-parametric | Parametric |
|---|---|---|
| 함수 형태 | 가정 없음 | 가정 있음 |
| 파라미터 수 | 데이터와 함께 증가 | 고정 |
| 예시 | k-NN, Decision Tree, Random Forest | Logistic Regression, LDA, Neural Networks |

## 9. k-Nearest Neighbors (k-NN)
- 가장 가까운 k개 이웃의 다수결로 분류
- 거리: Euclidean distance (sqrt 생략 가능 — 단조증가)
- Decision boundary: Voronoi diagram
- k 작으면 overfitting, k 크면 underfitting
- Rule of thumb: k < √n
- Hyperparameters: k, distance measure
- Normalization: zero mean, unit variance

## 10. Setting Hyperparameters
- Idea #1: training data 기준 → BAD (k=1이 항상 최고)
- Idea #2: test data 기준 → BAD (test에 overfitting)
- Idea #3: train/val/test 분리 → GOOD
- Idea #4: Cross-validation → 더 robust하지만 느림

## 11. Curse of Dimensionality
- 차원↑ → 공간 기하급수적 증가 → 데이터 sparse
- 필요한 데이터: O((1/ε)ᵈ)
- Saving grace: intrinsic dimension이 낮을 수 있음

## 12. Overfitting vs Underfitting
- Capacity: hypothesis space의 크기
- Overfitting: capacity 높음, train error↓ test error↑
- Underfitting: capacity 낮음, 둘 다 높음
- Bias-Variance Tradeoff: Generalization Error = Bias² + Variance + Noise
- U자 그래프: optimal capacity에서 generalization error 최소

## 13. Linear Classifier (Perceptron)
- f(x; w) = w₀ + wᵀx
- ŷ = sign(w₀ + wᵀx)
- Decision boundary: hyperplane (w₀ + wᵀx = 0)
  - 1D: 점(threshold), 2D: 직선, 3D: 평면, dD: 초평면
- w: boundary의 법선벡터 (수직)
- w₀: boundary의 위치 shift

## 14. Loss Functions for Classification
- 0-1 loss: 자연스럽지만 NP-hard, 불연속
- Squared Error: classification에 부적합 (맞아도 벌받음, outlier에 약함)
- Surrogate loss: 최적화 가능한 대리 loss 사용

## 15. Logistic Regression
- Sigmoid: σ(z) = 1/(1+e⁻ᶻ) — 출력을 [0,1]로 squash
- Cross-Entropy Loss: L_CE = -ylogŷ - (1-y)log(1-ŷ)
- Logistic-Cross-Entropy (수치 안정): L_LCE = ylog(1+e⁻ᶻ) + (1-y)log(1+eᶻ)
- MLE + Bernoulli → Cross-Entropy 유도
- Closed-form 없음 → gradient descent 필요

## 16. Parametric 접근법 (4단계)
1. Model: 변수 간 관계 설정
2. Loss function: fit 품질 측정
3. Regularizer: prior knowledge 인코딩
4. Optimization: 파라미터 학습

## 17. Hypothesis Class 선택 원칙
- Occam's razor: 단순한 모델 선호
- VC dimension theory
- Minimum description length
- Bias-variance tradeoff
- Curse of dimensionality

## 18. SVM (Support Vector Machine)

### Separating Hyperplanes
- 여러 hyperplane이 training data를 완벽 분리 가능
- Empirical loss는 같지만 generalization은 다름
- Margin이 큰 hyperplane이 더 좋음

### Math formulation
- Claim 1: w는 hyperplane에 수직 (법선벡터)
  - 증명: wᵀ(x₁-x₂)=0
- Claim 2: 점 x에서 hyperplane까지 거리 = |f(x)|/||w||
  - x = x⊥ + r·(w/||w||), f(x) = r·||w||

### Margin 유도
- 3개 평면: wᵀx+w₀ = +1, 0, -1
- x_pos = x_neg + λw
- λ = 2/(wᵀw)
- **Margin = 2/||w||**

### SVM 최적화 문제
```
min  ½||w||²
s.t. y⁽ⁱ⁾(wᵀx⁽ⁱ⁾ + w₀) ≥ 1,  ∀i
```
(Hard-margin SVM)
