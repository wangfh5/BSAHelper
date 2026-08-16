# Form-1 修正点（Y−Δ）不确定度估计 与 `bootstrap_bsa_analysis` 的已知缺陷

状态：待办（机制曾实现后于 2026-08 移除，本文档保留设计要点供将来重做）

## 背景与当前共识

form-1 的 data collapse 图（`correction_view=:subtracted`）绘制修正点 `Y−Δ` 对 `X1`，其中
`Δ = X2·F1(X1)`，`X2 = L^(−c3)`。当前约定：

- collapse 图是**固定参数可视化**：所有拟合参数（c3 及 GP 超参数 θ）固定为 bootstrap 均值后
  重建曲线与修正项，用途是验证均值参数的可靠性，**不**在图上表达拟合参数的不确定度。
- 因此图中修正点的 yerr 直接用原始观测误差 `E`，这是正确约定（类比：correlation ratio
  collapse 图的横轴也不画 ν 的 bootstrap 误差）。
- 本文档描述的机制**只服务于稳健性审计**（回答"修正点的不确定度到底多大"），不接进绘图管线。

## 一、Y−Δ 不确定度的正确估计：配对预测 bootstrap

`Δ` 的不确定度来自 c3 和 GP 超参数的拟合不确定度。正确做法是**配对 bootstrap**：

1. 每个 bootstrap 样本：重采样数据 → 完整重拟合（参数 + GP）→ 用该样本自己的拟合结果
   计算逐点预测（`mu_full`、`correction`、`mu_zero`）。
2. 逐点预测按 **sample index** 存储（而非线程完成顺序），保证同一点、同一样本的
   `Y_corrected`、`correction`、模型预测严格配对。
3. `Y−Δ` 的不确定度 = 配对样本间的标准差（`corrected_std`）；多点协方差跨样本计算，
   所有 corrected-point 的协方差共用同一个 outlier mask（`included`）。
4. `Y_corrected`、`correction`、模型预测来自同一样本，**不能**把各自标准差独立平方相加。

## 二、为什么不能用 √(E² + std_correction_latent²)

`point_predictions` 区段输出的 `std_correction_latent` 看似可以和观测误差 E 做误差传播，
但这样得到的 error bar 是错的，原因有三：

1. **latent std 是条件不确定度**：它是 GP 在"参数（c3、θ）固定"条件下的潜变量标准差，
   完全不含参数本身的 bootstrap 不确定度——而后者往往才是主导项。
2. **Y 与 Δ 不独立**：GP 是用同一份 Y 训练出来的，Y 涨则拟合出的 F1 也跟着涨，
   Var(Y−Δ) ≠ Var(Y) + Var(Δ)，交叉项不可忽略。
3. **训练点收缩**：在训练数据点上，GP 后验均值向观测值收缩，latent std 在数据点附近
   系统性偏小，直接套用在数据点上会低估不确定度。

## 三、老函数 `bootstrap_bsa_analysis` 的三个已知缺陷

1. **线程全局 RNG，结果不可复现**：重采样用线程共享的全局 RNG，样本序列依赖线程调度，
   同一 seed 跑两次结果不同。修法：每个 sample 的 RNG 由 sample index 派生
   （如 `MersenneTwister(seed + sample_idx)`），与调度无关。
2. **X、Y 独立重采样，丢失同源相关性**：m²–R 通道中 X=R 与 Y=m² 来自同一次蒙特卡洛测量，
   两者误差相关；独立抽样会高估 χ² 的离散度。修法：联合 (X,Y) 重采样
   （`xy_covariance`；m²–R 可传 `:m2R`，复用 `compute_cr_m2_covariance` 逐点算协方差后联合抽样）。
   注：G_AB–R 通道已有定量判据支持忽略 Cov(X,Y)（x 传播项主导 σ²_eff，r≈3.5，畸变 ~0.5ρ）。
3. **`y_sample_relative_error` 约定反直觉**：老函数 `y_sample_err = |y|·rel_err`，
   rel_err=0 时样本 y 误差为零（视为精确观测），默认行为反而丢弃了原始异方差误差。
   修法：默认保留每个点原始的异方差 `Y` 误差，只有显式设置 `y_sample_relative_error`
   才切换到相对误差模式。

## 四、已移除的参考实现要点（2026-08 删除）

曾实现 `bootstrap_bsa_analysis_with_predictions`（opt-in）+ `BootstrapPredictionResult`，
已按当前共识（collapse 图不需要该机制）移除。重做时的要点：

- RNG 由 sample index 决定，不受线程调度影响；
- 默认保留原始异方差 `y_err`（修复缺陷 3）；
- `xy_covariance` 用原始 X/Y 坐标，`:m2R` 走联合抽样（修复缺陷 2）；
- `successful` 记录实际成功样本，`included` 是所有 corrected-point covariance 共用的
  outlier mask；`corrected_std` 是配对样本间标准差，`interval_kind` 明确记录这一约定；
- 调用方需把 `context`、`result`、`predictions` 一起存入 JLD2。

## 验收标准（将来重做时）

- [ ] 同一 seed 多线程跑两次，逐点预测逐位一致（可复现性）。
- [ ] m²–R 通道配对 bootstrap 的 χ² 分布与独立抽样对比，量化 Cov(X,Y) 的影响。
- [ ] 修正点 `corrected_std` 显著大于 √(E²+std_latent²) 的朴素估计（验证第二节的三条理由）。
