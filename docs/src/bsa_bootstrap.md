# 模块二：`bsa_bootstrap.jl`

## Bootstrap FSS 分析流程

理解 Bootstrap 流程是掌握这些数据结构的关键。完整的 Bootstrap FSS 分析包含以下步骤：

```
1. 定义问题 (BSAProblem)
   ↓ 指定数据、列名
   
2. 配置 Bootstrap (BootstrapConfig)
   ↓ 设置样本数、jitter 参数
   
3. 重复 N 次迭代：
   ├─ 重采样数据：X_new ~ X ± σ_X, Y_new ~ Y ± σ_Y
   ├─ 随机化初值：Tc_init_new ~ Tc_init ± jitter
   ├─ 调用 bsa_core.jl（use_mc=false）拟合
   └─ 收集拟合参数
   
4. 统计汇总 (BootstrapResult)
   └─ 计算均值、标准差、成功率、衍生量（如 ν = 1/c_1）
```

**双重不确定性来源**：
- **数据不确定性**：X 和 Y 的测量误差（通过高斯噪声重采样处理）
- **拟合不确定性**：初始参数敏感性（通过 jitter 处理）

---

## 输出数据流总结

```
bootstrap_bsa_analysis()
  ↓ 每次迭代：BSACore.extract_parameter_dict(metadata)
  ↓ 收集无歧义参数：{"Tc", "c1", "c2", ...}
  ↓ 统计
BootstrapResult
  └─ param_means/param_stds: {"Tc", "c1", "c2", ...}  ← 存储的是无歧义参数
  
  ↓ extract_physical_params(result, ...)  ← 配套函数
  
physical_params: {"Uc", "nu", "eta_psi", ...}  ← 物理量（带 bootstrap 误差）
```

**关键设计**：
- ✅ Bootstrap 全流程使用无歧义参数（params 层级）
- ✅ 提供配套的 `extract_physical_params()` 进行物理诠释
- ✅ `save_bootstrap_summary()` 集成所有上下文，输出完整总结

---

## 核心数据结构

### `BSAProblem`：定义一个 FSS 问题

打包所有与**问题本身**相关的信息（数据、列名）。

```julia
Base.@kwdef struct BSAProblem
    name::String                        # 问题名称（日志标识）
    L_values::Vector{Int}               # 要分析的系统尺寸
    data::DataFrame                     # 原始数据（必须包含 L, X, Y, Y_err 列）
    x_col::Symbol                       # X 列名（如 :U, :R_AFM_correlation_ratio）
    y_col::Symbol                       # Y 列名（如 :m2, :R_AFM_correlation_ratio）
    y_err_col::Symbol                   # Y 误差列名（必需）
    x_err_col::Union{Symbol,Nothing}    # X 误差列名（可选，控制参数通常为 nothing）
end
```

**关键说明**：

- **X 有无误差**：
  - U-dependent FSS：`x_err_col = nothing`（U 是控制参数，无误差）
  - R-dependent FSS：`x_err_col = :R_AFM_err_correlation_ratio`（R 是观测量，有误差）

- **Y 标度处理**：
  - ❌ **不要**预先标度 Y（`Y' = Y * L^c2`），这会混淆参数物理意义
  - ✅ **正确做法**：通过 `BSAParameters` 的 `c2_init` 和 `c2_fixed` 控制标度
  - BSA 会自动处理 `Y = L^c2 * F(X, L^{-c3})` 的标度关系

- **R-dependent FSS 的简化方法**：
  - 设置 `c1 = 0` 和 `Tc = 0`，使得 BSA 直接使用 X = R，无需任何变换
  - 标度形式简化为：$Y(R, L) = L^{c_2} F(R)$

---

### `BootstrapConfig`：配置 Bootstrap 采样

控制 bootstrap 的**采样行为**和**参数随机化**策略。

```julia
Base.@kwdef struct BootstrapConfig
    n_samples::Int = 1000                          # Bootstrap 迭代次数（生产: 1000, 测试: 100）
    seed::Union{Int,Nothing} = nothing             # 随机种子（42 = 可重复）
    jitter_params::Dict{Symbol,Float64} = Dict()   # 多参数 jitter（见下文详解）
    y_sample_relative_error::Float64 = 0.0         # 每个 bootstrap 样本的 y 相对误差（0 = 视为精确观测）
    verbose::Bool = false                          # 每 100 次迭代打印进度
    keep_failed::Bool = false                      # 是否保留失败的拟合（调试用）
    tempdir::Union{Nothing,String} = nothing       # 临时文件目录（nothing = 自动清理）
end
```

**核心设计：多参数 Jitter**

`jitter_params` 是**最重要的成员**，决定如何随机化拟合参数的初值：

```julia
# 典型配置：同时对 Tc, c1, c2 进行 jitter
jitter_params = Dict(
    :Tc_init => 0.1,   # Tc 在 [原值 - 0.1, 原值 + 0.1] 内随机
    :c1_init => 0.05,  # 1/ν 在 ±0.05 范围内随机
    :c2_init => 0.05   # c2 在 ±0.05 范围内随机
)
```

**半径选择推荐**：临界点 ~3-5%，指数 ~0.05-0.1

---

### `BootstrapResult`：存储 Bootstrap 结果

```julia
struct BootstrapResult
    param_means::Dict{String,Float64}               # 无歧义参数均值（"Tc", "c1", "c2", ...）
    param_stds::Dict{String,Float64}                # 参数标准差
    n_success::Int                                  # 成功拟合次数
    n_trials::Int                                   # 总样本数
    bootstrap_samples::Vector{Dict{String,Float64}} # 所有成功样本（用于进一步分析）
end
```

**健康指标**：成功率 = `n_success / n_trials`
- \> 90%：正常
- 70-90%：可接受
- < 70%：初值或数据质量问题

## 使用示例

### 场景 1：U-dependent FSS（关联比 vs U）

- **`AFMCorrelationRatio-FSS-analysis.jl`**：U-dependent bootstrap FSS
  - X = U（无误差）
  - Y = R_AF（关联比）
  - 提取 Uc 和 ν

```julia
# BSA 参数配置
params = BSACore.BSAParameters(
    Tc_init = 3.8,
    c1_init = 1.0,    # 1/ν
    c2_init = 0.0,
    c2_fixed = true,  # 固定 c2 = 0（关联比无需标度）
    ...
)

# 问题定义
problem = BSAProblem(
    name = "AFM Correlation Ratio",
    data = data,  # 用户应在外部准备干净的数据
    x_col = :U,
    x_err_col = nothing,      # U 无误差（控制参数）
    y_col = :R_AFM,
    y_err_col = :R_AFM_err
)
```

### 场景 2：R-dependent FSS（m² vs R）

- **`AFMm2-CorrelationRatio-FSS-analysis.jl`**：R-dependent bootstrap FSS
  - X = R_AF（有误差！）
  - Y = m²（磁化强度平方）
  - 提取 η_φ（反常维度）

**关键技巧**：设置 `c1 = 0` 和 `Tc = 0`，使得标度形式简化为 $Y = L^{c_2} F(R)$，X 直接使用 R 无需变换。

```julia
# BSA 参数配置
params = BSACore.BSAParameters(
    Tc_init = 0.0,
    Tc_fixed = true,  # X_c = 0
    c1_init = 0.0,    # c1 = 0，使得 BSA 直接使用 X = R
    c1_fixed = true,
    c2_init = -1.7,   # c2 = -(1 + η_φ)
    c2_fixed = false,
    ...
)

# 问题定义
problem = BSAProblem(
    name = "m² vs R",
    data = merged_data,  # 用户应在外部准备干净的数据（如筛选 U > 3.7）
    x_col = :R_AFM,
    x_err_col = :R_AFM_err,  # R 有误差（观测量）
    y_col = :m2,
    y_err_col = :m2_err
)
```

---

## 临时文件管理

Bootstrap 为每次迭代生成临时文件：
- `bootstrap_N.dat`：FSS 输入文件
- `bootstrap_N.op`：BSA 输出文件
- `bootstrap_N.log`：BSA 日志文件

这些文件在每次迭代后自动清理。调试时可以通过设置 `BootstrapConfig.tempdir` 保留文件。

