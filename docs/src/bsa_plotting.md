# 模块三：`bsa_plotting.jl`（绘图扩展）

本文档专注介绍 `BSAHelper.BSAPlotting` 提供的绘图入口 `plot_bsa_data_collapse`：将 BSA/FSS 的坍缩数据绘制为（主图 + 可选残差图），并在需要时展示 Bootstrap 引入的误差与 X 方向误差棒。

作为背景，BSA 拟合的目标可以抽象为把不同系统尺寸的观测量写成统一的标度形式（示意）：

```
Y(X, L) = L^{c_2} × F[(X - X_c) L^{1/ν}, ...]
```

其中：
- `Y`：观测量（如磁化强度、关联比等）
- `X`：控制参数（如温度、相互作用强度）
- `L`：系统尺寸
- `X_c`（= `Tc` 或 `Uc`）：临界点
- `c_1 = 1/ν`：关联长度临界指数的倒数
- `c_2`：反常维度标度指数
- `F`：普适标度函数

本模块**仅提供唯一绘图函数**：`plot_bsa_data_collapse`, 用于绘制有限尺寸标度数据坍缩图（主图 + 残差图），支持**两种使用场景**：
1. **直接 BSA 拟合**：显示 MC 误差
2. **Bootstrap 分析**：显示 Bootstrap 误差

## 快速开始

绘图模块由 `BSAHelper` 的扩展提供，需要同时加载 `PyPlot` 与 `LaTeXStrings`：

```julia
using PyPlot
using LaTeXStrings
using BSAHelper

import BSAHelper: BSACore, BSABootstrap, BSAPlotting
```

## 数据流和控制流

### 数据流转换

**输入数据来源**： 以下两个输入数据来源, 也对应后面的两个使用场景. 
- `BSAHelper.BSACore.parse_bsa_output()`：解析 BSA 输出文件（仅包含 Y 误差）
- `BSAHelper.BSABootstrap.prepare_bootstrap_plot_data()`：准备 Bootstrap 绘图数据（注入 X 误差）

**核心数据流**：
```
输入: (metadata, data_sections)
  ├─ metadata 包含 p[0], p[1], p[2], ... (BSA 原始参数 + 误差)
  ├─ data_sections[1]: 标度后数据
  │   ├─ form=0 (无修正): [X, Y, E, L, T, A, dA] (7 列)
  │   ├─ form=1 (有修正): [X, Y, E, X2, L, T, A, dA] (8 列)
  │   └─ 扩展 (X 误差): 如果 BSAProblem 定义了 x_err_col，则追加 xerr 列
  └─ data_sections[2]: Scaling function (X_func, mu_func, sigma_func)
  
  ↓ 参数映射
metadata → BSACore.extract_parameter_dict() 
  → params (Tc, c1, c2, c3, ... 无歧义参数名)
  
  ↓ 物理诠释
params → BSACore.extract_physical_params()
  → phys (Uc/Tc/Jc, ν, η_φ/η_ψ)
  
  ↓ 绘图使用
phys → 标题和坐标轴标签
metadata (chi2, n_points, n_freeparams) → χ²_reduced 文本框
data_sections[1] 的 xerr 列 → 双向误差棒（主图和残差图）
```

---

### 绘图流程

**主图 (ax1)**：
1. 从 `data_sections[1]` 提取标度后数据点，按 L 值分组绘制
   - 自动检测是否包含 `xerr` 列
   - 如有 X 误差：绘制**双向误差棒**（X 和 Y 方向）
   - 如无 X 误差：仅绘制 Y 方向误差棒（传统模式）
2. 从 `data_sections[2]` 绘制 scaling function 曲线和置信区间
3. 使用 `phys` 构建标题（显示临界点、临界指数及误差）
4. 根据 `eta_type` 设置 Y 轴标签

**残差图 (ax2)** *(仅 plot_mode=:full)*：
1. 计算残差：`(Y - F(X)) / E`（F(X) 通过插值获得）
2. 绘制残差点（支持 X 误差棒），添加参考线 (y=0, y=±2)
3. 显示 χ²_reduced（自动从 `metadata` 提取 `n_freeparams`）

**拟合窗口指示器 (ax0)** *(仅 plot_mode=:simple)*：
1. 在主图上方显示拟合窗口，替代传统 legend
2. 使用 ticklabel 颜色区分窗口内外（黑色/浅灰色）
3. 绘制示意性双向误差棒（legend 风格，非定量）
4. 与主图 ticklabel 垂直对齐，无需重复 x 轴 ticks

## 使用场景

### 使用方式一：直接 BSA 拟合

```julia
# 运行 BSA 并解析输出
metadata, data_sections = BSACore.parse_bsa_output("output.op")

# 绘图（3 参数版本，自动从 metadata 生成 phys_fmt）
BSAPlotting.plot_bsa_data_collapse(
    metadata, data_sections, figs_dir;
    critical_param_name = "Uc",
    eta_type = :none
)
```

---

### 使用方式二：Bootstrap 分析

**关键步骤**：使用 `prepare_bootstrap_plot_data` 注入 Bootstrap 误差和 X 误差

```julia
# 1. Bootstrap 分析
bootstrap_result = BSABootstrap.bootstrap_bsa_analysis(problem, boot_cfg, bsa_cfg)

# 2. 准备绘图数据（注入 Bootstrap 误差 + X 误差，并生成 phys_fmt）
metadata, data_sections, phys_fmt = BSABootstrap.prepare_bootstrap_plot_data(
    problem, bootstrap_result, bsa_cfg;
    critical_param_name = "Uc",
    eta_type = :none
)

# 3. 绘图（显示 Bootstrap 误差 + 双向误差棒）
if !isempty(metadata) && !isempty(phys_fmt)
    BSAPlotting.plot_bsa_data_collapse(
        metadata, data_sections, phys_fmt, figs_dir;
        save_prefix = "bootstrap_fss",
        critical_param_name = "Uc",
        eta_type = :none
    )
end
```

**工作原理**：
- `prepare_bootstrap_plot_data` 会将 Bootstrap 均值和标准差写入 `metadata`，并基于 `BootstrapResult` 生成 `phys_fmt`
- `phys_fmt` 使用 `DataProcessforDQMC.statistics` 中的 `round_error` / `format_value_error`，根据 error-of-std 自洽确定有效数字，生成 `value_str` / `error_str`
- 如果 `problem.x_err_col !== nothing`，从原始数据提取 X 误差并追加到 `data_sections[1]`，绘图函数自动检测 `xerr` 列并绘制双向误差棒
- 绘图函数从 `data_sections` 读取数据散点（含 Bootstrap 误差）和平均拟合参数的标度函数，从 `metadata` 提取 χ² 信息，通过 `phys_fmt` 构建标题和 Y 轴缩放文字

## 两大绘图模式

`plot_bsa_data_collapse` 函数现在支持两种绘图模式：

1. **`:full` 模式**（默认）：全功能分析图
2. **`:simple` 模式**：论文发表用简洁图

### `:full` 模式 - 全功能分析图

**特点：**
- ✅ 散点数据（带误差棒）
- ✅ Scaling function 曲线 F(X)
- ✅ 残差图 (Residuals / σ)
- ✅ 完整标题（显示所有拟合参数：Uc, ν, c2/η, ω）
- ✅ χ²_reduced 信息框

**适用场景：**
- 数据分析验证
- 拟合质量检查
- 内部讨论汇报

**示例代码（单次 BSA 拟合，不使用 Bootstrap）**：
```julia
BSAPlotting.plot_bsa_data_collapse(
    metadata, data_sections, figs_dir;
    plot_mode=:full,  # 或省略（默认）
    critical_param_name="Uc",
    eta_type=:none,
    observable_label="G_{ab}"
)
```

**输出文件：**
- `bsa_fss_collapse.png`

---

### `:simple` 模式 - 论文发表图

**特点：**
- ✅ 散点数据（带误差棒）
- ❌ 无 Scaling function 曲线
- ❌ 无残差图
- ❌ 无标题
- ❌ 无 χ² 信息框

**适用场景：**
- 论文图片
- 演讲展示
- 视觉验证 data collapse 效果

**示例代码：**
```julia
BSAPlotting.plot_bsa_data_collapse(
    metadata, data_sections, figs_dir;
    plot_mode=:simple,
    critical_param_name="Uc",
    eta_type=:eta_phi,
    observable_label="G_{ab}"
)
```

**输出文件：**
- `bsa_fss_collapse_simple.png`

### 完整示例

```julia
using PyPlot
using LaTeXStrings
using BSAHelper
import BSAHelper: BSACore, BSAPlotting

# 读取 BSA 输出
metadata, data_sections = BSACore.parse_bsa_output("results/bsa_output.dat")

# 生成全功能图（用于内部分析）
BSAPlotting.plot_bsa_data_collapse(
    metadata, data_sections, "figs";
    plot_mode=:full,
    save_prefix="afm_analysis",
    critical_param_name="Uc",
    eta_type=:none,
    observable_label="G_{\\mathrm{AFM}}"
)

# 生成简洁图（用于论文）
BSAPlotting.plot_bsa_data_collapse(
    metadata, data_sections, "figs";
    plot_mode=:simple,
    save_prefix="afm_paper",
    critical_param_name="Uc",
    eta_type=:eta_phi,
    observable_label="G_{\\mathrm{AFM}}"
)
```

### 输出文件命名规则

- **`:full` 模式**: `{save_prefix}_collapse.png`
- **`:simple` 模式**: `{save_prefix}_collapse_simple.png`

示例：
- `:full` → `afm_analysis_collapse.png`
- `:simple` → `afm_analysis_collapse_simple.png`

---

### 注意事项

1. **`:simple` 模式专为论文设计**，省略了所有分析性元素，仅保留数据点
2. 在论文中使用 `:simple` 模式图片时，建议在 caption 中说明拟合参数
3. 两种模式可以同时生成，互不干扰
4. 图片尺寸暂时保持一致（8×9），后续可根据需要调整

## 主函数`plot_bsa_data_collapse`的参数说明

### 通用参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `plot_mode` | Symbol | `:full` | 绘图模式 (`:full` 或 `:simple`) |
| `save_prefix` | String | `"bsa_fss"` | 文件名前缀 |
| `critical_param_name` | String | `"Uc"` | 临界参数名称 |
| `eta_type` | Symbol | `:none` | η 类型 (`:none`, `:eta_phi`, `:eta_psi`) |
| `observable_label` | AbstractString | `"A"` | 可观测量标签（支持 String 或 LaTeXString） |
| `xlabel_custom` | AbstractString? | `nothing` | 自定义 x 轴标签（支持 String 或 LaTeXString） |

### 样式参数

| 参数 | 说明 |
|------|------|
| `tick_params` | 坐标轴刻度参数 |
| `font_legend` | 图例字体设置 |
| `mycolor` | 自定义颜色函数 |
| `mymarker` | 自定义标记函数 |
| `show_confidence` | 是否显示置信区间（仅 `:full` 模式） |
