# 模块一：`bsa_core.jl`

## 设计目标

封装外部 BSA（Bayesian Scaling Analysis）命令行工具，提供结构化的 Julia 接口。

BSA 命令行格式：
```bash
./bsa [Options] [Data file] [Parameters]
```

本模块将上述命令行参数封装为两个结构体：
- **`BSAConfig`**：对应 Options（命令行选项）
- **`BSAParameters`**：对应 Parameters（拟合参数序列）

---

## 核心数据结构

### `BSAConfig`：BSA 命令行选项

```julia
Base.@kwdef struct BSAConfig
    binary::String = "/Users/ssqc/bin/bsa"  # 可执行文件路径
    scaling_form::Int = 0                   # -f SCALING::FORM [0:standard, 1:with correction]
    use_mc::Bool = true                     # -c (estimate confidential intervals by MC)
    xscale::Float64 = 1.0                  # -w OUTPUT::XSCALE (xscale of outputted scaling function)
end
```

**成员变量映射到 BSA 选项**：

| 成员变量 | 对应 BSA 选项 | 说明 |
|---------|-------------|------|
| `binary` | 可执行文件路径 | 指定 BSA 程序位置 |
| `scaling_form` | `-f SCALING::FORM` | `0`: 标准形式, `1`: 带修正形式 |
| `use_mc` | `-c` | 是否启用蒙特卡洛误差估计（bootstrap 时设为 `false`） |
| `xscale` | `-w OUTPUT::XSCALE` | 标度函数输出范围（拟合通常用 `1.0`，绘图用 `1.8-2.0`） |

---

### `BSAParameters`：BSA 拟合参数序列

BSA 命令行参数格式：`mask initial_value`
- `mask`：`0` = 固定，`1` = 自由
- `initial_value`：初始猜测值

本结构体将这些参数打包：

```julia
Base.@kwdef struct BSAParameters
    Tc_init::Float64 = 3.85          # p[0]: Xc (临界点)
    Tc_fixed::Bool = false
    
    c1_init::Float64 = 0.9           # p[1]: 1/ν (关联长度指数倒数)
    c1_fixed::Bool = false
    
    c2_init::Float64 = 0.1           # p[2]: c2 (反常维度标度指数)
    c2_fixed::Bool = false
    
    c3_init::Float64 = 0.5           # p[3]: ω (修正指数, 仅 scaling_form=1)
    c3_fixed::Bool = false
    
    theta0_init::Float64 = 1.0       # theta[0-4]: 核函数超参数
    theta1_init::Float64 = 1.0
    theta2_init::Float64 = 1.0
    theta3_init::Float64 = 1.0
    theta4_init::Float64 = 1.0
    theta_fixed::Bool = false
end
```

**参数映射说明**：

Julia 结构体到 BSA 命令行参数的映射**取决于 `scaling_form`**：

**1. `scaling_form = 0`（标准形式）**

```
Params[0] = Tc      ← Tc_init
Params[1] = c1      ← c1_init (1/ν)
Params[2] = c2      ← c2_init
Params[3] = θ0      ← theta0_init
Params[4] = θ1      ← theta1_init
Params[5] = θ2      ← theta2_init
```

```bash
./bsa data.dat \
    1 0.42 \    # p[0]: Tc
    1 0.9  \    # p[1]: c1 (1/ν)
    1 0.1  \    # p[2]: c2
    1 1    \    # p[3]: θ0
    1 1    \    # p[4]: θ1
    1 1         # p[5]: θ2
```

**2. `scaling_form = 1`（带修正形式）**

⚠️ **关键差异**：`p[2]` 和 `p[3]` 的位置互换，theta 扩展到 5 个

```
Params[0] = Tc      ← Tc_init
Params[1] = c1      ← c1_init (1/ν)
Params[2] = c3      ← c3_init (ω, 修正指数)  ← 注意：这里是 c3！
Params[3] = c2      ← c2_init               ← 注意：这里是 c2！
Params[4] = θ0      ← theta0_init
Params[5] = θ1      ← theta1_init
Params[6] = θ2      ← theta2_init
Params[7] = θ3      ← theta3_init
Params[8] = θ4      ← theta4_init
```

```bash
./bsa -f 1 data.dat \
    1 0.42  \    # p[0]: Tc
    1 0.9   \    # p[1]: c1 (1/ν)
    1 0.5   \    # p[2]: c3 (ω, correction exponent)
    1 0.1   \    # p[3]: c2
    1 1     \    # p[4]: θ0
    1 1     \    # p[5]: θ1
    1 1     \    # p[6]: θ2
    1 1     \    # p[7]: θ3
    1 1          # p[8]: θ4
```

**总结**：
- 标准形式：`p[2] = c2`，3 个 theta 参数
- 带修正形式：`p[2] = c3`，`p[3] = c2`（位置互换），5 个 theta 参数
- `bsa_core.jl` 自动处理这种映射差异


## BSA 输入数据格式

BSA 输入文件必须遵循以下结构：

```
# L    X           Y           Error
  24   3.700000    0.435018    0.002429
  24   3.750000    0.464603    0.002417
  ...
  
  27   3.700000    0.425342    0.003315
  27   3.750000    0.456021    0.001053
  ...
```

- 不同 L 值之间用空行分隔
- **Y（观测量）**：必须有误差棒
- **X（控制参数）**：无误差棒（假设 X 被精确控制）


## 核心设计：Data Flow（数据处理流）

`bsa_core.jl` 的核心是实现 BSA 工具的**完整数据流封装**，包括输入（stdin）和输出（stdout）的双向桥接。

### 数据流全景

```
┌─────────────────────────────────────────────────────────────────┐
│  Input (stdin): Julia → BSA                                     │
├─────────────────────────────────────────────────────────────────┤
│  BSAParameters{Tc_init, c1_init, c2_init, c3_init, ...}        │
│         ↓ build_parameter_segment()                             │
│  BSA command line: {p[0], p[1], p[2], ...}                     │
│         ↓ (处理 scaling_form 差异)                               │
│  BSA executable runs                                            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  Output (stdout): BSA → Julia → Physics                         │
├─────────────────────────────────────────────────────────────────┤
│  Layer 1: *.op file                                             │
│         ↓ parse_bsa_output()                                    │
│  Layer 2: metadata {"p0", "p1", "p2", ...}  ← 歧义               │
│         ↓ extract_parameter_dict(metadata, scaling_form)        │
│  Layer 3: parameter_dict {"Tc", "c1", "c2", "c3", ...}  ← 无歧义 │
│         ↓ extract_physical_params(..., critical_param, eta_type)│
│  Layer 4: physical_params {"Uc", "nu", "eta_psi", ...}  ← 物理量 │
│         ↓ print_summary() (可选)                                 │
│  Layer 5: 格式化输出 [Parameters] + [Physical Quantities] + [Fit] │
└─────────────────────────────────────────────────────────────────┘
```

**关键设计点**：

| 转换步骤 | 输入 | 输出 | 核心任务 |
|---------|------|------|---------|
| `parse_bsa_output()` | `*.op` 文件 | `metadata{"p0", "p1", ...}` | 解析 BSA 原始输出 |
| `extract_parameter_dict()` | `metadata` | `{"Tc", "c1", "c2", ...}` | 消除 p[2] 歧义 ⚠️ |
| `extract_physical_params()` | `parameter_dict` | `{"Uc", "nu", "eta_psi", ...}` | 物理解释 |
| `print_summary()` | `metadata` | 格式化输出 | 三层混合展示 |

⚠️ **核心创新**：`extract_parameter_dict()` 统一处理 `scaling_form` 差异
- `scaling_form=0`: p[2]=c2
- `scaling_form=1`: p[2]=c3, p[3]=c2 (位置互换)
- `scaling_form` 自动从 `metadata["form"]` 提取，无需手动传递

### 数据层级的使用场景

不同场景下使用不同层级的数据：

| 使用场景 | 使用的数据层级 | 原因 |
|---------|--------------|------|
| **与 BSA 交互** | `metadata` | BSA 工具的原生接口 |
| **内部参数传递** | `params` (无歧义参数) | 消除了 `scaling_form` 依赖，与 `BSAParameters` 一致,便于模块间传递 |
| **给人看的输出** | `phys` (物理量) | 用户友好的命名 (Uc/Tc/Jc, ν, η) |
| **绘图、分析** | `phys` (物理量) | 物理含义明确，适合展示和报告 |

**设计哲学**：
- 📦 `params` 是内部中间层，确保参数命名无歧义
- 📊 `phys` 是外部展示层，物理含义清晰，适合给用户看

### 模块职责分离

| 模块 | 核心职责 | 关键函数 |
|------|---------|---------|
| **`bsa_core.jl`** | 单次 BSA 调用 + 参数解释 | `parse → extract_parameter_dict → extract_physical_params` |
| **`bsa_bootstrap.jl`** | Bootstrap 统计 + 批量复用 | 调用 `BSACore.extract_physical_params()` 处理 bootstrap 样本，并使用 `DataProcessforDQMC.statistics` 生成 `phys_fmt`（预格式化物理量） |

✅ 物理转换逻辑只在 `bsa_core.jl` 实现一次，`bsa_bootstrap.jl` 通过组合复用  
✅ 所有下游模块（绘图、分析）共享同一套参数映射和物理转换

### 快速上手

```julia
# 完整数据流（单次 BSA 结果）
metadata, _ = parse_bsa_output("result.op")
params = extract_parameter_dict(metadata)
phys = extract_physical_params(params, critical_param_name="Uc", eta_type=:eta_psi)
print_summary(metadata, critical_param_name="Uc", eta_type=:eta_psi)  # 一步到位
```

---

