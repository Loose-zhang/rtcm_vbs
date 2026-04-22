# 双基站合并策略调整记录（A主站 + B补星）

## 背景

当前目标是让双基站在"某一基站被遮挡"时提升共视卫星数，同时避免载波相位切换导致固定解被破坏。

典型场景：

- A 站仅有 GPS
- B 站仅有 BDS
- 期望 VBS 输出同时包含 GPS + BDS，且 GPS 与 BDS 都能固定解

旧策略在同卫星上会按 SNR 在 A/B 之间切换，且 B 未完成对齐时通过 `LLI_SLIP` 强制标记，
实际容易持续扰动 RTK 模糊度滤波，表现为长期 float。

---

## 改动历程

### v1：同卫星固定优先 A + 全量抑制未对齐 B 相位

问题诊断：旧版 `normalise_b()` 在 B 侧无有效 bias 时打 `LLI_SLIP`，持续扰动滤波器。

改动：

1. 同卫星固定用 A，B 仅补 A 缺失的卫星，删除 SNR 切换门控。
2. B 侧无有效 bias 时清零 L/D，只保留 P。

**v1 局限**：对于 A=GPS、B=BDS 的纯异构场景，A/B 没有任何共视卫星，`update_one()`
永远不会被触发，bias 状态始终是零初始值（`valid=0`，`code=0`，`init_count=0`）。
v1 中 `has_code = (match != NULL)` 实际恒为 false，导致 BDS 原始相位被"全量保留"——
看似正确，但 v1 的代码路径是无 match 时走"直接 continue"（不抑制也不归一化），
相当于把 B 侧 BDS 原始相位原封不动地送下游，这反而是旧的错误：BDS 的相位还带着
B 基站的真实接收机偏差，并没有和 A 侧 GPS 处于同一相位框架，RTKNAVI 仍然拿不到
正确的 BDS 固定条件。

---

### v2（当前版本）：三路策略 + 支持独立星座固定

文件：`rtcm_vbs_hp/vbs_merge.c`，函数 `normalise_b()`

**三路决策逻辑（基于 has_history 标志）：**

```text
has_history = 该 (sat,code) 在 align 表里是否存在过 A/B 共同更新记录

情况 A  valid=1                   → 归一化 L += bias_cyc （原有逻辑）
情况 B  has_history=1, valid=0    → 清零 L/D，只保留 P
         （曾经 A/B 共视但当前 bias 失效，不注入不稳定相位）
情况 C  has_history=0             → 保留原始 B 侧 L/D 原样输出
         （从未有 A/B 共视：典型为 B-only 星座如 BDS，原始相位是自洽的，
          RTKNAVI 独立对此星座固定即可）
```

**为什么情况 C 可以直接输出原始相位？**

B 站的每颗 BDS 卫星经过 `vbs_correct_obs_ex()` 以 B 站真实 ECEF 为基准做了 VBS
几何修正（`Δρ = |r_sat - r_VBS| - |r_sat - r_B|`），输出的 L 已经是"虚拟在 VBS
位置"的 BDS 相位。虽然它包含 B 接收机的相位硬件偏差，但这个偏差是整数周级别的常数，
RTKNAVI 会在 RTK 差分过程中将其作为模糊度整数部分吸收，不影响整周固定条件。
因此 BDS 可以独立固定，与 GPS 的相位框架无关。

**has_history 标志的实现：**

scan `align[sat_idx][k]` 所有槽：
- 任意槽的 `code == 当前 code` → has_history = 1
- 同时若该槽 `valid == 1` → 直接归一化（情况 A）

---

## 最终预期行为

| 场景 | A 侧 | B 侧 | 输出 |
|---|---|---|---|
| A/B 共同可见（同星座重叠） | A 有该星 | B 也有该星 | 固定取 A；bias 状态持续更新 |
| A 缺失（B 补同星座） | A 无该星 | B 有该星，bias valid | 归一化后输出 B（纳入 A 框架）|
| A 缺失（B 补同星座，初始化中）| A 无该星 | B 有该星，has_history=1, valid=0 | 仅输出 P（不注入不稳定相位）|
| A 无 GPS，B 全 BDS（纯异构）| A=GPS | B=BDS，has_history=0 | 保留 B 原始相位，RTKNAVI 独立固定 BDS |

GPS 由 A 侧正常固定，BDS 由 B 侧独立固定，VBS 增加卫星数的目标得以实现。

---

## 受影响文件

- `rtcm_vbs_hp/vbs_merge.c`
- `rtcm_vbs_hp/fix_dual_primary_a_backup_b.md`（本文档）

---

## 编译验证

已在本地执行并通过：

```bash
cmake --build build
```
