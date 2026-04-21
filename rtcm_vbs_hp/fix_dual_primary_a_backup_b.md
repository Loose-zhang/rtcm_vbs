# 双基站合并策略调整记录（A主站 + B补星）

## 背景

当前目标是让双基站在“某一基站被遮挡”时提升共视卫星数，同时避免载波相位切换导致固定解被破坏。

典型场景：

- A 站仅有 GPS
- B 站仅有 BDS
- 期望 VBS 输出同时包含 GPS + BDS

旧策略在同卫星上会按 SNR 在 A/B 之间切换，且 B 未完成对齐时通过 `LLI_SLIP` 强制标记，实际容易持续扰动 RTK 模糊度滤波，表现为长期 float。

---

## 本次改动

### 1) 同卫星固定优先 A（不再 A/B 来回切换）

文件：`rtcm_vbs_hp/vbs_merge.c`

- 在 `merger_align_and_merge()` 的 A/B 同时可见分支中，改为：
  - 同一颗卫星两侧都存在时，始终输出 A 侧观测。
  - 仅当卫星在 A 侧缺失时，才使用 B 侧补星。

效果：

- 保持相位参考一致，避免同卫星在 A/B 之间抖动切换。
- 满足“B 用于补充 A 缺失卫星”的设计目标。

### 2) B 补星未对齐时抑制相位（保留伪距）

文件：`rtcm_vbs_hp/vbs_merge.c`

- 在 `normalise_b()` 中，针对 B 侧信号：
  - 若已有有效 B→A 对齐偏差：继续做相位归一化（`L += bias_cyc`）。
  - 若无有效偏差：将该信号 `L`、`D` 清零，仅保留 `P`。
  - 不再对这类信号追加 `LLI_SLIP`，避免长期 slip 触发影响滤波稳定。

效果：

- B 侧未对齐相位不会注入到输出观测中。
- 仍可使用 B 侧伪距扩展星座和卫星数量。

### 3) 清理旧切换门控逻辑

文件：`rtcm_vbs_hp/vbs_merge.c`

- 删除不再使用的 `b_fully_aligned()` 路径。
- `n_chose_b` 仅统计 B 作为“补星输出”的卫星。

---

## 预期行为

1. A/B 同时有同卫星：使用 A。  
2. A 无、B 有同卫星：使用 B；若未对齐则输出伪距-only。  
3. 对于 A=GPS、B=BDS 场景：VBS 输出应形成 GPS+BDS 的互补组合。

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
