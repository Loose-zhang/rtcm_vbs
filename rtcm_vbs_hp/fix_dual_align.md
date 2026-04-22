# 双基站载波相位对齐修复记录

## 背景

单基站 VBS 偏移 `3000 3000 30 --trop` 在 RTKNAVI 中可以稳定固定解，说明几何/对流层修正链路正常。开启 `--dual` 后，即使基站使用中国移动 CORS，也始终拿不到固定解。

根因：旧版双基站流水线先在 `vbs_merge.c` 中按 SNR 在 A/B 之间逐卫星混选原始观测，再统一调用 `vbs_correct_obs_ex` 做 VBS 修正。几何距离虽被搬到同一个 VBS 坐标，但 A、B 两台真实接收机的载波相位各自带有独立的接收机相位偏差和模糊度。每次某颗卫星从一侧切换到另一侧，输出的 `L` 会出现几百周量级的台阶，但 `LLI` 不带任何标记，RTKNAVI 把它当作正常观测，模糊度滤波器持续被破坏，固定永远不收敛。

修复目标：把 A、B 两站观测分别改正到同一个 VBS 几何框架，再估计并补偿 B 相对于 A 的载波相位偏差，使最终输出的虚拟基站载波序列在 A/B 切换时保持连续。

## 修改总览

| 文件 | 类型 | 说明 |
|------|------|------|
| `rtcm_vbs_hp/vbs_merge.h` | 重写 | 新增 `align_state_t`，扩展 `merger_t`，替换合并接口 |
| `rtcm_vbs_hp/vbs_merge.c` | 重写 | 实现 per-(sat, code) 偏差估计、B 端载波归一化、带连续性约束的合并策略 |
| `rtcm_vbs_hp/main.c` | 局部 | 重排 dual 分支为「分别 VBS 改正 → 对齐 → 合并 → 输出」，扩展 `[stat]` 字段 |

## 1. 数据流改造（main.c dual 分支）

旧流程：

```
A/B 原始 obs → merger_poll(SNR 选边) → vbs_correct_obs_ex → emit_msm
```

新流程：

```
A/B 原始 obs
   → merger_poll_pair (返回两侧原始 obs，不做选边)
   → 分别 vbs_correct_obs_ex(side=0/1)
   → merger_align_and_merge (估计 B→A bias，归一化 B，按连续性约束合并)
   → emit_msm
```

`--iono` 仍在 A、B 同时存在时构建一份 `sat_extra_m`，两侧统一施加，避免 net 修正在两侧不一致。

## 2. 对齐状态结构

`merger_t` 内新增：

```c
align_state_t align[MAXSAT][NFREQ+NEXOBS];
```

每个槽位记录：

- `valid`：B→A 偏差是否可用于归一化
- `code`：本槽位绑定的 `CODE_???` 信号码
- `init_count` / `init_sum`：初始化阶段的样本计数与累加
- `bias_cyc`：valid 后的偏差，单位 cycles
- `last_diff` / `last_pair_t`：调试用与超时判定

key 使用 `(sat-1, slot)` 并以 `code` 字段二次校验，确保 A/B 两侧信号槽顺序不同时不会错配。

## 3. 偏差估计与失效规则

每次 A、B 同历元同卫星同 `code` 都可见时计算：

```
diff_cyc = L_A_vbs - L_B_vbs
```

规则：

| 阶段 | 条件 | 动作 |
|------|------|------|
| 初始化 | 连续 3 次 `|diff - mean| ≤ 0.15 cycle` | 锁定 `bias_cyc = mean`，置 `valid=1` |
| 跟踪 | `|diff - bias| ≤ 0.15 cycle` | IIR 更新 `bias = bias + 0.05*(diff-bias)` |
| reset | A 或 B `LLI & LLI_SLIP` | 清空 → `n_align_reset_slip++` |
| reset | 槽位 `code` 变化 | 清空 → `n_align_reset_code++` |
| reset | 配对中断 > 5 s | 清空 → `n_align_reset_gap++` |
| reset | `|diff - bias| > 0.75 cycle` | 清空 → `n_align_reset_resid++` |

reset 后该槽位重新进入初始化，期间任何输出 B 的载波都会被强制 `LLI |= LLI_SLIP`。

## 4. B 端载波归一化

`normalise_b()` 对 B 中每条非零载波：

- 在该卫星的对齐表中按 `code` 找匹配槽位
- 若 `valid`：`L[j] += bias_cyc`，统计 `n_b_norm_applied++`
- 若无效：`LLI[j] |= LLI_SLIP`，统计 `n_b_unaligned_lli++`

注意：B 的伪距 `P[j]` **不**叠加 B→A 载波偏差，只保留几何/对流层/网络改正；伪距和载波的接收机硬件偏差量级不同、性质不同，不该共用一个标量。

## 5. 合并策略（带连续性约束）

`merger_align_and_merge` 在 A/B 都齐时按卫星：

- 仅 A 有：直接用 A
- 仅 B 有：归一化后输出，未对齐信号置 LLI
- 两侧都有：
  - 计算两侧平均 SNR
  - SNR 差额需超过 2000（即 2 dB-Hz）才考虑切到 B（hysteresis 防抖动）
  - 即使 SNR 倾向 B，也要满足「B 的所有活跃载波都已 aligned」才真正切，否则保留 A 并 `n_switch_blocked++`

仅 A 或仅 B 历元（A-only / B-only）走 single-side 路径，B-only 仍走归一化逻辑。

## 6. LLI 策略

- 始终保留源 RTCM 的 LLI 位（包括 `LLI_HALFC` 等）
- 仅在以下情形额外置 `LLI_SLIP`：
  - 输出 B 的某条载波但该 (sat, code) 无有效 bias
  - bias 状态本历元发生 reset
- A↔B 切换若 bias 有效且残差正常，**不**额外置 slip，这正是修复后的关键预期

## 7. 统计字段扩展

dual 模式 `[stat]` 行新增字段：

```
align init=N reset(slip=N code=N gap=N resid=N)
Bnorm=N Bunal=N switch_blk=N
```

| 字段 | 含义 | 期望趋势 |
|------|------|----------|
| `align init` | 成功完成初始化的 (sat, code) 槽位累计次数 | 启动阶段快速增长，稳定后偶发 |
| `reset_slip/code/gap/resid` | 各类 reset 累计 | 短基线、稳定信号下应保持低位 |
| `Bnorm` | 实际归一化后输出 B 载波的次数 | 至少在 A/B 都活跃时持续增长 |
| `Bunal` | 因无对齐而被强制置 LLI 的 B 载波次数 | 主要出现在启动初期或 reset 后 |
| `switch_blk` | SNR 倾向 B 但因未对齐被回退到 A 的次数 | 启动后应趋于 0 |

## 8. 验证方法

1. 构建：`cmake --build build_mingw`，产物 `build_mingw/rtcm_vbs_hp.exe`
2. 单站回归：`--offset 3000 3000 30 --trop`，RTKNAVI 应继续可固定
3. 双站基础：`--dual --offset 3000 3000 30 --trop`
   - 启动日志应同时出现 `[vbs] base A ECEF` 和 `[vbs] base B ECEF`
   - `align init` 在前几十秒内快速增长
   - `Bnorm` 持续增长，`switch_blk` 趋零
   - RTKNAVI 进入固定或 fixed ratio 明显改善
4. 双站 + iono：`--dual --trop --iono`
   - `reset_resid` 不应飙升
5. 压力测试：手动断开 A 或 B 的 NTRIP 一段时间再恢复
   - 期间 `Bunal` 短暂增长，恢复后 `align init` 再次推进
   - RTKNAVI 在恢复后能重新固定

## 文件清单

- `rtcm_vbs_hp/vbs_merge.h`
- `rtcm_vbs_hp/vbs_merge.c`
- `rtcm_vbs_hp/main.c`
