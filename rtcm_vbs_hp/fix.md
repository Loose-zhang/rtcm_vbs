# net_interp.c 电离层改正修复记录

## Fix 1：对流层重复改正问题

### 问题描述

`--iono` 与 `--trop` 同时开启时，`net_interp_build_dual` 计算的残差

```
dres = (P_B - P_A) - (rho_B - rho_A)
```

中已包含 A、B 两站之间的对流层差分量 `trop(B) - trop(A)`。而 `vbs_core.c` 的 `vbs_correct_obs_ex` 在施加 `dnet`（即本模块输出的 `sat_extra_m`）的同时，又独立计算并叠加了 `dtrop = trop(VBS) - trop(base)`，导致对流层误差被改正两次。

### 修复方法

在 `dres` 计算后、进入插值前，若 `vbs->apply_trop == 1`，调用 `tropmodel` 计算 A/B 两站对该卫星的对流层延迟，并从 `dres` 中减去其差值，使残差只含电离层分量：

```c
if (vbs->apply_trop) {
    double ta = tropmodel(oa->time, pos_a, azel_a, vbs->humi);
    double tb = tropmodel(ob->time, pos_b, azel_b, vbs->humi);
    dres -= (tb - ta);
}
```

`azel_a`/`azel_b` 由 `satazel()` 在同一循环体内计算（Fix 3 中也需要），复用无额外开销。`vbs->humi` 直接读取 `vbs_ctx_t` 中已有的字段（默认 0.7）。

---

## Fix 2：一维线性投影改为 IDW 二维插值

### 问题描述

原 `baseline_fraction` 函数将 VBS 位置投影到 A→B 连线方向，得到标量系数 `alpha`，再用 `corr = alpha * dres` 插值。当 VBS 不在 A-B 连线上时，垂直于基线方向的大气误差完全未被内插，且 `alpha` 被人工限幅到 `[-1, 2]`，外推区域精度难以保证。

### 修复方法

删除 `baseline_fraction`，新增 `idw_weights` 函数，基于 VBS 到 A、B 两站的三维欧氏距离计算反距离加权系数 `wa`/`wb`：

```c
static void idw_weights(const double *vbs_ecef,
                        const double *a_ecef, const double *b_ecef,
                        double *wa, double *wb)
{
    // 计算 VBS 到 A、B 的距离 da、db
    // wa = (1/da) / (1/da + 1/db)
    // wb = (1/db) / (1/da + 1/db)
}
```

`wa + wb = 1`，自然限定在 `[0, 1]`，无需手动限幅。内插公式改为 `corr = wb * dres`（以 A 站为参考零点，B 站残差按距离权重 wb 分配到 VBS），同时覆盖 VBS 偏离基线垂直方向的场景。

---

## Fix 3：仰角倾斜因子归一化

### 问题描述

电离层延迟是沿信号路径的斜路径值（slant），其大小与卫星仰角相关（低仰角路径更长，延迟更大）。A、B 两站对同一颗卫星的仰角并不相同，直接对斜路径残差 `dres` 做线性/IDW 插值，等同于隐含地假设各站仰角相同，对低仰角卫星或长基线场景误差可达数十厘米。

### 修复方法

利用 RTKLIB 的 `ionmapf()` 计算各站的电离层倾斜映射因子，将 `dres` 折算到天顶（VTEC）域后再插值，最终用 VBS 侧的映射因子还原为斜路径延迟：

```c
double mf_a = ionmapf(pos_a, azel_a);
double mf_b = ionmapf(pos_b, azel_b);
double mf_mean = 0.5 * (mf_a + mf_b);

double dres_vtec = dres / mf_mean;          // 折算到天顶域
double mf_v = wa * mf_a + wb * mf_b;        // VBS 侧映射因子（IDW 近似）
double corr = wb * dres_vtec * mf_v;        // 还原为 VBS 斜路径
```

**退化条件**：若任一站的映射因子超出合理范围（`< 1e-3` 或 `>= 5`，对应仰角约低于 11.5°），跳过天顶化，直接使用 `corr = wb * dres`，防止低仰角噪声被放大。

当 `mf_a ≈ mf_b`（高仰角卫星）时，`mf_mean ≈ mf_v`，公式自动退化为 `wb * dres`，与修复前行为一致。

---

## 修改文件

| 文件 | 改动 |
|------|------|
| `rtcm_vbs_hp/net_interp.c` | 删除 `baseline_fraction`；新增 `idw_weights`；重写 `net_interp_build_dual` 内循环 |
| `rtcm_vbs_hp/net_interp.h` | 无改动（接口不变） |
| `rtcm_vbs_hp/vbs_core.c` | 无改动（`apply_trop` 字段直接读取） |

## 使用建议

| 启动参数组合 | 说明 |
|---|---|
| `--dual --iono` | 推荐：电离层插值，对流层不改正（对流层差由 dres 隐含吸收） |
| `--dual --iono --trop` | 修复后可安全使用：对流层在两处均被正确且不重复地施加 |
| 基线 < 50 km | IDW + 倾斜因子改正效果最佳 |
| 基线 > 50 km | 建议配合 `--trop` 同时开启，可进一步消除对流层差残余 |
