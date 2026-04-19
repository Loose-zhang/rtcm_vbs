# rtcm_vbs_hp — 高精度 RTCM3 虚拟基站

基于 **RTKLIB**（本仓库 `../src` + `../include`）重构的 C 实现，替代原 Python 版 `rtcm_vbs*.py`。完整解码 → 物理量修正 → 重新编码，彻底摆脱原版基于字节原位修改带来的量化与精度天花板。

---

## 目录

1. [特性概览](#1-特性概览)
2. [为什么精度更高](#2-为什么精度更高)
3. [工作原理](#3-工作原理)
4. [依赖与编译](#4-依赖与编译)
5. [运行模式](#5-运行模式)
   - [5.1 单基站](#51-单基站)
   - [5.2 独立星历流](#52-独立星历流)
   - [5.3 独立 SSR 流](#53-独立-ssr-流ppp--ppp-rtk-增强)
   - [5.4 双基站 CNR 合并](#54-双基站-cnr-合并模式)
6. [命令行参数](#6-命令行参数)
7. [消息处理规则](#7-消息处理规则)
8. [日志与统计](#8-日志与统计)
9. [与 Python 版精度对比](#9-与-python-版精度对比)
10. [常见问题](#10-常见问题)
11. [当前限制](#11-当前限制)
12. [代码结构](#12-代码结构)

---

## 1. 特性概览

| 能力 | 说明 |
|---|---|
| 全精度 VBS 修正 | `satposs()`（τ 迭代 + 相对论 + 钟差 + BDS GEO）+ `geodist()`（Sagnac）|
| 全星座 | GPS / GLONASS / Galileo / BeiDou / QZSS / IRNSS |
| 差分对流层 | Saastamoinen 差分（可选 `--trop`），天向偏移大时尤为重要 |
| 载波相位一致 | `L += Δρ/λ`，λ 由 `sat2freq()` 按信号逐一精确计算，保持整周模糊连续 |
| 相位重编码 | `gen_rtcm3()`，输出真 RTCM3 MSM4 / MSM7，无 fine_pr / fine_ph 量化溢出 |
| 双基站 CNR 合并 | `--dual`，每颗卫星从 SNR 更高的一侧取观测 |
| SSR 原样透传 | 1057–1068 / 1240–1270 / 11–14 等直接帧级转发，覆盖 HR-clock、phase-bias |
| 独立输入源 | 观测流 / 星历流 / SSR 流可在三个不同的 TCP 端口 |
| 多客户端广播 | 支持 N 个下游接收端同时连接同一输出端口 |

---

## 2. 为什么精度更高

原 Python 实现通过**字节原位修改** MSM 帧中的 `rough_int / rough_mod / fine_pr / fine_ph` 字段实现 VBS，存在以下不可突破的天花板；本 C 版全部解决。

| 问题 | 原 Python 版 | `rtcm_vbs_hp` |
|---|---|---|
| 卫星位置计算 | 手写简化 Kepler，无 τ 迭代 | `satposs()` 含 `τ = ρ/c` 迭代 + 相对论 + 钟差 + BDS GEO 旋转 |
| 几何距离 | 欧氏距离，无 Sagnac | `geodist()` 内置 Sagnac（消除 ~30 m 系统偏差）|
| TOW 计算 | 早期版本未取模一周（已用 fix 补） | 直接使用 `gtime_t`，无任何时间跳变 |
| MSM 字段量化 | `rough_mod` 步长 293 m，大偏移溢出 | 完整解码→物理量修正→`gen_rtcm3()` 重新编码，无溢出 |
| 对流层差分 | 未处理（VBS 高差时显著误差） | Saastamoinen 差分（可选） |
| 支持星座 | 仅 GPS / GAL / BDS | GPS / GLO / GAL / BDS / QZS / IRN |
| 相位连续性 | fine_ph 直接偏移，易破坏整周模糊 | `L_new = L_old + Δρ/λ`，全精度浮点累计 |
| SSR 消息 | 未处理 | 帧级透传（含 HR-clock、phase-bias） |
| 双基站合并 | 有（`dual` 模式） | 有（CNR 合并 + 历元对齐） |

---

## 3. 工作原理

对每个历元、每颗卫星：

```
r_sat_i   = satposs(t, obs_i, nav)        # 含 τ 迭代、钟差、相对论、BDS GEO
r_base    = 真实基站 ECEF（来自 1005/1006）
r_VBS     = f(r_base, N/E/U 偏移)   或   pos2ecef(LAT, LON, ALT)

Δρ_i      = geodist(r_sat_i, r_VBS)  -  geodist(r_sat_i, r_base)
Δtrop_i   = tropmodel(VBS, el) - tropmodel(base, el)   [可选, --trop]

P_new[j]  = P_old[j]  + Δρ  + Δtrop
L_new[j]  = L_old[j]  + (Δρ + Δtrop) / λ_j            # λ_j = c / sat2freq()
D_new[j]  = D_old[j]                                  # VBS 速度与基站相同
```

然后 `gen_rtcm3(MSM4 / MSM7)` 把修正后的 `obsd_t` 数组重编码为标准 RTCM3 消息，CRC-24Q 自动重算。

### 消息分发（`input_rtcm3` 返回值）

| 返回 | 含义 | 动作 |
|---:|---|---|
| `1` | 观测帧（MSM）| 单基站直接修正+编码；双基站先进合并器 |
| `2` | 广播星历（1019/1020/1042/1045/1046）| 镜像到 `rout.nav`；可选重新编码转发 |
| `5` | 测站参数（1005/1006）| 学习基站 ECEF → 计算 VBS ECEF → 改写并发送 |
| `10` | SSR | **字节级透传**（不解码再编码）|
| 其他 | 其他消息 | 丢弃（目前未用到）|

---

## 4. 依赖与编译

### 系统

- macOS / Linux（已验证 arm64 macOS）；Windows 下 MinGW + winsock2
- CMake ≥ 3.10
- `pthread`、`libm`

### 构建

```bash
cd rtcm_vbs_hp
cmake -B build
cmake --build build -j
```

产物：`build/rtcm_vbs_hp`。

CMake 只拉取 RTKLIB 中本程序需要的最小源集（`rtcm*.c`、`ephemeris.c`、`rtkcmn.c`、`trace.c`、`tides.c`、`ionex.c`、`sbas.c`、`preceph.c`、`options.c`、`datum.c`、`geoid.c`、`B2b.c` 以及 `f2c/*.c`），无需引入 GUI/接收机解码器等模块，编译约 10 秒。

---

## 5. 运行模式

约定：`str2str` 是 RTKLIB 自带工具。以下示例假设 NTRIP 源是 `ntrip.data.gnss.ga.gov.au:2101`。

### 5.1 单基站

```bash
# 主观测流 → :50001
str2str -in ntrip://user:pass@ntrip.data.gnss.ga.gov.au:2101/OBS_MOUNT \
        -out tcpsvr://:50001 &

# 生成一个在真实基站 北 1000 m、东 500 m、上 50 m 的虚拟基站
./build/rtcm_vbs_hp \
    --src-port 50001 --out-port 50002 \
    --offset  1000  500  50  --trop

# 接收端（rtkrcv / RTKLIB / u-blox RTKLIB 分支等）指向 tcp://127.0.0.1:50002
```

绝对经纬度模式：

```bash
./build/rtcm_vbs_hp \
    --src-port 50001 --out-port 50002 \
    --vbs  22.937764  113.206806  120.253  --trop
```

### 5.2 独立星历流

主挂载点不含 1019/1020/1042/1045/1046 时：

```bash
str2str -in ntrip://user:pass@ntrip.data.gnss.ga.gov.au:2101/BCEP00BKG0 \
        -out tcpsvr://:50010 &

./build/rtcm_vbs_hp \
    --src-port 50001 --eph-port 50010 --out-port 50002 \
    --offset 1000 500 50 --trop
```

### 5.3 独立 SSR 流（PPP / PPP-RTK 增强）

SSR 修正描述的是**卫星本身**的轨道/钟差/偏差误差，与测站位置无关，因此本程序对所有 SSR 帧**原样透传**（不解码再编码）——既避免 RTKLIB SSR 编码器覆盖不全（缺 HR-clock 1062/1068、phase-bias 11–14）的问题，又把延迟降到最低。

```bash
# SSR 源（如 GA 的某个 SSR 挂载点）
str2str -in ntrip://user:pass@ntrip.data.gnss.ga.gov.au:2101/SSR_MOUNT \
        -out tcpsvr://:50020 &

./build/rtcm_vbs_hp \
    --src-port 50001 --eph-port 50010 --ssr-port 50020 \
    --out-port 50002 --offset 1000 500 50 --trop
```

> SSR 也可以直接嵌在主 `--src-port` 流里——程序对所有通道都会识别并透传 SSR 帧。
> 关闭 SSR 转发：`--no-ssr`。

### 5.4 双基站 CNR 合并模式

两个物理基站合成一个 VBS。每颗卫星从平均 SNR 更高的一侧选取观测，**并使用该侧真实基站的 ECEF 做几何修正**——即 A 站选来的卫星按 `A → VBS` 修正，B 站选来的卫星按 `B → VBS` 修正。这样即使 A/B 基线达几公里～几十公里，合并后的观测也保持与 VBS 一致的几何，不会引入"错配 base"导致的伪距偏差。

下游只看到一条 VBS 流，因此外送的 1005/1006 帧只发一个坐标（VBS 的 ARP）；A/B 各自的 1005/1006 会被分别学习到 `vbs_ctx_t.base_ecef[0/1]` 供修正使用，但不直接广播。VBS 的锚点规则：偏移模式下以 A 站第一次学习到的位置为基准套 NED；绝对模式下直接使用 `--vbs` 给定的 LLH。

```bash
str2str -in ntrip://.../BASE_A -out tcpsvr://:50001 &
str2str -in ntrip://.../BASE_B -out tcpsvr://:50003 &

./build/rtcm_vbs_hp \
    --dual --src-port 50001 --b-port 50003 \
    --eph-port 50010 --ssr-port 50020 \
    --out-port 50002 \
    --vbs  22.937764  113.206806  120.253 \
    --window 40  --trop
```

**合并策略**：

| 情况 | 处理 |
|---|---|
| 同一历元仅 A 有某卫星 | 使用 A 的观测 |
| 同一历元仅 B 有某卫星 | 使用 B 的观测 |
| A、B 都有该卫星 | 取平均 SNR 更高的一侧 |
| A、B 历元时间差 > 10 ms | 先输出较旧的那一帧，不丢历元 |
| 仅一侧到达、等满 `--window MS` | 单侧输出 |

---

## 6. 命令行参数

| 参数 | 说明 |
|---|---|
| `--src-host HOST`                 | 主源 RTCM TCP 主机（默认 `127.0.0.1`）|
| `--src-port PORT`                 | 主源端口（默认 `50001`）|
| `--dual`                          | 启用双基站 CNR 合并模式 |
| `--b-host HOST` / `--b-port PORT` | 双基站模式下 B 站地址 |
| `--window MS`                     | 双基站历元对齐窗口（默认 `40` ms）|
| `--eph-host / --eph-port`         | 可选独立星历流 |
| `--ssr-host / --ssr-port`         | 可选独立 SSR 流（原样透传到输出）|
| `--out-port PORT`                 | VBS 输出端口（默认 `50002`）|
| `--offset N E U`                  | NED 偏移模式（米）|
| `--vbs LAT LON ALT`               | 绝对坐标模式（度, 度, 米）|
| `--trop`                          | 启用差分 Saastamoinen 对流层修正 |
| `--msm N`                         | 输出 MSM 子类型：`4` 或 `7`（默认 `7`）|
| `--no-eph`                        | 不再次发送星历消息 |
| `--no-ssr`                        | 不转发 SSR 消息 |
| `--trace FILE` / `--level N`      | 打开 RTKLIB trace 日志 |

---

## 7. 消息处理规则

| RTCM3 消息 | 处理 |
|---|---|
| **1005 / 1006** | 学习真实 ARP → 计算 VBS ARP → 以 **1006** 重新编码后发送 |
| **1019 / 1020 / 1042 / 1044 / 1045 / 1046** | 存入 `rout.nav`；默认重新编码转发（`--no-eph` 关闭）|
| **MSM4 / MSM7 (全星座)** | 解码 → VBS 修正 → `gen_rtcm3(MSM4 或 MSM7)` 按星座拆分编码后发送 |
| **MSM1 / MSM2 / MSM3 / MSM5 / MSM6** | 被 RTKLIB 解码为 obs；输出端统一以 MSM4 或 MSM7 重新编码 |
| **1057–1068 / 1240–1270 / 11–14 (SSR)** | **字节级原样透传**（不解码再编码）|
| **1007 / 1008 / 1033 (天线信息)** | 暂未同步重写（不影响主流 rover）|
| **其它** | 丢弃 |

---

## 8. 日志与统计

启动时打印配置摘要：

```
rtcm_vbs_hp: A=127.0.0.1:50001 -> out=:50002 | mode=NED-offset | trop=1 | MSM=7 | SSR fwd=1
        eph stream 127.0.0.1:50010
        ssr stream 127.0.0.1:50020
        offset N=1000.000 E=500.000 U=50.000 m
[tcp] VBS output listening on 0.0.0.0:50002
[src] connected 127.0.0.1:50001
[eph] connected 127.0.0.1:50010
[ssr] connected 127.0.0.1:50020
[vbs] base ECEF = (-2324876.317, 5387635.496, 2491675.078)  LLH = (23.14638978, 113.34119827, 103.504)
[vbs] VBS  ECEF = (-2325197.872, 5387118.693, 2492614.235)  LLH = (23.15541915, 113.34608074, 153.602)
```

每 15 秒打印一次统计：

**单基站**
```
[stat] 1005/6=11 MSM=104 | obs in=2993 out=2993 | sat corr=2993 no_eph=0 | SSR fwd=37 | clients=1
```

**双基站**
```
[stat] 1005/6=11 MSM=104 | obs in=2993 out=2993 | sat corr=2993 no_eph=0 \
       | epochs merged=98 A-only=2 B-only=1 \
       | CNR both_sat=780 chose_A=512 chose_B=268 \
       | SSR fwd=37 | clients=1
```

| 字段 | 含义 |
|---|---|
| `1005/6` | 已输出的 1005/1006 帧数（始终是 VBS 坐标） |
| `MSM` | 已输出的 MSM 帧数（多星座拆分后）|
| `obs in / out` | 进入修正的卫星数 / 成功输出的卫星数 |
| `sat corr` | 应用了 VBS 修正的累计卫星次数 |
| `no_eph` | 因缺星历（或对应侧 base 未学习）被丢弃的卫星次数 |
| `epochs merged` | A 和 B 都齐的历元数 |
| `epochs A-only / B-only` | 只有一侧抵达、在超时窗口后单独发出的历元数 |
| `CNR both_sat` | 同卫星两侧都有、做了 CNR 比较的累计次数 |
| `CNR chose_A / chose_B` | CNR 比较时选了 A 或 B 的次数 |
| `SSR fwd` | 已透传的 SSR 帧数 |
| `clients` | 当前连接到输出端口的下游客户端数 |

---

## 9. 与 Python 版精度对比

在 1 km NED 偏移、静态基准条件下的理论误差量级：

| 误差源 | Python 版 | `rtcm_vbs_hp` |
|---|---|---|
| 卫星位置（无 τ 迭代）         | 20–60 m      | < 1 cm       |
| Sagnac 未修正                 | ~30 m        | 0            |
| `rough_mod` 量化（293 m 步长）| 亚米–米级    | 无（重编码）|
| 对流层（VBS 高差 50 m）       | 未修正 ~2 dm | 可修正到 mm |
| **综合 Δρ 精度**             | **分米–米级** | **毫米级**   |

下游 RTK 求解精度随之提升，载波相位连续性保持，有利于整周模糊度 fix。

---

## 10. 常见问题

**Q：日志里 `no_eph` 持续增长？**
A：星历不足。检查主 `--src-port` 是否包含 1019/1020/1042/1045/1046，或为星历流单独指定 `--eph-port`。冷启动 1–2 分钟星历累计到位后即为 0。

**Q：下游 rover 刚连上就断开？**
A：试着用
```
str2str -in tcp://127.0.0.1:50002 -out /tmp/vbs.rtcm3
convbin -r rtcm3 /tmp/vbs.rtcm3 -o /tmp/vbs.obs -n /tmp/vbs.nav
grep "APPROX POSITION" /tmp/vbs.obs
```
确认输出流解码正常，且 `APPROX POSITION XYZ` 与日志里 `[vbs] VBS ECEF` 一致。

**Q：双基站下 `merge B=0` 或始终远小于 `merge A`？**
A：B 站数据迟到超过 `--window`，被强制单 A 输出。增大 `--window` 到 80–100 ms，或检查 B 站 str2str 是否健康。

**Q：上游 1005/1006 发得很慢（~7.5 s / 次）？**
A：由 NTRIP 挂载点播发策略决定，本程序不加工频率。对 RTK 求解无实质影响（位置准静态）。

**Q：如何只看 SSR 是否透传成功？**
A：观察 `[stat]` 里的 `SSR fwd` 每轮递增。若始终为 0，确认 `--ssr-port` 上游是 `tcpsvr` 方向、且确实有 SSR 消息（1057–1068 等）。

**Q：能和原 Python 版同时跑吗？**
A：可以。两者默认端口相同（50001 输入、50002 输出），要同时跑请把其中一个的 `--out-port` 改到其它端口（如 50012），以便 A/B 对比。

---

## 11. 当前限制

1. **天线描述**（1007/1008/1033）暂未同步重写；对绝大多数 rover 无影响。
2. **LLI / lock 计数**按输入原样透传；理论上 VBS 不改变失锁事件，正确。
3. **双基站**合并时外送的 1005/1006 始终是 VBS 的单一坐标；A/B 各自的 1005/1006 仅在内部分别维护为两个真实基站 ECEF，用于对每颗卫星按来源侧做 `base → VBS` 几何修正（baseline 为 km 级时仍保持精度）。
4. **SSR** 走帧级透传，不对 SSR 修正做 VBS 补偿；SSR 修正描述卫星本身误差，与测站位置无关，这样处理是正确的。
5. **大气延迟（电离层/对流层）** 在 A、B 两站本来就不同。几何差分已被正确处理；差分对流层可用 --trop 启用，但前提是 side 对应的 base LLH 正确——这点也已vbs_correct_obs 里按 side 选 base_pos_llh[side] 做 Saastamoinen 差分。
6. **电离层差分未做** 如果 A/B 相距几十公里、且不过 SSR 电离层修正，单频短基线下这部分残差需要 rover 侧 RTK 滤波去吸收。如果需要，可用 SSR VTEC（MT 1264）进一步校正；目前 SSR 是透传。

---

## 12. 代码结构

```
rtcm_vbs_hp/
├── CMakeLists.txt        # 最小化拉取 RTKLIB 源 + 本程序三个模块
├── README.md             # 本文件
├── main.c                # 事件循环 + CLI + channel 抽象 + SSR 识别
├── vbs_core.h / .c       # 核心 VBS 物理修正（位置 + 观测量）
├── vbs_merge.h / .c      # 双基站 CNR 合并 + 历元对齐
├── tcp_io.h  / .c        # TCP 客户端 + 多客户端广播服务器
└── rtklib_glue.c         # 提供 RTKLIB 未定义的 PPP_Glo 全局变量
```

### 主要 RTKLIB API 使用点

| RTKLIB 函数 | 使用目的 |
|---|---|
| `input_rtcm3()`           | 输入字节流的增量解码 |
| `gen_rtcm3()`             | 输出 MSM4/7、1005/1006、1019/1020/1042/1045/1046 编码 |
| `satposs()`               | 发射时刻卫星位置 / 钟差 / 相对论 / BDS GEO |
| `geodist()`               | 含 Sagnac 的几何距离 |
| `satazel()`               | 高度角（对流层模型用）|
| `tropmodel()`             | Saastamoinen 对流层延迟 |
| `ecef2pos() / pos2ecef()` | 大地坐标 ⇄ ECEF |
| `enu2ecef()`              | NED 偏移 → ECEF 平移向量 |
| `sat2freq()`              | 每颗卫星-每信号的精确载波频率（含 GLONASS FCN）|

---

*文档最后更新：2026-04-18*
