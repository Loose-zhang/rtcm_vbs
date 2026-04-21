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
   - [5.5 双基站远距增强（电离层/网络插值）](#55-双基站远距增强电离层网络插值)
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
| 远距电离层/网络插值 | `--iono`，基于 A/B 同卫星伪距残差减几何差的结果，沿基线向 VBS 做一维插值 |
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
| 长基线空间残差 | 无显式处理 | 新增 `--iono` 一维网络残差插值（双基站） |

---

## 3. 工作原理

对每个历元、每颗卫星：

```text
r_sat_i   = satposs(t, obs_i, nav)        # 含 τ 迭代、钟差、相对论、BDS GEO
r_base    = 真实基站 ECEF（来自 1005/1006）
r_VBS     = f(r_base, N/E/U 偏移)   或   pos2ecef(LAT, LON, ALT)

Δρ_i      = geodist(r_sat_i, r_VBS)  -  geodist(r_sat_i, r_base)
Δtrop_i   = tropmodel(VBS, el) - tropmodel(base, el)   [可选, --trop]
Δnet_i    = net_interp(A, B, sat, VBS)                 [可选, --iono, 双基站]

P_new[j]  = P_old[j]  + Δρ + Δtrop + Δnet
L_new[j]  = L_old[j]  + (Δρ + Δtrop + Δnet) / λ_j
D_new[j]  = D_old[j]
```

其中 `--iono` 的一阶段实现不是完整 VRS 网络模型，而是：

1. 找出 A/B 两站同历元共同观测到的卫星
2. 取该卫星的首个有效伪距 `P_A`、`P_B`
3. 计算站间残差去几何后的量：

```text
dres = (P_B - P_A) - (ρ_B - ρ_A)
```

4. 计算 VBS 在 A→B 基线方向上的投影比例 `α`
5. 得到该卫星的附加修正：

```text
Δnet = α * dres
```

这部分主要吸收长基线下的空间相关残差，尤其是电离层导致的差异。当前版本是**双基站、一维沿基线插值**，适合作为远距增强的第一阶段方案。

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

- macOS / Linux（已验证 arm64 macOS）；Windows 下 MSVC 或 MinGW + winsock2
- CMake ≥ 3.10
- `pthread`、`libm`

### 通用构建

```bash
cd rtcm_vbs_hp
cmake -B build
cmake --build build -j
```

产物：`build/rtcm_vbs_hp` 或 Windows 下的 `build-win/Release/rtcm_vbs_hp.exe`。

### Windows + MSVC

推荐在 **x64 Native Tools Command Prompt for VS** 或 **Developer PowerShell for VS** 中执行：

```bat
cd /d e:\RTKLIB_demo5_b34i\rtcm_vbs\rtcm_vbs_hp
cmake -S . -B build-win -G "Visual Studio 17 2022" -A x64
cmake --build build-win --config Release
```

如已安装更早版本的 Visual Studio，可把生成器替换为对应版本，例如：

```bat
-G "Visual Studio 16 2019"
```

CMake 只拉取 RTKLIB 中本程序需要的最小源集（`rtcm*.c`、`ephemeris.c`、`rtkcmn.c`、`trace.c`、`tides.c`、`ionex.c`、`sbas.c`、`preceph.c`、`options.c`、`datum.c`、`geoid.c`、`B2b.c`、`f2c/*.c`，以及本程序的 `net_interp.c`），无需引入 GUI/接收机解码器等模块。

---

## 5. 运行模式

约定：`str2str` 是 RTKLIB 自带工具。以下示例假设 NTRIP 源是 `ntrip.data.gnss.ga.gov.au:2101`。

### 5.1 单基站

```bash
str2str -in ntrip://user:pass@ntrip.data.gnss.ga.gov.au:2101/OBS_MOUNT \
        -out tcpsvr://:50001 &

./build/rtcm_vbs_hp \
    --src-port 50001 --out-port 50002 \
    --offset 1000 500 50 --trop
```

绝对经纬度模式：

```bash
./build/rtcm_vbs_hp \
    --src-port 50001 --out-port 50002 \
    --vbs 22.937764 113.206806 120.253 --trop
```

### 5.2 独立星历流

```bash
str2str -in ntrip://user:pass@ntrip.data.gnss.ga.gov.au:2101/BCEP00BKG0 \
        -out tcpsvr://:50010 &

./build/rtcm_vbs_hp \
    --src-port 50001 --eph-port 50010 --out-port 50002 \
    --offset 1000 500 50 --trop
```

### 5.3 独立 SSR 流（PPP / PPP-RTK 增强）

```bash
str2str -in ntrip://user:pass@ntrip.data.gnss.ga.gov.au:2101/SSR_MOUNT \
        -out tcpsvr://:50020 &

./build/rtcm_vbs_hp \
    --src-port 50001 --eph-port 50010 --ssr-port 50020 \
    --out-port 50002 --offset 1000 500 50 --trop
```

### 5.4 双基站 CNR 合并模式

```bash
str2str -in ntrip://.../BASE_A -out tcpsvr://:50001 &
str2str -in ntrip://.../BASE_B -out tcpsvr://:50003 &

./build/rtcm_vbs_hp \
    --dual --src-port 50001 --b-port 50003 \
    --eph-port 50010 --ssr-port 50020 \
    --out-port 50002 \
    --vbs 22.937764 113.206806 120.253 \
    --window 40 --trop
```

**合并策略**：

| 情况 | 处理 |
|---|---|
| 同一历元仅 A 有某卫星 | 使用 A 的观测 |
| 同一历元仅 B 有某卫星 | 使用 B 的观测 |
| A、B 都有该卫星 | 取平均 SNR 更高的一侧 |
| A、B 历元时间差 > 10 ms | 先输出较旧的那一帧，不丢历元 |
| 仅一侧到达、等满 `--window MS` | 单侧输出 |

### 5.5 双基站远距增强（电离层/网络插值）

在双基站模式下，额外启用：

```bash
./build/rtcm_vbs_hp \
    --dual --src-port 50001 --b-port 50003 \
    --out-port 50002 \
    --vbs 22.937764 113.206806 120.253 \
    --window 40 --trop --iono
```

适用场景：

- VBS 偏移较远（如 10–30 km）
- 原始几何 + 对流层修正后，仍存在明显空间相关残差
- 需要在现有双站结构上做第一阶段远距增强

注意：

- `--iono` 当前**必须搭配 `--dual`**
- 当前实现是**一维沿基线插值**，不是完整 3 站/4 站二维网络模型
- 当前使用的是**伪距残差代理模型**，不是严格双频电离层反演

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
| `--iono`                          | 启用双基站电离层/网络残差插值 |
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

启动时打印配置摘要，例如：

```text
rtcm_vbs_hp: A=127.0.0.1:50001 -> out=:50002 | mode=absolute-LLH | trop=1 | iono=1 | MSM=7 | SSR fwd=1
        B=127.0.0.1:50003  epoch window=40 ms
```

每 15 秒打印一次统计。

**单基站**

```text
[stat] 1005/6=11 MSM=104 | obs in=2993 out=2993 | sat corr=2993 no_eph=0 | SSR fwd=37 | clients=1
```

**双基站 + `--iono`**

```text
[stat] 1005/6=11 MSM=104 | obs in=2993 out=2993 | sat corr=2993 no_eph=0 net=846 \
       | epochs merged=98 A-only=2 B-only=1 \
       | CNR both_sat=780 chose_A=512 chose_B=268 \
       | iono epochs=96 pairs=742 interp=690 rej=52 \
       | SSR fwd=37 | clients=1
```

| 字段 | 含义 |
|---|---|
| `1005/6` | 已输出的 1005/1006 帧数（始终是 VBS 坐标） |
| `MSM` | 已输出的 MSM 帧数（多星座拆分后） |
| `obs in / out` | 进入修正的卫星数 / 成功输出的卫星数 |
| `sat corr` | 应用了 VBS 修正的累计卫星次数 |
| `no_eph` | 因缺星历（或对应侧 base 未学习）被丢弃的卫星次数 |
| `net` | 实际叠加了 `Δnet` 的累计卫星次数 |
| `epochs merged` | A 和 B 都齐的历元数 |
| `epochs A-only / B-only` | 单侧输出的历元数 |
| `CNR both_sat` | 同卫星两侧都有、做了 CNR 比较的累计次数 |
| `CNR chose_A / chose_B` | CNR 比较时选了 A 或 B 的次数 |
| `iono epochs` | 实际构建了网络/电离层插值的历元数 |
| `pairs` | A/B 共同卫星对总数 |
| `interp` | 生成了插值修正的卫星次数 |
| `rej` | 因残差异常过大被拒绝的次数 |
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
| 长基线空间残差                | 未处理       | 可由 `--iono` 部分吸收 |
| **综合 Δρ / Δmodel 精度**     | **分米–米级** | **毫米级 + 远距增强** |

---

## 10. 常见问题

**Q：日志里 `no_eph` 持续增长？**  
A：星历不足。检查主 `--src-port` 是否包含 1019/1020/1042/1045/1046，或为星历流单独指定 `--eph-port`。

**Q：双基站下固定率还是不理想？**  
A：先依次对比：`--dual`、`--dual --trop`、`--dual --trop --iono`。如果 20 km 以上仍难 fix，通常说明一维双站插值仍不够，需要升级到 3 站/4 站二维网络模型。

**Q：`--iono` 为什么要求 `--dual`？**  
A：当前实现只基于 A/B 两站的共同卫星残差做插值，没有单站版本也没有多站二维版本。

**Q：`rej` 持续增长说明什么？**  
A：A/B 某些共同卫星的残差异常大，可能由周跳、码噪声、多路径、时间不同步、或过长基线导致。当前实现会按阈值拒绝，以避免把异常值直接加到 VBS 上。

**Q：Windows 下怎么用 MSVC 编译？**  
A：请在已安装 Visual Studio 或 Build Tools 的开发者命令行中运行文档里的 CMake 命令。若普通 PowerShell 下提示找不到 `cl` / `msbuild` / `cmake`，说明开发环境尚未加入 PATH。

---

## 11. 当前限制

1. **天线描述**（1007/1008/1033）暂未同步重写。
2. **LLI / lock 计数**按输入原样透传。
3. **双基站**外送的 1005/1006 始终是 VBS 的单一坐标；A/B 仅在内部保留真实 ECEF。
4. **SSR** 走帧级透传，不对 SSR 修正做 VBS 补偿。
5. **`--trop`** 只做差分 Saastamoinen，对湿延迟的实时空间变化建模能力有限。
6. **`--iono` 当前不是完整网络 RTK/VRS**：
   - 仅支持双基站
   - 仅做沿 A→B 基线的一维插值
   - 当前以首个有效伪距残差做代理，不是严格双频电离层反演
   - 更适合作为 10–30 km 远距增强的第一阶段方案
7. **若要进一步提升远距 fix 率**，下一步应考虑：
   - 3 站/4 站二维平面插值
   - 双频电离层提取
   - 更强的周跳与鲁棒权重控制

---

## 12. 代码结构

```text
rtcm_vbs_hp/
├── CMakeLists.txt        # 最小化拉取 RTKLIB 源 + 本程序模块
├── README.md             # 本文件
├── main.c                # 主事件循环 + CLI + 双站插值调用入口
├── vbs_core.h / .c       # 核心 VBS 物理修正（位置 + 观测量）
├── vbs_merge.h / .c      # 双基站 CNR 合并 + 历元对齐
├── net_interp.h / .c     # 双基站电离层/网络残差一维插值
├── tcp_io.h / .c         # TCP 客户端 + 多客户端广播服务器
└── rtklib_glue.c         # 提供 RTKLIB 未定义的 PPP_Glo 全局变量
```

### 主要 RTKLIB API 使用点

| RTKLIB 函数 | 使用目的 |
|---|---|
| `input_rtcm3()` | 输入字节流的增量解码 |
| `gen_rtcm3()` | 输出 MSM4/7、1005/1006、1019/1020/1042/1045/1046 编码 |
| `satposs()` | 发射时刻卫星位置 / 钟差 / 相对论 / BDS GEO |
| `geodist()` | 含 Sagnac 的几何距离 |
| `satazel()` | 高度角（对流层模型用） |
| `tropmodel()` | Saastamoinen 对流层延迟 |
| `ecef2pos() / pos2ecef()` | 大地坐标 ⇄ ECEF |
| `enu2ecef()` | NED 偏移 → ECEF 平移向量 |
| `sat2freq()` | 每颗卫星-每信号的精确载波频率（含 GLONASS FCN） |

---

*文档最后更新：2026-04-20* 
