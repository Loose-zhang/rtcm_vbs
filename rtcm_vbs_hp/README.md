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
   - [5.4 双基站 A 优先合并](#54-双基站-a-优先合并模式)
   - [5.5 双基站远距增强（电离层/网络插值）](#55-双基站远距增强电离层网络插值)
6. [命令行参数](#6-命令行参数)
7. [消息处理规则](#7-消息处理规则)
8. [日志与统计](#8-日志与统计)
9. [与 Python 版精度对比](#9-与-python-版精度对比)
10. [常见问题](#10-常见问题)
11. [当前限制](#11-当前限制)
12. [代码结构](#12-代码结构)
13. [设计与修复记录](#13-设计与修复记录)

---

## 1. 特性概览

| 能力 | 说明 |
|---|---|
| 全精度 VBS 修正 | `satposs()`（τ 迭代 + 相对论 + 钟差 + BDS GEO）+ `geodist()`（Sagnac）|
| 全星座 | GPS / GLONASS / Galileo / BeiDou / QZSS / IRNSS |
| 差分对流层 | Saastamoinen 差分（可选 `--trop`），天向偏移大时尤为重要 |
| 载波相位一致 | `L += Δρ/λ`，λ 由 `sat2freq()` 按信号逐一精确计算，保持整周模糊连续 |
| 相位重编码 | `gen_rtcm3()`，输出真 RTCM3 MSM4 / MSM7，无 fine_pr / fine_ph 量化溢出 |
| 双基站 A 主站辅助改正 | `--dual`，A/B 历元对齐；外送 MSM 始终使用 A 观测，B 仅作内部辅助改正源 |
| 远距电离层/网络插值 | `--iono`，基于 A/B 同卫星伪距残差减几何差的结果，按 VBS 到 A/B 的距离权重插值 |
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
| 双基站合并 | 有（`dual` 模式） | 有（A 主站外送 + B 内部辅助改正 + 历元对齐） |
| 长基线空间残差 | 无显式处理 | 新增 `--iono` 双站距离加权网络残差插值 |

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

4. 计算 VBS 到 A/B 两站的 3D 距离反比权重 `wa/wb`
5. 可用时把残差按电离层映射函数归一到近似天顶域，再映射回 VBS 方向
6. 得到该卫星的附加修正：

```text
Δnet = wb * dres                  # 极端几何 fallback
Δnet = wb * (dres / mf_mean) * mf_v # 常规电离层映射修正路径
```

这部分主要吸收长基线下的空间相关残差，尤其是电离层导致的差异。当前版本是**双基站、距离加权残差插值**，不是完整网络 RTK/VRS 模型，适合作为远距增强的第一阶段方案。

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

### 5.4 双基站 A 优先合并模式

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
| 同一历元仅 B 有某卫星 | 不外送该观测，B 仅作为内部改正源 |
| A、B 都有该卫星 | 使用 A 的观测，保持 A 为共享卫星的相位参考 |
| A、B 时间差 <= `--window MS` | 作为同一合成历元输出，并把两侧观测时间统一到 A 历元 |
| A、B 时间差 > `--window MS` | 先输出较旧的那一帧，不丢历元 |
| 仅一侧到达、等满 `--window MS` | 单侧输出 |

B 站不再作为外送 MSM 的观测来源。双站模式的输出观测始终以 A 为主：A 有的卫星才进入输出，B 只用于共同卫星状态统计和可选 `--iono` 网络残差改正，避免 RTKNAVI 看到来自两个物理接收机的 MSM 观测集跳变。
双站配对不再使用 10 ms 固定阈值，而是统一使用 `--window MS`。只要 A/B 均已到达且时间差在窗口内，无论是否有共同卫星，都会输出为同一个合成历元，以保证输出 MSM 的历元时间一致。

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
- 当前实现是**双站距离加权残差插值**，不是完整 3 站/4 站二维网络模型
- 当前使用的是**伪距残差代理模型**，不是严格双频电离层反演
- 当前 A-master 策略下，B 不直接外送观测，只通过 `Δnet` 影响 A 的输出观测
- B 短时休眠或中断时，最后一次有效 `Δnet` 会保持 5 s，然后在 30 s 内线性衰减到 0，避免改正量突变导致 rover 从固定解掉到浮点解

---

## 6. 命令行参数

| 参数 | 说明 |
|---|---|
| `--src-host HOST`                 | 主源 RTCM TCP 主机（默认 `127.0.0.1`）|
| `--src-port PORT`                 | 主源端口（默认 `50001`）|
| `--dual`                          | 启用双基站 A 优先合并模式 |
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
| **1005 / 1006** | 学习真实 ARP → 计算 VBS ARP → 以 **1006** 重新编码后发送；双基站模式只外送 A 站号下的单一 VBS 1006，B 站 1005/1006 仅内部使用 |
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
| `epochs A-only / B-only` | 单侧到达或超时的历元数；当前 A-master 输出策略下 B-only 不外送 MSM |
| `CNR both_sat` | 日志字段名沿用旧称；含义是 A/B 同历元共同卫星累计次数 |
| `CNR chose_A / chose_B` | 日志字段名沿用旧称；当前实现中外送观测只计入 `chose_A`，`chose_B` 保留为旧策略字段 |
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
A：先依次对比：`--dual`、`--dual --trop`、`--dual --trop --iono`。如果 20 km 以上仍难 fix，通常说明双站距离加权残差插值仍不够，需要升级到 3 站/4 站二维网络模型。

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
3. **双基站**外送的 1005/1006 和 MSM 始终使用 A 的站号和 VBS 单一坐标；B 的站号、真实 ECEF 仅在内部保留，避免 rover/RTKNAVI 识别成两个参考站。
4. **SSR** 走帧级透传，不对 SSR 修正做 VBS 补偿。
5. **`--trop`** 只做差分 Saastamoinen，对湿延迟的实时空间变化建模能力有限。
6. **`--iono` 当前不是完整网络 RTK/VRS**：
   - 仅支持双基站
   - 仅做双站距离加权残差插值
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
├── vbs_merge.h / .c      # 双基站 A 优先合并 + B 相位归一 + 历元对齐
├── net_interp.h / .c     # 双基站电离层/网络残差距离加权插值
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

## 13. 设计与修复记录

本节合并原独立修复文档中的长期有用信息，便于后续只维护 `README.md`。

### 13.1 网络/电离层插值修复

`net_interp_build_dual()` 的输入是 A/B 同历元同卫星的首个有效伪距，基础残差为：

```text
dres = (P_B - P_A) - (rho_B - rho_A)
```

该残差用于估计 A→B 方向上的空间相关误差，主要面向长基线下的电离层差异。当前实现包含三个关键修正。

**对流层重复改正避免**

当 `--iono` 与 `--trop` 同时开启时，`dres` 中天然包含 A/B 站间对流层差 `trop(B) - trop(A)`。而 `vbs_core.c` 后续还会分别对 A 或 B 观测施加 `trop(VBS) - trop(base)`。为避免对流层项重复进入 `dnet`，当前实现会在插值前从 `dres` 中减去站间对流层差：

```c
if (vbs->apply_trop) {
    double ta = tropmodel(oa->time, pos_a, azel_a, vbs->humi);
    double tb = tropmodel(ob->time, pos_b, azel_b, vbs->humi);
    dres -= (tb - ta);
}
```

因此 `--dual --iono --trop` 可以一起使用：对流层由 `vbs_core.c` 负责，`net_interp.c` 输出更接近电离层/空间残差的补偿。

**距离反比权重替代一维基线投影**

旧思路按 VBS 在 A→B 连线上的投影比例插值。当前实现改为使用 VBS 到 A/B 两站的 3D 距离反比权重：

```text
wa = (1 / dist(VBS,A)) / ((1 / dist(VBS,A)) + (1 / dist(VBS,B)))
wb = (1 / dist(VBS,B)) / ((1 / dist(VBS,A)) + (1 / dist(VBS,B)))
```

`wa + wb = 1`，自然限定在 `[0, 1]`。以 A 站为参考零点时，VBS 侧附加修正的基础形式为 `corr = wb * dres`。

**电离层映射函数归一化**

电离层延迟是斜路径量，受卫星仰角影响。当前实现用 RTKLIB 的 `ionmapf()` 将 A/B 残差近似归一到天顶域，再映射回 VBS 方向：

```c
double mf_a = ionmapf(pos_a, azel_a);
double mf_b = ionmapf(pos_b, azel_b);
double mf_mean = 0.5 * (mf_a + mf_b);
double dres_vtec = dres / mf_mean;
double mf_v = wa * mf_a + wb * mf_b;
double corr = wb * dres_vtec * mf_v;
```

若映射因子异常或低仰角风险过高，代码退化为 `corr = wb * dres`，避免把低仰角噪声放大。

使用建议：

| 启动参数组合 | 说明 |
|---|---|
| `--dual --iono` | 使用双站残差插值，不额外做显式对流层差分 |
| `--dual --iono --trop` | 推荐用于较长基线或高差明显场景；当前实现已避免对流层重复改正 |
| 基线 < 50 km | IDW + 映射函数修正通常更稳 |
| 基线 > 50 km | 双站模型可能不够，应考虑 3/4 站二维网络模型 |

### 13.2 双基站载波相位对齐

双站 VBS 的主要风险不是几何改正，而是 A/B 两台真实接收机的载波相位各自包含不同接收机偏差和模糊度。若逐卫星在 A/B 之间直接切换，输出 `L` 可能出现几百周量级台阶，RTK 滤波会长期无法固定。

当前双站数据流为：

```text
A/B 原始 obs
   -> merger_poll_pair() 返回两侧原始历元
   -> 分别 vbs_correct_obs_ex(side=0/1) 修正到 VBS 几何位置
   -> merger_align_and_merge() 更新 B->A 载波偏置、归一化 B、合并输出
   -> emit_msm()
```

`merger_t` 内维护：

```c
align_state_t align[MAXSAT][NFREQ + NEXOBS];
```

状态按 `(sat-1, signal code)` 约束，避免 A/B 信号槽顺序不一致时错配。共同卫星、共同信号可见时计算：

```text
diff_cyc = L_A_vbs - L_B_vbs
```

偏置状态规则：

| 阶段 | 条件 | 动作 |
|---|---|---|
| 初始化 | 连续 3 次 `diff` 相对均值稳定在 0.15 cycle 内 | 锁定 `bias_cyc`，置 `valid=1` |
| 跟踪 | `diff - bias` 在 0.15 cycle 内 | 用 0.05 IIR 系数慢速更新 bias |
| reset | A/B 任一侧 LLI 标记周跳 | 清空状态 |
| reset | 信号码变化 | 清空状态 |
| reset | 配对中断超过 5 s | 清空状态 |
| reset | 残差超过 0.75 cycle | 清空状态 |

B 侧观测输出策略：

| 情况 | 处理 |
|---|---|
| A/B 共同卫星 | 输出 A，B 只参与内部状态和网络残差 |
| B-only 卫星 | 不外送 |
| B-only 历元 | 不外送 |

B 的伪距和载波都不直接进入外送 MSM。这样牺牲 B 补星能力，换取 RTKNAVI 端稳定的单接收机载波观测集。

### 13.3 A 主站 + B 辅助改正策略

当前双站合并的目标是：A 持续播发时保持 A 的相位连续性；B 休眠、断开或仅提供部分星座时，不阻塞 A；B 仅作为辅助改正源，不改变外送 MSM 的观测来源集合。

最终策略：

| 场景 | 输出行为 |
|---|---|
| A/B 共同可见同一卫星 | 固定使用 A，B 只用于内部统计和 `--iono` 残差 |
| A 有、B 无 | 输出 A |
| A 无、B 有 | 不输出该卫星 |
| 仅 B 历元到达 | 不输出 MSM |
| A/B 时间差 <= `--window` | 合成为同一历元，输出时间统一为 A |
| A/B 时间差 > `--window` | 先输出较旧的一侧，不强行等待 |
| B 休眠超过 `--window` | A-only 持续输出 |
| A 休眠超过 `--window` | B-only 不外送 |

这种策略刻意取消了共享卫星按 SNR 在 A/B 间切换，也取消了 B 补星外送，优先保护载波相位连续性和 RTKNAVI 固定解稳定性。

### 13.4 多输入流容错

上游 TCP client 当前为非阻塞读取：

| 情况 | 处理 |
|---|---|
| 某通道暂无数据 | `tcpc_read()` 返回 `TCPC_AGAIN`，主循环继续处理其他通道 |
| eph 中断 | 使用已经缓存到 `rout.nav` 的星历继续处理观测；新星历恢复后继续镜像 |
| ssr 中断 | 停止 SSR 透传，但观测流继续 |
| B 休眠或断开 | A-only 按窗口超时持续输出 |
| A 断开、B 正常 | B 可继续参与内部缓存，但双基站输出不会切换到 B 站号；若尚未学习 A 站号则暂不外送 MSM |
| 真实断线 | 关闭 socket，并按 1 秒退避重连 |

连接尝试使用短超时，避免 eph/ssr/B 长期不可达时拖慢主观测流。

#### 多通道读空与接收缓冲（实现说明）

当同时开启 A / B / 独立星历 / 独立 SSR 等多路 TCP 输入时，若主循环对**每个通道每轮只 `recv` 一次**（例如 4 KB），码率较高的一路可能在**轮询其他通道期间**持续堆积本机内核接收缓冲；对端发送缓冲被反压后，部分转发或中间件会**主动断开连接**， stderr 上表现为频繁的 `connected` / `disconnected`，并影响星历与观测的及时处理。

当前实现做了两点缓解（`main.c` / `tcp_io.c`）：

| 措施 | 说明 |
|---|---|
| **读至 `EAGAIN`** | 对每个已连接通道，在内层循环中反复 `tcpc_read()`，直到返回 `TCPC_AGAIN`（当前无更多可读数据），再处理下一通道，避免“每路每轮只啃一口”造成的慢消费。 |
| **套接字调优** | 在 `tcpc_connect` 成功后对客户端套接字设置较大的 `SO_RCVBUF`（约 4 MB）与 `SO_KEEPALIVE`，提高突发余量并便于发现半开连接。 |

若对端是**单客户端、新连接会踢掉旧连接**等策略，仍可能出现断线，需在**上游转发/服务**侧调整，而非仅靠本程序。

### 13.5 验证建议

基础验证：

```bash
cmake --build build -j
```

运行验证建议：

1. 单站回归：`--offset 3000 3000 30 --trop`，确认仍可稳定输出和固定。
2. 双站基础：`--dual --offset 3000 3000 30 --trop`，观察 `epochs merged / A-only / B-only`。
3. 双站 + 插值：`--dual --trop --iono`，观察 `iono epochs / pairs / interp / rej`，`rej` 不应持续飙升。
4. 断流测试：分别断开 eph、ssr、A、B，确认仍在线观测流持续推送，并在恢复后自动重连。
5. B 休眠测试：A 持续播发、B 间歇播发时，输出应在 merged / A-only 间自然切换。

---

*文档最后更新：2026-04-27*
