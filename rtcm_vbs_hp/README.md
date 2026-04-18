# rtcm_vbs_hp — 高精度 RTCM3 虚拟基站

基于 **RTKLIB** (`../src` + `../include`) 重构的 C 实现，替代原 Python 版本。

## 为什么精度更高

原 Python 版通过**字节原位修改** MSM 帧中的 `rough_int / rough_mod / fine_pr / fine_ph` 字段实现 VBS，存在以下天花板：

| 问题                                      | 原 Python 版                    | 本 C 版 (RTKLIB)                                    |
| ----------------------------------------- | ------------------------------- | --------------------------------------------------- |
| 卫星位置计算                              | 手写简化 Kepler，无 τ 迭代      | `satposs()` 含 `τ = ρ/c` 迭代 + 相对论 + 钟差 + BDS GEO |
| 几何距离                                  | 欧氏距离，无 Sagnac             | `geodist()` 内含 Sagnac（~30 m 系统偏差消除）       |
| TOW 计算                                  | 早期版本未取模一周（已用 fix 补）| 直接使用 RTKLIB `gtime_t`，无任何时间跳变           |
| MSM 字段量化                              | `rough_mod` 步长 293 m，大偏移溢出 | 完整解码→物理量修正→`gen_rtcm3()` 重新编码，无溢出 |
| 对流层差分 (高度差异影响)                 | 未处理                          | 可选 Saastamoinen 差分 (`--trop`)                   |
| 支持星座                                  | 仅 GPS / GAL / BDS              | GPS / GLO / GAL / BDS / QZS / IRN (RTKLIB 已启用)   |
| 相位连续性                                | fine_ph 直接偏移，易破坏整周模糊 | `L_new = L_old + Δρ/λ`，RTKLIB 全精度浮点累计        |

数学上，对每颗卫星：

```
Δρ_i  =  geodist(r_sat^i, r_VBS)  -  geodist(r_sat^i, r_base)
P_new = P_old + Δρ
L_new = L_old + Δρ / λ            (+ 可选 Δtrop 项)
```

所有卫星位置 `r_sat^i` 由 `satposs()` 统一按信号发射时刻算出。

## 编译

```bash
cd rtcm_vbs_hp
cmake -B build
cmake --build build -j
```

产物：`build/rtcm_vbs_hp`

依赖：POSIX (`pthread`)，macOS/Linux 原生编译即可；Windows 下 MinGW + winsock2。

## 使用

### 最小示例（NED 偏移模式）

```bash
# 真实基站 RTCM 流在 :50001
str2str -in ntrip://user:pass@host:2101/MOUNT -out tcpsvr://:50001 &

# 生成一个在真实基站北 1000 m、东 500 m、上 50 m 的虚拟基站，输出到 :50002
./build/rtcm_vbs_hp \
    --src-port 50001 \
    --out-port 50002 \
    --offset   1000.0  500.0  50.0 \
    --trop

# 接收端（如 rtkrcv）指向 tcp://127.0.0.1:50002
```

### 绝对经纬度模式

```bash
./build/rtcm_vbs_hp \
    --src-port 50001 \
    --out-port 50002 \
    --vbs  22.937764  113.206806  120.253
```

### 独立星历流

如果主挂载点不含 1019/1020/1042/1045/1046：

```bash
str2str -in ntrip://.../EPH_MOUNT -out tcpsvr://:50010 &

./build/rtcm_vbs_hp --src-port 50001 --eph-port 50010 --out-port 50002 \
    --offset 1000 500 50 --trop
```

## 命令行参数

| 参数                         | 说明                                           |
| ---------------------------- | ---------------------------------------------- |
| `--src-host HOST`            | 源 RTCM TCP 主机（默认 `127.0.0.1`）           |
| `--src-port PORT`            | 源 RTCM TCP 端口（默认 `50001`）               |
| `--eph-host/--eph-port`      | 可选的独立星历流                               |
| `--out-port PORT`            | VBS 输出端口（默认 `50002`）                   |
| `--offset N E U`             | NED 偏移模式（米）                             |
| `--vbs LAT LON ALT`          | 绝对坐标模式（度, 度, 米）                     |
| `--trop`                     | 启用差分 Saastamoinen 对流层修正               |
| `--msm N`                    | 输出 MSM 子类型：`4` 或 `7`（默认 `7`）        |
| `--no-eph`                   | 不转发星历消息                                 |
| `--trace FILE` / `--level N` | 打开 RTKLIB trace 日志                         |

## 与 Python 版本精度对比（理论）

在 1 km NED 偏移、静态基准条件下：

| 误差源                         | Python 版   | rtcm_vbs_hp |
| ------------------------------ | ----------- | ------------ |
| 卫星位置（无 τ 迭代）          | 20–60 m     | < 1 cm       |
| Sagnac 未修正                  | ~30 m       | 0            |
| rough_mod 量化（293 m 步长）   | 亚米~米级 | 无（重编码）|
| 对流层 (VBS 高差 50 m)         | 未修正 ~2 dm| 可修正到 mm |
| 综合 Δρ 精度                   | 分米~米级  | **毫米级**  |

下游 RTK 求解精度自然随之提升；载波相位连续性亦得以保持，有利于整周模糊度 fix。

## 当前限制

1. **SSR 消息**（1057–1068 等）暂未透传，仅转发 MSM + 1005/6 + 广播星历。
2. **天线描述**（1007/1008/1033）暂未同步重写；下游通常无影响。
3. 时间对齐假设上游观测帧之间时序正常；极端丢包下自动重连。
4. 整周 LLI / 锁定计数按输入原样透传，未做 VBS 额外修正（理论上 VBS 不改变失锁事件，这样处理正确）。
5. 双基站 CNR 合并模式（原 Python dual 模式）尚未移植如需继续扩展，下一步推荐先做 SSR 透传（只需多一个 gen_rtcm3 的 case），然后再做双基站合并（可复vbs_correct_obs 的输出，在 emit_msm 之前按星座合并 A/B 的 obsd_t 数组）。
