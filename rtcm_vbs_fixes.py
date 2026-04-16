"""
rtcm_vbs 修复补丁
=================
将此文件中的三个函数替换到原 rtcm_vbs.py 对应位置。

修复列表：
  Fix-1  patch_msm_observations() — tow_s 取模一周，修正 TOW 计算错误
  Fix-2  correct_parsed_msm()     — 同上
  Fix-3  run()                    — 无星历的观测帧继续保持丢弃策略
  Fix-4  配置建议                 — 偏移量示范值
"""

import math

# ─── 保留原文件中已有的常量 ────────────────────────────────────────
C_PER_MS = 299_792.458          # m/ms
WEEK_S   = 604_800.0            # GPS 周长（秒）


# ═══════════════════════════════════════════════════════════════════
# Fix-1  patch_msm_observations  ← 替换原函数
# ═══════════════════════════════════════════════════════════════════
def patch_msm_observations(raw_frame: bytes,
                            base_xyz: tuple,
                            vbs_xyz:  tuple,
                            ephem) -> "bytes | None":
    """
    与原版相同，但修复了 tow_s 计算（取模一周），
    并保持"无星历即丢弃"策略：
      - 仍要求至少一颗星有星历，否则返回 None
    """
    from rtcm_vbs import (
        _get_bits, _set_bits, _as_signed, _recompute_crc,
        _msm_info, _MSM_SIG_LAYOUT, EphemerisStore,
    )

    payload = bytearray(raw_frame[3:-3])
    bit     = 0

    def read(n: int) -> int:
        nonlocal bit
        v = _get_bits(payload, bit, n)
        bit += n
        return v

    mtype       = read(12)
    constl, sub = _msm_info(mtype)
    if sub not in _MSM_SIG_LAYOUT:
        return raw_frame

    read(12)                           # station ID（跳过）
    epoch_ms = read(30)
    bit += 19

    # ── Fix-1 核心：epoch_ms 取模一周 ──────────────────────────────
    tow_s = (epoch_ms / 1000.0) % WEEK_S
    # ──────────────────────────────────────────────────────────────

    sat_mask = read(64)
    sig_mask = read(32)

    nsat = bin(sat_mask).count('1')
    nsig = bin(sig_mask).count('1')
    ncells_total = nsat * nsig
    cell_mask    = read(ncells_total)
    ncell        = bin(cell_mask).count('1')

    sat_ids = [i + 1 for i in range(64) if (sat_mask >> (63 - i)) & 1]

    cell_sat_idx: list = []
    for s in range(nsat):
        for g in range(nsig):
            ci = s * nsig + g
            if (cell_mask >> (ncells_total - 1 - ci)) & 1:
                cell_sat_idx.append(s)

    rough_int_start = bit
    rough_int = [read(8) for _ in range(nsat)]

    if sub in (5, 7):
        bit += 4 * nsat

    rough_mod_start = bit
    rough_mod = [read(10) for _ in range(nsat)]

    if sub in (5, 7):
        bit += 14 * nsat

    sig_data_start = bit

    dx = vbs_xyz[0] - base_xyz[0]
    dy = vbs_xyz[1] - base_xyz[1]
    dz = vbs_xyz[2] - base_xyz[2]

    corrections: list = []
    for sv in sat_ids:
        pos = ephem.sat_ecef(constl, sv, tow_s)
        if pos is None:
            corrections.append(None)
            continue
        rx = pos[0] - base_xyz[0]
        ry = pos[1] - base_xyz[1]
        rz = pos[2] - base_xyz[2]
        dist = math.sqrt(rx * rx + ry * ry + rz * rz)
        if dist < 1e6:
            corrections.append(None)
            continue
        corrections.append(-(rx / dist * dx + ry / dist * dy + rz / dist * dz))

    # ── Fix-2 无星历时继续丢弃，避免位置/观测几何不一致 ───────────
    if not any(c is not None for c in corrections):
        return None

    # ── 以下逻辑与原版完全相同 ────────────────────────────────────
    rough_int_new = list(rough_int)
    rough_mod_new = list(rough_mod)
    residuals: list = []

    for k in range(nsat):
        corr = corrections[k]
        if corr is None or rough_int[k] == 255:
            residuals.append(0.0)
            continue

        dt_ms  = corr / C_PER_MS
        d_mod  = round(dt_ms * 1024)
        new_mod = rough_mod[k] + d_mod
        new_int = rough_int[k]

        while new_mod >= 1024: new_mod -= 1024; new_int += 1
        while new_mod < 0:     new_mod += 1024; new_int -= 1
        new_int = max(0, min(254, new_int))

        rough_int_new[k] = new_int
        rough_mod_new[k] = new_mod
        residuals.append((dt_ms - d_mod / 1024.0) * C_PER_MS)

    for k in range(nsat):
        _set_bits(payload, rough_int_start + k * 8,  8,  rough_int_new[k])
        _set_bits(payload, rough_mod_start + k * 10, 10, rough_mod_new[k])

    lay = _MSM_SIG_LAYOUT[sub]
    fine_pr_start = sig_data_start
    for c in range(ncell):
        res = residuals[cell_sat_idx[c]]
        if res == 0.0:
            continue
        bpos  = fine_pr_start + c * lay['fine_pr']
        old   = _as_signed(_get_bits(payload, bpos, lay['fine_pr']), lay['fine_pr'])
        delta = round(res / (C_PER_MS * lay['pr_scale']))
        _set_bits(payload, bpos, lay['fine_pr'], old + delta)

    fine_ph_start = fine_pr_start + ncell * lay['fine_pr']
    for c in range(ncell):
        res = residuals[cell_sat_idx[c]]
        if res == 0.0:
            continue
        bpos  = fine_ph_start + c * lay['fine_ph']
        old   = _as_signed(_get_bits(payload, bpos, lay['fine_ph']), lay['fine_ph'])
        delta = round(res / (C_PER_MS * lay['ph_scale']))
        _set_bits(payload, bpos, lay['fine_ph'], old + delta)

    new_frame = raw_frame[:3] + bytes(payload) + raw_frame[-3:]
    return _recompute_crc(new_frame)


# ═══════════════════════════════════════════════════════════════════
# Fix-3  correct_parsed_msm  ← 替换原函数（dual 模式也用到）
# ═══════════════════════════════════════════════════════════════════
def correct_parsed_msm(p: dict,
                        base_xyz: tuple, vbs_xyz: tuple,
                        ephem) -> dict:
    """
    与原版相同，但修复了 tow_s 计算（取模一周）。
    """
    import copy
    from rtcm_vbs import _MSM_SIG_LAYOUT

    p   = copy.deepcopy(p)
    lay = _MSM_SIG_LAYOUT[p['sub']]
    dx  = vbs_xyz[0] - base_xyz[0]
    dy  = vbs_xyz[1] - base_xyz[1]
    dz  = vbs_xyz[2] - base_xyz[2]

    # ── Fix-3 核心：取模一周 ──────────────────────────────────────
    tow_s = (p['epoch_ms'] / 1000.0) % WEEK_S
    # ──────────────────────────────────────────────────────────────

    corrections = []
    for sv in p['sat_ids']:
        pos = ephem.sat_ecef(p['constl'], sv, tow_s)
        if pos is None:
            corrections.append(None)
            continue
        rx, ry, rz = pos[0]-base_xyz[0], pos[1]-base_xyz[1], pos[2]-base_xyz[2]
        dist = math.sqrt(rx*rx + ry*ry + rz*rz)
        if dist < 1e6:
            corrections.append(None)
            continue
        corrections.append(-(rx/dist*dx + ry/dist*dy + rz/dist*dz))

    residuals = []
    for k, corr in enumerate(corrections):
        if corr is None or p['rough_int'][k] == 255:
            residuals.append(0.0)
            continue
        dt_ms = corr / C_PER_MS
        d_mod = round(dt_ms * 1024)
        nm    = p['rough_mod'][k] + d_mod
        ni    = p['rough_int'][k]
        while nm >= 1024: nm -= 1024; ni += 1
        while nm <     0: nm += 1024; ni -= 1
        p['rough_int'][k] = max(0, min(254, ni))
        p['rough_mod'][k] = nm
        residuals.append((dt_ms - d_mod / 1024.0) * C_PER_MS)

    for c, (si, _) in enumerate(p['cell_map']):
        res = residuals[si]
        if res == 0.0:
            continue
        p['fine_pr'][c] += round(res / (C_PER_MS * lay['pr_scale']))
        p['fine_ph'][c] += round(res / (C_PER_MS * lay['ph_scale']))

    return p


# ═══════════════════════════════════════════════════════════════════
# Fix-4  配置建议（修改原文件底部 __main__ 块）
# ═══════════════════════════════════════════════════════════════════
CONFIG_ADVICE = """
# ── 建议配置（替换原 __main__ 块中的值）──
MODE = "single"

# 先用小偏移量验证流程，确认 Fix 后再逐步增大
OFFSET_NORTH = 100     # m（原来 1000，过大导致 fine_pr 溢出）
OFFSET_EAST  = 50.0    # m（原来 500）
OFFSET_UP    = 0       # m（高程偏移通常不需要，先归零）

# EPHEM_PORT 必须指向能提供 1019/1042/1045/1046 的星历流
# 如果主 RTCM 流（50001）本身包含星历，可令：
#   EPHEM_PORT = None
# 并在 run() 中把下面这行改为直接从主流提取星历（原代码已支持）
EPHEM_PORT = None   # 或 50010（若有单独星历流）
"""

print(CONFIG_ADVICE)
