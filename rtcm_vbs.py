#!/usr/bin/env python3
"""
RTCM3 Virtual Base Station (VBS) Interpolator
==============================================

Reads an RTCM3 stream from a TCP port (e.g. str2str output on 50001),
applies a configurable North/East/Up offset to produce a geometrically
consistent Virtual Base Station stream, then re-serves it on another port.

Two layers of correction are applied:

1. **Position (MSG 1005/1006)**
   The ARP ECEF coordinates are shifted to the virtual location.

2. **Observations (MSM4/5/6/7)**
   Pseudoranges and carrier phases are corrected by the geometric range
   change  Δρ = −û·Δr  for each satellite, where û is the unit vector
   from the real base to the satellite (computed from broadcast ephemeris
   in MSG 1019 / 1042 / 1045 / 1046) and Δr is the ECEF displacement.
   The correction is split between the rough-range field (satellite-level)
   and the fine-PR / fine-phase fields (signal-level).

Only GPS / Galileo / BeiDou are forwarded. Unsupported constellations are
filtered out so the VBS geometry stays consistent.

Usage
-----
1. Start str2str to push the base-station RTCM to port 50001:
       str2str -in ntrip://... -out tcpsvr://:50001

2. Run this script:
       python rtcm_vbs.py

3. Point your rover / RTKLIB to   tcp://127.0.0.1:50002

Configuration
-------------
Edit the  ── Configuration ──  block at the bottom of this file.
"""

import copy
import math
import socket
import struct
import threading
import time
import logging

from pyrtcm import RTCMReader
from pyrtcm.exceptions import RTCMParseError

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("rtcm_vbs")

# ─────────────────────────────────────────────────────────────────────────────
# WGS-84
# ─────────────────────────────────────────────────────────────────────────────
WGS84_A  = 6_378_137.0              # semi-major axis (m)
WGS84_F  = 1.0 / 298.257_223_563    # flattening
WGS84_B  = WGS84_A * (1.0 - WGS84_F)
WGS84_E2 = 2.0 * WGS84_F - WGS84_F ** 2   # first eccentricity²

# ─────────────────────────────────────────────────────────────────────────────
# RTCM 1005 / 1006  bit layout (offsets inside the *payload*, MSB-first)
#
#  DF002  12 bits   [  0.. 11]  Message Number
#  DF003  12 bits   [ 12.. 23]  Reference Station ID
#  DF021   6 bits   [ 24.. 29]  ITRF Realisation Year
#  DF022   1 bit    [ 30]       GPS Indicator
#  DF023   1 bit    [ 31]       Secondary-system Indicator
#  DF024   1 bit    [ 32]       Galileo Indicator (reserved)
#  DF141   1 bit    [ 33]       Reference-Station Indicator
#  DF025  38 bits   [ 34.. 71]  ARP ECEF-X  (signed, ×0.0001 m)
#  DF142   1 bit    [ 72]       Single Receiver Oscillator Indicator
#  DF001   1 bit    [ 73]       Reserved
#  DF026  38 bits   [ 74..111]  ARP ECEF-Y  (signed, ×0.0001 m)
#  DF364   2 bits   [112..113]  Quarter Cycle Indicator
#  DF027  38 bits   [114..151]  ARP ECEF-Z  (signed, ×0.0001 m)
#
#  1006 appends:
#  DF028  16 bits   [152..167]  Antenna Height (unsigned, ×0.0001 m)
# ─────────────────────────────────────────────────────────────────────────────
_BIT_X     = 34
_BIT_Y     = 74
_BIT_Z     = 114
_BIT_GPS   = 30
_BIT_SYS   = 31
_BIT_GAL   = 32
_XYZ_BITS  = 38
_XYZ_SCALE = 1e-4        # raw unit → metres
_PREAMBLE  = 0xD3
ALLOWED_CONSTELLATIONS = frozenset(("GPS", "GAL", "BDS"))


# ─────────────────────────────────────────────────────────────────────────────
# Coordinate helpers
# ─────────────────────────────────────────────────────────────────────────────

def ecef_to_geodetic(x: float, y: float, z: float):
    """ECEF (m) → (lat °, lon °, alt m)  — Bowring iterative method."""
    lon = math.atan2(y, x)
    p   = math.hypot(x, y)
    # initial guess
    theta = math.atan2(z * WGS84_A, p * WGS84_B)
    lat = math.atan2(
        z + (WGS84_A ** 2 - WGS84_B ** 2) / WGS84_B * math.sin(theta) ** 3,
        p - WGS84_E2 * WGS84_A         * math.cos(theta) ** 3,
    )
    N   = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    alt = p / math.cos(lat) - N
    return math.degrees(lat), math.degrees(lon), alt


def geodetic_to_ecef(lat_d: float, lon_d: float, alt: float):
    """(lat °, lon °, alt m) → ECEF (m)."""
    lat = math.radians(lat_d)
    lon = math.radians(lon_d)
    N   = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    x   = (N + alt) * math.cos(lat) * math.cos(lon)
    y   = (N + alt) * math.cos(lat) * math.sin(lon)
    z   = (N * (1.0 - WGS84_E2) + alt) * math.sin(lat)
    return x, y, z


def apply_ned_offset(lat: float, lon: float, alt: float,
                     dn: float, de: float, du: float = 0.0):
    """
    Shift a geodetic position by (dn, de, du) metres in North / East / Up.
    Returns new (lat °, lon °, alt m).
    """
    R     = WGS84_A                     # local Earth radius (approx.)
    dlat  = math.degrees(dn / R)
    dlon  = math.degrees(de / (R * math.cos(math.radians(lat))))
    return lat + dlat, lon + dlon, alt + du


# ─────────────────────────────────────────────────────────────────────────────
# Bit-level helpers
# ─────────────────────────────────────────────────────────────────────────────

def _get_bits(data: bytes, start: int, n: int) -> int:
    """Extract n unsigned bits from data starting at bit offset start (MSB first)."""
    val = 0
    for i in range(n):
        byte_i = (start + i) >> 3
        bit_i  = 7 - ((start + i) & 7)
        val    = (val << 1) | ((data[byte_i] >> bit_i) & 1)
    return val


def _set_bits(ba: bytearray, start: int, n: int, value: int) -> None:
    """Write n bits of value into bytearray ba at bit offset start (MSB first, two's complement)."""
    if value < 0:
        value += (1 << n)           # convert signed → two's complement
    value &= (1 << n) - 1          # mask to n bits
    for i in range(n):
        byte_i = (start + i) >> 3
        bit_i  = 7 - ((start + i) & 7)
        if (value >> (n - 1 - i)) & 1:
            ba[byte_i] |=  (1 << bit_i)
        else:
            ba[byte_i] &= ~(1 << bit_i)


def _as_signed(val: int, bits: int) -> int:
    """Reinterpret unsigned val as a signed integer of given bit width."""
    if val >= (1 << (bits - 1)):
        val -= (1 << bits)
    return val


def _constellation_allowed(constl: str | None) -> bool:
    return constl in ALLOWED_CONSTELLATIONS


def _set_supported_system_flags(payload: bytearray) -> None:
    # RTCM 1005/1006 has flags for GPS / secondary-system / Galileo only.
    _set_bits(payload, _BIT_GPS, 1, 1)
    _set_bits(payload, _BIT_SYS, 1, 0)
    _set_bits(payload, _BIT_GAL, 1, 1)


# ─────────────────────────────────────────────────────────────────────────────
# CRC-24Q
# ─────────────────────────────────────────────────────────────────────────────

def crc24q(data: bytes) -> int:
    crc = 0
    for b in data:
        crc ^= b << 16
        for _ in range(8):
            crc <<= 1
            if crc & 0x1000000:
                crc ^= 0x1864CFB
    return crc & 0xFFFFFF


def _recompute_crc(frame: bytes) -> bytes:
    """Return a copy of an RTCM3 frame with its CRC-24Q recalculated."""
    crc  = crc24q(frame[:-3])
    tail = struct.pack(">I", crc)[1:]   # 3 bytes big-endian
    return frame[:-3] + tail


def _frame_crc(frame: bytes) -> int:
    return struct.unpack(">I", b"\x00" + frame[-3:])[0]


# ─────────────────────────────────────────────────────────────────────────────
# 1005 / 1006 position patch
# ─────────────────────────────────────────────────────────────────────────────

def patch_arp_position(raw_frame: bytes,
                       dn: float, de: float, du: float = 0.0,
                       verbose: bool = True) -> bytes:
    """
    Given a raw RTCM3 1005 or 1006 frame, return a new frame whose
    ARP ECEF position has been shifted by (dn, de, du) metres (N/E/U).
    The CRC is recomputed.
    """
    # Payload sits between the 3-byte header and the 3-byte CRC
    payload = bytearray(raw_frame[3:-3])

    # ── Read raw integer coords ──
    x_raw = _as_signed(_get_bits(payload, _BIT_X, _XYZ_BITS), _XYZ_BITS)
    y_raw = _as_signed(_get_bits(payload, _BIT_Y, _XYZ_BITS), _XYZ_BITS)
    z_raw = _as_signed(_get_bits(payload, _BIT_Z, _XYZ_BITS), _XYZ_BITS)

    x = x_raw * _XYZ_SCALE
    y = y_raw * _XYZ_SCALE
    z = z_raw * _XYZ_SCALE

    # ── ECEF → geodetic → offset → ECEF ──
    lat, lon, alt = ecef_to_geodetic(x, y, z)
    lat_v, lon_v, alt_v = apply_ned_offset(lat, lon, alt, dn, de, du)
    xv, yv, zv = geodetic_to_ecef(lat_v, lon_v, alt_v)

    if verbose:
        log.info(
            "Base ARP  : lat=%+.8f°  lon=%+.8f°  alt=%.3f m",
            lat, lon, alt,
        )
        log.info(
            "Base ECEF : X=%.4f  Y=%.4f  Z=%.4f m", x, y, z,
        )
        log.info(
            "Offset applied: N=%.1f m  E=%.1f m  U=%.1f m", dn, de, du,
        )
        log.info(
            "VBS  ARP  : lat=%+.8f°  lon=%+.8f°  alt=%.3f m",
            lat_v, lon_v, alt_v,
        )
        log.info(
            "VBS  ECEF : X=%.4f  Y=%.4f  Z=%.4f m", xv, yv, zv,
        )

    # ── Write new coords into payload ──
    _set_bits(payload, _BIT_X, _XYZ_BITS, round(xv / _XYZ_SCALE))
    _set_bits(payload, _BIT_Y, _XYZ_BITS, round(yv / _XYZ_SCALE))
    _set_bits(payload, _BIT_Z, _XYZ_BITS, round(zv / _XYZ_SCALE))
    _set_supported_system_flags(payload)

    new_frame = raw_frame[:3] + bytes(payload) + raw_frame[-3:]
    return _recompute_crc(new_frame)


def _msg_type(frame: bytes) -> int:
    """Return the 12-bit message type from a raw RTCM3 frame."""
    return _get_bits(frame, 24, 12)   # payload starts at byte 3 = frame bit 24


# ─────────────────────────────────────────────────────────────────────────────
# GNSS constants  (for satellite position computation)
# ─────────────────────────────────────────────────────────────────────────────
C_LIGHT  = 299_792_458.0          # m/s
C_PER_MS = C_LIGHT * 1e-3        # metres per millisecond  (299 792.458 m/ms)
GM_GPS   = 3.986_005e14           # WGS-84 gravitational parameter  (m³/s²)
OMEGA_E  = 7.292_115_1467e-5     # WGS-84 Earth rotation rate  (rad/s)

# MSM constellation bases — message number = base + subtype (subtype 1-7)
_MSM_CONSTL = (
    ('GPS', 1070), ('GAL', 1090),
    ('SBS', 1100), ('QZS', 1110), ('BDS', 1120), ('NAV', 1130),
)


def _msm_info(mtype: int):
    """Return (constellation_name, subtype 1-7) or (None, None)."""
    for name, base in _MSM_CONSTL:
        if base < mtype <= base + 7:
            return name, mtype - base
    return None, None


# ─────────────────────────────────────────────────────────────────────────────
# Keplerian satellite position  (GPS / Galileo / BeiDou)
# ─────────────────────────────────────────────────────────────────────────────

def _kepler_ecef(eph: dict, tow_s: float):
    """Return satellite ECEF (m) from broadcast Keplerian elements at TOW (s)."""
    dt = tow_s - eph['t_oe']
    if dt >  302_400: dt -= 604_800
    if dt < -302_400: dt += 604_800

    a = eph['sqrtA'] ** 2
    n = math.sqrt(GM_GPS / a ** 3) + eph['delta_n']
    M = eph['M0'] + n * dt

    e, E = eph['e'], float(M)
    for _ in range(12):                              # Newton–Raphson for E
        E -= (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))

    v   = math.atan2(math.sqrt(1.0 - e * e) * math.sin(E), math.cos(E) - e)
    phi = v + eph['omega']
    s2, c2 = math.sin(2 * phi), math.cos(2 * phi)

    u = phi              + eph['C_us'] * s2 + eph['C_uc'] * c2
    r = a * (1.0 - e * math.cos(E)) + eph['C_rs'] * s2 + eph['C_rc'] * c2
    i = eph['i0'] + eph['IDOT'] * dt  + eph['C_is'] * s2 + eph['C_ic'] * c2

    xp, yp = r * math.cos(u), r * math.sin(u)
    Om  = eph['Omega0'] + (eph['OmegaDot'] - OMEGA_E) * dt - OMEGA_E * eph['t_oe']
    co, so = math.cos(Om), math.sin(Om)
    ci, si = math.cos(i),  math.sin(i)

    return (xp * co - yp * ci * so,
            xp * so + yp * ci * co,
            yp * si)


# ─────────────────────────────────────────────────────────────────────────────
# Ephemeris store  (GPS 1019 / Galileo 1045-1046 / BeiDou 1042)
# ─────────────────────────────────────────────────────────────────────────────

class EphemerisStore:
    """Accumulates Keplerian broadcast ephemeris from RTCM navigation messages."""

    # (constellation, sv_id_field, {orbit_key: DF_attribute})
    _MAPS: dict = {
        1019: ('GPS', 'DF009', dict(
            sqrtA='DF092', e='DF090',    M0='DF088',    delta_n='DF087',
            omega='DF099', Omega0='DF095', OmegaDot='DF100',
            i0='DF097',    IDOT='DF079',
            C_uc='DF089',  C_us='DF091', C_rc='DF098',  C_rs='DF086',
            C_ic='DF094',  C_is='DF096', t_oe='DF093',
        )),
        1045: ('GAL', 'DF252', dict(       # Galileo F/NAV
            sqrtA='DF303', e='DF301',    M0='DF299',    delta_n='DF298',
            omega='DF310', Omega0='DF306', OmegaDot='DF311',
            i0='DF308',    IDOT='DF292',
            C_uc='DF300',  C_us='DF302', C_rc='DF309',  C_rs='DF297',
            C_ic='DF305',  C_is='DF307', t_oe='DF304',
        )),
        1046: ('GAL', 'DF252', dict(       # Galileo I/NAV  (same DF numbering)
            sqrtA='DF303', e='DF301',    M0='DF299',    delta_n='DF298',
            omega='DF310', Omega0='DF306', OmegaDot='DF311',
            i0='DF308',    IDOT='DF292',
            C_uc='DF300',  C_us='DF302', C_rc='DF309',  C_rs='DF297',
            C_ic='DF305',  C_is='DF307', t_oe='DF304',
        )),
        1042: ('BDS', 'DF488', dict(       # BeiDou
            sqrtA='DF504', e='DF502',    M0='DF500',    delta_n='DF499',
            omega='DF511', Omega0='DF507', OmegaDot='DF512',
            i0='DF509',    IDOT='DF491',
            C_uc='DF501',  C_us='DF503', C_rc='DF510',  C_rs='DF498',
            C_ic='DF506',  C_is='DF508', t_oe='DF505',
        )),
    }

    def __init__(self):
        self._store: dict = {'GPS': {}, 'GAL': {}, 'BDS': {}}

    def update(self, mtype: int, parsed) -> None:
        """Ingest a pyrtcm-parsed navigation message."""
        entry = self._MAPS.get(mtype)
        if entry is None or parsed is None:
            return
        constl, sv_df, field_map = entry
        try:
            sv  = int(getattr(parsed, sv_df))
            eph = {k: float(getattr(parsed, df)) for k, df in field_map.items()}
            self._store[constl][sv] = eph
            log.debug("Ephemeris: %s SV%02d  t_oe=%.0f s", constl, sv, eph['t_oe'])
        except (AttributeError, TypeError, ValueError) as exc:
            log.debug("Ephemeris skip MSG%d: %s", mtype, exc)

    def sat_ecef(self, constl: str, sv: int, tow_s: float):
        """Return satellite ECEF (m) or None when ephemeris is unavailable."""
        eph = self._store.get(constl, {}).get(sv)
        if eph is None:
            return None
        try:
            return _kepler_ecef(eph, tow_s)
        except Exception:
            return None

    @property
    def counts(self) -> dict:
        return {k: len(v) for k, v in self._store.items()}


# ─────────────────────────────────────────────────────────────────────────────
# Extract ARP ECEF from a 1005/1006 frame (before patching)
# ─────────────────────────────────────────────────────────────────────────────

def extract_arp_ecef(raw_frame: bytes):
    """Return (X, Y, Z) metres from the ARP fields of a 1005/1006 frame."""
    p = raw_frame[3:-3]
    x = _as_signed(_get_bits(p, _BIT_X, _XYZ_BITS), _XYZ_BITS) * _XYZ_SCALE
    y = _as_signed(_get_bits(p, _BIT_Y, _XYZ_BITS), _XYZ_BITS) * _XYZ_SCALE
    z = _as_signed(_get_bits(p, _BIT_Z, _XYZ_BITS), _XYZ_BITS) * _XYZ_SCALE
    return x, y, z


# ─────────────────────────────────────────────────────────────────────────────
# MSM4/5/6/7 observation patcher
# ─────────────────────────────────────────────────────────────────────────────

# Signal-data field widths (bits) and scaling (ms / unit) per MSM subtype
_MSM_SIG_LAYOUT: dict = {
    #        fine_pr  fine_ph  lock  hc  cnr  rate   pr_scale      ph_scale
    4: dict(fine_pr=15, fine_ph=22, lock=4,  hc=1, cnr=6,  rate=0,
            pr_scale=2.0 ** -24, ph_scale=2.0 ** -29),
    5: dict(fine_pr=15, fine_ph=22, lock=4,  hc=1, cnr=6,  rate=15,
            pr_scale=2.0 ** -24, ph_scale=2.0 ** -29),
    6: dict(fine_pr=20, fine_ph=24, lock=10, hc=1, cnr=10, rate=0,
            pr_scale=2.0 ** -29, ph_scale=2.0 ** -31),
    7: dict(fine_pr=20, fine_ph=24, lock=10, hc=1, cnr=10, rate=15,
            pr_scale=2.0 ** -29, ph_scale=2.0 ** -31),
}


def patch_msm_observations(raw_frame: bytes,
                            base_xyz: tuple,
                            vbs_xyz:  tuple,
                            ephem:    EphemerisStore) -> "bytes | None":
    """
    Correct MSM4/5/6/7 pseudoranges and phase ranges for a VBS position shift.

    For each satellite with available ephemeris, computes the geometric range
    change  Δρ = −û · Δr  (û = unit vector base→satellite, Δr = vbs−base),
    then distributes it across:
      • rough range modulo 1 ms  (satellite-level, DF398, resolution ~293 m)
      • fine pseudorange residual (signal-level, DF400/DF405)
      • fine phase range residual (signal-level, DF401/DF406)

    Satellites with no ephemeris are left uncorrected.
    Returns a new frame with recalculated CRC-24Q, or **None** when no
    ephemeris is available for any satellite in this constellation (the
    caller should drop the frame entirely to avoid geometric inconsistency
    with a patched 1005/1006 base position).
    """
    payload = bytearray(raw_frame[3:-3])
    bit     = 0

    def read(n: int) -> int:
        nonlocal bit
        v = _get_bits(payload, bit, n)
        bit += n
        return v

    # ── Parse MSM header ──
    mtype       = read(12)
    constl, sub = _msm_info(mtype)
    if sub not in _MSM_SIG_LAYOUT:
        return raw_frame               # MSM1/2/3 not handled
    if not _constellation_allowed(constl):
        return None

    read(12)                           # reference station ID  (skip)
    epoch_ms = read(30)
    bit += 19          # DF393(1) + DF409(3) + DF001(7) + DF411(2) + DF412(2) + DF417(1) + DF418(3)

    tow_s    = epoch_ms / 1000.0
    sat_mask = read(64)
    sig_mask = read(32)

    nsat = bin(sat_mask).count('1')
    nsig = bin(sig_mask).count('1')

    ncells_total = nsat * nsig
    cell_mask    = read(ncells_total)
    ncell        = bin(cell_mask).count('1')

    # Satellite IDs (1-based, MSB of sat_mask = SV 1)
    sat_ids = [i + 1 for i in range(64) if (sat_mask >> (63 - i)) & 1]

    # Which satellite index does each present cell belong to?
    cell_sat_idx: list = []
    for s in range(nsat):
        for g in range(nsig):
            ci = s * nsig + g
            if (cell_mask >> (ncells_total - 1 - ci)) & 1:
                cell_sat_idx.append(s)

    # ── Satellite data ──
    rough_int_start = bit
    rough_int = [read(8) for _ in range(nsat)]   # DF397, unsigned, 255 = invalid

    if sub in (5, 7):
        bit += 4 * nsat                           # DF399 extended info (skip)

    rough_mod_start = bit
    rough_mod = [read(10) for _ in range(nsat)]  # DF398, unsigned 0-1023

    if sub in (5, 7):
        bit += 14 * nsat                          # DF419 rough phase-rate (skip)

    sig_data_start = bit    # pointer to the first signal field (fine PR)

    # ── Compute per-satellite geometric range corrections (metres) ──
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
        if dist < 1e6:             # sanity: satellite must be > 1000 km away
            corrections.append(None)
            continue
        # Δρ = −û · Δr  (positive = range increases when base moves toward sat)
        corrections.append(-(rx / dist * dx + ry / dist * dy + rz / dist * dz))

    if not any(c is not None for c in corrections):
        return None                # no ephemeris for this constellation — drop frame

    # ── Distribute correction into rough_mod + fine residual ──
    rough_int_new = list(rough_int)
    rough_mod_new = list(rough_mod)
    residuals: list = []           # metres remaining after rough adjustment

    for k in range(nsat):
        corr = corrections[k]
        if corr is None or rough_int[k] == 255:
            residuals.append(0.0)
            continue

        dt_ms  = corr / C_PER_MS
        d_mod  = round(dt_ms * 1024)    # nearest rough_mod unit (~292.8 m each)
        new_mod = rough_mod[k] + d_mod
        new_int = rough_int[k]

        while new_mod >= 1024: new_mod -= 1024; new_int += 1
        while new_mod < 0:     new_mod += 1024; new_int -= 1
        new_int = max(0, min(254, new_int))

        rough_int_new[k] = new_int
        rough_mod_new[k] = new_mod
        residuals.append((dt_ms - d_mod / 1024.0) * C_PER_MS)

    # ── Write updated rough ranges ──
    for k in range(nsat):
        _set_bits(payload, rough_int_start + k * 8,  8,  rough_int_new[k])
        _set_bits(payload, rough_mod_start + k * 10, 10, rough_mod_new[k])

    # ── Apply residuals to fine pseudoranges (all cells packed first) ──
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

    # ── Apply residuals to fine phase ranges (packed after fine PRs) ──
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


# ─────────────────────────────────────────────────────────────────────────────
# TCP broadcast server  (output side)
# ─────────────────────────────────────────────────────────────────────────────

class TcpBroadcastServer:
    """
    Simple TCP server that accepts any number of clients and broadcasts
    every call to send() to all connected clients.
    """

    def __init__(self, host: str = "0.0.0.0", port: int = 50002):
        self.host = host
        self.port = port
        self._clients: list[socket.socket] = []
        self._lock    = threading.Lock()
        self._server  = None

    def start(self):
        self._server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self._server.bind((self.host, self.port))
        self._server.listen(10)
        log.info("VBS output server listening on %s:%d", self.host, self.port)
        threading.Thread(target=self._accept_loop, daemon=True, name="vbs-accept").start()

    def _accept_loop(self):
        while True:
            try:
                conn, addr = self._server.accept()
                conn.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
                log.info("Client connected: %s:%d", addr[0], addr[1])
                with self._lock:
                    self._clients.append(conn)
            except Exception as exc:
                log.error("Accept loop error: %s", exc)
                break

    def send(self, data: bytes):
        """Broadcast data to every connected client; drop dead connections."""
        dead = []
        with self._lock:
            for cli in self._clients:
                try:
                    cli.sendall(data)
                except (BrokenPipeError, ConnectionResetError, OSError):
                    dead.append(cli)
            for cli in dead:
                log.info("Client disconnected")
                self._clients.remove(cli)

    @property
    def client_count(self) -> int:
        with self._lock:
            return len(self._clients)


# ─────────────────────────────────────────────────────────────────────────────
# Ephemeris feeder  (optional dedicated ephemeris stream)
# ─────────────────────────────────────────────────────────────────────────────

def start_ephem_feeder(host: str, port: int, ephem: EphemerisStore,
                       label: str = "EPH") -> None:
    """
    Start a background thread that connects to a dedicated RTCM ephemeris
    stream (TCP) and feeds navigation messages (1019/1042/1045/1046) into
    the shared EphemerisStore.  Reconnects automatically on failure.

    Typical setup:
        str2str -in ntrip://user:pass@caster/MOUNT_EPH -out tcpsvr://:50010
    then pass host="127.0.0.1", port=50010.
    """
    def _feeder():
        while True:
            log.info("[%s] Connecting to ephemeris source %s:%d …", label, host, port)
            src = None
            try:
                src = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                src.connect((host, port))
                log.info("[%s] Ephemeris source connected.", label)
                rdr = RTCMReader(src.makefile("rb"), quitonerror=1)
                for _raw, parsed_msg in rdr:
                    if parsed_msg is None:
                        continue
                    try:
                        mtype = int(parsed_msg.identity)
                    except (AttributeError, ValueError):
                        continue
                    if mtype in (1019, 1042, 1045, 1046):
                        ephem.update(mtype, parsed_msg)
                log.warning("[%s] Ephemeris source closed connection.", label)
            except RTCMParseError as exc:
                log.warning("[%s] RTCM parse error in ephem stream (reconnecting): %s", label, exc)
            except (ConnectionRefusedError, OSError) as exc:
                log.warning("[%s] %s — retry in 5 s", label, exc)
            finally:
                if src:
                    src.close()
            time.sleep(5)

    threading.Thread(target=_feeder, daemon=True, name=f"ephem-{label}").start()


# ─────────────────────────────────────────────────────────────────────────────
# Main processing loop
# ─────────────────────────────────────────────────────────────────────────────

def run(source_host: str, source_port: int,
        output_port: int,
        offset_north: float, offset_east: float, offset_up: float,
        ephem_host: "str | None" = None, ephem_port: "int | None" = None):

    server = TcpBroadcastServer(port=output_port)
    server.start()

    ephem                  = EphemerisStore()
    base_xyz: tuple | None = None   # real base ARP in ECEF (m)
    vbs_xyz:  tuple | None = None   # virtual base ARP in ECEF (m)
    first_1005             = True
    total = patched_pos = patched_msm = 0

    if ephem_host and ephem_port:
        start_ephem_feeder(ephem_host, ephem_port, ephem)

    while True:
        # ── (Re)connect to source ──
        log.info("Connecting to RTCM source %s:%d …", source_host, source_port)
        src = None
        try:
            src = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            src.connect((source_host, source_port))
            log.info("Connected to source.")
        except (ConnectionRefusedError, OSError) as exc:
            log.warning("Cannot connect: %s  — retrying in 5 s …", exc)
            if src:
                src.close()
            time.sleep(5)
            continue

        try:
            src_file = src.makefile("rb")
            rdr      = RTCMReader(src_file, quitonerror=1)

            for raw_frame, parsed_msg in rdr:
                if raw_frame is None:
                    continue

                total += 1
                mtype = (int(parsed_msg.identity)
                         if parsed_msg is not None
                         else _msg_type(raw_frame))

                # ── 1005/1006: patch ARP position ──
                if mtype in (1005, 1006):
                    # Always record the real ARP; only rewrite if ephemeris is
                    # available so that MSM observations can also be corrected.
                    # Patching 1005 without correcting observations creates a
                    # geometric inconsistency that degrades RTK performance.
                    base_xyz = extract_arp_ecef(raw_frame)
                    vbs_lat, vbs_lon, vbs_alt = apply_ned_offset(
                        *ecef_to_geodetic(*base_xyz),
                        offset_north, offset_east, offset_up,
                    )
                    vbs_xyz = geodetic_to_ecef(vbs_lat, vbs_lon, vbs_alt)

                    ephem_ready = any(v > 0 for v in ephem.counts.values())
                    if ephem_ready:
                        verbose   = first_1005
                        out_frame = patch_arp_position(
                            raw_frame,
                            offset_north, offset_east, offset_up,
                            verbose=verbose,
                        )
                        patched_pos += 1
                        first_1005   = False
                        if verbose and parsed_msg is not None:
                            log.info("pyrtcm parsed (original) MSG%d: %s",
                                     mtype, parsed_msg)
                    else:
                        out_frame = raw_frame
                        if first_1005:
                            log.warning(
                                "MSG%d received but no ephemeris yet — "
                                "forwarding unchanged to avoid inconsistency. "
                                "GPS/GAL/BDS MSM observations will be dropped "
                                "until ephemeris arrives. Connect an ephemeris "
                                "stream (EPHEM_PORT).",
                                mtype,
                            )

                # ── Navigation messages: accumulate ephemeris ──
                elif mtype in (1019, 1042, 1045, 1046):
                    ephem.update(mtype, parsed_msg)
                    out_frame = raw_frame

                # ── Filter unsupported constellations / auxiliaries ──
                elif mtype in (1020, 1230):
                    continue

                # ── MSM4/5/6/7: correct observations for VBS position ──
                else:
                    constl, sub = _msm_info(mtype)
                    if constl is not None and not _constellation_allowed(constl):
                        continue
                    if sub in _MSM_SIG_LAYOUT and base_xyz is not None:
                        out_frame = patch_msm_observations(
                            raw_frame, base_xyz, vbs_xyz, ephem,
                        )
                        if out_frame is None:
                            # Forwarding uncorrected observations while 1005
                            # declares a shifted VBS position causes a
                            # geometric inconsistency.
                            log.debug("Dropping %s MSM (no ephemeris) to keep stream consistent", constl)
                            continue
                        patched_msm += 1
                    else:
                        out_frame = raw_frame

                server.send(out_frame)

                if total % 500 == 0:
                    log.info(
                        "Stats — total: %d  |  1005/6 patched: %d"
                        "  |  MSM patched: %d  |  ephem: %s  |  clients: %d",
                        total, patched_pos, patched_msm,
                        ephem.counts, server.client_count,
                    )

            log.warning("Source closed connection.")

        except RTCMParseError as exc:
            log.warning("RTCM parse error (bad frame skipped, reconnecting): %s", exc)
        except OSError as exc:
            log.error("Socket error: %s", exc)
        finally:
            src.close()

        log.info("Reconnecting in 3 s …")
        time.sleep(3)


# ─────────────────────────────────────────────────────────────────────────────
# MSM full parser / builder / merger  (dual-station VBS)
# ─────────────────────────────────────────────────────────────────────────────

def parse_msm_frame(raw_frame: bytes) -> "dict | None":
    """
    Fully parse an MSM4/5/6/7 frame into a structured dict.
    Returns None if the frame is not a supported MSM subtype.

    Dict keys: mtype, constl, sub, station_id, epoch_ms, flags_bits,
    sat_mask, sig_mask, sat_ids, sig_ids, nsat, nsig,
    cell_mask_raw, cell_map [(sat_idx, sig_local_idx)…], ncell,
    rough_int, ext_info, rough_mod, rough_rate,
    fine_pr, fine_ph, lock, hc, cnr, rate_sig
    """
    payload = raw_frame[3:-3]
    bit = 0

    def read(n: int) -> int:
        nonlocal bit
        v = _get_bits(payload, bit, n)
        bit += n
        return v

    mtype       = read(12)
    constl, sub = _msm_info(mtype)
    if sub not in _MSM_SIG_LAYOUT:
        return None

    lay        = _MSM_SIG_LAYOUT[sub]
    station_id = read(12)
    epoch_ms   = read(30)
    flags_bits = read(19)   # DF393+DF409+DF001+DF411+DF412+DF417+DF418
    sat_mask   = read(64)
    sig_mask   = read(32)

    sat_ids = [i + 1 for i in range(64) if (sat_mask >> (63 - i)) & 1]
    sig_ids = [i     for i in range(32) if (sig_mask >> (31 - i)) & 1]
    nsat    = len(sat_ids)
    nsig    = len(sig_ids)

    ncells_total  = nsat * nsig
    cell_mask_raw = read(ncells_total)

    cell_map = []
    for s in range(nsat):
        for g in range(nsig):
            ci = s * nsig + g
            if (cell_mask_raw >> (ncells_total - 1 - ci)) & 1:
                cell_map.append((s, g))
    ncell = len(cell_map)

    rough_int  = [read(8)  for _ in range(nsat)]
    ext_info   = [read(4)  for _ in range(nsat)] if sub in (5, 7) else None
    rough_mod  = [read(10) for _ in range(nsat)]
    rough_rate = [_as_signed(read(14), 14) for _ in range(nsat)] if sub in (5, 7) else None

    fine_pr  = [_as_signed(read(lay['fine_pr']), lay['fine_pr']) for _ in range(ncell)]
    fine_ph  = [_as_signed(read(lay['fine_ph']), lay['fine_ph']) for _ in range(ncell)]
    lock     = [read(lay['lock'])                                 for _ in range(ncell)]
    hc       = [read(1)                                           for _ in range(ncell)]
    cnr      = [read(lay['cnr'])                                  for _ in range(ncell)]
    rate_sig = ([_as_signed(read(lay['rate']), lay['rate']) for _ in range(ncell)]
                if sub in (5, 7) else None)

    return dict(
        mtype=mtype, constl=constl, sub=sub,
        station_id=station_id, epoch_ms=epoch_ms, flags_bits=flags_bits,
        sat_mask=sat_mask, sig_mask=sig_mask,
        sat_ids=sat_ids, sig_ids=sig_ids,
        nsat=nsat, nsig=nsig,
        cell_mask_raw=cell_mask_raw, cell_map=cell_map, ncell=ncell,
        rough_int=rough_int, ext_info=ext_info,
        rough_mod=rough_mod, rough_rate=rough_rate,
        fine_pr=fine_pr, fine_ph=fine_ph,
        lock=lock, hc=hc, cnr=cnr, rate_sig=rate_sig,
    )


def _sat_avg_cnr(parsed: dict, sat_idx: int) -> float:
    """Mean CNR across all active cells for a given satellite index."""
    vals = [parsed['cnr'][c]
            for c, (si, _) in enumerate(parsed['cell_map']) if si == sat_idx]
    return sum(vals) / len(vals) if vals else 0.0


def correct_parsed_msm(p: dict,
                        base_xyz: tuple, vbs_xyz: tuple,
                        ephem: EphemerisStore) -> dict:
    """
    Apply VBS position correction to a parsed MSM dict.
    Returns a deep-copied, corrected dict (rough ranges + fine fields updated).
    """
    p   = copy.deepcopy(p)
    lay = _MSM_SIG_LAYOUT[p['sub']]
    dx  = vbs_xyz[0] - base_xyz[0]
    dy  = vbs_xyz[1] - base_xyz[1]
    dz  = vbs_xyz[2] - base_xyz[2]
    tow_s = p['epoch_ms'] / 1000.0

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


def merge_parsed_msm(ca: dict, cb: dict) -> dict:
    """
    Merge two corrected (already translated to VBS position) parsed MSM dicts
    of the same constellation and subtype.

    Strategy
    --------
    * Satellites seen only by A  → use A's observation.
    * Satellites seen only by B  → use B's observation.
    * Satellites seen by both    → use the one with higher average CNR.
    * Signal mask → union of both stations' masks.
    """
    sub = ca['sub']

    # ── Decide source for each satellite ──
    sv_source: dict = {}   # sv_id -> ('A'|'B', local_idx_in_source)
    for i, sv in enumerate(ca['sat_ids']):
        if sv in cb['sat_ids']:
            j = cb['sat_ids'].index(sv)
            sv_source[sv] = ('A', i) if _sat_avg_cnr(ca, i) >= _sat_avg_cnr(cb, j) else ('B', j)
        else:
            sv_source[sv] = ('A', i)
    for j, sv in enumerate(cb['sat_ids']):
        if sv not in sv_source:
            sv_source[sv] = ('B', j)

    merged_svs  = sorted(sv_source)
    nsat_m      = len(merged_svs)

    # ── Union signal mask ──
    sig_mask_m = ca['sig_mask'] | cb['sig_mask']
    sig_ids_m  = [i for i in range(32) if (sig_mask_m >> (31 - i)) & 1]
    nsig_m     = len(sig_ids_m)

    rough_int_m, rough_mod_m = [], []
    ext_info_m   = [] if sub in (5, 7) else None
    rough_rate_m = [] if sub in (5, 7) else None
    cell_map_m   = []
    fine_pr_m, fine_ph_m = [], []
    lock_m, hc_m, cnr_m  = [], [], []
    rate_sig_m = [] if sub in (5, 7) else None

    for sv_idx, sv in enumerate(merged_svs):
        src, src_idx = sv_source[sv]
        p = ca if src == 'A' else cb

        rough_int_m.append(p['rough_int'][src_idx])
        rough_mod_m.append(p['rough_mod'][src_idx])
        if sub in (5, 7):
            ext_info_m.append((p['ext_info']   or [0]*p['nsat'])[src_idx])
            rough_rate_m.append((p['rough_rate'] or [0]*p['nsat'])[src_idx])

        for g_idx, gsig in enumerate(sig_ids_m):
            # Translate global signal ID to source station's local sig index
            try:
                local_g = p['sig_ids'].index(gsig)
            except ValueError:
                continue   # this station doesn't have this signal type

            # Find matching cell in source station
            cell_idx = next(
                (c for c, (si, gi) in enumerate(p['cell_map'])
                 if si == src_idx and gi == local_g),
                None
            )
            if cell_idx is None:
                continue   # cell not active for this satellite/signal

            cell_map_m.append((sv_idx, g_idx))
            fine_pr_m.append(p['fine_pr'][cell_idx])
            fine_ph_m.append(p['fine_ph'][cell_idx])
            lock_m.append(p['lock'][cell_idx])
            hc_m.append(p['hc'][cell_idx])
            cnr_m.append(p['cnr'][cell_idx])
            if sub in (5, 7):
                rate_sig_m.append((p['rate_sig'] or [0]*p['ncell'])[cell_idx])

    return dict(
        mtype=ca['mtype'], constl=ca['constl'], sub=sub,
        station_id=ca['station_id'], epoch_ms=ca['epoch_ms'],
        flags_bits=ca['flags_bits'],
        sig_mask=sig_mask_m,
        sat_ids=merged_svs, sig_ids=sig_ids_m,
        nsat=nsat_m, nsig=nsig_m,
        cell_map=cell_map_m, ncell=len(cell_map_m),
        rough_int=rough_int_m, ext_info=ext_info_m,
        rough_mod=rough_mod_m, rough_rate=rough_rate_m,
        fine_pr=fine_pr_m, fine_ph=fine_ph_m,
        lock=lock_m, hc=hc_m, cnr=cnr_m, rate_sig=rate_sig_m,
    )


def build_msm_frame(parsed: dict, out_station_id: "int | None" = None) -> bytes:
    """
    Build a complete raw RTCM3 MSM frame (header + payload + CRC-24Q)
    from a parsed/merged dict as returned by parse_msm_frame() or merge_parsed_msm().
    """
    sub   = parsed['sub']
    lay   = _MSM_SIG_LAYOUT[sub]
    nsat  = parsed['nsat']
    nsig  = parsed['nsig']
    ncell = parsed['ncell']
    ncells_total = nsat * nsig

    header_bits = 12 + 12 + 30 + 19 + 64 + 32 + ncells_total
    sat_bits    = nsat * (8
                         + (4  if sub in (5, 7) else 0)
                         + 10
                         + (14 if sub in (5, 7) else 0))
    sig_bits    = ncell * (lay['fine_pr'] + lay['fine_ph'] + lay['lock'] + 1
                           + lay['cnr']
                           + (lay['rate'] if sub in (5, 7) else 0))
    total_bits  = header_bits + sat_bits + sig_bits
    total_bytes = (total_bits + 7) // 8

    payload = bytearray(total_bytes)
    bit = 0

    def write(n: int, v: int) -> None:
        nonlocal bit
        _set_bits(payload, bit, n, v)
        bit += n

    # Rebuild satellite mask from sat_ids list
    sat_mask = 0
    for sv in parsed['sat_ids']:
        sat_mask |= 1 << (64 - sv)

    # Rebuild cell mask from cell_map
    cell_set  = set(parsed['cell_map'])
    cell_mask = 0
    for s in range(nsat):
        for g in range(nsig):
            ci = s * nsig + g
            if (s, g) in cell_set:
                cell_mask |= 1 << (ncells_total - 1 - ci)

    sid = out_station_id if out_station_id is not None else parsed['station_id']

    write(12, parsed['mtype'])
    write(12, sid)
    write(30, parsed['epoch_ms'])
    write(19, parsed['flags_bits'])
    write(64, sat_mask)
    write(32, parsed['sig_mask'])
    write(ncells_total, cell_mask)

    for k in range(nsat): write(8,  parsed['rough_int'][k])
    if sub in (5, 7):
        ei = parsed['ext_info'] or [0] * nsat
        for k in range(nsat): write(4, ei[k])
    for k in range(nsat): write(10, parsed['rough_mod'][k])
    if sub in (5, 7):
        rr = parsed['rough_rate'] or [0] * nsat
        for k in range(nsat): write(14, rr[k])

    for c in range(ncell): write(lay['fine_pr'], parsed['fine_pr'][c])
    for c in range(ncell): write(lay['fine_ph'], parsed['fine_ph'][c])
    for c in range(ncell): write(lay['lock'],    parsed['lock'][c])
    for c in range(ncell): write(1,              parsed['hc'][c])
    for c in range(ncell): write(lay['cnr'],     parsed['cnr'][c])
    if sub in (5, 7):
        rs = parsed['rate_sig'] or [0] * ncell
        for c in range(ncell): write(lay['rate'], rs[c])

    header = bytes([_PREAMBLE, (total_bytes >> 8) & 0x03, total_bytes & 0xFF])
    frame  = header + bytes(payload) + b'\x00\x00\x00'
    return _recompute_crc(frame)


def set_arp_ecef(raw_frame: bytes, x: float, y: float, z: float) -> bytes:
    """Return a copy of a 1005/1006 frame with ARP ECEF set to exact (x, y, z) m."""
    payload = bytearray(raw_frame[3:-3])
    _set_bits(payload, _BIT_X, _XYZ_BITS, round(x / _XYZ_SCALE))
    _set_bits(payload, _BIT_Y, _XYZ_BITS, round(y / _XYZ_SCALE))
    _set_bits(payload, _BIT_Z, _XYZ_BITS, round(z / _XYZ_SCALE))
    _set_supported_system_flags(payload)
    new_frame = raw_frame[:3] + bytes(payload) + raw_frame[-3:]
    return _recompute_crc(new_frame)


# ─────────────────────────────────────────────────────────────────────────────
# Per-station RTCM state cache
# ─────────────────────────────────────────────────────────────────────────────

class StationBuffer:
    """
    Accumulates MSM frames and ephemeris for one base station.
    Thread-safe epoch storage; base_xyz / ephem are read under the GIL.
    """
    _MAX_EPOCHS = 10

    def __init__(self, label: str):
        self.label     = label
        self.base_xyz: "tuple | None" = None
        self.ephem     = EphemerisStore()
        self._epochs: dict = {}   # epoch_ms → {constl: raw_frame}
        self._lock     = threading.Lock()

    def add_msm(self, epoch_ms: int, constl: str, raw_frame: bytes) -> None:
        with self._lock:
            self._epochs.setdefault(epoch_ms, {})[constl] = raw_frame
            while len(self._epochs) > self._MAX_EPOCHS:
                del self._epochs[min(self._epochs)]

    def pop_epoch(self, epoch_ms: int) -> "dict | None":
        with self._lock:
            return self._epochs.pop(epoch_ms, None)

    def latest_epochs(self) -> list:
        with self._lock:
            return sorted(self._epochs)


# ─────────────────────────────────────────────────────────────────────────────
# Dual-station main processing loop
# ─────────────────────────────────────────────────────────────────────────────

def run_dual(host_a: str, port_a: int,
             host_b: str, port_b: int,
             output_port: int,
             vbs_lat: float, vbs_lon: float, vbs_alt: float,
             epoch_window_ms: float = 200.0,
             ephem_host: "str | None" = None, ephem_port: "int | None" = None):
    """
    Dual base station VBS runner.

    Connects independently to two NTRIP/TCP RTCM sources (A = primary,
    B = supplementary).  For each GPS epoch:

    * MSG 1005/1006 from A → rewritten with exact VBS ECEF and forwarded.
    * Navigation msgs (1019/1042/1045/1046) from both → stored in per-station
      EphemerisStore.
    * MSM4-7 frames:
        - constellation observed only by A → corrected A→VBS, forwarded.
        - constellation observed only by B → corrected B→VBS, forwarded.
        - constellation observed by both   → merged (CNR-wins per satellite),
          corrected to VBS, forwarded as a single frame.

    Epoch alignment: wait up to `epoch_window_ms` ms for B after A's epoch
    first appears; output with A-only data if B doesn't arrive in time.
    """
    vbs_xyz = geodetic_to_ecef(vbs_lat, vbs_lon, vbs_alt)
    log.info("VBS target : lat=%+.8f°  lon=%+.8f°  alt=%.3f m",
             vbs_lat, vbs_lon, vbs_alt)
    log.info("VBS ECEF   : X=%.4f  Y=%.4f  Z=%.4f m", *vbs_xyz)

    server = TcpBroadcastServer(port=output_port)
    server.start()

    buf_a = StationBuffer("A")
    buf_b = StationBuffer("B")

    if ephem_host and ephem_port:
        start_ephem_feeder(ephem_host, ephem_port, buf_a.ephem, label="EPH-A")
        start_ephem_feeder(ephem_host, ephem_port, buf_b.ephem, label="EPH-B")

    # epoch_ms → monotonic time when A first reported this epoch
    epoch_first_seen: dict = {}
    output_done:      set  = set()
    epoch_window_s = epoch_window_ms / 1000.0

    # ── per-station reader thread ───────────────────────────────────────────
    def reader(host: str, port: int, buf: StationBuffer,
               label: str, is_primary: bool) -> None:
        while True:
            log.info("[%s] Connecting to %s:%d …", label, host, port)
            src = None
            try:
                src = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                src.connect((host, port))
                log.info("[%s] Connected.", label)
            except (ConnectionRefusedError, OSError) as exc:
                log.warning("[%s] %s — retry in 5 s", label, exc)
                if src:
                    src.close()
                time.sleep(5)
                continue

            try:
                rdr = RTCMReader(src.makefile("rb"), quitonerror=1)
                for raw_frame, parsed_msg in rdr:
                    if raw_frame is None:
                        continue
                    mtype = (int(parsed_msg.identity)
                             if parsed_msg is not None
                             else _msg_type(raw_frame))

                    if mtype in (1005, 1006):
                        buf.base_xyz = extract_arp_ecef(raw_frame)
                        if is_primary:
                            # Forward 1005 with VBS position
                            server.send(set_arp_ecef(raw_frame, *vbs_xyz))
                            real = ecef_to_geodetic(*buf.base_xyz)
                            log.info("[A] Base: lat=%+.8f° lon=%+.8f° alt=%.3f m",
                                     *real)

                    elif mtype in (1019, 1042, 1045, 1046):
                        buf.ephem.update(mtype, parsed_msg)

                    else:
                        constl, sub = _msm_info(mtype)
                        if constl is not None and not _constellation_allowed(constl):
                            continue
                        if mtype in (1020, 1230):
                            continue
                        if sub in _MSM_SIG_LAYOUT and buf.base_xyz is not None:
                            epoch_ms = _get_bits(raw_frame[3:-3], 24, 30)
                            buf.add_msm(epoch_ms, constl, raw_frame)
                            if is_primary:
                                epoch_first_seen.setdefault(
                                    epoch_ms, time.monotonic())

                log.warning("[%s] Source closed connection.", label)
            except RTCMParseError as exc:
                log.warning("[%s] RTCM parse error (reconnecting): %s", label, exc)
            except OSError as exc:
                log.error("[%s] Socket error: %s", label, exc)
            finally:
                if src:
                    src.close()
            log.info("[%s] Reconnecting in 3 s …", label)
            time.sleep(3)

    threading.Thread(target=reader,
                     args=(host_a, port_a, buf_a, "A", True),
                     daemon=True, name="reader-A").start()
    threading.Thread(target=reader,
                     args=(host_b, port_b, buf_b, "B", False),
                     daemon=True, name="reader-B").start()

    # ── merger / output loop ────────────────────────────────────────────────
    total_out = merged_out = 0

    while True:
        time.sleep(0.05)   # 50 ms polling cadence
        now      = time.monotonic()
        ep_a_set = set(buf_a.latest_epochs())
        ep_b_set = set(buf_b.latest_epochs())

        for ep in sorted(epoch_first_seen):
            if ep in output_done:
                continue
            if ep not in ep_a_set:
                continue   # A's data already popped or not yet arrived

            has_b = ep in ep_b_set
            age   = now - epoch_first_seen[ep]

            if not has_b and age < epoch_window_s:
                continue   # still waiting for B within the window

            # ── output this epoch ─────────────────────────────────────────
            output_done.add(ep)
            data_a = buf_a.pop_epoch(ep)
            data_b = buf_b.pop_epoch(ep)

            if data_a is None:
                continue

            all_constls = set(data_a)
            if data_b:
                all_constls |= set(data_b)

            for constl in sorted(all_constls):
                fa = data_a.get(constl)
                fb = data_b.get(constl) if data_b else None

                if fa is None:
                    # B-only constellation
                    if buf_b.base_xyz:
                        out = patch_msm_observations(
                            fb, buf_b.base_xyz, vbs_xyz, buf_b.ephem)
                        if out is not None:
                            server.send(out)
                            total_out += 1
                        else:
                            log.debug("Dropping %s MSM-B (no ephemeris)", constl)

                elif fb is None:
                    # A-only constellation
                    if buf_a.base_xyz:
                        out = patch_msm_observations(
                            fa, buf_a.base_xyz, vbs_xyz, buf_a.ephem)
                        if out is not None:
                            server.send(out)
                            total_out += 1
                        else:
                            log.debug("Dropping %s MSM-A (no ephemeris)", constl)

                else:
                    # Both have this constellation → parse, correct, merge
                    try:
                        pa = parse_msm_frame(fa)
                        pb = parse_msm_frame(fb)
                        if pa and pb and pa['sub'] == pb['sub'] and buf_a.base_xyz and buf_b.base_xyz:
                            ca = correct_parsed_msm(pa, buf_a.base_xyz, vbs_xyz, buf_a.ephem)
                            cb = correct_parsed_msm(pb, buf_b.base_xyz, vbs_xyz, buf_b.ephem)
                            merged = merge_parsed_msm(ca, cb)
                            server.send(build_msm_frame(merged))
                            total_out  += 1
                            merged_out += 1
                            log.debug(
                                "Merged %s ep=%d  A:%d+B:%d→%d sats",
                                constl, ep,
                                pa['nsat'], pb['nsat'], merged['nsat'],
                            )
                        else:
                            # Sub-type mismatch or missing position — fall back to A
                            if buf_a.base_xyz:
                                out = patch_msm_observations(
                                    fa, buf_a.base_xyz, vbs_xyz, buf_a.ephem)
                                if out is not None:
                                    server.send(out)
                                    total_out += 1
                    except Exception as exc:
                        log.debug("Merge error %s ep=%d: %s — using A", constl, ep, exc)
                        if buf_a.base_xyz:
                            out = patch_msm_observations(
                                fa, buf_a.base_xyz, vbs_xyz, buf_a.ephem)
                            if out is not None:
                                server.send(out)
                                total_out += 1

        # ── periodic stats ────────────────────────────────────────────────
        if total_out > 0 and total_out % 200 == 0:
            log.info(
                "Stats — MSM out: %d  |  merged: %d"
                "  |  eph A:%s  B:%s  |  clients: %d",
                total_out, merged_out,
                buf_a.ephem.counts, buf_b.ephem.counts,
                server.client_count,
            )

        # ── housekeeping ──────────────────────────────────────────────────
        if len(output_done) > 2000:
            cutoff = (min(ep_a_set) - 120_000) if ep_a_set else 0
            output_done        = {e for e in output_done if e >= cutoff}
            epoch_first_seen   = {e: t for e, t in epoch_first_seen.items()
                                  if e >= cutoff}


# ─────────────────────────────────────────────────────────────────────────────
# ── Configuration ──
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── Mode ────────────────────────────────────────────────────────────────
    # "single" : one base station, fixed N/E/U offset  (original behaviour)
    # "dual"   : two base stations merged, explicit VBS lat/lon/alt
    MODE = "single"

    # ── Common ──────────────────────────────────────────────────────────────
    OUTPUT_PORT = 50002

    if MODE == "single":
        # TCP address of str2str output (RTCM source)
        SOURCE_HOST  = "127.0.0.1"
        SOURCE_PORT  = 50001

        # VBS offset from the real base station (metres, positive = N / E / U)
        OFFSET_NORTH = 1000
        OFFSET_EAST  = 500.0
        OFFSET_UP    = 50

        # Optional dedicated ephemeris stream — set to None to disable.
        # Start str2str separately:
        #   str2str -in ntrip://user:pass@caster/MOUNT_EPH -out tcpsvr://:50010
        EPHEM_HOST = "127.0.0.1"
        EPHEM_PORT = 50010          # set to None to disable
        # EPHEM_PORT = None         # ← uncomment to disable

        run(
            source_host  = SOURCE_HOST,
            source_port  = SOURCE_PORT,
            output_port  = OUTPUT_PORT,
            offset_north = OFFSET_NORTH,
            offset_east  = OFFSET_EAST,
            offset_up    = OFFSET_UP,
            ephem_host   = EPHEM_HOST,
            ephem_port   = EPHEM_PORT,
        )

    else:  # "dual"
        # TCP addresses of the two str2str / NTRIP outputs
        # (each str2str connects to one mount point and pushes to a local TCP port)
        #   str2str -in ntrip://caster/MOUNT_A -out tcpsvr://:50001
        #   str2str -in ntrip://caster/MOUNT_B -out tcpsvr://:50003
        HOST_A = "127.0.0.1";  PORT_A = 50001
        HOST_B = "127.0.0.1";  PORT_B = 50003

        # Target VBS position — set to the rover's approximate location
        # (WGS-84 degrees / metres above ellipsoid)
        VBS_LAT =  31.0       # degrees North  ← edit to your area
        VBS_LON = 121.0       # degrees East   ← edit to your area
        VBS_ALT =  10.0       # metres

        # How long (ms) to wait for station B before outputting A-only data
        EPOCH_WINDOW_MS = 200.0

        # Optional dedicated ephemeris stream (shared by both stations)
        EPHEM_HOST = "127.0.0.1"
        EPHEM_PORT = 50010          # set to None to disable
        # EPHEM_PORT = None         # ← uncomment to disable

        run_dual(
            host_a          = HOST_A,
            port_a          = PORT_A,
            host_b          = HOST_B,
            port_b          = PORT_B,
            output_port     = OUTPUT_PORT,
            vbs_lat         = VBS_LAT,
            vbs_lon         = VBS_LON,
            vbs_alt         = VBS_ALT,
            epoch_window_ms = EPOCH_WINDOW_MS,
            ephem_host      = EPHEM_HOST,
            ephem_port      = EPHEM_PORT,
        )
