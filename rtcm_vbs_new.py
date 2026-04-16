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

_MSM_LOG_INTERVAL_S = 15.0
_msm_log_state: dict = {}
_frame_validation_count: int = 0
_FRAME_VALIDATION_LIMIT = 5
_sat_detail_logged: set = set()


def _validate_patched_frame(original: bytes, patched: bytes, constl: str) -> None:
    """One-time validation: check that the patched frame has valid RTCM3 structure."""
    global _frame_validation_count
    if _frame_validation_count >= _FRAME_VALIDATION_LIMIT:
        return
    _frame_validation_count += 1

    if len(patched) < 6:
        log.error("[VALIDATE] %s patched frame too short: %d bytes", constl, len(patched))
        return
    if patched[0] != 0xD3:
        log.error("[VALIDATE] %s patched frame preamble wrong: 0x%02X", constl, patched[0])
        return
    declared_len = ((patched[1] & 0x03) << 8) | patched[2]
    actual_len   = len(patched) - 6  # minus 3-byte header and 3-byte CRC
    if declared_len != actual_len:
        log.error("[VALIDATE] %s length mismatch: declared=%d actual=%d",
                  constl, declared_len, actual_len)
        return
    crc_check = crc24q(patched[:-3])
    crc_frame = struct.unpack(">I", b"\x00" + patched[-3:])[0]
    if crc_check != crc_frame:
        log.error("[VALIDATE] %s CRC mismatch: computed=0x%06X frame=0x%06X",
                  constl, crc_check, crc_frame)
        return
    if len(original) == len(patched):
        diff_bytes = sum(1 for a, b in zip(original, patched) if a != b)
        log.info("[VALIDATE] %s frame OK: %d bytes, %d bytes changed (in-place patch)",
                 constl, len(patched), diff_bytes)
    else:
        log.info("[VALIDATE] %s frame OK: orig=%d bytes, patched=%d bytes",
                 constl, len(original), len(patched))


def _has_any_ephemeris(ephem: "EphemerisStore") -> bool:
    """True once any supported constellation has at least one ephemeris."""
    return any(v > 0 for v in ephem.counts.values())


def _maybe_log_msm_correction_status(
        mtype: int,
        constl: str,
        sub: int,
        valid_corr: int,
        nsat_total: int,
        nsat_out: int,
        ncell_out: int,
        eph_cached: int) -> None:
    """Log MSM correction coverage with simple rate limiting."""
    key = (mtype, constl, sub)
    now = time.monotonic()
    state = (valid_corr, nsat_total, nsat_out, ncell_out, eph_cached)
    prev = _msm_log_state.get(key)

    if prev is not None:
        prev_state, prev_time = prev
        if prev_state == state and (now - prev_time) < _MSM_LOG_INTERVAL_S:
            return

    level = logging.INFO if valid_corr > 0 else logging.WARNING
    log.log(
        level,
        "MSG%d %s MSM%d: valid corrections %d/%d sats, output sats=%d, active cells=%d, cached ephem=%d",
        mtype, constl, sub, valid_corr, nsat_total, nsat_out, ncell_out, eph_cached,
    )
    _msm_log_state[key] = (state, now)

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


def _signed_max(bits: int) -> int:
    return (1 << (bits - 1)) - 1


def _invalid_msm_fine(bits: int) -> int:
    return -(1 << (bits - 1))


def _clamp_msm_fine(value: int, bits: int) -> int:
    # The most-negative code is the MSM invalid sentinel, so clamp above it.
    return max(_invalid_msm_fine(bits) + 1, min(_signed_max(bits), value))


def _set_supported_system_flags(payload: bytearray) -> None:
    # Preserve original DF022-DF024 flags from the base station.
    # DF022=GPS, DF023=GLONASS/secondary, DF024=Galileo.
    # Only clear GLONASS flag since we filter GLONASS MSM;
    # leave GPS and Galileo flags at their original values.
    _set_bits(payload, _BIT_SYS, 1, 0)  # DF023: no GLONASS


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
# GNSS constants  (matching RTKLIB ephemeris.c / rtklib.h exactly)
# ─────────────────────────────────────────────────────────────────────────────
C_LIGHT  = 299_792_458.0              # speed of light (m/s)  — rtklib.h
C_PER_MS = C_LIGHT * 1e-3            # = RANGE_MS in rtcm3.c  (299 792.458 m/ms)

# Gravitational parameters  — ephemeris.c
MU_GPS   = 3.9860050e14              # GPS   (m³/s²)
MU_GAL   = 3.986004418e14            # Galileo
MU_CMP   = 3.986004418e14            # BeiDou

# Earth rotation rates  — rtklib.h + ephemeris.c
OMGE_GPS = 7.2921151467e-5           # GPS / Galileo  (rad/s)
OMGE_CMP = 7.292115e-5               # BeiDou         (rad/s)

# Per-constellation lookup (mu, omge)  — mirrors ephemeris.c:258–261
_CONSTL_PARAMS: dict[str, tuple[float, float]] = {
    "GPS": (MU_GPS,  OMGE_GPS),
    "GAL": (MU_GAL,  OMGE_GPS),   # Galileo uses same OMGE as GPS
    "BDS": (MU_CMP,  OMGE_CMP),
}

# MSM constellation bases — message number = base + subtype (subtype 1-7)
_MSM_CONSTL = (
    ('GPS', 1070), ('GLO', 1080), ('GAL', 1090),
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
#   — mirrors RTKLIB ephemeris.c eph2pos() exactly
# ─────────────────────────────────────────────────────────────────────────────

# BDS GEO: prn<=5 || prn>=59  — ephemeris.c:285
_BDS_GEO_PRNS = frozenset(range(1, 6)) | frozenset(range(59, 65))
# sin/cos(-5°) constants from ephemeris.c:94-95
_SIN_5 = -0.0871557427476582
_COS_5 =  0.9961946980917456


def _kepler_ecef(eph: dict, tow_s: float, constl: str = "GPS", sv: int = 0):
    """Satellite ECEF (m) from broadcast Keplerian elements.

    Matches RTKLIB eph2pos() including:
    - Constellation-specific mu/omge (ephemeris.c:258-261)
    - BDS GEO special orbit (prn<=5||prn>=59) (ephemeris.c:284-301)
    """
    mu, omge = _CONSTL_PARAMS.get(constl, (MU_GPS, OMGE_GPS))

    dt = tow_s - eph['t_oe']
    if dt >  302_400: dt -= 604_800
    if dt < -302_400: dt += 604_800

    a = eph['sqrtA'] ** 2
    n = math.sqrt(mu / a ** 3) + eph['delta_n']
    M = eph['M0'] + n * dt

    e, E = eph['e'], float(M)
    for _ in range(12):
        E -= (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))

    sinE, cosE = math.sin(E), math.cos(E)
    v   = math.atan2(math.sqrt(1.0 - e * e) * sinE, cosE - e)
    phi = v + eph['omega']
    s2, c2 = math.sin(2 * phi), math.cos(2 * phi)

    u = phi + eph['C_us'] * s2 + eph['C_uc'] * c2
    r = a * (1.0 - e * cosE) + eph['C_rs'] * s2 + eph['C_rc'] * c2
    i = eph['i0'] + eph['IDOT'] * dt + eph['C_is'] * s2 + eph['C_ic'] * c2

    x, y = r * math.cos(u), r * math.sin(u)
    cosi = math.cos(i)

    is_bds_geo = (constl == "BDS" and sv in _BDS_GEO_PRNS)

    if is_bds_geo:
        # ephemeris.c:285-294  — BeiDou GEO
        O = eph['Omega0'] + eph['OmegaDot'] * dt - omge * eph['t_oe']
        sinO, cosO = math.sin(O), math.cos(O)
        xg = x * cosO - y * cosi * sinO
        yg = x * sinO + y * cosi * cosO
        zg = y * math.sin(i)
        sino = math.sin(omge * dt)
        coso = math.cos(omge * dt)
        return ( xg * coso + yg * sino * _COS_5 + zg * sino * _SIN_5,
                -xg * sino + yg * coso * _COS_5 + zg * coso * _SIN_5,
                            -yg * _SIN_5         + zg * _COS_5)
    else:
        # ephemeris.c:296-301  — GPS / GAL / BDS MEO-IGSO
        O = eph['Omega0'] + (eph['OmegaDot'] - omge) * dt - omge * eph['t_oe']
        sinO, cosO = math.sin(O), math.cos(O)
        return (x * cosO - y * cosi * sinO,
                x * sinO + y * cosi * cosO,
                y * math.sin(i))


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
            return _kepler_ecef(eph, tow_s, constl=constl, sv=sv)
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


def patch_msm_inplace(raw_frame: bytes,
                      base_xyz: tuple,
                      vbs_xyz: tuple,
                      ephem: EphemerisStore) -> "bytes | None":
    """
    In-place patch of MSM4/5/6/7 rough-range + fine-PR/phase fields.

    Instead of rebuilding the frame from scratch (which risks encoding
    mismatches with rtkrcv), this modifies only the range-related bit fields
    directly in a copy of the original frame bytes, preserving all other
    structure, padding, and reserved bits.

    Returns a new frame (bytes) with corrections applied and CRC recomputed,
    or None if the constellation is unsupported or has zero ephemeris coverage.
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
        return raw_frame
    if not _constellation_allowed(constl):
        return None

    lay        = _MSM_SIG_LAYOUT[sub]
    _station_id = read(12)
    epoch_ms    = read(30)
    _flags      = read(19)
    sat_mask    = read(64)
    sig_mask    = read(32)

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

    # Record bit offsets for satellite-level data
    rough_int_offsets = []
    rough_int_vals    = []
    for k in range(nsat):
        rough_int_offsets.append(bit)
        rough_int_vals.append(read(8))

    if sub in (5, 7):
        for _ in range(nsat): read(4)  # ext_info — skip

    rough_mod_offsets = []
    rough_mod_vals    = []
    for k in range(nsat):
        rough_mod_offsets.append(bit)
        rough_mod_vals.append(read(10))

    if sub in (5, 7):
        for _ in range(nsat): read(14)  # rough_rate — skip

    # Record bit offsets for signal-level data
    fine_pr_offsets, fine_pr_vals = [], []
    for c in range(ncell):
        fine_pr_offsets.append(bit)
        fine_pr_vals.append(_as_signed(read(lay['fine_pr']), lay['fine_pr']))

    fine_ph_offsets, fine_ph_vals = [], []
    for c in range(ncell):
        fine_ph_offsets.append(bit)
        fine_ph_vals.append(_as_signed(read(lay['fine_ph']), lay['fine_ph']))

    # ── Compute geometric corrections (exact range difference, no linear approx) ──
    tow_s = (epoch_ms / 1000.0) % 604_800.0

    corrections: list[float | None] = []
    valid_count = 0
    for sv in sat_ids:
        pos = ephem.sat_ecef(constl, sv, tow_s)
        if pos is None:
            corrections.append(None)
            continue
        # Exact range from base and VBS to satellite — no Taylor expansion error
        dist_base = math.sqrt((pos[0] - base_xyz[0]) ** 2 +
                              (pos[1] - base_xyz[1]) ** 2 +
                              (pos[2] - base_xyz[2]) ** 2)
        dist_vbs  = math.sqrt((pos[0] - vbs_xyz[0]) ** 2 +
                              (pos[1] - vbs_xyz[1]) ** 2 +
                              (pos[2] - vbs_xyz[2]) ** 2)
        if dist_base < 1e6:
            corrections.append(None)
            continue
        corrections.append(dist_vbs - dist_base)
        valid_count += 1

    _maybe_log_msm_correction_status(
        mtype, constl, sub,
        valid_count, nsat, nsat, ncell,
        ephem.counts.get(constl, 0),
    )

    if valid_count == 0:
        return None

    # One-time per-satellite detail log (first epoch per constellation)
    if constl not in _sat_detail_logged:
        _sat_detail_logged.add(constl)
        for k, sv in enumerate(sat_ids):
            corr = corrections[k]
            ri = rough_int_vals[k]
            rm = rough_mod_vals[k]
            rng_m = (ri + rm / 1024.0) * C_PER_MS if ri != 255 else 0.0
            if corr is not None:
                log.info("[DETAIL] %s SV%02d: rough_int=%3d rough_mod=%4d "
                         "range=%.3f m  corr=%+.4f m",
                         constl, sv, ri, rm, rng_m, corr)
            else:
                log.warning("[DETAIL] %s SV%02d: rough_int=%3d  NO EPHEMERIS → invalidated",
                            constl, sv, ri)

    # ── Apply corrections in-place ──
    out = bytearray(raw_frame)
    out_payload_start = 3
    def write_payload(offset: int, n: int, val: int) -> None:
        _set_bits(out, out_payload_start * 8 + offset, n, val)

    residuals: list[float] = []
    for k in range(nsat):
        corr = corrections[k]
        if rough_int_vals[k] == 255:
            residuals.append(0.0)
            continue
        if corr is None:
            # No ephemeris → mark satellite invalid so rtkrcv discards it.
            # rtcm3.c decode_msm4: rng==255 → r[j]=0.0 → obs not stored.
            write_payload(rough_int_offsets[k], 8, 255)
            residuals.append(0.0)
            continue
        dt_ms = corr / C_PER_MS
        d_mod = round(dt_ms * 1024)
        nm = rough_mod_vals[k] + d_mod
        ni = rough_int_vals[k]
        while nm >= 1024: nm -= 1024; ni += 1
        while nm <     0: nm += 1024; ni -= 1
        if ni < 0 or ni > 254:
            write_payload(rough_int_offsets[k], 8, 255)
            residuals.append(0.0)
            continue
        write_payload(rough_int_offsets[k], 8, ni)
        write_payload(rough_mod_offsets[k], 10, nm)
        residuals.append((dt_ms - d_mod / 1024.0) * C_PER_MS)

    for c, (si, _gi) in enumerate(cell_map):
        res = residuals[si]
        if res == 0.0:
            continue
        if fine_pr_vals[c] != _invalid_msm_fine(lay['fine_pr']):
            new_pr = _clamp_msm_fine(
                fine_pr_vals[c] + round(res / (C_PER_MS * lay['pr_scale'])),
                lay['fine_pr'],
            )
            write_payload(fine_pr_offsets[c], lay['fine_pr'], new_pr)
        if fine_ph_vals[c] != _invalid_msm_fine(lay['fine_ph']):
            new_ph = _clamp_msm_fine(
                fine_ph_vals[c] + round(res / (C_PER_MS * lay['ph_scale'])),
                lay['fine_ph'],
            )
            write_payload(fine_ph_offsets[c], lay['fine_ph'], new_ph)

    result = bytes(_recompute_crc(bytes(out)))

    # ── Round-trip verification: decode patched frame as RTKLIB would ──
    if _frame_validation_count < _FRAME_VALIDATION_LIMIT:
        _roundtrip_verify(
            raw_frame, result, payload, lay,
            rough_int_offsets, rough_mod_offsets,
            fine_pr_offsets, fine_ph_offsets,
            fine_pr_vals, fine_ph_vals,
            rough_int_vals, rough_mod_vals,
            corrections, sat_ids, cell_map,
            nsat, constl,
        )

    return result


def _roundtrip_verify(
    raw_frame, patched_frame, orig_payload, lay,
    rough_int_offsets, rough_mod_offsets,
    fine_pr_offsets, fine_ph_offsets,
    orig_fine_pr, orig_fine_ph,
    orig_rough_int, orig_rough_mod,
    corrections, sat_ids, cell_map,
    nsat, constl,
):
    """One-time round-trip verification: decode total observations from patched
    frame and compare with expected (original + correction).

    Only the TOTAL observation (rough + fine) is checked — comparing rough range
    alone is meaningless because it changes in ~293 m steps, with the fine fields
    absorbing the remainder.  Phase fields that carry the invalid sentinel in the
    original frame are reported once for awareness but are NOT counted as errors.
    """
    global _frame_validation_count
    if _frame_validation_count > _FRAME_VALIDATION_LIMIT:
        return
    patched_payload = patched_frame[3:-3]

    pr_sentinel  = _invalid_msm_fine(lay['fine_pr'])
    ph_sentinel  = _invalid_msm_fine(lay['fine_ph'])

    errors_found = False
    invalid_sats = 0
    skipped_ph   = 0

    for k in range(nsat):
        sv   = sat_ids[k]
        corr = corrections[k]
        new_ri = _get_bits(patched_payload, rough_int_offsets[k], 8)
        if corr is None and orig_rough_int[k] != 255:
            if new_ri != 255:
                log.error("[ROUNDTRIP] %s SV%02d: expected invalidated (255) "
                          "but got rough_int=%d", constl, sv, new_ri)
                errors_found = True
            else:
                invalid_sats += 1

    for c, (si, _gi) in enumerate(cell_map):
        sv   = sat_ids[si]
        corr = corrections[si]
        if corr is None or orig_rough_int[si] == 255:
            continue

        new_fpr = _as_signed(
            _get_bits(patched_payload, fine_pr_offsets[c], lay['fine_pr']),
            lay['fine_pr'])
        new_fph = _as_signed(
            _get_bits(patched_payload, fine_ph_offsets[c], lay['fine_ph']),
            lay['fine_ph'])

        ori_ri = orig_rough_int[si]
        ori_rm = orig_rough_mod[si]
        new_ri = _get_bits(patched_payload, rough_int_offsets[si], 8)
        new_rm = _get_bits(patched_payload, rough_mod_offsets[si], 10)

        # --- Pseudorange (total = rough + fine_pr) ---
        pr_valid = (orig_fine_pr[c] != pr_sentinel)
        if pr_valid:
            orig_total_pr = ((ori_ri + ori_rm / 1024.0) * C_PER_MS
                             + orig_fine_pr[c] * lay['pr_scale'] * C_PER_MS)
            new_total_pr  = ((new_ri + new_rm / 1024.0) * C_PER_MS
                             + new_fpr * lay['pr_scale'] * C_PER_MS)
            pr_err = new_total_pr - (orig_total_pr + corr)
            if abs(pr_err) > 0.02:
                log.error("[ROUNDTRIP] %s SV%02d sig%d: PR total err=%+.4f m "
                          "(orig=%.3f new=%.3f expected=%.3f)",
                          constl, sv, _gi, pr_err,
                          orig_total_pr, new_total_pr, orig_total_pr + corr)
                errors_found = True

        # --- Carrier phase (total = rough + fine_ph) ---
        ph_valid = (orig_fine_ph[c] != ph_sentinel)
        if ph_valid:
            orig_total_ph = ((ori_ri + ori_rm / 1024.0) * C_PER_MS
                             + orig_fine_ph[c] * lay['ph_scale'] * C_PER_MS)
            new_total_ph  = ((new_ri + new_rm / 1024.0) * C_PER_MS
                             + new_fph * lay['ph_scale'] * C_PER_MS)
            ph_err = new_total_ph - (orig_total_ph + corr)
            if abs(ph_err) > 0.002:
                log.error("[ROUNDTRIP] %s SV%02d sig%d: PH total err=%+.6f m "
                          "(orig=%.3f new=%.3f expected=%.3f)",
                          constl, sv, _gi, ph_err,
                          orig_total_ph, new_total_ph, orig_total_ph + corr)
                errors_found = True
        else:
            skipped_ph += 1

    summary_parts = []
    if invalid_sats:
        summary_parts.append("%d sats invalidated (no ephem)" % invalid_sats)
    if skipped_ph:
        summary_parts.append("%d signals with original invalid phase (sentinel)" % skipped_ph)
    extra = (" — " + ", ".join(summary_parts)) if summary_parts else ""

    if not errors_found:
        log.info("[ROUNDTRIP] %s: %d sats × %d cells verified OK "
                 "(PR<0.02m, PH<0.002m)%s",
                 constl, nsat, len(cell_map), extra)
    else:
        log.warning("[ROUNDTRIP] %s: encoding errors detected%s",
                    constl, extra)


def patch_msm_observations(raw_frame: bytes,
                            base_xyz: tuple,
                            vbs_xyz:  tuple,
                            ephem:    EphemerisStore) -> "bytes | None":
    """
    Correct MSM4/5/6/7 pseudoranges and phase ranges for a VBS position shift.

    Uses in-place patching to preserve the original frame structure exactly,
    only modifying rough-range and fine-PR/phase fields.  This avoids any
    rebuild-related encoding issues that could prevent rtkrcv from decoding
    the frame.

    Returns a new frame with recalculated CRC-24Q, or **None** when no
    ephemeris is available for any satellite in this constellation.
    """
    return patch_msm_inplace(raw_frame, base_xyz, vbs_xyz, ephem)


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
                src.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
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
    ephem_ready_once       = False
    position_patch_active  = False
    deferred_1005          = 0
    startup_msm_logged: set[str] = set()
    total = patched_pos = patched_msm = 0

    if ephem_host and ephem_port:
        start_ephem_feeder(ephem_host, ephem_port, ephem)

    while True:
        # ── (Re)connect to source ──
        log.info("Connecting to RTCM source %s:%d …", source_host, source_port)
        src = None
        try:
            src = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            src.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
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
                    # Always record the real ARP. During startup, defer the VBS
                    # position patch until ephemeris is available so that raw
                    # 1005/1006 and raw MSM stay geometrically consistent.
                    base_xyz = extract_arp_ecef(raw_frame)
                    vbs_lat, vbs_lon, vbs_alt = apply_ned_offset(
                        *ecef_to_geodetic(*base_xyz),
                        offset_north, offset_east, offset_up,
                    )
                    vbs_xyz = geodetic_to_ecef(vbs_lat, vbs_lon, vbs_alt)

                    if not ephem_ready_once and not _has_any_ephemeris(ephem):
                        out_frame = raw_frame
                        deferred_1005 += 1
                        if deferred_1005 == 1 or deferred_1005 % 20 == 0:
                            log.warning(
                                "MSG%d: startup phase, no ephemeris yet; deferring "
                                "1005/1006 patch (%d deferred, ephem=%s)",
                                mtype, deferred_1005, ephem.counts,
                            )
                    else:
                        ephem_ready_once = True
                        verbose   = first_1005
                        out_frame = patch_arp_position(
                            raw_frame,
                            offset_north, offset_east, offset_up,
                            verbose=verbose,
                        )
                        patched_pos += 1
                        first_1005  = False
                        position_patch_active = True
                        if verbose:
                            log.info(
                                "Ephemeris ready; enabling VBS 1005/1006 patch "
                                "(deferred startup frames: %d, ephem=%s)",
                                deferred_1005, ephem.counts,
                            )
                            if parsed_msg is not None:
                                log.info("pyrtcm parsed (original) MSG%d: %s",
                                         mtype, parsed_msg)

                # ── Navigation messages: accumulate ephemeris ──
                elif mtype in (1019, 1042, 1045, 1046):
                    ephem.update(mtype, parsed_msg)
                    out_frame = raw_frame
                    if not ephem_ready_once and _has_any_ephemeris(ephem):
                        ephem_ready_once = True
                        log.info(
                            "Ephemeris ready after startup; VBS patching will start "
                            "from the next 1005/1006 frame (deferred 1005/1006: %d, ephem=%s)",
                            deferred_1005, ephem.counts,
                        )

                # ── Filter unsupported constellations / auxiliaries ──
                # GLONASS ephemeris, bias, and legacy observation messages
                # must be blocked — they can't be geometrically corrected and
                # would be inconsistent with the patched VBS base position.
                elif mtype in (1020, 1230,
                               1001, 1002, 1003, 1004,
                               1009, 1010, 1011, 1012):
                    continue

                # ── MSM4/5/6/7: correct observations for VBS position ──
                else:
                    constl, sub = _msm_info(mtype)
                    if constl is not None and not _constellation_allowed(constl):
                        continue
                    if sub in _MSM_SIG_LAYOUT and base_xyz is not None:
                        if not position_patch_active:
                            out_frame = raw_frame
                            if constl not in startup_msm_logged:
                                log.info(
                                    "Startup phase: forwarding raw %s MSM until first patched 1005/1006 is sent",
                                    constl,
                                )
                                startup_msm_logged.add(constl)
                        else:
                            out_frame = patch_msm_observations(
                                raw_frame, base_xyz, vbs_xyz, ephem,
                            )
                            if out_frame is None:
                                log.debug("Dropping %s MSM (no ephemeris)", constl)
                                continue
                            _validate_patched_frame(raw_frame, out_frame, constl or "?")
                            patched_msm += 1
                    else:
                        out_frame = raw_frame

                server.send(out_frame)

                if total % 500 == 0:
                    log.info(
                        "Stats — total: %d  |  1005/6 patched: %d"
                        "  |  1005/6 deferred: %d  |  MSM patched: %d"
                        "  |  ephem ready: %s  |  ephem: %s  |  clients: %d",
                        total, patched_pos, deferred_1005, patched_msm,
                        ephem_ready_once,
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
    tow_s = (p['epoch_ms'] / 1000.0) % 604_800.0   # Fix-2: modulo one GPS week

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
        if p['fine_pr'][c] != _invalid_msm_fine(lay['fine_pr']):
            p['fine_pr'][c] = _clamp_msm_fine(
                p['fine_pr'][c] + round(res / (C_PER_MS * lay['pr_scale'])),
                lay['fine_pr'],
            )
        if p['fine_ph'][c] != _invalid_msm_fine(lay['fine_ph']):
            p['fine_ph'][c] = _clamp_msm_fine(
                p['fine_ph'][c] + round(res / (C_PER_MS * lay['ph_scale'])),
                lay['fine_ph'],
            )

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
             epoch_window_ms: float = 40.0,
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
                src.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
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
                        if mtype in (1020, 1230,
                                     1001, 1002, 1003, 1004,
                                     1009, 1010, 1011, 1012):
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
        time.sleep(0.005)   # 5 ms polling cadence to keep VBS latency low
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
        OFFSET_NORTH = 10
        OFFSET_EAST  = 10
        OFFSET_UP    = 10

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
        VBS_LAT = 22.937764       # degrees North  ← edit to your area
        VBS_LON = 113.206806       # degrees East   ← edit to your area
        VBS_ALT = 120.253       # metres

        # How long (ms) to wait for station B before outputting A-only data
        EPOCH_WINDOW_MS = 40.0

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
