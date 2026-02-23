#!/usr/bin/env python3
"""
BlackRoad Marine Biology Lab - Scientific Analysis Module
Real marine science algorithms for field research and lab analysis.
"""

import sqlite3
import math
import argparse
from dataclasses import dataclass, field
from typing import Optional
from datetime import datetime

# ANSI ocean colour palette
OCEAN_DEEP = '\033[38;5;21m'
OCEAN_BLUE = '\033[38;5;39m'
CORAL      = '\033[38;5;209m'
SEAFOAM    = '\033[38;5;122m'
SAND       = '\033[38;5;222m'
KELP       = '\033[38;5;64m'
BIOLUM     = '\033[38;5;51m'
BOLD       = '\033[1m'
RESET      = '\033[0m'

DB_PATH = "marine_lab.db"

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SpeciesObservation:
    species_name: str
    common_name: str
    count: int
    location_lat: float
    location_lon: float
    depth_m: float
    water_temp_c: float
    trophic_level: float = 2.0
    observer: str = "unknown"
    notes: str = ""
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())


@dataclass
class ReefHealthScore:
    reef_id: str
    coral_cover_pct: float        # 0–100 %
    bleaching_pct: float          # 0–100 % (penalty)
    fish_density_per_100m2: float
    water_clarity_m: float        # Secchi depth
    macroalgae_cover_pct: float   # 0–100 % (penalty)
    invertebrate_density: float
    score: float = 0.0
    grade: str = ""
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())


@dataclass
class WaterQualityIndex:
    station_id: str
    dissolved_oxygen_mgl: float
    ph: float
    turbidity_ntu: float
    temperature_c: float
    nitrate_mgl: float
    phosphate_mgl: float
    salinity_ppt: float
    wqi: float = 0.0
    category: str = ""
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())


@dataclass
class PopulationModel:
    species: str
    initial_population: float
    carrying_capacity: float
    growth_rate: float
    time_steps: int
    prey_species: Optional[str] = None
    predator_species: Optional[str] = None


@dataclass
class BiodiversityMetrics:
    site_id: str
    species_counts: dict
    shannon_h: float = 0.0
    simpson_d: float = 0.0
    chao1: float = 0.0
    ace: float = 0.0
    species_richness: int = 0
    evenness: float = 0.0
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())


@dataclass
class MigrationRecord:
    species: str
    tag_id: str
    origin_lat: float
    origin_lon: float
    destination_lat: float
    destination_lon: float
    duration_days: float
    distance_km: float = 0.0
    avg_speed_kmday: float = 0.0
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------

def init_db(db_path: str = DB_PATH) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS observations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            species_name TEXT, common_name TEXT, count INTEGER,
            location_lat REAL, location_lon REAL, depth_m REAL,
            water_temp_c REAL, trophic_level REAL, observer TEXT,
            notes TEXT, timestamp TEXT
        );
        CREATE TABLE IF NOT EXISTS reef_health (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            reef_id TEXT, coral_cover_pct REAL, bleaching_pct REAL,
            fish_density_per_100m2 REAL, water_clarity_m REAL,
            macroalgae_cover_pct REAL, invertebrate_density REAL,
            score REAL, grade TEXT, timestamp TEXT
        );
        CREATE TABLE IF NOT EXISTS water_quality (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            station_id TEXT, dissolved_oxygen_mgl REAL, ph REAL,
            turbidity_ntu REAL, temperature_c REAL,
            nitrate_mgl REAL, phosphate_mgl REAL,
            salinity_ppt REAL, wqi REAL, category TEXT, timestamp TEXT
        );
        CREATE TABLE IF NOT EXISTS migrations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            species TEXT, tag_id TEXT,
            origin_lat REAL, origin_lon REAL,
            destination_lat REAL, destination_lon REAL,
            duration_days REAL, distance_km REAL,
            avg_speed_kmday REAL, timestamp TEXT
        );
    """)
    conn.commit()
    return conn


# ---------------------------------------------------------------------------
# Scientific Algorithms
# ---------------------------------------------------------------------------

def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance (km) via Haversine formula."""
    R = 6371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi    = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return R * 2 * math.asin(math.sqrt(a))


def dissolved_oxygen_saturation_mgl(temp_c: float, salinity_ppt: float) -> float:
    """
    Garcia & Gordon (1992) dissolved-oxygen saturation equation.
    Returns DO saturation in mg/L at 1 atm.

      Ts = ln((298.15 - T) / (273.15 + T))
      ln(O2_sat) = A0 + A1·Ts + A2·Ts² + A3·Ts³ + A4·Ts⁴ + A5·Ts⁵
                 + S·(B0 + B1·Ts + B2·Ts² + B3·Ts³) + C0·S²
    """
    A0, A1, A2, A3, A4, A5 =  2.00907,  3.22014,  4.05010,  4.94457, -0.256847,  3.88767
    B0, B1, B2, B3          = -0.00624523, -0.00737614, -0.010341, -0.00817083
    C0                      = -4.88682e-7
    Ts = math.log((298.15 - temp_c) / (273.15 + temp_c))
    ln_do = (A0 + A1*Ts + A2*Ts**2 + A3*Ts**3 + A4*Ts**4 + A5*Ts**5
             + salinity_ppt * (B0 + B1*Ts + B2*Ts**2 + B3*Ts**3)
             + C0 * salinity_ppt**2)
    return math.exp(ln_do) * 1.4276   # ml/L → mg/L


def shannon_wiener_index(species_counts: dict) -> float:
    """Shannon-Wiener H' = -Σ pi·ln(pi).  Returns H' in nats."""
    total = sum(species_counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for n in species_counts.values():
        if n > 0:
            pi = n / total
            h -= pi * math.log(pi)
    return h


def simpsons_diversity_index(species_counts: dict) -> float:
    """Simpson's D = 1 - Σ ni(ni-1) / N(N-1).  Range [0, 1]."""
    counts = [v for v in species_counts.values() if v > 0]
    N = sum(counts)
    if N <= 1:
        return 0.0
    return 1.0 - sum(n * (n - 1) for n in counts) / (N * (N - 1))


def pielou_evenness(shannon_h: float, species_richness: int) -> float:
    """Pielou's J = H' / ln(S).  Range [0, 1]."""
    if species_richness <= 1:
        return 0.0
    return shannon_h / math.log(species_richness)


def chao1_estimator(species_counts: dict) -> float:
    """
    Chao1 non-parametric richness estimator.
      f1 = singletons, f2 = doubletons
      Ŝ = Sobs + f1²/(2·f2)   [if f2 > 0]
      Ŝ = Sobs + f1·(f1-1)/2  [if f2 = 0]
    """
    sobs = sum(1 for v in species_counts.values() if v > 0)
    f1   = sum(1 for v in species_counts.values() if v == 1)
    f2   = sum(1 for v in species_counts.values() if v == 2)
    if f2 > 0:
        return sobs + f1**2 / (2.0 * f2)
    if f1 > 0:
        return sobs + f1 * (f1 - 1) / 2.0
    return float(sobs)


def ace_estimator(species_counts: dict, rare_threshold: int = 10) -> float:
    """
    ACE (Abundance-based Coverage Estimator) species richness.
    Splits community into rare (ni ≤ threshold) and abundant sets.
      Ĉ_ACE = 1 - f1/n_rare
      Ŝ_ACE = S_abund + S_rare/Ĉ_ACE + (f1/Ĉ_ACE)·γ²_ACE
    """
    rare  = {k: v for k, v in species_counts.items() if 0 < v <= rare_threshold}
    abund = {k: v for k, v in species_counts.items() if v > rare_threshold}
    S_rare, S_abund = len(rare), len(abund)
    n_rare = sum(rare.values())
    f1 = sum(1 for v in rare.values() if v == 1)
    if n_rare == 0 or S_rare == 0:
        return float(len([v for v in species_counts.values() if v > 0]))
    C_ACE = max(1e-9, 1.0 - f1 / n_rare)
    gamma_sq = max(0.0, (S_rare / C_ACE) *
                   (sum(v * (v - 1) for v in rare.values()) / (n_rare * (n_rare - 1))) - 1.0)
    return S_abund + S_rare / C_ACE + (f1 / C_ACE) * gamma_sq


def compute_biodiversity(site_id: str, species_counts: dict) -> BiodiversityMetrics:
    """Compute full suite of biodiversity metrics."""
    valid = {k: v for k, v in species_counts.items() if v > 0}
    h = shannon_wiener_index(valid)
    s = len(valid)
    return BiodiversityMetrics(
        site_id=site_id, species_counts=valid,
        shannon_h=round(h, 4),
        simpson_d=round(simpsons_diversity_index(valid), 4),
        chao1=round(chao1_estimator(valid), 2),
        ace=round(ace_estimator(valid), 2),
        species_richness=s,
        evenness=round(pielou_evenness(h, s), 4),
    )


def score_reef_health(r: ReefHealthScore) -> ReefHealthScore:
    """
    Composite reef health score (0–100).
    Weights: coral cover 30 %, bleaching 25 % (penalty), fish 20 %,
             clarity 15 %, macroalgae 10 % (penalty).
    """
    coral   = min(r.coral_cover_pct / 50.0, 1.0) * 100
    bleach  = max(0.0, 1.0 - r.bleaching_pct / 100.0) * 100
    fish    = min(r.fish_density_per_100m2 / 50.0, 1.0) * 100
    clarity = min(r.water_clarity_m / 20.0, 1.0) * 100
    algae   = max(0.0, 1.0 - r.macroalgae_cover_pct / 100.0) * 100
    r.score = round(0.30*coral + 0.25*bleach + 0.20*fish + 0.15*clarity + 0.10*algae, 2)
    r.grade = (
        "Excellent" if r.score >= 80 else
        "Good"      if r.score >= 60 else
        "Fair"      if r.score >= 40 else
        "Poor"      if r.score >= 20 else "Critical"
    )
    return r


def compute_wqi(w: WaterQualityIndex) -> WaterQualityIndex:
    """
    Water Quality Index (0–100).  Weighted sub-indices inspired by NSF-WQI,
    tuned for open-ocean marine conditions.
    DO 17 %, pH 11 %, turbidity 8 %, temperature 10 %,
    nitrate 10 %, phosphate 10 %, salinity 34 %.
    """
    do_sat  = dissolved_oxygen_saturation_mgl(w.temperature_c, w.salinity_ppt)
    do_idx  = min(w.dissolved_oxygen_mgl / do_sat * 100.0, 100.0)
    ph_idx  = max(0.0, 100.0 - abs(w.ph - 8.1) * 50.0)
    trb_idx = max(0.0, 100.0 - w.turbidity_ntu * 5.0)
    tmp_idx = max(0.0, 100.0 - abs(w.temperature_c - 26.0) * 4.0)
    nit_idx = max(0.0, 100.0 - (w.nitrate_mgl   / 0.5) * 100.0)
    pho_idx = max(0.0, 100.0 - (w.phosphate_mgl / 0.1) * 100.0)
    sal_dev = max(0.0, 33 - w.salinity_ppt) + max(0.0, w.salinity_ppt - 37)
    sal_idx = max(0.0, 100.0 - sal_dev * 5.0)
    w.wqi = round(
        0.17*do_idx + 0.11*ph_idx + 0.08*trb_idx + 0.10*tmp_idx +
        0.10*nit_idx + 0.10*pho_idx + 0.34*sal_idx, 2)
    w.category = (
        "Excellent" if w.wqi >= 90 else
        "Good"      if w.wqi >= 70 else
        "Medium"    if w.wqi >= 50 else
        "Bad"       if w.wqi >= 25 else "Very Bad"
    )
    return w


def logistic_growth(N0: float, r: float, K: float,
                    t_max: int, dt: float = 1.0) -> list:
    """
    Logistic growth  dN/dt = r·N·(1 - N/K)  integrated with RK4.
    Returns list of (t, N) tuples.
    """
    results = [(0.0, N0)]
    N, t = N0, 0.0
    for _ in range(t_max):
        f  = lambda n: r * n * (1.0 - n / K)
        k1 = f(N)
        k2 = f(N + dt*k1/2)
        k3 = f(N + dt*k2/2)
        k4 = f(N + dt*k3)
        N  = max(0.0, N + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4))
        t += dt
        results.append((round(t, 2), round(N, 4)))
    return results


def lotka_volterra(N0: float, P0: float, alpha: float, beta: float,
                   delta: float, gamma: float,
                   t_max: int, dt: float = 0.1) -> list:
    """
    Lotka-Volterra predator-prey (RK4):
      dN/dt = α·N − β·N·P   (prey growth − predation)
      dP/dt = δ·N·P − γ·P   (predator growth − mortality)
    Analytical equilibrium: N* = γ/δ,  P* = α/β.
    Returns list of (t, N, P) tuples.
    """
    results = [(0.0, N0, P0)]
    N, P, t = N0, P0, 0.0
    steps = int(round(t_max / dt))
    for _ in range(steps):
        dN = lambda n, p: alpha*n - beta*n*p
        dP = lambda n, p: delta*n*p - gamma*p
        k1n, k1p = dN(N, P), dP(N, P)
        k2n, k2p = dN(N+dt*k1n/2, P+dt*k1p/2), dP(N+dt*k1n/2, P+dt*k1p/2)
        k3n, k3p = dN(N+dt*k2n/2, P+dt*k2p/2), dP(N+dt*k2n/2, P+dt*k2p/2)
        k4n, k4p = dN(N+dt*k3n,   P+dt*k3p),   dP(N+dt*k3n,   P+dt*k3p)
        N = max(0.0, N + (dt/6.0)*(k1n + 2*k2n + 2*k3n + k4n))
        P = max(0.0, P + (dt/6.0)*(k1p + 2*k2p + 2*k3p + k4p))
        t = round(t + dt, 6)
        results.append((round(t, 2), round(N, 4), round(P, 4)))
    return results


def trophic_transfer(producer_biomass: float, level: int,
                     efficiency: float = 0.10) -> float:
    """Lindeman's 10 % trophic efficiency rule."""
    return producer_biomass * (efficiency ** (level - 1))


# ---------------------------------------------------------------------------
# ASCII Visualisations
# ---------------------------------------------------------------------------

def ascii_population_chart(data: list, width: int = 58, height: int = 14,
                            label: str = "N") -> str:
    """Render a time-series population chart in ASCII art."""
    if not data:
        return ""
    has_pred = len(data[0]) == 3
    vn = [row[1] for row in data]
    vp = [row[2] for row in data] if has_pred else []
    all_v = vn + vp
    y_min, y_max = min(all_v), max(all_v)
    if y_max == y_min:
        y_max = y_min + 1
    step = max(1, len(data) // width)
    sn = vn[::step][:width]
    sp = vp[::step][:width] if has_pred else []
    cols = len(sn)
    grid = [[' '] * cols for _ in range(height)]

    def row_of(v):
        return max(0, min(height-1, int((1.0 - (v - y_min)/(y_max - y_min)) * (height-1))))

    for c, v in enumerate(sn):
        grid[row_of(v)][c] = f"{OCEAN_BLUE}●{RESET}"
    for c, v in enumerate(sp):
        r = row_of(v)
        if grid[r][c] == ' ':
            grid[r][c] = f"{CORAL}▲{RESET}"

    lines = [f"{BOLD}{OCEAN_DEEP}  {'─'*(cols+10)}{RESET}"]
    for i, row in enumerate(grid):
        yv = y_max - (y_max - y_min) * i / (height - 1)
        lines.append(f"{SAND}{yv:>9.1f} │{RESET}" + ''.join(row))
    lines.append(f"          └{'─'*cols}")
    lines.append(f"           t=0{' '*(cols-8)}t={data[-1][0]}")
    legend = f"  {OCEAN_BLUE}● {label}{RESET}"
    if has_pred:
        legend += f"   {CORAL}▲ Predator{RESET}"
    lines.append(legend)
    return '\n'.join(lines)


def ascii_gauge(score: float, label: str, grade: str) -> str:
    """Horizontal ASCII gauge for a 0–100 score."""
    w = 48
    filled = int(score / 100.0 * w)
    colour = (SEAFOAM if score >= 80 else
              OCEAN_BLUE if score >= 60 else
              SAND if score >= 40 else
              CORAL if score >= 20 else '\033[38;5;196m')
    bar = f"{colour}{'█'*filled}{RESET}{'░'*(w-filled)}"
    return (f"\n  {BOLD}{label}{RESET}\n"
            f"  [{bar}] {BOLD}{colour}{score:.1f}/100{RESET}\n"
            f"  Grade: {BOLD}{colour}{grade}{RESET}\n")


# ---------------------------------------------------------------------------
# CLI command handlers
# ---------------------------------------------------------------------------

def cmd_observe(args, conn):
    obs = SpeciesObservation(
        species_name=args.species, common_name=args.common_name or args.species,
        count=args.count, location_lat=args.lat, location_lon=args.lon,
        depth_m=args.depth, water_temp_c=args.temp,
        trophic_level=args.trophic, observer=args.observer or "field_team",
        notes=args.notes or "",
    )
    conn.execute(
        "INSERT INTO observations VALUES (NULL,?,?,?,?,?,?,?,?,?,?,?)",
        (obs.species_name, obs.common_name, obs.count, obs.location_lat,
         obs.location_lon, obs.depth_m, obs.water_temp_c, obs.trophic_level,
         obs.observer, obs.notes, obs.timestamp))
    conn.commit()
    biomass = trophic_transfer(10_000, int(obs.trophic_level))
    print(f"\n{SEAFOAM}{BOLD}✓ Observation recorded{RESET}")
    print(f"  {OCEAN_BLUE}Species:{RESET}  {obs.species_name}  (count={obs.count})")
    print(f"  {OCEAN_BLUE}Location:{RESET} {obs.location_lat}°N, {obs.location_lon}°E  depth={obs.depth_m}m")
    print(f"  {SAND}Available biomass at trophic level {int(obs.trophic_level)}: "
          f"{biomass:.1f} kg/km² (10 % rule){RESET}")


def cmd_model_population(args, conn):
    print(f"\n{BOLD}{OCEAN_DEEP}━━━ Population Dynamics Model ━━━{RESET}\n")
    if args.predator_pop is not None:
        data = lotka_volterra(args.population, args.predator_pop,
                              args.alpha, args.beta, args.delta, args.gamma,
                              args.time_steps)
        N_eq, P_eq = args.gamma / args.delta, args.alpha / args.beta
        print(f"  Model:       {CORAL}Lotka-Volterra{RESET}  "
              f"(α={args.alpha} β={args.beta} δ={args.delta} γ={args.gamma})")
        print(f"  Equilibrium: N*={N_eq:.2f}  P*={P_eq:.2f}")
        t, N, P = data[-1]
        print(f"  t={t}: Prey={N:.2f}  Predator={P:.2f}\n")
        print(ascii_population_chart(data, label="Prey"))
    else:
        data = logistic_growth(args.population, args.growth_rate,
                               args.carrying_capacity, args.time_steps)
        t, N = data[-1]
        pct = N / args.carrying_capacity * 100
        print(f"  Model:       {SEAFOAM}Logistic Growth{RESET}  "
              f"(r={args.growth_rate}  K={args.carrying_capacity})")
        print(f"  t={t}: N={N:.2f}  ({pct:.1f}% of K)\n")
        print(ascii_population_chart([(t, n) for t, n in data],
                                     label=args.species or "N"))


def cmd_reef_health(args, conn):
    r = ReefHealthScore(
        reef_id=args.reef_id, coral_cover_pct=args.coral,
        bleaching_pct=args.bleaching, fish_density_per_100m2=args.fish_density,
        water_clarity_m=args.clarity, macroalgae_cover_pct=args.macroalgae,
        invertebrate_density=args.invertebrates)
    r = score_reef_health(r)
    conn.execute(
        "INSERT INTO reef_health VALUES (NULL,?,?,?,?,?,?,?,?,?,?)",
        (r.reef_id, r.coral_cover_pct, r.bleaching_pct, r.fish_density_per_100m2,
         r.water_clarity_m, r.macroalgae_cover_pct, r.invertebrate_density,
         r.score, r.grade, r.timestamp))
    conn.commit()
    print(f"\n{BOLD}{OCEAN_DEEP}━━━ Reef Health: {r.reef_id} ━━━{RESET}")
    print(ascii_gauge(r.score, "Reef Health Score", r.grade))
    print(f"  {CORAL}Coral cover:{RESET}   {r.coral_cover_pct}%")
    print(f"  {CORAL}Bleaching:{RESET}     {r.bleaching_pct}%")
    print(f"  {OCEAN_BLUE}Fish density:{RESET}  {r.fish_density_per_100m2}/100m²")
    print(f"  {OCEAN_BLUE}Secchi depth:{RESET}  {r.water_clarity_m} m")
    print(f"  {KELP}Macroalgae:{RESET}    {r.macroalgae_cover_pct}%")


def cmd_water_quality(args, conn):
    w = WaterQualityIndex(
        station_id=args.station, dissolved_oxygen_mgl=args.do_mgl,
        ph=args.ph, turbidity_ntu=args.turbidity, temperature_c=args.temp,
        nitrate_mgl=args.nitrate, phosphate_mgl=args.phosphate,
        salinity_ppt=args.salinity)
    w = compute_wqi(w)
    do_sat = dissolved_oxygen_saturation_mgl(w.temperature_c, w.salinity_ppt)
    conn.execute(
        "INSERT INTO water_quality VALUES (NULL,?,?,?,?,?,?,?,?,?,?,?)",
        (w.station_id, w.dissolved_oxygen_mgl, w.ph, w.turbidity_ntu,
         w.temperature_c, w.nitrate_mgl, w.phosphate_mgl,
         w.salinity_ppt, w.wqi, w.category, w.timestamp))
    conn.commit()
    print(f"\n{BOLD}{OCEAN_DEEP}━━━ Water Quality: {w.station_id} ━━━{RESET}")
    print(ascii_gauge(w.wqi, "Water Quality Index", w.category))
    print(f"  {BIOLUM}DO:{RESET}          {w.dissolved_oxygen_mgl} mg/L "
          f"({w.dissolved_oxygen_mgl/do_sat*100:.1f}% sat.)")
    print(f"  {BIOLUM}DO sat (G&G):{RESET} {do_sat:.3f} mg/L @ {w.temperature_c}°C, {w.salinity_ppt} ppt")
    print(f"  {BIOLUM}pH:{RESET}          {w.ph}    {BIOLUM}Turbidity:{RESET} {w.turbidity_ntu} NTU")
    print(f"  {BIOLUM}Nitrate:{RESET}     {w.nitrate_mgl} mg/L    "
          f"{BIOLUM}Phosphate:{RESET} {w.phosphate_mgl} mg/L")


def cmd_biodiversity(args, conn):
    counts = {}
    for item in args.species_counts:
        sp, n = item.rsplit(":", 1)
        counts[sp.strip()] = int(n)
    bd = compute_biodiversity(args.site_id, counts)
    total = sum(counts.values())
    print(f"\n{BOLD}{OCEAN_DEEP}━━━ Biodiversity: {bd.site_id} ━━━{RESET}")
    print(f"  {SEAFOAM}Richness (obs):{RESET}  {bd.species_richness}   "
          f"{SEAFOAM}Chao1:{RESET} {bd.chao1}   {SEAFOAM}ACE:{RESET} {bd.ace}")
    print(f"  {OCEAN_BLUE}Shannon H':{RESET}     {bd.shannon_h:.4f} nats")
    print(f"  {OCEAN_BLUE}Simpson D:{RESET}      {bd.simpson_d:.4f}")
    print(f"  {OCEAN_BLUE}Pielou J:{RESET}       {bd.evenness:.4f}")
    print(f"\n  {SAND}Species (n={total}):{RESET}")
    for sp, n in sorted(counts.items(), key=lambda x: -x[1]):
        pi = n / total
        bar = f"{OCEAN_BLUE}{'█' * int(pi*30)}{RESET}"
        print(f"    {BIOLUM}{sp:<30}{RESET}  {n:>5}  ({pi*100:5.1f}%)  {bar}")


def cmd_migration(args, conn):
    dist  = haversine_km(args.orig_lat, args.orig_lon, args.dest_lat, args.dest_lon)
    speed = dist / args.days if args.days > 0 else 0.0
    rec   = MigrationRecord(
        species=args.species, tag_id=args.tag_id,
        origin_lat=args.orig_lat, origin_lon=args.orig_lon,
        destination_lat=args.dest_lat, destination_lon=args.dest_lon,
        duration_days=args.days, distance_km=round(dist, 2),
        avg_speed_kmday=round(speed, 2))
    conn.execute(
        "INSERT INTO migrations VALUES (NULL,?,?,?,?,?,?,?,?,?,?)",
        (rec.species, rec.tag_id, rec.origin_lat, rec.origin_lon,
         rec.destination_lat, rec.destination_lon, rec.duration_days,
         rec.distance_km, rec.avg_speed_kmday, rec.timestamp))
    conn.commit()
    print(f"\n{BOLD}{OCEAN_DEEP}━━━ Migration: {rec.species} ({rec.tag_id}) ━━━{RESET}")
    print(f"  {CORAL}Origin:{RESET}      {rec.origin_lat}°N, {rec.origin_lon}°E")
    print(f"  {CORAL}Destination:{RESET} {rec.destination_lat}°N, {rec.destination_lon}°E")
    print(f"  {SEAFOAM}Distance:{RESET}    {rec.distance_km} km (Haversine)")
    print(f"  {SEAFOAM}Duration:{RESET}    {rec.duration_days} days")
    print(f"  {SEAFOAM}Avg speed:{RESET}   {rec.avg_speed_kmday} km/day")


def cmd_report(args, conn):
    print(f"\n{BOLD}{OCEAN_DEEP}{'═'*52}")
    print(f"   BlackRoad Marine Biology Lab — Report")
    print(f"{'═'*52}{RESET}\n")
    r = conn.execute("SELECT COUNT(*), COUNT(DISTINCT species_name) FROM observations").fetchone()
    print(f"  {SEAFOAM}Observations:{RESET}  {r[0]} records, {r[1]} species")
    rows = conn.execute(
        "SELECT reef_id, score, grade FROM reef_health ORDER BY timestamp DESC LIMIT 5").fetchall()
    if rows:
        print(f"\n  {CORAL}Recent Reef Assessments:{RESET}")
        for rid, sc, gr in rows:
            bar = '█' * int(sc/5)
            print(f"    {rid:<22} {sc:5.1f}/100  [{bar:<20}] {gr}")
    rows = conn.execute(
        "SELECT station_id, wqi, category FROM water_quality ORDER BY timestamp DESC LIMIT 5").fetchall()
    if rows:
        print(f"\n  {BIOLUM}Recent Water Quality:{RESET}")
        for sid, wqi, cat in rows:
            print(f"    {sid:<22} WQI={wqi:5.1f}  {cat}")
    rows = conn.execute(
        "SELECT species, distance_km, avg_speed_kmday FROM migrations ORDER BY timestamp DESC LIMIT 5").fetchall()
    if rows:
        print(f"\n  {OCEAN_BLUE}Recent Migrations:{RESET}")
        for sp, dist, spd in rows:
            print(f"    {sp:<28} {dist:8.1f} km   {spd:.1f} km/day")
    print(f"\n{SAND}Report generated: {datetime.utcnow().isoformat()}Z{RESET}\n")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="marine",
        description=f"{BOLD}BlackRoad Marine Biology Lab{RESET} — Scientific Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    # observe
    p = sub.add_parser("observe", help="Record a species observation")
    p.add_argument("--species", required=True)
    p.add_argument("--common-name")
    p.add_argument("--count", type=int, default=1)
    p.add_argument("--lat", type=float, required=True)
    p.add_argument("--lon", type=float, required=True)
    p.add_argument("--depth", type=float, default=0.0)
    p.add_argument("--temp", type=float, default=25.0)
    p.add_argument("--trophic", type=float, default=2.0)
    p.add_argument("--observer")
    p.add_argument("--notes")

    # model-population
    p = sub.add_parser("model-population", help="Run logistic or Lotka-Volterra model")
    p.add_argument("--species")
    p.add_argument("--population", type=float, required=True)
    p.add_argument("--carrying-capacity", type=float, default=1000.0)
    p.add_argument("--growth-rate", type=float, default=0.3)
    p.add_argument("--time-steps", type=int, default=100)
    p.add_argument("--predator-pop", type=float, default=None)
    p.add_argument("--alpha", type=float, default=1.0,   help="Prey growth rate")
    p.add_argument("--beta",  type=float, default=0.1,   help="Predation rate")
    p.add_argument("--delta", type=float, default=0.075, help="Predator efficiency")
    p.add_argument("--gamma", type=float, default=1.5,   help="Predator mortality")

    # reef-health
    p = sub.add_parser("reef-health", help="Compute reef health composite score")
    p.add_argument("--reef-id", required=True)
    p.add_argument("--coral",        type=float, required=True)
    p.add_argument("--bleaching",    type=float, required=True)
    p.add_argument("--fish-density", type=float, required=True)
    p.add_argument("--clarity",      type=float, required=True, help="Secchi depth (m)")
    p.add_argument("--macroalgae",   type=float, required=True)
    p.add_argument("--invertebrates",type=float, default=10.0)

    # water-quality
    p = sub.add_parser("water-quality", help="Compute Water Quality Index")
    p.add_argument("--station",   required=True)
    p.add_argument("--do-mgl",    type=float, required=True)
    p.add_argument("--ph",        type=float, required=True)
    p.add_argument("--turbidity", type=float, required=True)
    p.add_argument("--temp",      type=float, required=True)
    p.add_argument("--nitrate",   type=float, required=True)
    p.add_argument("--phosphate", type=float, required=True)
    p.add_argument("--salinity",  type=float, required=True)

    # biodiversity
    p = sub.add_parser("biodiversity", help="Compute Shannon, Simpson, Chao1, ACE")
    p.add_argument("--site-id", required=True)
    p.add_argument("--species-counts", nargs="+", required=True, metavar="SPECIES:COUNT")

    # migration
    p = sub.add_parser("migration", help="Record migration track (Haversine distance)")
    p.add_argument("--species",  required=True)
    p.add_argument("--tag-id",   required=True)
    p.add_argument("--orig-lat", type=float, required=True)
    p.add_argument("--orig-lon", type=float, required=True)
    p.add_argument("--dest-lat", type=float, required=True)
    p.add_argument("--dest-lon", type=float, required=True)
    p.add_argument("--days",     type=float, required=True)

    # report
    sub.add_parser("report", help="Database summary report")
    return parser


def main():
    parser = build_parser()
    args   = parser.parse_args()
    conn   = init_db()
    {
        "observe":          cmd_observe,
        "model-population": cmd_model_population,
        "reef-health":      cmd_reef_health,
        "water-quality":    cmd_water_quality,
        "biodiversity":     cmd_biodiversity,
        "migration":        cmd_migration,
        "report":           cmd_report,
    }[args.command](args, conn)
    conn.close()


if __name__ == "__main__":
    main()
