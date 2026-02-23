"""
Tests for BlackRoad Marine Biology Lab — scientific algorithm validation.
All assertions use known analytical results or published reference values.
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from marine import (
    shannon_wiener_index,
    simpsons_diversity_index,
    pielou_evenness,
    chao1_estimator,
    ace_estimator,
    compute_biodiversity,
    dissolved_oxygen_saturation_mgl,
    lotka_volterra,
    logistic_growth,
    score_reef_health,
    compute_wqi,
    haversine_km,
    ReefHealthScore,
    WaterQualityIndex,
)


# ---------------------------------------------------------------------------
# Shannon-Wiener Index
# ---------------------------------------------------------------------------

class TestShannonWienerIndex:
    """H' = -Σ pi·ln(pi)"""

    def test_two_equal_species(self):
        # 2 equally abundant species → H' = ln(2)
        h = shannon_wiener_index({"A": 50, "B": 50})
        assert abs(h - math.log(2)) < 1e-9

    def test_four_equal_species(self):
        # 4 equally abundant → H' = ln(4)
        h = shannon_wiener_index({"A": 25, "B": 25, "C": 25, "D": 25})
        assert abs(h - math.log(4)) < 1e-9

    def test_single_dominant_species(self):
        # One dominant species → H' ≈ 0
        h = shannon_wiener_index({"Dom": 10000, "Rare": 1})
        assert h < 0.01

    def test_known_distribution(self):
        # sp1=10, sp2=20, sp3=70, N=100 → pi = 0.1, 0.2, 0.7
        counts = {"sp1": 10, "sp2": 20, "sp3": 70}
        expected = -(0.1*math.log(0.1) + 0.2*math.log(0.2) + 0.7*math.log(0.7))
        assert abs(shannon_wiener_index(counts) - expected) < 1e-9

    def test_empty_community(self):
        assert shannon_wiener_index({}) == 0.0

    def test_zero_count_species_ignored(self):
        # Zero-count species must not affect index
        h = shannon_wiener_index({"sp1": 10, "sp2": 0, "sp3": 10})
        assert abs(h - math.log(2)) < 1e-9

    def test_h_prime_increases_with_richness(self):
        # Adding an equally-common species must raise H'
        h2 = shannon_wiener_index({"A": 10, "B": 10})
        h3 = shannon_wiener_index({"A": 10, "B": 10, "C": 10})
        assert h3 > h2


# ---------------------------------------------------------------------------
# Simpson's Index
# ---------------------------------------------------------------------------

class TestSimpsonsIndex:

    def test_four_equal_species(self):
        # D = 1 - 4·(25·24)/(100·99) = 1 - 2400/9900 ≈ 0.7576
        counts = {"A": 25, "B": 25, "C": 25, "D": 25}
        d = simpsons_diversity_index(counts)
        expected = 1.0 - 4 * 25 * 24 / (100 * 99)
        assert abs(d - expected) < 1e-9

    def test_single_species_zero(self):
        assert simpsons_diversity_index({"A": 100}) == 0.0

    def test_known_small_community(self):
        # n=(3,2,1), N=6 → D = 1 - (6+2+0)/30 = 1 - 8/30
        counts = {"A": 3, "B": 2, "C": 1}
        expected = 1.0 - (3*2 + 2*1 + 1*0) / (6*5)
        assert abs(simpsons_diversity_index(counts) - expected) < 1e-9

    def test_range(self):
        for counts in [{"A": 1}, {"A": 5, "B": 3}, {"A": 2, "B": 2, "C": 2}]:
            d = simpsons_diversity_index(counts)
            assert 0.0 <= d <= 1.0


# ---------------------------------------------------------------------------
# Lotka-Volterra
# ---------------------------------------------------------------------------

class TestLotkaVolterra:

    def test_equilibrium_stability(self):
        """Starting exactly at N*=γ/δ, P*=α/β should remain stationary."""
        alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5
        N_eq, P_eq = gamma / delta, alpha / beta   # 20.0, 10.0
        data = lotka_volterra(N_eq, P_eq, alpha, beta, delta, gamma,
                              t_max=50, dt=0.05)
        N_vals = [r[1] for r in data]
        P_vals = [r[2] for r in data]
        assert max(abs(n - N_eq) for n in N_vals) < 0.5
        assert max(abs(p - P_eq) for p in P_vals) < 0.5

    def test_oscillation_off_equilibrium(self):
        """Starting away from equilibrium produces oscillations."""
        data = lotka_volterra(40, 9, 1.0, 0.1, 0.075, 1.5,
                              t_max=100, dt=0.05)
        N_vals = [r[1] for r in data]
        assert max(N_vals) > min(N_vals) * 1.5

    def test_populations_non_negative(self):
        data = lotka_volterra(10, 5, 0.8, 0.05, 0.04, 0.6,
                              t_max=200, dt=0.05)
        assert all(r[1] >= 0 for r in data)
        assert all(r[2] >= 0 for r in data)

    def test_analytical_equilibrium_values(self):
        alpha, beta, delta, gamma = 2.0, 0.2, 0.1, 1.0
        N_eq, P_eq = gamma / delta, alpha / beta  # 10, 10
        data = lotka_volterra(N_eq, P_eq, alpha, beta, delta, gamma,
                              t_max=30, dt=0.05)
        # Equilibrium should be conserved within 2 %
        tol = N_eq * 0.02
        assert all(abs(r[1] - N_eq) < tol for r in data)


# ---------------------------------------------------------------------------
# Water Quality Index
# ---------------------------------------------------------------------------

class TestWaterQualityIndex:

    def _pristine(self):
        return WaterQualityIndex(
            station_id="pristine",
            dissolved_oxygen_mgl=7.2, ph=8.1, turbidity_ntu=0.2,
            temperature_c=25.0, nitrate_mgl=0.01, phosphate_mgl=0.005,
            salinity_ppt=35.0)

    def test_pristine_seawater_high_score(self):
        w = compute_wqi(self._pristine())
        assert w.wqi >= 85.0
        assert w.category in ("Excellent", "Good")

    def test_polluted_water_low_score(self):
        w = WaterQualityIndex(
            station_id="polluted",
            dissolved_oxygen_mgl=1.0, ph=5.0, turbidity_ntu=50.0,
            temperature_c=35.0, nitrate_mgl=5.0, phosphate_mgl=2.0,
            salinity_ppt=10.0)
        w = compute_wqi(w)
        assert w.wqi < 50.0

    def test_wqi_range(self):
        w = compute_wqi(self._pristine())
        assert 0.0 <= w.wqi <= 100.0

    def test_category_assigned(self):
        w = compute_wqi(self._pristine())
        assert w.category != ""


# ---------------------------------------------------------------------------
# Reef Health Scoring
# ---------------------------------------------------------------------------

class TestReefHealthScoring:

    def test_pristine_reef_excellent(self):
        r = score_reef_health(ReefHealthScore(
            reef_id="GBR-pristine",
            coral_cover_pct=70, bleaching_pct=0,
            fish_density_per_100m2=80, water_clarity_m=25,
            macroalgae_cover_pct=0, invertebrate_density=30))
        assert r.score >= 90.0
        assert r.grade == "Excellent"

    def test_bleached_degraded_reef_critical(self):
        r = score_reef_health(ReefHealthScore(
            reef_id="bleached",
            coral_cover_pct=3, bleaching_pct=95,
            fish_density_per_100m2=2, water_clarity_m=2,
            macroalgae_cover_pct=85, invertebrate_density=1))
        assert r.score < 20.0
        assert r.grade == "Critical"

    def test_score_bounded_0_100(self):
        r = score_reef_health(ReefHealthScore(
            "test", 50, 0, 50, 20, 0, 20))
        assert 0.0 <= r.score <= 100.0

    def test_bleaching_reduces_score(self):
        base = ReefHealthScore("B", 60, 0, 50, 15, 5, 20)
        high = ReefHealthScore("H", 60, 80, 50, 15, 5, 20)
        assert score_reef_health(base).score > score_reef_health(high).score


# ---------------------------------------------------------------------------
# Logistic Growth
# ---------------------------------------------------------------------------

class TestLogisticGrowth:

    def test_approaches_carrying_capacity(self):
        """N must converge to within 1 % of K after 50 steps."""
        K = 1000.0
        data = logistic_growth(N0=10, r=0.5, K=K, t_max=50)
        assert abs(data[-1][1] - K) / K < 0.01

    def test_growth_below_k(self):
        data = logistic_growth(N0=100, r=0.3, K=1000, t_max=20)
        assert data[-1][1] > data[0][1]

    def test_decline_above_k(self):
        K = 500.0
        data = logistic_growth(N0=900, r=0.3, K=K, t_max=30)
        assert data[-1][1] < data[0][1]
        assert data[-1][1] > K * 0.90

    def test_output_length(self):
        data = logistic_growth(100, 0.3, 1000, 40)
        assert len(data) == 41   # t = 0 … 40

    def test_population_non_negative(self):
        data = logistic_growth(N0=5, r=2.0, K=100, t_max=30)
        assert all(row[1] >= 0 for row in data)


# ---------------------------------------------------------------------------
# Chao1 Estimator
# ---------------------------------------------------------------------------

class TestChao1Estimator:

    def test_known_f1_f2(self):
        """Sobs=5, f1=2, f2=1  →  Chao1 = 5 + 4/2 = 7.0"""
        counts = {"A": 1, "B": 1, "C": 2, "D": 5, "E": 10}
        assert abs(chao1_estimator(counts) - 7.0) < 1e-9

    def test_no_singletons_equals_sobs(self):
        """With f1=0 and f2=0, estimator falls back to Sobs."""
        counts = {"A": 5, "B": 8, "C": 12, "D": 3}
        assert chao1_estimator(counts) == 4.0

    def test_chao1_ge_sobs(self):
        counts = {"A": 1, "B": 2, "C": 3, "D": 1, "E": 1}
        sobs = len(counts)
        assert chao1_estimator(counts) >= sobs

    def test_all_singletons_no_doubletons(self):
        """f1=4, f2=0  →  Sobs + f1*(f1-1)/2 = 4 + 6 = 10"""
        counts = {"A": 1, "B": 1, "C": 1, "D": 1}
        assert abs(chao1_estimator(counts) - 10.0) < 1e-9


# ---------------------------------------------------------------------------
# Dissolved Oxygen Saturation (Garcia & Gordon 1992)
# ---------------------------------------------------------------------------

class TestDissolvedOxygenSaturation:

    def test_standard_seawater_25c_35ppt(self):
        """At 25 °C and 35 ppt DO_sat ≈ 6.9–7.5 mg/L (published tables)."""
        do_sat = dissolved_oxygen_saturation_mgl(25.0, 35.0)
        assert 6.5 <= do_sat <= 7.5

    def test_cold_water_higher_do(self):
        do_warm = dissolved_oxygen_saturation_mgl(30.0, 35.0)
        do_cold = dissolved_oxygen_saturation_mgl(5.0, 35.0)
        assert do_cold > do_warm

    def test_higher_salinity_lower_do(self):
        """Salting-out effect: DO decreases with salinity."""
        do_fresh = dissolved_oxygen_saturation_mgl(20.0, 0.0)
        do_salt  = dissolved_oxygen_saturation_mgl(20.0, 35.0)
        assert do_fresh > do_salt

    def test_tropical_surface_range(self):
        do_sat = dissolved_oxygen_saturation_mgl(28.0, 34.0)
        assert 5.0 <= do_sat <= 8.5

    def test_do_positive(self):
        for t in [0, 10, 20, 30]:
            assert dissolved_oxygen_saturation_mgl(t, 35.0) > 0


# ---------------------------------------------------------------------------
# Haversine Distance
# ---------------------------------------------------------------------------

class TestHaversineDistance:

    def test_same_point_zero(self):
        assert haversine_km(-18.0, 147.0, -18.0, 147.0) == 0.0

    def test_known_distance_cairns_to_sydney(self):
        """Cairns (~16.9°S, 145.8°E) to Sydney (~33.9°S, 151.2°E) ≈ 1900 km."""
        dist = haversine_km(-16.9, 145.8, -33.9, 151.2)
        assert 1700 < dist < 2100

    def test_antipodal_half_circumference(self):
        """(0,0) to (0,180) ≈ half Earth circumference = 20015 km."""
        dist = haversine_km(0.0, 0.0, 0.0, 180.0)
        assert abs(dist - 20015.09) < 20

    def test_symmetry(self):
        d1 = haversine_km(-10, 130, -30, 160)
        d2 = haversine_km(-30, 160, -10, 130)
        assert abs(d1 - d2) < 1e-6


# ---------------------------------------------------------------------------
# Full Biodiversity Metrics
# ---------------------------------------------------------------------------

class TestBiodiversityMetrics:

    def test_compute_full_metrics(self):
        counts = {"Parrotfish": 25, "Surgeonfish": 18, "Angelfish": 7,
                  "Damselfish": 42, "Moorishidol": 3}
        bd = compute_biodiversity("TestSite", counts)
        assert bd.species_richness == 5
        assert bd.shannon_h > 0
        assert 0 < bd.simpson_d < 1
        assert bd.chao1 >= 5
        assert 0 <= bd.evenness <= 1

    def test_single_species_minimum_diversity(self):
        bd = compute_biodiversity("Mono", {"OnlySpecies": 100})
        assert bd.shannon_h == 0.0
        assert bd.simpson_d == 0.0
        assert bd.evenness == 0.0

    def test_more_species_higher_h_prime(self):
        bd2 = compute_biodiversity("S2", {"A": 50, "B": 50})
        bd4 = compute_biodiversity("S4", {"A": 25, "B": 25, "C": 25, "D": 25})
        assert bd4.shannon_h > bd2.shannon_h
