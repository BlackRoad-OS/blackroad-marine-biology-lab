# 🌊 BlackRoad Marine Biology Lab

A scientific marine biology analysis toolkit implementing peer-reviewed ecological algorithms for field research, reef monitoring, and population dynamics modelling.  All calculations use only the Python standard library.

---

## Scientific Background

### Biodiversity Indices

**Shannon-Wiener Index H′**

Quantifies information entropy of a community:

$$H' = -\sum_{i=1}^{S} p_i \ln(p_i)$$

where $p_i = n_i / N$ is the relative abundance of species $i$.  $H' = 0$ for a monoculture; $H' = \ln S$ when all species are equally abundant.

**Simpson's Diversity Index D**

Probability that two randomly drawn individuals belong to *different* species:

$$D = 1 - \frac{\sum_{i=1}^{S} n_i(n_i - 1)}{N(N-1)}$$

Range $[0, 1]$; 1 = maximum diversity.

**Pielou's Evenness J**

Normalised Shannon index removing the species-richness effect:

$$J = \frac{H'}{\ln S}$$

**Species Richness Estimators**

*Chao1* (Chao, 1984) — non-parametric estimator using singletons $f_1$ and doubletons $f_2$:

$$\hat{S}_{Chao1} = S_{obs} + \frac{f_1^2}{2 f_2}$$

*ACE* (Chao & Lee, 1992) — abundance-based coverage estimator:

$$\hat{S}_{ACE} = S_{abund} + \frac{S_{rare}}{\hat{C}_{ACE}} + \frac{f_1}{\hat{C}_{ACE}} \hat{\gamma}^2_{ACE}, \quad \hat{C}_{ACE} = 1 - \frac{f_1}{n_{rare}}$$

---

### Population Dynamics

**Logistic Growth**

$$\frac{dN}{dt} = rN\!\left(1 - \frac{N}{K}\right)$$

Integrated numerically with 4th-order Runge-Kutta.  The population grows exponentially when $N \ll K$ and plateaus at carrying capacity $K$.

**Lotka-Volterra Predator-Prey**

$$\frac{dN}{dt} = \alpha N - \beta NP$$
$$\frac{dP}{dt} = \delta NP - \gamma P$$

| Symbol | Meaning |
|--------|---------|
| $\alpha$ | Prey intrinsic growth rate |
| $\beta$ | Predation rate coefficient |
| $\delta$ | Predator reproduction efficiency |
| $\gamma$ | Predator natural mortality rate |

Analytical equilibrium: $N^* = \gamma / \delta$, $P^* = \alpha / \beta$.  The system produces neutrally-stable oscillations around this fixed point.

**Trophic Transfer (Lindeman, 1942)**

$$\text{Biomass at level } L = B_0 \cdot e^L_{eff}^{(L-1)}$$

where $e^L_{eff} \approx 10\%$ per trophic step.

---

### Water Quality

**Dissolved Oxygen Saturation — Garcia & Gordon (1992)**

$$\ln\!\left[\frac{DO_{sat}}{\text{ml L}^{-1}}\right] = \sum_{k=0}^{5} A_k T_s^k + S\!\sum_{k=0}^{3} B_k T_s^k + C_0 S^2$$

where $T_s = \ln\!\bigl(\frac{298.15 - T}{273.15 + T}\bigr)$ and $S$ is salinity in g kg⁻¹.  Multiply by 1.4276 to convert ml/L → mg/L.

**Water Quality Index**

Weighted composite (0–100) of seven parameters tuned for open-ocean conditions:

| Parameter | Weight | Marine reference |
|-----------|--------|-----------------|
| DO saturation | 17 % | ≥ 90 % |
| pH | 11 % | 8.1 ± 0.2 |
| Turbidity | 8 % | < 1 NTU |
| Temperature | 10 % | 26 °C optimum |
| Nitrate | 10 % | < 0.07 mg/L |
| Phosphate | 10 % | < 0.01 mg/L |
| Salinity | 34 % | 33–37 ppt |

---

### Reef Health

Composite score (0–100) from weighted sub-indices:

| Parameter | Weight | Ideal |
|-----------|--------|-------|
| Coral cover | 30 % | ≥ 50 % |
| Bleaching (penalty) | 25 % | 0 % |
| Fish density | 20 % | ≥ 50 per 100 m² |
| Water clarity (Secchi) | 15 % | ≥ 20 m |
| Macroalgae (penalty) | 10 % | 0 % |

Grades: ≥ 80 Excellent · ≥ 60 Good · ≥ 40 Fair · ≥ 20 Poor · < 20 Critical.

---

### Migration

**Haversine Great-Circle Distance**

$$d = 2R \arcsin\!\left(\sqrt{\sin^2\!\tfrac{\Delta\phi}{2} + \cos\phi_1\cos\phi_2\sin^2\!\tfrac{\Delta\lambda}{2}}\right)$$

$R = 6371$ km.  Used to compute migration distance and average daily swim speed from GPS tag records.

---

## Installation

```bash
# No external dependencies — uses Python stdlib only
python --version   # 3.9+
```

---

## Usage

```bash
# Record a species observation
python src/marine.py observe \
  --species "Acropora millepora" --common-name "Staghorn Coral" \
  --count 47 --lat -18.28 --lon 147.70 --depth 6.0 --temp 27.1 --trophic 1.0

# Assess reef health (Great Barrier Reef station)
python src/marine.py reef-health \
  --reef-id GBR-S01 --coral 52 --bleaching 8 \
  --fish-density 38 --clarity 14 --macroalgae 4 --invertebrates 22

# Compute Water Quality Index
python src/marine.py water-quality \
  --station AIMS-01 --do-mgl 7.2 --ph 8.15 --turbidity 0.8 \
  --temp 26.5 --nitrate 0.03 --phosphate 0.008 --salinity 35.2

# Biodiversity analysis
python src/marine.py biodiversity --site-id Palau-North \
  --species-counts "Parrotfish:25" "Surgeonfish:18" "Angelfish:7" \
                   "Damselfish:42" "Moorish Idol:3" "Butterflyfish:11"

# Lotka-Volterra predator-prey model
python src/marine.py model-population \
  --species "Sardine" --population 10000 --predator-pop 200 \
  --alpha 1.0 --beta 0.005 --delta 0.002 --gamma 0.8 --time-steps 100

# Logistic growth
python src/marine.py model-population \
  --species "Bottlenose Dolphin" --population 120 \
  --carrying-capacity 500 --growth-rate 0.15 --time-steps 60

# Migration track
python src/marine.py migration \
  --species "Humpback Whale" --tag-id HW-2024-007 \
  --orig-lat -23.5 --orig-lon 152.0 \
  --dest-lat 18.0  --dest-lon 155.0 --days 47

# Full database report
python src/marine.py report
```

---

## Running Tests

```bash
# Run all tests
pytest tests/ -v --tb=short

# With coverage
pytest tests/ --cov=src --cov-report=term-missing
```

---

## Data Persistence

All observations, reef assessments, water-quality readings, and migration records are persisted to `marine_lab.db` (SQLite) in the working directory.  The database is created automatically on first run.

Schema tables: `observations`, `reef_health`, `water_quality`, `migrations`.

---

## References

- Garcia, H. E. & Gordon, L. I. (1992). Oxygen solubility in seawater: Better fitting equations. *Limnol. Oceanogr.* **37**(6), 1307–1312.
- Shannon, C. E. & Weaver, W. (1949). *The Mathematical Theory of Communication*. Univ. Illinois Press.
- Lotka, A. J. (1925). *Elements of Physical Biology*. Williams & Wilkins.
- Volterra, V. (1926). Fluctuations in the abundance of a species considered mathematically. *Nature* **118**, 558–560.
- Chao, A. (1984). Nonparametric estimation of the number of classes in a population. *Scand. J. Statist.* **11**, 265–270.
- Chao, A. & Lee, S.-M. (1992). Estimating the number of classes via sample coverage. *J. Am. Stat. Assoc.* **87**(417), 210–217.
- Pielou, E. C. (1966). The measurement of diversity in different types of biological collections. *J. Theor. Biol.* **13**, 131–144.
- Lindeman, R. L. (1942). The trophic-dynamic aspect of ecology. *Ecology* **23**(4), 399–417.
