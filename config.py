"""
Centralized configuration for BioStat Master.

Values here used to be scattered across logic.py, plotting.py, utils.py,
gui.py, dialogs.py, and reports.py. Changing a constant here should change
behavior project-wide; user-facing captions interpolate these values so
printed documentation stays in sync with actual computation.
"""
import re

# ============================================================
# EXPERIMENTAL CONSTANTS
# ============================================================

DISC_DIAMETER_MM: float = 6.0
"""Paper-disc diameter used in disc-diffusion assays. This is also the
target zone diameter the log-linear MIC extrapolation solves for."""

ALPHA: float = 0.05
"""Statistical significance threshold. Applied to Shapiro-Wilk normality
test, Levene variance homogeneity test, the ANOVA/Kruskal-Wallis main
test, and every pairwise post-hoc comparison."""

# ============================================================
# EXCEL SCHEMA
# ============================================================

COL_GROUP: str = 'Grupa'
"""Exact header name of the treatment-group column."""

COL_MEASUREMENT: str = 'Srednica_mm'
"""Exact header name of the zone-diameter measurement column."""

COL_BACT_SUBSTRING: str = 'Bakteri'
"""Case-sensitive substring used to auto-detect the bacteria column.
The first column whose name contains this substring is used."""

# ============================================================
# COHEN'S D EFFECT-SIZE THRESHOLDS (Cohen, 1988)
# ============================================================

COHENS_D_SMALL: float = 0.2
COHENS_D_MEDIUM: float = 0.5
COHENS_D_LARGE: float = 0.8
"""Boundaries for the qualitative interpretation bucket:
  |d| < 0.2     -> negligible ("znikomy")
  0.2 <= |d| < 0.5 -> small        ("mały")
  0.5 <= |d| < 0.8 -> medium       ("średni")
  |d| >= 0.8    -> large           ("DUŻY")
"""

# ============================================================
# FIGURE DPI
# ============================================================

SCREEN_DPI: int = 100
"""DPI used when building figures for on-screen display in the Tk canvas."""

PDF_DPI: int = 150
"""DPI used when rasterizing figures for embedding in the PDF report.
Balances file size against legibility at A4 print size."""

EXPORT_DPI: int = 300
"""DPI used when saving a standalone figure via the Zapisz Wykres button
(publication quality)."""

# ============================================================
# CONCENTRATION PARSING
# ============================================================

CONCENTRATION_UNITS: tuple = ('mg/ml', 'ug/ml', '%')
"""Unit strings recognized inside group labels such as 'Extract (50 mg/ml)'.
Used by utils.parse_concentration to extract the numeric concentration
and the unit for dose-response plots and MIC extrapolation."""

_UNITS_ALT = '|'.join(re.escape(u) for u in CONCENTRATION_UNITS)
CONCENTRATION_SEARCH_PATTERN: str = r"([\d,.]+)\s*(" + _UNITS_ALT + r")"
CONCENTRATION_STRIP_PATTERN: str = r"\s*\(?[\d,.]+\s*(" + _UNITS_ALT + r").*\)?"

# ============================================================
# OUTLIER DETECTION - Dixon Q-test critical values (alpha = 0.10)
# ============================================================

DIXON_Q90_CRITICAL: dict = {
    3: 0.941, 4: 0.765, 5: 0.642, 6: 0.560,
    7: 0.507, 8: 0.468, 9: 0.437, 10: 0.412,
}
"""Critical values for the Dixon Q test at the 90% confidence level,
for sample sizes n=3..10. Values above these indicate an outlier at
alpha = 0.10. Source: standard Dixon Q tables."""
