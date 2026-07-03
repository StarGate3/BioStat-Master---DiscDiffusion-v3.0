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
# WYKRYWANIE KONTROLI (grupa referencyjna)
# ============================================================

NEGATIVE_CONTROL_SIGNALS: tuple = ('(-)', 'woda', 'water', 'negatyw', 'vehicle', 'dmso')
"""Case-insensitive substrings identifying a NEGATIVE control group
(untreated/vehicle baseline). A group must match one of these AND none of
POSITIVE_CONTROL_SIGNALS / KNOWN_ANTIBIOTIC_SUBSTRINGS to auto-qualify as
the reference group for statistical comparisons."""

POSITIVE_CONTROL_SIGNALS: tuple = ('(+)', 'pozytyw', 'positive')
"""Case-insensitive substrings identifying a POSITIVE control group
(reference antibiotic). Groups matching these are NEVER auto-selected as
the reference group, even if they also contain a "kontrol"-like substring."""

KNOWN_ANTIBIOTIC_SUBSTRINGS: tuple = (
    'ampicylin', 'ampicillin', 'ampycylin', 'amoksycylin', 'amoxicillin',
    'penicylin', 'penicillin', 'streptomycyn', 'streptomycin',
    'tetracyklin', 'tetracycline', 'erytromycyn', 'erythromycin',
    'gentamycyn', 'gentamicin', 'wankomycyn', 'vancomycin',
    'cyprofloksacyn', 'ciprofloxacin', 'chloramfenikol', 'chloramphenicol',
    'kanamycyn', 'kanamycin',
)
"""Case-insensitive substrings of common reference-antibiotic names, incl.
the 'Ampycylina' typo seen in real lab spreadsheets. Treated the same as
POSITIVE_CONTROL_SIGNALS: never auto-selected as the reference group. Not
exhaustive - extend as needed for antibiotics used in your assays."""

REF_PLACEHOLDER: str = "-- Wybierz ręcznie (niejednoznaczne) --"
"""Sentinel shown in the reference-group dropdown when auto-detection could
not unambiguously identify a negative control. Never a real group name;
run_analysis refuses to proceed while this value is selected."""

# ============================================================
# MIC ESTIMATION - QUALITY GATE
# ============================================================

MIC_MIN_R2: float = 0.80
"""Minimum R^2 (log-linear fit quality) required to report a numeric MIC
estimate as reliable. Below this threshold the MIC is suppressed (not
shown as a number) and flagged as unreliable instead. This is a pragmatic
screening-tool bar, not a validated bioanalytical-method-validation
threshold: strict enough to reject poor fits, but not so strict that a
3-point biological dose-response curve (the practical minimum this
regression accepts) gets rejected outright. Extrapolation beyond the
tested concentration range is flagged separately (see estimate_mic),
independent of this threshold."""

# ============================================================
# NOWY FORMAT WEJŚCIOWY (arkusz "Dane" z powtórzeniami bio/tech)
# ============================================================

NEW_FORMAT_SHEET_NAME: str = "Dane"
"""Nazwa arkusza nowego szablonu wejściowego. Jeśli obecna w pliku, ten
arkusz jest czytany zamiast domyślnego pierwszego arkusza; w przeciwnym
razie loader cofa się do dotychczasowego zachowania (pierwszy/domyślny
arkusz) - wsteczna zgodność ze starymi plikami bez tego arkusza."""

COL_SUBSTANCE: str = 'Substancja'
COL_CONCENTRATION: str = 'Stezenie'
COL_UNIT: str = 'Jednostka'
COL_TYPE: str = 'Typ'
COL_REP_BIO: str = 'Rep_biologiczna'
COL_REP_TECH: str = 'Rep_techniczna'
"""Nagłówki kolumn nowego formatu wejściowego (arkusz 'Dane'). Każda z nich
jest wykrywana i obsługiwana NIEZALEŻNIE - plik może mieć część, ale nie
wszystkie z nich (np. Typ bez kolumn powtórzeń); braki są uzupełniane
rozsądnymi wartościami domyślnymi, patrz utils.build_internal_representation."""

TYPE_TESTED: str = 'Badana'
TYPE_NEG_CONTROL: str = 'Kontrola negatywna'
TYPE_POS_CONTROL: str = 'Kontrola pozytywna'
VALID_TYPES: tuple = (TYPE_TESTED, TYPE_NEG_CONTROL, TYPE_POS_CONTROL)
"""Dozwolone wartości kolumny Typ w nowym formacie."""

INTERNAL_TYPE_COL: str = '_Typ'
INTERNAL_SUBSTANCE_COL: str = '_Substancja'
INTERNAL_CONC_COL: str = '_Stezenie'
INTERNAL_UNIT_COL: str = '_Jednostka'
INTERNAL_REP_BIO_COL: str = '_Rep_biologiczna'
INTERNAL_REP_TECH_COL: str = '_Rep_techniczna'
"""Kolumny WEWNĘTRZNE, nie pochodzą bezpośrednio z pliku wejściowego -
zawsze obecne po przejściu przez utils.build_internal_representation,
niezależnie od formatu źródłowego. Wypełnione NaN (Typ/Substancja/
Stezenie/Jednostka) albo wartościami domyślnymi sygnalizującymi "cały
wiersz danej grupy to jedno powtórzenie biologiczne" (Rep_biologiczna=1,
Rep_techniczna=kolejny numer w obrębie grupy) gdy dany aspekt nie
występuje w pliku wejściowym (stary format)."""

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
