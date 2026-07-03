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
"""Paper-disc diameter used in disc-diffusion assays. Drawn as a
reference line on the main comparison plot."""

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
and the unit for the dose-response trend plot."""

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
# NOWY FORMAT WEJŚCIOWY (arkusz "Dane" z powtórzeniami bio/tech)
# ============================================================

NEW_FORMAT_SHEET_NAME: str = "Dane_dyfuzja"
"""Nazwa arkusza z danymi dyfuzji krążkowej w wieloarkuszowym szablonie
(patrz sekcja "ROUTER WIELOARKUSZOWY" niżej). Jeśli obecna w pliku, ten
arkusz jest czytany zamiast domyślnego pierwszego arkusza; w przeciwnym
razie loader cofa się do dotychczasowego zachowania (pierwszy/domyślny
arkusz) - wsteczna zgodność ze starymi, jednoarkuszowymi plikami."""

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
# ROUTER WIELOARKUSZOWY (MIC/MBC, Faza 1: wykrywanie i wybór - BEZ analizy)
# ============================================================

SHEET_DIFFUSION: str = NEW_FORMAT_SHEET_NAME  # "Dane_dyfuzja" - alias dla czytelności w kontekście routera
SHEET_MIC_VISUAL: str = "MIC_wizualny"
SHEET_MIC_OD: str = "MIC_OD"
SHEET_MBC: str = "MBC_posiew"
SHEET_CONTROLS: str = "Kontrole"
SHEET_INSTRUCTIONS: str = "Instrukcja"
SHEET_SETTINGS: str = "Ustawienia"
"""Nazwy wszystkich arkuszy docelowego, wieloarkuszowego szablonu wejściowego
(maks. 7 arkuszy). SHEET_INSTRUCTIONS i SHEET_SETTINGS są zawsze ignorowane
przez router - nie są ani arkuszem danych, ani błędem."""

KNOWN_SHEET_NAMES: tuple = (
    SHEET_DIFFUSION, SHEET_MIC_VISUAL, SHEET_MIC_OD, SHEET_MBC,
    SHEET_CONTROLS, SHEET_INSTRUCTIONS, SHEET_SETTINGS,
)
DATA_SHEET_NAMES: tuple = (SHEET_DIFFUSION, SHEET_MIC_VISUAL, SHEET_MIC_OD, SHEET_MBC)
"""DATA_SHEET_NAMES to arkusze niosące dane analizowalne per szczep (mają
kolumnę Bakteria). SHEET_CONTROLS ma inny schemat (Przebieg/Kontrola_wzrostu/
Kontrola_jalowosci/Inokulum_CFU_t0, bez kolumny Bakteria) i nie wchodzi do
mapy dostępności per szczep w tej fazie."""

ANALYSIS_DIFFUSION: str = "dyfuzja"
ANALYSIS_MIC: str = "mic"
ANALYSIS_MBC: str = "mbc"
ANALYSIS_TYPES: tuple = (ANALYSIS_DIFFUSION, ANALYSIS_MIC, ANALYSIS_MBC)
"""Typy analiz w mapie dostępności per szczep (utils.route_workbook).
ANALYSIS_MIC pokrywa łącznie SHEET_MIC_VISUAL i SHEET_MIC_OD - w tej fazie
routera oba tylko zasilają tę samą flagę dostępności; rozstrzygnięcie który
tryb ma pierwszeństwo, gdy oba są obecne dla tego samego szczepu+substancji,
jest odłożone do kolejnej fazy (patrz utils.route_workbook 'warnings')."""

# ============================================================
# MIC ZE STUDZIENEK (Faza 2: silnik per-wiersz, arkusze MIC_wizualny/MIC_OD)
# ============================================================

COL_RUN: str = "Przebieg"
COL_STEZ_S1: str = "Stez_S1"
COL_DILUTION_FACTOR: str = "Wsp_rozc"
"""Kolumny meta wspólne dla MIC_wizualny/MIC_OD/MBC_posiew. Bakteria/
Substancja/Typ/Jednostka/Rep_biologiczna/Rep_techniczna używają dokładnie
tych samych nazw i stałych co arkusz Dane_dyfuzja (COL_SUBSTANCE, COL_TYPE,
COL_UNIT, COL_REP_BIO, COL_REP_TECH) - te schematy dzielą nazewnictwo."""

WELL_COUNT: int = 10
WELL_COLUMNS: tuple = tuple(f"S{i}" for i in range(1, WELL_COUNT + 1))
"""Nazwy kolumn studzienek S1..S10. S1 = najwyższe testowane stężenie;
stężenie studzienki n = Stez_S1 / (Wsp_rozc ** (n-1)) - malejący szereg
rozcieńczeń w kolejności S1->S10."""

WELL_STATUS_GROWTH: str = "wzrost"
WELL_STATUS_NO_GROWTH: str = "brak"
"""Rozpoznawane (case-insensitive, po przycięciu białych znaków) wartości
tekstowe studzienek w arkuszu MIC_wizualny."""

MIC_OD_GROWTH_THRESHOLD: float = 0.10
"""Próg względny odróżniający "wzrost" od "brak" w MIC_OD:
procent_wzrostu = (OD_studzienki - tło) / (OD_kontroli_wzrostu - tło);
studzienka = "brak", gdy procent_wzrostu < próg. 10% to standardowy,
konserwatywny próg spotykany w kolorymetrycznych/OD-owych odczytach
wzrostu (odcina szum odczytu i nieswoiste zmętnienie, nie odcina
częściowo zahamowanego, ale realnie rosnącego inokulum)."""

MIC_OD_MIN_GROWTH_SIGNAL: float = 0.10
"""Minimalna wymagana różnica (Kontrola_wzrostu - Kontrola_jalowosci) w
jednostkach OD, żeby uznać przebieg za ważny. Typowy odczyt OD pustej/
jałowej studzienki mieści się w szumie rzędu 0.02-0.05 (błąd czytnika
płytek); różnica poniżej 0.10 OD oznacza, że albo kontrola wzrostu
realnie nie urosła, albo kontrola jałowości jest podejrzanie wysoka
(prawdopodobne skażenie) - w obu przypadkach różnicowanie wzrost/brak
w tym przebiegu nie jest wiarygodne. To pragmatyczny, konfigurowalny
próg bezpieczeństwa, nie wartość z walidowanego protokołu klinicznego."""

COL_GROWTH_CONTROL: str = "Kontrola_wzrostu"
COL_STERILITY_CONTROL: str = "Kontrola_jalowosci"
COL_INOCULUM: str = "Inokulum_CFU_t0"
"""Kolumny arkusza Kontrole (SHEET_CONTROLS), łączonego po COL_RUN. Dla MIC
istotne są COL_GROWTH_CONTROL (OD kontroli wzrostu) i COL_STERILITY_CONTROL
(OD kontroli jałowości, pełni rolę "tła"). COL_INOCULUM nie jest jeszcze
używana w Fazie 2."""

MIC_STATUS_OK: str = "ok"
MIC_STATUS_INVALID_RUN: str = "nieważny"
MIC_STATUS_NEEDS_REVIEW: str = "wymaga_weryfikacji"
MIC_STATUS_CENSORED_LOW: str = "cenzurowane_dol"
MIC_STATUS_CENSORED_HIGH: str = "cenzurowane_gora"
MIC_STATUS_MISSING_CONTROLS: str = "brak_kontroli"
"""Kody statusu wyniku MIC per wiersz (mic_logic.py). Każdy niepowodzenie/
zastrzeżenie ma odpowiadający, czytelny "reason" - żaden status nie jest
ciche (patrz mic_logic.compute_mic_for_row*)."""

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
