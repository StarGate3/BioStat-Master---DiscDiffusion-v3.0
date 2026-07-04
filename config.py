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
jednostkach OD, żeby uznać KONTROLĘ WZROSTU za wiarygodną (bakteria
realnie urosła). To TYLKO jeden z dwóch niezależnych warunków ważności
przebiegu - patrz MIC_OD_STERILITY_MAX dla drugiego (czystość podłoża).
Rozdzielenie jest celowe: wysoka kontrola wzrostu nie może "maskować"
skażonej kontroli jałowości, bo różnica mogłaby wciąż wyjść wystarczająco
duża. Typowy odczyt OD pustej/jałowej studzienki mieści się w szumie
rzędu 0.02-0.05 (błąd czytnika płytek), więc 0.10 to pragmatyczny,
konfigurowalny próg bezpieczeństwa, nie wartość z walidowanego protokołu
klinicznego."""

MIC_OD_STERILITY_MAX: float = 0.20
"""Maksymalna dopuszczalna wartość OD kontroli jałowości (Kontrola_jalowosci),
NIEZALEŻNIE od poziomu kontroli wzrostu - powyżej tego progu podłoże uznaje
się za skażone. Czysta/jałowa studzienka (samo podłoże, bez inokulum) daje
zwykle OD rzędu 0.05-0.15 (mętność podłoża + szum czytnika); 0.20 zostawia
margines nad typową zmiennością czystego blanku, jednocześnie wciąż łapiąc
realne skażenie, które zwykle podnosi OD dużo wyżej (często porównywalnie
do samej kontroli wzrostu, rzędu 0.3-1.0+). Podobnie jak
MIC_OD_MIN_GROWTH_SIGNAL - to pragmatyczny próg narzędzia przesiewowego,
nie wartość z walidowanego protokołu klinicznego."""

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
# MIC: AGREGACJA POWTÓRZEŃ (Faza 3)
# ============================================================

MIC_STATUS_NO_DATA: str = "brak_danych"
"""Status na poziomie AGREGATU (powtórzenie biologiczne albo grupa
Bakteria x Substancja) - gdy po wykluczeniu nieważnych/bezwartościowych
powtórzeń technicznych (albo biologicznych) nie zostaje ani jedna wartość
do zagregowania. Nie występuje na poziomie pojedynczego wiersza (Faza 2),
tylko po agregacji (mic_logic.aggregate_technical_to_biological /
summarize_mic_group)."""

MIC_LOW_N_BIO_WARNING: str = (
    "Brak replikacji biologicznej (n_bio=1) - wynik orientacyjny, nie "
    "potwierdzony niezależnymi powtórzeniami."
)
"""Ostrzeżenie dołączane do podsumowania grupy MIC (Bakteria x Substancja),
gdy n_bio<2, spójne w duchu z ostrzeżeniem n_bio=1 z modułu dyfuzji
(logika/treść analogiczna, ale osobna stała - moduły MIC i dyfuzja nie
współdzielą kodu prezentacji)."""

# ============================================================
# MIC: PORÓWNANIA MIĘDZY GRUPAMI (Faza 4)
# ============================================================

MIC_MEANINGFUL_DILUTION_DIFF: float = 2.0
"""Minimalna różnica median (w krokach rozcieńczenia = log2 stężenia)
uznawana za metodologicznie sensowną ("meaningful"). 2 kroki = różnica
4-krotna. Uzasadnienie: pojedynczy krok rozcieńczenia (2-krotny) mieści
się w typowej zmienności metody mikrorozcieńczeń (assay reproducibility
jest zwykle podawana jako ±1 rozcieńczenie w standardowych protokołach
mikrobiologicznych), więc różnica 1 kroku nie jest odróżnialna od zwykłego
szumu pomiarowego. Różnica ≥2 kroków (4-krotna) przekracza tę typową
zmienność i jest powszechnie przyjmowanym progiem "realnej" różnicy
aktywności w literaturze porównań MIC. Próg czysto metodologiczny, nie
kliniczny/regulacyjny."""

MIC_MAX_CENSORED_FRACTION: float = 0.50
"""Maksymalny dopuszczalny odsetek wartości cenzurowanych (≤min lub >max)
w KTÓREJKOLWIEK porównywanej grupie, powyżej którego test istotności
(Warstwa 2) jest wygaszany (zostaje tylko opis, Warstwa 1). Uzasadnienie:
gdy połowa lub więcej odczytów w grupie to otwarte granice, a nie
dokładne wartości, obserwowana różnica rang jest w dużej mierze efektem
KONWENCJI wiązania rang cenzurowanych wartości (patrz
mic_logic._censoring_surrogate), a nie realnego sygnału biologicznego -
p-value w takiej sytuacji sugerowałby pewność, jakiej dane nie
uzasadniają. 50% to prosty, łatwy do uzasadnienia próg ("większość
danych w grupie to nie są liczby") - nie wartość z walidowanego
protokołu."""

MIN_N_BIO_FOR_COMPARISON: int = 1
"""Minimalna liczba powtórzeń biologicznych z wartością, żeby grupa mogła
w ogóle wejść do porównania (choćby tylko opisowego). Poniżej tego progu
(n_bio=0, czyli grupa bez żadnej wartości MIC) nie da się policzyć nawet
mediany/zakresu - porównanie z udziałem takiej grupy jest blokowane z
jawnym powodem. n_bio=1 SPEŁNIA ten minimalny próg (bo 1 >= 1), więc nie
jest tu blokowany - dostaje własne, osobne traktowanie: p-value jest
liczone normalnie, ale z banerem ostrzegawczym (patrz MIC_LOW_N_BIO_WARNING
i mic_logic.compare_mic_groups), spójnie z resztą programu (moduł
dyfuzji)."""

# ============================================================
# MBC (Faza MBC): odczyt z posiewu + iloraz MBC/MIC
# ============================================================

MBC_REDUCTION_THRESHOLD: float = 0.999
"""Próg redukcji CFU (1 - CFU_studzienki/CFU_t0) uznawany za spełnienie
kryterium bójczego. 99.9% = redukcja 3-log10 - to NIE jest pragmatyczny
próg autorski jak inne stałe w tym pliku, tylko STANDARDOWA, powszechnie
przyjęta w mikrobiologii klinicznej definicja "bakteriobójczości" (MBC)
w testach mikrorozcieńczeń - patrz standardowe protokoły oznaczania MBC."""

MBC_MIC_BACTERICIDAL_MAX_D: int = 2
MBC_MIC_BACTERIOSTATIC_MIN_D: int = 3
"""Progi klasyfikacji ilorazu MBC/MIC na SKALI RÓŻNICY INDEKSÓW ROZCIEŃCZEŃ
d = log2(MBC) - log2(MIC) (iloraz do wyświetlenia = 2**d): d<=2 (iloraz
<=4) -> "bakteriobójcze"; d>=3 (iloraz >=8) -> "bakteriostatyczne". Progi
podane wprost w specyfikacji zadania - odpowiadają powszechnie cytowanej
klasycznej konwencji farmakologicznej (MBC/MIC <=4 = bakteriobójcze).

UWAGA (audyt 1.4): te progi na skali log2 są zdefiniowane w literaturze
DLA SERII DWUKROTNYCH ROZCIEŃCZEŃ (Wsp_rozc=2) - tylko wtedy każdy możliwy
krok d jest liczbą całkowitą, więc próg "d<=2 albo d>=3" nigdy nie trafia
w "martwą strefę" między nimi. Dla innego współczynnika (np. 5-krotnego,
d=log2(5)=2.32) próg ten wypadałby w tej martwej strefie mimo w pełni
dokładnego, niecenzurowanego pomiaru - dlatego mic_logic.compute_mbc_mic_ratio
i mic_logic.summarize_mbc_mic_ratio liczą MIC/MBC/iloraz normalnie dla
DOWOLNEGO Wsp_rozc, ale samą KLASYFIKACJĘ bakteriobójcze/bakteriostatyczne
podają tylko, gdy Wsp_rozc == MBC_MIC_CLASSIFICATION_DILUTION_FACTOR."""

MBC_MIC_CLASSIFICATION_DILUTION_FACTOR: int = 2
"""Jedyny współczynnik rozcieńczenia (Wsp_rozc), dla którego klasyfikacja
bakteriobójcze/bakteriostatyczne (progi wyżej) jest w ogóle podawana -
patrz uzasadnienie w komentarzu do MBC_MIC_BACTERICIDAL_MAX_D. Dla innego
Wsp_rozc MIC, MBC i sam iloraz są liczone i pokazywane jak zawsze -
niedostępna jest WYŁĄCZNIE etykieta bakteriobójcze/bakteriostatyczne,
z jawnym powodem zamiast zgadywania."""

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
