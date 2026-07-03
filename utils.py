import re
import numpy as np
import pandas as pd
from scipy import stats
from config import (
    COL_GROUP, COL_MEASUREMENT, COL_BACT_SUBSTRING,
    COHENS_D_SMALL, COHENS_D_MEDIUM, COHENS_D_LARGE,
    DIXON_Q90_CRITICAL,
    CONCENTRATION_SEARCH_PATTERN, CONCENTRATION_STRIP_PATTERN,
    NEGATIVE_CONTROL_SIGNALS, POSITIVE_CONTROL_SIGNALS, KNOWN_ANTIBIOTIC_SUBSTRINGS,
    NEW_FORMAT_SHEET_NAME, COL_SUBSTANCE, COL_CONCENTRATION, COL_UNIT, COL_TYPE,
    COL_REP_BIO, COL_REP_TECH, TYPE_NEG_CONTROL, TYPE_POS_CONTROL, VALID_TYPES,
    INTERNAL_TYPE_COL, INTERNAL_SUBSTANCE_COL, INTERNAL_CONC_COL, INTERNAL_UNIT_COL,
    INTERNAL_REP_BIO_COL, INTERNAL_REP_TECH_COL,
)

# --- SORTOWANIE I PARSOWANIE ---
def smart_sort_key(group_name):
    """Sortowanie naturalne dla nazw grup z liczbami."""
    match = re.match(r"(.+?)\s*\(([\d,.]+)\s*(.+?)\)", group_name)
    if match:
        name = match.group(1).strip()
        val_str = match.group(2).replace(',', '.')
        try: val = float(val_str)
        except ValueError: val = 0.0
        return (name, val)
    else: return (group_name, 0.0)

def find_bacteria_column(df):
    """
    Zwraca nazwę pierwszej kolumny zawierającej COL_BACT_SUBSTRING,
    albo None, jeśli żadna kolumna nie pasuje. Dopasowanie jest
    wrażliwe na wielkość liter (case-sensitive).
    """
    return next((c for c in df.columns if COL_BACT_SUBSTRING in c), None)

def parse_concentration(group_name):
    """Wyciąganie stężenia i jednostki z nazwy grupy."""
    match = re.search(CONCENTRATION_SEARCH_PATTERN, group_name)
    if match:
        val_str = match.group(1).replace(',', '.')
        try:
            conc = float(val_str)
            unit = match.group(2)
            substance = re.sub(CONCENTRATION_STRIP_PATTERN, "", group_name).strip()
            return substance, conc, unit
        except ValueError: return None, None, None
    return None, None, None

# --- WYKRYWANIE KONTROLI NEGATYWNEJ (grupa referencyjna) ---
def select_negative_control(groups):
    """
    Próbuje jednoznacznie wskazać grupę kontroli NEGATYWNEJ spośród `groups`.

    Zwraca (grupa, ambiguous):
        - (nazwa_grupy, False) - dokładnie jedna jednoznaczna grupa negatywna.
        - (None, True)         - brak jednoznacznej grupy (0 lub >1 kandydatów),
                                  wymaga ręcznego wyboru przez użytkownika.

    Grupa pasująca do POSITIVE_CONTROL_SIGNALS lub KNOWN_ANTIBIOTIC_SUBSTRINGS
    nigdy nie jest zwracana, nawet jeśli zawiera też sygnał negatywny (np.
    "Kontrola (+) Ampicylina" zawiera "kontrol" ale jest kontrolą pozytywną).
    """
    def is_positive(name):
        low = name.lower()
        return (any(s in low for s in POSITIVE_CONTROL_SIGNALS)
                or any(s in low for s in KNOWN_ANTIBIOTIC_SUBSTRINGS))

    def is_negative(name):
        low = name.lower()
        return any(s in low for s in NEGATIVE_CONTROL_SIGNALS)

    candidates = [g for g in groups if is_negative(g) and not is_positive(g)]

    if len(candidates) == 1:
        return candidates[0], False
    return None, True

def select_reference_group(df_subset, group_col=COL_GROUP, type_col=INTERNAL_TYPE_COL):
    """
    Wybiera grupę referencyjną (kontrolę negatywną) dla wierszy `df_subset`
    (typowo: dane jednego szczepu bakterii).

    Gdy kolumna `type_col` (nowy format - wartości Typ z pliku) jest obecna
    i ma jakiekolwiek nie-puste wartości, referencja jest wyznaczana
    JEDNOZNACZNIE z Typu ("Kontrola negatywna"), NIGDY z nazwy grupy -
    "Kontrola pozytywna" nigdy nie jest zwracana jako referencja.

    W przeciwnym razie (stary format, brak Typu) używa dotychczasowej,
    opartej na nazwie grupy logiki select_negative_control.

    Zwraca (grupa, ambiguous) - patrz select_negative_control.
    """
    if type_col in df_subset.columns and df_subset[type_col].notna().any():
        neg_groups = sorted(df_subset.loc[df_subset[type_col] == TYPE_NEG_CONTROL, group_col].unique())
        if len(neg_groups) == 1:
            return neg_groups[0], False
        return None, True

    groups = sorted(df_subset[group_col].unique(), key=smart_sort_key)
    return select_negative_control(groups)

# --- STATYSTYKA: EFFECT SIZE ---
def calculate_cohens_d(group1_data, group2_data):
    n1, n2 = len(group1_data), len(group2_data)
    # Wariancja wewnątrzgrupowa (a więc i pooled SD) jest niepoliczalna przy
    # n<2 (np. brak replikacji biologicznej, n_bio=1) - NaN, nie 0.0: d=0
    # oznaczałoby fałszywie "brak różnicy", podczas gdy w rzeczywistości po
    # prostu nie da się tego policzyć z tych danych.
    if n1 < 2 or n2 < 2: return np.nan

    var1 = np.var(group1_data, ddof=1)
    var2 = np.var(group2_data, ddof=1)
    
    s_pooled = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    mean1 = np.mean(group1_data)
    mean2 = np.mean(group2_data)

    if s_pooled == 0:
         return 0.0 
    
    return (mean1 - mean2) / s_pooled

def get_effect_size_interpretation(d):
    if pd.isna(d): return "nieokreślony (brak replikacji)"
    d = abs(d)
    if d < COHENS_D_SMALL: return "znikomy"
    elif d < COHENS_D_MEDIUM: return "mały"
    elif d < COHENS_D_LARGE: return "średni"
    else: return "DUŻY"

# --- STATYSTYKA: OUTLIERS (DIXON LOGIC) ---
def find_outliers_dixon(df):
    """Zwraca listę wykrytych outlierów (logika bez GUI)."""
    detected = []

    for group in df[COL_GROUP].unique():
        values = df[df[COL_GROUP] == group][COL_MEASUREMENT].values
        values = sorted(values)
        n = len(values)
        if n < 3 or n > 10: continue
        r = values[-1] - values[0]
        if r == 0: continue

        gap_low = values[1] - values[0]
        q_calc_low = gap_low / r
        gap_high = values[-1] - values[-2]
        q_calc_high = gap_high / r
        q_crit = DIXON_Q90_CRITICAL.get(n, 0.941)
        
        if q_calc_low > q_crit: 
            detected.append({'group': group, 'value': values[0], 'others': str(values[1:])})
        if q_calc_high > q_crit:
            detected.append({'group': group, 'value': values[-1], 'others': str(values[:-1])})

    return detected

# --- WALIDACJA STRUKTURY PLIKU EXCEL ---
def validate_excel_structure(df):
    """
    Strukturalna walidacja wczytanego DataFrame (bez sprawdzania typów / wartości).
    Obsługuje zarówno stary format (kolumna 'Grupa') jak i nowy format
    (kolumna 'Substancja' + 'Stezenie' + 'Jednostka' + 'Typ') - wykrywane
    niezależnie kolumna-po-kolumnie, patrz build_internal_representation.
    Zwraca (is_valid, error_messages) — error_messages to lista polskich komunikatów.
    """
    errors = []

    if find_bacteria_column(df) is None:
        errors.append(f"Brak wymaganej kolumny z nazwą bakterii (kolumny zawierającej w nazwie {COL_BACT_SUBSTRING!r}).")

    has_grupa = COL_GROUP in df.columns
    has_substancja = COL_SUBSTANCE in df.columns
    if not has_grupa and not has_substancja:
        errors.append(
            f"Brak kolumny {COL_GROUP!r} (stary format) ani {COL_SUBSTANCE!r} (nowy format) "
            "- nie można zbudować grup do porównania."
        )

    if COL_MEASUREMENT not in df.columns:
        errors.append(f"Brak wymaganej kolumny {COL_MEASUREMENT!r}.")

    if COL_CONCENTRATION in df.columns:
        if COL_SUBSTANCE not in df.columns:
            errors.append(f"Kolumna {COL_CONCENTRATION!r} wymaga też kolumny {COL_SUBSTANCE!r}.")
        if COL_UNIT not in df.columns:
            errors.append(f"Kolumna {COL_CONCENTRATION!r} wymaga też kolumny {COL_UNIT!r}.")

    if len(df) == 0:
        errors.append("Plik nie zawiera żadnych wierszy danych.")

    return (len(errors) == 0, errors)

# --- WALIDACJA WARTOŚCI KOMÓREK ---
def validate_excel_data(df, bacteria_col):
    """
    Walidacja wartości komórek (po walidacji strukturalnej). Sprawdza kolumny
    Typ/Stezenie tylko jeśli są obecne (nowy format) - stary format bez tych
    kolumn zachowuje się dokładnie jak dotychczas.
    Zwraca (cleaned_df, rejected_rows_info), gdzie rejected_rows_info to lista
    tupli (excel_row_number, reason_polish). Numer wiersza = indeks pandas + 2
    (zero-indexing + wiersz nagłówka), aby odpowiadać numeracji w Excelu.
    """
    rejected = []
    valid_indices = []

    has_grupa = COL_GROUP in df.columns
    has_substancja = COL_SUBSTANCE in df.columns
    has_type = COL_TYPE in df.columns
    has_conc = COL_CONCENTRATION in df.columns

    for idx in df.index:
        excel_row = idx + 2
        reason = None

        val = df.at[idx, COL_MEASUREMENT]
        if pd.isna(val):
            reason = f"Brak wartości w kolumnie {COL_MEASUREMENT!r}."
        else:
            try:
                fval = float(val)
                if not np.isfinite(fval):
                    reason = f"Nieprawidłowa wartość {COL_MEASUREMENT!r} ({val!r}) — wartość nieskończona."
                elif fval < 0:
                    reason = f"Nieprawidłowa wartość {COL_MEASUREMENT!r} ({val!r}) — nie może być wartością ujemną (0 oznacza brak strefy zahamowania i jest dopuszczalne)."
                elif fval > 100:
                    reason = f"Podejrzana wartość {COL_MEASUREMENT!r} ({val!r}) — większa niż 100 mm (prawdopodobny błąd pomiaru)."
            except (ValueError, TypeError):
                reason = f"Nieprawidłowa wartość {COL_MEASUREMENT!r} ({val!r}) — nie można przekształcić na liczbę."

        if reason is None:
            if has_grupa:
                grp = df.at[idx, COL_GROUP]
                if pd.isna(grp) or (isinstance(grp, str) and grp.strip() == ""):
                    reason = f"Pusta wartość w kolumnie {COL_GROUP!r}."
            elif has_substancja:
                sub = df.at[idx, COL_SUBSTANCE]
                if pd.isna(sub) or (isinstance(sub, str) and sub.strip() == ""):
                    reason = f"Pusta wartość w kolumnie {COL_SUBSTANCE!r}."

        if reason is None:
            bact_val = df.at[idx, bacteria_col]
            if pd.isna(bact_val) or (isinstance(bact_val, str) and bact_val.strip() == ""):
                reason = f"Pusta wartość w kolumnie bakterii ('{bacteria_col}')."

        if reason is None and has_type:
            typ_val = df.at[idx, COL_TYPE]
            if pd.isna(typ_val) or str(typ_val).strip() not in VALID_TYPES:
                reason = (
                    f"Nieprawidłowa wartość {COL_TYPE!r} ({typ_val!r}) — dozwolone: "
                    + ", ".join(repr(t) for t in VALID_TYPES) + "."
                )

        if reason is None and has_conc:
            conc_val = df.at[idx, COL_CONCENTRATION]
            if pd.notna(conc_val):
                try:
                    fconc = float(conc_val)
                    if not np.isfinite(fconc) or fconc < 0:
                        reason = f"Nieprawidłowa wartość {COL_CONCENTRATION!r} ({conc_val!r}) — musi być liczbą nieujemną."
                except (ValueError, TypeError):
                    reason = f"Nieprawidłowa wartość {COL_CONCENTRATION!r} ({conc_val!r}) — nie można przekształcić na liczbę."
            # Puste Stężenie jest dopuszczalne (np. dla kontroli bez określonego stężenia).

        if reason is None:
            valid_indices.append(idx)
        else:
            rejected.append((excel_row, reason))

    cleaned_df = df.loc[valid_indices].reset_index(drop=True)
    return cleaned_df, rejected

# --- FORMAT NORMALIZATION (stary <-> nowy format wejściowy) ---
def _fmt_num(value):
    """Czytelny format liczby w etykiecie grupy: 50 zamiast 50.0, 0.5 zostaje 0.5."""
    try:
        f = float(value)
    except (ValueError, TypeError):
        return str(value)
    if f == int(f):
        return str(int(f))
    return f"{f:g}"

def build_internal_representation(df):
    """
    Buduje wewnętrzną (dotychczasową) reprezentację danych z surowego,
    zwalidowanego DataFrame - niezależnie czy pochodzi ze starego czy
    nowego formatu wejściowego.

    Zwraca (df2, format_info):
      - df2 ma zawsze kolumnę COL_GROUP (istniejącą albo zsyntetyzowaną
        z Substancja/Stezenie/Jednostka/Typ) oraz kolumny wewnętrzne
        INTERNAL_TYPE_COL/INTERNAL_SUBSTANCE_COL/INTERNAL_CONC_COL/
        INTERNAL_UNIT_COL/INTERNAL_REP_BIO_COL/INTERNAL_REP_TECH_COL,
        wypełnione z nowego formatu gdy dostępne, inaczej NaN / wartościami
        domyślnymi (patrz niżej).
      - format_info: dict {'has_type', 'has_conc', 'has_reps'} do celów
        logowania/UI (np. "wykryto nowy format").
    """
    df = df.copy()
    has_type = COL_TYPE in df.columns
    has_conc = COL_CONCENTRATION in df.columns and COL_SUBSTANCE in df.columns and COL_UNIT in df.columns
    has_reps = COL_REP_BIO in df.columns

    df[INTERNAL_TYPE_COL] = df[COL_TYPE] if has_type else np.nan

    if has_conc:
        df[INTERNAL_SUBSTANCE_COL] = df[COL_SUBSTANCE]
        df[INTERNAL_CONC_COL] = pd.to_numeric(df[COL_CONCENTRATION], errors='coerce')
        df[INTERNAL_UNIT_COL] = df[COL_UNIT]
    elif COL_SUBSTANCE in df.columns:
        # Stezenie/Jednostka brakuje, ale Substancja jest - traktujemy jak
        # substancję bez określonego stężenia (np. sama nazwa kontroli).
        df[INTERNAL_SUBSTANCE_COL] = df[COL_SUBSTANCE]
        df[INTERNAL_CONC_COL] = np.nan
        df[INTERNAL_UNIT_COL] = np.nan
    else:
        df[INTERNAL_SUBSTANCE_COL] = np.nan
        df[INTERNAL_CONC_COL] = np.nan
        df[INTERNAL_UNIT_COL] = np.nan

    if COL_GROUP not in df.columns:
        def _label(row):
            typ = row[INTERNAL_TYPE_COL] if has_type else None
            sub = row[INTERNAL_SUBSTANCE_COL] if pd.notna(row[INTERNAL_SUBSTANCE_COL]) else None
            conc = row[INTERNAL_CONC_COL]
            unit = row[INTERNAL_UNIT_COL]

            if typ in (TYPE_NEG_CONTROL, TYPE_POS_CONTROL):
                return f"{typ} ({sub})" if sub else typ
            if sub and pd.notna(conc) and pd.notna(unit) and str(unit).strip() != "":
                return f"{sub} ({_fmt_num(conc)} {unit})"
            return sub or "?"
        df[COL_GROUP] = df.apply(_label, axis=1)

    bact_col = find_bacteria_column(df)
    if has_reps:
        df[INTERNAL_REP_BIO_COL] = df[COL_REP_BIO]
        df[INTERNAL_REP_TECH_COL] = df[COL_REP_TECH] if COL_REP_TECH in df.columns else 1
    else:
        # Brak kolumn powtórzeń: cały wiersz danej grupy to JEDNO powtórzenie
        # biologiczne (n_bio=1); kolejne wiersze w obrębie (bakteria, grupa)
        # traktujemy jako powtórzenia TECHNICZNE tego jedynego powtórzenia
        # biologicznego (indeksowane 1..k w kolejności występowania).
        df[INTERNAL_REP_BIO_COL] = 1
        df[INTERNAL_REP_TECH_COL] = df.groupby([bact_col, COL_GROUP]).cumcount() + 1

    format_info = {"has_type": has_type, "has_conc": has_conc, "has_reps": has_reps}
    return df, format_info

def read_excel_any_format(path):
    """
    Wczytuje plik Excel: jeśli zawiera arkusz NEW_FORMAT_SHEET_NAME ('Dane'),
    czyta ten arkusz; w przeciwnym razie czyta domyślny/pierwszy arkusz
    (zachowanie identyczne jak dotychczas dla starych plików).
    Przycina białe znaki w nazwach kolumn i w wartościach tekstowych.
    Może rzucić te same wyjątki co pandas.read_excel (FileNotFoundError,
    PermissionError, ValueError, itp.) - wywołujący łapie je jak dotychczas.
    """
    xls = pd.ExcelFile(path)
    sheet = NEW_FORMAT_SHEET_NAME if NEW_FORMAT_SHEET_NAME in xls.sheet_names else 0
    df = pd.read_excel(xls, sheet_name=sheet)
    df.columns = df.columns.str.strip()
    for col in df.select_dtypes(['object']).columns:
        df[col] = df[col].str.strip()
    return df

def validate_and_normalize(df):
    """
    Waliduje strukturalnie i wartościowo (stary LUB nowy format, wykrywane
    niezależnie kolumna-po-kolumnie) i buduje wewnętrzną reprezentację.

    Zwraca (df_internal, bacteria_col, rejected, struct_errors, format_info).
    Gdy struct_errors jest niepuste, df_internal/bacteria_col są None i
    wywołujący powinien przerwać wczytywanie i pokazać struct_errors.
    """
    is_valid, struct_errors = validate_excel_structure(df)
    if not is_valid:
        return None, None, [], struct_errors, {}

    bacteria_col = find_bacteria_column(df)
    cleaned_df, rejected = validate_excel_data(df, bacteria_col)
    df_internal, format_info = build_internal_representation(cleaned_df)
    return df_internal, bacteria_col, rejected, [], format_info

# --- POWTÓRZENIA: AGREGACJA TECHNICZNYCH, PODSUMOWANIE BIOLOGICZNYCH ---
def aggregate_technical_replicates(df, bacteria_col):
    """
    Uśrednia powtórzenia TECHNICZNE do jednej wartości na każde
    (bacteria_col, COL_GROUP, INTERNAL_REP_BIO_COL) - PRZED testami
    istotności, żeby n grupy w teście statystycznym = liczba powtórzeń
    BIOLOGICZNYCH (n_bio), a nie liczba surowych wierszy. Liczenie
    powtórzeń technicznych jako niezależnych obserwacji byłoby
    pseudoreplikacją (sztucznie zawyżałoby moc testu / zaniżało p-value).

    Dla starego formatu (brak kolumn powtórzeń) każda grupa ma dokładnie
    jedno powtórzenie biologiczne obejmujące wszystkie jej surowe wiersze
    (patrz build_internal_representation) - więc tutaj kolapsuje do
    JEDNEGO wiersza na grupę (n_bio=1), co jest zamierzone: to właśnie
    uwidacznia, że stare dane nie mają udokumentowanej replikacji
    biologicznej.

    Zwraca DataFrame z jednym wierszem na (bacteria_col, COL_GROUP,
    INTERNAL_REP_BIO_COL): COL_MEASUREMENT to średnia techniczna, kolumna
    'n_tech' (liczba uśrednionych powtórzeń technicznych w tym wierszu),
    oraz metadane (_Typ/_Substancja/_Stezenie/_Jednostka) wzięte jako
    pierwsza wartość (są stałe w obrębie klucza z założenia).
    """
    key = [bacteria_col, COL_GROUP, INTERNAL_REP_BIO_COL]
    agg_spec = {COL_MEASUREMENT: 'mean'}
    for meta_col in (INTERNAL_TYPE_COL, INTERNAL_SUBSTANCE_COL, INTERNAL_CONC_COL, INTERNAL_UNIT_COL):
        if meta_col in df.columns:
            agg_spec[meta_col] = 'first'

    df_bio = df.groupby(key, as_index=False).agg(agg_spec)
    n_tech = df.groupby(key).size().reset_index(name='n_tech')
    df_bio = df_bio.merge(n_tech, on=key)
    return df_bio

def build_group_summary(df_bio, group_col=COL_GROUP):
    """
    Tabela opisowa per grupa na podstawie DANYCH JUŻ ZAGREGOWANYCH do
    powtórzeń biologicznych (df_bio, patrz aggregate_technical_replicates):
      - n_bio: liczba powtórzeń biologicznych (liczba wierszy df_bio per grupa)
      - n_tech: łączna liczba surowych pomiarów technicznych złożonych na n_bio
      - mean: średnia z powtórzeń biologicznych
      - sd_bio: odchylenie standardowe MIĘDZY powtórzeniami biologicznymi -
        NaN gdy n_bio<2 (nie ma czego liczyć - pandas .std(ddof=1) na
        jednym punkcie zwraca NaN). NIGDY nie jest to rozproszenie
        techniczne - te dwa rodzaje zmienności nie są tu mieszane.
    """
    means = df_bio.groupby(group_col)[COL_MEASUREMENT].mean().reset_index(name='mean')
    n_bio = df_bio.groupby(group_col)[COL_MEASUREMENT].size().reset_index(name='n_bio')
    sd_bio = df_bio.groupby(group_col)[COL_MEASUREMENT].std().reset_index(name='sd_bio')
    n_tech = df_bio.groupby(group_col)['n_tech'].sum().reset_index(name='n_tech')
    summary = means.merge(n_bio, on=group_col).merge(sd_bio, on=group_col).merge(n_tech, on=group_col)
    return summary

def has_low_n_bio(df_bio, group_col=COL_GROUP):
    """
    True jeśli KTÓRAKOLWIEK grupa w df_bio (już zagregowanym do powtórzeń
    biologicznych) ma n_bio<2 - tzn. porównanie obejmuje przynajmniej jedną
    grupę bez replikacji biologicznej, więc istotność wyniku dla całej
    analizy jest tylko orientacyjna.
    """
    return bool((df_bio.groupby(group_col)[COL_MEASUREMENT].size() < 2).any())