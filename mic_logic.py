"""
Silnik odczytu MIC ze studzienek (Faza 2 modułu MIC/MBC).

Zakres tej fazy: dla KAŻDEGO pojedynczego wiersza (powtórzenia) arkusza
MIC_wizualny lub MIC_OD, wyznacz wartość MIC z serii studzienek S1..S10.
Nie ma tu jeszcze agregacji powtórzeń, statystyki międzygrupowej ani
wykresów - to kolejne fazy.

Podstawowy algorytm (compute_mic_from_wells) jest wspólny dla obu trybów;
różni je tylko sposób sprowadzenia surowej wartości studzienki do
booleana "wzrost/brak" (patrz classify_wizualny_well / classify_od_well).
"""
import pandas as pd
import utils
from config import (
    MIC_STATUS_OK, MIC_STATUS_NEEDS_REVIEW, MIC_STATUS_CENSORED_LOW, MIC_STATUS_CENSORED_HIGH,
    MIC_STATUS_INVALID_RUN, MIC_STATUS_MISSING_CONTROLS,
    COL_RUN, COL_STEZ_S1, COL_DILUTION_FACTOR, WELL_COLUMNS,
    WELL_STATUS_GROWTH, WELL_STATUS_NO_GROWTH, MIC_OD_GROWTH_THRESHOLD, MIC_OD_MIN_GROWTH_SIGNAL,
    COL_SUBSTANCE, COL_TYPE, COL_UNIT, COL_REP_BIO, COL_REP_TECH,
    COL_GROWTH_CONTROL, COL_STERILITY_CONTROL,
)


def compute_well_concentration(stez_s1, wsp_rozc, well_number):
    """
    Stężenie studzienki S{well_number} (well_number=1 dla S1, najwyższe
    stężenie, malejąco dla kolejnych numerów): Stez_S1 / (Wsp_rozc ** (n-1)).
    """
    return stez_s1 / (wsp_rozc ** (well_number - 1))


def compute_mic_from_wells(conc_growth_pairs):
    """
    Rdzeń algorytmu - niezależny od trybu (wizualny/OD) i od arkusza.

    conc_growth_pairs: lista (stężenie: float, wzrost: bool) TYLKO dla
    studzienek z rozpoznaną wartością (studzienki puste/nierozpoznane mają
    być odfiltrowane PRZED wywołaniem). Kolejność wejściowa dowolna -
    funkcja sama sortuje malejąco po stężeniu (S1 = najwyższe stężenie).

    wzrost=True  -> studzienka urosła (nieshamowana)
    wzrost=False -> studzienka "brak" (zahamowana)

    Zwraca dict:
      - mic_value (float | None): wartość liczbowa MIC, albo granica przy
        cenzurze (patrz "censored"), albo None gdy nie da się wyznaczyć.
      - status: jeden z MIC_STATUS_* (config.py).
      - reason (str): zawsze czytelny, nawet dla status="ok".
      - censored: "low" | "high" | None. Gdy nie None, mic_value to GRANICA,
        NIE dokładna wartość MIC - wyświetlaj zawsze z prefiksem ≤/> (patrz
        format_mic_display), nigdy jako zwykłą liczbę.

    Zasady:
      1. Cała seria "brak"  -> MIC <= najniższe testowane stężenie (censored="low").
      2. Cała seria "wzrost" -> MIC > najwyższe testowane stężenie (censored="high").
      3. Wzrost widoczny już w najwyższym testowanym stężeniu, a mimo to
         jakieś niższe stężenie pokazuje "brak" -> wynik niejednoznaczny,
         status wymaga_weryfikacji, mic_value=None (nie ma wiarygodnej
         wartości do podania).
      4. W pozostałych przypadkach: MIC = stężenie ostatniej studzienki
         "brak" WYSTĘPUJĄCEJ PRZED pierwszym napotkanym "wzrost" (skanując
         od najwyższego stężenia w dół). To jest CELOWO ta sama formuła
         niezależnie od tego, czy seria jest potem monotoniczna - wybrano
         ją jako regułę zachowawczą: nigdy nie zgłasza MIC niższego niż
         stężenie, przy którym już raz zaobserwowano wzrost. Gdy po tym
         pierwszym "wzrost" pojawi się jednak jeszcze jakieś "brak"
         (naruszenie monotoniczności - możliwy "skip well"), wartość MIC
         zostaje taka sama, ale status zmienia się na wymaga_weryfikacji
         z opisem, które studzienki są anomalią.
    """
    if not conc_growth_pairs:
        return {
            "mic_value": None,
            "status": MIC_STATUS_NEEDS_REVIEW,
            "reason": "Brak odczytanych (rozpoznanych) studzienek - nie można wyznaczyć MIC.",
            "censored": None,
        }

    ordered = sorted(conc_growth_pairs, key=lambda pair: pair[0], reverse=True)
    concs = [c for c, _ in ordered]
    growths = [g for _, g in ordered]

    if all(not g for g in growths):
        lowest = min(concs)
        return {
            "mic_value": lowest,
            "status": MIC_STATUS_CENSORED_LOW,
            "reason": (
                f"Brak wzrostu we wszystkich testowanych studzienkach - prawdziwe MIC może być "
                f"niższe niż najniższe testowane stężenie ({lowest:g})."
            ),
            "censored": "low",
        }

    if all(growths):
        highest = max(concs)
        return {
            "mic_value": highest,
            "status": MIC_STATUS_CENSORED_HIGH,
            "reason": (
                f"Wzrost we wszystkich testowanych studzienkach - prawdziwe MIC może być wyższe "
                f"niż najwyższe testowane stężenie ({highest:g})."
            ),
            "censored": "high",
        }

    first_growth_idx = growths.index(True)

    if first_growth_idx == 0:
        return {
            "mic_value": None,
            "status": MIC_STATUS_NEEDS_REVIEW,
            "reason": (
                "Wzrost zaobserwowano już przy najwyższym testowanym stężeniu, mimo że niższe "
                "stężenie(a) wykazały brak wzrostu - wynik niejednoznaczny, nie można wyznaczyć "
                "MIC standardową metodą (brak studzienki 'powyżej' punktu odcięcia)."
            ),
            "censored": None,
        }

    mic_value = concs[first_growth_idx - 1]
    is_monotonic = all(growths[first_growth_idx:])

    if is_monotonic:
        return {
            "mic_value": mic_value,
            "status": MIC_STATUS_OK,
            "reason": "Seria monotoniczna - MIC odczytane jednoznacznie.",
            "censored": None,
        }

    skip_concs = [concs[i] for i in range(first_growth_idx, len(growths)) if not growths[i]]
    skip_desc = ", ".join(f"{c:g}" for c in skip_concs)
    reason = (
        f"Niemonotoniczna seria: wzrost pojawił się już przy stężeniu {concs[first_growth_idx]:g}, "
        f"ale brak wzrostu wystąpił ponownie przy niższym stężeniu/stężeniach ({skip_desc}) - "
        f"możliwy błąd studzienki (skip well) lub skażenie. MIC odczytane ZACHOWAWCZO jako ostatnia "
        f"studzienka 'brak' PRZED pierwszym napotkanym 'wzrost' ({mic_value:g}), a NIE jako najniższe "
        f"'brak' w całej serii - wymaga ręcznej weryfikacji."
    )
    return {
        "mic_value": mic_value,
        "status": MIC_STATUS_NEEDS_REVIEW,
        "reason": reason,
        "censored": None,
    }


def format_mic_display(mic_value, censored):
    """
    Czytelna reprezentacja tekstowa mic_value/censored - NIGDY nie zwraca
    samej liczby granicznej bez prefiksu ≤/>, żeby cenzura nie została
    pomylona ze zwykłą wartością.
    """
    if mic_value is None:
        return "n/d"
    if censored == "low":
        return f"≤{mic_value:g}"
    if censored == "high":
        return f">{mic_value:g}"
    return f"{mic_value:g}"


# ============================================================
# KLASYFIKACJA POJEDYNCZEJ STUDZIENKI: surowa wartość -> wzrost (bool) | None
# ============================================================

def classify_wizualny_well(raw_value):
    """
    MIC_wizualny: studzienka to tekst "wzrost"/"brak" (case-insensitive, po
    przycięciu białych znaków). Zwraca True/False, albo None gdy pusta LUB
    nierozpoznana (wywołujący musi to odnotować, nie ignorować po cichu).
    """
    if pd.isna(raw_value):
        return None
    v = str(raw_value).strip().lower()
    if v == WELL_STATUS_GROWTH:
        return True
    if v == WELL_STATUS_NO_GROWTH:
        return False
    return None


def classify_od_well(od_value, kontrola_wzrostu, kontrola_jalowosci, threshold=MIC_OD_GROWTH_THRESHOLD):
    """
    MIC_OD: studzienka to liczba OD. Sprowadzana do wzrost/brak progiem
    względnym: procent_wzrostu = (OD - tło) / (Kontrola_wzrostu - tło),
    gdzie tło = Kontrola_jalowosci. "brak" gdy procent_wzrostu < threshold.
    Zwraca None gdy od_value brakuje, albo gdy mianownik jest niedodatni
    (nie powinno się zdarzyć po przejściu validate_run, ale zabezpieczenie
    przed dzieleniem przez zero/ujemną liczbę zamiast wyjątku).
    """
    if pd.isna(od_value):
        return None
    denom = kontrola_wzrostu - kontrola_jalowosci
    if denom <= 0:
        return None
    percent = (od_value - kontrola_jalowosci) / denom
    return percent >= threshold


# ============================================================
# KONTROLE (arkusz Kontrole, łączone po Przebieg)
# ============================================================

def lookup_controls(controls_df, przebieg):
    """
    Zwraca dict {Kontrola_wzrostu, Kontrola_jalowosci} dla danego Przebiegu,
    albo None gdy controls_df jest None, Przebieg jest pusty, albo nie ma
    dopasowania w arkuszu Kontrole.
    """
    if controls_df is None or pd.isna(przebieg):
        return None
    matches = controls_df[controls_df[COL_RUN] == przebieg]
    if matches.empty:
        return None
    row = matches.iloc[0]
    return {
        COL_GROWTH_CONTROL: row.get(COL_GROWTH_CONTROL),
        COL_STERILITY_CONTROL: row.get(COL_STERILITY_CONTROL),
    }


def validate_run(kontrola_wzrostu, kontrola_jalowosci):
    """
    Sprawdza, czy przebieg jest wiarygodny na podstawie kontroli OD:
    różnica (Kontrola_wzrostu - Kontrola_jalowosci) musi wynosić co najmniej
    MIC_OD_MIN_GROWTH_SIGNAL. Ten jeden test celowo pokrywa oba przypadki
    z zadania: kontrola wzrostu, która realnie nie urosła, ORAZ kontrola
    jałowości sugerująca skażenie - obie pchają tę różnicę w stronę zera.

    Zwraca (is_valid: bool, reason: str | None).
    """
    signal = kontrola_wzrostu - kontrola_jalowosci
    if signal < MIC_OD_MIN_GROWTH_SIGNAL:
        if kontrola_jalowosci >= kontrola_wzrostu:
            reason = (
                f"Kontrola jałowości ({kontrola_jalowosci:g} OD) >= kontrola wzrostu "
                f"({kontrola_wzrostu:g} OD) - podejrzenie skażenia podłoża. Przebieg nieważny."
            )
        else:
            reason = (
                f"Sygnał kontroli wzrostu po odjęciu tła ({signal:g} OD) poniżej minimum "
                f"({MIC_OD_MIN_GROWTH_SIGNAL:g} OD) - kontrola wzrostu nie urosła wystarczająco. "
                f"Przebieg nieważny."
            )
        return False, reason
    return True, None


# ============================================================
# PRZETWARZANIE WIERSZY ARKUSZA (MIC_wizualny / MIC_OD)
# ============================================================

def _base_fields(row, bact_col):
    return {
        "Przebieg": row.get(COL_RUN),
        "Bakteria": row.get(bact_col),
        "Substancja": row.get(COL_SUBSTANCE),
        "Typ": row.get(COL_TYPE),
        "Jednostka": row.get(COL_UNIT),
        "Rep_biologiczna": row.get(COL_REP_BIO),
        "Rep_techniczna": row.get(COL_REP_TECH),
    }


def _result(base, mic_value, status, reason, censored, wells):
    return {
        **base,
        "mic_value": mic_value,
        "status": status,
        "reason": reason,
        "censored": censored,
        "wells": wells,
    }


def _process_row(row, bact_col, well_cols, controls_df, mode):
    """
    Wspólna ścieżka dla jednego wiersza MIC_wizualny (mode='wizualny') albo
    MIC_OD (mode='od'). Zwraca dict wyniku (patrz _result / moduł docstring).
    """
    base = _base_fields(row, bact_col)
    przebieg = row.get(COL_RUN)
    ctrl = lookup_controls(controls_df, przebieg)
    has_numeric_ctrl = (
        ctrl is not None
        and pd.notna(ctrl.get(COL_GROWTH_CONTROL))
        and pd.notna(ctrl.get(COL_STERILITY_CONTROL))
    )

    if mode == "od" and not has_numeric_ctrl:
        return _result(
            base, None, MIC_STATUS_MISSING_CONTROLS,
            f"Brak wpisu kontroli (Kontrola_wzrostu/Kontrola_jalowosci) w arkuszu Kontrole dla "
            f"Przebiegu '{przebieg}' - wymagane do wyznaczenia tła i progu wzrostu OD.",
            None, [],
        )

    if has_numeric_ctrl:
        is_valid, invalid_reason = validate_run(ctrl[COL_GROWTH_CONTROL], ctrl[COL_STERILITY_CONTROL])
        if not is_valid:
            return _result(base, None, MIC_STATUS_INVALID_RUN, invalid_reason, None, [])

    stez_s1 = row.get(COL_STEZ_S1)
    wsp_rozc = row.get(COL_DILUTION_FACTOR)
    if pd.isna(stez_s1) or pd.isna(wsp_rozc):
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"Brak {COL_STEZ_S1!r} lub {COL_DILUTION_FACTOR!r} - nie można wyliczyć stężeń studzienek.",
            None, [],
        )
    if wsp_rozc <= 1:
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"{COL_DILUTION_FACTOR!r}={wsp_rozc:g} musi być > 1 (malejący szereg rozcieńczeń).",
            None, [],
        )

    conc_growth_pairs = []
    wells_detail = []
    parse_notes = []
    for i, col in enumerate(well_cols, start=1):
        raw_value = row.get(col)
        conc = compute_well_concentration(stez_s1, wsp_rozc, i)
        if mode == "wizualny":
            growth = classify_wizualny_well(raw_value)
            if growth is None and pd.notna(raw_value):
                parse_notes.append(f"{col}: nierozpoznana wartość {raw_value!r} - pominięto.")
        else:
            growth = classify_od_well(raw_value, ctrl[COL_GROWTH_CONTROL], ctrl[COL_STERILITY_CONTROL])
        wells_detail.append({"well": col, "conc": conc, "raw": raw_value, "growth": growth})
        if growth is not None:
            conc_growth_pairs.append((conc, growth))

    mic_result = compute_mic_from_wells(conc_growth_pairs)
    reason = mic_result["reason"]
    if parse_notes:
        reason = reason + " | " + "; ".join(parse_notes)

    return _result(
        base, mic_result["mic_value"], mic_result["status"], reason,
        mic_result["censored"], wells_detail,
    )


def process_mic_wizualny(df, controls_df):
    """Przetwarza KAŻDY wiersz arkusza MIC_wizualny. Zwraca listę dictów wyniku."""
    bact_col = utils.find_bacteria_column(df)
    well_cols = [c for c in WELL_COLUMNS if c in df.columns]
    return [_process_row(row, bact_col, well_cols, controls_df, "wizualny") for _, row in df.iterrows()]


def process_mic_od(df, controls_df):
    """Przetwarza KAŻDY wiersz arkusza MIC_OD. Zwraca listę dictów wyniku."""
    bact_col = utils.find_bacteria_column(df)
    well_cols = [c for c in WELL_COLUMNS if c in df.columns]
    return [_process_row(row, bact_col, well_cols, controls_df, "od") for _, row in df.iterrows()]
