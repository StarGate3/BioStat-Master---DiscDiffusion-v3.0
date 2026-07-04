"""
Silnik odczytu MIC ze studzienek (Faza 2 modułu MIC/MBC).

Zakres tej fazy: dla KAŻDEGO pojedynczego wiersza (powtórzenia) arkusza
MIC_wizualny lub MIC_OD, wyznacz wartość MIC z serii studzienek S1..S10.
Nie ma tu jeszcze agregacji powtórzeń, statystyki międzygrupowej ani
wykresów - to kolejne fazy.

Podstawowy algorytm (compute_endpoint_from_wells) jest w pełni generyczny -
wspólny dla obu trybów odczytu MIC i (od Fazy MBC) dla MBC; różni je tylko
sposób sprowadzenia surowej wartości studzienki do booleana "kryterium
niespełnione" (patrz classify_wizualny_well / classify_od_well / classify_mbc_well).
"""
import math
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
import scikit_posthocs as sp
import utils
from config import (
    MIC_STATUS_OK, MIC_STATUS_NEEDS_REVIEW, MIC_STATUS_CENSORED_LOW, MIC_STATUS_CENSORED_HIGH,
    MIC_STATUS_INVALID_RUN, MIC_STATUS_MISSING_CONTROLS, MIC_STATUS_NO_DATA, MIC_LOW_N_BIO_WARNING,
    COL_RUN, COL_STEZ_S1, COL_DILUTION_FACTOR, WELL_COLUMNS,
    WELL_STATUS_GROWTH, WELL_STATUS_NO_GROWTH, MIC_OD_GROWTH_THRESHOLD, MIC_OD_MIN_GROWTH_SIGNAL,
    MIC_OD_STERILITY_MAX,
    COL_SUBSTANCE, COL_TYPE, COL_UNIT, COL_REP_BIO, COL_REP_TECH,
    COL_GROWTH_CONTROL, COL_STERILITY_CONTROL, COL_INOCULUM,
    MIC_MEANINGFUL_DILUTION_DIFF, MIC_MAX_CENSORED_FRACTION, MIN_N_BIO_FOR_COMPARISON,
    ALPHA, MBC_REDUCTION_THRESHOLD, MBC_MIC_BACTERICIDAL_MAX_D, MBC_MIC_BACTERIOSTATIC_MIN_D,
    MBC_MIC_CLASSIFICATION_DILUTION_FACTOR,
)


def compute_well_concentration(stez_s1, wsp_rozc, well_number):
    """
    Stężenie studzienki S{well_number} (well_number=1 dla S1, najwyższe
    stężenie, malejąco dla kolejnych numerów): Stez_S1 / (Wsp_rozc ** (n-1)).
    """
    return stez_s1 / (wsp_rozc ** (well_number - 1))


def compute_endpoint_from_wells(conc_growth_pairs, endpoint_name="MIC"):
    """
    Rdzeń algorytmu - w pełni GENERYCZNY: niezależny od trybu odczytu
    (wizualny/OD/CFU), od arkusza, i od tego, czy wyznacza się MIC czy MBC
    (endpoint_name steruje tylko treścią komunikatów). Reużywane bez
    duplikacji przez _process_row (MIC) i _process_mbc_row (MBC).

    conc_growth_pairs: lista (stężenie: float, kryterium_niespełnione: bool)
    TYLKO dla studzienek z rozpoznaną wartością (studzienki puste/
    nierozpoznane mają być odfiltrowane PRZED wywołaniem). Kolejność
    wejściowa dowolna - funkcja sama sortuje malejąco po stężeniu (S1 =
    najwyższe stężenie).

    Drugi element pary (nazwany tu ogólnie "kryterium niespełnione", bo
    znaczenie zależy od kontekstu):
      - MIC: True = studzienka urosła ("wzrost", NIE zahamowana),
             False = "brak" wzrostu (zahamowana - kryterium inhibicji spełnione).
      - MBC: True = studzienka NIE osiągnęła progu redukcji CFU (bakterie
             przeżyły w istotnej liczbie - kryterium bójcze NIE spełnione),
             False = próg redukcji osiągnięty (kryterium bójcze spełnione).
    W obu przypadkach oczekiwany wzorzec jest identyczny: kryterium
    spełnione (False) przy wysokich stężeniach, przechodzące w niespełnione
    (True) przy niskich - stąd wspólny algorytm.

    Zwraca dict:
      - mic_value (float | None): wartość liczbowa (MIC albo MBC), albo
        granica przy cenzurze (patrz "censored"), albo None gdy nie da się
        wyznaczyć. Nazwa klucza zachowana jako "mic_value" dla obu
        zastosowań (patrz uzasadnienie w module mic_logic - to jeden,
        współdzielony kształt wyniku używany dalej przez Fazy 3/4).
      - status: jeden z MIC_STATUS_* (config.py) - te same stałe dla MIC i MBC.
      - reason (str): zawsze czytelny, nawet dla status="ok".
      - censored: "low" | "high" | None. Gdy nie None, mic_value to GRANICA,
        NIE dokładna wartość - wyświetlaj zawsze z prefiksem ≤/> (patrz
        format_mic_display), nigdy jako zwykłą liczbę.

    Zasady:
      1. Kryterium spełnione we WSZYSTKICH studzienkach -> wartość <=
         najniższe testowane stężenie (censored="low").
      2. Kryterium NIE spełnione w ŻADNEJ studzience -> wartość > najwyższe
         testowane stężenie (censored="high").
      3. Kryterium niespełnione widoczne już w najwyższym testowanym
         stężeniu, a mimo to jakieś niższe stężenie JE spełnia -> wynik
         niejednoznaczny, status wymaga_weryfikacji, mic_value=None (nie ma
         wiarygodnej wartości do podania).
      4. W pozostałych przypadkach: wartość = stężenie ostatniej studzienki
         SPEŁNIAJĄCEJ kryterium, WYSTĘPUJĄCEJ PRZED pierwszym napotkanym
         niespełnieniem (skanując od najwyższego stężenia w dół). To jest
         CELOWO ta sama formuła niezależnie od tego, czy seria jest potem
         monotoniczna - wybrano ją jako regułę zachowawczą: nigdy nie
         zgłasza wartości niższej niż stężenie, przy którym już raz
         zaobserwowano niespełnienie kryterium. Gdy po tym pierwszym
         niespełnieniu pojawi się jednak jeszcze jakieś spełnienie
         (naruszenie monotoniczności - możliwy "skip well"), wartość
         zostaje taka sama, ale status zmienia się na wymaga_weryfikacji
         z opisem, które studzienki są anomalią.
    """
    if not conc_growth_pairs:
        return {
            "mic_value": None,
            "status": MIC_STATUS_NEEDS_REVIEW,
            "reason": f"Brak odczytanych (rozpoznanych) studzienek - nie można wyznaczyć {endpoint_name}.",
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
                f"Kryterium spełnione we wszystkich testowanych studzienkach - prawdziwe {endpoint_name} "
                f"może być niższe niż najniższe testowane stężenie ({lowest:g})."
            ),
            "censored": "low",
        }

    if all(growths):
        highest = max(concs)
        return {
            "mic_value": highest,
            "status": MIC_STATUS_CENSORED_HIGH,
            "reason": (
                f"Kryterium NIE spełnione w żadnej testowanej studzience - prawdziwe {endpoint_name} "
                f"może być wyższe niż najwyższe testowane stężenie ({highest:g})."
            ),
            "censored": "high",
        }

    first_growth_idx = growths.index(True)

    if first_growth_idx == 0:
        return {
            "mic_value": None,
            "status": MIC_STATUS_NEEDS_REVIEW,
            "reason": (
                f"Kryterium NIE spełnione już przy najwyższym testowanym stężeniu, mimo że niższe "
                f"stężenie(a) je spełniły - wynik niejednoznaczny, nie można wyznaczyć {endpoint_name} "
                f"standardową metodą (brak studzienki 'powyżej' punktu odcięcia)."
            ),
            "censored": None,
        }

    mic_value = concs[first_growth_idx - 1]
    is_monotonic = all(growths[first_growth_idx:])

    if is_monotonic:
        return {
            "mic_value": mic_value,
            "status": MIC_STATUS_OK,
            "reason": f"Seria monotoniczna - {endpoint_name} odczytane jednoznacznie.",
            "censored": None,
        }

    skip_concs = [concs[i] for i in range(first_growth_idx, len(growths)) if not growths[i]]
    skip_desc = ", ".join(f"{c:g}" for c in skip_concs)
    reason = (
        f"Niemonotoniczna seria: kryterium przestało być spełnione już przy stężeniu "
        f"{concs[first_growth_idx]:g}, ale zostało spełnione ponownie przy niższym "
        f"stężeniu/stężeniach ({skip_desc}) - możliwy błąd studzienki (skip well) lub skażenie. "
        f"{endpoint_name} odczytane ZACHOWAWCZO jako ostatnia studzienka spełniająca kryterium PRZED "
        f"pierwszym napotkanym niespełnieniem ({mic_value:g}), a NIE jako najniższe spełnienie w całej "
        f"serii - wymaga ręcznej weryfikacji."
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
    Zwraca (ctrl, ambiguous_reason):
      - (dict, None): dokładnie jedno dopasowanie w arkuszu Kontrole -
        {Kontrola_wzrostu, Kontrola_jalowosci, Inokulum_CFU_t0} (wartości
        SUROWE - Kontrola_wzrostu/jalowosci mogą być liczbą OD albo tekstem
        "wzrost"/"brak", patrz validate_run; Inokulum_CFU_t0 używane tylko
        przez MBC).
      - (None, None): brak dopasowania - controls_df jest None, Przebieg
        jest pusty, albo ta wartość Przebiegu w ogóle nie występuje w
        arkuszu Kontrole.
      - (None, str) (audyt 1.8): NIEJEDNOZNACZNE - arkusz Kontrole zawiera
        WIĘCEJ NIŻ JEDEN wiersz dla tego samego Przebiegu. Zamiast po cichu
        brać pierwszy napotkany wiersz (poprzednie zachowanie - który wiersz
        jest "pierwszy" zależy od przypadkowej kolejności w pliku, więc
        wybór byłby niedeterministyczny w duchu, nawet jeśli technicznie
        powtarzalny), traktujemy to jak BRAK użytecznej kontroli, ale ze
        zwróconym, jawnym powodem - spójnie z filozofią "nigdy nie zgaduj"
        stosowaną gdzie indziej w programie (np. REF_PLACEHOLDER wymusza
        ręczny wybór grupy referencyjnej zamiast zgadywać, gdy jest
        niejednoznaczna).

    Dopasowanie jest PO SAMEJ WARTOŚCI `przebieg`, GLOBALNIE względem
    całego pliku - nie rozróżnia, z którego arkusza (MIC_wizualny/MIC_OD/
    MBC_posiew) pochodzi wywołujący wiersz. To zamierzone - patrz pełne
    uzasadnienie i wymagany zakres unikalności w docstringu config.COL_RUN
    (audyt 1.7).
    """
    if controls_df is None or pd.isna(przebieg):
        return None, None
    matches = controls_df[controls_df[COL_RUN] == przebieg]
    if matches.empty:
        return None, None
    if len(matches) > 1:
        return None, (
            f"Arkusz Kontrole zawiera {len(matches)} wpisy dla Przebiegu {przebieg!r} - "
            f"niejednoznaczne, żaden z nich nie został użyty. Popraw arkusz Kontrole, tak by "
            f"każdy Przebieg występował w nim dokładnie raz."
        )
    row = matches.iloc[0]
    return {
        COL_GROWTH_CONTROL: row.get(COL_GROWTH_CONTROL),
        COL_STERILITY_CONTROL: row.get(COL_STERILITY_CONTROL),
        COL_INOCULUM: row.get(COL_INOCULUM),
    }, None


def _try_parse_numeric_control(value):
    """Zwraca float(value) jeśli to rzeczywiście liczba (OD), inaczej None
    (obejmuje NaN/puste i tekst typu "wzrost"/"brak")."""
    if pd.isna(value) or isinstance(value, str):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _try_parse_well_param(value):
    """
    Bezpiecznie sprowadza Stez_S1/Wsp_rozc do float. Zwraca None, gdy
    wartość jest pusta/NaN albo nie da się jej przekształcić na liczbę
    (np. literówka "5o") - NIGDY nie rzuca wyjątku, w przeciwieństwie do
    bezpośredniego użycia surowej wartości w dzieleniu/porównaniu (naprawa
    audytu 1.2 - jeden zły wiersz nie może wywalać całego arkusza). Liczbowy
    tekst (np. "1024" jako komórka tekstowa) jest akceptowany, tak samo jak
    utils.validate_excel_data robi to dla Srednica_mm w module dyfuzji.
    """
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _excel_row_number(row):
    """
    Numer wiersza W EXCELU (nie w pandas) do komunikatów o odrzuceniu wiersza
    - TA SAMA konwencja co utils.validate_excel_data: indeks pandas + 2
    (zero-indexing + wiersz nagłówka), żeby wskazanie "który wiersz" było
    spójne w całej aplikacji.
    """
    return row.name + 2


def validate_run(kontrola_wzrostu_raw, kontrola_jalowosci_raw):
    """
    Sprawdza WAŻNOŚĆ przebiegu DWOMA NIEZALEŻNYMI warunkami - żaden nie może
    zamaskować drugiego (np. bardzo wysoka kontrola wzrostu nie "odkupuje"
    skażonej kontroli jałowości):

      a) Kontrola wzrostu musi realnie urosnąć:
         liczbowo: (Kontrola_wzrostu - Kontrola_jalowosci) >= MIC_OD_MIN_GROWTH_SIGNAL
         tekstowo: wartość != "brak" (patrz WELL_STATUS_NO_GROWTH)
      b) Kontrola jałowości musi pozostać czysta:
         liczbowo: Kontrola_jalowosci <= MIC_OD_STERILITY_MAX (próg BEZWZGLĘDNY,
                   niezależny od poziomu kontroli wzrostu)
         tekstowo: wartość != "wzrost" (patrz WELL_STATUS_GROWTH)

    Działa identycznie niezależnie od trybu odczytu MIC (wizualny/OD) - liczy
    się tylko to, czy KONKRETNA kontrola dla tego Przebiegu jest podana
    liczbowo czy tekstowo, ocenianie każdej z dwóch kontroli jest niezależne.

    Gdy dana kontrola nie jest ani liczbą, ani rozpoznanym tekstem (brak
    wpisu/nierozpoznana wartość), TEN warunek jest pomijany (nie da się
    zweryfikować - nie blokujemy z tego powodu, zgodnie z dotychczasowym
    zachowaniem dla brakujących kontroli).

    Zwraca (is_valid: bool, reason: str | None); reason jasno wskazuje, który
    warunek (albo oba) zawiódł.
    """
    problems = []

    kw_num = _try_parse_numeric_control(kontrola_wzrostu_raw)
    kj_num = _try_parse_numeric_control(kontrola_jalowosci_raw)

    # --- a) Kontrola wzrostu ---
    if kw_num is not None and kj_num is not None:
        signal = kw_num - kj_num
        if signal < MIC_OD_MIN_GROWTH_SIGNAL:
            problems.append(
                f"kontrola wzrostu nie urosła wystarczająco (sygnał po odjęciu tła = {signal:g} OD, "
                f"wymagane >= {MIC_OD_MIN_GROWTH_SIGNAL:g} OD)"
            )
    else:
        kw_bool = classify_wizualny_well(kontrola_wzrostu_raw)
        if kw_bool is False:
            problems.append("kontrola wzrostu nie urosła (odczyt wizualny: 'brak')")
        # kw_bool is True ("wzrost") -> warunek spełniony; None -> nie da się zweryfikować, pomijamy.

    # --- b) Kontrola jałowości (próg BEZWZGLĘDNY, niezależny od a) ---
    if kj_num is not None:
        if kj_num > MIC_OD_STERILITY_MAX:
            problems.append(
                f"kontrola jałowości wskazuje skażenie podłoża ({kj_num:g} OD, dopuszczalne "
                f"<= {MIC_OD_STERILITY_MAX:g} OD)"
            )
    else:
        kj_bool = classify_wizualny_well(kontrola_jalowosci_raw)
        if kj_bool is True:
            problems.append("kontrola jałowości wskazuje skażenie podłoża (odczyt wizualny: 'wzrost')")
        # kj_bool is False ("brak") -> warunek spełniony; None -> nie da się zweryfikować, pomijamy.

    if problems:
        return False, "Przebieg nieważny: " + " ORAZ ".join(problems) + "."
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
        # Wsp_rozc: None tutaj (jeszcze niesparsowane w tym miejscu wywołania) -
        # _process_row/_process_mbc_row nadpisuje wartością SPARSOWANĄ (float)
        # zaraz po udanej walidacji, żeby compute_mbc_mic_ratio (audyt 1.4)
        # wiedziało, czy seria była dwukrotna (jedyny współczynnik, dla którego
        # klasyfikacja bakteriobójcze/bakteriostatyczne jest zdefiniowana).
        "Wsp_rozc": None,
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
    ctrl, ctrl_ambiguous_reason = lookup_controls(controls_df, przebieg)
    if ctrl_ambiguous_reason:
        return _result(base, None, MIC_STATUS_MISSING_CONTROLS, ctrl_ambiguous_reason, None, [])
    kw_num = _try_parse_numeric_control(ctrl.get(COL_GROWTH_CONTROL)) if ctrl is not None else None
    kj_num = _try_parse_numeric_control(ctrl.get(COL_STERILITY_CONTROL)) if ctrl is not None else None
    has_numeric_ctrl = kw_num is not None and kj_num is not None

    if mode == "od" and not has_numeric_ctrl:
        return _result(
            base, None, MIC_STATUS_MISSING_CONTROLS,
            f"Brak liczbowych kontroli (Kontrola_wzrostu/Kontrola_jalowosci) w arkuszu Kontrole dla "
            f"Przebiegu '{przebieg}' - wymagane do wyznaczenia tła i progu wzrostu OD.",
            None, [],
        )

    # Ważność przebiegu jest egzekwowana NIEZALEŻNIE OD TRYBU, zawsze gdy dla
    # tego Przebiegu w ogóle podano jakiekolwiek kontrole (liczbowo lub
    # tekstowo) - patrz validate_run. Gdy kontrole w ogóle nie zostały
    # podane (ctrl is None), nie blokujemy - jak dotychczas.
    if ctrl is not None:
        is_valid, invalid_reason = validate_run(ctrl.get(COL_GROWTH_CONTROL), ctrl.get(COL_STERILITY_CONTROL))
        if not is_valid:
            return _result(base, None, MIC_STATUS_INVALID_RUN, invalid_reason, None, [])

    stez_s1_raw = row.get(COL_STEZ_S1)
    wsp_rozc_raw = row.get(COL_DILUTION_FACTOR)
    stez_s1 = _try_parse_well_param(stez_s1_raw)
    wsp_rozc = _try_parse_well_param(wsp_rozc_raw)
    excel_row = _excel_row_number(row)

    if stez_s1 is None:
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"Wiersz {excel_row}: nieprawidłowa wartość {COL_STEZ_S1!r} ({stez_s1_raw!r}) - nie można "
            f"przekształcić na liczbę. Wiersz pominięty.",
            None, [],
        )
    if wsp_rozc is None:
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"Wiersz {excel_row}: nieprawidłowa wartość {COL_DILUTION_FACTOR!r} ({wsp_rozc_raw!r}) - nie "
            f"można przekształcić na liczbę. Wiersz pominięty.",
            None, [],
        )
    if wsp_rozc <= 1:
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"Wiersz {excel_row}: {COL_DILUTION_FACTOR!r}={wsp_rozc:g} musi być > 1 (malejący szereg "
            f"rozcieńczeń). Wiersz pominięty.",
            None, [],
        )
    base["Wsp_rozc"] = wsp_rozc

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
            growth = classify_od_well(raw_value, kw_num, kj_num)
        wells_detail.append({"well": col, "conc": conc, "raw": raw_value, "growth": growth})
        if growth is not None:
            conc_growth_pairs.append((conc, growth))

    mic_result = compute_endpoint_from_wells(conc_growth_pairs, endpoint_name="MIC")
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


# ============================================================
# SILNIK ODCZYTU MBC (Faza MBC, arkusz MBC_posiew)
# ============================================================
#
# Reużywa maksymalnie maszynerię MIC: compute_well_concentration,
# compute_endpoint_from_wells (rdzeń algorytmu, bez zmian), lookup_controls
# i validate_run (ważność przebiegu - dokładnie ta sama logika co dla MIC).
# Jedyny naprawdę NOWY element to classify_mbc_well (redukcja CFU zamiast
# wzrost/brak albo progu OD) i wymóg dodatkowej danej: CFU_t0 z arkusza
# Kontrole.
#
# ZAŁOŻENIE (odnotowane tu jawnie, patrz też komunikat w wyniku): CFU
# studzienki i CFU_t0 MUSZĄ być na tej samej podstawie objętościowej (to
# samo rozcieńczenie/objętość posiewu) - redukcja jest liczona jako prosty
# stosunek 1 - CFU_studzienki/CFU_t0, funkcja nie ma jak zweryfikować tej
# spójności, dysponując wyłącznie liczbami.

def classify_mbc_well(cfu_value, cfu_t0, threshold=MBC_REDUCTION_THRESHOLD):
    """
    Zwraca True gdy studzienka NIE spełnia kryterium bójczego (redukcja <
    próg - odpowiednik "wzrostu" dla wspólnego algorytmu), False gdy
    kryterium jest spełnione (redukcja >= próg - odpowiednik "brak"), albo
    None gdy cfu_value/cfu_t0 nie da się odczytać jako liczby.

    redukcja = 1 - (CFU_studzienki / CFU_t0). Patrz założenie o wspólnej
    podstawie objętościowej w nagłówku sekcji.
    """
    if pd.isna(cfu_value) or cfu_t0 is None or pd.isna(cfu_t0) or cfu_t0 <= 0:
        return None
    try:
        cfu_value = float(cfu_value)
    except (TypeError, ValueError):
        return None
    redukcja = 1.0 - (cfu_value / cfu_t0)
    return redukcja < threshold


def _process_mbc_row(row, bact_col, well_cols, controls_df):
    """
    Odpowiednik _process_row dla MBC_posiew. Reużywa lookup_controls i
    validate_run identycznie jak MIC (ta sama ważność przebiegu), z jednym
    dodatkowym wymogiem specyficznym dla MBC: liczbowe CFU_t0 w arkuszu
    Kontrole (bez tego nie da się policzyć redukcji w ogóle).
    """
    base = _base_fields(row, bact_col)
    przebieg = row.get(COL_RUN)
    ctrl, ctrl_ambiguous_reason = lookup_controls(controls_df, przebieg)
    if ctrl_ambiguous_reason:
        return _result(base, None, MIC_STATUS_MISSING_CONTROLS, ctrl_ambiguous_reason, None, [])

    if ctrl is not None:
        is_valid, invalid_reason = validate_run(ctrl.get(COL_GROWTH_CONTROL), ctrl.get(COL_STERILITY_CONTROL))
        if not is_valid:
            return _result(base, None, MIC_STATUS_INVALID_RUN, invalid_reason, None, [])

    cfu_t0_raw = ctrl.get(COL_INOCULUM) if ctrl is not None else None
    cfu_t0 = _try_parse_numeric_control(cfu_t0_raw)
    if cfu_t0 is None or cfu_t0 <= 0:
        return _result(
            base, None, MIC_STATUS_MISSING_CONTROLS,
            f"Brak liczbowej wartości {COL_INOCULUM!r} (CFU t0) w arkuszu Kontrole dla Przebiegu "
            f"'{przebieg}' - wymagana do policzenia redukcji CFU. Nie liczę MBC dla tego wiersza.",
            None, [],
        )

    stez_s1_raw = row.get(COL_STEZ_S1)
    wsp_rozc_raw = row.get(COL_DILUTION_FACTOR)
    stez_s1 = _try_parse_well_param(stez_s1_raw)
    wsp_rozc = _try_parse_well_param(wsp_rozc_raw)
    excel_row = _excel_row_number(row)

    if stez_s1 is None:
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"Wiersz {excel_row}: nieprawidłowa wartość {COL_STEZ_S1!r} ({stez_s1_raw!r}) - nie można "
            f"przekształcić na liczbę. Wiersz pominięty.",
            None, [],
        )
    if wsp_rozc is None:
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"Wiersz {excel_row}: nieprawidłowa wartość {COL_DILUTION_FACTOR!r} ({wsp_rozc_raw!r}) - nie "
            f"można przekształcić na liczbę. Wiersz pominięty.",
            None, [],
        )
    if wsp_rozc <= 1:
        return _result(
            base, None, MIC_STATUS_NEEDS_REVIEW,
            f"Wiersz {excel_row}: {COL_DILUTION_FACTOR!r}={wsp_rozc:g} musi być > 1 (malejący szereg "
            f"rozcieńczeń). Wiersz pominięty.",
            None, [],
        )
    base["Wsp_rozc"] = wsp_rozc

    conc_growth_pairs = []
    wells_detail = []
    parse_notes = []
    for i, col in enumerate(well_cols, start=1):
        raw_value = row.get(col)
        conc = compute_well_concentration(stez_s1, wsp_rozc, i)
        does_not_meet = classify_mbc_well(raw_value, cfu_t0)
        if does_not_meet is None and pd.notna(raw_value):
            parse_notes.append(f"{col}: nierozpoznana wartość CFU {raw_value!r} - pominięto.")
        wells_detail.append({"well": col, "conc": conc, "raw": raw_value, "growth": does_not_meet})
        if does_not_meet is not None:
            conc_growth_pairs.append((conc, does_not_meet))

    mbc_result = compute_endpoint_from_wells(conc_growth_pairs, endpoint_name="MBC")
    reason = (
        mbc_result["reason"] + " [Założenie: CFU studzienki i CFU_t0 są na tej samej podstawie "
        "objętościowej (to samo rozcieńczenie/objętość posiewu) - redukcja liczona jako prosty stosunek.]"
    )
    if parse_notes:
        reason = reason + " | " + "; ".join(parse_notes)

    return _result(
        base, mbc_result["mic_value"], mbc_result["status"], reason,
        mbc_result["censored"], wells_detail,
    )


def process_mbc(df, controls_df):
    """Przetwarza KAŻDY wiersz arkusza MBC_posiew. Zwraca listę dictów wyniku
    (dokładnie ten sam kształt co process_mic_wizualny/process_mic_od)."""
    bact_col = utils.find_bacteria_column(df)
    well_cols = [c for c in WELL_COLUMNS if c in df.columns]
    return [_process_mbc_row(row, bact_col, well_cols, controls_df) for _, row in df.iterrows()]


# ============================================================
# AGREGACJA POWTÓRZEŃ (Faza 3): techniczne -> biologiczne -> grupa
# ============================================================
#
# WAŻNE: wszystkie operacje poniżej działają na SKALI ROZCIEŃCZEŃ (log2
# stężenia), nigdy na surowym stężeniu jako liczbie ciągłej. log2(stężenie)
# jest niezależny od tego, jakie Stez_S1/Wsp_rozc miało konkretne
# powtórzenie - dzięki temu powtórzenia o różnych schematach rozcieńczeń
# wciąż da się sensownie agregować na wspólnej skali.
#
# REGUŁA MEDIANY ("high median"): dla N wartości posortowanych rosnąco na
# skali log2, wynikiem jest element o indeksie N//2 (indeksowanie od 0).
# Dla N nieparzystego to zwykły środkowy element. Dla N parzystego to
# WYŻSZY z dwóch środkowych elementów - CELOWO nie uśredniamy dwóch
# środkowych (interpolacja dałaby stężenie, które nigdy nie było testowane,
# czyli "spoza skali rozcieńczeń"). Wybieramy zawsze ISTNIEJĄCY wpis, nigdy
# nowo obliczoną liczbę - to jest też MECHANIZM propagacji cenzury: skoro
# wynik to zawsze jeden z rzeczywistych wpisów wejściowych, jego flaga
# censored automatycznie staje się flagą wyniku, bez żadnej dodatkowej
# logiki "co zrobić, gdy mediana wypada w cenzurze".

def _censor_rank(censored):
    """
    Porządek pomocniczy do sortowania przy remisach na tej samej wartości
    log2: '<=X' sortuje TUŻ PRZED dokładnym '=X' (bo prawdziwa wartość może
    być jeszcze niższa), a '>X' TUŻ PO (bo prawdziwa wartość może być
    jeszcze wyższa). Rzadki przypadek (dokładny remis liczbowy), ale
    zapewnia jednoznaczny, przewidywalny porządek zamiast losowego.
    """
    return {"low": -1, None: 0, "high": 1}[censored]


def _high_median_entry(entries):
    """
    entries: niepusta lista dictów zawierających co najmniej klucze
    'log2' (float) i 'censored' (None|'low'|'high'). Dowolne dodatkowe
    klucze są zachowywane bez zmian.

    Zwraca WYBRANY wpis (ten sam obiekt co na wejściu, nie nową wartość)
    wg reguły "high median" opisanej w nagłówku sekcji.
    """
    ordered = sorted(entries, key=lambda e: (e["log2"], _censor_rank(e["censored"])))
    return ordered[len(ordered) // 2]


def _distinct_units(rows):
    """Zbiera niepuste, przycięte, unikalne wartości Jednostka z listy
    wierszy (posortowane, dla przewidywalnego wyświetlania) - budulec do
    wykrywania niespójności jednostek (audyt 1.5, patrz summarize_mic_group)."""
    units = set()
    for r in rows:
        u = r.get("Jednostka")
        if u is None:
            continue
        if isinstance(u, float) and math.isnan(u):
            continue
        u_txt = str(u).strip()
        if u_txt:
            units.add(u_txt)
    return sorted(units)


def aggregate_technical_to_biological(row_results):
    """
    Redukuje powtórzenia TECHNICZNE (wiersze z process_mic_wizualny/
    process_mic_od dla JEDNEJ kombinacji Bakteria+Substancja+Rep_biologiczna)
    do JEDNEJ wartości MIC tego powtórzenia biologicznego.

    Kroki:
      1. Wykluczamy wiersze bez wyznaczonej wartości (mic_value=None -
         obejmuje status nieważny/brak_kontroli, oraz brzegowy przypadek
         wymaga_weryfikacji "wzrost już w S1"). Zliczamy i opisujemy powody.
      2. Jeśli po wykluczeniu nic nie zostaje: wynik = MIC_STATUS_NO_DATA,
         z jawnym powodem wymieniającym wykluczone powtórzenia.
      3. W przeciwnym razie: konwertujemy pozostałe mic_value na log2 i
         wybieramy _high_median_entry - wynik (wartość + flaga cenzury)
         jest DOKŁADNIE jednym z wejściowych powtórzeń technicznych.

    Wiersze o statusie wymaga_weryfikacji NADAL uczestniczą w agregacji
    (mają realną, tylko zachowawczo wyznaczoną wartość) - liczymy je w
    n_flagged, żeby niepewność była widoczna w wyniku, a nie ukryta.

    "unit" to Jednostka WYBRANEGO powtórzenia (pierwszy wiersz - jak
    dotychczas, do samego wyświetlania). "tech_units" (audyt 1.5) to
    NIEZALEŻNIE od tego - zbiór WSZYSTKICH odrębnych wartości Jednostka
    zaobserwowanych wśród powtórzeń technicznych tego jednego powtórzenia
    biologicznego (może mieć >1 element, gdy dane są niespójne - to
    summarize_mic_group przenosi dalej do jawnego ostrzeżenia grupy, patrz
    tam; ta funkcja tylko zbiera fakty, nie ostrzega sama).

    Zwraca dict: Bakteria, Substancja, Rep_biologiczna, mic_value, unit,
    tech_units, status, reason, censored, Wsp_rozc (audyt 1.4 - współczynnik
    rozcieńczenia WYBRANEGO powtórzenia, potrzebny compute_mbc_mic_ratio do
    ustalenia, czy klasyfikacja bakteriobójcze/bakteriostatyczne w ogóle ma
    zastosowanie - patrz tam), n_tech_total, n_tech_used, n_tech_excluded,
    n_flagged. Klucz "rep_bio_fallback_warning" jest dopisywany PRZEZ
    WYWOŁUJĄCEGO (aggregate_all), nie przez tę funkcję - patrz tam.
    """
    if not row_results:
        raise ValueError("aggregate_technical_to_biological: pusta lista wejściowa.")

    first = row_results[0]
    base = {
        "Bakteria": first.get("Bakteria"),
        "Substancja": first.get("Substancja"),
        "Rep_biologiczna": first.get("Rep_biologiczna"),
    }
    unit = first.get("Jednostka")
    tech_units = _distinct_units(row_results)

    usable = []
    excluded_reasons = []
    for r in row_results:
        if r.get("mic_value") is None:
            excluded_reasons.append(
                f"Rep_techniczna={r.get('Rep_techniczna')}: wykluczone ({r['status']}) - {r['reason']}"
            )
        else:
            usable.append(r)

    n_tech_total = len(row_results)
    n_tech_excluded = len(excluded_reasons)

    if not usable:
        reason = (
            f"Brak danych: wszystkie {n_tech_total} powtórzenia(a) techniczne wykluczone. "
            + " | ".join(excluded_reasons)
        )
        return {
            **base, "mic_value": None, "unit": unit, "tech_units": tech_units, "status": MIC_STATUS_NO_DATA,
            "reason": reason, "censored": None, "Wsp_rozc": None,
            "n_tech_total": n_tech_total, "n_tech_used": 0,
            "n_tech_excluded": n_tech_excluded, "n_flagged": 0,
        }

    entries = [
        {"log2": math.log2(r["mic_value"]), "censored": r["censored"], "source": r}
        for r in usable
    ]
    picked = _high_median_entry(entries)
    picked_row = picked["source"]
    n_flagged = sum(1 for r in usable if r["status"] == MIC_STATUS_NEEDS_REVIEW)

    reason = (
        f"Mediana (high-median, n_tech={len(usable)}/{n_tech_total}) powtórzeń technicznych - "
        f"wybrane powtórzenie Rep_techniczna={picked_row.get('Rep_techniczna')} "
        f"(status źródłowy: {picked_row['status']})."
    )
    if excluded_reasons:
        reason += " Wykluczone: " + " | ".join(excluded_reasons)
    if n_flagged:
        reason += f" UWAGA: {n_flagged}/{len(usable)} użytych odczytów miało status wymaga_weryfikacji."

    return {
        **base, "mic_value": picked_row["mic_value"], "unit": unit, "tech_units": tech_units,
        "status": picked_row["status"], "reason": reason, "censored": picked_row["censored"],
        "Wsp_rozc": picked_row.get("Wsp_rozc"),
        "n_tech_total": n_tech_total, "n_tech_used": len(usable),
        "n_tech_excluded": n_tech_excluded, "n_flagged": n_flagged,
    }


def _value_display(mic_value, censored, unit):
    if mic_value is None:
        return {"mic_value": None, "censored": None, "unit": unit, "display": "n/d"}
    return {
        "mic_value": mic_value, "censored": censored, "unit": unit,
        "display": f"{format_mic_display(mic_value, censored)} {unit}".strip(),
    }


def summarize_mic_group(bio_results):
    """
    Podsumowanie OPISOWE jednej grupy (Bakteria x Substancja) na podstawie
    powtórzeń BIOLOGICZNYCH (wyniki aggregate_technical_to_biological, jeden
    na Rep_biologiczna). Nie porównuje grup między sobą - to kolejna faza
    (testy istotności).

    n_bio = liczba powtórzeń biologicznych Z WARTOŚCIĄ (bio-repy o statusie
    MIC_STATUS_NO_DATA są wykluczone i zliczone osobno w n_bio_excluded).

    Zwraca dict:
      - Bakteria, Substancja, n_bio, n_bio_excluded, n_flagged (bio-repy
        zbudowane na podstawie flagowanych odczytów technicznych)
      - median: {mic_value, censored, unit, display} - ta sama reguła
        high-median co przy agregacji technicznej->biologicznej (patrz jej
        docstring), więc cenzura propaguje się tak samo automatycznie.
      - range_min / range_max: {..., display} - skrajne wartości grupy,
        KAŻDA z WŁASNĄ flagą cenzury (np. jeśli najniższa wartość w grupie
        jest akurat "<=X", dolna granica zakresu to "<=X", nie "X").
      - mode: lista {..., display} najczęstszych wartości (>=2 wystąpienia
        po zaokrągleniu log2 do 6 miejsc), albo None gdy wszystkie różne.
      - geo_mean: {mic_value, display} albo None - liczona TYLKO gdy ŻADNA
        użyta wartość nie jest cenzurowana (średnia geometryczna wymaga
        dokładnych liczb, nie da się jej sensownie policzyć z granicy typu
        "przynajmniej"/"co najwyżej"); geo_mean_reason wyjaśnia brak.
        NIGDY nie jest to średnia arytmetyczna - liczona na log2, potem 2**x.
      - low_n_bio_warning (bool), warning (str | None) = MIC_LOW_N_BIO_WARNING
        gdy n_bio<2 - wynik nadal jest liczony i pokazany, tylko oznaczony.
      - unit_warning (str | None) - patrz niżej.

    ZAKRES SPRAWDZANIA SPÓJNOŚCI JEDNOSTEK (audyt 1.5) - uzasadnienie: cała
    mediana/zakres/moda/średnia geometryczna tutaj oraz porównania między
    grupami (compare_mic_groups) liczą się na log2(mic_value) z CAŁEJ puli
    powtórzeń tej (Bakteria, Substancja) na raz - więc niespójność jednostek
    GDZIEKOLWIEK w tej puli (czy to między powtórzeniami technicznymi
    jednego powtórzenia biologicznego - patrz "tech_units" z
    aggregate_technical_to_biological, czy między różnymi powtórzeniami
    biologicznymi tej samej substancji - "unit" per bio-rep tutaj) jednakowo
    unieważnia sens dalszej arytmetyki. Sprawdzanie WĘŻSZEGO zakresu (tylko
    wewnątrz jednego powtórzenia biologicznego) przeoczyłoby przypadek, gdy
    CAŁE powtórzenie nr 1 jest w mg/ml, a CAŁE powtórzenie nr 2 w µg/ml -
    każde z osobna spójne, ale razem tak samo błędne jak zademonstrowany w
    audycie przykład (50 mg/ml + 50 µg/ml). Dlatego sprawdzamy na poziomie
    CAŁEJ GRUPY, łącząc jednostki z obu źródeł. To WYŁĄCZNIE ostrzeżenie -
    nigdy nie blokuje ani nie przelicza jednostek (zbyt ryzykowne bez
    jawnie podanego przelicznika).
    """
    if not bio_results:
        raise ValueError("summarize_mic_group: pusta lista wejściowa.")

    first = bio_results[0]
    bakteria, substancja, unit = first.get("Bakteria"), first.get("Substancja"), first.get("unit")

    usable = [b for b in bio_results if b["mic_value"] is not None]
    n_bio_excluded = len(bio_results) - len(usable)
    n_flagged = sum(b.get("n_flagged", 0) for b in bio_results)
    rep_bio_fallback_warnings = [b["rep_bio_fallback_warning"] for b in bio_results if b.get("rep_bio_fallback_warning")]
    rep_bio_fallback_warning = "; ".join(rep_bio_fallback_warnings) if rep_bio_fallback_warnings else None

    all_units = set()
    for b in bio_results:
        if b.get("unit") is not None and str(b["unit"]).strip():
            all_units.add(str(b["unit"]).strip())
        all_units.update(b.get("tech_units") or [])
    unit_warning = None
    if len(all_units) > 1:
        units_txt = ", ".join(sorted(all_units))
        unit_warning = (
            f"Substancja {substancja!r} ({bakteria}): wykryto niespójne jednostki ({units_txt}). "
            f"Wyniki mogą być nieporównywalne - ujednolić jednostki."
        )

    if not usable:
        return {
            "Bakteria": bakteria, "Substancja": substancja,
            "n_bio": 0, "n_bio_excluded": n_bio_excluded, "n_flagged": n_flagged,
            "median": _value_display(None, None, unit), "range_min": _value_display(None, None, unit),
            "range_max": _value_display(None, None, unit), "mode": None,
            "geo_mean": None, "geo_mean_reason": "Brak danych.",
            "low_n_bio_warning": True, "warning": MIC_LOW_N_BIO_WARNING,
            "rep_bio_fallback_warning": rep_bio_fallback_warning,
            "unit_warning": unit_warning,
        }

    entries = [{"log2": math.log2(b["mic_value"]), "censored": b["censored"], "source": b} for b in usable]
    n_bio = len(usable)

    median_entry = _high_median_entry(entries)
    min_entry = min(entries, key=lambda e: e["log2"])
    max_entry = max(entries, key=lambda e: e["log2"])

    rounded = [round(e["log2"], 6) for e in entries]
    counts = {}
    for v in rounded:
        counts[v] = counts.get(v, 0) + 1
    max_count = max(counts.values())
    if max_count >= 2:
        mode_keys = [v for v, c in counts.items() if c == max_count]
        mode = [
            _value_display(m["source"]["mic_value"], m["censored"], unit)
            for m in (next(e for e in entries if round(e["log2"], 6) == k) for k in mode_keys)
        ]
    else:
        mode = None

    any_censored = any(e["censored"] is not None for e in entries)
    if any_censored:
        geo_mean = None
        geo_mean_reason = (
            "Nie liczono - grupa zawiera wartości cenzurowane (średnia geometryczna wymaga "
            "wyłącznie dokładnych wartości, nie da się jej sensownie wyznaczyć z granicy typu "
            "'≤'/'>'). Użyj mediany jako wartości centralnej."
        )
    else:
        geo_log2_mean = sum(e["log2"] for e in entries) / len(entries)
        geo_mean_value = 2 ** geo_log2_mean
        geo_mean = {
            "mic_value": geo_mean_value, "unit": unit,
            "display": f"{format_mic_display(geo_mean_value, None)} {unit} (średnia geometryczna)".strip(),
        }
        geo_mean_reason = None

    low_n_bio_warning = n_bio < 2

    return {
        "Bakteria": bakteria, "Substancja": substancja,
        "n_bio": n_bio, "n_bio_excluded": n_bio_excluded, "n_flagged": n_flagged,
        "median": _value_display(median_entry["source"]["mic_value"], median_entry["censored"], unit),
        "range_min": _value_display(min_entry["source"]["mic_value"], min_entry["censored"], unit),
        "range_max": _value_display(max_entry["source"]["mic_value"], max_entry["censored"], unit),
        "mode": mode,
        "geo_mean": geo_mean, "geo_mean_reason": geo_mean_reason,
        "low_n_bio_warning": low_n_bio_warning,
        "warning": MIC_LOW_N_BIO_WARNING if low_n_bio_warning else None,
        "rep_bio_fallback_warning": rep_bio_fallback_warning,
        "unit_warning": unit_warning,
    }


def _normalize_rep_bio(value):
    """
    Sprowadza surową wartość Rep_biologiczna do porównywalnej postaci przed
    grupowaniem, żeby MIESZANE TYPY tej samej kolumny (typowy artefakt
    Excela - część komórek sformatowana jako liczba, część jako tekst, np.
    int 1 i str '1' dla tego samego zamierzonego powtórzenia) nie tworzyły
    PO CICHU dwóch różnych kluczy grupowania (co po cichu rozbijałoby jedno
    prawdziwe powtórzenie biologiczne na dwa - naprawa 1.3 audytu).

    None/NaN/pusty string -> None (WSPÓLNY klucz "brak deklaracji" - patrz
    aggregate_all: to jest właśnie fallback "cały wiersz = n_bio=1" znany z
    modułu dyfuzji, celowo zachowany, ale teraz głośny gdy scala realnie
    różne wartości, patrz _detect_rep_bio_fallback_warning).
    """
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text = str(value).strip()
    return text if text else None


def _detect_rep_bio_fallback_warning(tech_rows):
    """
    tech_rows: wiersze techniczne, które trafiły do WSPÓLNEGO "fallback"
    powtórzenia biologicznego (Rep_biologiczna brakowała/była pusta dla
    wszystkich - patrz _normalize_rep_bio). Fallback sam w sobie jest
    zamierzony (patrz utils.build_internal_representation - ta sama
    konwencja co stary format dyfuzji: "brak deklaracji = n_bio=1"), ale
    STAJE SIĘ MYLĄCY, gdy scalone wiersze mają WIĘCEJ NIŻ JEDNĄ odrębną
    wartość endpointu (MIC albo MBC) - wtedy realna zmienność (biologiczna
    albo błąd danych) znika bez śladu w raporcie (patrz audyt, przypadek
    512/32/8 mg/ml -> "32-32" bez ostrzeżenia).

    Zwraca czytelne ostrzeżenie z WYPISANYMI wykrytymi wartościami, albo
    None gdy scalona grupa ma 0 lub 1 odrębną wartość (wtedy to zwykłe,
    niegroźne n_bio=1 - baner MIC_LOW_N_BIO_WARNING już to pokrywa, nie ma
    co dublować ostrzeżeniem bez powodu).
    """
    usable = [r for r in tech_rows if r.get("mic_value") is not None]
    if len(usable) < 2:
        return None

    unit = usable[0].get("Jednostka")
    distinct_display = {}
    for r in usable:
        disp = format_mic_display(r["mic_value"], r["censored"])
        distinct_display.setdefault(disp, True)

    if len(distinct_display) <= 1:
        return None

    values_txt = ", ".join(distinct_display.keys())
    unit_txt = f" {unit}" if unit else ""
    return (
        f"Wykryto {len(distinct_display)} różne wartości scalone jako jedno powtórzenie "
        f"biologiczne z powodu braku/niespójnej kolumny {COL_REP_BIO!r}: {values_txt}{unit_txt}. "
        f"Uzupełnij {COL_REP_BIO!r}, aby analiza traktowała je jako osobne powtórzenia."
    )


def aggregate_all(row_results):
    """
    Pełny potok Fazy 3: powtórzenia techniczne -> biologiczne -> podsumowanie
    grupy (Bakteria x Substancja). Przyjmuje listę wyników z
    process_mic_wizualny i/lub process_mic_od (można połączyć obie listy -
    grupowanie jest po Bakteria+Substancja+Rep_biologiczna, więc oba tryby
    współistnieją naturalnie jako kolejne "powtórzenia techniczne", gdyby
    kiedyś tak wystąpiły w danych; rozstrzygnięcie pierwszeństwa
    MIC_wizualny/MIC_OD zostało już zasygnalizowane w Fazie 1).

    Rep_biologiczna jest znormalizowana PRZED grupowaniem (_normalize_rep_bio)
    - naprawia to ciche rozbicie jednego powtórzenia na dwa przy mieszanych
    typach (int 1 vs str '1'). Gdy Rep_biologiczna brakuje/jest pusta dla
    wielu wierszy tej samej (Bakteria, Substancja), trafiają one do JEDNEGO
    wspólnego "fallback" powtórzenia (zamierzone, spójne ze starym formatem
    dyfuzji) - ale jeśli to scalenie ukrywa więcej niż jedną odrębną
    wartość endpointu, bio-rep dostaje jawne pole "rep_bio_fallback_warning"
    (patrz _detect_rep_bio_fallback_warning), które summarize_mic_group
    przenosi dalej do "summary", i stamtąd trafia wszędzie tam, gdzie inne
    ostrzeżenia (wykres, tabela zbiorcza PDF/Excel).

    Podobnie "summary['unit_warning']" sygnalizuje niespójne Jednostka w
    obrębie całej grupy (Bakteria, Substancja) - patrz
    aggregate_technical_to_biological/summarize_mic_group. Też WYŁĄCZNIE
    ostrzeżenie, analiza liczy się dalej.

    NIE porównuje grup między sobą (brak testów istotności) - tylko opisuje
    każdą z osobna.

    Zwraca dict {(Bakteria, Substancja): {"bio_results": [...], "summary": {...}}}.
    """
    by_tech_group = defaultdict(list)
    for r in row_results:
        norm_rep_bio = _normalize_rep_bio(r.get("Rep_biologiczna"))
        key = (r.get("Bakteria"), r.get("Substancja"), norm_rep_bio)
        by_tech_group[key].append(r)

    bio_by_group = defaultdict(list)
    for (bakteria, substancja, norm_rep_bio), tech_rows in by_tech_group.items():
        bio_result = aggregate_technical_to_biological(tech_rows)
        bio_result["rep_bio_fallback_warning"] = (
            _detect_rep_bio_fallback_warning(tech_rows) if norm_rep_bio is None else None
        )
        bio_by_group[(bakteria, substancja)].append(bio_result)

    return {
        key: {"bio_results": bio_results, "summary": summarize_mic_group(bio_results)}
        for key, bio_results in bio_by_group.items()
    }


# ============================================================
# PORÓWNANIA MIĘDZY GRUPAMI (Faza 4) - WARSTWA 1: OPIS
# ============================================================
#
# Warstwa 1 działa ZAWSZE, niezależnie od cenzury i n_bio - to najbardziej
# podstawowy, zawsze-interpretowalny wynik: mediana, zakres, i różnica w
# LICZBIE ROZCIEŃCZEŃ (kroków log2) między medianami porównywanych grup.

def _bio_results_to_entries(bio_results):
    """Konwertuje wyniki biologiczne (Faza 3) na {'log2','censored','source'},
    pomijając te bez wartości (status=brak_danych)."""
    return [
        {"log2": math.log2(b["mic_value"]), "censored": b["censored"], "source": b}
        for b in bio_results if b["mic_value"] is not None
    ]


def describe_mic_group(bio_results):
    """
    Opis Warstwy 1 dla JEDNEJ grupy: n_bio, mediana (high-median, patrz
    Faza 3), zakres min-max (każdy koniec z WŁASNĄ flagą cenzury). To jest
    dokładnie to samo co summarize_mic_group bez mody/geo_mean/ostrzeżenia -
    powtórzone tutaj jako samodzielna, minimalna jednostka do porównań,
    żeby compare_mic_groups nie musiało zależeć od pełnego kształtu
    summarize_mic_group.
    """
    entries = _bio_results_to_entries(bio_results)
    unit = bio_results[0].get("unit") if bio_results else None
    if not entries:
        return {
            "n_bio": 0, "unit": unit,
            "median": _value_display(None, None, unit),
            "range_min": _value_display(None, None, unit),
            "range_max": _value_display(None, None, unit),
            "median_log2": None,
        }
    median_entry = _high_median_entry(entries)
    min_entry = min(entries, key=lambda e: e["log2"])
    max_entry = max(entries, key=lambda e: e["log2"])
    return {
        "n_bio": len(entries), "unit": unit,
        "median": _value_display(median_entry["source"]["mic_value"], median_entry["censored"], unit),
        "range_min": _value_display(min_entry["source"]["mic_value"], min_entry["censored"], unit),
        "range_max": _value_display(max_entry["source"]["mic_value"], max_entry["censored"], unit),
        "median_log2": median_entry["log2"],
    }


def dilution_difference(desc_a, desc_b):
    """
    Różnica median DWÓCH grup (już opisanych przez describe_mic_group) w
    LICZBIE ROZCIEŃCZEŃ (kroki log2) - NIGDY na surowym stężeniu. Dodatnia
    wartość diff_dilutions oznacza, że mediana A > mediana B (A "bardziej
    oporny"/wyższe MIC).

    meaningful=True gdy |diff_dilutions| >= MIC_MEANINGFUL_DILUTION_DIFF
    (domyślnie 2 kroki = różnica 4-krotna - patrz uzasadnienie w config.py).

    Zwraca None gdy którejś z grup brakuje mediany (n_bio=0).
    """
    if desc_a["median_log2"] is None or desc_b["median_log2"] is None:
        return None
    # Zaokrąglenie przed porównaniem z progiem: log2 z podzielenia stężeń
    # które "na papierze" różnią się o dokładną liczbę rozcieńczeń (np.
    # 50/12.5=4=2^2) potrafi dać 1.9999999999999996 zamiast 2.0 przez błąd
    # zmiennoprzecinkowy - bez zaokrąglenia próg ">=" mógłby błędnie odrzucić
    # dokładnie granicznie sensowną różnicę.
    diff = round(desc_a["median_log2"] - desc_b["median_log2"], 9)
    meaningful = abs(diff) >= MIC_MEANINGFUL_DILUTION_DIFF
    if diff > 0:
        direction = "A > B"
    elif diff < 0:
        direction = "A < B"
    else:
        direction = "A = B"
    return {
        "diff_dilutions": diff,
        "meaningful": meaningful,
        "direction": direction,
        "threshold": MIC_MEANINGFUL_DILUTION_DIFF,
    }


# ============================================================
# WARSTWA 2: TEST INFERENCYJNY (Mann-Whitney / Kruskal-Wallis + Dunn)
# ============================================================
#
# REGUŁA CENZURY DLA RANGOWANIA (opisana tu raz, obowiązuje wszędzie w tej
# warstwie): testy rangowe (Mann-Whitney, Kruskal-Wallis, Dunn) liczą rangi
# na POŁĄCZONEJ próbie wszystkich porównywanych grup. Wartość "≤X" jest z
# założenia nie większa niż jakakolwiek wartość dokładna w zbiorze (bo
# dokładne MIC z Fazy 2 zawsze mieszczą się w przetestowanym zakresie
# rozcieńczeń), a ">X" nie mniejsza niż jakakolwiek wartość dokładna -
# dlatego WSZYSTKIE "≤" w połączonej próbie dostają JEDNĄ WSPÓLNĄ,
# NAJNIŻSZĄ wiązaną rangę (są sobie nawzajem "remisowane" niezależnie od
# konkretnego X), a WSZYSTKIE ">" dostają JEDNĄ WSPÓLNĄ, NAJWYŻSZĄ wiązaną
# rangę - wartości dokładne dostają zwykłe rangi (ze standardowym
# remisowaniem average-tie) POMIĘDZY tymi dwoma blokami.
#
# Implementacja: zamiast ręcznie przeliczać wzory testów na gotowych
# rangach (ryzykowne, łatwo o subtelny błąd we wzorze korekty na remisy),
# każdej cenzurowanej wartości przypisujemy WSPÓLNĄ liczbę zastępczą
# (surogat) gwarantowanie mniejszą/większą od każdej wartości dokładnej w
# połączonej próbie, i przekazujemy surogaty do standardowych funkcji
# scipy/scikit-posthocs - ich WŁASNE, przetestowane rangowanie (average-tie)
# automatycznie zwiąże wszystkie cenzurowane-nisko razem na dole i
# cenzurowane-wysoko razem na górze, dokładnie wg powyższej reguły.

def _censoring_surrogate(entries):
    """
    Zwraca listę wartości liczbowych (ta sama długość/kolejność co entries)
    gotową do bezpośredniego podania do scipy.stats.mannwhitneyu/kruskal
    albo scikit_posthocs.posthoc_dunn - patrz reguła cenzury w nagłówku
    sekcji.
    """
    exact = [e["log2"] for e in entries if e["censored"] is None]
    if exact:
        lo, hi = min(exact), max(exact)
    else:
        lo = hi = 0.0
    span = max(hi - lo, 1.0)
    surrogate_low = lo - span - 1.0
    surrogate_high = hi + span + 1.0

    out = []
    for e in entries:
        if e["censored"] == "low":
            out.append(surrogate_low)
        elif e["censored"] == "high":
            out.append(surrogate_high)
        else:
            out.append(e["log2"])
    return out


# ============================================================
# WARSTWA 3: OSŁONY (gaszą p-value, zostaje opis z Warstwy 1)
# ============================================================

def _censored_fraction(entries):
    if not entries:
        return 0.0
    return sum(1 for e in entries if e["censored"] is not None) / len(entries)


def _check_layer3_guards(groups_entries):
    """
    groups_entries: dict {label: lista entries (z _bio_results_to_entries)}.

    Sprawdza warunki Warstwy 3. Zwraca dict:
      - suppressed (bool): czy Warstwa 2 (test) ma być CAŁKOWICIE pominięta.
      - reason (str | None): jawny powód WYGASZENIA (gdy suppressed=True) -
        albo zbyt duża cenzura, albo mniej niż 2 grupy z danymi po
        wykluczeniu pustych.
      - low_n_bio_warning (bool): czy KTÓRAKOLWIEK z UŻYTYCH grup ma n_bio=1.
      - excluded_empty_groups (list[str]): etykiety grup bez ŻADNEJ wartości
        MIC (n_bio=0).

    Audyt 1.6 (spójność z modułem dyfuzji): grupa bez żadnej wartości NIGDY
    nie znika po cichu z wyniku - jest zawsze wymieniona w
    excluded_empty_groups, a compare_mic_groups dokłada z tego jawną notatkę
    do layer2, WIDOCZNĄ niezależnie od tego, czy test w ogóle się wykonał.
    Sama pusta grupa NIE blokuje już całego porównania (jak wcześniej) -
    jest WYKLUCZANA z testu Warstwy 2 (nie ma z niej czego liczyć), a test
    biegnie na POZOSTAŁYCH grupach, jeśli zostaje ich >=2 - dokładnie tak,
    jak moduł dyfuzji (gui.py.run_analysis) teraz też robi: wyklucz + jawnie
    ostrzeż, zamiast całkiem blokować albo cicho pomijać.
    """
    excluded_empty_groups = [label for label, entries in groups_entries.items() if len(entries) == 0]
    remaining_labels = [label for label in groups_entries if label not in excluded_empty_groups]

    if len(remaining_labels) < 2:
        parts = []
        if excluded_empty_groups:
            parts.append(f"grupa(y) bez żadnej wartości MIC (n_bio=0): {', '.join(excluded_empty_groups)}")
        parts.append(
            f"po ich wykluczeniu zostaje mniej niż 2 grupy z danymi "
            f"({', '.join(remaining_labels) if remaining_labels else 'brak'})"
        )
        return {
            "suppressed": True, "excluded_empty_groups": excluded_empty_groups, "low_n_bio_warning": False,
            "reason": "Porównanie (Warstwa 2) niemożliwe: " + "; ".join(parts) + ".",
        }

    low_n_bio_warning = any(len(groups_entries[label]) == 1 for label in remaining_labels)

    # (a)/(b) odsetek cenzury (100% to szczegolny przypadek tego samego testu) -
    # liczony TYLKO na grupach UŻYTYCH w teście (bez pustych, wykluczonych wyżej).
    over_threshold = {
        label: _censored_fraction(groups_entries[label])
        for label in remaining_labels
        if _censored_fraction(groups_entries[label]) > MIC_MAX_CENSORED_FRACTION
    }
    if over_threshold:
        parts = []
        for label, frac in over_threshold.items():
            if frac >= 1.0:
                parts.append(f"'{label}' jest CAŁKOWICIE cenzurowana (100%)")
            else:
                parts.append(f"'{label}' ma {frac:.0%} wartości cenzurowanych")
        reason = (
            f"Test istotności wygaszony: " + "; ".join(parts) +
            f" (próg: >{MIC_MAX_CENSORED_FRACTION:.0%}). Obserwowana różnica byłaby w dużej mierze "
            f"efektem konwencji wiązania rang wartości cenzurowanych, nie realnego sygnału. "
            f"Zostaje opis z Warstwy 1 (mediana/zakres/różnica rozcieńczeń)."
        )
        return {
            "suppressed": True, "excluded_empty_groups": excluded_empty_groups,
            "low_n_bio_warning": low_n_bio_warning, "reason": reason,
        }

    return {
        "suppressed": False, "excluded_empty_groups": excluded_empty_groups,
        "low_n_bio_warning": low_n_bio_warning, "reason": None,
    }


# ============================================================
# ORKIESTRACJA: WARSTWY 1+2+3 RAZEM
# ============================================================

def compare_mic_groups(groups, method="holm"):
    """
    Porównuje >=2 grup MIC (np. różnych substancji dla jednego szczepu).

    groups: dict {etykieta: bio_results} - bio_results to lista wyników
    biologicznych z Fazy 3 (aggregate_technical_to_biological), jedna
    lista na etykietę porównywanej grupy.
    method: korekta wielokrotnych porównań dla >2 grup przy Dunn's test
    ("holm" domyślnie, spójnie z modułem dyfuzji; "fdr_bh"/"bonferroni").

    Zwraca dict:
      - descriptions: {etykieta: describe_mic_group(...)}                    [Warstwa 1,
        dla WSZYSTKICH grup przekazanych w `groups`, w tym pustych (n_bio=0,
        opisywane jako "n/d") - pusta grupa NIGDY nie znika z opisu, nawet
        jeśli jest wykluczona z testu Warstwy 2 (patrz layer2 niżej)]
      - pairwise: {(etykieta_a, etykieta_b): dilution_difference(...)}       [Warstwa 1,
        dla WSZYSTKICH par, zawsze liczone niezależnie od tego, czy test
        w Warstwie 2 się wykonał]
      - layer2: dict z kluczami attempted/suppressed_reason/low_n_bio_warning/
        excluded_empty_groups/exclusion_note/test/statistic/p_value/posthoc/
        correction_method. [Warstwa 2/3] excluded_empty_groups (audyt 1.6)
        to grupy z n_bio=0 wykluczone z TEGO testu - exclusion_note opisuje
        to jawnie i jest ustawiona NIEZALEŻNIE od tego, czy test w ogóle się
        wykonał (attempted True albo False), żeby wykluczenie nigdy nie było
        ciche.
    """
    if len(groups) < 2:
        raise ValueError("compare_mic_groups: potrzeba co najmniej 2 grup do porównania.")

    labels = list(groups.keys())
    entries_by_group = {label: _bio_results_to_entries(bio) for label, bio in groups.items()}
    descriptions = {label: describe_mic_group(bio) for label, bio in groups.items()}

    pairwise = {
        (a, b): dilution_difference(descriptions[a], descriptions[b])
        for a, b in combinations(labels, 2)
    }

    guard = _check_layer3_guards(entries_by_group)
    excluded_empty_groups = guard["excluded_empty_groups"]
    remaining_labels = [label for label in labels if label not in excluded_empty_groups]
    exclusion_note = (
        f"Wykluczono z testu Warstwy 2 (brak żadnej wartości MIC, n_bio=0): "
        f"{', '.join(excluded_empty_groups)}. Warstwa 1 (opis) nadal obejmuje te grupy jako 'n/d'."
        if excluded_empty_groups else None
    )
    layer2 = {
        "attempted": not guard["suppressed"],
        "suppressed_reason": guard["reason"],
        "low_n_bio_warning": guard["low_n_bio_warning"],
        "excluded_empty_groups": excluded_empty_groups,
        "exclusion_note": exclusion_note,
        "test": None, "statistic": None, "p_value": None,
        "posthoc": None, "correction_method": None,
    }

    if not guard["suppressed"]:
        if len(remaining_labels) == 2:
            a, b = remaining_labels
            # Surogat MUSI być liczony na PULI OBU grup naraz (nie osobno na
            # każdej) - inaczej wysoki surogat cenzury z grupy A mógłby
            # wypaść NIŻEJ niż realne wartości grupy B, łamiąc regułę
            # "cenzurowane zawsze na skraju połączonej próby".
            pooled = entries_by_group[a] + entries_by_group[b]
            surrogate = _censoring_surrogate(pooled)
            na = len(entries_by_group[a])
            sa, sb = surrogate[:na], surrogate[na:]
            stat, p = stats.mannwhitneyu(sa, sb, alternative="two-sided")
            layer2["test"] = "Mann-Whitney"
            layer2["statistic"] = float(stat)
            layer2["p_value"] = float(p)
        else:
            pooled_entries, pooled_labels = [], []
            for label in remaining_labels:
                for e in entries_by_group[label]:
                    pooled_entries.append(e)
                    pooled_labels.append(label)
            surrogate = _censoring_surrogate(pooled_entries)
            groups_for_kw = [
                [surrogate[i] for i in range(len(surrogate)) if pooled_labels[i] == label]
                for label in remaining_labels
            ]
            stat, p = stats.kruskal(*groups_for_kw)
            layer2["test"] = "Kruskal-Wallis"
            layer2["statistic"] = float(stat)
            layer2["p_value"] = float(p)
            if p < ALPHA:
                df_posthoc = pd.DataFrame({"value": surrogate, "group": pooled_labels})
                posthoc_df = sp.posthoc_dunn(df_posthoc, val_col="value", group_col="group", p_adjust=method)
                layer2["posthoc"] = posthoc_df
                layer2["correction_method"] = method

    return {"descriptions": descriptions, "pairwise": pairwise, "layer2": layer2}


# ============================================================
# ILORAZ MBC/MIC + KLASYFIKACJA (Faza MBC)
# ============================================================
#
# Wszystko na RÓŻNICY INDEKSÓW ROZCIEŃCZEŃ d = log2(MBC) - log2(MIC)
# (iloraz do wyświetlenia = 2**d). Cenzura obu wielkości jest reprezentowana
# jako PRZEDZIAŁ [lo, hi] w log2 (None = nieskończoność w danym kierunku),
# nie jako punkt - dzięki temu jeden spójny mechanizm obsługuje: obie
# wartości dokładne, tylko MIC cenzurowane, tylko MBC cenzurowane, i obie
# cenzurowane naraz, bez osobnych gałęzi kodu dla każdej kombinacji.

def _endpoint_bounds(result):
    """
    (lo, hi) w log2 dla POJEDYNCZEGO wyniku MIC albo MBC (ten sam kształt
    dict: 'mic_value'/'censored'). None oznacza nieskończoność w danym
    kierunku:
      - dokładny:        (log2(v), log2(v))
      - censored='low'  (<=v): (None, log2(v))   - "co najwyżej v"
      - censored='high' (>v):  (log2(v), None)   - "więcej niż v"
    """
    log2_val = math.log2(result["mic_value"])
    if result["censored"] == "low":
        return (None, log2_val)
    if result["censored"] == "high":
        return (log2_val, None)
    return (log2_val, log2_val)


def _format_ratio(d_exact, d_lo, d_hi):
    """Reprezentacja tekstowa ilorazu 2**d - NIGDY liczba bez prefiksu ≥/≤
    gdy d nie jest dokładne (patrz format_mic_display - ta sama zasada)."""
    if d_exact is not None:
        return f"{2 ** d_exact:g}"
    lo_str = f"{2 ** d_lo:g}"
    if d_hi is None:
        return f"≥{lo_str}"
    hi_str = f"{2 ** d_hi:g}"
    if d_lo <= 0:
        return f"≤{hi_str}"
    return f"{lo_str}–{hi_str}"


def compute_mbc_mic_ratio(mic_result, mbc_result):
    """
    Liczy d = log2(MBC) - log2(MIC) dla JEDNEGO sparowanego powtórzenia
    (mic_result/mbc_result - dowolny poziom: wynik per-wiersz z Fazy 2 albo
    per-bio-rep z Fazy 3, oba mają klucze 'mic_value'/'censored').

    KONTROLA SPÓJNOŚCI (strukturalna, zawsze wykonywana jako pierwsza):
    MBC nie może być mniejsze niż MIC. Wyrażona na PRZEDZIAŁACH: jeśli
    górna granica MBC (mbc_hi) jest mniejsza niż dolna granica MIC (mic_lo),
    to MBC<MIC jest PEWNE niezależnie od dokładnych wartości - błąd danych.
    Ten JEDEN test pokrywa wszystkie wymagane w zadaniu "niemożliwe
    kombinacje cenzury" jednym mechanizmem:
      - MIC dokładne > MBC dokładne (prosty przypadek)
      - MIC ">max_MIC" (mic_lo=max_MIC) razem z jakimkolwiek MBC, którego
        górna granica < max_MIC (MIC ">max" "wymusza" MBC ">max" - inny
        wynik MBC jest niespójny)
      - MBC "<=min_MBC" (mbc_hi=min_MBC) razem z dokładnym MIC > min_MBC

    Gdy spójne: d_lo/d_hi to najciaśniejsza możliwa granica na d, z
    DOMYŚLNYM dolnym ograniczeniem d>=0 (bo MBC>=MIC jest wymuszone -
    nawet gdy surowa cenzura nie daje żadnej dolnej granicy, biologia
    i tak ją daje). MIC, MBC i sam iloraz (d/d_lo/d_hi/ratio_display) są
    liczone TAK SAMO niezależnie od współczynnika rozcieńczenia (Wsp_rozc)
    użytego w danym powtórzeniu - to się NIE zmienia.

    KLASYFIKACJA (audyt 1.4) - NATOMIAST jest podawana TYLKO, gdy zarówno
    MIC, jak i MBC tego powtórzenia pochodzą z serii DWUKROTNYCH rozcieńczeń
    (Wsp_rozc == MBC_MIC_CLASSIFICATION_DILUTION_FACTOR). Uzasadnienie: progi
    d<=2/d>=3 są zdefiniowane w literaturze właśnie dla serii dwukrotnych,
    gdzie każdy możliwy krok d jest liczbą całkowitą - próg nigdy nie trafia
    "między" 2 a 3. Dla innego współczynnika (np. 5-krotnego, gdzie jeden
    krok to log2(5)=2.32) próg wypadałby w tej martwej strefie mimo w pełni
    dokładnego pomiaru, co dawałoby fałszywe "nieoznaczalny" (patrz audyt).
    Gdy współczynnik nie jest dwukrotny: classification="niedostępna" z
    jawnym powodem - MIC/MBC/iloraz nadal są zwracane i pokazywane.

    Gdy klasyfikacja MA zastosowanie (oba Wsp_rozc == 2):
      - d_hi <= MBC_MIC_BACTERICIDAL_MAX_D (2)   -> "bakteriobójcze" PEWNE
        (prawdziwe d nie może przekroczyć d_hi).
      - d_lo >= MBC_MIC_BACTERIOSTATIC_MIN_D (3) -> "bakteriostatyczne" PEWNE
        (prawdziwe d nie może być niższe niż d_lo).
      - w przeciwnym razie: "nieoznaczalny" - granica [d_lo, d_hi] obejmuje
        próg klasyfikacji, NIE zgadujemy który wariant jest prawdziwy.

    Zwraca dict: status (ok/blad_spojnosci/nieoznaczalny/brak_danych),
    d (float|None, tylko gdy dokładne), d_lo, d_hi (float|None),
    classification (bakteriobójcze/bakteriostatyczne/nieoznaczalny/niedostępna/None),
    ratio_display (str), reason (str).
    """
    if mic_result.get("mic_value") is None or mbc_result.get("mic_value") is None:
        return {
            "status": "brak_danych", "d": None, "d_lo": None, "d_hi": None,
            "classification": None, "ratio_display": None,
            "reason": "Brak wartości MIC i/lub MBC dla tego powtórzenia - iloraz nie policzony.",
        }

    mic_lo, mic_hi = _endpoint_bounds(mic_result)
    mbc_lo, mbc_hi = _endpoint_bounds(mbc_result)

    if mbc_hi is not None and mic_lo is not None and mbc_hi < mic_lo:
        return {
            "status": "blad_spojnosci", "d": None, "d_lo": None, "d_hi": None,
            "classification": None, "ratio_display": None,
            "reason": (
                f"Niemożliwa kombinacja: MBC "
                f"({format_mic_display(mbc_result['mic_value'], mbc_result['censored'])}) nie może być "
                f"mniejsze niż MIC ({format_mic_display(mic_result['mic_value'], mic_result['censored'])}) "
                f"dla tego samego powtórzenia - błąd danych. Iloraz nie policzony."
            ),
        }

    d_lo = 0.0
    if mbc_lo is not None and mic_hi is not None:
        d_lo = max(d_lo, mbc_lo - mic_hi)
    d_hi = None if (mbc_hi is None or mic_lo is None) else (mbc_hi - mic_lo)

    d_lo = round(d_lo, 9)
    if d_hi is not None:
        d_hi = round(d_hi, 9)
    d_exact = d_lo if (d_hi is not None and d_lo == d_hi) else None
    ratio_display = _format_ratio(d_exact, d_lo, d_hi)

    # --- Audyt 1.4: klasyfikacja wymaga serii DWUKROTNEJ dla OBU wielkości ---
    mic_wsp = mic_result.get("Wsp_rozc")
    mbc_wsp = mbc_result.get("Wsp_rozc")
    classification_applicable = (
        mic_wsp == MBC_MIC_CLASSIFICATION_DILUTION_FACTOR and mbc_wsp == MBC_MIC_CLASSIFICATION_DILUTION_FACTOR
    )

    if not classification_applicable:
        classification = "niedostępna"
        status = "ok" if d_exact is not None else "nieoznaczalny"
        if mic_wsp is None or mbc_wsp is None:
            wsp_desc = "nieznany"
        elif mic_wsp == mbc_wsp:
            wsp_desc = f"{mic_wsp:g}"
        else:
            wsp_desc = f"MIC={mic_wsp:g}, MBC={mbc_wsp:g}"
        reason = (
            f"Iloraz MBC/MIC = {ratio_display}. Klasyfikacja bakteriobójcze/bakteriostatyczne jest "
            f"zdefiniowana dla serii dwukrotnych rozcieńczeń; przy Wsp_rozc={wsp_desc} nie jest podawana."
        )
    elif d_hi is not None and d_hi <= MBC_MIC_BACTERICIDAL_MAX_D:
        classification = "bakteriobójcze"
        status = "ok"
        reason = (
            f"d <= {d_hi:g} <= {MBC_MIC_BACTERICIDAL_MAX_D:g} (pewne, niezależnie od cenzury) "
            f"-> bakteriobójcze."
        )
    elif d_lo >= MBC_MIC_BACTERIOSTATIC_MIN_D:
        classification = "bakteriostatyczne"
        status = "ok"
        reason = (
            f"d >= {d_lo:g} >= {MBC_MIC_BACTERIOSTATIC_MIN_D:g} (pewne, niezależnie od cenzury) "
            f"-> bakteriostatyczne."
        )
    else:
        classification = "nieoznaczalny"
        status = "nieoznaczalny"
        reason = (
            f"Granica d obejmuje próg klasyfikacji ({MBC_MIC_BACTERICIDAL_MAX_D:g}-"
            f"{MBC_MIC_BACTERIOSTATIC_MIN_D:g}) - nie można jednoznacznie sklasyfikować bez zgadywania. "
            f"Zebrane dane pozwalają jedynie ograniczyć d do przedziału pokazanego w ratio_display."
        )

    return {
        "status": status, "d": d_exact, "d_lo": d_lo, "d_hi": d_hi,
        "classification": classification, "ratio_display": ratio_display,
        "reason": reason,
    }


def pair_mic_mbc_by_bio_rep(mic_bio_results, mbc_bio_results):
    """
    Paruje wyniki biologiczne MIC i MBC (z aggregate_technical_to_biological,
    osobno dla każdego typu) po (Bakteria, Substancja, Rep_biologiczna).

    Zwraca (paired, unmatched):
      - paired: lista {Bakteria, Substancja, Rep_biologiczna, mic, mbc, ratio}
        dla powtórzeń obecnych w OBU zbiorach (ratio = compute_mbc_mic_ratio(...)).
      - unmatched: lista (klucz, powód) dla powtórzeń biologicznych obecnych
        tylko w jednym z dwóch zbiorów - nie da się policzyć ilorazu bez obu
        wartości, więc są pominięte z jawnym powodem, nie po cichu.
    """
    mic_by_key = {(b["Bakteria"], b["Substancja"], b["Rep_biologiczna"]): b for b in mic_bio_results}
    mbc_by_key = {(b["Bakteria"], b["Substancja"], b["Rep_biologiczna"]): b for b in mbc_bio_results}

    paired, unmatched = [], []
    for key in sorted(set(mic_by_key) | set(mbc_by_key)):
        mic_b, mbc_b = mic_by_key.get(key), mbc_by_key.get(key)
        if mic_b is None or mbc_b is None:
            unmatched.append((key, "brak wyniku MIC dla tego powtórzenia" if mic_b is None else "brak wyniku MBC dla tego powtórzenia"))
            continue
        paired.append({
            "Bakteria": key[0], "Substancja": key[1], "Rep_biologiczna": key[2],
            "mic": mic_b, "mbc": mbc_b, "ratio": compute_mbc_mic_ratio(mic_b, mbc_b),
        })
    return paired, unmatched


def _representative_d(ratio_result):
    """
    Punkt reprezentatywny d dla jednego sparowanego powtórzenia biologicznego,
    do dalszej agregacji (mediana grupowa). Uproszczenie świadome: gdy
    klasyfikacja per-powtórzenie rozstrzygnęła się dzięki granicy (nie
    dokładnej wartości), używamy TEJ granicy jako reprezentanta (d_hi dla
    "bakteriobójcze", d_lo dla "bakteriostatyczne") - to nadal wartość, co do
    której wiemy PEWNIE, po której stronie progu leży. "nieoznaczalny" i
    "blad_spojnosci"/"brak_danych" NIE mają defensywnej reprezentacji -
    zwracają None i są wykluczone z mediany grupowej (zliczone osobno).

    "niedostępna" (audyt 1.4 - Wsp_rozc != 2) NADAL zwraca d, gdy jest ono
    dokładne (ratio_result["d"] is not None) - to nadal poprawna, policzona
    liczba, warta zachowania w medianie d dla RAPORTOWANIA (median_d/
    ratio_display). To summarize_mbc_mic_ratio, nie ta funkcja, decyduje,
    czy WYNIKOWA klasyfikacja grupy ma być podana, sprawdzając, czy
    KTÓRYKOLWIEK z użytych wkładów miał classification="niedostępna".
    """
    if ratio_result["status"] in ("blad_spojnosci", "brak_danych"):
        return None
    if ratio_result["d"] is not None:
        return ratio_result["d"]
    if ratio_result["classification"] == "bakteriobójcze":
        return ratio_result["d_hi"]
    if ratio_result["classification"] == "bakteriostatyczne":
        return ratio_result["d_lo"]
    return None


def summarize_mbc_mic_ratio(paired_results):
    """
    Podsumowanie ZBIORCZE ilorazu MBC/MIC dla jednej grupy (Bakteria x
    Substancja): mediana d PO POWTÓRZENIACH BIOLOGICZNYCH (ta sama reguła
    high-median co w Fazie 3/Warstwie 1 Fazy 4 - reużyta wprost,
    _high_median_entry), i klasyfikacja na tej medianie.

    paired_results: lista z pair_mic_mbc_by_bio_rep dla JEDNEJ grupy.

    Zwraca dict: n_bio (liczba powtórzeń z rozstrzygalnym d),
    n_bio_excluded (nieoznaczalne/błędne, wykluczone z mediany),
    median_d, classification, ratio_display, warning (n_bio<2, spójnie z
    resztą programu).
    """
    reps_with_d = [(p, _representative_d(p["ratio"])) for p in paired_results]
    usable = [(p, d) for p, d in reps_with_d if d is not None]
    n_bio_excluded = len(paired_results) - len(usable)

    if not usable:
        return {
            "n_bio": 0, "n_bio_excluded": n_bio_excluded, "median_d": None,
            "classification": None, "ratio_display": None,
            "warning": "Brak żadnego powtórzenia z rozstrzygalną wartością d - nie można podać klasyfikacji zbiorczej.",
        }

    entries = [{"log2": d, "censored": None, "source": p} for p, d in usable]
    median_entry = _high_median_entry(entries)
    median_d = round(median_entry["log2"], 9)

    # Audyt 1.4: jeśli KTÓREKOLWIEK z użytych powtórzeń miało klasyfikację
    # "niedostępna" (Wsp_rozc != 2 dla MIC i/lub MBC tego powtórzenia), cała
    # klasyfikacja GRUPOWA jest też "niedostępna" - median_d/ratio_display
    # nadal są liczone i pokazywane normalnie (to poprawne liczby niezależnie
    # od współczynnika rozcieńczenia), tylko etykieta bakteriobójcze/
    # bakteriostatyczne nie jest zgadywana z mediany zbudowanej częściowo
    # z niesklasyfikowalnych wkładów.
    non_2fold = [p for p, d in usable if p["ratio"].get("classification") == "niedostępna"]

    if non_2fold:
        classification = "niedostępna"
    elif median_d <= MBC_MIC_BACTERICIDAL_MAX_D:
        classification = "bakteriobójcze"
    elif median_d >= MBC_MIC_BACTERIOSTATIC_MIN_D:
        classification = "bakteriostatyczne"
    else:
        classification = "nieoznaczalny"

    n_bio = len(usable)
    warning = MIC_LOW_N_BIO_WARNING if n_bio < 2 else None
    if non_2fold:
        cls_warning = (
            f"Klasyfikacja bakteriobójcze/bakteriostatyczne niedostępna: {len(non_2fold)}/{n_bio} "
            f"użytych powtórzeń nie pochodziło z serii dwukrotnych rozcieńczeń (klasyfikacja "
            f"zdefiniowana tylko dla Wsp_rozc=={MBC_MIC_CLASSIFICATION_DILUTION_FACTOR:g})."
        )
        warning = f"{warning} {cls_warning}" if warning else cls_warning

    return {
        "n_bio": n_bio, "n_bio_excluded": n_bio_excluded, "median_d": median_d,
        "classification": classification, "ratio_display": f"{2 ** median_d:g}",
        "warning": warning,
    }


# ============================================================
# TABELA ZBIORCZA + WIERSZE SZCZEGÓŁOWE (Faza finalna - eksport)
# ============================================================
#
# Czysta budowa danych (bez matplotlib/reportlab) do użytku przez
# mic_plotting.py (opcjonalnie) i reports.py (PDF/Excel) - nic tu nie jest
# liczone od nowa, tylko spinane w tabelaryczny kształt z gotowych wyników
# Faz 2-MBC.

def build_mic_summary_rows(mic_grouped, mbc_grouped=None):
    """
    Tabela zbiorcza "do publikacji": jeden wiersz per (Bakteria, Substancja).
    mic_grouped (wymagane) i mbc_grouped (opcjonalne) to wyniki
    aggregate_all - gdy dana substancja nie ma MBC, kolumny MBC/iloraz są
    puste ("n/d"), nie błędem.

    Dla par obecnych w OBU zbiorach liczy iloraz przez
    pair_mic_mbc_by_bio_rep + summarize_mbc_mic_ratio (Faza MBC, bez zmian).

    Kolumna "Uwagi" zbiera WSZYSTKIE statusy/ostrzeżenia dotyczące danej
    grupy w jeden czytelny tekst (n_bio=1 z MIC i/lub MBC, wykluczone
    powtórzenia biologiczne, odczyty wymaga_weryfikacji, błędy spójności
    MBC<MIC, klasyfikacje nieoznaczalne, niesparowane powtórzenia) - żaden
    z tych sygnałów nie może być widoczny TYLKO w logu.
    """
    mbc_grouped = mbc_grouped or {}
    rows = []
    all_keys = sorted(set(mic_grouped.keys()) | set(mbc_grouped.keys()))
    for key in all_keys:
        bakteria, substancja = key
        mic_entry = mic_grouped.get(key)
        mbc_entry = mbc_grouped.get(key)
        mic_summary = mic_entry["summary"] if mic_entry else None
        mbc_summary = mbc_entry["summary"] if mbc_entry else None

        warnings = []
        if mic_summary and mic_summary.get("warning"):
            warnings.append(f"MIC: {mic_summary['warning']}")
        if mic_summary and mic_summary.get("n_bio_excluded"):
            warnings.append(f"MIC: {mic_summary['n_bio_excluded']} powtórzenie(a) biologiczne bez wartości (wykluczone).")
        if mic_summary and mic_summary.get("n_flagged"):
            warnings.append(f"MIC: {mic_summary['n_flagged']} odczyt(y) techniczne oznaczone jako wymaga_weryfikacji.")
        if mic_summary and mic_summary.get("rep_bio_fallback_warning"):
            warnings.append(f"MIC: {mic_summary['rep_bio_fallback_warning']}")
        if mic_summary and mic_summary.get("unit_warning"):
            warnings.append(f"MIC: {mic_summary['unit_warning']}")
        if mbc_summary and mbc_summary.get("warning"):
            warnings.append(f"MBC: {mbc_summary['warning']}")
        if mbc_summary and mbc_summary.get("n_bio_excluded"):
            warnings.append(f"MBC: {mbc_summary['n_bio_excluded']} powtórzenie(a) biologiczne bez wartości (wykluczone).")
        if mbc_summary and mbc_summary.get("n_flagged"):
            warnings.append(f"MBC: {mbc_summary['n_flagged']} odczyt(y) techniczne oznaczone jako wymaga_weryfikacji.")
        if mbc_summary and mbc_summary.get("rep_bio_fallback_warning"):
            warnings.append(f"MBC: {mbc_summary['rep_bio_fallback_warning']}")
        if mbc_summary and mbc_summary.get("unit_warning"):
            warnings.append(f"MBC: {mbc_summary['unit_warning']}")

        ratio_display, classification = None, None
        if mic_entry and mbc_entry:
            paired, unmatched = pair_mic_mbc_by_bio_rep(mic_entry["bio_results"], mbc_entry["bio_results"])
            errors = [p for p in paired if p["ratio"]["status"] == "blad_spojnosci"]
            if errors:
                warnings.append(
                    f"BŁĄD SPÓJNOŚCI w {len(errors)} powtórzeniu/powtórzeniach biologicznych "
                    f"(MBC mniejsze niż MIC) - iloraz dla nich nie policzony, patrz szczegóły."
                )
            if unmatched:
                warnings.append(f"{len(unmatched)} powtórzenie(a) biologiczne bez pary MIC/MBC (pominięte w ilorazie).")
            if paired:
                ratio_summary = summarize_mbc_mic_ratio(paired)
                ratio_display = ratio_summary["ratio_display"]
                classification = ratio_summary["classification"]
                if classification == "nieoznaczalny":
                    warnings.append("Klasyfikacja MBC/MIC nieoznaczalna - granica obejmuje próg, patrz iloraz.")
                if ratio_summary.get("warning"):
                    warnings.append(f"Iloraz MBC/MIC: {ratio_summary['warning']}")

        rows.append({
            "Bakteria": bakteria, "Substancja": substancja,
            "n_bio_MIC": mic_summary["n_bio"] if mic_summary else 0,
            "MIC": mic_summary["median"]["display"] if mic_summary else "n/d",
            "n_bio_MBC": mbc_summary["n_bio"] if mbc_summary else 0,
            "MBC": mbc_summary["median"]["display"] if mbc_summary else "n/d",
            "Iloraz_MBC_MIC": ratio_display or "n/d",
            "Klasyfikacja": classification or "n/d",
            "Uwagi": "; ".join(warnings) if warnings else "",
        })
    return rows


def build_mic_replicate_rows(mic_grouped, mbc_grouped=None):
    """
    Wiersze SZCZEGÓŁOWE (jeden na powtórzenie biologiczne) do arkusza
    Excela - odpowiednik "Dane Surowe" z modułu dyfuzji, ale już na
    poziomie powtórzeń biologicznych (po agregacji technicznej, Faza 3),
    z osobną kolumną Typ (MIC/MBC) i statusem/powodem widocznym wprost.
    """
    mbc_grouped = mbc_grouped or {}
    rows = []
    for typ, grouped in (("MIC", mic_grouped), ("MBC", mbc_grouped)):
        for (bakteria, substancja), entry in grouped.items():
            for b in entry["bio_results"]:
                wartosc = format_mic_display(b["mic_value"], b["censored"]) if b["mic_value"] is not None else "n/d"
                rows.append({
                    "Bakteria": bakteria, "Substancja": substancja, "Typ": typ,
                    "Rep_biologiczna": b.get("Rep_biologiczna"),
                    "Wartosc": wartosc, "Status": b.get("status"), "Powod": b.get("reason", ""),
                })
    return rows
