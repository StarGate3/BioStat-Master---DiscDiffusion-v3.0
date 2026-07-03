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
from config import (
    MIC_STATUS_OK, MIC_STATUS_NEEDS_REVIEW, MIC_STATUS_CENSORED_LOW, MIC_STATUS_CENSORED_HIGH,
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
