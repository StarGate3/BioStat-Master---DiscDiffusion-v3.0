"""
Wykresy MIC/MBC (Faza finalna modułu MIC/MBC).

Odpowiednik plotting.py dla danych MIC/MBC - ale zamiast surowych wierszy
DataFrame (jak Plotter), operuje na WYNIKACH z mic_logic.py (bio_results z
Fazy 3, opisy/porównania z Fazy 4, ilorazy z Fazy MBC). Żadna z funkcji tu
NIE liczy niczego od nowa - tylko wizualizuje gotowe wyniki.

ZASADY NADRZĘDNE (patrz też zadanie):
- Skala MIC/MBC jest zawsze logarytmiczna (log2) - oś pokazuje rzeczywiste
  stężenia, nigdy log2 jako liczbę. Patrz _setup_log2_axis.
- Wartości cenzurowane są rysowane jako trójkąt-strzałka (skierowany w
  stronę, w którą "ucieka" nieznana prawdziwa wartość), NIGDY jako zwykły
  punkt na granicy - patrz MARKER_BY_CENSOR.
- Statusy/ostrzeżenia z poprzednich faz (n_bio=1, wymaga_weryfikacji,
  nieoznaczalny, nieważny przebieg) muszą być widoczne NA figurze (podpis/
  adnotacja), nie tylko w logu - ten sam standard co plotting.draw_bar_plot
  (baner n_bio=1 wypalony na figurze).
"""
import math

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

import mic_logic
from config import SCREEN_DPI

MIC_COLOR = "#1F6AA5"
MBC_COLOR = "#8B0000"
CLASSIFICATION_COLORS = {
    "bakteriobójcze": "#2E7D32",
    "bakteriostatyczne": "#EF6C00",
    "nieoznaczalny": "#9E9E9E",
}
MARKER_BY_CENSOR = {None: "o", "low": "<", "high": ">"}


def _setup_log2_axis(ax, values):
    """
    Oś logarytmiczna o podstawie 2 z etykietami RZECZYWISTYCH stężeń (nie
    log2!) - wymóg nadrzędny zadania. Wartości cenzurowane mają swoje
    WŁASNE punkty na tej samej osi (patrz MARKER_BY_CENSOR) - to jedyna oś
    używana przez cały moduł MIC/MBC, nigdy skala liniowa.
    """
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda x, pos: ""))
    vmin, vmax = min(values), max(values)
    ax.set_xlim(vmin / 2.5, vmax * 2.5)


def _status_annotation_text(bio_results):
    """
    Zbiera z listy wyników biologicznych/wierszowych WSZYSTKIE statusy inne
    niż "ok", żeby można je było wypisać wprost na figurze - wymóg
    nadrzędny: żaden status/ostrzeżenie nie może być widoczny TYLKO w logu.
    """
    problems = [b for b in bio_results if b.get("status") not in (mic_logic.MIC_STATUS_OK, None)]
    if not problems:
        return None
    counts = {}
    for b in problems:
        counts[b["status"]] = counts.get(b["status"], 0) + 1
    return "Statusy: " + ", ".join(f"{status} (n={n})" for status, n in counts.items())


def draw_mic_mbc_distribution(bact, substancja, mic_bio_results, mbc_bio_results=None, config=None):
    """
    WYJŚCIE 1: rozkład wartości powtórzeń biologicznych MIC (i opcjonalnie
    MBC, obok siebie) na skali log2, per substancja x szczep. Mediana
    (◆) i zakres (linia) zaznaczone przez reużycie describe_mic_group -
    bez przeliczania niczego od nowa. Cenzura pokazana jako trójkąt.
    """
    config = config or {}
    f_lbl = config.get("font_labels", 10)
    f_ttl = config.get("font_title", 12)

    has_mbc = bool(mbc_bio_results)
    categories = ["MIC", "MBC"] if has_mbc else ["MIC"]
    y_positions = {cat: i for i, cat in enumerate(categories)}

    fig = plt.Figure(figsize=(8, 1.8 + 1.1 * len(categories)), dpi=SCREEN_DPI)
    ax = fig.add_subplot(111)

    all_vals = [b["mic_value"] for b in mic_bio_results if b.get("mic_value") is not None]
    if has_mbc:
        all_vals += [b["mic_value"] for b in mbc_bio_results if b.get("mic_value") is not None]

    if not all_vals:
        ax.text(0.5, 0.5, "Brak danych do wyświetlenia.", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"{bact} / {substancja}: rozkład MIC" + (" i MBC" if has_mbc else ""), fontsize=f_ttl + 2)
        fig.tight_layout()
        return fig

    annotation_lines = []

    def _plot_series(bio_results, y_pos, color, endpoint_name):
        for b in bio_results:
            if b.get("mic_value") is None:
                continue
            jitter = (hash((b.get("Rep_biologiczna"), endpoint_name)) % 1000 / 1000 - 0.5) * 0.24
            marker = MARKER_BY_CENSOR.get(b["censored"], "o")
            ax.scatter(
                b["mic_value"], y_pos + jitter, marker=marker, s=70, color=color,
                edgecolor="black", linewidth=0.6, zorder=5, alpha=0.85,
            )

        desc = mic_logic.describe_mic_group(bio_results)
        if desc["n_bio"] == 0:
            return
        min_val = desc["range_min"]["mic_value"]
        max_val = desc["range_max"]["mic_value"]
        med_val = desc["median"]["mic_value"]
        med_marker = MARKER_BY_CENSOR.get(desc["median"]["censored"], "D")

        ax.plot([min_val, max_val], [y_pos, y_pos], color=color, linewidth=2, alpha=0.5, zorder=2)
        ax.scatter([med_val], [y_pos], marker="D" if med_marker == "o" else med_marker,
                   s=150, color=color, edgecolor="black", linewidth=1.3, zorder=6)
        ax.text(
            med_val, y_pos - 0.32,
            f"mediana: {desc['median']['display']}  (n_bio={desc['n_bio']}, zakres "
            f"{desc['range_min']['display']}–{desc['range_max']['display']})",
            fontsize=f_lbl - 1, ha="center", va="top", color=color, fontweight="bold",
        )

        status_txt = _status_annotation_text(bio_results)
        if status_txt:
            annotation_lines.append(f"{endpoint_name}: {status_txt}")
        summary = mic_logic.summarize_mic_group(bio_results) if bio_results else None
        if summary and summary.get("warning"):
            annotation_lines.append(f"{endpoint_name}: {summary['warning']}")

    _plot_series(mic_bio_results, y_positions["MIC"], MIC_COLOR, "MIC")
    if has_mbc:
        _plot_series(mbc_bio_results, y_positions["MBC"], MBC_COLOR, "MBC")

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()), fontsize=f_lbl + 1, fontweight="bold")
    ax.set_ylim(-0.9, len(categories) - 0.1)

    _setup_log2_axis(ax, all_vals)
    ax.set_xlabel("Stężenie", fontsize=f_ttl)
    ax.set_title(
        f"{bact} / {substancja}: rozkład MIC" + (" i MBC" if has_mbc else "")
        + "  (●/▲/▼ = powtórzenie biologiczne, ◆ = mediana)",
        fontsize=f_ttl, fontweight="bold",
    )

    if annotation_lines:
        fig.subplots_adjust(bottom=fig.subplotpars.bottom + 0.10 * len(annotation_lines))
        fig.text(
            0.5, 0.01, "\n".join(annotation_lines), ha="center", va="bottom",
            fontsize=max(8, f_lbl - 1), color="darkred", fontweight="bold", wrap=True,
        )
    else:
        fig.tight_layout()
    return fig
