import io
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from config import COL_GROUP, PDF_DPI


def _add_plot_to_pdf(elements, styles, fig, title):
    """
    Dodaje wykres (matplotlib Figure) do listy elementów PDF - WSPÓLNA dla
    sekcji dyfuzji i sekcji MIC/MBC (Faza finalna reużywa ją wprost, bez
    duplikowania logiki skalowania/osadzania obrazka), żeby styl osadzania
    (skalowanie do szerokości strony, limit wysokości dla bardzo wysokich
    wykresów typu heatmapa) był spójny w całym raporcie.
    """
    if not fig:
        return
    elements.append(Paragraph(title, styles['Heading2']))
    elements.append(Spacer(1, 6))

    img_buf = io.BytesIO()
    fig.savefig(img_buf, format='png', dpi=PDF_DPI, bbox_inches='tight')
    img_buf.seek(0)

    img = Image(img_buf)
    aspect = img.imageHeight / img.imageWidth
    target_width = 6 * inch
    target_height = target_width * aspect

    if target_height > 9 * inch:
        ratio = (9 * inch) / target_height
        target_height = 9 * inch
        target_width = target_width * ratio

    img.drawWidth = target_width
    img.drawHeight = target_height

    elements.append(img)
    elements.append(Spacer(1, 12))


def _build_mic_mbc_elements(styles, f_norm, f_bold, mic_mbc_data):
    """
    Sekcja MIC/MBC raportu PDF (Faza finalna modułu MIC/MBC). Reużywa
    _add_plot_to_pdf (ta sama funkcja co sekcja dyfuzji) i styl tabel/
    ostrzeżeń modułu dyfuzji (styl 'Warning' zdefiniowany w generate_pdf).

    mic_mbc_data: dict {
        'bact': str,
        'table_rows': lista z mic_logic.build_mic_summary_rows (tabela
            zbiorcza - WYJŚCIE 4),
        'figures': lista (tytuł, matplotlib Figure) - WYJŚCIA 1-3, w
            kolejności do wstawienia,
    }

    Kolumna "Uwagi" tabeli zbiorczej jest wypisywana OSOBNO jako lista
    ostrzeżeń pod tabelą (styl 'Warning') - wymóg nadrzędny: żaden status/
    ostrzeżenie nie może być widoczny TYLKO w logu ani zgnieciony do
    nieczytelnej komórki tabeli.
    """
    elements = []
    elements.append(Paragraph(f"Sekcja MIC/MBC: {mic_mbc_data.get('bact', '')}", styles['Heading2']))
    elements.append(Spacer(1, 12))

    table_rows = mic_mbc_data.get('table_rows') or []
    if table_rows:
        header = ['Substancja', 'n_bio\nMIC', 'MIC', 'n_bio\nMBC', 'MBC', 'Iloraz\nMBC/MIC', 'Klasyfikacja']
        table_data = [header]
        warning_lines = []
        for row in table_rows:
            table_data.append([
                row['Substancja'], str(row['n_bio_MIC']), row['MIC'],
                str(row['n_bio_MBC']), row['MBC'], row['Iloraz_MBC_MIC'], row['Klasyfikacja'],
            ])
            if row.get('Uwagi'):
                warning_lines.append(f"• <b>{row['Substancja']}</b>: {row['Uwagi']}")

        t = Table(table_data, colWidths=[85, 38, 62, 38, 62, 58, 85])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), f_bold),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 10),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTNAME', (0, 1), (-1, -1), f_norm),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
        ]))
        elements.append(t)
        elements.append(Spacer(1, 16))

        if warning_lines:
            elements.append(Paragraph("Statusy i ostrzeżenia (per substancja):", styles['Heading2']))
            elements.append(Spacer(1, 6))
            for line in warning_lines:
                elements.append(Paragraph(line, styles['Warning']))
                elements.append(Spacer(1, 4))
            elements.append(Spacer(1, 12))
    else:
        elements.append(Paragraph("Brak danych MIC/MBC do zestawienia dla tego szczepu.", styles['Normal']))
        elements.append(Spacer(1, 12))

    for title, fig in (mic_mbc_data.get('figures') or []):
        _add_plot_to_pdf(elements, styles, fig, title)

    return elements


def generate_pdf(file_path, metadata, stats_summary=None, figures=None, detailed_results=None, mic_mbc_data=None):
    """
    Generuje raport PDF.

    Args:
        file_path (str): Ścieżka do zapisu pliku.
        metadata (dict): Dane o dacie, bakterii, ref group (dot. sekcji dyfuzji).
        stats_summary (DataFrame | None): Tabela statystyk opisowych (dyfuzja).
            None -> pomija całą sekcję dyfuzji (raport MIC/MBC-only).
        figures (dict | None): Słownik obiektów matplotlib Figure (dyfuzja).
        detailed_results (list | None): Wyniki post-hoc/effect size (dyfuzja).
        mic_mbc_data (dict | None): Patrz _build_mic_mbc_elements. None -> brak
            sekcji MIC/MBC. Gdy podane RAZEM z danymi dyfuzji, raport zawiera
            OBIE sekcje w jednym pliku PDF (jak dla szczepu z danymi obu metod).
    """
    try:
        doc = SimpleDocTemplate(file_path, pagesize=A4)

        # Konfiguracja czcionek
        try:
            pdfmetrics.registerFont(TTFont('Arial', 'arial.ttf'))
            pdfmetrics.registerFont(TTFont('Arial-Bold', 'arialbd.ttf'))
            f_norm = 'Arial'
            f_bold = 'Arial-Bold'
        except Exception:
            f_norm = 'Helvetica'
            f_bold = 'Helvetica-Bold'

        styles = getSampleStyleSheet()
        styles['Normal'].fontName = f_norm
        styles['Heading2'].fontName = f_bold
        styles['Title'].fontName = f_bold
        styles.add(ParagraphStyle(
            name='Warning', parent=styles['Normal'], fontName=f_bold,
            textColor=colors.HexColor('#8B0000'), borderColor=colors.HexColor('#8B0000'),
            borderWidth=1, borderPadding=6, backColor=colors.HexColor('#FFF0F0'),
        ))

        elements = []
        has_diffusion = stats_summary is not None
        has_mic_mbc = mic_mbc_data is not None

        # 1. Tytuł - odzwierciedla, co faktycznie jest w raporcie
        if has_diffusion and has_mic_mbc:
            title = "Raport z analizy: Dyfuzja krążkowa + MIC/MBC"
        elif has_diffusion:
            title = "Raport z analizy Disk Diffusion"
        else:
            title = "Raport z analizy MIC/MBC"
        elements.append(Paragraph(title, styles['Title']))
        elements.append(Spacer(1, 12))

        if has_diffusion:
            # 2. Metryczka
            meta_text = (
                f"<b>Data:</b> {metadata['date']}<br/>"
                f"<b>Bakteria:</b> {metadata['bact']}<br/>"
                f"<b>Post-hoc:</b> {metadata['method']}<br/>"
                f"<b>Grupa referencyjna (kontrola, do której porównywano istotność):</b> {metadata['ref']}"
            )
            elements.append(Paragraph(meta_text, styles['Normal']))
            elements.append(Spacer(1, 12))

            if metadata.get('test_used') == "ANOVA":
                elements.append(Paragraph(
                    "UWAGA: dla testu ANOVA porównania parami wykonano testem Tukey HSD, który ma "
                    f"wbudowaną własną korektę wielokrotnych porównań. Wybrana metoda post-hoc "
                    f"('{metadata['method']}') NIE ma tu zastosowania - dotyczy wyłącznie ścieżki "
                    "Kruskal-Wallis/Dunn.",
                    styles['Normal']
                ))
                elements.append(Spacer(1, 12))

            elements.append(Paragraph(
                "Niniejsza analiza porównuje średnice stref zahamowania metodami statystycznymi "
                "(test istotności + wielkość efektu) i NIE wylicza klinicznych kategorii S/I/R "
                "(Susceptible/Intermediate/Resistant) wg CLSI (M100) ani EUCAST.",
                styles['Normal']
            ))
            elements.append(Spacer(1, 24))

            if metadata.get('low_n_bio_warning'):
                elements.append(Paragraph(
                    "⚠ WYNIKI ORIENTACYJNE: co najmniej jedna porównywana grupa nie ma replikacji "
                    "biologicznej (n_bio&lt;2 - patrz kolumna n_bio w tabeli poniżej). P-value i istotność "
                    "statystyczna poniżej NIE są potwierdzone niezależnymi powtórzeniami biologicznymi; "
                    "traktuj je wyłącznie jako wskazówkę, nie dowód.",
                    styles['Warning']
                ))
                elements.append(Spacer(1, 18))

            # 3. Tabela Statystyk
            if stats_summary is not None:
                elements.append(Paragraph("Statystyki Opisowe", styles['Heading2']))
                table_data = [[COL_GROUP, 'Średnia (mm)', 'SD (między-biologiczne)', 'n_bio', 'n_tech']]
                for index, row in stats_summary.iterrows():
                    sd_txt = f"{row['sd_bio']:.2f}" if pd.notna(row['sd_bio']) else "—"
                    table_data.append([
                        row[COL_GROUP], f"{row['mean']:.2f}", sd_txt,
                        f"{int(row['n_bio'])}", f"{int(row['n_tech'])}",
                    ])
                t = Table(table_data, colWidths=[170, 70, 110, 45, 45])
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                    ('FONTNAME', (0, 0), (-1, 0), f_bold),
                    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                    ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                    ('GRID', (0, 0), (-1, -1), 1, colors.black),
                    ('FONTNAME', (0, 1), (-1, -1), f_norm)
                ]))
                elements.append(t)
                elements.append(Spacer(1, 24))

            # 4. Wykresy
            figures = figures or {}
            if figures.get('bar'):
                _add_plot_to_pdf(elements, styles, figures['bar'], "Wykres Porównawczy (Główny)")
            if figures.get('effect'):
                _add_plot_to_pdf(elements, styles, figures['effect'], "Analiza Wielkości Efektu (Cohen's d)")
            if figures.get('heat'):
                _add_plot_to_pdf(elements, styles, figures['heat'], "Mapa Ciepła (Aktywność)")
            if figures.get('pvalue'):
                _add_plot_to_pdf(elements, styles, figures['pvalue'], "Mapa Istotności Statystycznych (P-value)")
            if figures.get('trend'):
                _add_plot_to_pdf(elements, styles, figures['trend'], "Trend Zależności od Dawki")
            if figures.get('cross'):
                _add_plot_to_pdf(elements, styles, figures['cross'], "Porównanie Międzygatunkowe")
            if figures.get('pca'):
                _add_plot_to_pdf(elements, styles, figures['pca'], "Analiza PCA (Główne Składowe)")

            # 5. Werdykt
            elements.append(Paragraph("Werdykt Statystyczny (Istotne różnice)", styles['Heading2']))
            verdicts = []

            if detailed_results:
                for row in detailed_results:
                    if row['Significant']:
                        d_val = row["Cohen's d"]
                        interp = row["Effect Size"]
                        # d_val jest NaN gdy n<2 w ktorejs grupie (brak replikacji
                        # biologicznej) - "d=nan" wygladaloby jak blad, wiec
                        # pokazujemy samo slowne wyjasnienie zamiast liczby.
                        if d_val != d_val:  # NaN check bez importu math/pandas
                            v_text = f"• Istotna różnica: <b>{row['Group 1']}</b> vs <b>{row['Group 2']}</b> (p={row['P-adj']:.4f}). Wielkość efektu: {interp}."
                        else:
                            v_text = f"• Istotna różnica: <b>{row['Group 1']}</b> vs <b>{row['Group 2']}</b> (p={row['P-adj']:.4f}). Wielkość efektu d={d_val:.2f} ({interp})."
                        verdicts.append(v_text)

            if not verdicts:
                elements.append(Paragraph("Nie stwierdzono różnic istotnych statystycznie.", styles['Normal']))
            else:
                for v in verdicts:
                    elements.append(Paragraph(v, styles['Normal']))
            elements.append(Spacer(1, 12))

        # 6. Sekcja MIC/MBC (Faza finalna)
        if has_mic_mbc:
            elements.extend(_build_mic_mbc_elements(styles, f_norm, f_bold, mic_mbc_data))

        doc.build(elements)
        return True, "Raport PDF został wygenerowany!"

    except Exception as e:

        return False, str(e)


def export_mic_mbc_excel(file_path, table_rows, replicate_rows=None, warnings=None):
    """
    Eksport MIC/MBC do Excela - ten sam wzorzec co gui.py.export_to_excel
    dla dyfuzji: arkusz tabeli zbiorczej, arkusz powtórzeń biologicznych
    (odpowiednik "Dane Surowe"), i osobny JAWNY arkusz "UWAGA - MIC_MBC"
    z ostrzeżeniami - analogicznie do "UWAGA - n_bio" w module dyfuzji,
    żeby żadne ostrzeżenie nie było widoczne tylko wewnątrz komórek tabeli.
    """
    with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
        pd.DataFrame(table_rows).to_excel(writer, sheet_name="Tabela zbiorcza MIC-MBC", index=False)
        if replicate_rows:
            pd.DataFrame(replicate_rows).to_excel(writer, sheet_name="Powtorzenia biologiczne", index=False)
        if warnings:
            pd.DataFrame({"Uwaga": warnings}).to_excel(writer, sheet_name="UWAGA - MIC_MBC", index=False)
