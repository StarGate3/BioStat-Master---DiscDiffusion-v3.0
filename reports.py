import io
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from config import COL_GROUP, PDF_DPI

def generate_pdf(file_path, metadata, stats_summary, figures, detailed_results):
    """
    Generuje raport PDF.
    
    Args:
        file_path (str): Ścieżka do zapisu pliku.
        metadata (dict): Dane o dacie, bakterii, ref group.
        stats_summary (DataFrame): Tabela statystyk opisowych.
        figures (dict): Słownik obiektów matplotlib Figure.
        detailed_results (list): Lista słowników z wynikami post-hoc/effect size.
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
        
        elements = []

        # 1. Metryczka
        elements.append(Paragraph("Raport z analizy Disk Diffusion", styles['Title']))
        elements.append(Spacer(1, 12))
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

        # 2. Tabela Statystyk
        if stats_summary is not None:
            elements.append(Paragraph("Statystyki Opisowe", styles['Heading2']))
            table_data = [[COL_GROUP, 'Średnia (mm)', 'SD (±)', 'N']]
            for index, row in stats_summary.iterrows():
                table_data.append([row[COL_GROUP], f"{row['mean']:.2f}", f"{row['std']:.2f}", f"{int(row['count'])}"])
            t = Table(table_data, colWidths=[200, 80, 80, 50])
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

        # 3. Wykresy
        def add_plot_to_pdf(fig, title):
            """Pomocnicza funkcja dodająca wykres do listy elementów PDF"""
            if fig:
                elements.append(Paragraph(title, styles['Heading2']))
                elements.append(Spacer(1, 6))
                
                img_buf = io.BytesIO()
                # Zapisujemy wykres do bufora pamięci
                fig.savefig(img_buf, format='png', dpi=PDF_DPI, bbox_inches='tight')
                img_buf.seek(0)
                
                img = Image(img_buf)
                # Skalowanie obrazka, aby mieścił się na stronie A4 (szerokość ok. 6 cali)
                aspect = img.imageHeight / img.imageWidth
                target_width = 6 * inch
                target_height = target_width * aspect
                
                # Zabezpieczenie przed zbyt wysokimi wykresami (np. heatmapa)
                if target_height > 9 * inch:
                    ratio = (9 * inch) / target_height
                    target_height = 9 * inch
                    target_width = target_width * ratio

                img.drawWidth = target_width
                img.drawHeight = target_height
                
                elements.append(img)
                elements.append(Spacer(1, 12))

        # --- DODAWANIE WSZYSTKICH WYKRESÓW ---
        if figures.get('bar'): 
            add_plot_to_pdf(figures['bar'], "Wykres Porównawczy (Główny)")
            
        if figures.get('effect'): 
            add_plot_to_pdf(figures['effect'], "Analiza Wielkości Efektu (Cohen's d)")
            
        if figures.get('heat'): 
            add_plot_to_pdf(figures['heat'], "Mapa Ciepła (Aktywność)")
            
        if figures.get('pvalue'): 
            add_plot_to_pdf(figures['pvalue'], "Mapa Istotności Statystycznych (P-value)")
            
        if figures.get('trend'): 
            add_plot_to_pdf(figures['trend'], "Trend Zależności od Dawki")
            
        if figures.get('cross'):
            add_plot_to_pdf(figures['cross'], "Porównanie Międzygatunkowe")

        if figures.get('pca'):
            add_plot_to_pdf(figures['pca'], "Analiza PCA (Główne Składowe)")

        # 4. Werdykt
        elements.append(Paragraph("Werdykt Statystyczny (Istotne różnice)", styles['Heading2']))
        verdicts = []
        
        if detailed_results:
            for row in detailed_results:
                if row['Significant']:
                    d_val = row["Cohen's d"]
                    interp = row["Effect Size"]
                    v_text = f"• Istotna różnica: <b>{row['Group 1']}</b> vs <b>{row['Group 2']}</b> (p={row['P-adj']:.4f}). Wielkość efektu d={d_val:.2f} ({interp})."
                    verdicts.append(v_text)

        if not verdicts: 
            elements.append(Paragraph("Nie stwierdzono różnic istotnych statystycznie.", styles['Normal']))
        else:
            for v in verdicts: 
                elements.append(Paragraph(v, styles['Normal']))

        doc.build(elements)
        return True, "Raport PDF został wygenerowany!"
        
    except Exception as e:

        return False, str(e)
