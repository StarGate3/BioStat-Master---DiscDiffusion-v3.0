import customtkinter as ctk

# ======================================================
# OKNO DIALOGOWE - OUTLIERY (DIXON)
# ======================================================
class OutlierDialog(ctk.CTkToplevel):
    def __init__(self, parent, outlier_data):
        super().__init__(parent)
        self.title("Wykryto wartości odstające")
        self.geometry("500x400")
        
        # Okno zawsze na wierzchu
        self.lift()
        self.attributes("-topmost", True)
        
        self.result = [] 

        ctk.CTkLabel(self, text="Wykryto potencjalne błędy pomiarowe (Test Dixona).\nZaznacz wartości, które chcesz WYKLUCZYĆ z analizy:", 
                      font=ctk.CTkFont(size=14, weight="bold"), wraplength=450).pack(pady=10)

        self.scroll = ctk.CTkScrollableFrame(self, width=450, height=250)
        self.scroll.pack(pady=5, padx=10, fill="both", expand=True)

        self.check_vars = {}

        for item in outlier_data:
            group = item['group']
            value = item['value']
            others = item['others']
            desc = f"Grupa: {group}\nOdstający: {value} mm (Pozostałe: {others})"
            
            var = ctk.IntVar(value=1)
            chk = ctk.CTkCheckBox(self.scroll, text=desc, variable=var, font=ctk.CTkFont(size=12))
            chk.pack(anchor="w", pady=5, padx=5)
            self.check_vars[(group, value)] = var

        btn_frame = ctk.CTkFrame(self, fg_color="transparent")
        btn_frame.pack(pady=10, fill="x")

        ctk.CTkButton(btn_frame, text="Potwierdź i Analizuj", fg_color="green", command=self.confirm).pack(side="right", padx=20)
        ctk.CTkButton(btn_frame, text="Ignoruj wszystkie", fg_color="gray", command=self.cancel).pack(side="right", padx=10)

    def confirm(self):
        for (group, val), var in self.check_vars.items():
            if var.get() == 1:
                self.result.append({'Group': group, 'Srednica_mm': val})
        self.destroy()

    def cancel(self):
        self.result = [] 
        self.destroy()


# ======================================================
# OKNO POMOCY / PODRĘCZNIK (FINAL)
# ======================================================
class HelpDialog(ctk.CTkToplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Podręcznik Metodyczny BioStat")
        self.geometry("800x900")
        
        # --- OKNO ZAWSZE NA WIERZCHU ---
        self.lift()
        self.attributes("-topmost", True)
        
        # Kontener z przewijaniem
        self.scroll = ctk.CTkScrollableFrame(self, fg_color="transparent")
        self.scroll.pack(fill="both", expand=True, padx=15, pady=15)

        # Tytuł główny
        ctk.CTkLabel(self.scroll, text="Przewodnik po Analizie i Interpretacji", 
                     font=ctk.CTkFont(size=22, weight="bold")).pack(pady=(10, 5))
        
        ctk.CTkLabel(self.scroll, text="Jak czytać wyniki i rozumieć wykresy generowane przez aplikację.", 
                     text_color="gray").pack(pady=(0, 20))

        # --- SEKCJA 1: ALGORYTM ---
        self.add_section("1. ALGORYTM DZIAŁANIA PROGRAMU")
        self.add_text("Program automatycznie dobiera odpowiedni test statystyczny, podążając za poniższą logiką (drzewo decyzyjne):")
        
        self.add_entry("Krok 1: Normalność (Shapiro-Wilk)", 
                       "Sprawdzamy, czy dane w każdej grupie układają się w 'krzywą dzwonową'.\n"
                       "• p > 0.05: Rozkład normalny (OK).\n"
                       "• p < 0.05: Rozkład inny niż normalny (częste w małych próbach n=3).")
        
        self.add_entry("Krok 2: Wariancja (Levene)", 
                       "Sprawdzamy, czy grupy mają podobny 'rozrzut' wyników.\n"
                       "Jeśli wariancje są różne, testy parametryczne mogą dawać błędne wyniki.")
        
        self.add_entry("Krok 3: Wybór Testu Głównego", 
                       "• ANOVA: Wybierana, gdy dane są normalne i mają równą wariancję (największa moc).\n"
                       "• Kruskal-Wallis: Wybierany, gdy założenia ANOVA nie są spełnione (bezpieczniejszy dla danych mikrobiologicznych).")

        # --- SEKCJA 2: KOREKTY POST-HOC ---
        self.add_section("2. KOREKTY POST-HOC (Którą wybrać?)")
        self.add_text("Gdy porównujesz wiele próbek naraz, rośnie ryzyko, że przypadkowo znajdziesz 'istotną' różnicę (Błąd I rodzaju). Korekty temu zapobiegają.")
        
        self.add_entry("Holm (Zalecana)", 
                       "Najlepszy balans. Jest silniejsza niż brak korekty, ale nie tak 'brutalna' jak Bonferroni. "
                       "Dobra do większości standardowych badań.")
        
        self.add_entry("Bonferroni", 
                       "Bardzo konserwatywna. Bardzo trudno uzyskać p < 0.05. "
                       "Stosuj tylko, gdy musisz mieć absolutną pewność i chcesz uniknąć fałszywych alarmów za wszelką cenę.")
        
        self.add_entry("FDR (Benjamini-Hochberg)", 
                       "Najmniej rygorystyczna. Dopuszcza pewien odsetek fałszywych odkryć. "
                       "Idealna do 'screeningu' (przesiewu) setek substancji, gdy nie chcesz przegapić niczego potencjalnie ciekawego.")

        # --- SEKCJA 3: INTERPRETACJA WYKRESÓW ---
        self.add_section("3. ANATOMIA WYKRESÓW")

        self.add_entry("Wykres Słupkowy (Barplot)", 
                       "• Wysokość słupka: Średnia arytmetyczna strefy zahamowania.\n"
                       "• Antenka (Słupek błędu): Pokazuje zmienność (SD - Odchylenie Standardowe). Im krótsza antenka, tym bardziej powtarzalne były wyniki.")

        self.add_entry("Wykres Pudełkowy (Boxplot)", 
                       "Bardziej szczegółowy niż słupkowy:\n"
                       "• Linia w środku pudełka: Mediana (wartość środkowa).\n"
                       "• Pudełko: Obejmuje 50% środkowych wyników (od 25. do 75. percentyla).\n"
                       "• Wąsy: Zasięg danych (min-max), z wyłączeniem wartości odstających.")

        self.add_entry("Wykres 'Lollipop' (Wielkość Efektu)", 
                       "Najważniejszy wykres do oceny 'siły' działania.\n"
                       "• Oś pozioma (Cohen's d): Mówi, ile 'odchyleń standardowych' dzieli dwie grupy.\n"
                       "• Kropka ZIELONA (W prawo): Grupa badana jest lepsza/silniejsza.\n"
                       "• Kropka CZERWONA (W lewo): Grupa badana jest gorsza/słabsza.\n"
                       "UWAGA: Wykres prezentuje wyłącznie pary różniące się istotnie statystycznie (p < 0.05), "
                       "aby zachować czytelność. Pełne wyniki dla wszystkich par (również nieistotnych) "
                       "znajdziesz w raporcie Excel (zakładka 'Post-hoc Details').")

        self.add_entry("Mapa Ciepła (Heatmap)", 
                       "Wizualizacja macierzy. Kolory ułatwiają szybkie wyłapanie liderów.\n"
                       "• Jasne/Ciepłe pola: Duże strefy zahamowania (aktywność).\n"
                       "• Ciemne pola: Brak aktywności.")
        
        self.add_entry("Trend (Dawka-Odpowiedź)", 
                       "Jeśli w nazwach grup są stężenia (np. 5mg, 10mg), ten wykres pokaże linię.\n"
                       "• Zacieniony obszar: Przedział ufności (95% CI).\n"
                       "• r (korelacja): Mówi, jak mocno stężenie wpływa na wynik (blisko 1.0 = idealna zależność).")

        self.add_entry("Szacowanie MIC (Minimalne Stężenie Hamujące)", 
                       "Program wyznacza teoretyczne MIC na podstawie punktu przecięcia linii trendu z osią średnicy krążka (6 mm). "
                       "Jest to model matematyczny logarytmiczno-liniowy. Wynik jest szacunkowy i służy do porównania siły substancji.")

        self.add_entry("Porównanie Międzygatunkowe", 
                       "Zestawienie działania wybranych substancji na wszystkie badane szczepy bakterii jednocześnie. "
                       "Pozwala ocenić spektrum działania (czy substancja działa na wszystko, czy tylko wybiórczo).")

        self.add_entry("Analiza PCA (Główne Składowe)", 
                       "Zaawansowana metoda wizualizacji podobieństwa między szczepami bakterii.\n"
                       "• Punkty blisko siebie: Szczepy o bardzo podobnym profilu wrażliwości na badane substancje.\n"
                       "• Punkty daleko od siebie: Szczepy reagujące odmiennie.\n"
                       "• Osie PC1 i PC2: Reprezentują główne kierunki zmienności w danych. Procent w nawiasie mówi, jak dużo informacji o różnicach widać na wykresie.")

        # --- SEKCJA 4: AUTOMATYCZNY OPIS METOD ---
        self.add_section("4. AUTOMATYCZNY OPIS (Materials and Methods)")
        
        # Logika dynamicznego tekstu
        generated_text = ""
        if parent.export_stats_main:
            test_info = parent.export_stats_main[0]
            used_test = test_info.get("Test", "")
            used_correction_raw = parent.combo_method.get()
            
            corr_map = {
                "holm": "Holm-Bonferroni correction",
                "bonferroni": "Bonferroni correction",
                "fdr_bh": "Benjamini-Hochberg (FDR) procedure",
                "None": "no correction"
            }
            correction_desc = corr_map.get(used_correction_raw, used_correction_raw)

            if "ANOVA" in used_test:
                generated_text = (
                    f"\"Statistical analysis was performed using Python (scipy, statsmodels). "
                    f"Normality was confirmed using the Shapiro-Wilk test. "
                    f"Differences between groups were analyzed using one-way ANOVA, followed by Tukey's HSD post-hoc test for multiple comparisons. "
                    f"Effect sizes were calculated using Cohen’s d estimator. "
                    f"A p-value < 0.05 was considered statistically significant.\""
                )
            elif "Kruskal" in used_test:
                generated_text = (
                    f"\"Statistical analysis was performed using Python (scipy, scikit-posthocs). "
                    f"Due to the non-normal distribution of data (Shapiro-Wilk test, p < 0.05), "
                    f"differences between groups were analyzed using the Kruskal-Wallis test. "
                    f"Pairwise comparisons were performed using Dunn's post-hoc test with {correction_desc}. "
                    f"Effect sizes were estimated using Cohen’s d. "
                    f"A p-value < 0.05 was considered statistically significant.\""
                )
            else:
                generated_text = "Analysis performed, but test type unrecognized."
            info_label = "Poniższy tekst został wygenerowany na podstawie Twoich OSTATNICH WYNIKÓW:"
        else:
            generated_text = (
                "\"Statistical analysis was performed using Python. "
                "Normality was assessed using the Shapiro-Wilk test. "
                "Differences between groups were analyzed using one-way ANOVA (for normal data) "
                "or Kruskal-Wallis test (for non-normal data), followed by post-hoc analysis "
                "with appropriate correction for multiple comparisons. "
                "Effect sizes were calculated using Cohen’s d estimator.\""
            )
            info_label = "To jest ogólny szablon. Uruchom analizę, aby uzyskać tekst dopasowany do Twoich danych."

        self.add_text(info_label)
        
        textbox = ctk.CTkTextbox(self.scroll, height=120, font=("Consolas", 11))
        textbox.pack(fill="x", pady=5)
        textbox.insert("0.0", generated_text)
        textbox.configure(state="disabled")

        # --- SEKCJA 5: UKRYTE MECHANIZMY ---
        self.add_section("5. WAŻNE UWAGI TECHNICZNE")

        self.add_entry("Wykres Trendu - Nazewnictwo", 
                       "Aby wykres trendu (dawka-odpowiedź) zadziałał, nazwa grupy w Excelu MUSI zawierać liczbę i jednostkę.\n"
                       "• Poprawnie: 'Ekstrakt (50 mg/ml)', 'Ojek (10%)', 'Próbka 0.5 ug/ml'.\n"
                       "• Źle: 'Próbka A', 'Stężenie wysokie'.\n"
                       "Bez tego program nie rozpozna osi X.")

        self.add_entry("Wykrywanie Outlierów (Dixon)", 
                       "Test Dixona uruchamia się automatycznie, ale tylko dla prób o liczebności N = 3 do 10. "
                       "Dla bardzo dużych prób (N > 10) test nie jest wykonywany, aby uniknąć błędów statystycznych.")

        self.add_entry("Minimalna liczebność próby", 
                       "Grupy posiadające mniej niż 2 wyniki (N < 2) są automatycznie pomijane w analizie statystycznej, "
                       "ponieważ niemożliwe jest obliczenie dla nich odchylenia standardowego.")

        self.add_entry("Inteligentne Sortowanie", 
                       "Program stosuje tzw. 'Natural Sort Order'. Oznacza to, że grupy 'Próbka 2' i 'Próbka 10' "
                       "ułożą się w kolejności 2 -> 10, a nie 10 -> 2 (jak w zwykłym sortowaniu alfabetycznym).")

        self.add_entry("Autokorekta Nazw", 
                       "Program automatycznie usuwa zbędne spacje z nazw w Excelu (np. zamienia 'E. coli ' na 'E. coli'). "
                       "Dzięki temu błędy typu 'spacja na końcu' nie są traktowane jako osobne grupy.")

        # Przycisk zamknięcia
        ctk.CTkButton(self.scroll, text="Zamknij Pomoc", fg_color="#333333", hover_color="#555555", command=self.destroy).pack(pady=30)

    def add_section(self, title):
        ctk.CTkLabel(self.scroll, text=title, anchor="w", 
                     font=ctk.CTkFont(size=16, weight="bold"), text_color="#1F6AA5").pack(fill="x", pady=(25, 10))
        ctk.CTkFrame(self.scroll, height=2, fg_color="gray").pack(fill="x", pady=(0, 10))

    def add_text(self, text):
        # Kolor tekstu: ("czarny", "jasny szary") - poprawiona czytelność
        ctk.CTkLabel(self.scroll, text=text, font=ctk.CTkFont(size=13), 
                     wraplength=650, justify="left", anchor="w", text_color=("black", "gray80")).pack(fill="x", pady=(0, 10))

    def add_entry(self, subtitle, description):
        frame = ctk.CTkFrame(self.scroll, fg_color="transparent")
        frame.pack(fill="x", pady=8)
        
        # Tytuł wpisu - auto-kolor (czarny w light, biały w dark)
        ctk.CTkLabel(frame, text=f"• {subtitle}", font=ctk.CTkFont(size=13, weight="bold"), anchor="w", text_color=("black", "#E0E0E0")).pack(fill="x")
        
        # Treść wpisu - czytelny szary w obu trybach
        content_frame = ctk.CTkFrame(frame, fg_color="transparent")
        content_frame.pack(fill="x", padx=(15, 0))
        ctk.CTkLabel(content_frame, text=description, font=ctk.CTkFont(size=12), wraplength=650, justify="left", anchor="w", text_color=("gray30", "gray80")).pack(fill="x")


# ======================================================
# OKNO O TWÓRCY (ABOUT)
# ======================================================
class AboutDialog(ctk.CTkToplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("O twórcy")
        self.geometry("400x300")
        
        # Okno na wierzchu
        self.lift()
        self.attributes("-topmost", True)
        
        # Layout
        self.grid_columnconfigure(0, weight=1)
        
        ctk.CTkLabel(self, text="BioStat Master v3.0", font=ctk.CTkFont(size=20, weight="bold")).pack(pady=(30, 10))
        
        ctk.CTkLabel(self, text="Modular Analysis Tool", font=ctk.CTkFont(size=12, slant="italic")).pack(pady=(0, 20))
        
        ctk.CTkLabel(self, text="Created by:", font=ctk.CTkFont(size=14)).pack(pady=(10, 5))
        ctk.CTkLabel(self, text="StarGate3", font=ctk.CTkFont(size=16, weight="bold"), text_color="#1F6AA5").pack(pady=0)
        
        # Link (jako readonly entry dla łatwego kopiowania lub label)
        link = "https://github.com/StarGate3"
        self.link_entry = ctk.CTkEntry(self, width=250, justify="center", fg_color="transparent", border_width=0)
        self.link_entry.insert(0, link)
        self.link_entry.configure(state="readonly") # Aby można było chociaż zaznaczyć i skopiować
        self.link_entry.pack(pady=10)
        
        ctk.CTkButton(self, text="Zamknij", command=self.destroy).pack(pady=20)
