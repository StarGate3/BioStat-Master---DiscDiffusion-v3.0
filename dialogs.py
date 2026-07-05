import customtkinter as ctk
from config import ALPHA

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

        ctk.CTkLabel(
            self,
            text=(
                "Uwaga: przy n=3 test Dixona jest bardzo czuły - punkt różniący się "
                "choćby o 1 mm od pozostałych dwóch często zostaje wykryty, mimo że to "
                "naturalna zmienność biologiczna, a nie błąd pomiaru. Nic nie jest "
                "domyślnie zaznaczone do usunięcia - decyzję podejmij świadomie."
            ),
            font=ctk.CTkFont(size=11, slant="italic"), text_color="gray", wraplength=450, justify="left",
        ).pack(pady=(0, 10), padx=10)

        self.scroll = ctk.CTkScrollableFrame(self, width=450, height=250)
        self.scroll.pack(pady=5, padx=10, fill="both", expand=True)

        self.check_vars = {}

        for item in outlier_data:
            group = item['group']
            value = item['value']
            others = item['others']
            desc = f"Grupa: {group}\nOdstający: {value} mm (Pozostałe: {others})"

            # Domyślnie ODZNACZONE: usunięcie punktu ma być świadomą decyzją
            # użytkownika (opt-in), nie automatycznym domyślnym zachowaniem.
            var = ctk.IntVar(value=0)
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
                       f"Sprawdzamy, czy dane w każdej grupie układają się w 'krzywą dzwonową'.\n"
                       f"• p > {ALPHA}: Rozkład normalny (OK).\n"
                       f"• p < {ALPHA}: Rozkład inny niż normalny (częste w małych próbach n=3).")
        
        self.add_entry("Krok 2: Wariancja (Levene)", 
                       "Sprawdzamy, czy grupy mają podobny 'rozrzut' wyników.\n"
                       "Jeśli wariancje są różne, testy parametryczne mogą dawać błędne wyniki.")
        
        self.add_entry("Krok 3: Wybór Testu Głównego",
                       "• ANOVA: Wybierana, gdy dane są normalne i mają równą wariancję (największa moc).\n"
                       "• Kruskal-Wallis: Wybierany, gdy założenia ANOVA nie są spełnione (bezpieczniejszy dla danych mikrobiologicznych).")

        self.add_entry("Uwaga praktyczna: dane bez replikacji biologicznej",
                       "Normalność (Krok 1) sprawdzana jest na powtórzeniach BIOLOGICZNYCH, nie na surowych "
                       "pomiarach technicznych. Jeśli plik nie rozróżnia powtórzeń biologicznych/technicznych "
                       "(stary, jednoarkuszowy format - każda grupa ma wtedy dokładnie jedno powtórzenie "
                       "biologiczne, n_bio=1), testu Shapiro-Wilka nie da się w ogóle wykonać na pojedynczej "
                       "wartości, więc program w praktyce zawsze przechodzi do testu nieparametrycznego "
                       "(Kruskal-Wallis) - ANOVA w takim przypadku nigdy nie zostanie wybrana, niezależnie od "
                       "tego, jak bardzo 'normalne' były surowe dane techniczne. ANOVA staje się osiągalna "
                       "dopiero przy realnych powtórzeniach biologicznych (kolumna Rep_biologiczna, nowy format).")

        # --- SEKCJA 2: KOREKTY POST-HOC ---
        self.add_section("2. KOREKTY POST-HOC (Którą wybrać?)")
        self.add_text("Gdy porównujesz wiele próbek naraz, rośnie ryzyko, że przypadkowo znajdziesz 'istotną' różnicę (Błąd I rodzaju). Korekty temu zapobiegają.")
        
        self.add_entry("Holm (Zalecana)", 
                       "Najlepszy balans. Jest silniejsza niż brak korekty, ale nie tak 'brutalna' jak Bonferroni. "
                       "Dobra do większości standardowych badań.")
        
        self.add_entry("Bonferroni",
                       f"Bardzo konserwatywna. Bardzo trudno uzyskać p < {ALPHA}. "
                       f"Stosuj tylko, gdy musisz mieć absolutną pewność i chcesz uniknąć fałszywych alarmów za wszelką cenę.")
        
        self.add_entry("FDR (Benjamini-Hochberg)",
                       "Najmniej rygorystyczna. Dopuszcza pewien odsetek fałszywych odkryć. "
                       "Idealna do 'screeningu' (przesiewu) setek substancji, gdy nie chcesz przegapić niczego potencjalnie ciekawego.")

        self.add_entry("Brak korekty (None)",
                       "Surowe p-value z każdego porównania, bez żadnej poprawki na wielokrotne testowanie. "
                       "RYZYKO: przy wielu porównaniach naraz mocno zawyża odsetek fałszywie 'istotnych' wyników "
                       "(Błąd I rodzaju) - to dokładnie to, przed czym mają chronić korekty opisane wyżej. "
                       "Sensowne najwyżej przy pojedynczym, z góry zaplanowanym porównaniu (np. tylko badana "
                       "substancja vs kontrola, bez żadnych innych par) - nie do rutynowego używania przy "
                       "porównaniu wielu grup.")

        # --- SEKCJA 3: INTERPRETACJA WYKRESÓW ---
        self.add_section("3. ANATOMIA WYKRESÓW")

        self.add_entry("Wykres Słupkowy (Barplot)", 
                       "• Wysokość słupka: Średnia arytmetyczna strefy zahamowania.\n"
                       "• Antenka (Słupek błędu): Pokazuje zmienność (SD - Odchylenie Standardowe). Im krótsza antenka, tym bardziej powtarzalne były wyniki.")

        self.add_entry("Wykres Pudełkowy (Boxplot)",
                       "Bardziej szczegółowy niż słupkowy - POD WARUNKIEM, że grupa ma więcej niż jedno "
                       "powtórzenie biologiczne (n_bio≥2):\n"
                       "• Linia w środku pudełka: Mediana (wartość środkowa).\n"
                       "• Pudełko: Obejmuje 50% środkowych wyników (od 25. do 75. percentyla).\n"
                       "• Wąsy: Zasięg danych (min-max), z wyłączeniem wartości odstających.\n"
                       "Przy n_bio=1 (np. stary, jednoarkuszowy format - patrz 'Minimalna liczebność próby' "
                       "niżej) grupa ma tylko JEDNĄ wartość biologiczną, więc pudełko degeneruje się do "
                       "pojedynczego punktu (brak IQR, brak wąsów do pokazania) - to nie błąd wykresu, tylko "
                       "wierne odwzorowanie braku replikacji biologicznej w danych źródłowych.")

        self.add_entry("Wykres 'Lollipop' (Wielkość Efektu)",
                       f"Najważniejszy wykres do oceny 'siły' działania.\n"
                       f"• Oś pozioma (Cohen's d): Mówi, ile 'odchyleń standardowych' dzieli dwie grupy.\n"
                       f"• Kropka ZIELONA (W prawo): Grupa badana jest lepsza/silniejsza.\n"
                       f"• Kropka CZERWONA (W lewo): Grupa badana jest gorsza/słabsza.\n"
                       f"UWAGA: Domyślnie wykres pokazuje wyłącznie istotne statystycznie (p < {ALPHA}) pary "
                       f"WZGLĘDEM GRUPY ODNIESIENIA (spójnie z wykresem głównym), a nie wszystkie istotne pary "
                       f"między wszystkimi grupami - przy wielu grupach lista wszystkich par bywa zbyt długa, "
                       f"żeby dało się ją czytelnie podpisać. Pełny widok 'wszystkie pary' włączysz "
                       f"przełącznikiem 'Wielkość efektu: wszystkie pary' w Opcjach Wykresu. Pełne wyniki dla "
                       f"wszystkich par (również nieistotnych) niezależnie od tego przełącznika znajdziesz w "
                       f"raporcie Excel (zakładka 'Post-hoc (Details)').")

        self.add_entry("Mapa Ciepła (Heatmap)", 
                       "Wizualizacja macierzy. Kolory ułatwiają szybkie wyłapanie liderów.\n"
                       "• Jasne/Ciepłe pola: Duże strefy zahamowania (aktywność).\n"
                       "• Ciemne pola: Brak aktywności.")
        
        self.add_entry("Trend (Dawka-Odpowiedź)", 
                       "Jeśli w nazwach grup są stężenia (np. 5mg, 10mg), ten wykres pokaże linię.\n"
                       "• Zacieniony obszar: Przedział ufności (95% CI).\n"
                       "• r (korelacja): Mówi, jak mocno stężenie wpływa na wynik (blisko 1.0 = idealna zależność).")

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
                    f"A p-value < {ALPHA} was considered statistically significant.\""
                )
            elif "Kruskal" in used_test:
                generated_text = (
                    f"\"Statistical analysis was performed using Python (scipy, scikit-posthocs). "
                    f"Due to the non-normal distribution of data (Shapiro-Wilk test, p < {ALPHA}), "
                    f"differences between groups were analyzed using the Kruskal-Wallis test. "
                    f"Pairwise comparisons were performed using Dunn's post-hoc test with {correction_desc}. "
                    f"Effect sizes were estimated using Cohen’s d. "
                    f"A p-value < {ALPHA} was considered statistically significant.\""
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
                       "Z analizy statystycznej pomijane są wyłącznie grupy CAŁKOWICIE puste (N=0). "
                       "Grupa z jednym wynikiem (N=1 - np. brak replikacji biologicznej) NIE jest pomijana: "
                       "test wciąż się wykonuje, ale taki wynik jest oznaczony jako orientacyjny, z jawnym "
                       "ostrzeżeniem 'brak replikacji biologicznej' (widocznym m.in. na wykresie głównym) - "
                       "bo dla N=1 nie da się policzyć odchylenia standardowego ani ocenić normalności rozkładu.")

        self.add_entry("Inteligentne Sortowanie", 
                       "Program stosuje tzw. 'Natural Sort Order'. Oznacza to, że grupy 'Próbka 2' i 'Próbka 10' "
                       "ułożą się w kolejności 2 -> 10, a nie 10 -> 2 (jak w zwykłym sortowaniu alfabetycznym).")

        self.add_entry("Autokorekta Nazw",
                       "Program automatycznie usuwa zbędne spacje z nazw w Excelu (np. zamienia 'E. coli ' na 'E. coli'). "
                       "Dzięki temu błędy typu 'spacja na końcu' nie są traktowane jako osobne grupy.")

        # --- SEKCJA 6: MODUŁ MIC/MBC ---
        self.add_section("6. MODUŁ MIC/MBC")

        self.add_entry("Czym jest MIC i MBC",
                       "MIC (Minimalne Stężenie Hamujące) to najniższe testowane stężenie substancji, przy "
                       "którym NIE obserwuje się wzrostu bakterii. MBC (Minimalne Stężenie Bójcze) to "
                       "najniższe stężenie, przy którym bakterie zostają faktycznie ZABITE (nie tylko "
                       "zahamowane) - zawsze wyższe lub równe MIC dla tej samej substancji/szczepu.")

        self.add_entry("Jak program wyznacza MIC",
                       "• Odczyt wizualny (arkusz MIC_wizualny): każda studzienka to 'wzrost' albo 'brak' - "
                       "wpisane wprost przez osobę wykonującą test.\n"
                       "• Odczyt OD (arkusz MIC_OD): każda studzienka to zmierzona gęstość optyczna (OD). "
                       "Program przelicza ją na 'wzrost'/'brak' progiem WZGLĘDNYM: procent wzrostu = "
                       "(OD studzienki - OD kontroli jałowości) / (OD kontroli wzrostu - OD kontroli "
                       "jałowości); poniżej progu (domyślnie 10%) uznaje się, że wzrostu nie ma. Próg jest "
                       "rozsądną, konfigurowalną wartością inżynierską, nie liczbą z konkretnego standardu "
                       "CLSI/EUCAST.\n"
                       "MIC to najniższe testowane stężenie, przy którym seria (skanowana od najwyższego "
                       "stężenia w dół) po raz pierwszy pokazuje 'brak wzrostu' i pozostaje przy tym "
                       "konsekwentnie aż do najniższego testowanego stężenia; przy zaburzeniu tej kolejności "
                       "('skip well') wynik jest zwracany zachowawczo, ale oznaczony jako wymagający ręcznej "
                       "weryfikacji. Wynik jest też uzależniony od tego, czy Przebieg w ogóle przeszedł "
                       "walidację kontroli (patrz niżej).")

        self.add_entry("Jak program wyznacza MBC",
                       "Arkusz MBC_posiew: każda studzienka to liczba kolonii (CFU) po posiewie na czyste "
                       "podłoże. Studzienkę uznaje się za 'zabójczą', gdy redukcja względem inokulum "
                       "wyjściowego (Inokulum_CFU_t0) wynosi co najmniej 99,9% (3-log10) - to faktycznie "
                       "standardowa, klinicznie przyjęta definicja działania bakteriobójczego, nie autorski "
                       "wybór. MBC to najniższe testowane stężenie spełniające ten warunek, wyznaczane tą "
                       "samą metodą skanowania serii co MIC.")

        self.add_entry("Cenzura (wartości '≤'/'>')",
                       "Jeśli NAJNIŻSZE testowane stężenie już daje 'brak'/'zabójcze', prawdziwe MIC/MBC "
                       "może być jeszcze niższe - program zwraca to jako wartość CENZUROWANĄ DOLNIE "
                       "('≤ najniższe stężenie'), nigdy jako zwykłą liczbę. Analogicznie, jeśli NAJWYŻSZE "
                       "testowane stężenie wciąż nie daje 'brak'/'zabójcze', MIC/MBC jest cenzurowane "
                       "GÓRNIE ('> najwyższe stężenie') - prawdziwa wartość może być wyższa niż to, co "
                       "przetestowano. Cenzura jest przenoszona przez wszystkie kolejne etapy (agregację do "
                       "powtórzenia biologicznego, medianę grupową, iloraz MBC/MIC) - nigdy nie jest po "
                       "cichu zamieniana na zwykłą liczbę.")

        self.add_entry("Rola kontroli (arkusz Kontrole)",
                       "Każdy Przebieg musi mieć dokładnie jeden wpis w arkuszu Kontrole. Dwa warunki "
                       "ważności przebiegu sprawdzane są NIEZALEŻNIE i identycznie dla obu trybów odczytu "
                       "MIC (wizualnego i OD) - różni się tylko FORMA samej kontroli:\n"
                       "• Kontrola wzrostu: potwierdza, że bakteria w ogóle urosła w tym przebiegu (bez "
                       "substancji). W trybie OD to warunek liczbowy (różnica OD kontroli wzrostu i "
                       "jałowości musi przekroczyć próg); w trybie wizualnym to po prostu słowo 'wzrost' w "
                       "tej kolumnie - 'brak' odrzuca cały przebieg jako nieważny.\n"
                       "• Kontrola jałowości: potwierdza, że samo podłoże (bez inokulum) jest czyste. W "
                       "trybie OD to próg liczbowy BEZWZGLĘDNY (niezależny od poziomu kontroli wzrostu); w "
                       "trybie wizualnym - słowo 'brak' ('wzrost' oznacza skażenie i odrzuca przebieg). "
                       "Wysoka kontrola wzrostu NIE maskuje skażonej kontroli jałowości - oba warunki muszą "
                       "być spełnione niezależnie, w obu trybach.\n"
                       "• Inokulum (CFU w chwili t0): punkt odniesienia do liczenia % redukcji CFU dla MBC "
                       "- używane tylko przy odczycie MBC, nie przy MIC.\n"
                       "Przebieg BEZ pasującego wpisu w Kontrole (albo z więcej niż jednym wpisem dla tego "
                       "samego Przebiegu - to błąd danych, nie zgadywane) jest odrzucany z jawnym powodem, "
                       "nigdy po cichu pomijany.")

        self.add_entry("Iloraz MBC/MIC i klasyfikacja bakteriobójcze/bakteriostatyczne",
                       "Iloraz liczony jest jako różnica LICZBY dwukrotnych rozcieńczeń (kroków log2) "
                       "między MBC a MIC, nie jako proste dzielenie stężeń w jednostkach - dzięki temu "
                       "wartości cenzurowane (≤/>) też dają sensowny, choć czasem tylko częściowo "
                       "rozstrzygalny wynik. Klasyfikacja jest podawana WYŁĄCZNIE dla serii dwukrotnych "
                       "rozcieńczeń (współczynnik rozcieńczenia = 2):\n"
                       "• MBC/MIC ≤ 4 (różnica ≤2 rozcieńczeń) → bakteriobójcze.\n"
                       "• MBC/MIC ≥ 8 (różnica ≥3 rozcieńczeń) → bakteriostatyczne.\n"
                       "• Gdy cenzura nie pozwala jednoznacznie rozstrzygnąć, po której stronie progu leży "
                       "prawdziwa wartość → 'nieoznaczalny'.\n"
                       "Dla dowolnego innego współczynnika rozcieńczenia MIC, MBC i sam iloraz są liczone i "
                       "pokazywane jak zawsze, ale etykieta klasyfikacji to jawne 'niedostępna' z podanym "
                       "powodem - program nigdy nie zgaduje kategorii, dla której nie ma metodologicznych "
                       "podstaw. Wynik MBC niższy niż MIC dla tej samej pary jest zgłaszany jako BŁĄD "
                       "SPÓJNOŚCI danych (nie jako 'nieoznaczalny') i wyłączony z klasyfikacji.")

        self.add_entry("Powtórzenia i ostrzeżenie n_bio=1",
                       "Tak jak w module dyfuzji, wynik z wielu powtórzeń TECHNICZNYCH tego samego "
                       "powtórzenia biologicznego jest agregowany regułą 'wysokiej mediany' do JEDNEJ "
                       "wartości na powtórzenie biologiczne, zanim cokolwiek zostanie porównane między "
                       "grupami. Grupa oparta na tylko jednym powtórzeniu biologicznym (n_bio=1) nie jest "
                       "odrzucana, ale każdy wynik z niej jest oznaczony jako orientacyjny - ostrzeżenie "
                       "widoczne wprost na wykresach, w tabeli zbiorczej i w eksporcie, nie tylko w logu.")

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
