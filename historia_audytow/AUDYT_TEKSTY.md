# AUDYT_TEKSTY.md — recenzja aktualności i poprawności merytorycznej dwóch generatorów tekstu

**Zakres:** zadanie wyłącznie diagnostyczne, wyłącznie odczytujące. Nic w repozytorium nie zostało zmienione ani zacommitowane. Praca wykonana na `main` (potwierdzone czystym `git status` przed i po).

**Metoda:** dla obu funkcji przeczytano pełny, aktualny kod źródłowy (nie interpretację/pamięć z wcześniejszych sesji), a następnie **uruchomiono obie funkcje na żywej instancji `App()`** (bez pętli mainloop, ten sam wzorzec testowy używany w całym tym projekcie) w kilku scenariuszach — łącznie z celowo skonstruowanymi danymi syntetycznymi, żeby wymusić konkretne ścieżki kodu (ANOVA vs Kruskal-Wallis, referencja = kontrola dodatnia, brak analizy) i zobaczyć DOSŁOWNY wygenerowany tekst, a nie tylko przeczytać kod.

---

## A. "Generuj Opisy Rycin" — lokalizacja i pełna treść

**Plik:** `gui.py`, metoda `open_caption_window` (linie 418–485), wywoływana przez przycisk `📝 Generuj Opisy Rycin` (linia 150).

Funkcja generuje **statyczny szablon 6 podpisów rycin** (Rycina 1–6, odpowiadające 6 zakładkom głównego okna: Wykres Główny, Mapa Ciepła, Wielkość Efektu, Mapa P-value, Trend, Porównanie Szczepów), z podstawieniem kilku zmiennych (`bact`, `ref_group`, `post_hoc`, `test_name`, `err_desc`, `viz_desc`, `ALPHA`, `DISC_DIAMETER_MM`). Pełna treść szablonu (fragment kodu, linie 452–484):

```python
captions = f"""--- OPISY RYCIN (Scientific Captions) ---\n
Możesz skopiować poniższe opisy bezpośrednio do manuskryptu (Word/LaTeX).

=== Rycina 1: Wykres Główny ===
Figure 1. Antibacterial activity of tested samples against {bact}.
{viz_desc}. Error bars indicate the {err_desc} of independent replicates.
Statistical significance was determined using {test_name} followed by {post_hoc} for multiple comparisons.
Asterisks (*) indicate a statistically significant difference (p < {ALPHA}) compared to the negative control ({ref_group}).
Red dashed line represents the diameter of the disk ({DISC_DIAMETER_MM:g} mm).

=== Rycina 2: Mapa Ciepła ===
Figure 2. Heatmap visualizing the magnitude of growth inhibition zones (mm) for {bact} treated with various substances.
Color intensity corresponds to the mean diameter of the inhibition zone. Warmer colors indicate higher antibacterial activity.

=== Rycina 3: Mapa Wielkości Efektu (Effect Size) ===
Figure 3. Lollipop chart displaying the standardized effect size (Cohen's d) for statistically significant pairwise comparisons.
Dots represent the magnitude of the difference between groups. Green dots indicate a positive difference (Group 1 > Group 2), while red dots indicate a negative difference.

=== Rycina 4: Mapa Istotności (P-value Matrix) ===
Figure 4. Pairwise comparison significance matrix (P-values).
The heatmap displays adjusted p-values for all pairwise comparisons. Blue shades indicate statistical significance (p < {ALPHA}), while red/white shades indicate non-significant differences.
P-values were adjusted for multiple comparisons using the {post_hoc} method.

=== Rycina 5: Trend Dawka-Odpowiedź ===
Figure 5. Dose-response relationship of antibacterial activity.
Lines represent the trend of inhibition zone diameter (mm) across increasing concentrations of tested substances.
Shaded areas indicate the confidence interval. Spearman correlation coefficients (r) are provided for each substance to quantify the strength of the monotonic relationship.

=== Rycina 6: Porównanie Szczepów ===
Figure 6. Cross-species comparison of antibacterial activity.
Bar chart summarizing the mean inhibition zone diameters for selected substances across different bacterial strains.
Error bars represent standard deviation. This overview highlights the differential susceptibility of tested pathogens to the antimicrobial agents.
"""
```

Brak jakiejkolwiek wzmianki o MIC/MBC — funkcja generuje wyłącznie te 6 podpisów, niezależnie od tego, czy wczytany plik ma dane MIC/MBC, czy nie. Przycisk `btn_captions` nie jest nigdy blokowany/odblokowywany warunkowo — jest zawsze klikalny, niezależnie od tego, czy plik jest wczytany, czy analiza została uruchomiona.

---

## B. "Podręcznik / Pomoc" — lokalizacja i pełna treść

**Plik:** `dialogs.py`, klasa `HelpDialog` (linie 71–272), otwierana przez `open_help_window` (`gui.py:408`) z przycisku `❓ Podręcznik / Pomoc`.

Statyczny tekst w 5 sekcjach + jeden dynamicznie generowany akapit (Sekcja 4). Pełna treść (parafraza 1:1 wywołań `add_section`/`add_text`/`add_entry`, patrz kod dla dosłownego brzmienia):

**1. ALGORYTM DZIAŁANIA PROGRAMU**
- Krok 1: Normalność (Shapiro-Wilk) — "p > α: rozkład normalny; p < α: rozkład inny niż normalny (częste w małych próbach n=3)".
- Krok 2: Wariancja (Levene).
- Krok 3: Wybór testu — ANOVA (normalne + równa wariancja) vs Kruskal-Wallis (w przeciwnym razie).

**2. KOREKTY POST-HOC** — Holm (Zalecana), Bonferroni (bardzo konserwatywna), FDR/Benjamini-Hochberg (najmniej rygorystyczna, do przesiewu).

**3. ANATOMIA WYKRESÓW** — Barplot, Boxplot, Lollipop (Wielkość Efektu — "prezentuje wyłącznie pary różniące się istotnie statystycznie (p < α)... Pełne wyniki... w raporcie Excel (zakładka 'Post-hoc Details')"), Heatmap, Trend Dawka-Odpowiedź (wymaga liczby+jednostki w nazwie grupy), Porównanie Międzygatunkowe, PCA. **Brak jakiejkolwiek sekcji o wykresach MIC/MBC.**

**4. AUTOMATYCZNY OPIS (Materials and Methods)** — generowany dynamicznie z `parent.export_stats_main`/`parent.combo_method`:
```python
if "ANOVA" in used_test:
    generated_text = ("... one-way ANOVA, followed by Tukey's HSD post-hoc test ...")
elif "Kruskal" in used_test:
    generated_text = ("... Kruskal-Wallis test. Pairwise comparisons were performed using "
                       "Dunn's post-hoc test with {correction_desc} ...")
```

**5. WAŻNE UWAGI TECHNICZNE** — nazewnictwo dla wykresu trendu, zakres Dixona (N=3–10), "Grupy posiadające mniej niż 2 wyniki (N<2) są automatycznie pomijane", naturalne sortowanie, autokorekta spacji w nazwach. **Brak jakiejkolwiek wzmianki o nowym formacie wieloarkuszowym, kolumnach Typ/Rep_biologiczna/Rep_techniczna, arkuszu Kontrole, czy o module MIC/MBC.**

---

## Znaleziska (uszeregowane wg powagi)

### KRYTYCZNE

---

**Znalezisko 1 — Rycina 1: etykieta "kontrola negatywna" jest bezwarunkowa, nawet gdy referencja to kontrola dodatnia**

(a) Cytat: `Asterisks (*) indicate a statistically significant difference (p < {ALPHA}) compared to the negative control ({ref_group}).`

(b) Co robi kod naprawdę: `combo_ref` (lista wyboru grupy referencyjnej) jest wypełniana **wszystkimi** grupami danego szczepu (`gui.py:377`, `self.combo_ref.configure(values=grupy_bact)`), nie tylko autodetekowaną kontrolą negatywną. Użytkownik może świadomie wybrać dowolną grupę jako referencję (np. kontrolę dodatnią, żeby porównać wszystko względem znanego antybiotyku). **Zweryfikowano na żywo** na `dane_disk.xlsx`: po ręcznym ustawieniu `combo_ref` na `"Kontrola (+) Ampycylina"` i uruchomieniu analizy, wygenerowany podpis brzmiał dosłownie:
```
Asterisks (*) indicate a statistically significant difference (p < 0.05) compared to the negative control (Kontrola (+) Ampycylina).
```

(c) Rozbieżność: tekst nazywa dowolnie wybraną grupę referencyjną "kontrolą negatywną" bez sprawdzenia, czy nią rzeczywiście jest. Zdanie jest wewnętrznie sprzeczne ("negative control (Kontrola (+)...)").

(d) Powaga: **KRYTYCZNE** — zdanie wprost sprzeczne ze sobą, trafiające dosłownie do manuskryptu.

---

**Znalezisko 2 — Rycina 1/4: nazwa korekty post-hoc ignoruje to, że ścieżka ANOVA zawsze używa Tukeya, niezależnie od wyboru użytkownika**

(a) Cytat: `Statistical significance was determined using {test_name} followed by {post_hoc} for multiple comparisons.` (gdzie `post_hoc` pochodzi wyłącznie z `self.combo_method.get()`, zmapowanego 1:1 na opis tekstowy, niezależnie od `test_name`).

(b) Co robi kod naprawdę: w `logic.py::run_statistics`, gałąź ANOVA **zawsze** wywołuje `pairwise_tukeyhsd(...)` (Tukey's HSD) — parametr `method` (wybór Holm/Bonferroni/FDR z paska bocznego) jest używany **wyłącznie** w gałęzi Kruskal-Wallis (`sp.posthoc_dunn(..., p_adjust=method)`). Wybór użytkownika w `combo_method` nie ma żadnego wpływu, gdy dane trafią na ścieżkę ANOVA. **Zweryfikowano na żywo**: zbudowano syntetyczny plik w nowym formacie z realnymi powtórzeniami biologicznymi (n_bio=6/grupę, Shapiro p=0.6–0.95, Levene p=0.77), wymuszając rzeczywiste użycie ANOVA (`test_used: ANOVA`, `Statistic=364.77`). Z `combo_method="holm"` wygenerowany podpis brzmiał:
```
Statistical significance was determined using One-way ANOVA followed by Holm-Bonferroni correction for multiple comparisons.
```
— podczas gdy kod faktycznie wykonał Tukey's HSD, nie korektę Holma. Dla kontrastu: **generator opisu metod w Podręczniku (Sekcja 4, dialogs.py) robi to poprawnie** — jawnie sprawdza `if "ANOVA" in used_test` i wtedy zawsze pisze "Tukey's HSD", ignorując wybór z dropdowna tylko dla tej gałęzi. Poprawna logika istnieje więc w kodzie, ale nie została przeniesiona do generatora podpisów rycin.

(c) Rozbieżność: podpis podaje fałszywą nazwę metody korekty (Holm-Bonferroni zamiast Tukey's HSD), ilekroć aktywna jest ścieżka ANOVA i użytkownik ma ustawioną dowolną inną korektę niż domyślna.

(d) Powaga: **KRYTYCZNE** — błędne stwierdzenie metodologiczne trafiające wprost do sekcji Methods/podpisu.

---

**Znalezisko 3 — Rycina 3 (Wielkość Efektu): podpis generowany bezwarunkowo, mimo że wykres jest PUSTY dla całej klasy plików**

(a) Cytat: `Figure 3. Lollipop chart displaying the standardized effect size (Cohen's d) for statistically significant pairwise comparisons.`

(b) Co robi kod naprawdę: `calculate_cohens_d` zwraca `NaN`, gdy którakolwiek z grup ma n<2 (`utils.py`: `if n1 < 2 or n2 < 2: return np.nan`). `draw_effect_plot` odfiltrowuje wszystkie wpisy z `NaN` Cohen's d i zwraca `None`, jeśli nic nie zostaje. **Zweryfikowano na żywo** na `dane_disk.xlsx` (plik w starym formacie, flagowy plik testowy tego projektu): `app.figures['effect'] is None` → **True**. Powód: patrz Znalezisko 6 — każda grupa w starym formacie kolapsuje do n_bio=1, więc Cohen's d jest zawsze `NaN` dla każdej pary.

(c) Rozbieżność: funkcja generuje opis ryciny, która w praktyce (dla każdego pliku w starym formacie, i dla nowego formatu z n_bio<2) nie istnieje — zakładka jest pusta, a podpis i tak zostaje wygenerowany, opisując "lollipop chart" którego nie ma.

(d) Powaga: **KRYTYCZNE** — dla bardzo dużej, praktycznej klasy plików podpis opisuje nieistniejącą rycinę.

---

**Znalezisko 4 — Rycina 1 i 6: "error bars = SD niezależnych powtórzeń" jest fałszywe dla każdego pliku w starym formacie (słupki błędu = zero)**

(a) Cytat: `Error bars indicate the {err_desc} of independent replicates.` (Rycina 1) oraz `Error bars represent standard deviation.` (Rycina 6, Porównanie Szczepów).

(b) Co robi kod naprawdę: dane do obu wykresów (`df_bio` w `run_analysis`, linia 678; `df_bio_all` dla Rycina 6, linia 775) pochodzą z `utils.aggregate_technical_replicates`. Dla plików w **starym formacie** (brak kolumn `Rep_biologiczna`/`Rep_techniczna`) funkcja ta — celowo, zgodnie z własnym docstringiem — **kolapsuje każdą grupę do dokładnie JEDNEGO wiersza (n_bio=1)**, uśredniając wszystkie surowe powtórzenia techniczne. **Zweryfikowano na żywo** na `dane_disk.xlsx`/E. coli: wszystkie 30 grup mają `n_bio=1` (`df_bio.groupby(Grupa).size().unique() == [1]`), a co za tym idzie SD w danych zasilających wykres wynosi **dokładnie 0.0 dla każdej z 30 grup** (`std()` jednej liczby = NaN → `.fillna(0)` → 0). Słupek błędu na Rycinie 1 i 6 jest więc zerowej wysokości dla KAŻDEGO pliku w starym formacie — w tym `dane_disk.xlsx`, głównego pliku testowego całego tego projektu.

(c) Rozbieżność: (i) słupek błędu faktycznie nie pokazuje żadnej zmienności (wysokość 0) — sprzeczne z "error bars indicate SD"; (ii) nawet gdyby SD było niezerowe, "independent replicates" jest mylące — n_bio=1 oznacza JEDNO powtórzenie biologiczne (uśrednione technicznie), nie wiele niezależnych powtórzeń.

(d) Powaga: **KRYTYCZNE** — opis wprost przeciwny temu, co widać na wykresie, dla plików reprezentujących większość dotychczasowych danych użytkownika.

---

**Znalezisko 5 — Całkowity brak modułu MIC/MBC w obu tekstach**

(a) Cytat: cały blok "Generuj Opisy Rycin" (6 rycin) i cała Sekcja 3 Podręcznika ("ANATOMIA WYKRESÓW") — żadne wystąpienie słów "MIC", "MBC" w żadnym z dwóch plików źródłowych (potwierdzone `grep`).

(b) Co robi kod naprawdę: moduł MIC/MBC (od dawna wpięty do GUI) generuje 4 własne wykresy (`mic_distribution`, `mic_pairs`, `mic_comparison`, `mbc_comparison`) w oknie z zakładkami `Rozkład MIC/MBC`, `Pary MIC↔MBC`, `Porównanie substancji`, `Tabela zbiorcza i ostrzeżenia`, oraz cały aparat pojęciowy: iloraz MBC/MIC, klasyfikacja bakteriobójcze/bakteriostatyczne (dostępna tylko dla serii dwukrotnych, `Wsp_rozc=2` — naprawa audytu 1.4), cenzura (`≤`/`≥`), ostrzeżenia o braku replikacji biologicznej.

(c) Rozbieżność: generator podpisów nie tworzy ŻADNEGO podpisu dla żadnego z 4 wykresów MIC/MBC — użytkownik analizujący MIC/MBC nie dostaje żadnego wsparcia z narzędzia, którego głównym zadaniem jest właśnie to. Podręcznik nie tłumaczy użytkownikowi ani jednego z tych pojęć (iloraz, klasyfikacja, cenzura, `Wsp_rozc=2`), mimo że są one widoczne wprost w interfejsie.

(d) Powaga: **KRYTYCZNE** dla generatora podpisów (blokuje legalne wykorzystanie modułu do publikacji), **ISTOTNE–KRYTYCZNE** dla podręcznika (całe pasmo funkcjonalności bez żadnej dokumentacji w aplikacji).

---

### ISTOTNE

---

**Znalezisko 6 — "Krok 1/Krok 3": ANOVA jest de facto nieosiągalna dla każdego pliku w starym formacie, bez żadnej wzmianki o tym w tekście**

(a) Cytat: *"Krok 1: Normalność (Shapiro-Wilk)... p < α: Rozkład inny niż normalny (częste w małych próbach n=3)."* / *"Krok 3:... ANOVA: Wybierana, gdy dane są normalne i mają równą wariancję."*

(b) Co robi kod naprawdę: test statystyczny jest uruchamiany na `df_bio` (dane PO agregacji do powtórzeń biologicznych), nie na surowych wierszach. Shapiro-Wilk wymaga `len(vals) >= 3` (`logic.py:46`), inaczej `is_norm` pozostaje `False` domyślnie. Dla starego formatu każda grupa ma **zawsze** n_bio=1 (patrz Znalezisko 4) → Shapiro nigdy się nie wykonuje → `all_normal` jest zawsze `False` → **ANOVA jest nieosiągalna, zawsze wybierany jest Kruskal-Wallis**, niezależnie od tego, jak bardzo normalne i jednorodne były surowe dane techniczne. **Zweryfikowano na żywo**: syntetyczny plik w starym formacie z jawnie normalnymi, homogenicznymi danymi surowymi (Shapiro p=0.6–0.8 na surowych wartościach) nadal dał `test_used: Kruskal-Wallis` po agregacji do n_bio=1. Dopiero plik w NOWYM formacie z realnym `Rep_biologiczna` ≥3/grupę faktycznie uruchomił ANOVA (zweryfikowano, patrz Znalezisko 2).

(c) Rozbieżność: drzewo decyzyjne opisane w tekście jest poprawne jako abstrakcyjny algorytm, ale dla ogromnej, dotychczas dominującej klasy plików (stary format) jedna z dwóch gałęzi jest martwa — a tekst nigdzie o tym nie informuje.

(d) Powaga: **ISTOTNE**.

---

**Znalezisko 7 — "Grupy z N<2 są automatycznie pomijane" — nieaktualne po naprawie audytu 1.6**

(a) Cytat: *"Minimalna liczebność próby: Grupy posiadające mniej niż 2 wyniki (N < 2) są automatycznie pomijane w analizie statystycznej, ponieważ niemożliwe jest obliczenie dla nich odchylenia standardowego."*

(b) Co robi kod naprawdę: `logic.py::run_statistics`, komentarz i kod wprost mówią co innego: *"Filtrujemy tylko grupy calkowicie puste. n=1... jest dopuszczalne — test wciaz sie wykonuje, ale wynik jest orientacyjny"* — filtrowane są wyłącznie grupy **całkowicie puste** (n=0); grupy z n=1 są WŁĄCZANE do testu, z osobnym ostrzeżeniem "brak replikacji biologicznej" (widocznym m.in. jako baner na wykresie).

(c) Rozbieżność: tekst opisuje starsze zachowanie (twarde odrzucanie N<2), zastąpione świadomie w ramach napraw audytu (1.6) modelem "włącz i jawnie ostrzeż" zamiast "po cichu pomiń".

(d) Powaga: **ISTOTNE**.

---

**Znalezisko 8 — Opis wykresu Lollipop nie uwzględnia nowego domyślnego trybu "vs referencja" (poprawka UX-2)**

(a) Cytat: *"UWAGA: Wykres prezentuje wyłącznie pary różniące się istotnie statystycznie (p < α), aby zachować czytelność."*

(b) Co robi kod naprawdę: od niedawnej poprawki UX (`plotting.py::draw_effect_plot`, `show_all_pairs=False` domyślnie), wykres domyślnie pokazuje jeszcze węższy podzbiór — **tylko istotne pary z udziałem grupy referencyjnej**, nie wszystkie istotne pary. Pełny widok "wszystkie pary" jest teraz opcją w "Opcje Wykresu" (`switch_effect_all_pairs`), a nie zachowaniem domyślnym.

(c) Rozbieżność: opis w Podręczniku jest nieaktualny o jeden poziom szczegółowości — nie tłumaczy, że domyślny widok jest teraz ograniczony do porównań z referencją, ani że istnieje przełącznik "wszystkie pary".

(d) Powaga: **ISTOTNE**.

---

**Znalezisko 9 — Brak jakiejkolwiek instrukcji przygotowania nowego, wieloarkuszowego formatu wejścia**

(a) Cytat: jedyna wzmianka o formacie danych to *"Wykres Trendu - Nazewnictwo: nazwa grupy w Excelu MUSI zawierać liczbę i jednostkę... Poprawnie: 'Ekstrakt (50 mg/ml)'..."*.

(b) Co robi kod naprawdę: aplikacja obsługuje od dawna wieloarkuszowy format (`Dane_dyfuzja`, `MIC_wizualny`, `MIC_OD`, `MBC_posiew`, `Kontrole`) z kolumnami `Substancja`, `Stezenie`, `Jednostka`, `Typ`, `Rep_biologiczna`, `Rep_techniczna`. Sama konwencja nazewnicza opisana w tym wpisie pozostaje technicznie trafna (auto-budowana etykieta grupy dla badanych substancji w nowym formacie faktycznie przyjmuje postać `"Substancja (Stężenie Jednostka)"` — zweryfikowano w `utils.py`, funkcja `_label`), ale to jedyne zdanie w całym Podręczniku dotyczące struktury pliku wejściowego. Nie ma ani słowa o tym, jakie arkusze/kolumny trzeba przygotować, ani o arkuszu `Kontrole`, ani o różnicy powtórzenie biologiczne/techniczne, ani o wciąż-widocznym na każdym wykresie starego-formatu banerze "brak replikacji biologicznej".

(c) Rozbieżność: nie błąd, ale dotkliwa luka — użytkownik przygotowujący dane od zera (zwłaszcza pod MIC/MBC) nie znajdzie w aplikacji żadnej wskazówki, jak zbudować plik.

(d) Powaga: **ISTOTNE**.

---

**Znalezisko 10 — Opis Boxplotu zakłada niezdegenerowany rozkład, który nie istnieje dla starego formatu**

(a) Cytat: *"Pudełko: Obejmuje 50% środkowych wyników (od 25. do 75. percentyla). Wąsy: Zasięg danych (min-max), z wyłączeniem wartości odstających."*

(b) Co robi kod naprawdę: ten sam mechanizm co w Znalezisku 4/6 — dla starego formatu `df_bio` ma dokładnie 1 wiersz na grupę, więc boxplot rysowany jest na pojedynczym punkcie na grupę (brak IQR, brak wąsów do pokazania).

(c) Rozbieżność: opis zakłada wielopunktowy rozkład w obrębie grupy, którego dla starego formatu po prostu nie ma na wejściu do wykresu.

(d) Powaga: **ISTOTNE** (ta sama przyczyna źródłowa co Znalezisko 4/6, ale osobna, konkretna rycina/opis).

---

### DROBNE

---

**Znalezisko 11 — Nazwa arkusza Excel nie zgadza się dosłownie**

(a) Cytat: *"Pełne wyniki... znajdziesz w raporcie Excel (zakładka 'Post-hoc Details')."*

(b) Co robi kod naprawdę: `gui.py:856`, faktyczna nazwa arkusza to `"Post-hoc (Details)"` (z nawiasami), nie `"Post-hoc Details"`.

(c) Rozbieżność: kosmetyczna literówka w nazwie, łatwa do znalezienia mimo to, ale dosłowne wyszukanie frazy z Podręcznika w Excelu jej nie znajdzie.

(d) Powaga: **drobne**.

---

**Znalezisko 12 — Generator podpisów nie ma żadnego zabezpieczenia przed uruchomieniem bez analizy/z niejednoznaczną referencją**

(a) Zweryfikowane na żywo, dwa scenariusze:
- Przed wczytaniem jakiegokolwiek pliku: `Statistical significance was determined using Statistical test followed by Holm-Bonferroni correction...` / `...compared to the negative control (...).`
- Z ręcznie/domyślnie nierozwiązaną niejednoznaczną referencją: `...compared to the negative control (-- Wybierz ręcznie (niejednoznaczne) --).`

(b) Co robi kod naprawdę: `btn_captions` nie jest nigdy blokowany warunkowo (jedyne wystąpienie w kodzie to `gui.py:150-151`, bez późniejszego `.configure(state=...)`). `open_caption_window` nie sprawdza `self.df`, `self.export_stats_main`, ani czy `ref_group` różni się od `REF_PLACEHOLDER`/`"..."`, zanim wstawi je do tekstu.

(c) Rozbieżność: nie literalny błąd merytoryczny (przy normalnym użyciu — wczytaj plik, uruchom analizę, dopiero potem generuj podpisy — tekst jest sensowny), ale funkcja nie chroni przed łatwym do popełnienia, przedwczesnym kliknięciem, które wstawia bełkot wprost do pola tekstowego przeznaczonego do kopiowania do manuskryptu.

(d) Powaga: **drobne** (nie występuje przy prawidłowym użyciu, ale ryzyko realne i tanie do naprawienia).

---

**Znalezisko 13 — Sekcja "Korekty post-hoc" nie wspomina o czwartej, realnie istniejącej opcji "None"**

(a) Cytat: Sekcja 2 opisuje tylko Holm, Bonferroni, FDR.

(b) Co robi kod naprawdę: `combo_method` ma cztery wartości: `["holm", "fdr_bh", "bonferroni", "None"]` (`gui.py:114`) — "None" (brak korekty) to prawidłowa, wybieralna opcja.

(c) Rozbieżność: jedna z czterech opcji w interfejsie nie ma żadnego opisu w Podręczniku.

(d) Powaga: **drobne**.

---

## Co jest aktualne i poprawne

- **Poprawność techniczna (punkt 6 zlecenia):** obie funkcje wykonują się bez wyjątku we wszystkich przetestowanych scenariuszach (przed wczytaniem pliku, z niejednoznaczną referencją, z referencją = kontrola dodatnia, ze ścieżką ANOVA i Kruskal-Wallis) — nie odwołują się do usuniętych zmiennych/funkcji. Problemy leżą wyłącznie w TREŚCI generowanego tekstu, nie w stabilności kodu.
- **MIC z dyfuzji (punkt 1 zlecenia):** żaden z dwóch tekstów nigdzie nie wspomina o szacowaniu MIC z testu krążkowego (funkcja usunięta) — to konkretne ryzyko się nie zmaterializowało.
- Nazwy 6 rycin/zakładek w generatorze podpisów odpowiadają dokładnie nazwom zakładek w GUI (`Wykres Główny`, `Mapa Ciepła`, `Mapa P-value`, `Trend (Dawka)`, `Wielkość Efektu`, `Porównanie Szczepów`).
- Opis kolorystyki mapy p-value ("niebieski = istotne, czerwony/biały = nieistotne") zgadza się dokładnie z użytą paletą (`cmap="RdBu_r", center=ALPHA`).
- Zakres N=3–10 dla testu Dixona zgadza się dokładnie z kodem (`find_outliers_dixon`).
- Opis autokorekty spacji w nazwach (`_strip_string_cells`) i naturalnego sortowania (`smart_sort_key`) zgadza się z kodem.
- Konwencja nazewnicza wykresu trendu ("Substancja (Stężenie Jednostka)") pozostaje trafna również dla nowego formatu — etykieta grupy jest budowana automatycznie w tej samej konwencji.
- Ogólny opis logiki Shapiro-Wilk → Levene → ANOVA/Kruskal-Wallis jest poprawny jako ABSTRAKCYJNY algorytm (błąd dotyczy praktycznej osiągalności gałęzi ANOVA dla starego formatu — patrz Znalezisko 6 — nie samej logiki).
- Generator opisu metod w Podręczniku (Sekcja 4) poprawnie rozróżnia Tukey (ANOVA) od Dunn+korekta (Kruskal-Wallis) — w przeciwieństwie do generatora podpisów rycin (Znalezisko 2).

## Co wymaga aktualizacji

1. Generator podpisów: uzależnić etykietę "kontrola negatywna" od faktycznego typu wybranej grupy referencyjnej, nie zakładać tego bezwarunkowo (Znalezisko 1).
2. Generator podpisów: nazwać metodę post-hoc zgodnie z faktycznie użytym testem (Tukey dla ANOVA, wybrana korekta tylko dla Kruskal-Wallis) — skopiować logikę już poprawnie zaimplementowaną w Sekcji 4 Podręcznika (Znalezisko 2).
3. Generator podpisów: nie generować opisu Ryciny 3, gdy wykres wielkości efektu jest pusty (Znalezisko 3).
4. Oba teksty: doprecyzować, co faktycznie reprezentują słupki błędu/wąsy pudełka dla plików w starym formacie (n_bio=1 → brak zmienności do pokazania) (Znaleziska 4, 10).
5. Oba teksty: dodać sekcję o module MIC/MBC — podpisy dla 4 wykresów MIC/MBC, wyjaśnienie ilorazu/klasyfikacji/cenzury w Podręczniku (Znalezisko 5).
6. Podręcznik: zaktualizować opis drzewa decyzyjnego normalności o praktyczną nieosiągalność ANOVA dla starego formatu (Znalezisko 6).
7. Podręcznik: poprawić opis obsługi grup N<2 (włączane z ostrzeżeniem, nie pomijane) (Znalezisko 7).
8. Podręcznik: zaktualizować opis domyślnego zakresu wykresu Wielkości Efektu (vs referencja, nie wszystkie pary) (Znalezisko 8).
9. Podręcznik: dodać sekcję o przygotowaniu nowego, wieloarkuszowego formatu wejścia (Znalezisko 9).
10. Podręcznik: poprawić literalną nazwę arkusza Excel ("Post-hoc (Details)") (Znalezisko 11).
11. Generator podpisów: rozważyć blokadę przed uruchomieniem bez wczytanego pliku/zakończonej analizy/rozwiązanej referencji (Znalezisko 12).
12. Podręcznik: dodać opis opcji "None" (brak korekty) w sekcji o korektach post-hoc (Znalezisko 13).

*(Diagnoza wyłącznie — żadna z powyższych zmian nie została wprowadzona.)*
