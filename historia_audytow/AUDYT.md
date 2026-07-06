# AUDYT PROJEKTU: BioStat Master (v3.0)

**Folder audytowany:** `C:\Users\Ariva\Desktop\BioStat-Master_AKTUALNY` (klon z GitHub, ostatni commit 2026-04-16)
**Dane testowe:** `C:\Users\Ariva\Desktop\dane_disk.xlsx`
**Data audytu:** 2026-07-03
**Charakter audytu:** wyłącznie odczytujący — żaden plik projektu ani danych nie został zmieniony. Wygenerowane artefakty testowe zapisano w podfolderze `audyt_wyniki/`.

---

## 1. CEL I KONTEKST

**BioStat Master** to desktopowa aplikacja (Python + customtkinter) automatyzująca analizę statystyczną wyników **testu dyfuzji krążkowej (disk diffusion assay)** — standardowej metody laboratoryjnej oceny aktywności przeciwdrobnoustrojowej substancji (ekstraktów roślinnych, związków chemicznych, antybiotyków referencyjnych) na podstawie średnicy strefy zahamowania wzrostu bakterii wokół krążka nasączonego badaną substancją.

Program przyjmuje surowe dane w formacie Excel (bakteria / grupa-substancja / średnica strefy w mm), automatycznie dobiera właściwy test statystyczny (parametryczny lub nie), porównuje badane substancje względem grupy referencyjnej (zwykle kontroli negatywnej), oblicza wielkość efektu, wykrywa potencjalne wartości odstające, szacuje orientacyjne MIC (minimalne stężenie hamujące) tam, gdzie dane na to pozwalają, oraz generuje gotowe do publikacji wykresy i raport PDF.

Odbiorcą jest najpewniej mikrobiolog/biotechnolog (albo student/doktorant) prowadzący przesiewowe badania aktywności przeciwbakteryjnej nowych substancji i potrzebujący szybkiej, ustandaryzowanej analizy bez ręcznego liczenia testów w R/SPSS — na co wskazują polskie etykiety GUI, wbudowany "Generator Opisów do Publikacji" (gotowe angielskie captions do manuskryptu) oraz podręcznik metodyczny tłumaczący dobór testu (`dialogs.py`, klasa `HelpDialog`).

---

## 2. ARCHITEKTURA

### Struktura plików (14 plików źródłowych/zasobów, bez `.git`/`__pycache__`/`.venv`)

| Plik | Linie | Rola |
|---|---|---|
| `main.py` | 4 | Punkt wejścia |
| `gui.py` | 675 | Widok/Kontroler (customtkinter) |
| `logic.py` | 237 | Model — `StatsEngine` (statystyka, PCA, MIC) |
| `plotting.py` | 369 | Widok — `Plotter` (wszystkie wykresy matplotlib/seaborn) |
| `utils.py` | 163 | Walidacja danych, parsowanie, Dixon Q, Cohen's d |
| `dialogs.py` | 295 | Okna dialogowe (outliery, pomoc, o autorze) |
| `reports.py` | 144 | Generowanie raportu PDF (reportlab) |
| `config.py` | 90 | Scentralizowane stałe |
| **Razem** | **1977** | |
| `README.md`, `LICENSE`, `requirements.txt`, `ikona.ico/.png` | — | Dokumentacja/zasoby |

### Punkt wejścia i uruchomienie
```python
# main.py
from gui import App
if __name__ == "__main__":
    app = App()
    app.mainloop()
```
Uruchomienie: `pip install -r requirements.txt` → `python main.py`. Weryfikacja: **potwierdzone empirycznie w tym audycie** — środowisko `.venv` (Python 3.11.9) zbudowane od zera, instalacja zależności bezbłędna, aplikacja startuje i wchodzi w pętlę zdarzeń bez wyjątków (zob. wcześniejszy smoke-test tej sesji).

### Przepływ danych (Excel → PDF)
```
dane_disk.xlsx
   │  pd.read_excel + strip nagłówków/wartości                       [gui.load_file]
   ▼
walidacja struktury (kolumny Bakterie/Grupa/Srednica_mm)              [utils.validate_excel_structure]
   ▼
walidacja wartości komórek (NaN, <0, >100, puste etykiety)            [utils.validate_excel_data]
   ▼
DataFrame oczyszczony → wybór szczepu + zaznaczonych grup w GUI
   ▼
detekcja outlierów (Dixon Q, n=3..10) → opcjonalne usunięcie          [utils.find_outliers_dixon, dialogs.OutlierDialog]
   ▼
StatsEngine.run_statistics: Shapiro → Levene → ANOVA/Kruskal-Wallis   [logic.py]
   ▼
StatsEngine.process_detailed_results: post-hoc + Cohen's d            [logic.py]
   ▼
Plotter: 7 wykresów (bar/heat/pvalue/trend+MIC/effect/cross/PCA)      [plotting.py]
   ▼
reports.generate_pdf: tabela opisowa + wykresy + werdykt tekstowy     [reports.py]
   ▼
raport.pdf / eksport Excel / PNG HQ
```

### Rozdzielenie logiki od GUI
Rozdział jest **rzeczywisty i konsekwentny** — to mocna strona tego projektu. `logic.py` (`StatsEngine`) nie importuje `tkinter`/`customtkinter` i nie ma żadnej zależności od GUI; przyjmuje i zwraca czyste `DataFrame`/`dict`. `plotting.py` (`Plotter`) buduje obiekty `matplotlib.Figure` i nigdy nie wywołuje `plt.show()` ani nie odwołuje się do widgetów Tk — dzięki temu dało się w tym audycie uruchomić cały pipeline analityczny **bez GUI** (zob. sekcja 6), co jest dobrym testem architektury. Jedyne miejsce łączące oba światy to `gui.py`, które orkiestruje wywołania.

### Wykorzystanie 11 zależności z `requirements.txt`

| Zależność | Użyta? | Gdzie |
|---|---|---|
| customtkinter | ✅ | `gui.py`, `dialogs.py` |
| matplotlib | ✅ | `plotting.py`, `gui.py` (canvas) |
| numpy | ✅ | `logic.py`, `utils.py`, `plotting.py` |
| pandas | ✅ | wszędzie |
| scipy | ✅ | `logic.py` (shapiro/levene/f_oneway/kruskal/linregress), `plotting.py` (spearmanr) |
| statsmodels | ✅ | `logic.py` (`pairwise_tukeyhsd`) |
| scikit-posthocs | ✅ | `logic.py` (`posthoc_dunn`) |
| scikit-learn | ✅ | `logic.py` (`PCA`, `StandardScaler`) |
| seaborn | ✅ | `plotting.py` |
| reportlab | ✅ | `reports.py` |
| openpyxl | ✅ pośrednio | silnik `pd.read_excel`/`pd.ExcelWriter` w `gui.py` — nigdy nie importowany wprost, ale wymagany transytywnie |

**Żadna z 11 zależności nie jest martwa.** Brak nadmiarowych bibliotek.

---

## 3. FUNKCJONALNOŚĆ (co użytkownik może zrobić)

**Wejście:** plik Excel (`.xlsx`/`.xls`) z kolumnami zawierającymi w nazwie "Bakteri..." (dowolna), `Grupa`, `Srednica_mm`.

**Konfiguracja analizy:**
- wybór szczepu bakterii do analizy,
- zaznaczanie/odznaczanie dowolnego podzbioru grup (checkboxy, "zaznacz/odznacz wszystko"),
- wybór metody korekty post-hoc: Holm (domyślna), Benjamini-Hochberg (FDR), Bonferroni, brak,
- wybór grupy referencyjnej (kontroli) z listy rozwijanej,
- orientacja wykresu głównego (pionowa/pozioma).

**Analizy wykonywane po kliknięciu "URUCHOM ANALIZĘ":**
1. Test normalności Shapiro-Wilka dla każdej grupy,
2. Test jednorodności wariancji Levene'a,
3. Automatyczny wybór: ANOVA (parametryczny) albo Kruskal-Wallis (nieparametryczny),
4. Post-hoc: Tukey HSD (po ANOVA) albo Dunn's Test z wybraną korekcją (po Kruskal-Wallis),
5. Cohen's d dla każdej istotnej pary grup + interpretacja słowna (znikomy/mały/średni/DUŻY),
6. Detekcja outlierów testem Q Dixona (automatycznie dla n=3..10) z możliwością wykluczenia,
7. Szacowanie MIC metodą regresji log-liniowej (gdy nazwa grupy zawiera stężenie),
8. Korelacja Spearmana dla trendu dawka-odpowiedź,
9. PCA (grupowanie szczepów bakterii wg profilu wrażliwości — wymaga ≥3 szczepów i ≥2 substancji).

**Wyjścia:**
- 7 zakładek wykresów: Barplot/Boxplot/Violinplot główny (z gwiazdkami istotności), Mapa Ciepła aktywności, Mapa P-value, Trend Dawka-Odpowiedź (z adnotacją MIC), Wielkość Efektu (lollipop), Porównanie Międzygatunkowe, PCA,
- log tekstowy z pełnym raportem statystycznym,
- eksport pojedynczego wykresu do PNG (300 DPI) lub PDF,
- eksport pełnych danych do Excela (4 arkusze: dane surowe, normalność, test główny, post-hoc szczegółowy),
- pełny raport PDF (metryczka, tabela opisowa, wszystkie wykresy, werdykt tekstowy istotnych różnic),
- generator gotowych angielskich opisów rycin ("Figure 1. Antibacterial activity...") do wklejenia w manuskrypt,
- wbudowany podręcznik metodyczny i generator automatycznego opisu metod ("Materials and Methods") dopasowanego do faktycznie użytego testu.

---

## 4. METODOLOGIA STATYSTYCZNA I NAUKOWA

To najważniejsza część audytu — poniżej każdy punkt z dokładnym cytatem kodu.

### 4.1 Drzewo decyzyjne testów i sprawdzanie założeń

`logic.py:37-57`:
```python
# 1. Normalność
all_normal = True
normality_results = []
for g in valid_groups:
    vals = df_run[df_run[COL_GROUP] == g][COL_MEASUREMENT]
    p_shapiro = 0
    is_norm = False
    if len(vals) >= 3 and vals.std() > 0:
        s, p_shapiro = stats.shapiro(vals)
        if p_shapiro >= ALPHA: is_norm = True
    if not is_norm: all_normal = False
    normality_results.append(...)

# 2. Levene
stat, p_levene = stats.levene(*dane_list)
use_parametric = all_normal and p_levene > ALPHA
```
Test główny: ANOVA (`stats.f_oneway`) jeśli **wszystkie** grupy normalne i wariancje jednorodne; w przeciwnym razie Kruskal-Wallis (`stats.kruskal`) — logika poprawna i standardowa. **Ale:**

- **Reguła "wszystko albo nic".** Jedna niespełniająca założeń grupa (np. jedna z 30) wymusza test nieparametryczny dla **całego** zestawu porównań, nawet jeśli pozostałe 29 grup jest w pełni normalnych. To podejście konserwatywne (bezpieczne), ale w praktyce z małymi próbami (n=3, patrz niżej) niemal zawsze ląduje się w Kruskal-Wallis.
- **Brak realnej mocy testu Shapiro-Wilka przy n=3.** Kod wymaga `len(vals) >= 3` do w ogóle podjęcia próby (minimalna liczba dla `scipy.stats.shapiro`), ale test normalności przy n=3 ma znikomą moc statystyczną — praktycznie nie da się nim wykryć odstępstwa od normalności ani go potwierdzić. W realnych danych (`dane_disk.xlsx`, sekcja 6) grupy mają dokładnie n=3 — czyli wynik "normalny"/"nienormalny" dla pojedynczej grupy jest w dużej mierze przypadkowy.
- **Grupy o zerowej wariancji (`vals.std() == 0`) automatycznie oznaczane jako nienormalne** (`p_shapiro=0`, `is_norm=False`), bez żadnego komunikatu dla użytkownika wyjaśniającego, że test w ogóle się nie wykonał (a nie: "wykonał się i wykazał brak normalności"). W realnych danych aż **15/30 grup (E. coli)** i **9/30 grup (S. aureus)** miało dokładnie taki przypadek (identyczne 3 pomiary) — więcej niż jedna trzecia wyniku "nienormalny" pochodzi z tego mechanizmu zastępczego, a nie z faktycznego testu.

### 4.2 Interpretacja średnic wg CLSI/EUCAST — **BRAK**

Sprawdzono cały kod (`grep -i "CLSI|EUCAST|breakpoint|susceptib|resistan|intermediate|M100"`) — **jedyne trafienie** to fragment angielskiego opisu rycin w `gui.py:346` ("differential *susceptibility* of tested pathogens"), czysto opisowy zwrot językowy, nie logika programu.

**Program nie zawiera żadnych tabel progów CLSI (M100) ani EUCAST i nie przypisuje kategorii S/I/R (Susceptible/Intermediate/Resistant).** Cała interpretacja opiera się wyłącznie na:
- porównaniu średniej średnicy strefy między grupami (test statystyczny + Cohen's d),
- linii referencyjnej na wykresie oznaczającej **średnicę krążka** (`DISC_DIAMETER_MM = 6.0`, `config.py:15`), a nie próg kliniczny.

To fundamentalna różnica względem standardu klinicznej mikrobiologii: prawdziwe kategorie S/I/R są **specyficzne dla pary bakteria-antybiotyk** (inny próg dla *E. coli*/ampicylina niż dla *S. aureus*/ampicylina) i pochodzą z tabel CLSI/EUCAST kalibrowanych względem MIC z bulionowej metody referencyjnej. Program tego nie robi — i szczerze mówiąc **nie powinien tego udawać** dla ekstraktów roślinnych/związków niebędących zarejestrowanymi antybiotykami (dla nich nie istnieją oficjalne progi CLSI/EUCAST). Jest to więc uzasadnione ograniczenie zakresu, ale musi być jasno zakomunikowane użytkownikowi — obecnie **nie jest** (README/GUI nie wspominają wprost, że S/I/R nie jest tu wyliczane).

### 4.3 Szacowanie MIC ze średnicy strefy — dokładny kod

`logic.py:175-237` (`estimate_mic`):
```python
def estimate_mic(self, df, selected_substances, target_diameter=DISC_DIAMETER_MM):
    """
    Estimates MIC for each substance using Log-Linear Regression.
    Model: Diameter = a + b * ln(Concentration)
    MIC = exp((Target - a) / b)
    """
    ...
    if len(set(x_concs)) < 3:
        continue          # <-- pomija substancję BEZ ŻADNEGO komunikatu

    log_x = np.log(x_concs)
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, y_diams)

    if slope > 0:
        ln_mic = (target_diameter - intercept) / slope
        mic = np.exp(ln_mic)
    else:
        mic = None
    results[sub] = {"MIC": mic, "Unit": valid_unit, "R2": r_value**2, "Slope": slope, "Intercept": intercept}
```

Kluczowe fakty:
- **To jeden uniwersalny wzór, nie kalibracja per antybiotyk/gatunek.** Każda substancja dostaje własną regresję log-liniową na podstawie *tylko swoich* punktów (stężenie, średnica), ale sam **model** (log-liniowy, cel = średnica krążka 6 mm) jest identyczny dla wszystkich — to **nie jest** metoda kalibracyjna CLSI (regresja MIC-vs-średnica budowana na wielu szczepach referencyjnych o znanym MIC z mikrorozcieńczeń bulionowych), tylko doraźna ekstrapolacja z własnych danych dawka-odpowiedź.
- **Cel ekstrapolacji (`target_diameter = 6.0 mm = średnica krążka`) jest założeniem, nie zweryfikowanym punktem biologicznym.** Zakłada, że "MIC" to stężenie, przy którym strefa zahamowania kurczy się dokładnie do rozmiaru samego krążka (brak widocznej strefy poza krążkiem). To sensowna heurystyka, ale nie jest to definicja MIC uznawana w mikrobiologii klinicznej (MIC pochodzi z mikrorozcieńczeń bulionowych, nie z ekstrapolacji strefy).
- **Wymaga tylko ≥3 unikalnych stężeń** — bez żadnego progu jakości dopasowania (R²). W realnych danych zaobserwowano dopasowanie R²=0.60 (S. aureus, substancja "53-084") wciąż zwracające konkretną liczbową wartość MIC bez ostrzeżenia.
- **Brak sprawdzenia, czy ekstrapolacja wykracza poza zakres przetestowanych stężeń.** W realnym przebiegu (sekcja 6) dla "53-084" (testowane 0.5 / 5 / 50 mg/ml) program zwrócił MIC = **0.0011 mg/ml** — wartość ok. 450× niższą niż najniższe przetestowane stężenie. To matematycznie poprawny wynik regresji, ale biologicznie ekstremalna ekstrapolacja poza zbadany zakres, bez ostrzeżenia w interfejsie.
- Substancje z <3 stężeniami lub błędem regresji są **pomijane bez jakiegokolwiek komunikatu** (`continue` / `except: pass`, `logic.py:208` i `234-235`) — użytkownik nie dowie się z GUI *dlaczego* dana substancja nie ma wyliczonego MIC.

### 4.4 Replikaty — ryzyko pseudoreplikacji

Schemat danych Excela to płasko: `Bakterie | Grupa | Srednica_mm` — **nie ma żadnej kolumny identyfikującej powtórzenie jako biologiczne vs. techniczne** (np. numer płytki/eksperymentu). Kod (`logic.py:24-32`) traktuje każdy wiersz danej grupy jako niezależną obserwację wchodzącą wprost do `n` w ANOVA/Kruskal-Wallis/Shapiro/Levene:
```python
for g in df_run[COL_GROUP].unique():
     data = df_run[df_run[COL_GROUP] == g][COL_MEASUREMENT].values
     if len(data) >= 2:
        valid_groups.append(g)
        dane_list.append(data)
```
Jeśli 3 pomiary w grupie to 3 techniczne powtórzenia z tej samej płytki/tego samego eksperymentu (co jest bardzo prawdopodobne w typowym teście dyfuzji krążkowej — patrz sekcja 6, gdzie wiele trójek ma identyczne wartości), to traktowanie ich jako 3 niezależne obserwacje **zawyża efektywne n i sztucznie zawyża moc statystyczną / zaniża p-value** — klasyczna pseudoreplikacja. Program nie ma mechanizmu rozróżnienia i nie ostrzega o tym ryzyku nigdzie (ani w kodzie, ani w podręczniku pomocy).

### 4.5 Korekta na wielokrotne porównania

Zaimplementowana i **skonfigurowana poprawnie w większości ścieżek**, ale z istotną niespójnością:

- Ścieżka Kruskal-Wallis → Dunn: `logic.py:79` `sp.posthoc_dunn(df_run, COL_MEASUREMENT, COL_GROUP, p_adjust=method)` — `method` pochodzi wprost z wyboru użytkownika (holm/fdr_bh/bonferroni/None). **Działa zgodnie z oczekiwaniem.**
- Ścieżka ANOVA → Tukey: `logic.py:70` `pairwise_tukeyhsd(df_run[COL_MEASUREMENT], df_run[COL_GROUP], ALPHA)` — **funkcja `pairwise_tukeyhsd` nie przyjmuje parametru `method` w ogóle.** Tukey HSD ma wbudowaną własną korektę (rozkład studentyzowanego rozstępu), niezależną od wyboru w GUI. Oznacza to, że **wybór metody post-hoc w interfejsie ("Holm"/"Bonferroni"/"FDR") jest po cichu ignorowany, gdy dane trafiają na ścieżkę parametryczną (ANOVA)** — użytkownik może wybrać np. "Bonferroni", a mimo to dostanie korektę Tukeya. Nie jest to błąd w sensie matematycznym (Tukey HSD to poprawna, ugruntowana korekta), ale jest to rozbieżność między tym, co GUI sugeruje, że się stanie, a tym, co faktycznie się dzieje — warto to udokumentować w interfejsie.

### 4.6 Powtarzalność (seed)

- **PCA** (`logic.py:166`: `PCA(n_components=2)`) — brak jawnego `random_state`. Dla małych macierzy (jak tutaj) `sklearn` domyślnie używa pełnego SVD (deterministycznego), więc w praktyce wynik jest powtarzalny przy tych rozmiarach danych, ale nie jest to **gwarantowane** przez kod — dla większych zbiorów `sklearn` może przełączyć się na randomizowany solver i wtedy wynik przestałby być identyczny między uruchomieniami bez ustawionego seeda.
- Brak bootstrapu w kodzie.
- Brak innych losowych elementów.

**Zalecenie:** dodać `PCA(n_components=2, random_state=42)` — to jednowierszowa, tania zmiana eliminująca całe ryzyko.

---

## 5. JAKOŚĆ I POPRAWNOŚĆ KODU

### Mocne strony
- Walidacja wejścia jest **rzeczywiście solidna**: strukturalna (`utils.validate_excel_structure`) i wartościowa (`utils.validate_excel_data`) z precyzyjnymi komunikatami zawierającymi numer wiersza Excela. Pomiar `0` (brak strefy) jest **prawidłowo dopuszczony** (`utils.py:141-142`), wartości ujemne/nieskończone/>100mm/puste — odrzucane z podaniem przyczyny.
- Obsługa błędów wczytania pliku (`gui.py:176-188`) rozróżnia `FileNotFoundError`/`PermissionError`/`ValueError` z osobnymi komunikatami — dobra praktyka.
- Stałe naukowe scentralizowane w `config.py` z opisowymi docstringami — bardzo dobra higiena wobec wersji sprzed refaktoryzacji.

### Realny błąd wykryty dynamicznie: awaria "Mapy P-value" na aktualnym środowisku

Podczas testu na prawdziwych danych (sekcja 6) wykres "Mapa P-value" **rzucił wyjątek za każdym razem**, gdy użyto ścieżki Kruskal-Wallis/Dunn (czyli w praktyce zawsze na tych danych). Zlokalizowano i potwierdzono przyczynę:

`plotting.py:159-161`:
```python
else: # Dunn
    p_matrix = export_stats_posthoc.copy()
    np.fill_diagonal(p_matrix.values, 1.0)
```
Zweryfikowano bezpośrednio: `p_matrix.values.flags.writeable` → `False`, mimo jawnego `.copy()`. Przyczyna: `pandas 3.0.3` (zainstalowany, bo `requirements.txt` **nie pinuje wersji**) wprowadza domyślne Copy-on-Write, przez co bufor numpy pod `.values` bywa oznaczony jako tylko-do-odczytu nawet po `.copy()`, dopóki modyfikacja nie przejdzie przez interfejs pandas. `np.fill_diagonal` pisze bezpośrednio do bufora numpy → `ValueError: assignment destination is read-only / underlying array is read-only`.

**Skutek w praktyce:** w żywym GUI wyjątek jest łapany (`gui.py: display_plot` → `_show_plot_error`), więc aplikacja się nie wywala, ale zakładka "Mapa P-value" pokazuje błąd zamiast wykresu, a w raporcie PDF ta strona **po cichu w ogóle się nie pojawia** (bez żadnej wzmianki, że coś zostało pominięte) — potwierdzone empirycznie w wygenerowanych plikach `raport_E_coli.pdf` / `raport_S_aureus.pdf` w `audyt_wyniki/`.

**Naprawa** (nie wykonana w ramach tego audytu, bo miał być tylko odczytujący): zbudować macierz jako świeżą tablicę numpy przed `fill_diagonal`, np.:
```python
arr = np.asarray(export_stats_posthoc, dtype=float).copy()
np.fill_diagonal(arr, 1.0)
p_matrix = pd.DataFrame(arr, index=export_stats_posthoc.index, columns=export_stats_posthoc.columns)
```
oraz przypięcie wersji w `requirements.txt` (np. `pandas>=2,<3`), żeby uniknąć podobnych niespodzianek przy przyszłych instalacjach.

### Realny błąd logiczny wykryty dynamicznie: auto-wybór grupy referencyjnej

`gui.py:253`:
```python
woda = next((g for g in grupy_bact if "woda" in g.lower() or "kontrol" in g.lower()), None)
if woda: self.combo_ref.set(woda)
```
Zweryfikowano na realnych danych: gdy istnieją obie grupy **"Kontrola (+) Ampycylina"** (kontrola pozytywna — antybiotyk) i **"Kontrola (-) Woda"** (kontrola negatywna — woda), program automatycznie wybiera jako referencję **"Kontrola (+) Ampycylina"**, ponieważ:
1. dopasowanie jest jednym wspólnym warunkiem OR ("woda" **lub** "kontrol") bez rozróżnienia `(+)`/`(-)`,
2. `next()` bierze pierwszy pasujący element z listy posortowanej przez `smart_sort_key`, a `"Kontrola (+)..."` sortuje się alfabetycznie przed `"Kontrola (-)..."` (znak `+` ma niższy kod ASCII niż `-`).

**To realny problem naukowy, nie tylko kosmetyczny:** gwiazdki istotności na wykresie głównym i logika "grupa referencyjna" w całej analizie odnoszą się wtedy do "różni się istotnie od antybiotyku", a nie standardowego "różni się istotnie od nieleczonej kontroli" — co jest zupełnie inną hipotezą naukową. Użytkownik musi **ręcznie poprawić** wybór grupy referencyjnej w GUI (rozwijana lista pozwala to zmienić), ale domyślne zachowanie jest mylące i łatwe do przeoczenia.

### Test Q Dixona przy n=3 — nadmierna czułość, potwierdzona na realnych danych

Krytyczne wartości Dixona (`config.py:84-90`, alpha=0.10) w połączeniu z n=3 są ekstremalnie czułe: przy trzech pomiarach, jeśli dwa są identyczne a trzeci różni się choćby o 1 mm, obliczone Q często przekracza próg krytyczny (0.941 dla n=3). Potwierdzono empirycznie: na realnych danych aż **22 z 30 grup (73%) dla E. coli** i **19 z 30 (63%) dla S. aureus** zostało oznaczonych jako zawierające "outlier" — w większości przypadków to pojedynczy pomiar różniący się o 1 mm od dwóch identycznych. To wskazuje, że test Dixona przy n=3 **praktycznie nie nadaje się do sensownego różnicowania błędu pomiaru od naturalnej zmienności biologicznej** w tych warunkach.

Dodatkowo: w oknie dialogowym wykluczania outlierów (`dialogs.py:33`) każdy wykryty punkt jest **domyślnie zaznaczony do usunięcia** (`var = ctk.IntVar(value=1)`), czyli działanie jest "opt-out", nie "opt-in". W połączeniu z powyższą nadczułością, nieuważne kliknięcie "Potwierdź i Analizuj" może usunąć większość punktów, które wcale nie są błędami pomiaru, tylko zwykłą zmiennością — sztucznie zawyżając pozorną precyzję/istotność wyników.

### Inne obserwacje
- **Case-sensitive wykrywanie kolumny bakterii** (`config.py:34-36`, `COL_BACT_SUBSTRING = 'Bakteri'`) — jeśli użytkownik nazwie kolumnę inaczej niż oczekiwany wzorzec, program **zgłasza to jawnie** jako błąd walidacji (nie jest to cichy błąd), ale wymaga dokładnej pisowni.
- Brak normalizacji wielkości liter w nazwach grup (np. "Kontrola" i "kontrola" byłyby dwiema różnymi grupami) — program już usuwa białe znaki (`gui.py:198-199`), ale nie różnice wielkości liter.
- Magiczne liczby rozrzucone poza `config.py`: próg "mała próba" (`total < 10`, `gui.py:457`), próg "mało obserwacji w grupie" (`n < 3`, `gui.py:446`), minimalna liczba szczepów do PCA (`< 3`, `logic.py:151`) — działają poprawnie, ale nie są scentralizowane jak reszta stałych w `config.py`, co utrudnia ich spójną zmianę w jednym miejscu.
- Generowanie PDF: obsługa braku `arial.ttf`/`arialbd.ttf` z łagodnym fallbackiem do Helvetiki (`reports.py:26-33`) — dobra praktyka defensywna (choć wtedy polskie znaki diakrytyczne mogą się źle wyświetlić w PDF — nie testowano tego konkretnie w tym audycie, bo generowane tu treści PDF były głównie angielskie/liczby).

---

## 6. DYNAMICZNY TEST NA REALNYCH DANYCH

### Struktura `dane_disk.xlsx`
- **180 wierszy**, 3 kolumny: `Bakterie` (str), `Grupa` (str), `Srednica_mm` (float64).
- **2 szczepy bakterii:** *E. coli*, *S. aureus*.
- **30 unikalnych grup na każdy szczep** (60 łącznie, każda z dokładnie **3 powtórzeniami** = 180 wierszy). Grupy obejmują: kontrolę negatywną (woda), kontrolę pozytywną (ampicylina), oraz kilkanaście kodowanych substancji (np. `55-156`, `C1`–`C4`) testowanych przy 1-3 różnych stężeniach.
- **Brak braków danych (NaN), brak wartości ujemnych/>100mm.** 6 pomiarów o wartości `0` (brak strefy zahamowania — poprawnie zaakceptowane przez walidację).

### Przebieg pipeline'u (nieinteraktywnie, bez GUI, przez bezpośrednie wywołania `StatsEngine`/`Plotter`/`reports.generate_pdf`)

Dla obu szczepów uruchomiono pełną analizę na wszystkich 30 grupach z korektą Holm, referencją auto-wybraną przez program:

| | E. coli | S. aureus |
|---|---|---|
| Liczba grup / obserwacji | 30 / 90 | 30 / 90 |
| Grupa referencyjna (auto) | **Kontrola (+) Ampycylina** ⚠ | **Kontrola (+) Ampycylina** ⚠ |
| Grup normalnych (Shapiro) | 1/30 | 5/30 |
| Test główny | Kruskal-Wallis | Kruskal-Wallis |
| Statystyka / p-value | H=87.69, p=8.2e-08 | H=84.70, p=2.3e-07 |
| Porównań post-hoc (Holm) | 435 | 435 |
| Istotnych statystycznie | 4 | 3 |
| Wykryte "outliery" Dixona | 22/30 grup (73%) | 19/30 grup (63%) |
| Substancji kwalifikujących się do MIC (≥3 stężenia) | 2 z 15 | 2 z 19 |
| PCA | odrzucone: "Za mało danych do PCA (wymagane min. 3 szczepy)" — **poprawne zachowanie**, w danych są tylko 2 szczepy | jw. |
| Mapa P-value | ❌ awaria (`ValueError: underlying array is read-only`) | ❌ awaria (jw.) |
| Raport PDF | ✅ wygenerowany, 829 KB | ✅ wygenerowany, 842 KB |

**Przykładowe wyniki MIC** (jedyne substancje z ≥3 stężeniami):
- E. coli / "55-156": MIC=0.34 mg/ml, R²=0.87
- E. coli / "C3": MIC=0.06%, R²=0.75
- S. aureus / "C3": MIC=6.17%, R²=0.76
- S. aureus / "53-084": MIC=**0.0011 mg/ml**, R²=0.60 — ekstrapolacja daleko poza przetestowany zakres (0.5–50 mg/ml), przy słabym dopasowaniu; wynik liczbowo "wygląda" precyzyjnie, ale nie ma pokrycia w faktycznie zbadanych stężeniach.

### Ocena sensowności wyników
- **Brak NaN-ów tam, gdzie nie powinno ich być** — wszystkie 60 grup dało kompletne statystyki opisowe.
- **Brak przypadków "Cohen's d" niepoliczalnego** (0 NaN/inf na 435+435 porównań).
- Wysoce istotne p-value testu głównego (p<1e-6) jest **oczekiwane i wiarygodne** przy 30 bardzo różnych grupach (kontrola negatywna 0mm vs. kontrola pozytywna ~27mm w tym samym zbiorze gwarantuje ogromną różnicę) — nie jest to sygnał błędu, tylko odzwierciedlenie realnie zróżnicowanych danych.
- **Jedyne "absurdalne" wartości to opisany wyżej MIC=0.0011 mg/ml** — matematycznie poprawny produkt regresji, ale biologicznie wątpliwy z powodu głębokiej ekstrapolacji.
- Oba wykryte błędy (mapa p-value, auto-wybór referencji) są **odtwarzalne w 100% przypadków** na tym zbiorze danych i tym środowisku — nie są to incydentalne/rzadkie usterki.

### Wygenerowane artefakty (`audyt_wyniki/`)
```
audyt_wyniki/
├── pipeline_log.txt                 (pełny log tekstowy przebiegu)
├── E_coli/
│   ├── 1_bar.png, 2_heat.png, 4_trend.png, 5_effect.png, 6_cross.png
│   ├── opisowe.csv, posthoc_details.csv
│   └── raport_E_coli.pdf            (829 KB, poprawny nagłówek %PDF-1.4)
└── S_aureus/
    └── (analogicznie, raport_S_aureus.pdf, 842 KB)
```
(Brak `3_pvalue.png` i `7_pca.png` — odpowiednio z powodu opisanej awarii i poprawnie zadziałanego zabezpieczenia "za mało szczepów".)

---

## PODSUMOWANIE KOŃCOWE

**(a) Czy podejście jest merytorycznie poprawne naukowo?**
W dużej mierze **tak, dla zastosowania jako narzędzie przesiewowe/eksploracyjne**: właściwy dobór testu parametryczny/nieparametryczny na podstawie realnie sprawdzanych założeń, poprawnie zaimplementowana korekta wielokrotnych porównań (przynajmniej na ścieżce Kruskal-Wallis), rzetelna walidacja danych wejściowych i przejrzysta, powtarzalna metodologia opisana w kodzie. To nie jest "fałszywa nauka" — to solidnie napisany, transparentny silnik statystyczny bez ukrytych manipulacji.

**(b) Najsłabsze punkty metodologiczne:**
1. **Pseudoreplikacja** — brak rozróżnienia powtórzeń technicznych od biologicznych w schemacie danych i w analizie (sekcja 4.4) — to prawdopodobnie najpoważniejsza luka naukowa.
2. **Test Shapiro-Wilka przy n=3** ma znikomą moc, a mimo to steruje wyborem całej ścieżki analizy dla wszystkich grup naraz.
3. **Szacowanie MIC** to uproszczona ekstrapolacja bez progu jakości dopasowania (R²) i bez ostrzeżenia przy ekstrapolacji poza zbadany zakres — potwierdzone konkretnym absurdalnym wynikiem na realnych danych.
4. **Auto-wybór grupy referencyjnej** może po cichu wskazać kontrolę pozytywną zamiast negatywnej (potwierdzone empirycznie).
5. **Test Dixona przy n=3** jest nadmiernie czuły (73%/63% grup oznaczonych jako zawierające outlier), a domyślne UI ułatwia niezamierzone, nadmierne usuwanie danych.
6. **Brak kategoryzacji CLSI/EUCAST** — akceptowalne jako świadome ograniczenie zakresu, ale niekomunikowane wprost użytkownikowi.

**(c) Czego brakuje, żeby to było narzędzie publikowalne/laboratoryjne:**
- Kolumna/mechanizm rozróżniający powtórzenia biologiczne od technicznych (np. model mieszany albo agregacja do średniej biologicznej przed testem głównym).
- Próg R² (i/lub ostrzeżenie przy ekstrapolacji poza zakres) w `estimate_mic`.
- Poprawka logiki auto-wyboru referencji (rozróżnienie `(+)`/`(-)` albo jawne pytanie do użytkownika zamiast zgadywania).
- Naprawa awarii Mapy P-value (przypięcie wersji pandas w `requirements.txt` + poprawka `plotting.py`).
- Jawna informacja w README/GUI, że wyniki nie są kategoriami klinicznymi S/I/R i nie zastępują interpretacji CLSI/EUCAST dla substancji, dla których takie progi istnieją.
- Podstawowy zestaw testów jednostkowych (obecnie: zero) — szczególnie dla `logic.py` i `utils.py`, gdzie regresje przy przyszłych zmianach (jak pokazał przypadek pandas 3.0) mogą przechodzić niezauważone.
- Rozważenie nieco większego progu minimalnej liczebności próby (obecnie N=2 wystarcza do przeprowadzenia testu, tylko z ostrzeżeniem, które można zignorować) dla wyników mających trafić do publikacji.

**Werdykt ogólny:** to dojrzałe, dobrze zaprojektowane narzędzie robocze do szybkiej, ustandaryzowanej analizy przesiewowej — nadaje się do wspierania własnych badań i wstępnej orientacji w wynikach, ale wymaga wymienionych poprawek (szczególnie punktów 1-4 z sekcji b) zanim wyniki z niego trafią bezpośrednio do publikacji recenzowanej bez dodatkowej, ręcznej weryfikacji statystycznej.
