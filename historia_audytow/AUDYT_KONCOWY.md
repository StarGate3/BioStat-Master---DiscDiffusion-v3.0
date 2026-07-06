# AUDYT KOŃCOWY — BioStat Master (gałąź `main`)

**Data audytu:** 2026-07-04
**Zakres:** moduł dyfuzji krążkowej + moduł MIC/MBC (broth microdilution), stan po scaleniu gałęzi `mic-mbc` do `main` (commit `13c61ae`).
**Charakter audytu:** WYŁĄCZNIE odczytujący. Nic nie zmieniono, nic nie naprawiono, nic nie zacommitowano. Wszystkie liczby w tym raporcie pochodzą z niezależnych skryptów weryfikacyjnych uruchomionych podczas audytu (nie z wcześniejszych deklaracji w komentarzach kodu ani z wcześniejszych sesji).
**Nastawienie:** recenzent zewnętrzny, zadanie = znaleźć słabości, nie potwierdzić poprawność.

---

## Streszczenie

Rdzeń obliczeniowy (wzory na stężenia studzienek, próg OD, redukcja CFU, mediana log2, test rangowy z cenzurą, zaokrąglenia zmiennoprzecinkowe) jest **poprawny i zweryfikowany niezależnym rachunkiem** — nie znalazłem błędu w żadnym z tych miejsc. Wsteczna zgodność ze starym formatem jest **potwierdzona bit-w-bit** (identyczny hash średnich, identyczna statystyka Kruskal-Wallisa i p-value, przed i po całym module MIC/MBC).

Natomiast poza tym rdzeniem znalazłem **11 konkretnych zastrzeżeń**, w tym **3 krytyczne**: moduł MIC/MBC nie jest wpięty do GUI (funkcjonalnie niedostępny dla użytkownika), jeden nienumeryczny wpis w kolumnie `Stez_S1`/`Wsp_rozc` wywala przetwarzanie **całego arkusza** (nie tylko jednego wiersza), a brak lub niespójny typ `Rep_biologiczna` potrafi po cichu scalić rzeczywiście różne powtórzenia biologiczne w jedno — zademonstrowane liczbowo (3 wartości MIC różniące się 64-krotnie zlewają się w wynik "32–32 mg/ml", zero zakresu).

---

## 1. Znaleziska (najpoważniejsze najpierw)

### 1.1 — Moduł MIC/MBC nie jest wpięty do GUI
**Czego dotyczy:** `gui.py` (całość); porównanie z `mic_logic.py`, `mic_plotting.py`.
**Co jest nie tak:** `gui.py` importuje i wywołuje wyłącznie `utils.route_workbook` (Faza 1 — wykrywanie/dostępność). Żadna funkcja z `mic_logic.py` (`process_mic_wizualny`, `process_mic_od`, `process_mbc`, `aggregate_all`, `compare_mic_groups`, `compute_mbc_mic_ratio`, `build_mic_summary_rows`...) ani z `mic_plotting.py` nie jest nigdzie w `gui.py` wywoływana. Zweryfikowałem to przez `grep` całego pliku — zero trafień poza `route_workbook`. Interfejs pokazuje tylko znacznik dostępności (`✓/✗` per szczep) w logu/etykiecie; gdy szczep ma tylko dane MIC/MBC (brak dyfuzji), GUI czyści panel i wypisuje "Analiza dyfuzji jest dla niego niedostępna" — **nie oferuje żadnej alternatywnej ścieżki**. Nie ma przycisku, zakładki ani menu uruchamiającego analizę MIC/MBC, wykresy z `mic_plotting.py` ani sekcję MIC/MBC z `reports.py`.
**Waga:** **KRYTYCZNE** (z punktu widzenia kompletności produktu — nie jest to błąd w istniejącym kodzie MIC/MBC, tylko całkowity brak dostępu do niego).
**Typ:** odpornościowy / kompletności funkcji, nie liczbowy.

### 1.2 — Jeden nienumeryczny wpis w `Stez_S1`/`Wsp_rozc` wywala przetwarzanie CAŁEGO arkusza
**Czego dotyczy:** `mic_logic.py::_process_row` / `_process_mbc_row` (dzielenie `stez_s1 / wsp_rozc**...` bez wcześniejszej walidacji typu), wywoływane przez `process_mic_wizualny` / `process_mic_od` / `process_mbc` (list comprehension po `df.iterrows()` bez obsługi wyjątków per-wiersz).
**Co jest nie tak:** `pd.isna(stez_s1)` nie wykrywa tekstu (np. literówki „nieznane”) — to nie jest NaN, tylko string. Kod przechodzi dalej do `compute_well_concentration`, gdzie `'nieznane' / 2` rzuca `TypeError`. Ponieważ `process_mic_wizualny` woła `_process_row` w pętli **bez `try/except`**, wyjątek z JEDNEGO złego wiersza przerywa całą listę — **wiersze przed i po nim, nawet całkowicie poprawne, nigdy nie zostają przetworzone**.
**Zweryfikowane bezpośrednio:** 3-wierszowy arkusz (wiersz 1 poprawny, wiersz 2 z literówką w `Stez_S1`, wiersz 3 poprawny) → `process_mic_wizualny` rzuca `TypeError: unsupported operand type(s) for /: 'str' and 'int'` i **nic** nie zostaje zwrócone, łącznie z wierszami 1 i 3.
**Kontrast z modułem dyfuzji:** ten sam rodzaj błędu (tekst w kolumnie liczbowej, `Srednica_mm='literowka'`) w `utils.validate_and_normalize` jest obsłużony **poprawnie i bezpiecznie** — wiersz trafia na listę `rejected` z czytelnym powodem, a pozostałe 2 poprawne wiersze są przetworzone normalnie. Zweryfikowałem to bezpośrednio (patrz sekcja 4).
**Waga:** **KRYTYCZNE**.
**Typ:** odpornościowy (crash), nie cichy błąd liczbowy — ale efekt jest gorszy niż cichy błąd: cała analiza się nie wykonuje.

### 1.3 — Brak lub niespójny typ `Rep_biologiczna` cicho scala (albo cicho dzieli) prawdziwe powtórzenia biologiczne
**Czego dotyczy:** `mic_logic.py::aggregate_all` (grupowanie `defaultdict(list)` po kluczu `(Bakteria, Substancja, r.get("Rep_biologiczna"))`), brak jakiejkolwiek walidacji/domyślnej wartości tej kolumny w ścieżce MIC/MBC (w przeciwieństwie do `utils.build_internal_representation`, patrz 4.1).
**Co jest nie tak:**
- Gdy kolumna `Rep_biologiczna` **nie istnieje w ogóle** w arkuszu MIC_wizualny/MIC_OD/MBC_posiew, `row.get(COL_REP_BIO)` zwraca `None` dla **każdego** wiersza. Ponieważ `None is None` jest zawsze prawdą, WSZYSTKIE wiersze danej pary (Bakteria, Substancja) trafiają do jednego klucza grupowania — traktowane jak powtórzenia TECHNICZNE jednego powtórzenia biologicznego, nie jak niezależne powtórzenia biologiczne.
- Gdy kolumna istnieje, ale zawiera **mieszane typy** (część komórek sformatowana jako liczba, część jako tekst — typowy artefakt Excela), np. `1` (int) i `"1"` (str) dla tego samego zamierzonego powtórzenia, tworzą **dwa różne** klucze krotki (`hash(1) != hash("1")`) — jedno prawdziwe powtórzenie biologiczne zostaje sztucznie rozbite na dwa.
**Zweryfikowane bezpośrednio (scenariusz „brak kolumny”):** 3 wiersze z MIC = 512, 32, 8 mg/ml (rozstęp 64×) BEZ kolumny `Rep_biologiczna` → `aggregate_all` zwraca **`n_bio=1`**, mediana **„32 mg/ml”**, zakres **„32 mg/ml – 32 mg/ml”**. Prawdziwy 64-krotny rozstęp znika całkowicie z raportu — wygląda, jakby nie było żadnej zmienności, podczas gdy w rzeczywistości mogła to być ogromna zmienność biologiczna albo zupełnie nieświadomie połączone różne substancje/powtórzenia.
**Zweryfikowane bezpośrednio (scenariusz „mieszany typ”):** `('E. coli','DrugX',1)` i `('E. coli','DrugX','1')` to dwa różne klucze słownika — potwierdzone wprost.
**Waga:** **KRYTYCZNE**.
**Typ:** liczbowy + odpornościowy (cichy błąd — zero wyjątku, zero ostrzeżenia, wynik wygląda na w pełni poprawny).

### 1.4 — Klasyfikacja MBC/MIC (d≤2 / d≥3) zakłada rozcieńczenia 2-krotne; dla innego `Wsp_rozc` dokładny pomiar może wypaść jako „nieoznaczalny”
**Czego dotyczy:** `mic_logic.py::compute_mbc_mic_ratio`, `config.py::MBC_MIC_BACTERICIDAL_MAX_D`/`MBC_MIC_BACTERIOSTATIC_MIN_D`. Kontrast: `mic_logic.py::_process_row`/`_process_mbc_row` sprawdzają tylko `wsp_rozc > 1` — **nie ma żadnego wymogu, by `Wsp_rozc` wynosiło akurat 2**.
**Co jest nie tak:** klasyfikacja bakteriobójcze/bakteriostatyczne jest oparta o **stałe progi na skali log2** (d≤2 ⇒ bakteriobójcze, d≥3 ⇒ bakteriostatyczne), które są bezpieczne WYŁĄCZNIE, gdy MIC i MBC pochodzą z serii 2-krotnych rozcieńczeń (wtedy `d` jest zawsze liczbą całkowitą, więc nigdy nie wpada między 2 a 3). Kod jednak **jawnie dopuszcza dowolny współczynnik rozcieńczenia > 1** (3-krotny, 5-krotny, dziesiętny itd.), a `d = log2(MBC/MIC)` dla takich serii **nie jest liczbą całkowitą** — może wypaść dokładnie w „martwej strefie” (2,3), mimo że pomiar jest w 100% dokładny, bez cenzury.
**Zweryfikowane bezpośrednio:** MIC=100 (dokładne), MBC=500 (dokładne) — jeden krok serii 5-krotnej. `d = log2(5) = 2.321928...`. Program zwraca `classification='nieoznaczalny'`, `status='nieoznaczalny'`, mimo że **nie ma tu żadnej niepewności pomiarowej** — oba wejścia są dokładnymi liczbami. Dla porównania: seria 3-krotna z MIC=100/MBC=300 (`d=log2(3)=1.585`) poprawnie klasyfikuje się jako bakteriobójcze, a MIC=100/MBC=900 (`d=log2(9)=3.170`) poprawnie jako bakteriostatyczne — błąd ujawnia się tylko w konkretnym przedziale współczynników rozcieńczenia (np. w okolicach 5-krotnego), więc łatwo go przeoczyć testując wyłącznie serie 2- i 3-krotne.
**Waga:** **ISTOTNE** (nie krytyczne, bo w praktyce mikrorozcieńczeń 2-krotne serie są zdecydowanie dominującym standardem — ale kod NIE wymusza tego założenia, mimo że cała klasyfikacja go zakłada).
**Typ:** liczbowy + metodologiczny.

### 1.5 — Zero walidacji spójności jednostek (`Jednostka`) na jakimkolwiek etapie
**Czego dotyczy:** `mic_logic.py::aggregate_technical_to_biological`, `summarize_mic_group`, cała ścieżka log2 (mediana, zakres, iloraz MBC/MIC, testy rangowe).
**Co jest nie tak:** wszystkie te funkcje liczą na `math.log2(mic_value)`, gdzie `mic_value` to surowa liczba — jednostka (`Jednostka`, np. „mg/ml” vs „µg/ml”, różnica 1000×) jest zapisywana WYŁĄCZNIE do wyświetlenia (`unit = first.get("Jednostka")` — bierze wartość z PIERWSZEGO wiersza, bez porównania z resztą) i nigdy nie jest sprawdzana pod kątem spójności.
**Zweryfikowane bezpośrednio:** dwa „powtórzenia techniczne” tego samego powtórzenia biologicznego, jedno zadeklarowane jako `50 mg/ml`, drugie jako `50 µg/ml` (rzeczywista różnica 1000×) → `aggregate_technical_to_biological` zwraca `mic_value=50.0`, `unit='mg/ml'` — **bez jakiegokolwiek ostrzeżenia o niespójności**, a raportowana jednostka („mg/ml”, wzięta z PIERWSZEGO wiersza) niekoniecznie nawet odpowiada wybranej wartości (w tym przykładzie wybrana została technicznie wartość z DRUGIEGO wiersza, który był w µg/ml).
**Waga:** **ISTOTNE**.
**Typ:** liczbowy + odpornościowy (cichy błąd).

### 1.6 — Niespójność międzymodułowa: puste grupy w porównaniu — cicho pominięte (dyfuzja) vs jawnie zablokowane (MIC)
**Czego dotyczy:** `logic.py::StatsEngine.run_statistics` (pętla `if len(data) >= 1: valid_groups.append(g)` — brak gałęzi `else`, brak listy wykluczonych grup) vs `mic_logic.py::_check_layer3_guards` (`blocked_groups` z jawnym, nazwanym powodem).
**Co jest nie tak:** gdy w module dyfuzji jedna z porównywanych grup ma `n=0` (pusta), zostaje ona po prostu wykluczona z `dane_list` bez żadnej wzmianki w zwracanym wyniku, KOMU dotyczy wykluczenie — jeśli pozostaje ≥2 grupy, analiza biegnie dalej tak, jakby ta grupa nigdy nie istniała, bez ostrzeżenia w raporcie/PDF/logu. Moduł MIC w analogicznej sytuacji (grupa o `n_bio=0`) **blokuje całe porównanie** i zwraca jawny, nazwany powód (`"Porównanie zablokowane: grupa(y) bez żadnej wartości MIC (n_bio=0): ..."`). To dwa różne standardy przejrzystości dla koncepcyjnie tej samej sytuacji w tej samej aplikacji.
**Waga:** **UMIARKOWANE**.
**Typ:** spójność międzymodułowa (nie jest to błąd liczbowy sam w sobie, ale różna filozofia zgłaszania tego samego problemu).

### 1.7 — `Przebieg` dopasowywane globalnie między arkuszami — brak dokumentacji zakresu unikalności
**Czego dotyczy:** `mic_logic.py::lookup_controls` (`controls_df[controls_df[COL_RUN] == przebieg]` — dopasowanie po samej wartości, bez informacji z jakiego arkusza pochodzi wiersz wywołujący); `config.py` (opis `COL_RUN` nie precyzuje zakresu unikalności).
**Co jest nie tak:** `Przebieg` z arkusza `MIC_wizualny` i `Przebieg` z arkusza `MBC_posiew` są dopasowywane do TEGO SAMEGO wiersza w `Kontrole`, jeśli mają tę samą wartość tekstową — co jest prawdopodobnie ZAMIERZONE (MBC zwykle pochodzi z posiewu tych samych studzienek MIC, więc współdzielenie identyfikatora przebiegu ma sens), ale nigdzie nie jest to jawnie udokumentowane jako wymóg. Użytkownik, który (rozsądnie, ale błędnie) numeruje `Przebieg` osobno w obrębie każdego arkusza (np. zaczynając znowu od "1" w MBC_posiew), **cicho** dopasuje kontrolę z zupełnie niepowiązanego przebiegu.
**Waga:** **UMIARKOWANE**.
**Typ:** dokumentacyjny + odpornościowy (potencjalnie cichy błąd, zależnie od konwencji użytkownika).

### 1.8 — `lookup_controls` cicho bierze pierwszy wiersz przy zduplikowanym `Przebieg` w arkuszu Kontrole
**Czego dotyczy:** `mic_logic.py::lookup_controls` — `row = matches.iloc[0]` bez sprawdzenia `len(matches) > 1`.
**Co jest nie tak:** jeśli arkusz `Kontrole` ma dwa wiersze z tym samym `Przebieg` (błąd danych — np. wklejone dwa razy), funkcja cicho używa pierwszego, ignorując drugi, bez żadnego ostrzeżenia o duplikacie.
**Waga:** **DROBNE**.
**Typ:** odpornościowy (cichy błąd, ale wymaga już błędu w danych źródłowych, żeby się ujawnić).

### 1.9 — Brak zbiorczego, "z góry" ostrzeżenia przy całkowitym braku arkusza Kontrole dla MIC_OD
**Czego dotyczy:** `utils.py::route_workbook` (lista `warnings` — nie zawiera nic o brakującym `Kontrole` przy obecności `MIC_OD`).
**Co jest nie tak:** każdy POJEDYNCZY wiersz `MIC_OD` poprawnie i jawnie zgłasza `MIC_STATUS_MISSING_CONTROLS`, gdy `Kontrole` nie istnieje — to jest bezpieczne. Ale nie ma ŻADNEGO zbiorczego ostrzeżenia na poziomie całego pliku ("masz arkusz MIC_OD, ale nie masz w ogóle arkusza Kontrole — WSZYSTKIE wiersze OD zawiodą"), więc użytkownik dowiaduje się o tym dopiero po uruchomieniu pełnej analizy, wiersz po wierszu, zamiast od razu przy wczytaniu pliku.
**Waga:** **DROBNE** (UX, nie poprawność).
**Typ:** odpornościowy / użyteczności.

### 1.10 — `classify_od_well`'s zabezpieczenie `denom<=0` to martwy kod przy obecnym porządku wywołań
**Czego dotyczy:** `mic_logic.py::classify_od_well` (`if denom <= 0: return None`).
**Co jest nie tak:** w obecnym przepływie (`_process_row`, tryb "od") `validate_run` ZAWSZE wykonuje się przed dotarciem do `classify_od_well`, a `validate_run` odrzuca przebieg, gdy `Kontrola_wzrostu - Kontrola_jalowosci < MIC_OD_MIN_GROWTH_SIGNAL (0.10)` — więc w praktyce `denom` nigdy nie może być ≤0 w momencie wywołania `classify_od_well`. To nie jest błąd (zabezpieczenie jest nieszkodliwe), ale sugeruje niepewność autora co do własnych niezmienników, a nie świadome, zweryfikowane podwójne zabezpieczenie.
**Waga:** **KOSMETYCZNE**.
**Typ:** brak (martwy kod, nieszkodliwy).

### 1.11 — Progi OD (10% wzrostu, 0.10/0.20 OD kontroli) opisane jako „standardowe”, choć nie są przywiązane do konkretnego cytowanego protokołu
**Czego dotyczy:** `config.py::MIC_OD_GROWTH_THRESHOLD`, `MIC_OD_MIN_GROWTH_SIGNAL`, `MIC_OD_STERILITY_MAX`.
**Co jest nie tak:** komentarze nazywają 10% „standardowym, konserwatywnym progiem” dla odczytu MIC metodą OD. CLSI M07 (główny standard mikrorozcieńczeń) definiuje MIC przez **odczyt wizualny/mętność**, nie przez skodyfikowany liczbowy próg procentowy OD — nie ma jednego, powszechnie cytowanego numeru „10%” w tym standardzie. Wartości 10-20% pojawiają się w różnych publikowanych metodach OD-owych, więc wybór 10% jest **rozsądny i mieści się w praktyce**, ale stwierdzenie „standardowy” w komentarzu sugeruje większy autorytet formalny, niż to uzasadnione. To samo dotyczy `MIC_OD_STERILITY_MAX=0.20` i `MIC_OD_MIN_GROWTH_SIGNAL=0.10` — kod sam przyznaje w dalszej części komentarza, że to „pragmatyczny próg narzędzia przesiewowego, nie wartość z walidowanego protokołu klinicznego” (czyli częściowo sam sobie zaprzecza w niuansie sformułowania).
**Waga:** **KOSMETYCZNE** (dotyczy języka dokumentacji, nie liczb).
**Typ:** metodologiczny (dokumentacyjny).

---

## 2. Rzeczy, które wyglądają poprawnie i POTWIERDZIŁEM to niezależnym rachunkiem

| # | Co sprawdzono | Metoda weryfikacji | Wynik |
|---|---|---|---|
| 1 | `compute_well_concentration` (Stez_S1/Wsp_rozc^(n-1)) | Ręczny rachunek dla S1..S10, Stez_S1=1024, Wsp_rozc=2 | Identyczne z programem (1024,512,...,2) |
| 2 | Próg względny OD (`classify_od_well`) | Ręczny rachunek na granicy progu (OD=0.17 vs 0.169, próg=0.10, Kw=0.8, Kj=0.1) | Zgodne z programem dokładnie na granicy (0.10000 → True, 0.09857 → False) |
| 3 | Redukcja CFU i próg 99,9% (`classify_mbc_well`) | Ręczny rachunek na granicy (CFU=199 vs 201, CFU_t0=200000) | Zgodne z programem dokładnie na granicy (0.999005 → spełnia, 0.998995 → nie spełnia) |
| 4 | Mediana log2, N nieparzyste (5) | Ręczny rachunek + przetasowana lista wejściowa | Zwraca poprawny środkowy element (3) niezależnie od kolejności wejścia |
| 5 | Mediana log2, N parzyste (6), bez remisów | Ręczny rachunek | Zwraca WYŻSZY z dwóch środkowych elementów (4), zgodnie z deklarowaną regułą |
| 6 | Remis na tej samej wartości log2 (dokładna vs cenzurowana „≤”) | Ręczna konstrukcja remisu | Program poprawnie preferuje wartość DOKŁADNĄ nad cenzurowaną przy remisie |
| 7 | Test rangowy z cenzurą (Mann-Whitney), w tym przypadek adwersarialny (własna granica cenzury MNIEJSZA niż realne wartości innej grupy) | Ręczne wyliczenie sum rang i U-statystyki od zera, niezależnie od kodu programu | U=9.0 identyczne z programem; potwierdzone, że wpis cenzurowany-górnie poprawnie ląduje na SZCZYCIE połączonej rangi mimo małej własnej granicy |
| 8 | Zaokrąglenie zmiennoprzecinkowe przy granicy klasyfikacji (log2(100)-log2(25)=2.0 dokładnie) | Rachunek bezpośredni w Pythonie | Zgodne, mechanizm zaokrąglania do 9 miejsc działa poprawnie tam, gdzie jest potrzebny |
| 9 | Wsteczna zgodność — `dane_disk.xlsz`, kod SPRZED całego modułu MIC/MBC vs kod PO scaleniu | Wyodrębniłem historyczne wersje `config.py`/`utils.py`/`logic.py` z commita `65a2a83~1` (bezpośrednio przed rozpoczęciem gałęzi mic-mbc) do izolowanego katalogu i uruchomiłem identyczny skrypt względem OBU wersji | **Identyczny hash MD5 średnich per grupa** (E. coli: `d8d28b60...`, S. aureus: `181053f0...`) i **identyczna statystyka/p-value Kruskala-Wallisa** (87.6948545564.../8.205387760068e-08 oraz 84.7044934112.../2.324498350703e-07) przed i po |
| 10 | Wartości krytyczne testu Q Dixona (n=3..10, α=0.10) | Porównanie z powszechnie cytowaną tabelą (Rorabacher 1991) | Zgodne |
| 11 | Progi Cohen's d (0.2/0.5/0.8) | Porównanie z Cohen (1988) | Zgodne ze standardową konwencją |
| 12 | Próg MBC 99,9% (redukcja 3-log10) | Porównanie ze standardową definicją kliniczną MBC (mikrobiologia kliniczna, protokoły mikrorozcieńczeń) | Zgodne — to jedna z niewielu stałych w projekcie, która JEST faktycznym standardem, nie autorską konwencją |
| 13 | Próg „sensowności metodologicznej” (≥2 rozcieńczenia = 4×) dla różnicy median MIC | Porównanie z powszechnie przywoływaną powtarzalnością testu mikrorozcieńczeń (±1 rozcieńczenie) | Rozsądny, broniony wybór, spójny z praktyką — pod warunkiem 2-krotnej serii rozcieńczeń (patrz jednak 1.4 dla NIESPÓJNEGO zastosowania tej samej logiki w klasyfikacji MBC/MIC) |
| 14 | Pula surogatów cenzury liczona RAZEM dla obu/wszystkich porównywanych grup (nie osobno per grupa) | Odczyt kodu `compare_mic_groups` + niezależny test adwersarialny (patrz #7 wyżej) | Potwierdzone - `_censoring_surrogate` jest wywoływane na połączonej próbie, dokładnie jak deklarowane |

---

## 3. Rzeczy wątpliwe lub wymagające Twojej decyzji

1. **Czy wpięcie modułu MIC/MBC do GUI jest zaplanowane jako osobny etap, czy to przeoczenie?** Obecnie cała praca 5 faz jest dostępna tylko programistycznie (import + wywołanie funkcji), nie z poziomu aplikacji. Warto świadomie zdecydować, czy to zamierzony podział pracy (silnik najpierw, UI później) czy luka do zamknięcia przed uznaniem modułu za "gotowy".
2. **Czy `Wsp_rozc` powinno być ograniczone do wartości 2 (lub próg klasyfikacji d powinien skalować się z faktycznym współczynnikiem rozcieńczenia)?** Obecny kod jawnie wspiera dowolne `Wsp_rozc>1`, ale klasyfikacja MBC/MIC zakłada 2-krotność. To wymaga decyzji projektowej, nie tylko poprawki technicznej.
3. **Czy `Rep_biologiczna` powinno być kolumną WYMAGANĄ (z twardym błędem walidacji przy braku/niespójności), zamiast cichego domyślnego zachowania odziedziczonego po starym formacie dyfuzji?** Obecne zachowanie (cichy fallback do `None`/mieszanych typów) jest szczególnie ryzykowne właśnie w formacie MIC/MBC, którego sensem jest rozróżnienie powtórzeń technicznych od biologicznych.
4. **Czy warto dodać walidację spójności jednostek** (przynajmniej ostrzeżenie, gdy technika/biologiczne powtórzenia w jednej grupie mają różne `Jednostka`)? To stosunkowo tania zmiana o dużym znaczeniu dla wiarygodności wyników.
5. **Czy `Przebieg` ma być dokumentowany jako identyfikator GLOBALNY** (spójny między MIC_wizualny/MIC_OD/MBC_posiew), czy per-arkuszowy? Obecne zachowanie kodu (globalne dopasowanie) powinno zostać albo udokumentowane wprost, albo zmienione — na razie jest to milcząca, niejasna konwencja.
6. **Czy próg 10% wzrostu / 0.10 / 0.20 OD powinny zostać przeformułowane w komentarzach jako jawnie autorskie/konfigurowalne** (a nie "standardowe"), żeby przyszły czytelnik kodu nie szukał nieistniejącego źródła w CLSI/EUCAST?

---

## 4. Log weryfikacji (skrócone dowody)

Poniższe polecenia zostały uruchomione podczas audytu (nie zmieniają repozytorium):

```
# 1.2 - crash całego arkusza
process_mic_wizualny([wiersz_OK, wiersz_z_literowka_w_Stez_S1, wiersz_OK])
→ TypeError: unsupported operand type(s) for /: 'str' and 'int'
   (żaden z 3 wierszy, łącznie z 2 poprawnymi, nie zostaje zwrócony)

# kontrast: to samo w module dyfuzji
utils.validate_and_normalize(df z Srednica_mm='literowka' w 1 z 3 wierszy)
→ brak wyjątku; 2 wiersze przetworzone, 1 odrzucony z jawnym powodem

# 1.3 - zniknięcie zmienności przy braku Rep_biologiczna
3 wiersze MIC = 512, 32, 8 mg/ml, BRAK kolumny Rep_biologiczna
→ n_bio=1, mediana="32 mg/ml", zakres="32 mg/ml - 32 mg/ml"

# 1.4 - luka klasyfikacji dla rozcieńczeń niebędących 2-krotnością
MIC=100 (dokładne), MBC=500 (dokładne), d=log2(5)=2.321928...
→ classification='nieoznaczalny' (mimo braku jakiejkolwiek cenzury/niepewności)

# 1.5 - brak walidacji jednostek
2 powtórzenia techniczne: 50 mg/ml i 50 µg/ml (te same liczby, różne jednostki)
→ wybrana wartość 50.0, unit='mg/ml' (z PIERWSZEGO wiersza), zero ostrzeżenia

# wsteczna zgodność - dane_disk.xlsx, kod sprzed vs po module MIC/MBC
hash średnich E. coli:   PRZED=d8d28b60f9495eed5460a83c43aabf8d  PO=d8d28b60f9495eed5460a83c43aabf8d
hash średnich S. aureus: PRZED=181053f0b160e44c1b1eae2c78133683 PO=181053f0b160e44c1b1eae2c78133683
Kruskal-Wallis E. coli:  PRZED H=87.69485455645736 p=8.205387760068333e-08
                         PO    H=87.69485455645736 p=8.205387760068333e-08
Kruskal-Wallis S. aureus: identycznie (84.70449341128518 / 2.3244983507029469e-07)
```
