import os
import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from datetime import datetime

# Impornty modułów
import utils
from dialogs import OutlierDialog, HelpDialog, AboutDialog
import reports
from logic import StatsEngine
from plotting import Plotter
import mic_logic
import mic_plotting
from config import DISC_DIAMETER_MM, ALPHA, COL_GROUP, COL_MEASUREMENT, EXPORT_DPI, REF_PLACEHOLDER

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")

class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("BioStat Master v3.0 - Modular")
        self.after(0, lambda: self.state('zoomed'))
        
        # --- ZMIENNE DANYCH ---
        self.df = None           
        self.col_bact_name = None
        self.checkboxes = []     
        self.sample_vars = {}    
        
        # --- ZMIENNE EXPORTU ---
        self.export_data_raw = None
        self.export_stats_normality = []
        self.export_stats_main = [] 
        self.export_stats_posthoc = None
        self.posthoc_detailed_results = []
        self.stats_summary = None
        self.low_n_bio_warning = False

        # --- ROUTER WIELOARKUSZOWY (MIC/MBC) ---
        self.route = None
        self.availability = {}

        # --- MIC/MBC (Faza integracji GUI - audyt 1.1) ---
        # Wyniki mic_logic.aggregate_all (Faza 3), policzone RAZ przy wczytaniu
        # pliku dla WSZYSTKICH szczepów naraz (tak jak router) - klucz to
        # (Bakteria, Substancja) -> {"bio_results":[...], "summary": {...}}.
        # Całkowicie osobne od self.df/self.stats_summary itd. (dyfuzja) -
        # nic tu nie może wpłynąć na ścieżkę dyfuzji.
        self.mic_grouped = {}
        self.mbc_grouped = {}
        self.mic_mbc_window = None
        self.mic_mbc_figures = {}

        # --- FIGURY ---
        self.figures = {
            'bar': None, 'heat': None, 'pvalue': None,
            'trend': None, 'effect': None, 'cross': None, 'pca': None
        }
        
        # --- KONFIGURACJA ---
        self.plot_config = {
            "font_labels": 10, "font_title": 12, "axis_max": 0, "star_offset": 0.03, "bar_width": 0.8,
            "show_disk_line": True, "palette": "viridis", "transparent_background": True,
            "plot_type": "Barplot (Słupkowy)", "error_bar": "SD (Odchylenie Std.)", "show_points": False,
            # UX: domyślnie tylko porównania vs grupa odniesienia na wykresie
            # wielkości efektu (spójnie z wykresem głównym) - przy wielu
            # grupach wszystkie pary (setki) robiły z niego nieczytelną masę.
            "effect_all_pairs": False,
        }
        
        self.available_palettes = ["viridis", "magma", "plasma", "inferno", "Blues", "Reds", "Greens", "Spectral", "coolwarm", "gray", "tab10"]
        self.available_plot_types = ["Barplot (Słupkowy)", "Boxplot (Pudełkowy)", "Violinplot (Skrzypcowy)"]
        self.available_error_bars = ["SD (Odchylenie Std.)", "SEM (Błąd Std.)", "95% CI (Przedział Ufności)"]

        # --- MODUŁY ---
        self.stats_engine = StatsEngine()
        self.plotter = Plotter(self.plot_config)

        # --- LAYOUT ---
        self._setup_layout()
        self.log("Witaj w wersji 3.0 (Modularnej)! Wczytaj plik Excel.")

    def _setup_layout(self):
        self.grid_columnconfigure(1, weight=1) 
        self.grid_columnconfigure(2, weight=0) 
        self.grid_rowconfigure(0, weight=1)

        # Lewy Panel
        self.sidebar = ctk.CTkFrame(self, width=220, corner_radius=0)
        self.sidebar.grid(row=0, column=0, sticky="nsew")
        self.sidebar.grid_rowconfigure(18, weight=1) 

        self.logo = ctk.CTkLabel(self.sidebar, text="Panel Sterowania", font=ctk.CTkFont(size=18, weight="bold"))
        self.logo.grid(row=0, column=0, padx=20, pady=(20, 10))

        self.btn_load = ctk.CTkButton(self.sidebar, text="1. Wczytaj Excel", command=self.load_file)
        self.btn_load.grid(row=1, column=0, padx=20, pady=10)
        self.lbl_file = ctk.CTkLabel(self.sidebar, text="Brak pliku", text_color="gray", font=("Arial", 10))
        self.lbl_file.grid(row=2, column=0, padx=20, pady=(0, 10))

        self.lbl_bact = ctk.CTkLabel(self.sidebar, text="2. Wybierz szczep:", anchor="w")
        self.lbl_bact.grid(row=3, column=0, padx=20, pady=(10, 0), sticky="w")
        self.combo_bact = ctk.CTkOptionMenu(self.sidebar, values=["..."], command=self.on_bacteria_change)
        self.combo_bact.grid(row=4, column=0, padx=20, pady=(5, 10))

        self.lbl_method = ctk.CTkLabel(self.sidebar, text="3. Korekta Post-hoc:", anchor="w")
        self.lbl_method.grid(row=5, column=0, padx=20, pady=(10, 0), sticky="w")
        self.combo_method = ctk.CTkOptionMenu(self.sidebar, values=["holm", "fdr_bh", "bonferroni", "None"])
        self.combo_method.grid(row=6, column=0, padx=20, pady=(5, 10))
        self.combo_method.set("holm")

        self.lbl_ref = ctk.CTkLabel(self.sidebar, text="4. Grupa odniesienia (*):", anchor="w")
        self.lbl_ref.grid(row=7, column=0, padx=20, pady=(10, 0), sticky="w")
        self.combo_ref = ctk.CTkOptionMenu(self.sidebar, values=["..."])
        self.combo_ref.grid(row=8, column=0, padx=20, pady=(5, 10))
        
        # Orientację przenoszę do configu w przyszłości, na razie zostawiam UI tutaj, ale logika w plotter
        self.lbl_orient = ctk.CTkLabel(self.sidebar, text="5. Orientacja wykresu:", anchor="w")
        self.lbl_orient.grid(row=9, column=0, padx=20, pady=(10, 0), sticky="w")
        self.seg_orient = ctk.CTkSegmentedButton(self.sidebar, values=["Pionowa", "Pozioma"], command=self.update_orientation)
        self.seg_orient.grid(row=10, column=0, padx=20, pady=(5, 10))
        self.seg_orient.set("Pozioma")
        # Inicjalizacja w plotterze
        self.plot_config["orientation"] = "Pozioma"

        self.btn_settings = ctk.CTkButton(self.sidebar, text="⚙ Opcje Wykresu", fg_color="#3B8ED0", command=self.open_plot_settings)
        self.btn_settings.grid(row=11, column=0, padx=20, pady=(20, 10))

        self.btn_run = ctk.CTkButton(self.sidebar, text="URUCHOM ANALIZĘ", fg_color="green", hover_color="darkgreen", 
                                     height=40, font=ctk.CTkFont(size=14, weight="bold"), command=self.run_analysis)
        self.btn_run.grid(row=12, column=0, padx=20, pady=(20, 10), sticky="s")

        ctk.CTkFrame(self.sidebar, height=2, fg_color="gray").grid(row=13, column=0, sticky="ew", padx=10, pady=10)
        
        self.btn_save_plot = ctk.CTkButton(self.sidebar, text="📷 Zapisz Wykres (HQ)", fg_color="#E59400", hover_color="#B37400", command=self.save_plot_image)
        self.btn_save_plot.grid(row=14, column=0, padx=20, pady=5)

        self.btn_export_excel = ctk.CTkButton(self.sidebar, text="💾 Eksportuj do Excela", fg_color="#1F6AA5", command=self.export_to_excel)
        self.btn_export_excel.grid(row=15, column=0, padx=20, pady=(5, 5))

        self.btn_export_pdf = ctk.CTkButton(self.sidebar, text="📄 Generuj Raport PDF", fg_color="#8B0000", hover_color="#600000", command=self.generate_pdf_report)
        self.btn_export_pdf.grid(row=16, column=0, padx=20, pady=(5, 5))
        
        self.btn_captions = ctk.CTkButton(self.sidebar, text="📝 Generuj Opisy Rycin", fg_color="#555555", hover_color="#333333", command=self.open_caption_window)
        self.btn_captions.grid(row=17, column=0, padx=20, pady=(5, 5))

        self.btn_help = ctk.CTkButton(self.sidebar, text="❓ Podręcznik / Pomoc", fg_color="transparent", border_width=1, text_color=("gray10", "gray90"), command=self.open_help_window)
        self.btn_help.grid(row=18, column=0, padx=20, pady=(10, 5), sticky="s")
        
        self.btn_about = ctk.CTkButton(self.sidebar, text="ℹ O twórcy", fg_color="transparent", text_color="gray", font=("Arial", 10), hover_color="#EEE", command=self.open_about_window)
        self.btn_about.grid(row=19, column=0, padx=20, pady=(0, 20), sticky="s")

        # Środkowy Panel
        self.main_view = ctk.CTkTabview(self)
        self.main_view.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        self.tab_plot = self.main_view.add("Wykres Główny")
        self.tab_heatmap = self.main_view.add("Mapa Ciepła")
        self.tab_pvalue = self.main_view.add("Mapa P-value")
        self.tab_trend = self.main_view.add("Trend (Dawka)")
        self.tab_effect = self.main_view.add("Wielkość Efektu")
        self.tab_cross = self.main_view.add("Porównanie Szczepów") 
        self.tab_pca = self.main_view.add("Analiza PCA")
        self.tab_log = self.main_view.add("Raport Statystyczny")
        
        self.textbox = ctk.CTkTextbox(self.tab_log, font=("Consolas", 12))
        self.textbox.pack(expand=True, fill="both", padx=5, pady=5)

        # Prawy Panel
        self.right_frame = ctk.CTkFrame(self, width=200, corner_radius=0)
        self.right_frame.grid(row=0, column=2, sticky="nsew")
        self.right_frame.grid_rowconfigure(3, weight=1)

        # Dostępność analiz per szczep (nad "Wybór próbek")
        ctk.CTkLabel(self.right_frame, text="Dostępność analiz:", font=ctk.CTkFont(weight="bold")).grid(row=0, column=0, padx=10, pady=(10, 0))
        self.lbl_availability = ctk.CTkLabel(
            self.right_frame, text="Wczytaj plik, aby zobaczyć dostępność.",
            justify="left", anchor="w", font=("Consolas", 11), wraplength=180,
        )
        self.lbl_availability.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="w")

        ctk.CTkLabel(self.right_frame, text="Wybór próbek:", font=ctk.CTkFont(weight="bold")).grid(row=2, column=0, padx=10, pady=(10, 5))
        self.scroll_samples = ctk.CTkScrollableFrame(self.right_frame, label_text="Dostępne grupy", height=600)
        self.scroll_samples.grid(row=3, column=0, padx=10, pady=5, sticky="nsew")
        self.btn_select_all = ctk.CTkButton(self.right_frame, text="Zaznacz wszystko", width=100, command=self.select_all)
        self.btn_select_all.grid(row=4, column=0, padx=10, pady=5)
        self.btn_deselect_all = ctk.CTkButton(self.right_frame, text="Odznacz wszystko", width=100, fg_color="gray", command=self.deselect_all)
        self.btn_deselect_all.grid(row=5, column=0, padx=10, pady=(5, 20))

        # --- Sekcja MIC/MBC (osobna od "Wybór próbek" powyżej - dyfuzja) ---
        ctk.CTkFrame(self.right_frame, height=2, fg_color="gray").grid(row=6, column=0, sticky="ew", padx=10, pady=(0, 10))
        ctk.CTkLabel(self.right_frame, text="Analiza MIC/MBC:", font=ctk.CTkFont(weight="bold")).grid(row=7, column=0, padx=10, pady=(0, 5))
        self.lbl_mic_mbc_status = ctk.CTkLabel(
            self.right_frame, text="Wczytaj plik, aby zobaczyć dostępność.",
            justify="left", anchor="w", font=("Consolas", 11), wraplength=180,
        )
        self.lbl_mic_mbc_status.grid(row=8, column=0, padx=10, pady=(0, 5), sticky="w")
        self.btn_run_mic_mbc = ctk.CTkButton(
            self.right_frame, text="🧫 Uruchom analizę MIC/MBC", fg_color="#2E7D32", hover_color="#1B5E20",
            state="disabled", command=self.open_mic_mbc_window,
        )
        self.btn_run_mic_mbc.grid(row=9, column=0, padx=10, pady=(0, 20))

    # ==================== LOGIKA POMOCNICZA ====================
    def log(self, text):
        self.textbox.insert("end", text + "\n")
        self.textbox.see("end")
    def clear_log(self): self.textbox.delete("1.0", "end")

    def load_file(self):
        path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx *.xls")])
        if not path:
            return
        self._load_from_path(path)

    def _load_from_path(self, path):
        """Rdzeń wczytywania wydzielony z load_file, żeby dało się go wywołać
        (np. w testach) bez przechodzenia przez natywne okno dialogowe."""
        try:
            route = utils.route_workbook(path)
        except FileNotFoundError:
            messagebox.showerror("Błąd", "Plik nie istnieje lub został usunięty.")
            return
        except PermissionError:
            messagebox.showerror("Błąd", "Nie można otworzyć pliku. Sprawdź, czy nie jest otwarty w Excelu.")
            return
        except ValueError:
            messagebox.showerror("Błąd", "Nieobsługiwany format pliku. Wymagany: .xlsx lub .xls.")
            return
        except Exception as e:
            messagebox.showerror("Błąd", f"Nie udało się wczytać pliku: {e}")
            return

        if route["errors"]:
            messagebox.showerror("Brak danych do analizy", "\n".join(route["errors"]))
            return

        # --- Ścieżka dyfuzji (dotychczasowa logika, bez zmian) - tylko gdy
        # router wykrył obecny i niepusty arkusz Dane_dyfuzja / stary format ---
        df_internal, bacteria_col, format_info = None, None, {}
        if route["diffusion_raw_df"] is not None:
            df_internal, bacteria_col, rejected, struct_errors, format_info = utils.validate_and_normalize(route["diffusion_raw_df"])

            if struct_errors:
                body = "Arkusz dyfuzji nie może zostać wczytany z powodu następujących problemów:\n\n" + "\n".join(struct_errors)
                messagebox.showerror("Błąd walidacji pliku", body)
                return

            total_rows = len(route["diffusion_raw_df"])
            if rejected:
                if len(df_internal) == 0:
                    messagebox.showerror(
                        "Brak poprawnych danych",
                        f"Wczytano {total_rows} wierszy dyfuzji. Wszystkie zostały odrzucone podczas walidacji — brak danych do analizy."
                    )
                    return
                max_show = 15
                shown_lines = [f"Wiersz {n}: {reason}" for n, reason in rejected[:max_show]]
                if len(rejected) > max_show:
                    shown_lines.append(f"... oraz {len(rejected) - max_show} więcej.")
                body = (
                    f"Wczytano {total_rows} wierszy dyfuzji. Odrzucono {len(rejected)}:\n\n"
                    + "\n".join(shown_lines)
                    + f"\n\nPozostało {len(df_internal)} wierszy do analizy."
                )
                messagebox.showwarning("Odrzucono nieprawidłowe wiersze", body)

        try:
            self.route = route
            self.availability = route["availability"]
            self.df = df_internal
            self.col_bact_name = bacteria_col
            self.lbl_file.configure(text=os.path.basename(path), text_color="white")

            # --- Przetwarzanie MIC/MBC (audyt 1.1 - wpięcie do GUI) ---
            # Całkowicie OSOBNE od ścieżki dyfuzji powyżej: reużywa wprost
            # mic_logic.py (Fazy 2-3, bez zmian logiki), na WSZYSTKICH
            # szczepach naraz (tak jak router) - filtrowanie per szczep
            # dzieje się dopiero przy otwieraniu okna MIC/MBC. Błąd tutaj
            # nigdy nie przerywa wczytania pliku dyfuzji - tylko ostrzega i
            # zostawia MIC/MBC niedostępne.
            self.mic_grouped = {}
            self.mbc_grouped = {}
            try:
                mic_row_results = []
                if route["mic_wizualny_raw_df"] is not None:
                    mic_row_results += mic_logic.process_mic_wizualny(route["mic_wizualny_raw_df"], route["controls_raw_df"])
                if route["mic_od_raw_df"] is not None:
                    mic_row_results += mic_logic.process_mic_od(route["mic_od_raw_df"], route["controls_raw_df"])
                if mic_row_results:
                    self.mic_grouped = mic_logic.aggregate_all(mic_row_results)

                if route["mbc_raw_df"] is not None:
                    mbc_row_results = mic_logic.process_mbc(route["mbc_raw_df"], route["controls_raw_df"])
                    if mbc_row_results:
                        self.mbc_grouped = mic_logic.aggregate_all(mbc_row_results)
            except Exception as e:
                self.mic_grouped, self.mbc_grouped = {}, {}
                messagebox.showwarning(
                    "MIC/MBC",
                    f"Nie udało się przetworzyć danych MIC/MBC: {e}\nAnaliza dyfuzji nie jest tym dotknięta."
                )

            bacts = sorted(self.availability.keys())
            self.combo_bact.configure(values=bacts)
            self.combo_bact.set(bacts[0])
            self.on_bacteria_change(bacts[0])

            self.clear_log()
            sheets_found = [s for s, present in route["sheets_present"].items() if present]
            self.log(f"Wczytano plik. Arkusze danych znalezione: {sheets_found or '(brak - stary jednoarkuszowy format)'}")
            self._log_availability_map()
            for w in route["warnings"]:
                self.log(f"UWAGA: {w}")
            if route["diffusion_raw_df"] is not None:
                is_new_format = format_info.get('has_type') or format_info.get('has_conc') or format_info.get('has_reps')
                fmt_label = "NOWY (Typ/Stężenie/Powtórzenia)" if is_new_format else "STARY"
                self.log(f"Dane dyfuzji: format {fmt_label}.")
        except Exception as e: messagebox.showerror("Błąd", f"Nie udało się wczytać: {e}")

    def _format_availability_line(self, strain):
        a = self.availability.get(strain, {})
        mark = lambda ok: "✓" if ok else "✗"
        return f"{strain}: dyfuzja {mark(a.get('dyfuzja'))}, MIC {mark(a.get('mic'))}, MBC {mark(a.get('mbc'))}"

    def _log_availability_map(self):
        self.log("Dostępność analiz per szczep:")
        for strain in sorted(self.availability.keys()):
            self.log(f"  {self._format_availability_line(strain)}")
        lines = [self._format_availability_line(s) for s in sorted(self.availability.keys())]
        self.lbl_availability.configure(text="\n".join(lines) if lines else "Brak danych.")

    def on_bacteria_change(self, selected_bact):
        avail = self.availability.get(selected_bact, {})
        self._update_mic_mbc_panel(selected_bact, avail)
        if not avail.get('dyfuzja', False):
            # Ten szczep ma dane tylko dla MIC i/lub MBC (albo router go w ogóle
            # nie widział w arkuszu dyfuzji) - analiza dyfuzji jest dla niego
            # niedostępna, więc czyścimy panel zamiast pokazywać mylące
            # "niejednoznaczna referencja".
            for cb in self.checkboxes: cb.destroy()
            self.checkboxes = []
            self.sample_vars = {}
            self.combo_ref.configure(values=["..."])
            self.combo_ref.set("...")
            self.log(
                f"Szczep '{selected_bact}': brak danych z testu dyfuzji krążkowej w tym pliku "
                f"(dostępne: MIC {'✓' if avail.get('mic') else '✗'}, MBC {'✓' if avail.get('mbc') else '✗'}). "
                f"Analiza dyfuzji jest dla niego niedostępna."
            )
            return

        if self.df is None: return
        try:
            all_groups = sorted(self.df[COL_GROUP].unique(), key=utils.smart_sort_key)
            
            for cb in self.checkboxes: cb.destroy()
            self.checkboxes = []
            self.sample_vars = {}
            for g in all_groups:
                var = ctk.IntVar(value=1)
                self.sample_vars[g] = var
                cb = ctk.CTkCheckBox(self.scroll_samples, text=g, variable=var)
                cb.pack(anchor="w", padx=5, pady=2)
                self.checkboxes.append(cb)
            
            df_temp = self.df[self.df[self.col_bact_name] == selected_bact]
            grupy_bact = sorted(df_temp[COL_GROUP].unique(), key=utils.smart_sort_key)

            ref_candidate, ambiguous = utils.select_reference_group(df_temp)
            if not ambiguous:
                self.combo_ref.configure(values=grupy_bact)
                self.combo_ref.set(ref_candidate)
            else:
                # Nigdy nie zgadujemy (np. wybierając kontrolę pozytywną) - użytkownik
                # musi wybrać ręcznie, dopóki nie ma jednoznacznej kontroli negatywnej.
                self.combo_ref.configure(values=[REF_PLACEHOLDER] + grupy_bact)
                self.combo_ref.set(REF_PLACEHOLDER)
                messagebox.showwarning(
                    "Wybierz grupę referencyjną ręcznie",
                    "Nie znaleziono jednoznacznej kontroli negatywnej wśród grup dla tego szczepu "
                    "(0 lub więcej niż 1 pasująca grupa).\n\n"
                    "Wybierz RĘCZNIE grupę referencyjną (kontrolę) w polu '4. Grupa odniesienia' "
                    "przed uruchomieniem analizy. Program nigdy nie wybiera automatycznie kontroli "
                    "pozytywnej (antybiotyku) jako referencji."
                )
                self.log("!!! UWAGA: nie wybrano automatycznie grupy referencyjnej - wybierz ją ręcznie w polu '4. Grupa odniesienia'. !!!")
        except Exception as e: self.log(f"Błąd zmiany bakterii: {e}")

    def select_all(self):
        for var in self.sample_vars.values(): var.set(1)
    def deselect_all(self):
        for var in self.sample_vars.values(): var.set(0)
    def get_selected_groups(self):
        return [g for g, var in self.sample_vars.items() if var.get() == 1]
    
    def update_orientation(self, value):
        self.plot_config["orientation"] = value
        self.plotter.update_config(self.plot_config)
        self.log(f"Zmieniono orientację na: {value}. Odśwież wykres.")

    # ==================== OKNA ====================
    def open_help_window(self):
        if hasattr(self, 'help_window') and self.help_window is not None and self.help_window.winfo_exists():
            self.help_window.lift()
        else: self.help_window = HelpDialog(self)

    def open_about_window(self):
        if hasattr(self, 'about_window') and self.about_window is not None and self.about_window.winfo_exists():
            self.about_window.lift()
        else: self.about_window = AboutDialog(self)

    def open_caption_window(self):
        win = ctk.CTkToplevel(self)
        win.title("Generator Opisów do Publikacji")
        win.geometry("700x600")
        win.attributes("-topmost", True)

        text_area = ctk.CTkTextbox(win, font=("Arial", 12), wrap="word")
        text_area.pack(fill="both", expand=True, padx=10, pady=10)

        bact = self.combo_bact.get()
        ref_group = self.combo_ref.get()
        used_correction_raw = self.combo_method.get()
        if used_correction_raw == "None": correction_desc = "no correction"
        elif used_correction_raw == "fdr_bh": correction_desc = "Benjamini-Hochberg (FDR) correction"
        elif used_correction_raw == "holm": correction_desc = "Holm-Bonferroni correction"
        else: correction_desc = "Bonferroni correction"

        # Audyt (Znalezisko 2): nazwa testu post-hoc musi odpowiadać temu, co
        # logic.py FAKTYCZNIE wykonuje, nie ślepo temu, co wybrano w
        # combo_method - gałąź ANOVA w run_statistics ZAWSZE używa Tukey's
        # HSD (pairwise_tukeyhsd), niezależnie od wyboru w tym dropdownie;
        # `method`/korekta z combo_method dotyczy WYŁĄCZNIE gałęzi
        # Kruskal-Wallis (posthoc_dunn). Ten sam podział testu ANOVA/Kruskal
        # jest już poprawnie zaimplementowany w dialogs.py (HelpDialog,
        # sekcja "Automatyczny opis") - tu używamy tego samego źródła prawdy
        # (export_stats_main[0]["Test"]), żeby nie rozjeżdżać się z nim.
        test_name = "Statistical test"
        posthoc_clause = ", followed by an appropriate post-hoc test"
        posthoc_matrix_desc = "an appropriate post-hoc"
        if self.export_stats_main:
            used_test = self.export_stats_main[0].get("Test", "")
            if "ANOVA" in used_test:
                test_name = "One-way ANOVA"
                posthoc_clause = ", followed by Tukey's HSD post-hoc test"
                posthoc_matrix_desc = "Tukey's HSD"
            elif "Kruskal" in used_test:
                test_name = "Kruskal-Wallis test"
                posthoc_clause = f", followed by Dunn's post-hoc test with {correction_desc}"
                posthoc_matrix_desc = f"Dunn's post-hoc test with {correction_desc}"

        err_conf = self.plot_config["error_bar"]
        if "SD" in err_conf: err_desc = "standard deviation (SD)"
        elif "SEM" in err_conf: err_desc = "standard error of the mean (SEM)"
        else: err_desc = "95% confidence interval (95% CI)"

        plot_type = self.plot_config["plot_type"]
        if "Barplot" in plot_type: 
            viz_desc = "Bars represent the mean inhibition zone diameter"
        else: 
            viz_desc = "Boxplots represent the median and interquartile range (IQR), with whiskers extending to the minimum and maximum values"

        captions = f"""--- OPISY RYCIN (Scientific Captions) ---\n
Możesz skopiować poniższe opisy bezpośrednio do manuskryptu (Word/LaTeX).

=== Rycina 1: Wykres Główny ===
Figure 1. Antibacterial activity of tested samples against {bact}.
{viz_desc}. Error bars indicate the {err_desc} of independent replicates.
Statistical significance was determined using {test_name}{posthoc_clause} for multiple comparisons.
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
P-values were adjusted for multiple comparisons using {posthoc_matrix_desc}.

=== Rycina 5: Trend Dawka-Odpowiedź ===
Figure 5. Dose-response relationship of antibacterial activity.
Lines represent the trend of inhibition zone diameter (mm) across increasing concentrations of tested substances.
Shaded areas indicate the confidence interval. Spearman correlation coefficients (r) are provided for each substance to quantify the strength of the monotonic relationship.

=== Rycina 6: Porównanie Szczepów ===
Figure 6. Cross-species comparison of antibacterial activity.
Bar chart summarizing the mean inhibition zone diameters for selected substances across different bacterial strains.
Error bars represent standard deviation. This overview highlights the differential susceptibility of tested pathogens to the antimicrobial agents.
"""
        text_area.insert("0.0", captions)

    def open_plot_settings(self):
        self.settings_win = ctk.CTkToplevel(self)
        self.settings_win.title("Ustawienia Wykresu")
        self.settings_win.geometry("400x800")
        self.settings_win.attributes("-topmost", True) 
        
        ctk.CTkLabel(self.settings_win, text="Typ wykresu:").pack(pady=(10,5))
        self.option_plot_type = ctk.CTkOptionMenu(self.settings_win, values=self.available_plot_types)
        self.option_plot_type.set(self.plot_config["plot_type"])
        self.option_plot_type.pack(pady=5)

        ctk.CTkLabel(self.settings_win, text="Rodzaj słupka błędu:").pack(pady=(10,5))
        self.option_error = ctk.CTkOptionMenu(self.settings_win, values=self.available_error_bars)
        self.option_error.set(self.plot_config["error_bar"])
        self.option_error.pack(pady=5)

        ctk.CTkLabel(self.settings_win, text="Styl kolorystyczny:").pack(pady=(10,5))
        self.option_palette = ctk.CTkOptionMenu(self.settings_win, values=self.available_palettes)
        self.option_palette.set(self.plot_config["palette"])
        self.option_palette.pack(pady=5)

        self.switch_points = ctk.CTkSwitch(self.settings_win, text="Pokaż punkty pomiarowe")
        if self.plot_config["show_points"]: self.switch_points.select()
        else: self.switch_points.deselect()
        self.switch_points.pack(pady=10)

        self.switch_line = ctk.CTkSwitch(self.settings_win, text=f"Pokaż linię krążka ({DISC_DIAMETER_MM:g}mm)")
        if self.plot_config["show_disk_line"]: self.switch_line.select()
        else: self.switch_line.deselect()
        self.switch_line.pack(pady=10)

        self.switch_trans = ctk.CTkSwitch(self.settings_win, text="Zapisz z przezroczystym tłem")
        if self.plot_config["transparent_background"]: self.switch_trans.select()
        else: self.switch_trans.deselect()
        self.switch_trans.pack(pady=10)

        self.switch_effect_all_pairs = ctk.CTkSwitch(
            self.settings_win, text="Wielkość efektu: wszystkie pary (zamiast vs referencja)"
        )
        if self.plot_config["effect_all_pairs"]: self.switch_effect_all_pairs.select()
        else: self.switch_effect_all_pairs.deselect()
        self.switch_effect_all_pairs.pack(pady=10)

        ctk.CTkLabel(self.settings_win, text="Wielkość etykiet osi:").pack(pady=(5,5))
        self.slider_font_labels = ctk.CTkSlider(self.settings_win, from_=6, to=20, number_of_steps=14)
        self.slider_font_labels.set(self.plot_config["font_labels"])
        self.slider_font_labels.pack(pady=5)

        ctk.CTkLabel(self.settings_win, text="Wielkość tytułów:").pack(pady=(5,5))
        self.slider_font_title = ctk.CTkSlider(self.settings_win, from_=8, to=24, number_of_steps=16)
        self.slider_font_title.set(self.plot_config["font_title"])
        self.slider_font_title.pack(pady=5)

        ctk.CTkLabel(self.settings_win, text="Maks zakres osi (0=auto):").pack(pady=(5,5))
        self.entry_axis_max = ctk.CTkEntry(self.settings_win)
        self.entry_axis_max.insert(0, str(self.plot_config["axis_max"]))
        self.entry_axis_max.pack(pady=5)

        ctk.CTkLabel(self.settings_win, text="Odległość gwiazdki:").pack(pady=(5,5))
        self.slider_star_offset = ctk.CTkSlider(self.settings_win, from_=0.01, to=0.2)
        self.slider_star_offset.set(self.plot_config["star_offset"])
        self.slider_star_offset.pack(pady=5)

        ctk.CTkButton(self.settings_win, text="Odśwież Wykres", fg_color="green", command=self.apply_settings).pack(pady=30)

    def apply_settings(self):
        self.plot_config["plot_type"] = self.option_plot_type.get()
        self.plot_config["error_bar"] = self.option_error.get() 
        self.plot_config["palette"] = self.option_palette.get()
        self.plot_config["show_disk_line"] = bool(self.switch_line.get())
        self.plot_config["show_points"] = bool(self.switch_points.get()) 
        self.plot_config["transparent_background"] = bool(self.switch_trans.get())
        self.plot_config["effect_all_pairs"] = bool(self.switch_effect_all_pairs.get())
        self.plot_config["font_labels"] = int(self.slider_font_labels.get())
        self.plot_config["font_title"] = int(self.slider_font_title.get())
        self.plot_config["star_offset"] = float(self.slider_star_offset.get())
        
        try: 
            val_text = self.entry_axis_max.get()
            if val_text.strip():
                self.plot_config["axis_max"] = float(val_text)
            else:
                 self.plot_config["axis_max"] = 0
        except ValueError: 
            self.plot_config["axis_max"] = 0 
            messagebox.showwarning("Ustawienia", "Nieprawidłowa wartość dla osi (musi być liczbą). Przyjęto auto.") 
        
        self.plotter.update_config(self.plot_config)
        if self.df is not None: self.run_analysis()

    # ==================== GŁÓWNA ANALIZA (REFACTORED) ====================
    def _check_analysis_preconditions(self, df_run, selected_groups):
        """
        Sanity-check pre-analizy. Zwraca True -> kontynuuj, False -> przerwij
        (z powodu twardego bloku albo anulowania przez użytkownika).
        """
        if len(selected_groups) < 2:
            messagebox.showerror(
                "Niewystarczająca liczba grup",
                "Wybrano mniej niż 2 grupy — nie można wykonać porównania statystycznego."
            )
            return False

        sparse = [(g, int((df_run[COL_GROUP] == g).sum())) for g in selected_groups]
        sparse = [(g, n) for g, n in sparse if n < 3]
        if sparse:
            sparse_lines = "\n".join(f"• '{g}' — {n} obserwacji" for g, n in sparse)
            msg = (
                f"Następujące grupy mają mniej niż 3 obserwacje:\n\n{sparse_lines}\n\n"
                "Testy statystyczne będą miały niską moc. Czy kontynuować?"
            )
            if not messagebox.askokcancel("Niska moc statystyczna", msg):
                return False

        total = len(df_run)
        if total < 10:
            msg = (
                f"Łącznie wybrano tylko {total} obserwacji. "
                "Wyniki mogą być niewiarygodne. Czy kontynuować?"
            )
            if not messagebox.askokcancel("Mała próba", msg):
                return False

        return True

    def run_analysis(self):
        if self.df is None: return
        bact = self.combo_bact.get()

        if not self.availability.get(bact, {}).get('dyfuzja', False):
            messagebox.showerror(
                "Brak danych dyfuzji",
                f"Szczep '{bact}' nie ma danych z testu dyfuzji krążkowej w tym pliku.\n\n"
                "Wybierz szczep, dla którego ta analiza jest dostępna (patrz panel 'Dostępność analiz')."
            )
            return

        method = self.combo_method.get()
        ref_group = self.combo_ref.get()
        if method == "None": method = None

        if ref_group in (REF_PLACEHOLDER, "...", ""):
            messagebox.showerror(
                "Brak grupy referencyjnej",
                "Nie wybrano prawidłowej grupy referencyjnej (kontroli).\n\n"
                "Program nie mógł jednoznacznie rozpoznać kontroli negatywnej automatycznie. "
                "Wybierz ją ręcznie w polu '4. Grupa odniesienia' przed uruchomieniem analizy."
            )
            return

        wybrane = self.get_selected_groups()
        if not wybrane:
            messagebox.showwarning("Stop", "Nie wybrano próbek!")
            return

        # 1. Filtrowanie wstępne
        df_run = self.df[
            (self.df[self.col_bact_name] == bact) & 
            (self.df[COL_GROUP].isin(wybrane))
        ].copy()

        if df_run.empty:
            messagebox.showerror(
                "Brak danych dla wybranych grup",
                f"Żadna z wybranych grup nie ma danych dla szczepu '{bact}'.\n\n"
                "Sprawdź, czy wybrane grupy (panel 'Wybór próbek') rzeczywiście dotyczą tego szczepu."
            )
            return

        if not self._check_analysis_preconditions(df_run, wybrane):
            return

        # 2. Outliery (UI Logic)
        outliers_data = utils.find_outliers_dixon(df_run)
        if outliers_data:
            dialog = OutlierDialog(self, outliers_data)
            self.wait_window(dialog) 
            if dialog.result:
                for item in dialog.result:
                    mask = (df_run[COL_GROUP] == item['Group']) & (df_run[COL_MEASUREMENT] == item['Srednica_mm'])
                    idx = df_run[mask].first_valid_index()
                    if idx is not None: df_run = df_run.drop(idx)
                self.log(f"!!! USUNIĘTO {len(dialog.result)} WARTOŚCI ODSTAJĄCYCH !!!")

        self.export_data_raw = df_run

        # 2b. AGREGACJA POWTÓRZEŃ TECHNICZNYCH (Delegacja)
        # Testy istotności, post-hoc i wykresy statystyczne liczą się na
        # średnich BIOLOGICZNYCH (n_bio), nie na surowych wierszach - inaczej
        # powtórzenia techniczne liczyłyby się jako niezależne obserwacje
        # (pseudoreplikacja, zawyżona moc / zaniżone p-value).
        df_bio = utils.aggregate_technical_replicates(df_run, self.col_bact_name)

        # Spójność "puste grupy" (audyt 1.6) - ujednolicone z modułem MIC/MBC
        # (mic_logic._check_layer3_guards): grupa bez ŻADNEJ wartości nigdy
        # nie znika po cichu z porównania - jest jawnie nazwana w komunikacie,
        # nawet jeśli w praktyce zostaje wykluczona (nie ma z niej czego
        # liczyć). NIE blokujemy całej analizy z tego powodu - checkboxy w
        # 'Wybór próbek' są budowane ze WSZYSTKICH grup w całym pliku,
        # niezależnie od wybranego szczepu (żeby zmiana szczepu nie gubiła
        # zaznaczeń), więc "wybrana, ale nieistniejąca dla TEGO szczepu"
        # grupa jest normalnym, częstym stanem, nie błędem - twarda blokada
        # zmuszałaby do ręcznego odznaczania checkboxów przy każdej zmianie
        # szczepu. Grupa może "zniknąć" też z drugiego powodu: usunięcie
        # wartości odstających (Dixon) skasowało WSZYSTKIE jej wiersze.
        present_groups = set(df_bio[COL_GROUP].unique())
        missing_groups = [g for g in wybrane if g not in present_groups]

        self.stats_summary = utils.build_group_summary(df_bio)
        self.low_n_bio_warning = utils.has_low_n_bio(df_bio)

        # 3. STAT ENGINE (Delegacja)
        summary_res, posthoc_df, error = self.stats_engine.run_statistics(df_bio, method, ref_group)

        if error:
            self.log(f"Blad Statystyki: {error}")
            return

        # Logowanie wyników
        self.clear_log()
        self.log(f"=== RAPORT v3: {bact} ===")
        self.log(f">>> GRUPA REFERENCYJNA (kontrola do porównań): {ref_group} <<<")
        self.log(
            "UWAGA: program porównuje średnice stref zahamowania statystycznie i NIE wylicza "
            "klinicznych kategorii S/I/R (Susceptible/Intermediate/Resistant) wg CLSI (M100) ani EUCAST."
        )
        if missing_groups:
            self.log(
                f"UWAGA: pominięto w porównaniu (brak jakichkolwiek danych dla '{bact}', albo "
                f"wszystkie wartości usunięte jako odstające): {', '.join(missing_groups)}"
            )
        if self.low_n_bio_warning:
            self.log("!" * 66)
            self.log("UWAGA: BRAK REPLIKACJI BIOLOGICZNEJ (n_bio<2) dla co najmniej jednej")
            self.log("porównywanej grupy. Poniższe p-value / istotności są WYŁĄCZNIE")
            self.log("orientacyjne - nie są potwierdzone niezależnymi powtórzeniami")
            self.log("biologicznymi (patrz kolumny n_bio/n_tech w tabeli opisowej).")
            self.log("!" * 66)
        self.export_stats_normality = summary_res['normality']
        self.export_stats_main = summary_res['main_stats']
        self.export_stats_posthoc = posthoc_df
        
        # Raportowanie Normalności
        for res in summary_res['normality']:
             # Opcjonalnie loguj wszystko lub tylko problemy
             pass
        if summary_res['all_normal']: self.log(">> Rozkład normalny: TAK")
        else: self.log(">> Rozkład normalny: NIE (użyto testów nieparametrycznych)")
        
        self.log(f"Test Główny: {summary_res['test_used']}")
        if summary_res['main_stats']:
            s = summary_res['main_stats'][0]
            self.log(f"Stat: {s['Statistic']:.2f}, p={s['p-value']:.6f}")

        if summary_res['test_used'] == "ANOVA":
            self.log(
                f"UWAGA: Tukey HSD ma wbudowaną własną korektę wielokrotnych porównań. "
                f"Wybrana metoda post-hoc ('{method}') NIE ma tu zastosowania - dotyczy wyłącznie "
                f"ścieżki Kruskal-Wallis/Dunn."
            )

        # 4. POST HOC DETALE (Delegacja)
        detailed, sig_set = self.stats_engine.process_detailed_results(posthoc_df, df_bio, ref_group, summary_res['test_used'])
        self.posthoc_detailed_results = detailed
        
        if detailed:
            self.log("\n[3] WYNIKI SZCZEGÓŁOWE (Effect Size):")
            for d in detailed:
                if d['Significant']:
                    metrics_d = d["Cohen's d"]
                    if metrics_d != metrics_d:  # NaN - brak replikacji biologicznej w jednej z grup
                        self.log(f"{d['Group 1']} vs {d['Group 2']} | p={d['P-adj']:.4f} | d=nieokreślony (brak replikacji)")
                    else:
                        self.log(f"{d['Group 1']} vs {d['Group 2']} | p={d['P-adj']:.4f} | d={metrics_d:.2f}")

        # 5. RYSOWANIE (Delegacja)
        self.display_plot(lambda: self.plotter.draw_bar_plot(df_bio, bact, ref_group, sig_set, low_n_bio_warning=self.low_n_bio_warning), self.tab_plot, 'bar')
        self.display_plot(lambda: self.plotter.draw_heatmap(df_bio, bact), self.tab_heatmap, 'heat')
        self.display_plot(lambda: self.plotter.draw_pvalue_heatmap(self.export_stats_posthoc, bact), self.tab_pvalue, 'pvalue')

        # Wykres trendu dawka-odpowiedź (na średnich biologicznych)
        fig_trend, err = self.plotter.draw_trend(df_bio, bact)
        if fig_trend:
             self.display_figure(fig_trend, self.tab_trend, 'trend')
        elif err:
             self._show_plot_error(self.tab_trend, err)

        # Porównanie międzygatunkowe - też na średnich biologicznych (cały self.df, wszystkie szczepy)
        df_bio_all = utils.aggregate_technical_replicates(self.df, self.col_bact_name)
        self.display_plot(lambda: self.plotter.draw_cross_species(df_bio_all, self.col_bact_name, wybrane), self.tab_cross, 'cross')
        self.display_plot(
            lambda: self.plotter.draw_effect_plot(
                self.posthoc_detailed_results, ref_group=ref_group,
                show_all_pairs=self.plot_config.get("effect_all_pairs", False),
            ),
            self.tab_effect, 'effect'
        )

        pca_res, pca_err = self.stats_engine.run_pca(self.df, self.col_bact_name, wybrane)
        if pca_res:
             self.display_plot(lambda: self.plotter.draw_pca(pca_res), self.tab_pca, 'pca')
        elif pca_err:
             self._show_plot_error(self.tab_pca, pca_err)

    # ==================== WSPARCIE UI DO RYSOWANIA ====================
    def display_plot(self, draw_func, tab_widget, fig_key):
        """Helper to clear tab, run draw function, and pack canvas."""
        try:
            fig = draw_func()
            self.display_figure(fig, tab_widget, fig_key)
        except Exception as e:
            self._show_plot_error(tab_widget, str(e))

    def display_figure(self, fig, tab_widget, fig_key):
        # 1. Wyczyść tab
        for w in tab_widget.winfo_children(): w.destroy()
        
        if fig is None: return

        # 2. Zapisz ref
        self.figures[fig_key] = fig
        
        # 3. Osadź
        canvas = FigureCanvasTkAgg(fig, master=tab_widget)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    def _show_plot_error(self, tab, msg):
        for w in tab.winfo_children(): w.destroy()
        ctk.CTkLabel(tab, text=f"Błąd wykresu: {msg}").pack(pady=20)

    # ==================== EXPORTY ====================
    def save_plot_image(self):
        current_tab = self.main_view.get()
        fig_to_save = None
        
        if current_tab == "Wykres Główny": fig_to_save = self.figures['bar']
        elif current_tab == "Mapa Ciepła": fig_to_save = self.figures['heat']
        elif current_tab == "Mapa P-value": fig_to_save = self.figures['pvalue']
        elif current_tab == "Trend (Dawka)": fig_to_save = self.figures['trend']
        elif current_tab == "Porównanie Szczepów": fig_to_save = self.figures['cross'] 
        elif current_tab == "Wielkość Efektu": fig_to_save = self.figures['effect']
        elif current_tab == "Analiza PCA": fig_to_save = self.figures['pca']
        
        if fig_to_save is None:
            messagebox.showwarning("Uwaga", "Brak wykresu do zapisania.")
            return

        file_path = filedialog.asksaveasfilename(defaultextension=".png", 
                                                 filetypes=[("PNG Image", "*.png"), ("PDF Document", "*.pdf")])
        if not file_path: return
        try:
            is_transparent = self.plot_config["transparent_background"]
            fig_to_save.savefig(file_path, dpi=EXPORT_DPI, bbox_inches='tight', transparent=is_transparent)
            messagebox.showinfo("Sukces", "Wykres zapisany!")
        except Exception as e: messagebox.showerror("Błąd Zapisu", str(e))

    def export_to_excel(self):
        if self.export_data_raw is None:
            messagebox.showwarning("Uwaga", "Najpierw przeprowadź analizę!")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel file", "*.xlsx")])
        if not file_path: return
        try:
            with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                self.export_data_raw.to_excel(writer, sheet_name="Dane Surowe", index=False)
                if self.stats_summary is not None: self.stats_summary.to_excel(writer, sheet_name="Statystyki Opisowe (n_bio)", index=False)
                if self.export_stats_normality: pd.DataFrame(self.export_stats_normality).to_excel(writer, sheet_name="Normalnosc", index=False)
                if self.export_stats_main: pd.DataFrame(self.export_stats_main).to_excel(writer, sheet_name="Test Glowny", index=False)
                if self.posthoc_detailed_results: pd.DataFrame(self.posthoc_detailed_results).to_excel(writer, sheet_name="Post-hoc (Details)", index=False)
                if self.low_n_bio_warning:
                    pd.DataFrame({"Uwaga": [
                        "BRAK REPLIKACJI BIOLOGICZNEJ (n_bio<2) dla co najmniej jednej porownywanej grupy.",
                        "Wyniki (p-value, istotnosc, Cohen's d) sa WYLACZNIE orientacyjne - nie sa",
                        "potwierdzone niezaleznymi powtorzeniami biologicznymi.",
                        "Zobacz kolumny n_bio/n_tech w arkuszu 'Statystyki Opisowe (n_bio)'.",
                    ]}).to_excel(writer, sheet_name="UWAGA - n_bio", index=False)
            messagebox.showinfo("Sukces", f"Zapisano wyniki w:\n{file_path}")
        except Exception as e: messagebox.showerror("Błąd Zapisu", str(e))

    def generate_pdf_report(self):
        if self.export_data_raw is None:
            messagebox.showwarning("Uwaga", "Najpierw przeprowadź analizę!")
            return
        
        file_path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF Document", "*.pdf")])
        if not file_path: return

        meta = {
            'date': datetime.now().strftime('%Y-%m-%d %H:%M'),
            'bact': self.combo_bact.get(),
            'method': self.combo_method.get(),
            'ref': self.combo_ref.get(),
            'test_used': self.export_stats_main[0]['Test'] if self.export_stats_main else "",
            'low_n_bio_warning': self.low_n_bio_warning,
        }
        
        success, msg = reports.generate_pdf(
            file_path,
            meta,
            self.stats_summary,
            self.figures,
            self.posthoc_detailed_results
        )

        if success:
            messagebox.showinfo("Sukces", msg)
        else:
            messagebox.showerror("Błąd PDF", msg)

    # ==================== MIC/MBC (integracja GUI - audyt 1.1) ====================
    #
    # Ta sekcja WYŁĄCZNIE podłącza istniejącą, przetestowaną logikę
    # (mic_logic.py / mic_plotting.py / reports.py) do interfejsu - nic tu
    # nie liczy niczego od nowa. Całkowicie osobna od sekcji dyfuzji powyżej:
    # własne okno (CTkToplevel), własne zmienne stanu (self.mic_grouped/
    # self.mbc_grouped/self._mic_mbc_state), żadna metoda dyfuzji nie jest
    # tu modyfikowana.

    def _update_mic_mbc_panel(self, selected_bact, avail=None):
        """
        Aktualizuje etykietę/przycisk MIC/MBC w prawym panelu dla wybranego
        szczepu - reużywa WPROST mapę dostępności z routera (Faza 1),
        niczego nie dedukuje samodzielnie. Wołana z on_bacteria_change.
        """
        if avail is None:
            avail = self.availability.get(selected_bact, {})
        has_mic = avail.get('mic', False)
        has_mbc = avail.get('mbc', False)
        mark = lambda ok: "✓" if ok else "✗"
        self.lbl_mic_mbc_status.configure(text=f"{selected_bact}:\nMIC {mark(has_mic)}   MBC {mark(has_mbc)}")
        self.btn_run_mic_mbc.configure(state="normal" if (has_mic or has_mbc) else "disabled")

    def _mic_mbc_method(self):
        """
        Reużywa TĘ SAMĄ korektę post-hoc wybraną dla dyfuzji (combo_method),
        żeby ustawienie było spójne między modułami - z zamianą "None"
        (dyfuzja: brak korekty) na "holm", bo compare_mic_groups (Dunn's
        test) nie ma opcji "brak korekty".
        """
        m = self.combo_method.get()
        return m if m in ("holm", "fdr_bh", "bonferroni") else "holm"

    @staticmethod
    def _filter_grouped_by_strain(grouped, bact):
        return {key: val for key, val in grouped.items() if key[0] == bact}

    def open_mic_mbc_window(self):
        """Punkt wejścia z przycisku '🧫 Uruchom analizę MIC/MBC'."""
        bact = self.combo_bact.get()
        avail = self.availability.get(bact, {})
        if not (avail.get('mic') or avail.get('mbc')):
            messagebox.showerror(
                "Brak danych MIC/MBC",
                f"Szczep '{bact}' nie ma danych MIC ani MBC w tym pliku.\n\n"
                "Wybierz szczep, dla którego ta analiza jest dostępna (patrz panel 'Analiza MIC/MBC')."
            )
            return

        mic_bact = self._filter_grouped_by_strain(self.mic_grouped, bact)
        mbc_bact = self._filter_grouped_by_strain(self.mbc_grouped, bact)
        substances = sorted({s for (_, s) in mic_bact} | {s for (_, s) in mbc_bact})
        if not substances:
            messagebox.showerror("Brak danych MIC/MBC", f"Nie znaleziono żadnej substancji MIC/MBC dla '{bact}'.")
            return

        if self.mic_mbc_window is not None and self.mic_mbc_window.winfo_exists():
            self.mic_mbc_window.destroy()

        win = ctk.CTkToplevel(self)
        win.title(f"Analiza MIC/MBC: {bact}")
        win.geometry("1100x750")
        self.mic_mbc_window = win
        self._bring_window_to_front(win)

        self._mic_mbc_state = {"bact": bact, "mic_bact": mic_bact, "mbc_bact": mbc_bact, "substances": substances}
        self.mic_mbc_figures = {}

        top = ctk.CTkFrame(win)
        top.pack(fill="x", padx=10, pady=(10, 5))
        ctk.CTkLabel(top, text="Substancja:").pack(side="left", padx=(5, 5))
        self.combo_mic_substance = ctk.CTkOptionMenu(
            top, values=substances, command=lambda _v: self._refresh_mic_mbc_distribution()
        )
        self.combo_mic_substance.set(substances[0])
        self.combo_mic_substance.pack(side="left", padx=(0, 20))

        ctk.CTkButton(top, text="💾 Eksportuj Excel (MIC/MBC)", fg_color="#1F6AA5",
                      command=self._export_mic_mbc_excel).pack(side="left", padx=5)
        ctk.CTkButton(top, text="📄 Generuj PDF (MIC/MBC)", fg_color="#8B0000", hover_color="#600000",
                      command=self._export_mic_mbc_pdf).pack(side="left", padx=5)

        tabview = ctk.CTkTabview(win)
        tabview.pack(fill="both", expand=True, padx=10, pady=(5, 10))
        self.mic_mbc_tabview = tabview
        self.tab_mic_dist = tabview.add("Rozkład MIC/MBC")
        self.tab_mic_pairs = tabview.add("Pary MIC↔MBC")
        self.tab_mic_compare = tabview.add("Porównanie substancji")
        self.tab_mic_table = tabview.add("Tabela zbiorcza i ostrzeżenia")

        self._refresh_mic_mbc_distribution()
        self._render_mic_mbc_pairs()
        self._render_mic_mbc_comparison()
        self._render_mic_mbc_table()

        # Budowa zakładek (wykresy) trwa chwilę - po jej zakończeniu okno
        # bywa znowu za głównym oknem na niektórych menedżerach okien
        # Windows, więc podnosimy je jeszcze raz na koniec.
        self._bring_window_to_front(win)

    def _bring_window_to_front(self, win):
        """
        UX (test na realnych danych): okno MIC/MBC otwierało się POD głównym
        oknem aplikacji, nie na wierzchu. Windows/Tk czasem ignoruje samo
        `.lift()` wywołane z callbacku przycisku, jeśli okno wywołujące ma
        fokus - chwilowe wymuszenie "-topmost" i jego natychmiastowe
        wyłączenie to standardowy, niezawodny sposób na wypchnięcie okna na
        wierzch BEZ trwałego przypinania go tam na stałe (to nie ma być
        modalne ani wiecznie "always on top" - użytkownik ma móc normalnie
        przełączyć się z powrotem na główne okno później).
        """
        win.lift()
        win.attributes("-topmost", True)
        win.after(50, lambda: win.attributes("-topmost", False))
        win.focus_force()

    def _clear_tab(self, tab):
        for w in tab.winfo_children():
            w.destroy()

    def _embed_mic_mbc_figure(self, fig, tab, fig_key, scroll_parent=None):
        parent = scroll_parent if scroll_parent is not None else tab
        if fig is None:
            ctk.CTkLabel(parent, text="Brak danych do wyświetlenia.").pack(pady=20)
            return
        self.mic_mbc_figures[fig_key] = fig
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=5, pady=5)

    def _refresh_mic_mbc_distribution(self):
        self._clear_tab(self.tab_mic_dist)
        state = self._mic_mbc_state
        bact = state["bact"]
        sub = self.combo_mic_substance.get()
        key = (bact, sub)
        mic_bio = state["mic_bact"].get(key, {}).get("bio_results", [])
        mbc_bio = state["mbc_bact"].get(key, {}).get("bio_results", [])
        try:
            fig = mic_plotting.draw_mic_mbc_distribution(bact, sub, mic_bio, mbc_bio, config=self.plot_config)
        except Exception as e:
            self._show_plot_error(self.tab_mic_dist, str(e))
            return
        self._embed_mic_mbc_figure(fig, self.tab_mic_dist, "mic_distribution")

    def _render_mic_mbc_pairs(self):
        self._clear_tab(self.tab_mic_pairs)
        state = self._mic_mbc_state
        bact = state["bact"]
        pair_rows = []
        for sub in state["substances"]:
            key = (bact, sub)
            mic_entry = state["mic_bact"].get(key)
            mbc_entry = state["mbc_bact"].get(key)
            if not mic_entry:
                continue
            mic_summary = mic_entry["summary"]
            mbc_summary = mbc_entry["summary"] if mbc_entry else None
            error_reason, classification, ratio_display = None, None, None
            if mic_entry and mbc_entry:
                paired, _unmatched = mic_logic.pair_mic_mbc_by_bio_rep(mic_entry["bio_results"], mbc_entry["bio_results"])
                for p in paired:
                    if p["ratio"]["status"] == "blad_spojnosci":
                        error_reason = p["ratio"]["reason"]
                        break
                if paired and error_reason is None:
                    ratio_summary = mic_logic.summarize_mbc_mic_ratio(paired)
                    classification = ratio_summary["classification"]
                    ratio_display = ratio_summary["ratio_display"]
            pair_rows.append({
                "Substancja": sub, "mic": mic_summary["median"],
                "mbc": mbc_summary["median"] if mbc_summary else {"mic_value": None, "censored": None},
                "classification": classification, "ratio_display": ratio_display, "error_reason": error_reason,
            })
        try:
            fig = mic_plotting.draw_mic_mbc_pairs(bact, pair_rows, config=self.plot_config)
        except Exception as e:
            self._show_plot_error(self.tab_mic_pairs, str(e))
            return
        self._embed_mic_mbc_figure(fig, self.tab_mic_pairs, "mic_pairs")

    def _render_mic_mbc_comparison(self):
        self._clear_tab(self.tab_mic_compare)
        scroll = ctk.CTkScrollableFrame(self.tab_mic_compare)
        scroll.pack(fill="both", expand=True, padx=5, pady=5)
        state = self._mic_mbc_state
        bact = state["bact"]
        method = self._mic_mbc_method()

        any_plot = False
        for endpoint_name, grouped, fig_key in (
            ("MIC", state["mic_bact"], "mic_comparison"),
            ("MBC", state["mbc_bact"], "mbc_comparison"),
        ):
            groups = {sub: grouped[(bact, sub)]["bio_results"] for sub in state["substances"] if (bact, sub) in grouped}
            groups = {sub: bio for sub, bio in groups.items() if any(b.get("mic_value") is not None for b in bio)}
            if len(groups) < 2:
                continue
            try:
                comparison = mic_logic.compare_mic_groups(groups, method=method)
                fig = mic_plotting.draw_mic_group_comparison(endpoint_name, comparison, label=bact, config=self.plot_config)
            except Exception as e:
                ctk.CTkLabel(scroll, text=f"Błąd porównania {endpoint_name}: {e}").pack(pady=10)
                continue
            self._embed_mic_mbc_figure(fig, self.tab_mic_compare, fig_key, scroll_parent=scroll)
            any_plot = True

        if not any_plot:
            ctk.CTkLabel(
                scroll,
                text="Za mało substancji z danymi (potrzeba >=2) dla MIC lub MBC, żeby zbudować porównanie."
            ).pack(pady=20)

    def _render_mic_mbc_table(self):
        self._clear_tab(self.tab_mic_table)
        state = self._mic_mbc_state
        table_rows = mic_logic.build_mic_summary_rows(state["mic_bact"], state["mbc_bact"])
        replicate_rows = mic_logic.build_mic_replicate_rows(state["mic_bact"], state["mbc_bact"])
        state["table_rows"] = table_rows
        state["replicate_rows"] = replicate_rows

        box = ctk.CTkTextbox(self.tab_mic_table, font=("Consolas", 12))
        box.pack(fill="both", expand=True, padx=5, pady=5)

        header = f"{'Substancja':<15}{'n_bio MIC':<11}{'MIC':<14}{'n_bio MBC':<11}{'MBC':<14}{'Iloraz':<9}{'Klasyfikacja':<18}\n"
        box.insert("end", header)
        box.insert("end", "-" * len(header) + "\n")
        any_warning = False
        for row in table_rows:
            box.insert(
                "end",
                f"{row['Substancja']:<15}{row['n_bio_MIC']:<11}{row['MIC']:<14}{row['n_bio_MBC']:<11}"
                f"{row['MBC']:<14}{row['Iloraz_MBC_MIC']:<9}{row['Klasyfikacja']:<18}\n"
            )
            if row["Uwagi"]:
                any_warning = True
                box.insert("end", f"    UWAGA: {row['Uwagi']}\n")
        if not any_warning:
            box.insert("end", "\n(Brak ostrzeżeń dla tego szczepu.)\n")

        # Też do głównego logu aplikacji - widoczność "wszędzie", nie tylko w tym oknie.
        self.log(f"=== MIC/MBC: {state['bact']} ===")
        for row in table_rows:
            self.log(
                f"  {row['Substancja']}: MIC={row['MIC']} MBC={row['MBC']} "
                f"iloraz={row['Iloraz_MBC_MIC']} klasyfikacja={row['Klasyfikacja']}"
            )
            if row["Uwagi"]:
                self.log(f"    UWAGA: {row['Uwagi']}")

    def _export_mic_mbc_excel(self):
        state = getattr(self, "_mic_mbc_state", None)
        if not state or not state.get("table_rows"):
            messagebox.showwarning("Uwaga", "Najpierw otwórz analizę MIC/MBC.")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel file", "*.xlsx")])
        if not file_path:
            return
        warnings_list = [f"{r['Bakteria']} / {r['Substancja']}: {r['Uwagi']}" for r in state["table_rows"] if r["Uwagi"]]
        try:
            reports.export_mic_mbc_excel(file_path, state["table_rows"], state["replicate_rows"], warnings_list)
            messagebox.showinfo("Sukces", f"Zapisano wyniki MIC/MBC w:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Błąd Zapisu", str(e))

    def _export_mic_mbc_pdf(self):
        state = getattr(self, "_mic_mbc_state", None)
        if not state or not state.get("table_rows"):
            messagebox.showwarning("Uwaga", "Najpierw otwórz analizę MIC/MBC.")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF Document", "*.pdf")])
        if not file_path:
            return
        figures = [
            (title, self.mic_mbc_figures[key]) for key, title in (
                ("mic_distribution", f"Rozkład MIC/MBC: {self.combo_mic_substance.get()}"),
                ("mic_pairs", f"Odstęp MIC↔MBC: {state['bact']}"),
                ("mic_comparison", f"Porównanie substancji (MIC): {state['bact']}"),
                ("mbc_comparison", f"Porównanie substancji (MBC): {state['bact']}"),
            ) if key in self.mic_mbc_figures
        ]
        meta = {"date": datetime.now().strftime('%Y-%m-%d %H:%M'), "bact": state["bact"]}
        mic_mbc_data = {"bact": state["bact"], "table_rows": state["table_rows"], "figures": figures}
        try:
            success, msg = reports.generate_pdf(file_path, meta, mic_mbc_data=mic_mbc_data)
            if success:
                messagebox.showinfo("Sukces", msg)
            else:
                messagebox.showerror("Błąd PDF", msg)
        except Exception as e:
            messagebox.showerror("Błąd PDF", str(e))


