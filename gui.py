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
        
        # --- FIGURY ---
        self.figures = {
            'bar': None, 'heat': None, 'pvalue': None,
            'trend': None, 'effect': None, 'cross': None, 'pca': None
        }
        
        # --- KONFIGURACJA ---
        self.plot_config = {
            "font_labels": 10, "font_title": 12, "axis_max": 0, "star_offset": 0.03, "bar_width": 0.8,
            "show_disk_line": True, "palette": "viridis", "transparent_background": True,
            "plot_type": "Barplot (Słupkowy)", "error_bar": "SD (Odchylenie Std.)", "show_points": False
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
        self.right_frame.grid_rowconfigure(1, weight=1)
        ctk.CTkLabel(self.right_frame, text="Wybór próbek:", font=ctk.CTkFont(weight="bold")).grid(row=0, column=0, padx=10, pady=(10, 5))
        self.scroll_samples = ctk.CTkScrollableFrame(self.right_frame, label_text="Dostępne grupy", height=600)
        self.scroll_samples.grid(row=1, column=0, padx=10, pady=5, sticky="nsew")
        self.btn_select_all = ctk.CTkButton(self.right_frame, text="Zaznacz wszystko", width=100, command=self.select_all)
        self.btn_select_all.grid(row=2, column=0, padx=10, pady=5)
        self.btn_deselect_all = ctk.CTkButton(self.right_frame, text="Odznacz wszystko", width=100, fg_color="gray", command=self.deselect_all)
        self.btn_deselect_all.grid(row=3, column=0, padx=10, pady=(5, 20))

    # ==================== LOGIKA POMOCNICZA ====================
    def log(self, text):
        self.textbox.insert("end", text + "\n")
        self.textbox.see("end")
    def clear_log(self): self.textbox.delete("1.0", "end")

    def load_file(self):
        path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx *.xls")])
        if path:
            try:
                self.df = pd.read_excel(path)
                self.df.columns = self.df.columns.str.strip()
                for col in self.df.select_dtypes(['object']).columns:
                    self.df[col] = self.df[col].str.strip()
                
                self.lbl_file.configure(text=path.split("/")[-1], text_color="white")
                col = next((c for c in self.df.columns if 'Bakteri' in c), None)
                if col:
                    self.col_bact_name = col
                    bacts = list(self.df[col].unique())
                    self.combo_bact.configure(values=bacts)
                    self.combo_bact.set(bacts[0])
                    self.on_bacteria_change(bacts[0])
                    self.log(f"Wczytano plik. Znaleziono szczepy: {bacts}")
                else: messagebox.showerror("Błąd", "Brak kolumny 'Bakterie'.")
            except Exception as e: messagebox.showerror("Błąd", f"Nie udało się wczytać: {e}")

    def on_bacteria_change(self, selected_bact):
        if self.df is None: return
        try:
            all_groups = sorted(self.df['Grupa'].unique(), key=utils.smart_sort_key)
            
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
            grupy_bact = sorted(df_temp['Grupa'].unique(), key=utils.smart_sort_key)
            self.combo_ref.configure(values=grupy_bact)
            woda = next((g for g in grupy_bact if "woda" in g.lower() or "kontrol" in g.lower()), None)
            if woda: self.combo_ref.set(woda)
            elif grupy_bact: self.combo_ref.set(grupy_bact[0])
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
        post_hoc = self.combo_method.get()
        if post_hoc == "None": post_hoc = "no correction"
        elif post_hoc == "fdr_bh": post_hoc = "Benjamini-Hochberg (FDR) correction"
        elif post_hoc == "holm": post_hoc = "Holm-Bonferroni correction"
        else: post_hoc = "Bonferroni correction"

        test_name = "Statistical test" 
        if self.export_stats_main:
            used_test = self.export_stats_main[0].get("Test", "")
            if "ANOVA" in used_test: test_name = "One-way ANOVA"
            elif "Kruskal" in used_test: test_name = "Kruskal-Wallis test"

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
Statistical significance was determined using {test_name} followed by {post_hoc} for multiple comparisons.
Asterisks (*) indicate a statistically significant difference (p < 0.05) compared to the negative control ({ref_group}).
Red dashed line represents the diameter of the disk (6 mm).

=== Rycina 2: Mapa Ciepła ===
Figure 2. Heatmap visualizing the magnitude of growth inhibition zones (mm) for {bact} treated with various substances.
Color intensity corresponds to the mean diameter of the inhibition zone. Warmer colors indicate higher antibacterial activity.

=== Rycina 3: Mapa Wielkości Efektu (Effect Size) ===
Figure 3. Lollipop chart displaying the standardized effect size (Cohen's d) for statistically significant pairwise comparisons.
Dots represent the magnitude of the difference between groups. Green dots indicate a positive difference (Group 1 > Group 2), while red dots indicate a negative difference.

=== Rycina 4: Mapa Istotności (P-value Matrix) ===
Figure 4. Pairwise comparison significance matrix (P-values).
The heatmap displays adjusted p-values for all pairwise comparisons. Blue shades indicate statistical significance (p < 0.05), while red/white shades indicate non-significant differences.
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

        self.switch_line = ctk.CTkSwitch(self.settings_win, text="Pokaż linię krążka (6mm)")
        if self.plot_config["show_disk_line"]: self.switch_line.select()
        else: self.switch_line.deselect()
        self.switch_line.pack(pady=10)

        self.switch_trans = ctk.CTkSwitch(self.settings_win, text="Zapisz z przezroczystym tłem")
        if self.plot_config["transparent_background"]: self.switch_trans.select()
        else: self.switch_trans.deselect()
        self.switch_trans.pack(pady=10)

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
    def run_analysis(self):
        if self.df is None: return
        bact = self.combo_bact.get()
        method = self.combo_method.get()
        ref_group = self.combo_ref.get()
        if method == "None": method = None

        wybrane = self.get_selected_groups()
        if not wybrane:
            messagebox.showwarning("Stop", "Nie wybrano próbek!")
            return

        # 1. Filtrowanie wstępne
        df_run = self.df[
            (self.df[self.col_bact_name] == bact) & 
            (self.df['Grupa'].isin(wybrane))
        ].copy()

        if df_run.empty: return

        # 2. Outliery (UI Logic)
        outliers_data = utils.find_outliers_dixon(df_run)
        if outliers_data:
            dialog = OutlierDialog(self, outliers_data)
            self.wait_window(dialog) 
            if dialog.result:
                for item in dialog.result:
                    mask = (df_run['Grupa'] == item['Group']) & (df_run['Srednica_mm'] == item['Srednica_mm'])
                    idx = df_run[mask].first_valid_index()
                    if idx is not None: df_run = df_run.drop(idx)
                self.log(f"!!! USUNIĘTO {len(dialog.result)} WARTOŚCI ODSTAJĄCYCH !!!")

        self.export_data_raw = df_run

        # 3. STAT ENGINE (Delegacja)
        summary_res, posthoc_df, error = self.stats_engine.run_statistics(df_run, method, ref_group)
        
        if error:
            self.log(f"Blad Statystyki: {error}")
            return
        
        # Logowanie wyników
        self.clear_log()
        self.log(f"=== RAPORT v3: {bact} ===")
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

        # 4. POST HOC DETALE (Delegacja)
        detailed, sig_set = self.stats_engine.process_detailed_results(posthoc_df, df_run, ref_group, summary_res['test_used'])
        self.posthoc_detailed_results = detailed
        
        if detailed:
            self.log("\n[3] WYNIKI SZCZEGÓŁOWE (Effect Size):")
            for d in detailed:
                if d['Significant']:
                    metrics_d = d["Cohen's d"]
                    self.log(f"{d['Group 1']} vs {d['Group 2']} | p={d['P-adj']:.4f} | d={metrics_d:.2f}")

        # 5. RYSOWANIE (Delegacja)
        # 5. RYSOWANIE (Delegacja)
        self.display_plot(lambda: self.plotter.draw_bar_plot(df_run, bact, ref_group, sig_set), self.tab_plot, 'bar')
        self.display_plot(lambda: self.plotter.draw_heatmap(df_run, bact), self.tab_heatmap, 'heat')
        self.display_plot(lambda: self.plotter.draw_pvalue_heatmap(self.export_stats_posthoc, bact), self.tab_pvalue, 'pvalue')
        
        # MIC ESTIMATION
        mic_results = self.stats_engine.estimate_mic(df_run, wybrane) # wybrane -> trzeba przekazac nazwy substancji, a nie grupy?
        # estimate_mic oczekuje selected_substances, ale w logice wyciągamy substancje z nazw grup. 
        # Przekażemy unikalne substancje z nazw grup.
        
        # W logic.py estimate_mic(df, selected_substances)
        # selected_substances to lista stringow (np. ['Amikacyna', 'Gentamycyna'])
        # A 'wybrane' to lista ['Amikacyna 10ug', 'Gentamycyna 5ug']
        # Musimy wyciągnąć unikalne nazwy
        
        unique_subs = set()
        for g in wybrane:
            s, _, _ = utils.parse_concentration(g)
            if s: unique_subs.add(s)
            
        mic_results = self.stats_engine.estimate_mic(df_run, list(unique_subs))
        
        if mic_results:
            self.log("\n[4] Oszacowane MIC (Theoretical):")
            for sub, res in mic_results.items():
                if res['MIC']:
                    self.log(f"{sub}: {res['MIC']:.3f} {res['Unit']} (R2={res['R2']:.2f})")
                else:
                    self.log(f"{sub}: Nie można wyznaczyć (<0 slope)")

        # Pass mic_results to draw_trend
        fig_trend, err = self.plotter.draw_trend(df_run, bact, mic_data=mic_results)
        if fig_trend: 
             self.display_figure(fig_trend, self.tab_trend, 'trend')
        elif err:
             self._show_plot_error(self.tab_trend, err)

        self.display_plot(lambda: self.plotter.draw_cross_species(self.df, self.col_bact_name, wybrane), self.tab_cross, 'cross')
        self.display_plot(lambda: self.plotter.draw_effect_plot(self.posthoc_detailed_results), self.tab_effect, 'effect')

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
            fig_to_save.savefig(file_path, dpi=300, bbox_inches='tight', transparent=is_transparent)
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
                if self.export_stats_normality: pd.DataFrame(self.export_stats_normality).to_excel(writer, sheet_name="Normalnosc", index=False)
                if self.export_stats_main: pd.DataFrame(self.export_stats_main).to_excel(writer, sheet_name="Test Glowny", index=False)
                if self.posthoc_detailed_results: pd.DataFrame(self.posthoc_detailed_results).to_excel(writer, sheet_name="Post-hoc (Details)", index=False)
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
            'ref': self.combo_ref.get()
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

