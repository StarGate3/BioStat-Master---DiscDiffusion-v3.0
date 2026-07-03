import pandas as pd
import numpy as np
from scipy import stats
import scikit_posthocs as sp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import utils
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from config import DISC_DIAMETER_MM, ALPHA, COL_GROUP, COL_MEASUREMENT, MIC_MIN_R2

class StatsEngine:
    def __init__(self):
        pass

    def run_statistics(self, df_run, method, ref_group):
        """
        Calculates main statistics (ANOVA/Kruskal) and Post-hoc.
        Returns:
            - results_summary (dict): 'test_name', 'p_value', 'statistic'
            - posthoc_df (DataFrame or None)
            - error_msg (str or None)
        """
        # Przygotowanie danych
        valid_groups = []
        dane_list = []
        
        # Filtrujemy grupy z < 2 pomiarami
        for g in df_run[COL_GROUP].unique():
             data = df_run[df_run[COL_GROUP] == g][COL_MEASUREMENT].values
             if len(data) >= 2:
                valid_groups.append(g)
                dane_list.append(data)

        if len(dane_list) < 2:
            return None, None, "Za mało ważnych grup do przeprowadzenia testów statystycznych."

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
            normality_results.append({"Grupa": g, "Shapiro p-value": p_shapiro, "Rozkład Normalny?": "TAK" if is_norm else "NIE"})

        # 2. Levene
        try: 
            stat, p_levene = stats.levene(*dane_list)
        except Exception as e: 
            print(f"Warning: Levene error: {e}")
            p_levene = 0
        
        use_parametric = all_normal and p_levene > ALPHA
        
        stats_main = []
        posthoc_df = None
        test_used = ""

        # 3. Testy Główne + Post-hoc
        if use_parametric:
            test_used = "ANOVA"
            try:
                f, p = stats.f_oneway(*dane_list)
                stats_main = [{"Test": "ANOVA", "Statistic": f, "p-value": p}]
                if p < ALPHA:
                    tukey = pairwise_tukeyhsd(df_run[COL_MEASUREMENT], df_run[COL_GROUP], ALPHA)
                    posthoc_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
            except Exception as e: return None, None, f"Błąd ANOVA: {e}"
        else:
            test_used = "Kruskal-Wallis"
            try:
                h, p = stats.kruskal(*dane_list)
                stats_main = [{"Test": "Kruskal-Wallis", "Statistic": h, "p-value": p}]
                if p < ALPHA:
                    posthoc_df = sp.posthoc_dunn(df_run, COL_MEASUREMENT, COL_GROUP, p_adjust=method)
            except Exception as e: return None, None, f"Błąd Kruskal: {e}"

        return {
            "normality": normality_results,
            "main_stats": stats_main,
            "test_used": test_used,
            "is_parametric": use_parametric,
            "all_normal": all_normal
        }, posthoc_df, None

    def process_detailed_results(self, posthoc_df, df_data, ref_group, test_type):
        """
        Przetwarza wyniki post-hoc na listę detali z Effect Size.
        Zwraca: (detailed_list, significant_set)
        """
        detailed_results = []
        sig_set = set()
        
        if posthoc_df is None: return [], set()

        # TUKEY
        if test_type == "ANOVA": 
            for i, r in posthoc_df.iterrows():
                g1, g2 = r['group1'], r['group2']
                is_sig = r['reject']
                p_adj = r['p-adj']
                
                self._add_detail(g1, g2, p_adj, is_sig, df_data, ref_group, detailed_results, sig_set)

        # DUNN (Kruskal)
        elif test_type == "Kruskal-Wallis":
            seen = set()
            for c in posthoc_df.columns:
                for r in posthoc_df.index:
                    if c != r:
                        pair = tuple(sorted((str(r), str(c))))
                        if pair not in seen:
                            pval = posthoc_df.loc[r, c]
                            is_sig = pval < ALPHA
                            self._add_detail(r, c, pval, is_sig, df_data, ref_group, detailed_results, sig_set)
                            seen.add(pair)
        
        return detailed_results, sig_set

    def _add_detail(self, g1, g2, p_val, is_sig, df_data, ref, results_list, sig_set):
        data1 = df_data[df_data[COL_GROUP] == g1][COL_MEASUREMENT].values
        data2 = df_data[df_data[COL_GROUP] == g2][COL_MEASUREMENT].values
        
        d_val = utils.calculate_cohens_d(data1, data2)
        d_interp = utils.get_effect_size_interpretation(d_val)
        
        results_list.append({
            "Group 1": g1, "Group 2": g2, "P-adj": p_val, 
            "Significant": is_sig, "Cohen's d": d_val, "Effect Size": d_interp
        })

        if is_sig:
            if g1 == ref: sig_set.add(g2)
            if g2 == ref: sig_set.add(g1)

    def run_pca(self, df, col_bact, selected_substances):
        """
        Runs PCA on the dataframe to visualize bacterial similarity based on sensitivity.
        Rows: Bacteria, Columns: Substances, Values: Mean Zone Diameter.
        """
        # 1. Filtrujemy dane tylko dla wybranych substancji
        df_filtered = df[df[COL_GROUP].isin(selected_substances)].copy()

        # 2. Pivot Table: Wiersze=Bakterie, Kolumny=Substancje
        df_pivot = df_filtered.pivot_table(index=col_bact, columns=COL_GROUP, values=COL_MEASUREMENT, aggfunc='mean')
        
        if df_pivot.empty or len(df_pivot) < 3:
            return None, "Za mało danych do PCA (wymagane min. 3 szczepy)."
            
        # 3. Uzupełnianie braków (jeśli jakaś bakteria nie ma pomiaru dla danej substancji -> 0)
        df_pivot = df_pivot.fillna(0)
        
        # Wymagane min 2 kolumny (cechy)
        if df_pivot.shape[1] < 2:
            return None, "Za mało cech do PCA (wymagane min. 2 substancje)."

        # 4. Skalowanie
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df_pivot)
        
        # 5. PCA
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(scaled_data)
        
        pca_df = pd.DataFrame(data=pcs, columns=['PC1', 'PC2'])
        pca_df['Bakteria'] = df_pivot.index.values # use values to avoid index issues
        
        explained_variance = pca.explained_variance_ratio_
        return (pca_df, explained_variance), None

    def estimate_mic(self, df, selected_substances, target_diameter=DISC_DIAMETER_MM):
        """
        Estimates MIC for each substance using Log-Linear Regression.
        Model: Diameter = a + b * ln(Concentration)
        MIC = exp((Target - a) / b)

        Every substance in `selected_substances` gets an entry in the
        returned dict, even when no reliable MIC could be computed -
        "Status"/"Reason" always explain what happened (never a silent
        omission). A numeric "MIC" is only populated when the log-linear
        fit meets MIC_MIN_R2; "Extrapolated" flags an estimate that falls
        outside the actually-tested concentration range.
        """
        results = {}

        for sub in selected_substances:
            # 1. Pobierz dane tylko dla tej substancji
            # Musimy wyciągnąć stężenia z nazw grup.
            # Używamy utils.parse_concentration

            sub_df = df[df[COL_GROUP].str.contains(sub, regex=False)].copy() # Wstępne filtrowanie, ale dokładne parsowanie niżej

            x_concs = []
            y_diams = []
            valid_unit = ""

            for g in sub_df[COL_GROUP].unique():
                parsed_sub, conc, unit = utils.parse_concentration(g)
                # Sprawdź czy to ta substancja (bo contains jest luźne)
                if parsed_sub and sub in parsed_sub and conc is not None:
                    # Pobierz wszystkie pomiary dla tej grupy
                    measurements = sub_df[sub_df[COL_GROUP] == g][COL_MEASUREMENT].values
                    for m in measurements:
                        if conc > 0:
                            x_concs.append(conc)
                            y_diams.append(m)
                            valid_unit = unit

            n_points = len(set(x_concs))
            if n_points < 3:
                results[sub] = {
                    "MIC": None, "Unit": valid_unit, "R2": None, "Slope": None, "Intercept": None,
                    "Status": "za_malo_stezen", "Extrapolated": False, "ConcRange": None,
                    "Reason": f"Za mało unikalnych stężeń do regresji (znaleziono {n_points}, wymagane min. 3).",
                }
                continue

            # 2. Regresja Liniowa na logarytmach
            try:
                log_x = np.log(x_concs)
                slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, y_diams)
                r2 = r_value ** 2
                conc_range = (min(x_concs), max(x_concs))
            # Skip this substance on degenerate regression input: linregress
            # raises ValueError on non-finite / zero-variance data, TypeError
            # guards non-numeric slipping past parse_concentration, LinAlgError
            # for any underlying solver failure. KeyboardInterrupt / SystemExit
            # still propagate correctly.
            except (ValueError, TypeError, np.linalg.LinAlgError) as e:
                results[sub] = {
                    "MIC": None, "Unit": valid_unit, "R2": None, "Slope": None, "Intercept": None,
                    "Status": "blad_regresji", "Extrapolated": False, "ConcRange": None,
                    "Reason": f"Regresja nie powiodła się ({type(e).__name__}: {e}).",
                }
                continue

            if slope <= 0:
                # Oczekujemy że strefa rośnie ze stężeniem - ujemny/zerowy
                # współczynnik kierunkowy nie ma sensu biologicznego.
                results[sub] = {
                    "MIC": None, "Unit": valid_unit, "R2": r2, "Slope": slope, "Intercept": intercept,
                    "Status": "ujemny_slope", "Extrapolated": False, "ConcRange": conc_range,
                    "Reason": f"Strefa nie rośnie ze stężeniem (nachylenie={slope:.3f} <= 0) - brak sensu biologicznego.",
                }
                continue

            # 3. Oblicz MIC: target = a + b * ln(MIC) => ln(MIC) = (target - a) / b
            ln_mic = (target_diameter - intercept) / slope
            mic = np.exp(ln_mic)
            extrapolated = not (conc_range[0] <= mic <= conc_range[1])

            if r2 < MIC_MIN_R2:
                results[sub] = {
                    "MIC": None, "Unit": valid_unit, "R2": r2, "Slope": slope, "Intercept": intercept,
                    "Status": "slabe_dopasowanie", "Extrapolated": extrapolated, "ConcRange": conc_range,
                    "Reason": (
                        f"MIC odrzucone: słabe dopasowanie regresji (R²={r2:.2f} < próg {MIC_MIN_R2:g}). "
                        f"Wartość {mic:.4g} {valid_unit} nie jest wiarygodna i nie jest raportowana."
                    ),
                }
                continue

            if extrapolated:
                reason = (
                    f"UWAGA: MIC ({mic:.4g} {valid_unit}) leży POZA przetestowanym zakresem stężeń "
                    f"({conc_range[0]:g}-{conc_range[1]:g} {valid_unit}) - to ekstrapolacja poza zbadane dane, "
                    f"traktuj wyłącznie orientacyjnie."
                )
            else:
                reason = "Dopasowanie w normie, MIC mieści się w przetestowanym zakresie stężeń."

            results[sub] = {
                "MIC": mic,
                "Unit": valid_unit,
                "R2": r2,
                "Slope": slope,
                "Intercept": intercept,
                "Status": "ekstrapolacja" if extrapolated else "ok",
                "Extrapolated": extrapolated,
                "ConcRange": conc_range,
                "Reason": reason,
            }

        return results
