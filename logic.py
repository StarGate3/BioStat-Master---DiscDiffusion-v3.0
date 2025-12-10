import pandas as pd
import numpy as np
from scipy import stats
import scikit_posthocs as sp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import utils
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

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
        for g in df_run['Grupa'].unique():
             data = df_run[df_run['Grupa'] == g]['Srednica_mm'].values
             if len(data) >= 2:
                valid_groups.append(g)
                dane_list.append(data)

        if len(dane_list) < 2:
            return None, None, "Za mało ważnych grup do przeprowadzenia testów statystycznych."

        # 1. Normalność
        all_normal = True
        normality_results = []
        for g in valid_groups:
            vals = df_run[df_run['Grupa'] == g]['Srednica_mm']
            p_shapiro = 0
            is_norm = False
            if len(vals) >= 3 and vals.std() > 0:
                s, p_shapiro = stats.shapiro(vals)
                if p_shapiro >= 0.05: is_norm = True
            if not is_norm: all_normal = False
            normality_results.append({"Grupa": g, "Shapiro p-value": p_shapiro, "Rozkład Normalny?": "TAK" if is_norm else "NIE"})

        # 2. Levene
        try: 
            stat, p_levene = stats.levene(*dane_list)
        except Exception as e: 
            print(f"Warning: Levene error: {e}")
            p_levene = 0
        
        use_parametric = all_normal and p_levene > 0.05
        
        stats_main = []
        posthoc_df = None
        test_used = ""

        # 3. Testy Główne + Post-hoc
        if use_parametric:
            test_used = "ANOVA"
            try:
                f, p = stats.f_oneway(*dane_list)
                stats_main = [{"Test": "ANOVA", "Statistic": f, "p-value": p}]
                if p < 0.05:
                    tukey = pairwise_tukeyhsd(df_run['Srednica_mm'], df_run['Grupa'], 0.05)
                    posthoc_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
            except Exception as e: return None, None, f"Błąd ANOVA: {e}"
        else:
            test_used = "Kruskal-Wallis"
            try:
                h, p = stats.kruskal(*dane_list)
                stats_main = [{"Test": "Kruskal-Wallis", "Statistic": h, "p-value": p}]
                if p < 0.05:
                    posthoc_df = sp.posthoc_dunn(df_run, 'Srednica_mm', 'Grupa', p_adjust=method)
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
                            is_sig = pval < 0.05
                            self._add_detail(r, c, pval, is_sig, df_data, ref_group, detailed_results, sig_set)
                            seen.add(pair)
        
        return detailed_results, sig_set

    def _add_detail(self, g1, g2, p_val, is_sig, df_data, ref, results_list, sig_set):
        data1 = df_data[df_data['Grupa'] == g1]['Srednica_mm'].values
        data2 = df_data[df_data['Grupa'] == g2]['Srednica_mm'].values
        
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
        df_filtered = df[df['Grupa'].isin(selected_substances)].copy()
        
        # 2. Pivot Table: Wiersze=Bakterie, Kolumny=Substancje
        df_pivot = df_filtered.pivot_table(index=col_bact, columns='Grupa', values='Srednica_mm', aggfunc='mean')
        
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

    def estimate_mic(self, df, selected_substances, target_diameter=6.0):
        """
        Estimates MIC for each substance using Log-Linear Regression.
        Model: Diameter = a + b * ln(Concentration)
        MIC = exp((Target - a) / b)
        """
        results = {}
        
        for sub in selected_substances:
            # 1. Pobierz dane tylko dla tej substancji
            # Musimy wyciągnąć stężenia z nazw grup.
            # Używamy utils.parse_concentration
            
            sub_df = df[df['Grupa'].str.contains(sub, regex=False)].copy() # Wstępne filtrowanie, ale dokładne parsowanie niżej
            
            x_concs = []
            y_diams = []
            valid_unit = ""
            
            for g in sub_df['Grupa'].unique():
                parsed_sub, conc, unit = utils.parse_concentration(g)
                # Sprawdź czy to ta substancja (bo contains jest luźne)
                if parsed_sub and sub in parsed_sub and conc is not None:
                    # Pobierz wszystkie pomiary dla tej grupy
                    measurements = sub_df[sub_df['Grupa'] == g]['Srednica_mm'].values
                    for m in measurements:
                        if conc > 0:
                            x_concs.append(conc)
                            y_diams.append(m)
                            valid_unit = unit

            if len(set(x_concs)) < 3:
                # Za mało punktów stężeń do regresji (min 3 lepie)
                continue
                
            # 2. Regresja Liniowa na logarytmach
            try:
                log_x = np.log(x_concs)
                slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, y_diams)
                
                # 3. Oblicz MIC: 6 = a + b * ln(MIC) => ln(MIC) = (6 - a) / b
                if slope > 0: # Oczekujemy że strefa rośnie ze stężeniem
                    ln_mic = (target_diameter - intercept) / slope
                    mic = np.exp(ln_mic)
                else:
                    mic = None # Ujemny lub zerowy współczynnik kierunkowy - brak sensu biol.
                
                results[sub] = {
                    "MIC": mic, 
                    "Unit": valid_unit, 
                    "R2": r_value**2,
                    "Slope": slope,
                    "Intercept": intercept
                }
            except:
                pass
                
        return results
