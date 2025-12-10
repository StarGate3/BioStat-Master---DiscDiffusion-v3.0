import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
import utils

class Plotter:
    def __init__(self, config):
        """
        Inicjalizacja z konfiguracją wykresów.
        config: dict z ustawieniami (font_labels, palette, etc.)
        """
        self.config = config
    
    def update_config(self, new_config):
        self.config = new_config

    def draw_bar_plot(self, df, bact, ref, sig_set):
        plt.close('all') 
        is_horiz = False 
        
        is_horiz = (self.config.get("orientation", "Pozioma") == "Pozioma")

        order = sorted(df['Grupa'].unique(), key=utils.smart_sort_key)
        if ref in order:
            order.remove(ref)
            order.insert(0, ref)

        f_lbl = self.config["font_labels"]
        f_ttl = self.config["font_title"]
        ax_max = self.config["axis_max"]
        s_off = self.config["star_offset"]
        pal = self.config["palette"] 
        show_line = self.config["show_disk_line"] 
        plot_type = self.config["plot_type"]
        error_bar_choice = self.config["error_bar"]
        show_points = self.config["show_points"]

        if "SD" in error_bar_choice: sb_error = "sd"
        elif "SEM" in error_bar_choice: sb_error = "se"
        elif "95% CI" in error_bar_choice: sb_error = ("ci", 95)
        else: sb_error = "sd"

        h = max(6, len(order)*0.4) if is_horiz else 6
        w = 8 if is_horiz else max(8, len(order)*0.3)
        
        fig = plt.Figure(figsize=(w, h), dpi=100)
        ax = fig.add_subplot(111)
        
        means = df.groupby('Grupa')['Srednica_mm'].mean()
        sds = df.groupby('Grupa')['Srednica_mm'].std().fillna(0)
        sems = df.groupby('Grupa')['Srednica_mm'].sem().fillna(0)
        maxs = df.groupby('Grupa')['Srednica_mm'].max()
        
        if "Barplot" in plot_type:
            if is_horiz:
                sns.barplot(x='Srednica_mm', y='Grupa', data=df, order=order, ax=ax, capsize=0.2, errorbar=sb_error, palette=pal, orient='h', edgecolor='black', hue='Grupa', legend=False)
                if show_points: sns.stripplot(x='Srednica_mm', y='Grupa', data=df, order=order, ax=ax, color='black', alpha=0.6, jitter=True, size=4)
            else:
                sns.barplot(x='Grupa', y='Srednica_mm', data=df, order=order, ax=ax, capsize=0.2, errorbar=sb_error, palette=pal, orient='v', edgecolor='black', hue='Grupa', legend=False)
                if show_points: sns.stripplot(x='Grupa', y='Srednica_mm', data=df, order=order, ax=ax, color='black', alpha=0.6, jitter=True, size=4)
            
            if "SD" in error_bar_choice: ref_points = means + sds
            elif "SEM" in error_bar_choice: ref_points = means + sems
            elif "95% CI" in error_bar_choice: ref_points = means + (1.96 * sems) 
            else: ref_points = means + sds
            max_val_data = ref_points.max()

        elif "Boxplot" in plot_type:
            if is_horiz:
                sns.boxplot(x='Srednica_mm', y='Grupa', data=df, order=order, ax=ax, palette=pal, orient='h', hue='Grupa', legend=False)
                if show_points: sns.stripplot(x='Srednica_mm', y='Grupa', data=df, order=order, ax=ax, color='black', alpha=0.6, jitter=True, size=4)
            else:
                sns.boxplot(x='Grupa', y='Srednica_mm', data=df, order=order, ax=ax, palette=pal, orient='v', hue='Grupa', legend=False)
                if show_points: sns.stripplot(x='Grupa', y='Srednica_mm', data=df, order=order, ax=ax, color='black', alpha=0.6, jitter=True, size=4)
            max_val_data = maxs.max()
            ref_points = maxs 

        elif "Violinplot" in plot_type:
            if is_horiz:
                sns.violinplot(x='Srednica_mm', y='Grupa', data=df, order=order, ax=ax, palette=pal, orient='h', inner="stick", hue='Grupa', legend=False)
                if show_points: sns.stripplot(x='Srednica_mm', y='Grupa', data=df, order=order, ax=ax, color='black', alpha=0.6, jitter=True, size=4)
            else:
                sns.violinplot(x='Grupa', y='Srednica_mm', data=df, order=order, ax=ax, palette=pal, orient='v', inner="stick", hue='Grupa', legend=False)
                if show_points: sns.stripplot(x='Grupa', y='Srednica_mm', data=df, order=order, ax=ax, color='black', alpha=0.6, jitter=True, size=4)
            max_val_data = maxs.max()
            ref_points = maxs

        if ax_max > 0: final_limit = ax_max
        else: final_limit = max_val_data * 1.15

        if is_horiz:
            if show_line: ax.axvline(x=6, color='red', linestyle='--', alpha=0.5, label='Krążek (6mm)')
            ax.set_xlim(0, final_limit)
            ax.tick_params(axis='y', labelsize=f_lbl) 
            ax.tick_params(axis='x', labelsize=f_lbl)
            ax.set_xlabel("Średnica strefy (mm)", fontsize=f_ttl)
            ax.set_ylabel("", fontsize=f_ttl)
        else:
            if show_line: ax.axhline(y=6, color='red', linestyle='--', alpha=0.5, label='Krążek (6mm)')
            ax.set_ylim(0, final_limit)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=f_lbl) 
            ax.tick_params(axis='y', labelsize=f_lbl)
            ax.set_ylabel("Średnica strefy (mm)", fontsize=f_ttl)
            ax.set_xlabel("", fontsize=f_ttl)

        if show_line: ax.legend(loc='upper right')

        offset_val = df['Srednica_mm'].max() * s_off
        for i, g in enumerate(order):
            if g in sig_set:
                try:
                    base_pos = ref_points[g]
                    pos = base_pos + offset_val
                    if is_horiz: ax.text(pos, i, "*", va='center', fontweight='bold', fontsize=f_ttl+2)
                    else: ax.text(i, pos, "*", ha='center', fontweight='bold', fontsize=f_ttl+2)
                except: pass

        ax.set_title(f"{bact} vs {ref}", fontsize=f_ttl+2, fontweight='bold')
        fig.tight_layout()
        return fig

    def draw_heatmap(self, df, bact):
        df_mean = df.groupby('Grupa')['Srednica_mm'].mean().sort_values(ascending=False)
        data = df_mean.to_frame(name="Średnica (mm)")
        h = max(6, len(data) * 0.4) 
        fig = plt.Figure(figsize=(8, h), dpi=100) 
        ax = fig.add_subplot(111)
        
        pal = self.config["palette"]
        try: 
            sns.heatmap(data, annot=True, fmt=".1f", cmap=pal, ax=ax, linewidths=1, linecolor='white')
        except ValueError as e: 
            print(f"Warning: Palette '{pal}' error: {e}. Using magma.")
            sns.heatmap(data, annot=True, fmt=".1f", cmap="magma", ax=ax, linewidths=1, linecolor='white')
        
        ax.set_title(f"Mapa aktywności: {bact}", fontsize=self.config["font_title"]+2)
        ax.tick_params(axis='both', labelsize=self.config["font_labels"])
        ax.set_ylabel("")
        fig.tight_layout()
        return fig

    def draw_pvalue_heatmap(self, export_stats_posthoc, bact):
        if export_stats_posthoc is None:
            return None

        p_matrix = None
        if 'group1' in export_stats_posthoc.columns: # Tukey
            df_res = export_stats_posthoc
            groups = sorted(list(set(df_res['group1']) | set(df_res['group2'])))
            p_matrix = pd.DataFrame(index=groups, columns=groups, dtype=float)
            np.fill_diagonal(p_matrix.values, 1.0) 
            for _, row in df_res.iterrows():
                g1, g2 = row['group1'], row['group2']
                pval = row['p-adj']
                p_matrix.at[g1, g2] = pval
                p_matrix.at[g2, g1] = pval
        else: # Dunn
            p_matrix = export_stats_posthoc.copy()
            np.fill_diagonal(p_matrix.values, 1.0)

        if p_matrix is None: return None

        h = max(6, len(p_matrix) * 0.5)
        w = max(8, len(p_matrix) * 0.5)
        fig = plt.Figure(figsize=(w, h), dpi=100)
        ax = fig.add_subplot(111)
        
        mask = np.triu(np.ones_like(p_matrix, dtype=bool))
        sns.heatmap(p_matrix, mask=mask, annot=True, fmt=".3f", 
                    cmap="RdBu_r", center=0.05, vmin=0, vmax=1,
                    ax=ax, linewidths=1, linecolor='white',
                    cbar_kws={'label': 'P-value (Istotność)'})

        ax.set_title(f"Mapa Istotności (P-value): {bact}", fontsize=self.config["font_title"]+2)
        ax.tick_params(axis='both', labelsize=self.config["font_labels"])
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
        fig.tight_layout()
        return fig

    def draw_trend(self, df, bact, mic_data=None):
        f_lbl = self.config["font_labels"]
        f_ttl = self.config["font_title"]
        pal = self.config["palette"]
        ax_max = self.config["axis_max"]
        
        trend_data = []
        for g in df['Grupa'].unique():
            sub, conc, unit = utils.parse_concentration(g)
            if sub is not None:
                measurements = df[df['Grupa'] == g]['Srednica_mm'].values
                for m in measurements:
                    trend_data.append({"Substancja": sub, "Stężenie": conc, "Jednostka": unit, "Średnica": m})
        if not trend_data:
            return None, "Nie wykryto stężeń w nazwach grup."
        
        df_trend = pd.DataFrame(trend_data)
        fig = plt.Figure(figsize=(8, 6), dpi=100)
        ax = fig.add_subplot(111)
        
        sns.lineplot(data=df_trend, x="Stężenie", y="Średnica", hue="Substancja", style="Substancja", markers=True, dashes=False, ax=ax, err_style="band", palette=pal)
        ax.set_title(f"Zależność Dawka-Odpowiedź: {bact}", fontsize=f_ttl+2)
        ax.set_ylabel("Średnica strefy (mm)", fontsize=f_ttl)
        unit_label = df_trend['Jednostka'].iloc[0] if not df_trend.empty else ""
        ax.set_xlabel(f"Stężenie ({unit_label})", fontsize=f_ttl)
        ax.tick_params(axis='both', labelsize=f_lbl)
        if ax_max > 0: ax.set_ylim(0, ax_max)
        
        correlations = []
        for sub in df_trend['Substancja'].unique():
            sub_df = df_trend[df_trend['Substancja'] == sub]
            if len(sub_df) > 2:
                means = sub_df.groupby("Stężenie")["Średnica"].mean().reset_index()
                if np.std(means["Stężenie"]) > 0 and np.std(means["Średnica"]) > 0:
                    try:
                        r, p = stats.spearmanr(means["Stężenie"], means["Średnica"])
                        correlations.append(f"{sub}: r={r:.2f} (p={p:.3f})")
                    except Exception as e: 
                        correlations.append(f"{sub}: błąd ({str(e)})")
                else: correlations.append(f"{sub}: brak zmienności")
        
        if correlations:
            # Dodaj MIC do boxa
            if mic_data:
                correlations.append("\n[MIC Estimates (D=6mm)]")
                for sub, res in mic_data.items():
                    if res and res['MIC']:
                        correlations.append(f"{sub}: {res['MIC']:.2f} {res['Unit']}")
                    else:
                        correlations.append(f"{sub}: > max conc?")

            ax.text(0.02, 0.95, "\n".join(correlations), transform=ax.transAxes, fontsize=f_lbl, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
        fig.tight_layout()
        return fig, None

    def draw_effect_plot(self, posthoc_detailed_results):
        if not posthoc_detailed_results: return None
        
        sig_results = [r for r in posthoc_detailed_results if r['Significant']]
        
        if not sig_results: return None

        sig_results.sort(key=lambda x: abs(x["Cohen's d"]), reverse=False)

        labels = [f"{r['Group 1']}\nvs {r['Group 2']}" for r in sig_results]
        values = [r["Cohen's d"] for r in sig_results]
        colors_list = ['red' if v < 0 else 'green' for v in values]

        h = max(6, len(sig_results) * 0.45)
        fig = plt.Figure(figsize=(9, h), dpi=100)
        ax = fig.add_subplot(111)
        
        y_pos = np.arange(len(labels))
        ax.hlines(y=y_pos, xmin=0, xmax=values, color=colors_list, alpha=0.6, linewidth=2)
        ax.scatter(values, y_pos, color=colors_list, s=80, zorder=3)
        
        ax.axvline(0, color='black', linestyle='--', alpha=0.7)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=self.config["font_labels"])
        ax.set_xlabel("Wielkość Efektu (Cohen's d)", fontsize=self.config["font_title"])
        ax.set_title("Siła różnic między grupami (Istotne statystycznie)", fontsize=self.config["font_title"]+2)
        
        for i, v in enumerate(values):
            offset = max(1, abs(v)*0.05) if v >= 0 else -max(1, abs(v)*0.05)
            ha_align = 'left' if v >= 0 else 'right'
            ax.text(v + offset, i, f"{v:.1f}", va='center', ha=ha_align, fontsize=9, fontweight='bold')

        fig.tight_layout()
        return fig

    def draw_cross_species(self, df, col_bact_name, selected_substances):
        # Walidacja
        if not selected_substances: return None

        df_cross = df[df['Grupa'].isin(selected_substances)].copy()
        if df_cross.empty: return None

        f_lbl = self.config["font_labels"]
        f_ttl = self.config["font_title"]
        pal = self.config["palette"]
        ax_max = self.config["axis_max"]
        
        num_bact = len(df_cross[col_bact_name].unique())
        num_substances = len(selected_substances)
        
        calc_width = max(10, num_bact * 2.5)
        calc_height = 11 
        
        fig = plt.Figure(figsize=(calc_width, calc_height), dpi=100) 
        ax = fig.add_subplot(111)
        
        sns.barplot(
            data=df_cross, 
            x=col_bact_name, 
            y='Srednica_mm', 
            hue='Grupa', 
            ax=ax, 
            palette=pal, 
            capsize=0.04, 
            errorbar='sd', 
            edgecolor='black', 
            linewidth=0.8,
            width=0.85
        )
        
        ax.set_title("Porównanie Międzygatunkowe", fontsize=f_ttl+6, fontweight='bold', pad=25)
        ax.set_xlabel("Szczep bakterii", fontsize=f_ttl+2, labelpad=15)
        ax.set_ylabel("Średnica strefy (mm)", fontsize=f_ttl+2)
        ax.tick_params(axis='x', labelsize=f_lbl+2)
        ax.tick_params(axis='y', labelsize=f_lbl)
        ax.yaxis.grid(True, linestyle='--', alpha=0.5, color='gray', zorder=0)
        ax.set_axisbelow(True)
        sns.despine(ax=ax)

        if ax_max > 0: ax.set_ylim(0, ax_max)
        
        leg_cols = 4 if num_substances > 15 else (3 if num_substances > 6 else 2)
        
        ax.legend(
            loc='upper center', 
            bbox_to_anchor=(0.5, -0.12),
            ncol=leg_cols, 
            frameon=True, 
            fontsize=11,          
            title="Badana substancja / Grupa",
            title_fontsize=12,
            framealpha=1.0,
            edgecolor='#cccccc',
            borderpad=1.5,
            labelspacing=0.6,
            columnspacing=2.0
        )

        fig.subplots_adjust(bottom=0.35, top=0.93, left=0.08, right=0.98)
        return fig

    def draw_pca(self, pca_data):
        (pca_df, explained_variance) = pca_data
        
        fig = plt.Figure(figsize=(8, 6), dpi=100)
        ax = fig.add_subplot(111)
        
        # Scatter plot
        sns.scatterplot(
            data=pca_df, x='PC1', y='PC2', 
            hue='Bakteria', style='Bakteria', 
            s=120, ax=ax, palette=self.config["palette"]
        )
        
        # Labels
        for i, row in pca_df.iterrows():
            ax.text(row['PC1']+0.05, row['PC2']+0.05, row['Bakteria'], fontsize=9, alpha=0.8)
            
        # Axes
        ax.set_xlabel(f"PC1 ({explained_variance[0]:.1%})", fontsize=self.config["font_title"])
        ax.set_ylabel(f"PC2 ({explained_variance[1]:.1%})", fontsize=self.config["font_title"])
        ax.set_title("PCA: Podobieństwo profili wrażliwości", fontsize=self.config["font_title"]+2)
        
        ax.grid(True, linestyle='--', alpha=0.5)
        # Add zero lines
        ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(0, color='gray', linestyle=':', alpha=0.5)
        
        fig.tight_layout()
        return fig

