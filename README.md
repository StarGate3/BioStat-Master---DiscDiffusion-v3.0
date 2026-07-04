# BioStat Master v3.0

**BioStat Master** is a specialized, modular Python application designed for the statistical analysis and visualization of **Disk Diffusion Assays** (Zone of Inhibition data). It automates the entire workflow from raw Excel data to publication-ready figures and reports, ensuring statistical rigor using standard biostatistical methods.

> **Disclaimer:** this tool compares inhibition-zone diameters statistically (significance tests + effect size) and does **not** compute clinical S/I/R (Susceptible/Intermediate/Resistant) categories per CLSI (M100) or EUCAST breakpoints.

---

## 🚀 Key Features

### 📊 Statistical Analysis
*   **Automated Decision Tree**: Automatically selects between Parametric (ANOVA) and Non-Parametric (Kruskal-Wallis) tests based on Normality (Shapiro-Wilk) and Homogeneity of Variance (Levene’s test).
*   **Post-hoc Corrections**: Supports **Tukey HSD** (for ANOVA) and **Dunn’s Test** (for Kruskal-Wallis) with multiple correction methods:
    *   Holm-Bonferroni (Default)
    *   Benjamini-Hochberg (FDR)
    *   Bonferroni
*   **Effect Size**: Calculates **Cohen’s *d*** for all pairwise comparisons to determine the magnitude of differences.
*   **Outlier Detection**: Implements **Dixon’s Q Test** to identify and suggest removal of technical outliers in small sample sizes ($3 \le n \le 10$).

### 🎨 Scientific Visualization
Generates high-resolution, publication-quality figures using `Matplotlib` and `Seaborn`:
1.  **Main Comparison Plot**: Barplots, Boxplots, or Violinplots with significance asterisks.
2.  **Heatmaps**: Activity heatmaps and P-value significance matrices.
3.  **Dose-Response Trends**: Line plots with Spearman correlation coefficients.
4.  **Effect Size Plot**: Lollipop charts visualizing the strength of differences (Cohen's d).
5.  **Cross-Species Comparison**: Summary view of activity across multiple bacterial strains.
6.  **PCA Analysis**: Principal Component Analysis to cluster bacterial strains based on their sensitivity profiles.

### 📝 Reporting
*   **PDF Reports**: detailed summary including methodology, statistical results, and embedded figures.
*   **Caption Generator**: Automatically generates scientific figure captions (e.g., "Figure 1. Antibacterial activity...") ready for copy-pasting into manuscripts.
*   **Excel Export**: Exports raw data, statistical summaries, and detailed post-hoc results.

---

## 🧫 MIC/MBC Module (Broth Microdilution)

In addition to disk diffusion, BioStat Master supports **MIC (Minimum Inhibitory Concentration)** and **MBC (Minimum Bactericidal Concentration)** analysis from broth microdilution assays, loaded from the same multi-sheet Excel workbook as the diffusion data (sheets `Dane_dyfuzja`, `MIC_wizualny`, `MIC_OD`, `MBC_posiew`, `Kontrole`, `Instrukcja`, `Ustawienia`). **The legacy single-sheet format is still fully supported** — existing diffusion-only files load and analyze exactly as before, with no changes required and no new sheets to add.

*   **Reading modes**:
    *   *Visual* (`MIC_wizualny`): well-by-well growth/no-growth calls ("wzrost"/"brak").
    *   *OD* (`MIC_OD`): optical-density readings interpreted via a **relative threshold** against that run's own growth and sterility controls (not a fixed absolute OD cutoff), so plate-to-plate baseline differences don't bias the call.
    *   *MBC* (`MBC_posiew`): colony-count (CFU) readings after subculture, converted to a % reduction vs. the inoculum (`Inokulum_CFU_t0`); a well is scored "killed" once **reduction ≥ 99.9% (3-log10)** — the standard clinical definition of bactericidal activity, not an author-chosen threshold.
*   **Run validity**: every run is checked against BOTH its growth control (must show adequate growth) AND its sterility control (must show no contamination), **independently** — a strong growth control can never mask a contaminated sterility control.
*   **Run identifiers (`Przebieg`)**: this identifier is **global to the whole workbook**, not scoped to a single sheet — a row in `MIC_wizualny`/`MIC_OD`/`MBC_posiew` is matched to its controls in `Kontrole` purely by this value, regardless of which data sheet it came from. This is intentional: MBC is typically read from a subculture of the *same* physical run as its MIC counterpart, so sharing one `Przebieg` (and therefore one set of controls) between a MIC row and its corresponding MBC row is correct. It does mean `Przebieg` values must be unique across the *entire* file for unrelated runs — restarting the numbering independently in each sheet (e.g. "1", "2"... reused for unconnected runs) will silently match the wrong controls.
*   **Log2 scale**: all MIC/MBC arithmetic (replicate aggregation, medians, group comparisons, the MBC/MIC ratio) happens on the **log2 dilution-index scale**, never on raw concentration — this is what makes "median MIC" and "difference in dilutions" statistically meaningful for a two-fold dilution series.
*   **Replicates**: technical replicates within one biological replicate collapse to a single value via a **high-median rule** (for an even count, the *higher* of the two middle values is kept — never an interpolated number), so censored calls propagate automatically instead of being silently "fixed" to a number. Biological replicates are then summarized per bacterium × substance: median, range, mode, and geometric mean (**never** the arithmetic mean).
*   **Censoring**: values at the edge of the tested range are always reported as `≤`/`>` bounds, never as a bare number, and are ranked accordingly (always lowest/highest) in any statistical comparison.
*   **Group comparisons**: the same three-layer approach as diffusion — Layer 1 (dilution difference + a configurable "meaningful difference" threshold, default ≥2 dilutions) is always available; Layer 2 (Mann-Whitney, or Kruskal-Wallis + Dunn's post-hoc for >2 groups) is computed when the data support it; Layer 3 guards suppress the p-value (with an explicit stated reason) when censoring is too heavy or too few biological replicates exist — a p-value is never fabricated.
*   **MBC/MIC ratio & classification**: the ratio (`d = log2(MBC) − log2(MIC)`, displayed as `2^d`) is always computed and shown, for any dilution factor. The **bactericidal/bacteriostatic label** (`d ≤ 2` → bactericidal; `d ≥ 3` → bacteriostatic), however, is only given when the run used a **two-fold dilution series** — these thresholds are only mathematically sound on a two-fold scale (every possible step is then a whole number, so it never falls between the two cutoffs); for any other dilution factor the ratio is still reported but the label is marked "unavailable" with an explicit reason instead of guessing. When either MIC or MBC is censored, the label is only given if the resulting bound is unambiguous — otherwise the result is reported as **undetermined**, with the tightest ratio bound shown (e.g. "≥8"), never guessed. A structural consistency check flags any case where MBC would have to be smaller than MIC.
*   **n_bio = 1 warning**: exactly like the diffusion module, a result based on a single biological replicate is still computed and shown, but carries an explicit "no biological replication" banner on every affected plot, table, and report — consistent across the whole application.
*   **Outputs**: MIC/MBC distribution plots, MIC↔MBC pair/gap plots (color-coded by classification), cross-substance comparison plots, and a publication-ready summary table — exported into the same PDF report and Excel workbook as the diffusion results (one combined report when a strain has data for both methods).

---

## 🛠️ Installation & Requirements

Ensure you have **Python 3.8+** installed.

1.  **Clone/Download** the repository.
2.  **Install Dependencies**:
    ```bash
    pip install -r requirements.txt
    ```
    *Key libraries: `customtkinter`, `pandas`, `scipy`, `statsmodels`, `seaborn`, `scikit-posthocs`, `reportlab`.*

---

## 💻 Usage

1.  **Run the Application**:
    ```bash
    python main.py
    ```
2.  **Load Data**: Click "1. Wczytaj Excel". The Excel file should be formatted with columns for 'Bakterie' (Bacteria), 'Grupa' (Group/Substance), and 'Srednica_mm' (Zone Diameter).
3.  **Configure**:
    *   Select the bacterial strain to analyze.
    *   Choose a **Reference Group** (Negative Control) for comparisons.
    *   Select a Post-hoc correction method.
4.  **Run Analysis**: Click "URUCHOM ANALIZĘ".
5.  **Explore Results**: Switch between tabs to view different plots and the statistical log.
6.  **Export**: Save figures as high-res PNGs or generate a full PDF report.

---

## 🏗️ Architecture (v3.0 Modular)

The application follows a clean, modular Model-View-Controller (MVC) pattern for maintainability:

*   **`gui.py` (View/Controller)**: Handles the user interface using `customtkinter`. Orchestrates the application flow.
*   **`logic.py` (Model)**: Contains the `StatsEngine`. Pure Python class responsible for all disk-diffusion statistical calculations (Shapiro, Levene, ANOVA/KW, Post-hoc). Independent of the GUI.
*   **`plotting.py` (View)**: Contains the `Plotter` class. Encapsulates all disk-diffusion `matplotlib` figure generation logic.
*   **`mic_logic.py` (Model)**: MIC/MBC well-reading engine, technical→biological replicate aggregation, cross-group comparisons, and the MBC/MIC ratio/classification. Independent of the GUI, mirrors `logic.py`'s role for the MIC/MBC module.
*   **`mic_plotting.py` (View)**: MIC/MBC figure generation (distribution, MIC↔MBC pairs, group comparison) — mirrors `plotting.py`'s role.
*   **`reports.py`**: Builds the PDF report (diffusion and/or MIC/MBC sections) and the MIC/MBC Excel export.
*   **`utils.py`**: Helper functions for outlier detection (Dixon), robust sorting, string parsing, and multi-sheet workbook routing (`route_workbook`).

---

## 🔬 Statistical Methods

The software adheres to standard biostatistical practices for small-sample biological data:
1.  **Normality Check**: Shapiro-Wilk test ($\alpha=0.05$).
2.  **Variance Check**: Levene’s test ($\alpha=0.05$).
3.  **Test Selection**:
    *   **Parametric**: If Normal AND Homogeneous Variances $\rightarrow$ **One-way ANOVA**.
    *   **Non-Parametric**: If Non-Normal OR Unequal Variances $\rightarrow$ **Kruskal-Wallis**.
4.  **Pairwise Comparison**:
    *   ANOVA $\rightarrow$ **Tukey HSD**.
    *   Kruskal-Wallis $\rightarrow$ **Dunn’s Test** (corrected).
