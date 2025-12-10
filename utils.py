import re
import numpy as np
from scipy import stats

# --- SORTOWANIE I PARSOWANIE ---
def smart_sort_key(group_name):
    """Sortowanie naturalne dla nazw grup z liczbami."""
    match = re.match(r"(.+?)\s*\(([\d,.]+)\s*(.+?)\)", group_name)
    if match:
        name = match.group(1).strip()
        val_str = match.group(2).replace(',', '.')
        try: val = float(val_str)
        except ValueError: val = 0.0
        return (name, val)
    else: return (group_name, 0.0)

def parse_concentration(group_name):
    """Wyciąganie stężenia i jednostki z nazwy grupy."""
    match = re.search(r"([\d,.]+)\s*(mg\/ml|ug\/ml|%)", group_name)
    if match:
        val_str = match.group(1).replace(',', '.')
        try:
            conc = float(val_str)
            unit = match.group(2)
            substance = re.sub(r"\s*\(?[\d,.]+\s*(mg\/ml|ug\/ml|%).*\)?", "", group_name).strip()
            return substance, conc, unit
        except ValueError: return None, None, None
    return None, None, None

# --- STATYSTYKA: EFFECT SIZE ---
def calculate_cohens_d(group1_data, group2_data):
    n1, n2 = len(group1_data), len(group2_data)
    if n1 < 2 or n2 < 2: return 0.0
    
    var1 = np.var(group1_data, ddof=1)
    var2 = np.var(group2_data, ddof=1)
    
    s_pooled = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    mean1 = np.mean(group1_data)
    mean2 = np.mean(group2_data)

    if s_pooled == 0:
         return 0.0 
    
    return (mean1 - mean2) / s_pooled

def get_effect_size_interpretation(d):
    d = abs(d)
    if d < 0.2: return "znikomy"
    elif d < 0.5: return "mały"
    elif d < 0.8: return "średni"
    else: return "DUŻY"

# --- STATYSTYKA: OUTLIERS (DIXON LOGIC) ---
def find_outliers_dixon(df):
    """Zwraca listę wykrytych outlierów (logika bez GUI)."""
    dixon_q90 = {3: 0.941, 4: 0.765, 5: 0.642, 6: 0.560, 7: 0.507, 8: 0.468, 9: 0.437, 10: 0.412}
    detected = []
    
    for group in df['Grupa'].unique():
        values = df[df['Grupa'] == group]['Srednica_mm'].values
        values = sorted(values)
        n = len(values)
        if n < 3 or n > 10: continue 
        r = values[-1] - values[0]
        if r == 0: continue
        
        gap_low = values[1] - values[0]
        q_calc_low = gap_low / r
        gap_high = values[-1] - values[-2]
        q_calc_high = gap_high / r
        q_crit = dixon_q90.get(n, 0.941)
        
        if q_calc_low > q_crit: 
            detected.append({'group': group, 'value': values[0], 'others': str(values[1:])})
        if q_calc_high > q_crit: 
            detected.append({'group': group, 'value': values[-1], 'others': str(values[:-1])})
            
    return detected