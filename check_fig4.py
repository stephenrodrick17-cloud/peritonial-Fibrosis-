"""
check_fig4.py - independent re-check of Figure 4 against the raw data.
Run from the project root (d:/Peritoneal Project):   python check_fig4.py
Nothing here reads the plotting script; every number is recomputed from raw inputs
and compared with what is printed on the figure. Lines starting FLAG need fixing.
"""
import sys, warnings, traceback
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.feature_selection import RFECV
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ---------------- paths (edit if yours differ) ----------------
P_EXPR = "results/tables/GSE125498_full_expression_matrix.csv"
P_META = "results/tables/GSE125498_sample_metadata.csv"
P_POOL = "results/revision/wgcna_paper_figures/matrisome_deg_module_intersection.csv"
P_MESO_CNT = "data/processed/pseudobulk/pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv"
P_MESO_META = "data/processed/pseudobulk/pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv"
P_ALL_LOG = "results/revision/effluent_validation/GSE248762/analysis_20261009/edger_logCPM.csv"
P_ALL_META = "results/revision/effluent_validation/GSE248762/analysis_20261009/all_cell_pseudobulk_metadata.csv"

# ---------------- what the figure currently shows ----------------
HUB = ["TNFSF15", "LTB", "EBI3", "CCBE1", "NAT6", "BMP7", "SFRP4", "TGM3", "COL8A1", "CLEC2D"]
ENS = {"TNFSF15": "ENSG00000181634", "LTB": "ENSG00000227507", "CCBE1": "ENSG00000183287",
       "NAT6": "ENSG00000243477", "EBI3": "ENSG00000105246", "CLEC2D": "ENSG00000104921",
       "TGM3": "ENSG00000125798", "BMP7": "ENSG00000153162", "SFRP4": "ENSG00000106483",
       "COL8A1": "ENSG00000144810"}
FIG_DISC_AUC = {"TNFSF15": .904, "LTB": .804, "EBI3": .792, "CCBE1": .788, "NAT6": .788,
                "BMP7": .727, "SFRP4": .700, "TGM3": .673, "COL8A1": .581, "CLEC2D": .465}
FIG_VER = {"TNFSF15": (.317, .37, .90), "LTB": (.650, .37, .90), "EBI3": (.450, .55, 1.00),
           "CCBE1": (.017, .00, .25), "NAT6": (.500, .22, .80)}          # (AUC, CI lo, CI hi) from Fig 4F legend
FIG_STARS_G = {"TNFSF15": "***", "LTB": "**", "EBI3": "**", "CCBE1": "**", "NAT6": "**",
               "BMP7": "*", "SFRP4": "ns", "TGM3": "ns", "COL8A1": "ns", "CLEC2D": "ns"}
FIG_STARS_H = {"TNFSF15": "ns", "LTB": "ns", "EBI3": "ns", "CCBE1": "**", "BMP7": "**",
               "SFRP4": "ns", "TGM3": "ns", "COL8A1": "ns", "CLEC2D": "ns"}
FIG_RF10 = ["FLT3LG", "TNFSF15", "ADAM19", "SERPINA10", "LTB", "SEMA3E", "CST7", "C1QB", "HAPLN3", "LTA"]

D = {}


def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < .05 else "ns"


def sec(t):
    print(f"\n{'=' * 10} {t} {'=' * 10}")


def flag(ok, msg):
    print(("  OK    " if ok else "  FLAG  ") + msg)


def auc_ci(y, s, n=2000, seed=1):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) < 2:
            continue
        out.append(roc_auc_score(y[i], s[i]))
    return np.percentile(out, [2.5, 97.5])


def stats2(df, g, a, b):
    x = df.loc[g, a].astype(float).values
    z = df.loc[g, b].astype(float).values
    p = mannwhitneyu(x, z, alternative="two-sided").pvalue
    yy = np.r_[np.ones(len(x)), np.zeros(len(z))]
    ss = np.r_[x, z]
    return x.mean() - z.mean(), np.median(x) - np.median(z), p, roc_auc_score(yy, ss), yy, ss


# ---------------- loading ----------------
def load():
    expr = pd.read_csv(P_EXPR).set_index("Gene")
    expr.columns = expr.columns.astype(str)
    meta = pd.read_csv(P_META)
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta = meta.set_index("sample_id")
    ids = [s for s in meta.index if s in expr.columns]
    grp = meta.loc[ids, "group"]
    D.update(expr=expr, ids=ids, y=(grp == "LPD").astype(int).values,
             lpd=[s for s in ids if grp[s] == "LPD"], spd=[s for s in ids if grp[s] == "SPD"])

    cnt = pd.read_csv(P_MESO_CNT, index_col=0).T          # genes x samples
    mm = pd.read_csv(P_MESO_META)
    lcpm = np.log2(cnt.divide(cnt.sum(axis=0), axis=1) * 1e6 + 1)
    D.update(meso=lcpm,
             meso_lv=[r.donor_id for r in mm.itertuples() if "LV" in str(r.group)],
             meso_sv=[r.donor_id for r in mm.itertuples() if str(r.group) == "SV"])

    al = pd.read_csv(P_ALL_LOG, index_col=0)
    am = pd.read_csv(P_ALL_META).drop_duplicates("sample_id").set_index("sample_id")
    D.update(allc=al,
             allc_lv=[s for s in am.index if am.loc[s, "primary_group"] == "LV"],
             allc_sv=[s for s in am.index if am.loc[s, "primary_group"] == "SV"])


def c_groups():
    flag(len(D["lpd"]) == 13 and len(D["spd"]) == 20,
         f"discovery LPD={len(D['lpd'])}, SPD={len(D['spd'])} (figure: 13/20)")
    flag(len(D["meso_lv"]) == 10 and len(D["meso_sv"]) == 6,
         f"mesothelial pseudobulk LV={len(D['meso_lv'])}, SV={len(D['meso_sv'])} (figure: 10/6)")
    flag(len(D["allc_lv"]) == 10 and len(D["allc_sv"]) == 6,
         f"all-cell pseudobulk LV={len(D['allc_lv'])}, SV={len(D['allc_sv'])}")


# ---------------- Panels E and G ----------------
def c_discovery():
    rows = []
    for g in HUB:
        md, mdn, p, auc, *_ = stats2(D["expr"], g, D["lpd"], D["spd"])
        rows.append(dict(gene=g, mean_diff=md, median_diff=mdn, p=p, raw_auc_LPDhigher=auc,
                         expected_fig_auc=auc if md > 0 else 1 - auc))
    df = pd.DataFrame(rows)
    df["q_BH"] = multipletests(df["p"], method="fdr_bh")[1]
    df["stars_raw_p"] = df["p"].map(stars)
    df["stars_fig"] = df["gene"].map(FIG_STARS_G)
    df["fig_auc"] = df["gene"].map(FIG_DISC_AUC)
    print(df.round(4).to_string(index=False))
    for r in df.itertuples():
        flag(abs(r.fig_auc - r.expected_fig_auc) < .0015,
             f"E legend {r.gene}: figure {r.fig_auc} vs recomputed {r.expected_fig_auc:.3f}")
        flag(r.stars_raw_p == r.stars_fig, f"G stars {r.gene}: figure {r.stars_fig} vs raw-p {r.stars_raw_p}")
        if np.sign(r.mean_diff) != np.sign(r.median_diff):
            flag(False, f"{r.gene}: mean and median differences disagree in sign -> direction is ambiguous")
        if r.stars_raw_p != "ns" and r.q_BH >= .05:
            flag(False, f"{r.gene}: starred on raw p={r.p:.4f} but BH q={r.q_BH:.3f} (not significant after FDR)")


# ---------------- Panels F and H ----------------
def ver_table(df, lv, sv, key):
    up = {g: D["expr"].loc[g, D["lpd"]].astype(float).mean() > D["expr"].loc[g, D["spd"]].astype(float).mean()
          for g in HUB}
    rows = []
    for g in HUB:
        k = key(g)
        if k not in df.index:
            rows.append(dict(gene=g, note=f"absent ({k})"))
            continue
        md, mdn, p, auc, yy, ss = stats2(df, k, lv, sv)
        a = auc if up[g] else 1 - auc
        s = ss if up[g] else -ss
        lo, hi = auc_ci(yy, s)
        rows.append(dict(gene=g, mean_diff_LV_SV=md, p=p, disc_consistent_auc=a, ci_lo=lo, ci_hi=hi,
                         frac_zero=float((ss == 0).mean())))
    out = pd.DataFrame(rows)
    ok = out["p"].notna() if "p" in out else pd.Series(False, index=out.index)
    out["q_BH"] = np.nan
    if ok.any():
        out.loc[ok, "q_BH"] = multipletests(out.loc[ok, "p"], method="fdr_bh")[1]
        out["stars_raw_p"] = [stars(p) if p == p else "" for p in out["p"]]
    return out


def c_verification():
    meso = ver_table(D["meso"], D["meso_lv"], D["meso_sv"], lambda g: g)
    allc = ver_table(D["allc"], D["allc_lv"], D["allc_sv"], lambda g: ENS[g])
    print("\n-- mesothelial pseudobulk (used in Fig 4H) --")
    print(meso.round(4).to_string(index=False))
    print("\n-- all-cell pseudobulk (sensitivity) --")
    print(allc.round(4).to_string(index=False))

    print()
    for g, (a, lo, hi) in FIG_VER.items():
        flag(lo <= a <= hi, f"F legend {g}: AUC {a} lies inside its own CI [{lo}, {hi}]")
        row = meso[meso.gene == g]
        if "disc_consistent_auc" in row and row["disc_consistent_auc"].notna().any():
            v = float(row["disc_consistent_auc"].iloc[0])
            flag(abs(v - a) < .01, f"F legend {g}: figure {a} vs mesothelial recomputed {v:.3f}")
        else:
            flag(False, f"F legend {g}: gene is ABSENT from mesothelial data, so this curve must come from another dataset")
    n_f, n_h = len(FIG_VER), len(FIG_STARS_H)
    flag(n_f == n_h, f"genes in Panel F ({n_f}) vs Panel H ({n_h}); both should show every evaluable gene")

    for r in meso.dropna(subset=["p"]).itertuples():
        flag(r.stars_raw_p == FIG_STARS_H.get(r.gene), f"H stars {r.gene}: figure {FIG_STARS_H.get(r.gene)} vs raw-p {r.stars_raw_p}")


def c_nat6():
    for s in ("NAT6", "NAA80", "NAA50"):
        print(f"  mesothelial matrix has symbol {s}: {s in D['meso'].index}")
    print(f"  all-cell matrix has {ENS['NAT6']}: {ENS['NAT6'] in D['allc'].index}")
    print(f"  first mesothelial row names: {list(D['meso'].index[:5])}")
    print("  NAT6 is the former symbol of NAA80; NAA50 is a different gene.")


# ---------------- Panel A (LASSO) ----------------
def c_lasso():
    pool = pd.read_csv(P_POOL)
    gcol = next(c for c in pool.columns if c.lower() in ("gene", "symbol", "gene_symbol", "genes", "gene.symbol", "annotated gene"))
    mcol = next((c for c in pool.columns if "module" in c.lower()), None)
    if mcol:
        pool = pool[pool[mcol].notna()]
    genes = [g for g in pool[gcol] if g in D["expr"].index]
    flag(len(genes) == 43, f"candidate pool size = {len(genes)} (report says 43)")
    X = StandardScaler().fit_transform(D["expr"].loc[genes, D["ids"]].T.values.astype(float))
    D.update(genes=genes, X=X)
    y, n = D["y"], len(D["y"])
    Cs = np.logspace(-2, 2.5, 40)
    ll = -np.log(Cs * n)                 # glmnet-style log(lambda) for sklearn C (mean-loss scaling)

    def fit(C):
        return LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=5000).fit(X, y)

    for K in (5, 10):
        errs = np.zeros((K, len(Cs)))
        for k, (tr, te) in enumerate(StratifiedKFold(K, shuffle=True, random_state=42).split(X, y)):
            for j, C in enumerate(Cs):
                m = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=5000).fit(X[tr], y[tr])
                errs[k, j] = 1 - m.score(X[te], y[te])
        mean, se = errs.mean(0), errs.std(0, ddof=1) / np.sqrt(K)
        jmin = int(mean.argmin())
        j1 = int(np.where(mean <= mean[jmin] + se[jmin])[0][0])      # smallest C = most regularised
        nz = lambda j: int((fit(Cs[j]).coef_ != 0).sum())
        print(f"  {K}-fold: lambda.min log={ll[jmin]:.2f} err={mean[jmin]:.3f} nonzero={nz(jmin)} | "
              f"lambda.1se log={ll[j1]:.2f} nonzero={nz(j1)}")
        flag(ll[j1] >= ll[jmin], f"{K}-fold: lambda.1se must be >= lambda.min on the log scale")

    print("  nonzero coefficients at the log(lambda) values printed on Fig 4A:")
    for L in (-8, -7.0, -5.4, -4, -2):
        print(f"    log(lambda)={L}: {int((fit(np.exp(-L) / n).coef_ != 0).sum())}")
    print(f"  nonzero at sklearn C=5.0 (the setting the report cites): {int((fit(5.0).coef_ != 0).sum())}")
    print("  Fig 4A claims: 17 genes at the marker labelled lambda.1se (-7.0), 11 at lambda.min (-5.4),"
          " and the panel title says 10-fold CV while the report says cv=5.")


# ---------------- Panels B and C ----------------
def gran(name, val, n):
    k = val * n
    flag(abs(k - round(k)) < .1, f"{name}: {val} x {n} = {k:.2f} events (must be an integer for one CV/OOB run)")


def c_rf_svm():
    X, y, genes = D["X"], D["y"], np.array(D["genes"])
    rf = RandomForestClassifier(n_estimators=20, warm_start=True, oob_score=True,
                                class_weight="balanced", max_features="sqrt", random_state=42)
    oob, e0, e1 = [], [], []
    for t in range(20, 501, 10):
        rf.set_params(n_estimators=t).fit(X, y)
        pr = np.argmax(np.nan_to_num(rf.oob_decision_function_), axis=1)
        oob.append((pr != y).mean())
        e0.append(int((pr[y == 0] != 0).sum()))
        e1.append(int((pr[y == 1] != 1).sum()))
    d = np.diff(oob[10:])
    print(f"  RF OOB at 500 trees: {oob[-1]:.3f} | SPD errors {e0[-1]}/20 | LPD errors {e1[-1]}/13")
    print(f"  OOB step-to-step noise after 200 trees (std of diff): {d.std():.4f}; upward steps: {(d > 0).sum()}")
    print("  A real OOB curve is a noisy staircase. Fig 4C draws smooth exponential-looking curves.")
    gran("Fig 4C final SPD error (approx. read from plot)", 0.06, 20)
    gran("Fig 4C final LPD error (approx. read from plot)", 0.088, 13)
    gran("Fig 4B SVM-RFE accuracy label", 0.954, 33)
    gran("Fig 4B SVM-RFE error label", 0.046, 33)

    m = RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=42).fit(X, y)
    order = [str(g) for g in genes[np.argsort(m.feature_importances_)[::-1]]]
    top10 = order[:10]
    print(f"\n  recomputed RF top 10: {top10}")
    print(f"  figure RF top 10    : {FIG_RF10}")
    print(f"  overlap = {len(set(top10) & set(FIG_RF10))}/10")
    for g in ("CCBE1", "NAT6"):
        r = order.index(g) + 1 if g in order else None
        flag(r is not None and r <= 10,
             f"{g} is highlighted in Fig 4C as an RF-selected gene; recomputed RF rank = {r} (must be <= 10; "
             f"the Venn says RF and SVM-RFE share no genes)")
    cnt = pd.Series(0, index=genes)
    for s in range(20):
        mm = RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=s).fit(X, y)
        cnt.loc[genes[np.argsort(mm.feature_importances_)[::-1][:10]]] += 1
    print("\n  RF top-10 selection frequency over 20 seeds:")
    print(cnt.sort_values(ascending=False).head(15).to_string())

    rfe = RFECV(SVC(kernel="linear", C=1.0), step=1, min_features_to_select=1, scoring="accuracy",
                cv=StratifiedKFold(10, shuffle=True, random_state=42)).fit(X, y)
    try:
        acc = np.asarray(rfe.cv_results_["mean_test_score"])
    except Exception:
        acc = np.asarray(rfe.grid_scores_).mean(1)
    turns = int((np.diff(acc)[1:] * np.diff(acc)[:-1] < 0).sum())
    print(f"\n  SVM-RFE: optimal n={rfe.n_features_}, best CV acc={acc.max():.3f}, direction changes in curve={turns}")
    print(f"  SVM-RFE selected: {[str(g) for g in genes[rfe.support_]]}")
    print("  Fig 4B shows ONE smooth peak at 10 features, acc 0.954; compare with the numbers above.")


if __name__ == "__main__":
    for name, fn in [("Load data", load), ("Group sizes", c_groups),
                     ("Panels E and G: discovery AUC, stars, FDR", c_discovery),
                     ("Panels F and H: verification", c_verification),
                     ("NAT6 lookup", c_nat6),
                     ("Panel A: LASSO path / lambda markers", c_lasso),
                     ("Panels B and C: RF and SVM-RFE re-fit", c_rf_svm)]:
        sec(name)
        try:
            fn()
        except Exception:
            traceback.print_exc()
