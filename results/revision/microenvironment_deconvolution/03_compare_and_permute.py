"""Compare LPD and SPD scores and perform random fixed-size label permutations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


OUT = Path(__file__).resolve().parent
SEED = 20261015
LABEL_PERMUTATIONS = 10_000
ALPHA = 0.05


def compare_scores(scores: np.ndarray, labels: np.ndarray) -> tuple[float, float, float]:
    lpd = scores[labels == "LPD"]
    spd = scores[labels == "SPD"]
    test = mannwhitneyu(lpd, spd, alternative="two-sided", method="auto")
    rank_biserial = 2.0 * float(test.statistic) / (len(lpd) * len(spd)) - 1.0
    return float(test.statistic), float(test.pvalue), rank_biserial


def main() -> None:
    scores = pd.read_csv(OUT / "sample_population_scores_long.csv")
    if scores.empty:
        raise ValueError("No sample-population scores were found.")
    comparisons: list[dict[str, object]] = []
    permutation_results: list[dict[str, object]] = []
    familywise_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(SEED)

    for variant, frame in scores.groupby("signature_variant", sort=False):
        populations = list(frame["population"].drop_duplicates())
        samples = frame[["sample_id", "group"]].drop_duplicates()
        if len(samples) != 33 or (samples["group"] == "LPD").sum() != 13 or (samples["group"] == "SPD").sum() != 20:
            raise ValueError("Expected fixed-size group labels of 13 LPD and 20 SPD.")
        samples = samples.sort_values("sample_id")
        labels = samples["group"].to_numpy()
        observed_by_population: dict[str, dict[str, float]] = {}
        vectors: dict[str, np.ndarray] = {}
        for population in populations:
            cell = frame[frame["population"] == population].set_index("sample_id").loc[samples["sample_id"]]
            values = cell["score_mean_z"].to_numpy(dtype=float)
            if not np.isfinite(values).all():
                raise ValueError(f"Non-finite scores found for {population} / {variant}.")
            vectors[population] = values
            U, p_value, effect = compare_scores(values, labels)
            lpd = values[labels == "LPD"]
            spd = values[labels == "SPD"]
            direction = "higher in LPD" if np.median(lpd) > np.median(spd) else (
                "higher in SPD" if np.median(lpd) < np.median(spd) else "equal medians"
            )
            observed_by_population[population] = {
                "U": U,
                "P": p_value,
                "effect": effect,
                "abs_effect": abs(effect),
                "LPD_median": float(np.median(lpd)),
                "SPD_median": float(np.median(spd)),
            }
            comparisons.append(
                {
                    "signature_variant": variant,
                    "population": population,
                    "LPD_n": len(lpd),
                    "SPD_n": len(spd),
                    "LPD_median": np.median(lpd),
                    "SPD_median": np.median(spd),
                    "direction_by_median": direction,
                    "Mann_Whitney_U_LPD": U,
                    "Mann_Whitney_P_two_sided": p_value,
                    "rank_biserial_correlation": effect,
                    "nominal_P_lt_0_05": p_value < ALPHA,
                }
            )

        begin = len(permutation_results)
        any_nominal = np.zeros(LABEL_PERMUTATIONS, dtype=bool)
        for permutation in range(LABEL_PERMUTATIONS):
            permuted_labels = rng.permutation(labels)
            perm_p_values = []
            for population in populations:
                values = vectors[population]
                perm_U, perm_p, perm_effect = compare_scores(values, permuted_labels)
                observed = observed_by_population[population]
                permutation_results.append(
                    {
                        "signature_variant": variant,
                        "permutation": permutation + 1,
                        "population": population,
                        "rank_biserial_correlation": perm_effect,
                        "mann_whitney_P_two_sided": perm_p,
                    }
                )
                perm_p_values.append(perm_p)
            any_nominal[permutation] = any(p < ALPHA for p in perm_p_values)
            if (permutation + 1) % 1000 == 0:
                print(f"{variant}: completed {permutation + 1}/{LABEL_PERMUTATIONS} permutations")

        variant_permutations = permutation_results[begin:]
        for comparison in comparisons:
            if comparison["signature_variant"] != variant:
                continue
            population = str(comparison["population"])
            observed = observed_by_population[population]
            extreme_count = sum(
                abs(float(row["rank_biserial_correlation"])) >= observed["abs_effect"]
                for row in variant_permutations
                if row["population"] == population
            )
            comparison["empirical_permutation_P"] = (
                extreme_count + 1
            ) / (LABEL_PERMUTATIONS + 1)
        familywise_rows.append(
            {
                "signature_variant": variant,
                "random_label_permutations": LABEL_PERMUTATIONS,
                "permutations_with_any_population_nominal_P_lt_0_05": int(any_nominal.sum()),
                "familywise_permutation_P": float(any_nominal.mean()),
                "percent_with_any_nominal_hit": 100 * float(any_nominal.mean()),
                "fixed_group_sizes": "13 LPD; 20 SPD",
                "seed": SEED,
            }
        )

    comparison_df = pd.DataFrame(comparisons)
    comparison_df["BH_FDR_within_signature_variant"] = comparison_df.groupby(
        "signature_variant"
    )["Mann_Whitney_P_two_sided"].transform(
        lambda p: multipletests(p, method="fdr_bh")[1]
    )
    comparison_df["FDR_lt_0_05"] = comparison_df["BH_FDR_within_signature_variant"] < ALPHA
    comparison_df["permutation_P_lt_0_05"] = comparison_df["empirical_permutation_P"] < ALPHA
    comparison_df.to_csv(OUT / "group_comparisons.csv", index=False)
    pd.DataFrame(familywise_rows).to_csv(OUT / "familywise_permutation_summary.csv", index=False)
    pd.DataFrame(
        permutation_results,
        columns=[
            "signature_variant", "permutation", "population",
            "rank_biserial_correlation", "mann_whitney_P_two_sided",
        ],
    ).to_csv(OUT / "sample_label_permutation_statistics.csv", index=False)
    comparison_df.loc[comparison_df["Mann_Whitney_P_two_sided"] < ALPHA].to_csv(
        OUT / "nominal_P_lt_0_05.csv", index=False
    )
    comparison_df.loc[comparison_df["BH_FDR_within_signature_variant"] < ALPHA].to_csv(
        OUT / "FDR_lt_0_05.csv", index=False
    )
    print(comparison_df.to_string(index=False))
    print(pd.DataFrame(familywise_rows).to_string(index=False))


if __name__ == "__main__":
    main()
