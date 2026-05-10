#Calculate Statistics
import pandas as pd
import glob
import os

POSITIVE_FAMILY = "hemoglobin"

# Load answer key
answer = pd.read_csv("answer_key.tsv", sep="\t")

positives = set(answer.loc[
    answer["family"] == POSITIVE_FAMILY,
    "subject_id"
])

negatives = set(answer.loc[
    answer["family"] != POSITIVE_FAMILY,
    "subject_id"
])

rows = []

# Loop through all BLAST result files
for file in glob.glob("results/*.tsv"):

    blast = pd.read_csv(
        file,
        sep="\t",
        names=[
            "qseqid",
            "sseqid",
            "evalue",
            "bitscore",
            "pident",
            "length"
        ]
    )

    # Unique hits
    hit_ids = set(blast["sseqid"].unique())

    # Statistics
    TP = len(hit_ids & positives)

    FP = len(hit_ids & negatives)

    FN = len(positives - hit_ids)

    TN = len(negatives - hit_ids)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0

    recall = TP / (TP + FN) if (TP + FN) > 0 else 0

    rows.append({
        "test": os.path.basename(file),
        "hits": len(hit_ids),
        "TP": TP,
        "FP": FP,
        "FN": FN,
        "TN": TN,
        "precision": precision,
        "recall": recall,
        "avg_evalue": blast["evalue"].mean(),
        "avg_bitscore": blast["bitscore"].mean(),
        "avg_percent_identity": blast["pident"].mean()
    })

summary = pd.DataFrame(rows)

summary = summary.sort_values("test")

summary.to_csv(
    "results/summary_stats.tsv",
    sep="\t",
    index=False
)

print(summary)
