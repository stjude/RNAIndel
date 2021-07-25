import os
import gzip
import pickle
import pandas as pd
from sklearn.ensemble import IsolationForest

HOMOPOLYEMRS = ["A/T insertion", "A/T deletion", "G/C insertion", "G/C deletion"]
PREFIXES = ["at_ins", "at_del", "gc_ins", "gc_del"]


def train_homolopolymer(df, model_dir):
    n = 4
    dfs = split(df, n)
    for _df, _polymer, _pf in zip(dfs, HOMOPOLYEMRS, PREFIXES):
        model = _model(_df, _polymer)
        if model:
            updater(model, _pf, _polymer, model_dir)


def split(df, n):

    (
        df["is_at_ins_{}_polymer".format(n)],
        df["is_at_del_{}_polymer".format(n)],
        df["is_gc_ins_{}_polymer".format(n)],
        df["is_gc_del_{}_polymer".format(n)],
    ) = zip(*df.apply(_wrapper, n=n, axis=1))

    df = df[df["truth"] == "artifact"].reset_index(drop=True)

    df["cov"], df["vaf"] = zip(*df.apply(cov_vaf, axis=1))

    df_at_ins_n = df[df["is_at_ins_{}_polymer".format(n)]].reset_index(drop=True)
    df_at_del_n = df[df["is_at_del_{}_polymer".format(n)]].reset_index(drop=True)
    df_gc_ins_n = df[df["is_gc_ins_{}_polymer".format(n)]].reset_index(drop=True)
    df_gc_del_n = df[df["is_gc_del_{}_polymer".format(n)]].reset_index(drop=True)

    features = ["alt_count", "vaf", "cov"]

    return (
        df_at_ins_n[features],
        df_at_del_n[features],
        df_gc_ins_n[features],
        df_gc_del_n[features],
    )


def cov_vaf(row):
    cov = row["alt_count"] + row["ref_count"]
    vaf = row["alt_count"] / cov if cov > 0 else 0

    return cov, vaf


def flag_target_polymers(row, n):
    at_ins_n = is_target_polymer(row, n, True, True)
    at_del_n = is_target_polymer(row, n, True, False)
    gc_ins_n = is_target_polymer(row, n, False, True)
    gc_del_n = is_target_polymer(row, n, False, False)

    return at_ins_n, at_del_n, gc_ins_n, gc_del_n


def is_target_polymer(row, n, for_at, for_ins):

    if row["ref_count"] + row["alt_count"] < 5:
        return False

    n = int(n)

    if for_at:
        if for_ins:
            if row["is_at_ins"] and row["repeat"] > n:
                return True
        else:
            if row["is_at_del"] and row["repeat"] > n:
                return True
    else:
        if for_ins:
            if row["is_gc_ins"] and row["repeat"] > n:
                return True
        else:
            if row["is_gc_del"] and row["repeat"] > n:
                return True
    return False


def _model(df, polymer):
    nrows = len(df.index)
    if nrows < 500:
        print(
            "Only {} artifact examples found for {}. At least 500 needed.".format(
                nrows, polymer
            )
        )
        return None

    iso = IsolationForest(random_state=42, contamination=0.05)

    return iso.fit(df)


def updater(model, prefix, polymer, model_dir):
    path = os.path.join(model_dir, prefix + ".pkl.gz")
    model_pkl = gzip.open(path, "wb")
    pickle.dump(model, model_pkl)
    model_pkl.close()

    print("Successfully completed for {} model training".format(polymer))
