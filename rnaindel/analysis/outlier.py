import os
import gzip
import warnings
import pandas as pd

warnings.simplefilter("ignore")
import pickle

def outlier_analysis(df, model_dir):

    _df = df[df["is_rescurable_homopolymer"]].reset_index(drop=True)

    if not len(_df):
        return df

    __df = df[~df["is_rescurable_homopolymer"]].reset_index(drop=True)

    at_ins_df = _df[_df["is_at_ins"] == 1].reset_index(drop=True)
    at_ins_df = find_outliers(at_ins_df, "at_ins", model_dir)

    at_del_df = _df[_df["is_at_del"] == 1].reset_index(drop=True)
    at_del_df = find_outliers(at_del_df, "at_del", model_dir)

    gc_ins_df = _df[_df["is_gc_ins"] == 1].reset_index(drop=True)
    gc_ins_df = find_outliers(gc_ins_df, "gc_ins", model_dir)

    gc_del_df = _df[_df["is_gc_del"] == 1].reset_index(drop=True)
    gc_del_df = find_outliers(gc_del_df, "gc_del", model_dir)

    return pd.concat([__df, at_ins_df, at_del_df, gc_ins_df, gc_del_df], axis=0)


def cov_vaf(row):
    cov = row["ref_count"] + row["alt_count"]
    vaf = row["alt_count"] / cov
    return cov, vaf


def reclassify_by_outlier_status(row):
    if row["outlying"] == -1:
        return "reclassifed_by_outlier_analysis", "somatic"
    else:
        return row["reclassified"], row["predicted_class"]


def find_outliers(df, homopolymer_type, model_dir):
    if not len(df):
        return df

    df["cov"], df["vaf"] = zip(*df.apply(cov_vaf, axis=1))

    saved_model = os.path.join(model_dir, "{}.pkl.gz".format(homopolymer_type))
    iso = pickle.load(gzip.open(saved_model, "rb"))

    _test = df[["alt_count", "vaf", "cov"]]

    pred = pd.DataFrame(data=iso.predict(_test), columns=["outlying"])

    _df = pd.concat([df, pred], axis=1)

    _df["reclassified"], _df["predicted_class"] = zip(
        *_df.apply(reclassify_by_outlier_status, axis=1)
    )

    return _df.drop(["cov", "vaf", "outlying"], axis=1)
