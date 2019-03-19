#!/usr/bin/env python3

import os


def reporter(
    indel_class,
    ds_beta,
    ds_f_beta,
    ds_precision,
    artifact_ratio,
    fs_beta,
    fs_f_beta,
    fs_precision,
    selected_features,
    pt_beta,
    pt_f_beta,
    pt_precision,
    max_features,
    log_dir,
):
    """Output training results

    Args:
        indel_class (str): s for single-nucleotide indels, m for multi-nucleotide indels
        *_beta (int): specify F beta scores optimized in each step: downsampling (ds), feature selection (fs), and parameter tuning (pt)
        *_f_beta (float): F beta scores optimized in each step
        *_precision (float): associated precisions in each step
        artifact_ratio (int): optimized downsampling ratio for artifact class
        selected_features (str): semicolon-delimited list of selected fetures
        max_features (int): optimized maximum num of features considered in sklearn random forest
        log_dir (str): path to ouput dir
    Returns:
        None
    """
    prefix = (
        "single_nucleotide_indel_" if indel_class == "s" else "multi_nucleotide_indel_"
    )
    report_file_name = prefix + "training_result.txt"
    report = os.path.join(log_dir, report_file_name)

    f = open(report, "w")
    ds_header = (
        "#downsampling step\nartifact ratio\t"
        + format_f_beta(ds_beta)
        + "\ttrp\tprecision\n"
    )
    ds_result = (
        "\t".join(
            [
                str(artifact_ratio),
                str(ds_f_beta),
                calculate_tpr(ds_beta, ds_f_beta, ds_precision),
                str(ds_precision),
            ]
        )
        + "\n"
    )

    fs_header = (
        "#feature selection step\nselected features\t"
        + format_f_beta(fs_beta)
        + "\ttrp\tprecision\n"
    )
    fs_result = (
        "\t".join(
            [
                selected_features,
                str(fs_f_beta),
                calculate_tpr(fs_beta, fs_f_beta, fs_precision),
                str(fs_precision),
            ]
        )
        + "\n"
    )

    pt_header = (
        "#parameter tuning step\nmax features\t"
        + format_f_beta(pt_beta)
        + "\ttrp\tprecision\n"
    )
    pt_result = (
        "\t".join(
            [
                str(max_features),
                str(pt_f_beta),
                calculate_tpr(pt_beta, pt_f_beta, pt_precision),
                str(pt_precision),
            ]
        )
        + "\n"
    )

    f.write(ds_header + ds_result + fs_header + fs_result + pt_header + pt_result)
    f.write(
        '\nthe metrics from "parameter tuning step" represent the final model performance.\n'
    )
    f.close()


def format_f_beta(beta):
    return "tpr" if beta > 100 else "f" + str(beta)


def calculate_tpr(beta, f_beta, precision):
    """Calculate TPR from F beta abd precision

    Args:
        beta (int): specify F beta score to be optimized
        f_beta (float): F beta score
        precision (float): associated presicion
    Returns:
        TPR (str): not float
    """
    denom = beta * beta * precision * f_beta
    numtr = (1 + beta * beta) * precision - f_beta

    if numtr == 0:
        return str(0)
    else:
        return str(denom / numtr)
