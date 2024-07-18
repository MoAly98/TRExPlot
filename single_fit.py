import argparse
import os
from collections import defaultdict
import numpy as np
import itertools
from tools.fit_results_reader  import read_pulls, read_corr, read_ranking, read_CovErrDecomp
from tools.config_reader import process_config
from tools.utils import rgx_match, extract_hash_substring, texify
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import atlasify
plt.rcParams['axes.linewidth'] = 1.3
#plt.rcParams['font.size'] = 12
plt.rcParams['xtick.major.pad']='10'
plt.rcParams['ytick.major.pad']='10'
plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

LINEWIDTH = 2
LABEL_SIZE = 16
MAJOR_TICK_LENGTH = 8
MINOR_TICK_LENGTH = 4

def get_args():
    '''
     Method to retreive arguments from command line
     Returns:
        An `ArgumentParser` object with the parsed values
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trex_cfg',                         help="TRExFitter configuration file")
    parser.add_argument('-o','--output',   required = True, help="Folder to dump plots to")
    parser.add_argument('-m','--mergecat',  nargs = '+',    help="Categories to merge")
    parser.add_argument('--impacts',        action = 'store_true',    help="Run impacts")
    parser.add_argument('--pulls',          action = 'store_true',    help="Run pulls")
    parser.add_argument('--corrmat',        action = 'store_true',    help="Run corrmat")
    parser.add_argument('--pullSignif',     action = 'store_true',    help="Run pull significances")
    parser.add_argument('-pt', '--pullmin', type=float, default = 3e-2,    help="Threshold of pull/cosntr to highlight in pull plot")
    return parser.parse_args()

def  main():

    # ============== Parse Args =================== #
    args = get_args()
    outdir = args.output
    os.makedirs(outdir, exist_ok = True)
    os.makedirs(f"{outdir}/pulls/", exist_ok = True)
    os.makedirs(f"{outdir}/correlations/", exist_ok = True)
    os.makedirs(f"{outdir}/rankings/", exist_ok = True)
    os.makedirs(f"{outdir}/pulls_signif/", exist_ok = True)
    os.makedirs(f"{outdir}/CovarianceRanking/", exist_ok = True)
    merge_syst_cats = args.mergecat
    run_pulls = args.pulls
    run_impacts = args.impacts
    run_corrmat = args.corrmat
    run_pull_signif = args.pullSignif
    run_cov_err_decomp = True

    # Threshold of pull/constr to highlight in pull plot
    pull_constr_threshold = args.pullmin
    # ======================================================================
    # TREx Config Reading
    # ======================================================================

    # Get TRExFitter config
    trex_config    = args.trex_cfg
    # Get Regions, Samples and OutputDirectory
    job, objs, attrs   = process_config(trex_config, **{'OutputDir': {},
                                                        'BlindingThreshold': {},
                                                        'POI': {},
                                                        'Region': {'attrs': ['Variable', 'Binning', 'Rebinning','VariableTitle', 'DropBins', 'AutomaticDropBins']},
                                                        'Sample': {'attrs': ['FillColor', 'Type', 'Title', 'Group']},
                                                        'Systematic': {'attrs': ['Title', 'Category'], 'split_by_sc':True},
                                                        'Fit': {'attrs': ['FitType', 'FitRegions']},
                                                        'NormFactor': {'attrs': ['Title']},
                                                        })

    trex_outdir = objs['OutputDir'][0]
    if trex_outdir == "":
        logger.error(f"Cannot find an output directory in the config file. This is required.")

    systematics = objs['Systematic']
    syst_attrs = {systematics[idx]: { 'label':    syst_attr['Title'], 'category': syst_attr['Category']} for idx, syst_attr in enumerate(attrs['Systematic'])}
    syst_cats  = set([v['category'] for v in syst_attrs.values()])

    nfs = objs['NormFactor']
    nfs_attrs = {nfs[idx]: { 'label':    nf_attr['Title']} for idx, nf_attr in enumerate(attrs['NormFactor'])}

    # ======================================================================
    #  Parsing TREx fit results
    # ======================================================================


    # =======================
    # Pulls and Covariance results
    # =======================
    fit_results_path = f"{trex_outdir}/{job}/Fits/{job}.txt"
    fit_pulls          = read_pulls(fit_results_path)
    fit_corrmat        = read_corr(fit_results_path)
    pois = objs['POI']
    assert len(pois) == 1, "Cannot handle more than 1 POI with this code.."
    nfs = objs['NormFactor']
    kfacts = [nf for nf in nfs if nf not in pois]
    np_pulls  = {k: v for k, v in fit_pulls.items() if k not in kfacts+pois}
    nf_pulls  = {k: v for k, v in fit_pulls.items() if k in kfacts}
    filt_syst_attrs = {k: v for k, v in syst_attrs.items() if k in np_pulls.keys()}

    # =======================
    # Ranking results
    # =======================
    if run_impacts:
        ranking_file_path = f"{trex_outdir}/{job}/Ranking_{pois[0]}.yaml"
        ranking = read_ranking(ranking_file_path)

    # =======================
    # Covariance-based Error decomposition results
    # =======================
    if run_cov_err_decomp:
        cov_err_decomp_results_path = f"{trex_outdir}/{job}/Fits/{job}_errDecomp_{pois[0]}.txt"
        total_impacts, nps_impacts = read_CovErrDecomp(cov_err_decomp_results_path)
        nps_impacts = {k: v for k, v in sorted(nps_impacts.items(), key=lambda item: item[1]['impact_symm'], reverse=True)}
        nps_impacts = dict(itertools.islice(nps_impacts.items(), 15))
        cov_err_decomp_pulls = {}
        for np_name in nps_impacts.keys():
            cov_err_decomp_pulls[np_name] = np_pulls[np_name] if np_name in np_pulls else nf_pulls[np_name] if np_name in nf_pulls else "None"

    # ======================================================================
    # PLOTTING
    # ======================================================================

    # ======================================================================
    #  Create ranking plot
    # ======================================================================
    if run_impacts:
        labels = [filt_syst_attrs[k]['label'] if k in  filt_syst_attrs else nfs_attrs[k]['label'] if k in  nfs_attrs else 'None' for k in ranking.keys()]
        labels = texify(labels)
        num_pars = len(list(ranking.keys()))
        legend_space = 1.0 / (num_pars + 3) + 1
        fig, ax_pulls = plt.subplots(figsize=(6, 2.5 + num_pars * 0.45), dpi=100, layout=None)
        ax_impact = ax_pulls.twiny()
        ax_pulls.set_zorder(1)
        ax_pulls.patch.set_visible(False)
        # lines through pulls of -1, 0, 1 for orientation
        ax_pulls.vlines(
                -1, -1, num_pars - 0.5, linestyles="dashed", color="black", linewidth=0.75
            )
        ax_pulls.vlines(
            0, -1, num_pars - 0.5, linestyles="dashed", color="black", linewidth=0.75
        )
        ax_pulls.vlines(
            1, -1, num_pars - 0.5, linestyles="dashed", color="black", linewidth=0.75
        )
        y_pos = np.arange(num_pars)[::-1]

        # pre-fit up
        impact_prefit_up = [v['impact_hi_pre'] for v in ranking.values()]
        pre_up = ax_impact.barh( y_pos, impact_prefit_up, fill=False, linewidth=1, edgecolor="C0")
        # pre-fit down
        impact_prefit_do= [v['impact_lo_pre'] for v in ranking.values()]
        pre_down = ax_impact.barh(y_pos, impact_prefit_do, fill=False, linewidth=1, edgecolor="C5")
        # post-fit up
        impact_postfit_up= [v['impact_hi_post'] for v in ranking.values()]
        post_up = ax_impact.barh(y_pos, impact_postfit_up, color="C0")
        # post-fit down
        impact_postfit_do= [v['impact_lo_post'] for v in ranking.values()]
        post_down = ax_impact.barh(y_pos, impact_postfit_do, color="C5")
        # pulls
        bestfit = [v['pull'] for  v in ranking.values()]
        up_err  = [v['unc_hi'] for v in ranking.values()]
        do_err  = [v['unc_lo'] for v in ranking.values()]

        pulls = ax_pulls.errorbar(bestfit, y_pos, xerr=[np.abs(do_err), up_err], fmt="o", color="k")

        ax_pulls.set_xlabel(r"$\left(\hat{\theta} - \theta_0\right) / \Delta \theta$", fontsize=LABEL_SIZE, labelpad=7)
        ax_impact.set_xlabel(r"$\Delta \mu$", fontsize=LABEL_SIZE, labelpad=7)
        ax_pulls.set_xlim([-2, 2])
        ax_impact.set_xlim([-5, 5])
        ax_pulls.set_ylim([-1, num_pars])
        impact_max = np.amax(
            np.fabs(
                np.hstack(
                    (
                        impact_prefit_up,
                        impact_prefit_do,
                        impact_postfit_up,
                        impact_postfit_do,
                    )
                )
            )
        )
        ax_impact.set_xlim([-impact_max * 1.1, impact_max * 1.1])
        for axis in [ax_pulls.xaxis, ax_impact.xaxis]:
            axis.set_minor_locator(AutoMinorLocator())

        ax_pulls.set_yticks(y_pos)
        ax_pulls.set_yticklabels(labels, fontsize=LABEL_SIZE)

        ax_pulls.tick_params(direction="in", which="both", labelsize=LABEL_SIZE)
        ax_impact.tick_params(direction="in", which="both", labelsize=LABEL_SIZE, pad=0.8)
        fig.legend(
            (pre_up, pre_down, post_up, post_down, pulls),
            (
                r"pre-fit impact: ($+ \Delta \theta$)", #$\theta = \hat{\theta} + \Delta \theta$",
                r"pre-fit impact: ($- \Delta \theta$)",
                r"post-fit impact: ($+ \Delta \hat{\theta}$)",
                r"post-fit impact: ($- \Delta \hat{\theta}$)",
                "pulls",
            ),
            frameon=False,
            loc="upper left",
            ncol=3,
            fontsize=LABEL_SIZE,
        )
        fig.savefig(f"{outdir}/rankings/Ranking.pdf", bbox_inches='tight')

    # ======================================================================
    #  Create ranking plot
    # ======================================================================
    if run_cov_err_decomp:
        nps_labels = [filt_syst_attrs[k]['label'] if k in  filt_syst_attrs else nfs_attrs[k]['label'] if k in nfs_attrs else 'None' for k in nps_impacts.keys()]

        nps_labels = texify(nps_labels)
        num_nps = 15
        num_totals = 4
        num_pars = num_nps + num_totals # Only plot top 20 parameters + 4 totals
        fig, ax_pulls = plt.subplots(figsize=(6, 2.5 + num_pars * 0.45), dpi=100, layout=None)
        ax_impact = ax_pulls.twiny()
        ax_pulls.set_zorder(1)
        ax_pulls.patch.set_visible(False)
        # lines through pulls of -1, 0, 1 for orientation
        ax_pulls.vlines(
                -1, -1, num_pars - 0.5, linestyles="dashed", color="black", linewidth=0.75
            )
        ax_pulls.vlines(
            0, -1, num_pars - 0.5, linestyles="dashed", color="black", linewidth=0.75
        )
        ax_pulls.vlines(
            1, -1, num_pars - 0.5, linestyles="dashed", color="black", linewidth=0.75
        )

        y_pos = np.arange(num_pars)[::-1]
        # Totals labels
        totals_labels = ["Total error", "Systematic error", "Data Stat error", "MC + Data Stat error"]
        # post-fit up
        impact_postfit_up= [v['impact_up'] for v in nps_impacts.values()]
        impact_postfit_up.insert(0, total_impacts["TOTSTAT_ERROR"]['impact_up'])
        impact_postfit_up.insert(0, total_impacts["STAT_ERROR"]['impact_up'])
        impact_postfit_up.insert(0, total_impacts["SYST_ERROR"]['impact_up'])
        impact_postfit_up.insert(0, total_impacts["TOT_ERROR"]['impact_up'])
        post_up = ax_impact.barh(y_pos, impact_postfit_up, color="C0")
        # post-fit down
        impact_postfit_do= [v['impact_do'] for v in nps_impacts.values()]
        impact_postfit_do.insert(0, total_impacts["TOTSTAT_ERROR"]['impact_do'])
        impact_postfit_do.insert(0, total_impacts["STAT_ERROR"]['impact_do'])
        impact_postfit_do.insert(0, total_impacts["SYST_ERROR"]['impact_do'])
        impact_postfit_do.insert(0, total_impacts["TOT_ERROR"]['impact_do'])

        post_down = ax_impact.barh(y_pos, impact_postfit_do, color="C5")
        # pulls
        bestfit =  [v['pull'] for  k, v  in cov_err_decomp_pulls.items()]
        up_err  =  [v['unc_hi'] for k, v in cov_err_decomp_pulls.items()]
        do_err  =  [v['unc_lo'] for k, v in cov_err_decomp_pulls.items()]

        pulls = ax_pulls.errorbar(bestfit, y_pos[num_totals:], xerr=[np.abs(do_err), up_err], fmt="o", color="k")

        ax_pulls.set_xlabel(r"$\left(\hat{\theta} - \theta_0\right) / \Delta \theta$", fontsize=LABEL_SIZE, labelpad=7)
        ax_impact.set_xlabel(r"$\Delta \mu$", fontsize=LABEL_SIZE, labelpad=7)
        ax_pulls.set_xlim([-2, 2])
        ax_impact.set_xlim([-5, 5])
        ax_pulls.set_ylim([-1, num_pars])
        impact_max = np.amax(
            np.fabs(
                np.hstack(
                    (
                        impact_postfit_up,
                        impact_postfit_do,
                    )
                )
            )
        )
        ax_impact.set_xlim([-impact_max * 1.1, impact_max * 1.1])
        for axis in [ax_pulls.xaxis, ax_impact.xaxis]:
            axis.set_minor_locator(AutoMinorLocator())

        ax_pulls.set_yticks(y_pos)
        ax_pulls.set_yticklabels(totals_labels+nps_labels, fontsize=LABEL_SIZE)

        ax_pulls.tick_params(direction="in", which="both", labelsize=LABEL_SIZE)
        ax_impact.tick_params(direction="in", which="both", labelsize=LABEL_SIZE, pad=0.8)
        fig.legend(
            ( post_up, post_down, pulls),
            (
                r"post-fit up   impact:",
                r"post-fit down impact:",
                "pulls",
            ),
            frameon=False,
            bbox_to_anchor=(-0.2, 0.97),
            loc="upper left",
            ncol=3,
            fontsize=LABEL_SIZE,
        )
        fig.savefig(f"{outdir}/CovarianceRanking/covImpacts.pdf", bbox_inches='tight')

    # ======================================================================
    #  Create pull and pul significance plots
    # ======================================================================
    if run_pulls:
        skip = []
        for syst_cat in syst_cats:
            if syst_cat in skip:    continue
            syst_cat_list = [syst_cat]
            syst_cat_name = syst_cat
            if merge_syst_cats is not None:
                for group in merge_syst_cats:
                    if syst_cat in group:
                        syst_cat_list = group.split(":")
                        syst_cat_name = "__".join(syst_cat_list)
                        skip.extend(syst_cat_list)
                        break


            cat_np_pulls = {k: v for k, v in np_pulls.items() if filt_syst_attrs[k]['category'] in syst_cat_list}
            pulls = [(v["pull"], v["constr"]) for k, v in cat_np_pulls.items()]
            bestfits, bestfits_unc = zip(*pulls)
            labels = [filt_syst_attrs[k]['label'] for k in cat_np_pulls.keys()]
            labels = texify(labels)
            num_nps = len(list(cat_np_pulls.keys()))


            # ======================================================================
            #  Create pulls
            # ======================================================================
            fig_pulls, ax_pulls = plt.subplots(figsize=(7, 1 + num_nps / 4), dpi=100,layout='tight')
            y_positions = np.arange(num_nps)[::-1]
            pull_colors = ['royalblue' if bf > pull_constr_threshold or abs(constr-1) > pull_constr_threshold else 'black'  for (bf,  constr) in zip(bestfits, bestfits_unc)]
            ax_pulls.errorbar(bestfits, y_positions, xerr=bestfits_unc, mfc='none', mec='none', ecolor=pull_colors, linestyle='None')
            ax_pulls.fill_between([-2, 2], -0.5, len(bestfits) - 0.5, color="#ffd166")#"#ebf5df")
            ax_pulls.fill_between([-1, 1], -0.5, len(bestfits) - 0.5, color="#06d6a0")#"#bad4aa")
            ax_pulls.vlines(0, -0.5, len(bestfits) - 0.5, linestyles="dotted", color="black")

            ax_pulls.scatter(bestfits, y_positions, marker="o", s=30, color=pull_colors, linestyle='None', alpha=0.8)
            # X-axis
            ax_pulls.set_xlim([-3, 3])
            ax_pulls.set_xlabel(r"$\left(\hat{\theta} - \theta_0\right) / \Delta \theta$", fontsize=LABEL_SIZE-2)
            ax_pulls.xaxis.set_minor_locator(AutoMinorLocator())

            # Y-axis
            ax_pulls.set_ylim([-0.5, num_nps - 0.5])
            ax_pulls.set_yticks(y_positions)
            ax_pulls.set_yticklabels(labels)

            # Tick params
            ax_pulls.tick_params(axis="both", which="major", pad=8, labelsize=LABEL_SIZE-4)
            ax_pulls.tick_params(direction="in", top=True, right=True, which="both")

            # # Save it!
            fig_pulls.savefig(f"{outdir}/pulls/pulls_{syst_cat_name}.pdf", dpi=300, bbox_inches='tight')

    # ======================================================================
    #  Create pull significance plot
    # ======================================================================
    if run_pull_signif:
        pulls = [(v["pull"], v["constr"]) for k, v in np_pulls.items()]
        bestfits, bestfits_unc = zip(*pulls)
        labels = [filt_syst_attrs[k]['label'] for k in np_pulls.keys()]
        labels = texify(labels)

        fig_signif, ax_signif = plt.subplots(figsize=(7, 1 + 20 / 4), dpi=100,layout='tight')
        y_positions = np.arange(20)[::-1]
        #bestfits bestfits_unc
        signif = [np.divide(np.abs(bf), np.sqrt(1-bf_unc**2))
                    if bf_unc != 1 else np.divide(np.abs(bf), np.sqrt(1-0.999**2))
                    for bf, bf_unc in zip(bestfits, bestfits_unc)]
        signif_ordered, labels_ordered = zip(*sorted(zip(signif, labels), reverse=True))
        signif_ordered = signif_ordered[:20]
        labels_ordered = labels_ordered[:20]
        ax_signif.barh( y_positions, signif_ordered, linewidth=1, color="midnightblue", alpha=0.9, edgecolor="grey",zorder=2)
        # X-axis
        ax_signif.set_xlim([0, 3])
        ax_signif.set_xlabel(r"Pull significance", fontsize=LABEL_SIZE-2)
        ax_signif.xaxis.set_minor_locator(AutoMinorLocator())

        # Y-axis
        ax_signif.set_ylim([-0.5, 20 - 0.5])
        ax_signif.set_yticks(y_positions)
        ax_signif.set_yticklabels(labels_ordered)

        # Grid
        ax_signif.grid(True, color='grey', alpha=0.2, axis="x",which='both',zorder=1)

        # Tick params
        ax_signif.tick_params(axis="both", which="major", pad=8, labelsize=LABEL_SIZE-4)
        ax_signif.tick_params(direction="in", top=True, right=True, which="both")
        fig_signif.savefig(f"{outdir}/pulls_signif/pull_signifs_ranked.pdf", dpi=300, bbox_inches='tight')

    # ======================================================================
    #  Create Correlation matrix
    # ======================================================================
    if run_corrmat:
        labels = [filt_syst_attrs[k]['label'] if k in  filt_syst_attrs else nfs_attrs[k]['label'] if k in  nfs_attrs else 'None' for k in fit_pulls.keys()]
        labels = texify(labels)
        corr_mat = np.zeros(shape=(len(labels), len(labels)))
        for i, row in enumerate(fit_corrmat):
            row_itr = row[::-1]
            for j, col in enumerate(row_itr):
                corr_mat[i][j] = col

        fit_corrmat = corr_mat
        below_threshold = np.where(np.abs(fit_corrmat) < 0.25, True, False)
        np.fill_diagonal(below_threshold, True)

        all_below_threshold = np.all(below_threshold, axis=0)
        fixed_parameter = np.all(np.equal(fit_corrmat, 0.0), axis=0)
        delete_indices =np.where(np.logical_or(all_below_threshold, fixed_parameter))

        fit_corrmat = np.delete(np.delete(fit_corrmat, delete_indices, axis=1), delete_indices, axis=0)
        labels = np.delete(labels, delete_indices)

        n_params = fit_corrmat.shape[0]

        fig, ax = plt.subplots(
            figsize=(round(5 + n_params / 1.6, 1), round(3 + n_params / 1.6, 1)),
            dpi=100,
            layout="constrained",
        )
        im = ax.imshow(fit_corrmat, vmin=-1, vmax=1, cmap="RdBu")
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticklabels(labels)

        # add correlation as text
        for (j, i), corr in np.ndenumerate(fit_corrmat):
            text_color = "white" if abs(fit_corrmat[j, i]) > 0.75 else "black"
            if abs(corr) > 0.005:
                if corr != 100:
                    ax.text(i, j, f"{corr*100:.1f}", ha="center", va="center", color=text_color)
                else:
                    ax.text(i, j, f"{corr*100:.0f}", ha="center", va="center", color=text_color)
        fig.colorbar(im, ax=ax)
        ax.set_aspect("auto")  # to get colorbar aligned with matrix
        fig.savefig(f"{outdir}/correlations/corrMat.pdf")

if __name__ == "__main__":
    main()
