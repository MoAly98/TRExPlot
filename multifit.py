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
import matplotlib.lines as mlines

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
    parser.add_argument('-m','--mergecat',  nargs = '+', default=["DataDriven:Instrumental_lumi:Instrumental:Instrumental_JESR:Instrumental_FTAG"], help="Categories to merge")
    parser.add_argument('-pt', '--pullmin', type=float, default = 3e-2,    help="Threshold of pull/cosntr to highlight in pull plot")
    return parser.parse_args()

def  main():

    # ============== Parse Args =================== #
    args = get_args()
    outdir = args.output
    os.makedirs(outdir, exist_ok = True)
    os.makedirs(f"{outdir}/pulls/", exist_ok = True)
    merge_syst_cats = args.mergecat

    # Threshold of pull/constr to highlight in pull plot
    pull_constr_threshold = args.pullmin
    # ======================================================================
    # TREx Config Reading
    # ======================================================================

    # Get TRExFitter config
    trex_config    = args.trex_cfg
    # Get Regions, Samples and OutputDirectory
    job, objs, attrs   = process_config(trex_config, **{'Fit': {'attrs': ['ConfigFile', 'Label']}})
    fits = {objs['Fit'][idx]: { 'config': fit_attr['ConfigFile'], 'label': fit_attr['Label']} for idx, fit_attr in enumerate(attrs['Fit'])}

    data_for_plot = defaultdict(list)
    for fit in fits:
        # Get Regions, Samples and OutputDirectory
        fit_job, fit_objs, fit_attrs   = process_config(fits[fit]['config'], **{'OutputDir': {},
                                                    'BlindingThreshold': {},
                                                    'POI': {},
                                                    'Region': {'attrs': ['Variable', 'Binning', 'Rebinning','VariableTitle', 'DropBins', 'AutomaticDropBins']},
                                                    'Sample': {'attrs': ['FillColor', 'Type', 'Title', 'Group']},
                                                    'Systematic': {'attrs': ['Title', 'Category'], 'split_by_sc':True},
                                                    'Fit': {'attrs': ['FitType', 'FitRegions']},
                                                    'NormFactor': {'attrs': ['Title']},
                                                    })

        trex_outdir = fit_objs['OutputDir'][0]
        if trex_outdir == "":
            logger.error(f"Cannot find an output directory in the config file. This is required.")

        systematics = fit_objs['Systematic']
        syst_attrs = {systematics[idx]: { 'label':    syst_attr['Title'], 'category': syst_attr['Category']} for idx, syst_attr in enumerate(fit_attrs['Systematic'])}
        syst_cats  = set([v['category'] for v in syst_attrs.values()])

        nfs = fit_objs['NormFactor']
        nfs_attrs = {nfs[idx]: { 'label':    nf_attr['Title']} for idx, nf_attr in enumerate(fit_attrs['NormFactor'])}

        # ========  # Pulls and Covariance results
        # =======================
        fit_results_path = f"{trex_outdir}/{fit_job}/Fits/{fit_job}.txt"
        fit_pulls          = read_pulls(fit_results_path)
        fit_corrmat        = read_corr(fit_results_path)
        pois = fit_objs['POI']
        assert len(pois) == 1, "Cannot handle more than 1 POI with this code.."
        nfs = fit_objs['NormFactor']
        kfacts = [nf for nf in nfs if nf not in pois]
        np_pulls  = {k: v for k, v in fit_pulls.items() if k not in kfacts+pois}
        nf_pulls  = {k: v for k, v in fit_pulls.items() if k in kfacts}
        filt_syst_attrs = {k: v for k, v in syst_attrs.items() if k in np_pulls.keys()}

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


            data_for_plot[syst_cat_name].append({'bestfits': bestfits, 'bestfits_unc': bestfits_unc, 'num_nps': num_nps, 'labels': labels, 'fit_label': fits[fit]['label']})

    # ======================================================================
    #  Figure creation
    # ======================================================================
    colours = ['black', 'red', 'blue']
    markers = ['o', 'X', 'v']
    for syst_cat_name, list_of_fits in data_for_plot.items():
        nfits = len(list_of_fits)
        assert nfits < 4, "ERROR: at most 3 fits can be compared"
        all_np_labels = []
        for fit_info in list_of_fits:
            for label in fit_info['labels']:
                if label not in all_np_labels:  all_np_labels.append(label)
        tot_num_nps = len(all_np_labels)

        extra_fig_height = 1 if tot_num_nps > 5 else 1.3
        fig_pulls, ax_pulls = plt.subplots(figsize=(7,  extra_fig_height + tot_num_nps / 4), dpi=100,layout='tight')
        ax_pulls.fill_between([-2, 2], -0.5,tot_num_nps - 0.5, color="#ffd166",zorder=0)#"#ebf5df")
        ax_pulls.fill_between([-1, 1], -0.5, tot_num_nps - 0.5, color="#06d6a0",zorder=1)#"#bad4aa")
        ax_pulls.vlines(0, -0.5, tot_num_nps - 0.5, linestyles="dotted", color="grey", alpha=0.7,zorder=3)

        # ======================================================================
        #  Figure plotting
        # ======================================================================

        fit_colours = {}
        for fit_idx, fit_info in enumerate(list_of_fits):
            bestfits     = fit_info['bestfits']
            bestfits_unc = fit_info['bestfits_unc']
            pull_color  = colours[fit_idx]
            pull_marker = markers[fit_idx]
            fit_label    = fit_info['fit_label']
            fit_colours[fit_label] = (pull_color,pull_marker)

            for np_global_idx, np_label in enumerate(all_np_labels):
                np_ypos_nominal = tot_num_nps - np_global_idx -1

                np_in_all_fits = all([np_label in fit['labels'] for fit in list_of_fits])

                if fit_idx == 0:
                    ax_pulls.hlines(np_ypos_nominal+0.5, -3, 3, color="grey", linewidth=0.3, linestyle='dashed',alpha=0.5)

                if np_label not in fit_info['labels']:  continue
                np_idx_in_this_fit = fit_info['labels'].index(np_label)

                if nfits % 2 == 0:
                    if fit_idx == 0:
                        np_ypos = np_ypos_nominal + 0.35*(fit_idx+1)/nfits
                    elif fit_idx % 2 == 0:
                        np_ypos = np_ypos_nominal + 0.35*(fit_idx)/nfits
                    else :
                        np_ypos = np_ypos_nominal - 0.35*(fit_idx)/nfits
                else:
                    # Calculate the offset, which increases with the distance from fit_idx 0
                    offset = 0.7 * ((fit_idx + 1) // 2) / (nfits)
                    if fit_idx % 2 == 0:
                        np_ypos = np_ypos_nominal + offset
                    else:
                        np_ypos = np_ypos_nominal - offset

                ax_pulls.errorbar(bestfits[np_idx_in_this_fit], np_ypos, xerr=bestfits_unc[np_idx_in_this_fit], mfc='none', mec='none', ecolor=pull_color, linestyle='None',zorder=4)
                ax_pulls.scatter(bestfits[np_idx_in_this_fit], np_ypos, marker=pull_marker, s=30*2/nfits, color=pull_color, linestyle='None', alpha=0.8, facecolors='none',zorder=5)


        # X-axis
        ax_pulls.set_xlim([-3, 3])
        ax_pulls.set_xlabel(r"$\left(\hat{\theta} - \theta_0\right) / \Delta \theta$", fontsize=LABEL_SIZE-2)
        ax_pulls.xaxis.set_minor_locator(AutoMinorLocator())

        # Y-axis
        extra_white_space = 1 if tot_num_nps > 5 else 1
        ax_pulls.set_ylim([-0.5, tot_num_nps+extra_white_space])
        ax_pulls.set_yticks(np.arange(tot_num_nps))
        ax_pulls.set_yticklabels(all_np_labels[::-1])

        # legend
        legend_lines = []
        for fit_label, fit_colour_marker in fit_colours.items():
            line = mlines.Line2D([], [], color=fit_colour_marker[0], marker=fit_colour_marker[1], markersize=6, label=fit_label, markerfacecolor='none')
            legend_lines.append(line)
        ax_pulls.legend(handles=legend_lines, loc="upper right", frameon=False, ncols=min(len(fit_colours), 5))

        # Tick params
        ax_pulls.tick_params(axis="both", which="major", pad=8, labelsize=LABEL_SIZE-4)
        ax_pulls.tick_params(direction="in", top=True, right=True, which="both")

        # # Save it!
        fig_pulls.savefig(f"{outdir}/pulls/pulls_{syst_cat_name}.pdf", dpi=300)

if __name__ == "__main__":
    main()
