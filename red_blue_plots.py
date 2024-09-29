import argparse
import os
from collections import defaultdict
import numpy as np
import itertools
from tools.fit_results_reader  import read_pulls, read_corr, read_ranking, read_CovErrDecomp
from tools.config_reader import process_config
from tools.utils import rgx_match, extract_hash_substring, texify
from tools.rebin import rebin
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
import atlasify
import uproot
import mplhep
import matplotlib.lines as mlines
import warnings
import gc
'''
Use warning package to suppres mplhep warnings, especially:

UserWarning: Integer weights indicate poissonian data. Will calculate Garwood interval if ``scipy`` is installed. Otherwise errors will be set to ``sqrt(w2)``.
'''

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
    return parser.parse_args()

def  main():

    # ============== Parse Args =================== #
    args = get_args()
    outdir = args.output + "/red_blue/"
    os.makedirs(outdir, exist_ok = True)

    # ======================================================================
    # TREx Config Reading
    # ======================================================================

    # Get TRExFitter config
    trex_config    = args.trex_cfg
    # Get Regions, Samples and OutputDirectory
    job, objs, attrs   = process_config(trex_config, **{'OutputDir': {},
                                                        'InputFolder': {},
                                                        'Region': {'attrs': ['Label', 'VariableTitle']},
                                                        'Sample': {'attrs': ['Title', 'Type']},
                                                        'Systematic': {'attrs': ['Title'], 'split_by_sc':True},
                                                        })

    histograms_folder = objs['InputFolder'][0]

    trex_outdir = objs['OutputDir'][0]
    if trex_outdir == "":
        logger.error(f"Cannot find an output directory in the config file. This is required.")

    systematics = objs['Systematic']
    syst_attrs = {systematics[idx]: { 'label':    syst_attr['Title']} for idx, syst_attr in enumerate(attrs['Systematic'])}

    regions = objs['Region']
    region_attrs = {regions[idx]: { 'label':    region_attr['Label'], 'varlabel': region_attr['VariableTitle']} for idx, region_attr in enumerate(attrs['Region'])}

    samples = objs['Sample']
    sample_attrs = {samples[idx]: { 'label':    sample_attr['Title'], 'type': sample_attr['Type']} for idx, sample_attr in enumerate(attrs['Sample'])}

    nmax = len(regions)*len(samples)*len(systematics)
    idx = 0
    for region in regions:
        histofile = f"{histograms_folder}/{job}_{region}_histos.root"
        with uproot.open(histofile) as file:
            print(histofile)
            for sample in samples:
                if sample_attrs[sample]['type'].lower() == "ghost":
                    idx += 1
                    continue
                for syst in systematics:
                    idx += 1
#                    if 'ttb_PS' not in syst:    continue
                    if idx % 500 == 0:
                        print(f"INFO:: {idx}/{nmax} done")

                    histopath_syst  = f"{region}/{sample}/{syst}"
                    histopath_nom  = f"{region}/{sample}/nominal"

                    nom_name   = f"{histopath_nom}/{region}_{sample}"
                    up_name    = f"{histopath_syst}/{region}_{sample}_{syst}_Up"
                    down_name  = f"{histopath_syst}/{region}_{sample}_{syst}_Down"
                    try:
                        up_hist         = file[up_name]
                        up_hist_orig    = file[up_name+"_orig"]
                        down_hist       = file[down_name]
                        down_hist_orig  = file[down_name+"_orig"]
                        nominal         = file[nom_name]
                        bin_edges       = nominal.to_numpy()[1]
                        bin_centers     = nominal.axes[0].centers()
                        bin_widths      = (bin_edges[1:] - bin_edges[:-1])
                        up_hist_orig    = rebin(up_hist_orig.to_hist(), up_hist_orig.to_hist().axes[0].name, edges=list(bin_edges))#, reltol=0.1, abstol=None)
                        down_hist_orig  = rebin(down_hist_orig.to_hist(), down_hist_orig.to_hist().axes[0].name, edges=list(bin_edges))#, reltol=0.1, abstol=None)

                    except uproot.exceptions.KeyInFileError:
                        #print(f"Warning:: {histopath} was not found in {histofile}")
                        continue


                    norm_up   = 100*(sum(up_hist.values())   - sum(nominal.values()))/sum(nominal.values())
                    norm_down = 100*(sum(down_hist.values()) - sum(nominal.values()))/sum(nominal.values())

                    def ratio_and_error(variation, nominal, doErr=False):
                        nominal_values = nominal.values()
                        nominal_vars   = nominal.variances()

                        variation_values = variation.values()
                        variation_vars   = variation.variances()

                        difference = variation_values - nominal_values
                        percent_variation =  100*np.divide(difference, nominal_values, out=np.full_like(nominal_values, 0), where=nominal_values>1e-5)

                        err = None
                        if doErr:
                            difference_vars = nominal_vars + variation_vars
                            err_part1 = np.divide(difference_vars, difference**2, out=np.full_like(difference_vars, 0), where=(abs(difference)>1e-5))
                            err_part2 = np.divide(nominal_vars**2 * nominal_values, nominal_values**4, out=np.full_like(difference_vars, 0), where=(nominal_values>1e-4))

                            err = percent_variation*np.sqrt(err_part1 + err_part2, out=np.full_like(difference_vars, 0), where=(abs(difference)>1e-5) & (nominal_values>1e-5))

                        return percent_variation, err

                    up_nom_difference, _    = ratio_and_error(up_hist, nominal, False)
                    up_orig_nom_difference, up_orig_nom_difference_err  = ratio_and_error(up_hist_orig, nominal, True)

                    down_nom_difference, _    = ratio_and_error(down_hist, nominal, False)
                    down_orig_nom_difference, down_orig_nom_difference_err  = ratio_and_error(down_hist_orig, nominal, True)

                    if max(max(up_nom_difference),max(down_nom_difference)) < 1:    continue


                    fig      = plt.figure(figsize=(8,6), dpi=100, layout="constrained")
                    gs       = fig.add_gridspec(ncols=1, nrows=2, height_ratios=(3,1), hspace=0.05)
                    ratio_ax = fig.add_subplot(gs[1, 0])
                    main_ax  = fig.add_subplot(gs[0, 0])
                    plt.setp(main_ax.get_xticklabels(), visible=False)
                    main_ax.set_xlabel('')


                    main_ax_histos = [nominal, up_hist, up_hist_orig, down_hist, down_hist_orig]
                    with warnings.catch_warnings(action="ignore"):
                        mplhep.histplot(main_ax_histos, ax=main_ax, color=['black', 'blue','blue','red','red'], linestyle=['solid', 'solid','dashed','solid','dashed'], linewidth=[1.2,1.2,1.2,1.2,1.2])
                    main_ax.set_xlim(bin_edges[0],bin_edges[-1])
                    main_ax.tick_params(axis="both", which="major", pad=8, labelsize=LABEL_SIZE-4)
                    main_ax.tick_params(direction="in", top=True, right=True, which="both")
                    main_ax.xaxis.set_minor_locator(AutoMinorLocator())
                    main_ax.yaxis.set_minor_locator(AutoMinorLocator())
                    maxes = []
                    for histo in main_ax_histos:
                        max_val = max(histo.values()+np.sqrt(histo.variances(), out=np.full_like(histo.values(), 0.), where=histo.variances()>1e-5))
                        maxes.append(max_val)
                    max_of_maxes = max(maxes)
                    main_ax.set_ylim(0, max_of_maxes*1.5)
                    main_ax.text(0.05, 0.92, texify([region_attrs[region]['label']])[0], horizontalalignment='left', verticalalignment='center', transform = main_ax.transAxes, fontsize=LABEL_SIZE-3)
                    main_ax.text(0.05, 0.86, texify([sample_attrs[sample]['label']])[0], horizontalalignment='left', verticalalignment='center', transform = main_ax.transAxes, fontsize=LABEL_SIZE-3)
                    main_ax.text(0.05, 0.80, texify([syst_attrs[syst]['label']])[0], horizontalalignment='left', verticalalignment='center', transform = main_ax.transAxes, fontsize=LABEL_SIZE-3)
                    red_line   = mlines.Line2D([], [], color='red', marker=None, linestyle='solid', label= fr'up variation ($+1\sigma$, modified) [ {norm_up:.1f}\% ]')
                    red_point  = mlines.Line2D([], [], color='red', marker='+', linestyle='dashed', label= r'up variation ($+1\sigma$, original)')
                    blue_line  = mlines.Line2D([], [], color='blue', marker=None, linestyle='solid', label= fr'down variation ($-1\sigma$, modified) [ {norm_down:.1f}\% ]')
                    blue_point = mlines.Line2D([], [], color='blue', marker='+', linestyle='dashed', label= r'down variation ($-1\sigma$, original)')
                    main_ax.legend(handles=[red_line,red_point,blue_line,blue_point], loc="upper right", frameon=False, ncols=1, fontsize=LABEL_SIZE-4)
                    main_ax.set_ylabel(r'Number of events', rotation='vertical', fontsize=LABEL_SIZE-2, labelpad=10)



                    with warnings.catch_warnings(action="ignore"):
                        mplhep.histplot(up_nom_difference, bins=bin_edges, ax=ratio_ax, color='red', linestyle='solid', linewidth=1.2)
                        mplhep.histplot(down_nom_difference, bins=bin_edges, ax=ratio_ax, color='blue', linestyle='solid', linewidth=1.2)
                    ratio_ax.hlines(y=0, xmin=bin_edges[0],xmax=bin_edges[-1], color='black', linestyle='solid', linewidth=1.2)

                    nominal_error      = np.sqrt(nominal.variances(), out = np.full_like(nominal.variances(), 0), where = (nominal.variances() >=0))
                    nominal_ratio_err =  100*np.divide(nominal_error, nominal.values(), out = np.full_like(nominal.values(), 1e3), where = (nominal.values() > 0))
                    ratio_ax.bar(bin_centers, 2*nominal_ratio_err, width= bin_widths, bottom=(-nominal_ratio_err), fill=False, linewidth=0, edgecolor="gray", hatch=3 * "/",)

                    def plot_orig(difference, difference_err, color):
                        eb = ratio_ax.errorbar(bin_centers, difference, yerr= abs(difference_err), xerr=bin_widths/2, linestyle='none', marker='+', color=color)
                        eb[-1][0].set_linestyle('--')
                        eb[-1][0].set_linewidth(1.2)

                    plot_orig(up_orig_nom_difference, up_orig_nom_difference_err, 'red')
                    plot_orig(down_orig_nom_difference, down_orig_nom_difference_err, 'blue')

                    max_ratio = max([max(abs(up_orig_nom_difference)),
                                    max(abs(down_orig_nom_difference)),
                                    max(abs(up_nom_difference)),
                                    max(abs(down_nom_difference)),
                                    ])

                    ratio_ax.set_ylim(-1*max_ratio*1.15,max_ratio*1.15)
                    ratio_ax.tick_params(axis="both", which="major", pad=8, labelsize=LABEL_SIZE-3)
                    ratio_ax.tick_params(direction="in", top=True, right=True, which="both")
                    ratio_ax.xaxis.set_minor_locator(AutoMinorLocator())
                    ratio_ax.yaxis.set_minor_locator(AutoMinorLocator())
                    ratio_ax.grid(True, which='both',axis='both',color='gray',alpha=0.1)
                    ratio_ax.set_xlim(bin_edges[0],bin_edges[-1])
                    ratio_ax.set_xlabel(texify([region_attrs[region]['varlabel']])[0], fontsize=LABEL_SIZE-2, loc='right', labelpad=10)
                    ratio_ax.set_ylabel(r'percentage variation', rotation='vertical', fontsize=LABEL_SIZE-2, labelpad=2)

                    os.makedirs(f"{outdir}/{region}/{syst}/", exist_ok=True)
                    fig.savefig(f"{outdir}/{region}/{syst}/{sample}.pdf", dpi=300)
                    plt.clf()
                    gc.collect()

if __name__ == "__main__":
    main()
