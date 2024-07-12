import uhepp
import argparse
import os
from glob import glob
from collections import defaultdict
import numpy as np
from tools.utils import rgx_match

def get_args():
    '''
     Method to retreive arguments from command line
     Returns:
        An `ArgumentParser` object with the parsed values
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', help="TRExFitter Plots folder")
    parser.add_argument('-o','--output',  required = True, help="Folder to dump plots to")
    parser.add_argument('--filter',         type =str, nargs='+',  help="Filter for the YAML files to process")
    parser.add_argument('--samples',        type =str, nargs='+',  help="Filter samples to plot")
    parser.add_argument('--exclude',        type =str, nargs='+',  help="Exclude lements to plot from a a given yaml")
    parser.add_argument('--noData',         action='store_true',   help="Do not plot data")
    parser.add_argument('--drawTotaledge',         action='store_true',   help="Draw the edge (no errorbar) of total from filtered samples")
    parser.add_argument('--regionName',         type=str,   help="Name of the region being plotted, if only one region is used for all plots.")
    parser.add_argument('--regionLabel',        type=str,   help="Label of the region being plotted, if only one region is used for all plots. Use single quotes in terminal, and entire string in math mode if any part needs math mode.")
    return parser.parse_args()

def swapPositions(lis, pos1, pos2):
    temp=lis[pos1]
    lis[pos1]=lis[pos2]
    lis[pos2]=temp
    return lis

ColorMap = {
  't#bar{t} + #geq1b': '#41764b',#488453',#386641', #'#1e81b0',#'#6666cc',
  't#bar{t} + #geq1c': '#6A994E', #'#76b5c5',#'#ccccff',
  't#bar{t} + light':  '#d4e09b',#A7C957',    #,'#154c79',#"#cc99ff",
  'Fakes': "#FF6F59", #"slategray",
  'single-top': "#BC4749",#"#eab676", #,"darkorange",
  'Signal': "chocolate", #"#873e23", # "darkred"
  'others': "slategray", #"darkgreen"
  #'tWH': "#660000",
}

LabelMap = {}
LabelMap["n_tophad_jets_CBT4"]                        = r"$n_{j,t_{\text{had}}}(\chi^2_{\text{min, ttl}},    \text{PCBT-bin } 4)$"
LabelMap["n_tophad_jets_CBT4_ttAll"]                  = r"$n_{j,t_{\text{had}}}(\chi^2_{\text{min, ttAll}},   \text{PCBT-bin } 4)$"
LabelMap["n_tophad_jets_CBT123_ttAll"]                = r"$n_{j, t_{\text{had}}}(\chi^2_{\text{min, ttAll}}, \text{PCBT-bin }\in\{1,2,3\})$"
LabelMap['leptons_ID_0']                               = r'$id(\ell_0)$'
LabelMap['fwdjets_pt_0']                               = r'$p^T_{\text{fwd},0}$'
LabelMap['sphericity']                                = 'Sphericity'
LabelMap["chi2_min_ttAll"]                            = r"$\chi^2_{\text{min,ttAll}}$"
LabelMap["n_nontophad_jets_CBT123"]                   = r"$n_{j, \text{not-}t_{\text{had}}}(\chi^2_{\text{min, ttl}}, \text{PCBT-bin}\in\{1,2,3\})$"
LabelMap['nonbjets_pt_1']                              = r'$p^T_{\text{light},1}$'
LabelMap['chi2_min']                                  = r'$\chi^2_{\text{top-higgs}}(M_{H}=111.5,M_{t}=168)$'
LabelMap["chi2_min_ttl"]                              = r"$\chi^2_{\text{min,ttl}}$"
LabelMap['nonbjets_tagWeightBin_DL1r_Continuous_0']    = r'$\text{DL1r}_{\text{light},0}$ '
LabelMap['nonbjets_pt_2']                              = r'$p^T_{\text{light},2}$ '
LabelMap['nonbjets_eta_0']                             = r'$\eta_{\text{light},0}$'
LabelMap['tophiggs_chi2_min']                         = r'$\chi^2_{\text{top-higgs}}(M_{H}=113,M_{t}=165)$'
LabelMap['nonbjets_pt_0']                             = r'$p^T_{\text{light},0}$ '
LabelMap['chi2_min_DeltaEta_tH']                      = r'$\Delta\eta_{\text{top-higgs}}$'
LabelMap["n_tophad_jets_CBT0"]                        = r"$n_{j,t_{\text{had}}}(\chi^2_{\text{min, ttl}}, \text{PCBT-bin } 0)$ "
LabelMap['njets']                                     = r'$n_j$'
LabelMap["chi2_min_tophad_m_ttAll"]                   = r"$m_{t_{\text{had}}}(\chi^2_{\text{min, ttAll}})$"
LabelMap['chi2_min_tophad_m']                         = r'$m_{t_{\text{had}}}(\chi^2_{\text{min, ttl}})$'
LabelMap['nonbjets_eta_2']                             = r'$\eta_{\text{light},2}$'
LabelMap['chi2_min_deltaRq1q2']                       = r'$\Delta R(q^{W}_{1},q^{W}_{2})$'
LabelMap["tagnonb_topb_m"]                            = r"$M(b_{\text{top}}, j_{\text{tag}})$"
LabelMap["tagnonb_eta"]                               = r"$\eta_{\text{tag}}$"


def main():
    args = get_args()
    trexFolder = args.folder
    outFolder  = args.output
    yamlFilter = args.filter
    samplesFilter = args.samples
    excludeFilter = args.exclude
    drawTotaledge = args.drawTotaledge
    region_label  = args.regionLabel
    region_name   = args.regionName
    os.makedirs(outFolder, exist_ok=True)
    uheppFiles = glob(f"{trexFolder}/*uhepp.yaml")
    if yamlFilter is not None:
        uheppFiles = [f for f in uheppFiles if any(rgx_match(yF, f.replace(trexFolder+"/","")) for yF in yamlFilter)]

    nFiles = len(uheppFiles)
    assert len(uheppFiles) != 0, "No Uhepp files found.."
    for idx, uheppFile in  enumerate(uheppFiles):
        print(f"INFO:: processing {idx+1}/{nFiles}")
        name = uheppFile.replace(trexFolder+"/","").replace(".uhepp.yaml","")
        if region_name is not  None:
            name = region_name
        region = "CR" if "CR_ttb" in name else "SR"
        if args.regionLabel is not None:
            region = region_label
        fit    = "post-fit" if "_postfit" in name else "pre-fit"
        loaded_hist_orig = uhepp.from_yaml(uheppFile)
        loaded_hist = loaded_hist_orig.clone()

        ## HACK!!
        actual_region =  name.replace("VR_SR_ALL_","").replace("VR_CR_ttb_ALL_","").replace("_postfit","").replace("_prefit","")
        if region_name is not  None:
            actual_region = region_name

        variable = next((label for label in LabelMap.keys() if label == actual_region), None)
        if variable is not None:
            xlabel = LabelMap[variable]
        else:
            xlabel = loaded_hist.symbol
            if "#" in xlabel and "$" not in xlabel:
                xlabel = rf'${xlabel}$'
            xlabel = xlabel.replace("#",'\\')
            if "\\" not in xlabel and "$" in xlabel:  xlabel = xlabel.replace("$","")
            xlabel = rf'{xlabel}'

        loaded_hist.symbol = xlabel
        simulation = True


        stacks = loaded_hist_orig.stacks
        idx_to_remove = None
        new_stack = None
        norm_signal_stack = None
        signal_norm_factor = 1
        new_total_si = None
        for idx, stack in enumerate(stacks):
            stackItems = stack.content
            if not ((stackItems[0].label == "Total") or (stackItems[0].label == "Data")):
                tot_bkgs_content = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                tot_bkgs_stats   = [0]*len(loaded_hist_orig.yields[stackItems[0].label])

                signal_content = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                signal_stats   = [0]*len(loaded_hist_orig.yields[stackItems[0].label])

                other_bkgs_content = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                other_bkgs_stats   = [0]*len(loaded_hist_orig.yields[stackItems[0].label])

                new_total_content = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                new_total_stats   = [1e-6]*len(loaded_hist_orig.yields[stackItems[0].label])

                other_bkgs_label = "other Bkgs"
                new_stack_items = []

                for stackItem in stackItems:
                    ProcLabel = stackItem.label
                    if samplesFilter is not None:
                        if not any(rgx_match(sF, ProcLabel) for sF in samplesFilter):   continue

                    if ProcLabel not in ["Signal", "Total", "Data"]:
                        for bin_index, bin_content in enumerate(loaded_hist_orig.yields[ProcLabel]):
                            tot_bkgs_content[bin_index] += bin_content
                    if ProcLabel == "Signal":
                        for bin_index, bin_content in enumerate(loaded_hist_orig.yields[ProcLabel]):
                            signal_content[bin_index] += bin_content

                    if ProcLabel not in (list(ColorMap.keys()) + ["Total", "Data"]):
                        for bin_index, bin_content in enumerate(loaded_hist_orig.yields[ProcLabel]):
                            other_bkgs_content[bin_index] += bin_content
                    else:
                        #stackItem.edgecolor = "black"
                        new_stack_items.append(stackItem)


                    if samplesFilter is not None and drawTotaledge:
                        for bin_index, bin_content in enumerate(loaded_hist_orig.yields[ProcLabel]):
                            new_total_content[bin_index] += bin_content

                if sum(signal_content) != 0:
                    # Now add a stack for normalised signal
                    tot_bkg         = sum(tot_bkgs_content)
                    signal_yields   = (np.array(signal_content)/sum(signal_content))*tot_bkg
                    signal_norm_factor = tot_bkg/sum(signal_content)
                    signal_yield    = uhepp.Yield(signal_yields, signal_stats)
                    loaded_hist.yields["NormSignal"] = signal_yield
                    norm_signal_si = uhepp.StackItem(["NormSignal"], r"Signal$\,^{\dagger}$", edgecolor=ColorMap['Signal'], linestyle='--')
                    norm_signal_stack = uhepp.Stack([norm_signal_si], bartype="step")


                other_bkgs_yield    = uhepp.Yield(other_bkgs_content, other_bkgs_stats)
                loaded_hist.yields["others"] = other_bkgs_yield

                if any(bin_content !=0 for bin_content in other_bkgs_content):
                    other_bkgs_si = uhepp.StackItem(["others"], "other Bkgs")
                    new_stack_items.insert(0, other_bkgs_si)

                if sum(signal_content) != 0:
                    signalItemIdx = next(idx for idx, item in enumerate(new_stack_items) if item.label == "Signal")
                    move_element = new_stack_items.pop(signalItemIdx)
                    new_stack_items.append(move_element)

                if samplesFilter is not None and drawTotaledge:
                    new_total_yield    = uhepp.Yield(new_total_content, new_total_stats)
                    loaded_hist.yields["trex_total"] = new_total_yield
                    new_total_si = uhepp.StackItem(["trex_total"], "Total")

                new_stack = uhepp.Stack(new_stack_items)
                idx_to_remove = idx

        data_stack_idx = next(idx for idx, stack in enumerate(loaded_hist.stacks) if stack.content[0].label == "Data")
        total_stack_idx = next(idx for idx, stack in enumerate(loaded_hist.stacks) if stack.content[0].label == "Total")
        if args.noData:
            loaded_hist.stacks.pop(data_stack_idx)
            loaded_hist.ratio = None
        if samplesFilter is not None:
            if not drawTotaledge:
                loaded_hist.stacks.pop(total_stack_idx)
            else:
                loaded_hist.stacks[total_stack_idx].error = "no"
                loaded_hist.stacks[total_stack_idx].content[0].edgecolor = "black"

        if idx_to_remove is not None and new_stack is not None:
            loaded_hist.stacks.pop(idx_to_remove)
            loaded_hist.stacks.append(new_stack)

        if norm_signal_stack is not None:
            loaded_hist.stacks.append(norm_signal_stack)


        stacks = loaded_hist.stacks
        for idx, stack in enumerate(stacks):
            for stackItem in stack.content:
                ProcLabel = stackItem.label
                if ProcLabel not in (list(ColorMap.keys()) + ["Total", "Data"]):
                    stackItem.color = ColorMap["others"]
                elif ProcLabel in ColorMap.keys():
                    stackItem.color = ColorMap[ProcLabel]
                if ProcLabel == "Total":
                    stackItem.color = "grey"
                if ProcLabel == "Data":
                    simulation = False
                ProcLabel = ProcLabel.strip('"')
                if "#" in ProcLabel and "$" not in ProcLabel:
                    ProcLabel = f'${ProcLabel}$'
                stackItem.label = ProcLabel.replace("#", "\\")

        outname = f"{outFolder}/{name}.pdf"
        loaded_hist.lumi = 140
        loaded_hist.subtext = (r'$\sqrt{s}=13\,\text{TeV}\,, \mathcal{L}=140\,\text{fb}^{-1}$'
                              + "\n"
                              + ('Simulation \n' if simulation else '')
                              + f"{region}, {fit}"
                            )


        fig, axes = loaded_hist.render()
        ratio_ax = None
        if not isinstance(axes, tuple):
            main_ax = axes
        else:
            main_ax = axes[0]
            ratio_ax = axes[1]

        handles, labels = main_ax.get_legend().legend_handles, [text_obj.get_text() for text_obj in fig.axes[0].get_legend().texts]
        main_ax.legend(handles=handles, labels=labels, ncols=2, loc="upper right", fontsize=8, bbox_to_anchor=(0.985, 0.975),frameon=False)

        if signal_norm_factor != 1:
            main_ax.text(0.45, 1.05, rf"$\dagger\,$ Normalised to total background ($\times {signal_norm_factor:.1f}$)", transform=main_ax.transAxes, fontsize=8, verticalalignment='top')

        if ratio_ax is not None:
            ratio_ax.axhline(1, color="black",linestyle='--',linewidth=0.5)
            if not simulation:
                ratio_ax.set_ylabel("Data/MC", fontsize=10)
            ratio_ax.set_ylim(0.89,1.11)
        fig.savefig(outname, dpi=300)

if __name__ == "__main__":
    main()