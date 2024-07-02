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

# ColorMap = {
#   't#bar{t} + #geq1b': '#9984D4', #'#1e81b0',#'#6666cc',
#   't#bar{t} + #geq1c': '#CAA8F5', #'#76b5c5',#'#ccccff',
#   't#bar{t} + light':  '#592E83',    #,'#154c79',#"#cc99ff",
#   'Fakes': "#F5811F", #"slategray",
#   'single-top': "#FEB702",#"#eab676", #,"darkorange",
#   'Signal': "#F2E8CF", #"#873e23", # "darkred"
#   'others': "slategray", #"darkgreen"
#   #'tWH': "#660000",
# }

# ColorMap = {
#   't#bar{t} + #geq1b': '#219EBC', #'#1e81b0',#'#6666cc',
#   't#bar{t} + #geq1c': '#7CC3DF', #'#76b5c5',#'#ccccff',
#   't#bar{t} + light':  '#14718C',    #,'#154c79',#"#cc99ff",
#   'Fakes': "#FFB703", #"slategray",
#   'single-top': "#FC9D02",#"#eab676", #,"darkorange",
#   'Signal': "#F98702", #"#873e23", # "darkred"
#   'others': "slategray", #"darkgreen"
#   #'tWH': "#660000",
# }

# ColorMap = {
#   't#bar{t} + #geq1b': '#005F60', #'#1e81b0',#'#6666cc',
#   't#bar{t} + #geq1c': '#008285',#008083', #'#76b5c5',#'#ccccff',
#   't#bar{t} + light':  '#00D9DD',#249EA0',    #,'#154c79',#"#cc99ff",
#   'Fakes': "#FAAB36", #"slategray",
#   'single-top': "#F78104",#"#eab676", #,"darkorange",
#   'Signal': "chocolate", #"#873e23", # "darkred"
#   'others': "slategray", #"darkgreen"
#   #'tWH': "#660000",
# }

def main():
    args = get_args()
    trexFolder = args.folder
    outFolder  = args.output
    yamlFilter = args.filter
    samplesFilter = args.samples
    excludeFilter = args.exclude
    drawTotaledge = args.drawTotaledge
    os.makedirs(outFolder, exist_ok=True)
    uheppFiles = glob(f"{trexFolder}/*uhepp.yaml")
    if yamlFilter is not None:
        uheppFiles = [f for f in uheppFiles if any(rgx_match(yF, f) for yF in yamlFilter)]

    assert len(uheppFiles) != 0, "No Uhepp files found.."
    for uheppFile in  uheppFiles:
        name = uheppFile.replace(trexFolder,"").replace(".uhepp.yaml","")

        region = "CR" if "CR_ttb" in name else "SR"
        fit    = "post-fit" if "_postfit" in name else "pre-fit"
        loaded_hist_orig = uhepp.from_yaml(uheppFile)
        loaded_hist = loaded_hist_orig.clone()
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
        fig.savefig(outname, dpi=300)

if __name__ == "__main__":
    main()