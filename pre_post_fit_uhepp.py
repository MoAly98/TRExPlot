import uhepp
import argparse
import os
from glob import glob
from collections import defaultdict
import numpy as np

def get_args():
    '''
     Method to retreive arguments from command line
     Returns:
        An `ArgumentParser` object with the parsed values
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', help="TRExFitter Plots folder")
    parser.add_argument('-o','--output',  required = True, help="Folder to dump plots to")

    return parser.parse_args()

def swapPositions(lis, pos1, pos2):
    temp=lis[pos1]
    lis[pos1]=lis[pos2]
    lis[pos2]=temp
    return lis

# ColorMap = {
#   't#bar{t} + #geq1b': '#386641', #'#1e81b0',#'#6666cc',
#   't#bar{t} + #geq1c': '#6A994E', #'#76b5c5',#'#ccccff',
#   't#bar{t} + light':  '#A7C957',    #,'#154c79',#"#cc99ff",
#   'Fakes': "#FF6F59", #"slategray",
#   'single-top': "#BC4749",#"#eab676", #,"darkorange",
#   'Signal': "#F2E8CF", #"#873e23", # "darkred"
#   'others': "slategray", #"darkgreen"
#   #'tWH': "#660000",
# }

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

ColorMap = {
  't#bar{t} + #geq1b': '#d4e09b', #'#1e81b0',#'#6666cc',
  't#bar{t} + #geq1c': '#99AD35', #'#76b5c5',#'#ccccff',
  't#bar{t} + light':  '#606c38',    #,'#154c79',#"#cc99ff",
  'Fakes': "#f19c79", #"slategray",
  'single-top': "#ffdab9",#"#eab676", #,"darkorange",
  'Signal': "chocolate", #"#873e23", # "darkred"
  'others': "slategray", #"darkgreen"
  #'tWH': "#660000",
}

def main():
    args = get_args()
    trexFolder = args.folder
    outFolder  = args.output
    os.makedirs(outFolder, exist_ok=True)
    uheppFiles = glob(f"{trexFolder}/*uhepp.yaml")
    for uheppFile in  uheppFiles:


        name = uheppFile.replace(trexFolder,"").replace(".uhepp.yaml","")
        region = "CR" if "CR_ttb" in name else "SR"
        fit    = "post-fit" if "_postfit" in name else "pre-fit"
        print(name)
        loaded_hist_orig = uhepp.from_yaml(uheppFile)
        loaded_hist = loaded_hist_orig.clone()
        xlabel = loaded_hist.symbol
        if "#" in xlabel and "$" not in xlabel:
            xlabel = rf'${xlabel}$'
        xlabel = rf'{xlabel}'
        loaded_hist.symbol = xlabel.replace("#",'\\')
        simulation = True


        stacks = loaded_hist_orig.stacks
        idx_to_remove = None
        new_stack = None
        norm_signal_stack = None
        signal_norm_factor = 1
        for idx, stack in enumerate(stacks):
            stackItems = stack.content
            if not ((stackItems[0].label == "Total") or (stackItems[0].label == "Data")):
                tot_bkgs_content = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                tot_bkgs_stats   = [0]*len(loaded_hist_orig.yields[stackItems[0].label])

                signal_content = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                signal_stats   = [0]*len(loaded_hist_orig.yields[stackItems[0].label])

                other_bkgs_content = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                other_bkgs_stats   = [0]*len(loaded_hist_orig.yields[stackItems[0].label])
                other_bkgs_label = "other Bkgs"
                new_stack_items = []

                for stackItem in stackItems:
                    ProcLabel = stackItem.label
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
                        new_stack_items.append(stackItem)

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

                other_bkgs_si = uhepp.StackItem(["others"], "other Bkgs")
                new_stack_items.insert(0, other_bkgs_si)

                signalItemIdx = next(idx for idx, item in enumerate(new_stack_items) if item.label == "Signal")
                move_element = new_stack_items.pop(signalItemIdx)
                new_stack_items.append(move_element)
                new_stack = uhepp.Stack(new_stack_items)
                idx_to_remove = idx

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
        #loaded_hist.brand = "ATLAS"
        loaded_hist.subtext = (r'$\sqrt{s}=13\,\text{TeV}\,, \mathcal{L}=140\,\text{fb}^{-1}$'
                              + "\n"
                              + ('Simulation \n' if simulation else '')
                              + f"{region}, {fit}"
                            )


        fig, axes = loaded_hist.render()
        handles, labels = axes[0].get_legend().legend_handles, [text_obj.get_text() for text_obj in fig.axes[0].get_legend().texts]
        axes[0].legend(handles=handles, labels=labels, ncols=2, loc="upper right", fontsize=8, bbox_to_anchor=(0.985, 0.975),frameon=False)

        axes[0].text(0.45, 1.05, rf"$\dagger\,$ Normalised to total background ($\times {signal_norm_factor:.1f}$)", transform=axes[0].transAxes, fontsize=8, verticalalignment='top')

        axes[1].axhline(1, color="black",linestyle='--',linewidth=0.5)
        if not simulation:
            axes[1].set_ylabel("Data/MC", fontsize=10)
        fig.savefig(outname, dpi=300)

if __name__ == "__main__":
    main()