import argparse
import os, re
from glob import glob
from collections import defaultdict
import yaml
from pprint import pprint
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
    parser.add_argument('-b', '--bins',   required = False, default=None, help="Which bins to plot composition for explicitly")
    parser.add_argument('--oneplot',        action = 'store_true', help="all pie charts on one plot? default is split by region+pre/postfit combos")
    parser.add_argument('--ncols',          default =3, type = int, help="number of pie charts per row. Default is 3")
    parser.add_argument('--preOnly',        action = 'store_true',  help="Prefit only?")
    parser.add_argument('--postOnly',       action = 'store_true',  help="Postfit only?")
    parser.add_argument('--noFitInTitle',   action = 'store_true',  help="Don't write pre/post fit on subplot titles")
    parser.add_argument('--filter',         type =str, nargs='+',  help="Filter for the YAML files to process")
    parser.add_argument('--samples',        type =str, nargs='+',  help="Filter samples to plot")
    return parser.parse_args()

def make_fig(yields_mapping, metadata, oneplot, max_columns):
    # Compute background fractions of total background
    n_plots = 0
    if oneplot:
        for yaml in yields_mapping:
            n_plots += len(list(yields_mapping[yaml].keys()))
    else:
        n_plots = len(yields_mapping.keys())

    rows = (n_plots + max_columns - 1) // max_columns
    columns = min(n_plots, max_columns)
    specs = []
    for i in range(rows):
        row_specs = [{'type': 'domain'} for _ in range(columns)]
        specs.append(row_specs)

    titles = []
    for element in metadata:
        if oneplot:
            for subplot in metadata[element]:
                titles.append(metadata[element][subplot]['title'])
        else:
            titles.append(metadata[element]['title'])

    fig = make_subplots(rows=rows, cols=columns,specs=specs, vertical_spacing=0.2, horizontal_spacing=0.01, subplot_titles=titles)
    return fig

# ColorMap = {
#   't#bar{t} + #geq1b': '#d4e09b', #'#1e81b0',#'#6666cc',
#   't#bar{t} + #geq1c': '#99AD35', #'#76b5c5',#'#ccccff',
#   't#bar{t} + light':  '#606c38',    #,'#154c79',#"#cc99ff",
#   'Fakes': "#f19c79", #"slategray",
#   'single-top': "#ffdab9",#"#eab676", #,"darkorange",
#   'Signal': "chocolate", #"#873e23", # "darkred"
#   'others': "slategray", #"darkgreen"
#   #'tWH': "#660000",
# }

ColorMap = {
  't#bar{t} + #geq1b': '#41764b',#488453',#386641', #'#1e81b0',#'#6666cc',
  't#bar{t} + #geq1c': '#6A994E', #'#76b5c5',#'#ccccff',
  't#bar{t} + light':  '#d4e09b',#A7C957',    #,'#154c79',#"#cc99ff",
  'Fakes': "#FF6F59", #"slategray",
  'single-top': "#BC4749",#"#eab676", #,"darkorange",
  'Signal': "chocolate", #"#873e23", # "darkred"
  'others': "slategray", #"darkgreen"
  't#bar{t} + 1b':  '#59a5d8',
  't#bar{t} + 2b':  '#386fa4',
  't#bar{t} + 1B':  '#133c55',

  #'tWH': "#660000",
}

def main():
    args = get_args()
    # Input folder and files
    trexFolder = args.folder
    regex = "*prefit.yaml" if args.preOnly else  "*postfit.yaml" if args.postOnly else "*.yaml"
    allFiles = glob(f"{trexFolder}/{regex}")
    samplesFilter = args.samples
    yamlFilter = args.filter
    notUheppFiles = [f for f in allFiles if "uhepp" not in f]
    if yamlFilter is not None:
        notUheppFiles = [f for f in notUheppFiles if any(rgx_match(yF, f) for yF in yamlFilter)]
    assert len(notUheppFiles) != 0, "ERROR:: No files found to read from"
    ordered_list = []
    for path in notUheppFiles:
        fname = path.replace(trexFolder,"")
        if "SR" in fname:
            no_prepost = fname.replace("_prefit.yaml","").replace("_postfit.yaml","")
            if args.preOnly:
                ordered_list.append(no_prepost+"_prefit.yaml")
            elif args.postOnly:
                ordered_list.append(no_prepost+"_postfit.yaml")
            else:
                ordered_list.extend([no_prepost+"_prefit.yaml", no_prepost+"_postfit.yaml"])
            break

    for path in notUheppFiles:
        fname = path.replace(trexFolder,"")
        if "CR_ttb" in fname:
            no_prepost = fname.replace("_prefit.yaml","").replace("_postfit.yaml","")
            if args.preOnly:
                ordered_list.append(no_prepost+"_prefit.yaml")
            elif args.postOnly:
                ordered_list.append(no_prepost+"_postfit.yaml")
            else:
                ordered_list.extend([no_prepost+"_prefit.yaml", no_prepost+"_postfit.yaml"])
            break

    notUheppFiles = [f"{trexFolder}/{ol}" for ol in ordered_list]

   # Output folder
    outFolder  = args.output
    os.makedirs(outFolder+"/BkgComp/", exist_ok=True)
    # Bins of interest
    bins_to_split = args.bins
    # Number of columns
    max_columns = args.ncols

    # ============================================== #
    # ========== Aggregate information ============ #
    # ============================================== #
    allRegions_yields_fracs_metadata = {}
    allRegions_yields_fracs = {}
    for yamlFile in  notUheppFiles:
        fname = yamlFile.replace(trexFolder,"").replace(".yaml","")
        region = "CR" if "CR_ttb" in fname else "SR"
        fit    = "post-fit" if "_postfit" in fname else "pre-fit"

        with open(yamlFile,'r') as f:
            yields = yaml.safe_load(f)


        per_bin_yields = defaultdict(lambda: defaultdict(float))
        for sample_data in yields['Samples']:
            name = sample_data['Name']
            if samplesFilter is  not None:
                if not any(rgx_match(sF, name) for sF in samplesFilter):    continue
            sample_yield = sample_data['Yield']
            for bin_idx, bin_yield in enumerate(sample_yield):
                per_bin_yields[bin_idx][name]       =  bin_yield
                if name != "Signal":
                    per_bin_yields[bin_idx]['totalBkg'] += bin_yield
            per_bin_yields['allBins'][name]       = sum(sample_yield)
            if name != "Signal":
                per_bin_yields['allBins']['totalBkg'] += sum(sample_yield)

        # ============================================== #
        # ========== Filter Bin information ============ #
        # ============================================== #
        filtered_per_bin_yields = defaultdict(lambda: defaultdict(float))
        filtered_per_bin_yields_fracs =  defaultdict(lambda: defaultdict(float))
        metadata = {}
        for bin_idx, bin_data in per_bin_yields.items():
            if bins_to_split is None and bin_idx != "allBins":   continue
            if bins_to_split is not None and bin_idx != "allBins" and bins_to_split.lower() != 'all':
                if not any(bin_idx == bin_to_split for bin_to_split in bins_to_split):    continue
            filtered_per_bin_yields[bin_idx] = bin_data


            for name, data in bin_data.items():

                if name == "Signal":    continue
                if name not in ColorMap.keys():
                    key = f"others"
                else:
                    key = f"{name}"
                if name == 'totalBkg':  continue

                filtered_per_bin_yields_fracs[bin_idx][key] += filtered_per_bin_yields[bin_idx][name]/filtered_per_bin_yields[bin_idx]['totalBkg']

            labels, colors = [], []
            for key in filtered_per_bin_yields_fracs[bin_idx]:
                label = key.replace('#','\\')
                label = fr'${label}$' if "\\" in label else label
                if label == "others":   label = "other Bkgs"
                labels.append(label)


                color = ColorMap[key]
                colors.append(color)


             # If per bin comp used, each subplot has title REGION + PRE/POSTFIT
            each_pie_title             = f"Bin {bin_idx}" if bin_idx != "allBins" else region
            if not args.noFitInTitle:
                each_pie_title += f", {fit}"

            metadata[bin_idx] = {'labels': labels, 'title': each_pie_title, 'colors': colors}

        # ================= Save for each region =============== #
        allRegions_yields_fracs_metadata[fname] = metadata
        allRegions_yields_fracs[fname] = filtered_per_bin_yields_fracs

    if args.oneplot:    fig = make_fig(allRegions_yields_fracs, allRegions_yields_fracs_metadata, args.oneplot, max_columns)

    layout = {
            "margin":   dict(t=40, b=40, l=60, r=60),
            "height":   400,
            "width":    750,
            "font":     dict(family="Serif"),  # Set the font family here
            "legend":   dict(font=dict(family="Serif"), orientation='h')
        }

    for idx_oneplot, (fname, filtered_per_bin_yields_fracs) in enumerate(allRegions_yields_fracs.items()):
        if not args.oneplot:    fig = make_fig(filtered_per_bin_yields_fracs, allRegions_yields_fracs_metadata[fname], args.oneplot, max_columns)
        for idx_manyplots, (bin_idx, bin_data) in enumerate(filtered_per_bin_yields_fracs.items()):
            labels   = allRegions_yields_fracs_metadata[fname][bin_idx]['labels']
            colors   = allRegions_yields_fracs_metadata[fname][bin_idx]['colors']
            title    = allRegions_yields_fracs_metadata[fname][bin_idx]['title']

            idx_of_plot = idx_oneplot if args.oneplot else idx_manyplots
            row = idx_of_plot // max_columns + 1
            col = idx_of_plot % max_columns + 1
            pie = go.Pie(labels=labels, values=list(filtered_per_bin_yields_fracs[bin_idx].values()), title="", name=bin_idx, marker=dict(colors=colors))
            fig.add_trace(pie,row, col)


        if not args.oneplot:
            fig.update_layout(layout)
            fig.write_image(f"{outFolder}/BkgComp/{fname}.png", scale=4)

    if args.oneplot:
        fig.update_layout(layout)
        fig.write_image(f"{outFolder}/BkgComp/summary.png", scale=4)

if __name__ == "__main__":
    main()