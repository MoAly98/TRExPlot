import os
import numpy
import yaml
def read_txtfile(file):
    # Hack when not all bootstraps are done
    if not os.path.exists(file):
        return []

    with open(file, 'r') as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines if l.strip() != '' and 'NUISANCE_PARAMETERS' not in l]
        corrmatrix_start = lines.index('CORRELATION_MATRIX')
        nps = lines[:corrmatrix_start]
        corr = lines[corrmatrix_start+2:-2]
        corrmat = []
        for line in corr:
            vals = line.split()
            corrmat.append(vals)

    return nps, corrmat

def read_ranking(file):
    # Hack when not all bootstraps are done
    assert os.path.exists(file), "Ranking file is not found"

    with open(file, 'r') as f:
        ranking = yaml.safe_load(f)

    return map_ranking(ranking)

def map_ranking(ranking_list):
    ranks = {}
    for param in ranking_list:
        # Name
        name         = param['Name']
        # The pull
        bestfit_pull  = param['NPhat']
        bestfit_errUp = param['NPerrHi']
        bestfit_errLow = param['NPerrLo']
        # The impact
        impactUp_postfit = param['POIup']
        impactDo_postfit = param['POIdown']

        impactUp_prefit = param['POIupPreFit']
        impactDo_prefit = param['POIdownPreFit']

        ranks[name] = {'pull': bestfit_pull, 'unc_hi': bestfit_errUp, 'unc_lo': bestfit_errLow,
                        'impact_hi_pre': impactUp_prefit, 'impact_lo_pre': impactDo_prefit,
                        'impact_hi_post': impactUp_postfit, 'impact_lo_post': impactDo_postfit,
                        }

    return ranks

def map_pulls(lines):
    pulls = {}
    for l in lines:
        nuisp, cen, left, right = l.split()
        pulls[nuisp] = {'pull': float(cen), 'constr': (abs(float(left))+abs(float(right)))/2, 'unc_hi': float(left), 'unc_lo': float(right)}
    return pulls

def read_pulls(file):
    return map_pulls(read_txtfile(file)[0])

def read_corr(file):
    nps, corr = read_txtfile(file)
    np_ver = numpy.array(corr, dtype='float64')

    return np_ver

def read_CovErrDecomp_txtfile(file):

    assert os.path.exists(file), "Error decomposition  by covariance file  is not found"
    with open(file, 'r') as f:
        lines = f.readlines()
        totals = lines[:7]
        nps = lines[7:]

    return map_CovErrDecomp(totals, nps)

def map_CovErrDecomp(totals, nps):
    total_impacts, component_impacts = {}, {}
    for total_line in totals:
        nuisp, impact_symm, impact_up, impact_do = total_line.split()
        # Hack for TREx bug swapping total syst and stat impacts
        if "SYST" in nuisp: nuisp = nuisp.replace("SYST", "STAT")
        elif "STAT" in nuisp and "MCSTAT" not in nuisp: nuisp = nuisp.replace("STAT", "SYST")
        total_impacts[nuisp] = {'impact_symm': float(impact_symm), 'impact_up': float(impact_up), 'impact_do': float(impact_do)}
    for np_line in nps:
        if np_line in ["\n", ""]:   continue
        nuisp, impact_symm, impact_up, impact_do = np_line.split()
        nuisp  = nuisp.replace("alpha_","")
        component_impacts[nuisp] =  {'impact_symm': float(impact_symm), 'impact_up': float(impact_up), 'impact_do': float(impact_do)}
    return total_impacts, component_impacts

def read_CovErrDecomp(file):
    totals, nps = read_CovErrDecomp_txtfile(file)
    return totals, nps

if __name__ == "__main__":
   read_CovErrDecomp("/eos/user/m/maly/thbb/thbb-fit/tHbb_v34_v3_FullSyst/AsimovResults/split_4v5FS_by_HF_no_ttb_split/tHbb_v34_v3_FullSyst_Rect_AltBDT_Asimov_CRSR_SPlusB_with4v5FS_decorrHF_withPtHard/tHbb_v34_v3_FullSyst_Rect_AltBDT/Fits/tHbb_v34_v3_FullSyst_Rect_AltBDT_errDecomp_tH_tWH_NORM.txt")