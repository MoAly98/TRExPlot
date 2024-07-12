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
    if not os.path.exists(file):
        return []

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
        pulls[nuisp] = {'pull': float(cen), 'constr': (abs(float(left))+abs(float(right)))/2}
    return pulls

def read_pulls(file):
    return map_pulls(read_txtfile(file)[0])

def read_corr(file):
    nps, corr = read_txtfile(file)
    np_ver = numpy.array(corr, dtype='float64')

    return np_ver