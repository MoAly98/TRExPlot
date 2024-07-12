'''
This script is responsible for parsing the configuration file
from trex-fitter, to return any groups of blocks. The hard-coded
blocks retrieved are Samples, Systematics and Regions
'''
# Python imports
import re
from collections import defaultdict
from tools.utils import rgx_match

# def remove_after_percent(s):
#     # Regular expression to match % inside quotes
#     pattern = re.compile(r'"[^"]*%"')

#     # Find all matches of % inside quotes and store them
#     quoted_percents = pattern.findall(s)

#     # Temporarily replace % inside quotes with a placeholder
#     temp_str = s
#     for i, match in enumerate(quoted_percents):
#         temp_str = temp_str.replace(match, f"__PLACEHOLDER_{i}__")

#     # Find the first unquoted %
#     percent_pos = temp_str.find('%')

#     # Remove everything after the first unquoted % if it exists
#     if percent_pos != -1:
#         temp_str = temp_str[:percent_pos]

#     # Restore the quoted % placeholders back to their original form
#     for i, match in enumerate(quoted_percents):
#         temp_str = temp_str.replace(f"__PLACEHOLDER_{i}__", match)

#     return temp_str

def remove_after_percent(s):
    # Split the string into parts that are either quoted segments or others
    parts = re.split(r'("[^"]*")', s)

    result = []
    for part in parts:
        if part.startswith('"') and part.endswith('"'):
            # If the part is quoted, keep it as is
            result.append(part)
        else:
            # If the part is not quoted, remove everything after the first %
            percent_pos = part.find('%')
            if percent_pos != -1:
                part = part[:percent_pos]
            result.append(part)

    return ''.join(result)


def open_config(cfgpath):
    '''
    Method to open the configuration file and read all lies
    Args:
        cfgpath (str): The path to the config file
    Returns:
        list of lines
    '''
    with open(cfgpath, 'r') as f:
        clines = f.readlines()
    return clines

def get_blocks_of(block, line, clines, *, attrs = None, split_by_sc = False):
    '''
    Method to look for the start of a block in a given line and
    return the name of the object found on that line as well as requested
    object attributes within the block. This assumes the configuration
    file is ordered, such that the definitions of Samples/Regions/etc..
    come in a row, not scattered across the config

    Block: A group of lines defining 1 Sample/Systematic/Region
    Object:  1x Sample/Systematic/Region
    Attribute: An option defined for an object, such as Type of sample

    Args:
        block (str): The block being retrieved (e.g. Sample, Systematic,..)
        line (str): The candidate first line of a block
        clines List(str): List of lines in the config file
        attrs (Optional[List[str]]): List of attributes to retreive from one block
        split_by_sc (bool): Decide if the block with multiple object definitions should be
                            split using the ";" delimeter
    Return:
        if not split_by_sc:
            object name and a dictionary of attributes
        else:
            list of object names and a list of dictionaries of attributes
    '''

    # Check if there is a block starting the definition of an object on the candidate line
    if re.search(f'^{block}:', line.strip()) is not None:

        # Save the object name found on the line
        objname = line.strip().replace(f"{block}:","").replace('"','').replace("'","").strip().rstrip()
        # If we are extractnig attributes
        if attrs is not None:
            # Turn the given attributes to a list, if they weren't already
            attrs = attrs if isinstance(attrs, list) else [attrs]

            # =============== Identify the end of the current block ====================
            # Find start of next block of same object type (e.g. next Sample)
            next_obj = next((l for l in clines[clines.index(line)+1:] if re.search(f'^{block}:', l.strip()) is not None), None)
            # Get index of the line where next block starts, unless there are no more, return None
            next_obj_idx = clines.index(next_obj) if next_obj is not None else None
            # Look for attributes up to next block start, or end of file if no new blocks found
            query_to = next_obj_idx if next_obj_idx is not None else -1

            # =============== Save object attributes to a dict ====================
            attr_data = {}
            for attr in attrs:
                # If the attribute is not found, an empty string
                attr_val = next((l for l in clines[clines.index(line)+1:query_to]
                                 if re.search(f'^{attr}:', l.strip()) is not None), '')

                #attr_val = attr_val.split('%', 1)[0] # Get rid of commented part, can be whole line.
                attr_val = remove_after_percent(attr_val)
                attr_data[attr] =   attr_val.replace(f"{attr}:","").replace('"','').replace("'","").strip()

            # If user asks to split block by ";"
            if split_by_sc:
                # Object name becomes a list of object names
                objname = objname.split(';')
                # Assume all objects have same attributes,
                # create a list of size = sub-objects with
                # elements being the same attributes dict
                new_attrs = []
                for idx, subobj in enumerate(objname):
                    new_data = {}
                    for attr, vals in attr_data.items():
                        attr_split = vals.split(';')
                        if len(attr_split) == len(objname):
                            attr_val = attr_split[idx]
                        elif len(attr_split) == 1:
                            attr_val = attr_split[0]
                        else:
                            print("ERROR:: Could not split by ;")
                            exit(1)
                        new_data[attr] = attr_val

                    new_attrs.append(new_data)

                return objname, new_attrs

            # If not splitting by ";", we are done
            else:
                return objname, attr_data

        # If no attributes to be retrieved
        else:
            if split_by_sc:
                objname = objname.split(';')
                attrs = [None]*len(objname)
            else:
                objname = objname
                attrs = None

            return objname, attrs

    # If no block-start for an object of the type
    # specified is found, return Nones.
    else:
        return None, None


def grab_from_config(clines, cfgpath, **kwargs):
    '''
    Method to decide what to grab from the config, and call get_blocks_of
    to grab the objects, along with attributes needed.

    TODO:: Implement a dictionary of Object type to attributes and split_by_sc
           to allow user to pick out anything from the config, rather than the
           hard-coded stuff

    Args:
        clines (List[str]): List of config file lines

    Returns:
        list of samples, list of sample attributes dicts
        list of regions, list of region attribute dicts
        list of systematics, list of systematic attributes dicts

    '''
    obj, obj_attrs = defaultdict(list), defaultdict(list)
    alllines = clines
    for i, line in enumerate(clines):
        #line = line.split('%', 1)[0] # Get rid of commented part, can be whole line.
        line = remove_after_percent(line)
        if "INCLUDE" in line:
            subconfig_relpath  = line.replace('INCLUDE:','').strip()
            if subconfig_relpath[0] == '/':
                subconfig_fullpath = subconfig_relpath
            else:
                subconfig_fullpath = "/".join(cfgpath.split('/')[:-1])+'/'+subconfig_relpath

            morelines = open_config(subconfig_fullpath)
            alllines += morelines

    for i, aline in enumerate(alllines):
        #line = aline.split('%', 1)[0] # Get rid of commented part, can be whole line.
        line = remove_after_percent(aline)
        for objflag, options in kwargs.items():
            split_by_sc = options.get('split_by_sc',False)
            attrs = options.get('attrs', None)
            objname, attrs = get_blocks_of(objflag, line, clines, attrs = attrs, split_by_sc = split_by_sc)
            if objname is not None:
                if not split_by_sc:
                    obj[objflag].append(objname)
                    obj_attrs[objflag].append(attrs)
                else:
                    obj[objflag].extend(objname)
                    obj_attrs[objflag].extend(attrs)

    return obj, obj_attrs

def process_config(cfgpath, **kwargs):
    '''
    Steering method for the config reader, which is called from other scripts
    Args:
        cfgpath (str): Path to configuration file

    Returns:
        job name,
        list of samples, list of sample attributes dicts
        list of regions, list of region attribute dicts
        list of systematics, list of systematic attributes dicts
    '''
    # Get the config lines
    clines = open_config(cfgpath)
    # Get the objects
    #samples, samp_attrs, regions, region_attrs, systematics, systematic_attrs = grab_from_config(clines, **kwargs)
    objs, objs_attrs = grab_from_config(clines, cfgpath, **kwargs)
    # Workout the job name
    job = [get_blocks_of("Job", line, clines)[0]  for line in clines if get_blocks_of("Job", line, clines)[0] is not None ]
    if job == []:
        job = [get_blocks_of("MultiFit", line, clines)[0]  for line in clines if get_blocks_of("MultiFit", line, clines)[0] is not None ]
    if job == []:
        print("ERROR:: Could not extract Job/Multifit name")
        exit(1)
    job = job[0]
    return job, objs, objs_attrs

def expected_histos(storage, reqregions = None, reqsamples = None, excluded = None):
    '''
    Method to workout which histograms are expected from the config
    Args:
        storage (Storage):   A `Storage` object, holding the information extracted
                            from the trex-fitter config and master root-files contents
    Returns:
        A list of maps specifying all expected region, sample, syst, template and histogram name combos
    '''
    in_reqregions = lambda x: (reqregions is None) or (any(rgx_match(reqreg, x) for reqreg in reqregions))
    in_reqsamples = lambda x: (reqsamples is None) or (any(rgx_match(reqsamp, x) for reqsamp in reqsamples))
    in_excluded   = lambda x: (excluded   is not None) and (any(rgx_match(excl, x) for excl in excluded))
    systs_maps = []
    nominals_maps = []
    for region in storage.regions:
        if not in_reqregions(region): continue
        if in_excluded(region): continue
        # Zip samples and their attributes to simplify loop
        samp_use_syst_zip = zip(storage.samples, storage.samples_attrs)
        for samp_item in samp_use_syst_zip:

            # Process sample name and attributes
            sample, sample_attributes = samp_item
            if not in_reqsamples(sample): continue
            if in_excluded(sample): continue

            use_syst, sample_type = sample_attributes['UseSystematics'], sample_attributes['Type']

            # Ghost samples don't get histograms
            if sample_type.lower() == 'ghost':  continue
            # Data samples only have nominal histograms
            if sample_type.lower() == 'data':  use_syst = 'false'

            # We are ready to define a nominal
            nominal = f'{region}_{sample}'
            nominals_maps.append(dict(region=region, sample=sample, hist = nominal))
            # If the sample uses systematics (either explicitly
            # or by not specifying UseSystematics option), loop over systs
            if use_syst.lower().strip() == 'true' or use_syst.lower().strip() == '':

                # Zip systematics and their attributes to simplify loop
                syst_exclude_samples_zip = zip(storage.systematics, storage.systematics_attrs)
                for sys_item in syst_exclude_samples_zip:

                    # Process sample name and attributes
                    systematic, syst_attributes = sys_item
                    exclude, keep_samples = syst_attributes['Exclude'], syst_attributes['Samples']
                    syst_symm = syst_attributes['Symmetrisation']
                    # Skip samples/regions excluded from this systematic
                    if any([rgx_match(excluded.replace('*','.*'), sample) and excluded != '' for excluded in exclude.split(',')]):  continue
                    if any([rgx_match(excluded.replace('*','.*'), region)  and excluded != '' for excluded in exclude.split(',')]):  continue
                     # if samples are specified for a systematic, then skip samples not in specified
                    if keep_samples.split(',') != ['']:
                        if all([not rgx_match(ks.replace('*','.*'), sample) for ks in keep_samples.split(',')]): continue


                    # Need (at least) 1x histogram for the up and 1x for down variations
                    for template in ['Up','Down']:
                        expected_hist = f'{region}_{sample}_{systematic}_{template}'
                        systs_maps.append(dict(region=region, sample=sample, systematic=systematic, template=template, hist=expected_hist, syst_symm=syst_symm))

    return nominals_maps+systs_maps
