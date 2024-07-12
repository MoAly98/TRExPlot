import re

def rgx_match(rgx, mystr):
    match = re.match(rgx, mystr)
    return match is not None and match.group(0) == mystr

# Function to extract the substring
def extract_hash_substring(s):
    matches = re.findall(r'#\S*', s)
    matches_with_space = []
    for match in matches:
        if match + " " in s:
            matches_with_space.append("$"+match+"$" + " ")
        else:
            matches_with_space.append("$"+match+"$")

    return matches_with_space if matches else None

def texify(labels_list):
    new_labels = []
    for label in labels_list:
        hash_strings = extract_hash_substring(label)
        new_label = label
        if hash_strings is not None:
            for hash_string in hash_strings:
                better_hash_string = hash_string.replace("#", '\\')
                new_label = new_label.replace(hash_string.replace("$",""), better_hash_string)
        new_labels.append(new_label.replace("%", "\\%"))
    return new_labels