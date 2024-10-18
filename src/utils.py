

def tag_substr(dest, tags, max_recursion=20):
    """ 'borrowed' from sotodlib because it's so useful. Do string substitution of all 
    our tags into dest (in-place if dest is a dict). Used to replace tags within yaml files.
    """
    assert(max_recursion > 0)  # Too deep this dictionary.
    if isinstance(dest, str):
        # Keep subbing until it doesn't change any more...
        new = dest.format(**tags)
        while dest != new:
            dest = new
            new = dest.format(**tags)
        return dest
    if isinstance(dest, list):
        return [tag_substr(x,tags) for x in dest]
    if isinstance(dest, tuple):
        return (tag_substr(x,tags) for x in dest)
    if isinstance(dest, dict):
        for k, v in dest.items():
            dest[k] = tag_substr(v,tags, max_recursion-1)
        return dest
    return dest