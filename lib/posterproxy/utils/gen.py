"""General utility functions."""


def set_default_kw(kw_in, kw_default={}):
    """Return new dictionary matching original dictionary with any missing keys
    set to default values. If kw_in is None, just return kw_default.
    """

    if kw_in is None:
        return kw_default
    else:
        kw_out = kw_in.copy()
        for key in kw_default.keys():
            if key not in kw_in.keys():
                kw_out[key] = kw_default[key]
        return kw_out

