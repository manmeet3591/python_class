# https://github.com/phydrus/pyet/blob/master/pyet/radiation.pyFrom 
def calc_lambda(tmean):
    """ Latent Heat of Vaporization [MJ kg-1].
    Parameters
    ----------
    tmean: pandas.Series/float, optional
        average day temperature [°C]
    Returns
    -------
    pandas.Series containing the calculated Latent Heat of Vaporization
        [MJ kg-1].
    Examples
    --------
    >>> lambd = calc_lambda(tmean)
    Notes
    -----
    Based on equation (3-1) in [allen_1998]_.
    """
    return 2.501 - 0.002361 * tmean
  
