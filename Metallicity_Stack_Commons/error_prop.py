from chun_codes import random_pdf, compute_onesig_pdf

def construct_pdf(values, RMS, seed_i=1, n_iter=1000):
    '''
    Constructs probability distribution function (PDF) based on input
    values and their associated uncertainty

    :param values: list or numpy array of values/parameters
    :param RMS: 1-sigma errors associated with values (same dimension)
    :param seed_i: integer value for initial seed for np.random. Default: 1
    :param n_iter: Number of iterations. Default: 1000

    :return pdf_arr: numpy array of size (size of values, n_iter)
    '''

    pdf_arr = random_pdf(values, RMS, seed_i=seed_i, n_iter=n_iter, silent=False)

    return pdf_arr
