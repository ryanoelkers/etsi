""" This serves as the configuration file for the ETSI pipeline. """


class Configuration:

    # Computer for reduction
    MACHINE = 'thoroslab'
    FILE_EXTENSION = '.fits'

    # update for different data products
    STAR = 'WASP3b'
    RA = 278.6333
    DEC = 35.6617
    TC = 2454143.851120
    PERIOD = 1.846835100
    PHOTOMETRY = 'APER'
    APERTURE_SHAPE = 'ellipse'

    # steps to skip
    CLEAN_SKIP = 'Y'
    WRITE_SKY = 'N'
    CALIBRATE_SKIP = 'Y'
    PHOTOMETRY_SKIP = 'Y'
    LIGHTCURVE_SKIP = 'Y'
    COLOR_SKIP = 'N'

    # get image information
    AXS_X = 1676
    AXS_Y = 1266
    AXS = 2048
    GAIN = 0.38  # in e-/ADU
    FOV = 5  # size of the image in arc minutes
    SEARCH_DIST = FOV / 60.0

    # rough size of PSF for star finding and clip functions, also avoidance for edge
    SIGMA_PSF = 7
    PSF_X = 41
    PSF_Y = 21
    PSF_EDGE_LIMIT = 100
    PSF_BUFFER = 0

    # update sky subtraction specific information
    PIX_BOX = 128
    PIX = 32

    # update the image axes to work for the a PIX_BOX setting
    X_CUT = int(AXS_X % PIX_BOX)  # work to make it divisible by a PIX x PIX box
    Y_CUT = int(AXS_Y % PIX_BOX)
    X_CENT = int(AXS_X / 2)  # get the center with respect to the old image
    Y_CENT = int(AXS_Y / 2)
    AXIS_X = int(AXS_X - X_CUT) # get the size of the new image
    AXIS_Y = int(AXS_Y - Y_CUT)

    # if the image is divisible by the box, then don't worry
    if (Y_CUT != 0) | (X_CUT != 0):
        CUT_IMAGE = 'Y'
    else:
        CUT_IMAGE = 'N'

    # a photometry configuration
    FWHM = 15.  # fwhm of the image
    THRESHOLD = 5.  # the threshold for a source above the background

    # aperture information
    CIRC_APER_SIZE = 32  # circular aperture
    ELIP_APER_A = 20
    ELIP_APER_B = 10
    RECT_APER_H = 20  # rectangular aperture y size
    RECT_APER_W = 40  # rectangular aperture x size

    # aperture annulus for the sky background automatically determined from the main aperture
    CIRC_ANNULI_INNER = CIRC_APER_SIZE + 2
    CIRC_ANNULI_OUTER = CIRC_APER_SIZE + 4

    RECT_ANNULI_H0 = RECT_APER_H + 2
    RECT_ANNULI_HN = RECT_APER_H + 4
    RECT_ANNULI_W0 = RECT_APER_W + 2
    RECT_ANNULI_WN = RECT_APER_W + 4

    ELIP_ANNULI_A0 = ELIP_APER_A + 2
    ELIP_ANNULI_AN = ELIP_APER_A + 4
    ELIP_ANNULI_B0 = ELIP_APER_B + 2
    ELIP_ANNULI_BN = ELIP_APER_B + 4

    # ETSI PSF separations
    FOOT_PRINT = 40
    NUM_PSF = 4
    ST1 = 0.
    ST2 = ST1 + 50.
    ST3 = ST1 + 120.
    ST4 = ST1 + 200.

    # output paths for logging, temporary files, figures etc
    WORKING_DIRECTORY = "C:\\Users\\astrolab\\Development\\etsi\\"
    ANALYSIS_DIRECTORY = WORKING_DIRECTORY + 'analysis\\'
    LOG_DIRECTORY = WORKING_DIRECTORY + 'logs\\'
    QUERIES_DIRECTORY = WORKING_DIRECTORY + 'queries\\'

    # input paths for data etc
    DATA_DIRECTORY = "F:\\WASP3b\\data\\"
    DARKS_DIRECTORY = DATA_DIRECTORY + "darks\\"
    BIAS_DIRECTORY = DATA_DIRECTORY + "bias\\"
    FLATS_DIRECTORY = DATA_DIRECTORY + "flats\\"
    RAW_DIRECTORY = DATA_DIRECTORY + "raw\\"
    CLEAN_DIRECTORY = DATA_DIRECTORY + "clean\\"
    MASTER_DIRECTORY = DATA_DIRECTORY + "master\\"
    CENTROID_DIRECTORY = MASTER_DIRECTORY + "centroids\\"
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + "lc\\"

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY,
                   QUERIES_DIRECTORY, CLEAN_DIRECTORY, MASTER_DIRECTORY, LIGHTCURVE_DIRECTORY,
                   CENTROID_DIRECTORY, RAW_DIRECTORY, FLATS_DIRECTORY, BIAS_DIRECTORY, DARKS_DIRECTORY]
