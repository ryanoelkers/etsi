""" This serves as the configuration file for the ETSI pipeline. """


class Configuration:

    # Computer for reduction
    MACHINE = 'oelkerrj'
    FILE_EXTENSION = '.fits'

    # update for different data products
    STAR = 'WASP3b'

    # steps to skip
    CLEAN_SKIP = 'Y'
    WRITE_SKY = 'N'
    PHOTOMETRY_SKIP = 'N'
    CALIBRATE_SKIP = 'Y'

    # get image information
    AXS_X = 1676
    AXS_Y = 1266
    AXS = 2048
    GAIN = 1  # in e-/ADU

    # rough size of PSF for star finding and clip functions, also avoidance for edge
    PSF_X = 40
    PSF_Y = 20
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
    FWHM = 10.  # fwhm of the image
    THRESHOLD = 7.  # the threshold for a source above the background
    APER_SIZE = 32  # size of the single aperture to use
    ANNULI_INNER = 34  # size of the inner annulus for the background (only used in local mode)
    ANNULI_OUTER = 36  # size of the outer annulus for the background (only used in local mode)

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

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY,
                   QUERIES_DIRECTORY, CLEAN_DIRECTORY, MASTER_DIRECTORY]
