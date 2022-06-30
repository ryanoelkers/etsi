""" This serves as the configuration file for the ETSI pipeline. """


class Configuration:

    # Computer for reduction
    MACHINE = 'barristan'

    # calibration steps to skip
    BIAS_SUBTRACT = 'Y'
    DARK_SUBTRACT = 'N'
    FLAT_DIVIDE = 'N'
    SKY_SUBTRACT = 'Y'
    ALIGNMENT = 'N'
    CUT_IMAGE = 'N'
    NUM_MASTER_FILES = 100

    # major steps to skip
    CLEAN_SKIP = 'Y'
    WRITE_SKY = 'N'
    PHOTOMETRY_SKIP = 'N'
    LIGHTCURVE_SKIP = 'Y'
    COLOR_SKIP = 'Y'

    # stellar parameters
    STAR = 'XO-1'
    RA = 213.1438604
    DEC = 47.0147826
    TC = 2455696.93695
    PERIOD = 4.301219

    # observation & image specific parameters
    EXPOSURE_TIME = 10.
    AXIS_X = 2048
    AXIS_Y = 2048
    AXS_X = 2048
    AXS_Y = 2048
    PIX = 128
    GAIN = 0.61  # in e-/ADU
    FILE_EXTENSION = '.fits'

    # photometry parameters
    PHOTOMETRY = 'PSF'
    APERTURE_SHAPE = 'ellipse'
    TIME = 'phase'
    BEAM_TYPE = 'transmission'

    # rough size of PSF for star finding and clip functions, also avoidance for edge
    SIGMA_PSF = 10
    PSF_THRESHOLD = 500
    PSF_X = 40
    PSF_Y = 25
    PSF_CUTOUT = 80

    # a photometry configuration
    FWHM = 10.  # fwhm of the image
    THRESHOLD = 5.  # the threshold for a source above the background

    # aperture information
    ELIP_APER_A = 40
    ELIP_APER_B = 25

    # aperture annulus for the sky background automatically determined from the main aperture
    ELIP_ANNULI_A0 = ELIP_APER_A + 2
    ELIP_ANNULI_AN = ELIP_APER_A + 4
    ELIP_ANNULI_B0 = ELIP_APER_B + 2
    ELIP_ANNULI_BN = ELIP_APER_B + 4

    # PSF information (transmitted band passes)
    NUM_PSF_TRANSMISSION = 8
    PSF1_TRANSMISSION = 0
    PSF2_TRANSMISSION = PSF1_TRANSMISSION + 80
    PSF3_TRANSMISSION = PSF1_TRANSMISSION + 160
    PSF4_TRANSMISSION = PSF1_TRANSMISSION + 250
    PSF5_TRANSMISSION = PSF1_TRANSMISSION + 340
    PSF6_TRANSMISSION = PSF1_TRANSMISSION + 430
    PSF7_TRANSMISSION = PSF1_TRANSMISSION + 530
    PSF8_TRANSMISSION = PSF1_TRANSMISSION + 630
    PSFS_TRANSMISSION = [PSF1_TRANSMISSION, PSF2_TRANSMISSION, PSF3_TRANSMISSION, PSF4_TRANSMISSION,
                         PSF5_TRANSMISSION, PSF6_TRANSMISSION, PSF7_TRANSMISSION, PSF8_TRANSMISSION]
    WAVELENGTHS_TRANSMISSION = ['937nm', '763nm', '660nm', '587nm', '533nm', '494nm', '467nm', '435nm']

    # PSF information (reflected band passes)
    NUM_PSF_REFLECTION = 7
    PSF1_REFLECTION = 0
    PSF2_REFLECTION = PSF1_REFLECTION + 85
    PSF3_REFLECTION = PSF1_REFLECTION + 170
    PSF4_REFLECTION = PSF1_REFLECTION + 255
    PSF5_REFLECTION = PSF1_REFLECTION + 340
    PSF6_REFLECTION = PSF1_REFLECTION + 425
    PSF7_REFLECTION = PSF1_REFLECTION + 510

    PSFS_REFLECTION = [PSF1_REFLECTION, PSF2_REFLECTION, PSF3_REFLECTION, PSF4_REFLECTION,
                       PSF5_REFLECTION, PSF6_REFLECTION, PSF7_REFLECTION]
    WAVELENGTHS_REFLECTION = ['448nm', '476nm', '512nm', '559nm', '620nm', '713nm', '873nm']

    if BEAM_TYPE == 'transmission':
        NUM_PSF = NUM_PSF_TRANSMISSION
    else:
        NUM_PSF = NUM_PSF_REFLECTION

    # output paths for logging, temporary files, figures etc
    WORKING_DIRECTORY = "C:\\Users\\barristan\\Development\\etsi\\"
    ANALYSIS_DIRECTORY = WORKING_DIRECTORY + 'analysis\\'
    LOG_DIRECTORY = WORKING_DIRECTORY + 'logs\\'
    QUERIES_DIRECTORY = WORKING_DIRECTORY + 'queries\\'

    # input paths for data etc
    DATA_DIRECTORY = "F:\\ETSI_June2022\\2022-06-22\\" + STAR + "\\"
    CALIBRATION_DIRECTORY = "F:\\ETSI_June2022\\2022-06-22\\"

    # directories to be generated in the reduction
    DARKS_DIRECTORY = CALIBRATION_DIRECTORY + "darks\\"
    BIAS_DIRECTORY = CALIBRATION_DIRECTORY + "Bias-12bit\\"
    FLATS_DIRECTORY = CALIBRATION_DIRECTORY + "flats\\"
    RAW_DIRECTORY = DATA_DIRECTORY + "raw\\" + BEAM_TYPE + "\\"
    CLEAN_DIRECTORY = DATA_DIRECTORY + "clean\\"
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + "lc\\"
    MASTER_DIRECTORY = DATA_DIRECTORY + "master\\"

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY,
                   QUERIES_DIRECTORY, CLEAN_DIRECTORY, LIGHTCURVE_DIRECTORY,
                   RAW_DIRECTORY, FLATS_DIRECTORY, BIAS_DIRECTORY, DARKS_DIRECTORY, MASTER_DIRECTORY]
