""" This serves as the configuration file for the ETSI pipeline. """


class Configuration:

    # Computer for reduction
    MACHINE = 'karl'

    # observatory information
    OBS_LAT = 30.6797
    OBS_LON = -104.0247
    OBS_HGT = 2077

    # calibration steps to skip
    BIAS_SUBTRACT = 'N'
    BIAS_TYPE = 'coadd'
    FLAT_DIVIDE = 'N'
    ALIGNMENT = 'N'
    CENTROID = 'Y'
    CUT_IMAGE = 'N'

    # pipeline choice
    COADD_IMAGE = 'N'
    COADDED = 'Y'
    CLEAN = 'N'
    MAKE_MASTER = 'N'
    PSF_PHOTOMETRY = 'N'
    PER_BAND = 'N'
    BAND = 'aper_2'

    # background subtraction
    SKY_SUBTRACT = 'none'

    # file information
    STAR = 'HAT-P-32b'  #
    DATE = '2022-10-04'  #
    MONTH = 'Oct'  #
    YEAR = 2022
    RA = 31.0428229
    DEC = 46.6878361
    EXPOSURE_TIME = 3
    EXPOSURE_COADD = 0
    BIN_NUM = 1.
    BIN_TIME = 60.
    DST = 'Y'  # is it daylight savings time? Also, make sure you are in Chicago time!

    # image information
    AXIS_X = 2048
    AXIS_Y = 2048
    GAIN = 0.61  # in e-/ADU
    FILE_EXTENSION = '.fits'
    BEAM_TYPE = 'transmission'

    # background subtraction
    PIX = 256

    # rough size of PSF for star finding and clip functions, also avoidance for edge
    SIGMA_PSF = 10
    PSF_THRESHOLD = 500
    if BEAM_TYPE == 'transmission':
        PSF_X = 35
        PSF_Y = 15
        PSF_CUTOUT = 80

        # aperture information
        ELIP_APER_A = 40
        ELIP_APER_B = 25
        SKY_POS_ABV = 100
        SKY_POS_BLW = 100

    else:
        PSF_X = 30
        PSF_Y = 10
        PSF_CUTOUT = 130

        # aperture information
        ELIP_APER_A = 65
        ELIP_APER_B = 30
        SKY_POS_ABV = 60
        SKY_POS_BLW = 60

    # a photometry configuration
    FWHM = 10.  # fwhm of the image
    THRESHOLD = 5.  # the threshold for a source above the background

    # aperture annulus for the sky background automatically determined from the main aperture
    ELIP_ANNULI_A0 = ELIP_APER_A + 2
    ELIP_ANNULI_AN = ELIP_APER_A + 4
    ELIP_ANNULI_B0 = ELIP_APER_B + 2
    ELIP_ANNULI_BN = ELIP_APER_B + 4

    # PSF information (transmitted band passes)
    NUM_PSF_TRANSMISSION = 8
    WAVELENGTHS_TRANSMISSION = ['937nm', '763nm', '660nm', '587nm', '553nm', '494nm', '467nm', '435nm']
    WAVELENGTHS_TRANSMISSION_NUMS = [937, 763, 660, 587, 553, 494, 467, 435]

    # PSF information (reflected band passes)
    NUM_PSF_REFLECTION = 7
    WAVELENGTHS_REFLECTION = ['873nm', '713nm', '620nm', '559nm', '512nm', '476nm', '448nm']
    WAVELENGTHS_REFLECTION_NUMS = [873, 713, 620, 559, 512, 476, 448]

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
    if MONTH == 'April':
        if YEAR == 2022:
            DATA_DIRECTORY = "F:\\ETSI_" + MONTH + "2022\\" + DATE + "\\" + STAR + "\\"
            CALIBRATION_DIRECTORY = "F:\\ETSI_" + MONTH + "2022\\" + DATE + "\\"
        else:
            DATA_DIRECTORY = "D:\\" + MONTH + "2023\\" + DATE + "\\" + STAR + "\\"
            CALIBRATION_DIRECTORY = "D:\\" + MONTH + "2023\\" + DATE + "\\"
    elif MONTH == 'May':
            DATA_DIRECTORY = "D:\\" + MONTH + "2023\\" + DATE + "\\" + STAR + "\\"
            CALIBRATION_DIRECTORY = "D:\\" + MONTH + "2023\\" + DATE + "\\"
    elif (MONTH == 'Sep') | (MONTH == 'Oct'):
        if YEAR == 2022:
            DATA_DIRECTORY = "D:\\" + MONTH + "2022\\" + DATE + "\\" + STAR + "\\"
            CALIBRATION_DIRECTORY = "D:\\" + MONTH + "2022\\" + DATE + "\\"
        else:
            DATA_DIRECTORY = "I:\\" + MONTH + "2023\\" + DATE + "\\" + STAR + "\\"
            CALIBRATION_DIRECTORY = "I:\\" + MONTH + "2023\\" + DATE + "\\"
    elif MONTH == 'June':
        if YEAR == 2022:
            DATA_DIRECTORY = "F:\\ETSI_" + MONTH + "2022\\" + DATE + "\\" + STAR + "\\"
            CALIBRATION_DIRECTORY = "F:\\ETSI_" + MONTH + "2022\\" + DATE + "\\"
        else:
            DATA_DIRECTORY = "G:\\" + MONTH + "2023\\" + DATE + "\\" + STAR + "\\"
            CALIBRATION_DIRECTORY = "G:\\" + MONTH + "2023\\" + DATE + "\\"
    elif MONTH == 'July':
        DATA_DIRECTORY = "G:\\" + MONTH + "2022\\" + DATE + "\\" + STAR + "\\"
        CALIBRATION_DIRECTORY = "G:\\" + MONTH + "2022\\" + DATE + "\\"
    else:
        print('Warning no month with data!')

    # directories to be generated in the reduction
    DARKS_DIRECTORY = CALIBRATION_DIRECTORY + "darks\\"
    BIAS_DIRECTORY = CALIBRATION_DIRECTORY + "Bias-12bit\\"
    FLATS_DIRECTORY = CALIBRATION_DIRECTORY + "flats\\"
    RAW_DIRECTORY = DATA_DIRECTORY + "raw\\" + BEAM_TYPE + "\\"
    COADD_DIRECTORY = RAW_DIRECTORY + "coadd\\"
    CLEAN_DIRECTORY = DATA_DIRECTORY + "clean\\"
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + "lc\\"
    LIGHTCURVE_BAND_DIRECTORY = LIGHTCURVE_DIRECTORY + "band\\"
    MASTER_DIRECTORY = DATA_DIRECTORY + "master\\"

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY,
                   QUERIES_DIRECTORY, CLEAN_DIRECTORY, LIGHTCURVE_DIRECTORY,
                   LIGHTCURVE_BAND_DIRECTORY, COADD_DIRECTORY, RAW_DIRECTORY,
                   FLATS_DIRECTORY, BIAS_DIRECTORY, DARKS_DIRECTORY, MASTER_DIRECTORY]
