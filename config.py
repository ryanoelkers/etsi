""" This serves as the configuration file for the ETSI pipeline. """


class Configuration:

    # McDonald Observatory information
    OBS_LAT = 30.6797
    OBS_LON = -104.0247
    OBS_HGT = 2077

    # pipeline choices
    ALIGNMENT = 'N'
    CENTROID = 'Y'
    COADD_IMAGE = 'N'

    # star / observation information
    STAR = 'HAT-P-44b'
    RA = 213.1440475
    DEC = 47.0147347

    # observation information
    DATE = '2022-04-23'
    MONTH = 'April'
    DST = 'Y'  # is it daylight savings time? Also, make sure you are in Chicago time!
    YEAR = 2022

    # exposure time and binning information
    EXPOSURE_TIME = 10  # the exposure time for a single image
    EXPOSURE_COADD = 1  # how many images were pre-added?
    BIN_TIME = 120  # how many seconds of integration do you want in your final co-add?

    # image information
    AXIS_X = 2048
    AXIS_Y = 2048
    GAIN = 0.61  # in e-/ADU
    FILE_EXTENSION = '.fits'
    BEAM_TYPE = 'transmission'

    if BEAM_TYPE == 'transmission':
        # aperture information
        ELIP_APER_A = 40
        ELIP_APER_B = 25
        SKY_POS_ABV = 150
        SKY_POS_BLW = 150
    else:
        # aperture information
        ELIP_APER_A = 65
        ELIP_APER_B = 30
        SKY_POS_ABV = 60
        SKY_POS_BLW = 60

    # PSF information (transmitted band passes)
    NUM_PSF_TRANSMISSION = 8
    WAVELENGTHS_TRANSMISSION = ['937nm', '763nm', '660nm', '587nm', '533nm', '494nm', '467nm', '435nm']

    # PSF information (reflected band passes)
    NUM_PSF_REFLECTION = 7
    WAVELENGTHS_REFLECTION = ['873nm', '713nm', '620nm', '559nm', '512nm', '476nm', '448nm']

    if BEAM_TYPE == 'transmission':
        NUM_PSF = NUM_PSF_TRANSMISSION
        WVE_STR = WAVELENGTHS_TRANSMISSION
    else:
        NUM_PSF = NUM_PSF_REFLECTION
        WVE_STR = WAVELENGTHS_REFLECTION

    # output paths for logging, temporary files, figures etc
    WORKING_DIRECTORY = "\\etsi\\"
    LOG_DIRECTORY = WORKING_DIRECTORY + 'logs\\'

    # input paths for data etc
    if MONTH == 'April':
        if YEAR == 2022:
            DATA_DIRECTORY = ""
        else:
            DATA_DIRECTORY = ""
    elif MONTH == 'May':
            DATA_DIRECTORY = ""
    elif (MONTH == 'Sep') | (MONTH == 'Oct'):
        if YEAR == 2022:
            DATA_DIRECTORY = ""
        else:
            DATA_DIRECTORY = ""
    elif MONTH == 'June':
        if YEAR == 2022:
            DATA_DIRECTORY = ""
        else:
            DATA_DIRECTORY = ""
    elif MONTH == 'July':
        DATA_DIRECTORY = ""
    else:
        print('Warning no month with data!')

    # directories to be generated in the reduction
    RAW_DIRECTORY = DATA_DIRECTORY + "raw\\" + BEAM_TYPE + "\\"
    COADD_DIRECTORY = DATA_DIRECTORY + "coadd\\"
    COADD_BEAM_DIRECTORY = COADD_DIRECTORY + BEAM_TYPE + "\\"
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + "lc\\"
    MISC_DIRECTORY = DATA_DIRECTORY + "misc\\"

    # directory_list
    DIRECTORIES = [DATA_DIRECTORY, LOG_DIRECTORY,
                   LIGHTCURVE_DIRECTORY, RAW_DIRECTORY, MISC_DIRECTORY,
                   COADD_DIRECTORY, COADD_BEAM_DIRECTORY]
