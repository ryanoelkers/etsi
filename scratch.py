@staticmethod
def mk_flat(image_directory):
    """ This function will make the master flat frame using the provided image list.
    :parameter image_directory - a directory where the images reside for combination

    :return - The bias frame is returned and written to the master directory
    """

    if os.path.isfile(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits") == 0:

        # read in the bias file
        bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY)

        # get the image list
        image_list = Utils.get_file_list(image_directory, '.fits')

        # determine the number of loops we need to move through for each image
        nfiles = len(image_list)
        nbulk = 20

        # get the integer and remainder for the combination
        full_bulk = nfiles // nbulk
        part_bulk = nfiles % nbulk
        if part_bulk > 0:
            hold_bulk = full_bulk + 1
        else:
            hold_bulk = full_bulk

        # here is the 'holder'
        hold_data = np.array(shape=(Configuration.AXS, Configuration.AXS))
        # np.ndarray(shape=(hold_bulk, Configuration.AXS, Configuration.AXS))

        # update the log
        Utils.log("Generating a master flat frame from multiple files in bulks of " + str(nbulk) +
                  " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                  str(hold_bulk) + " mini-files to combine.", "info")

        for kk in range(0, hold_bulk):

            # loop through the images in sets of nbulk
            if kk < full_bulk:
                # generate the image holder
                block_hold = np.ndarray(shape=(nbulk, Configuration.AXS, Configuration.AXS))

                # generate the max index
                mx_index = nbulk
            else:
                # generate the image holder
                block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS, Configuration.AXS))

                # generate the max index
                mx_index = part_bulk

            # make the starting index
            loop_start = kk * nbulk
            idx_cnt = 0

            Utils.log("Making mini-flat file " + str(kk) + ".", "info")

            # now loop through the images
            for jj in range(loop_start, mx_index + loop_start):
                # read in the flat file
                flat_tmp = fits.getdata(Configuration.FLAT_DIRECTORY + image_list[jj], 0)
                if np.ndim(flat_tmp) > 2:
                    flat_tmp = flat_tmp[0]

                # read in the image directly into the block_hold
                block_hold[idx_cnt] = flat_tmp - bias

                # increase the iteration
                idx_cnt += 1

            # median the data into a single file
            hold_data[kk] = np.median(block_hold, axis=0)

        # median the mini-images into one large image
        flat_image = np.median(hold_data, axis=0)
        nflat_image = flat_image / np.median(flat_image[512:1536, 512:1536])

        # pull the header information from the first file of the set
        flat_header = fits.getheader(Configuration.FLAT_DIRECTORY + image_list[0])
        flat_header['comb_typ'] = 'median'
        flat_header['norm_pix'] = '512to1536'
        flat_header['num_comb'] = nfiles

        # write the image out to the master directory
        fits.writeto(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits",
                     nflat_image, flat_header, overwrite=True)
    else:
        nflat_image = fits.getdata(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits", 0)

    return nflat_image