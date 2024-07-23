# etsi
This is the pipeline that was used to reduce the data from the Exoplanet Transmission Spectroscopy Imager (ETSI) that is
described in Oelkers et al. 2024. 

If you have ETSI data and would like to get this pipeline working, make sure you have a position file in a 'misc' 
directory and make sure your data is set up as '\date\star_name\raw\transmission\' or '\date\star_name\raw\reflection\'.

When making hand selected apertures, please remember to select the reddest bandpass for each star first! The code assumes
you are going red to blue for each star.
