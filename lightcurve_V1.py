import astropy
import astropy.io.fits as pyfits
import numpy
import numpy as np
import pandas
import csv
import matplotlib.pyplot as plt
import statistics
import sys

for fits_fn in sys.argv[1:]:
    #fits_fn = sys.argv[1]  # this is for running all data through command line
    #fits_fn = "4813000020.fits"
    hdulist = pyfits.open(fits_fn)
    print("loading fits file")
    hdulist.info()  # getting info of fits file


    # making table of lightcurve data, deleting columns with the formatting pandas doesn't like,
    # and creating a pandas dataframe of the lightcurve data

    lightcurve_tbl = astropy.table.Table.read(hdulist['LIGHTCURVE'])

    for col in ['date_cmd','date_avg','date_end','ingest_date']:
      del lightcurve_tbl[col]

    lightcurve = lightcurve_tbl.to_pandas()
    lightcurve.info()
    print("loading table")


    # turning the obsepoch data into a pandas dataframe

    obsepoch_tbl = astropy.table.Table.read(hdulist['OBS_EPOCH'])
    obsepoch = obsepoch_tbl.to_pandas()
    print("making table from obs_epoch")


    # creating a dataframe with the lightcurve data per epoch

    output_df = pandas.DataFrame()  # creating an empty dataframe
    epoch_lightcurves = []
    print("making columns in empty dataframe")

    for i, epoch in obsepoch.iterrows():

        is_in_this_epoch = (lightcurve['mjd'] >= epoch['mjd_start']) & (lightcurve['mjd'] <= epoch['mjd_end'])
        # creating a block of data within the start and end dates
        print("defining epoch")
        lightcurve_this_epoch = lightcurve[is_in_this_epoch]  # getting the lightcurve for the data block

    # saving lightcurve per epoch data into the empty pandas dataframe created above
        output_df.loc[i, 'mjd_start'] = epoch['mjd_start']
        output_df.loc[i, 'mjd_end'] = epoch['mjd_end']
        output_df.loc[i, 'mjd_mean'] = (epoch['mjd_start'] + epoch['mjd_end']) / 2.
        output_df.loc[i, 'n_expected'] = (epoch['mjd_end'] - epoch['mjd_start']) * 1440. / 40.
        # how many frames expected for the number of days, with exposure time of 40 minutes
        output_df.loc[i, 'n_actual'] = numpy.sum(is_in_this_epoch)
        outlier_threshold = 5 # number of actual data points
        output_df.loc[i, 'good_data'] = output_df.loc[i, 'n_actual'] > outlier_threshold
        single_epoch_lightcurve = {}

        for mag_name in ['mag_auto', 'mag_aper1', "mag_aper2", "phot_aper2", "magzero", "magzero_std"]:
            mags = lightcurve_this_epoch[mag_name]
            mjds = lightcurve_this_epoch['mjd']

            # save the straight average without outlier rejection
            avgmag = numpy.average(mags)
            out_name = "avg_" + mag_name
            output_df.loc[i, out_name] = avgmag

            # add outlier rejection
            good_data = output_df.loc[i, 'n_actual'] > outlier_threshold
            bad_values = (mags > 50)
            clean_avgmag = numpy.average(mags[~bad_values])
            output_df.loc[i, "cleanavg_" + mag_name] = clean_avgmag
            clean_sigma_raw = numpy.nanpercentile(mags[~bad_values], [16,50,84])
            output_df.loc[i, "cleanmedian_" + mag_name] = clean_sigma_raw[1]
            output_df.loc[i, "cleansigma_" + mag_name] = 0.5*(clean_sigma_raw[2]-clean_sigma_raw[0])

            single_epoch_lightcurve[mag_name] = (mags, output_df.loc[i, 'good_data'], mjds)

        # output_df.info()

        # for num in (output_df['avg_mag_aper2']):  # getting rid out outliers
        #     if num > 50:
        #         print("found num greater than 50")
        #         output_df['avg_mag_aper2'] = output_df['avg_mag_aper2'].replace([num], np.NaN)
        #         #print(num)
        #     print("removed outliers")
        #
        #     # getting variance of avg_mag_aper2 data
        #     variance_avgmag_aper2 = statistics.variance(range(output_df['avg_mag_aper2']))

        epoch_lightcurves.append(single_epoch_lightcurve)
        print("made single epoch lightcurve")



    # creating csv file of data frame
    # display(output_df)
    output_fn = fits_fn[:-5] + ".csv"
    print(output_fn)
    output_df.to_csv(output_fn)


    # plotting light curve
    fig = plt.figure()
    ax = fig.add_subplot(111)
    min_clean = numpy.min(output_df['cleanavg_mag_aper2'])
    max_clean = numpy.max(output_df['cleanavg_mag_aper2'])

    def mjd2year(mjd):
        return (mjd-54000)/365 + 1999.

    good_epoch = output_df['n_actual'] > outlier_threshold
    clean_output = output_df[good_epoch]
    ax.set_ylim((min_clean-0.5, max_clean+0.5))
    ax.set_title('mjd_mean   vs  cleanavg_mag_aper2')
    #bad_average = output_df['avg_mag_aper2'] > 50
    #ax.scatter(mjd2year(output_df['mjd_mean'][~bad_average]), output_df['avg_mag_aper2'][~bad_average])
    ax.scatter(mjd2year(clean_output['mjd_mean']), clean_output['cleanavg_mag_aper2'], c='red')
    #ax.errorbar(mjd2year(output_df['mjd_mean']), output_df['cleanavg_mag_aper2'], yerr=output_df['cleansigma_mag_aper2'])
    ax.plot(mjd2year(clean_output['mjd_mean']), clean_output['cleanavg_mag_aper2'])
    # saving graph
    plt.savefig(fits_fn[:-5] + ".png", dpi=300)
    plt.show()


