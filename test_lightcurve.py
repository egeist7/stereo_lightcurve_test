import astropy
import astropy.io.fits as pyfits
import numpy
import pandas
import csv
import matplotlib.pyplot as plt


fits_fn = "4813000025.fits"
hdulist = pyfits.open(fits_fn)
hdulist.info()  # getting info of fits file
print("loading fits file")



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
# obsepoch.info()
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
    # lightcurve_this_epoch.info()

    # saving lightcurve per epoch data into the empty pandas dataframe created above
    output_df.loc[i, 'mjd_start'] = epoch['mjd_start']
    output_df.loc[i, 'mjd_end'] = epoch['mjd_end']
    output_df.loc[i, 'mjd_mean'] = (epoch['mjd_start'] + epoch['mjd_end']) / 2.
    output_df.loc[i, 'n_expected'] = (epoch['mjd_end'] - epoch['mjd_start']) * 1440. / 40.
    # how many frames expected for the number of days, with exposure time of 40 minutes
    output_df.loc[i, 'n_actual'] = numpy.sum(is_in_this_epoch)
    # output_df.loc[i, 'local_zp'] = epoch['local_zp']



    single_epoch_lightcurve = {}

    for mag_name in ['mag_auto', 'mag_aper1', "mag_aper2", "phot_aper2", "magzero", "magzero_std"]:

        outlier_threshold = 5
        mags = lightcurve_this_epoch[mag_name]
        mjds = lightcurve_this_epoch['mjd']
        # start with all data that has good values
        good_values = numpy.isfinite(mags)

        # for iteration in range(3):
        #     # get percentiles
        #     # 16th --> -1 sigma
        #     # 50th --> median
        #     # 84th --> +1 sigma
        #     _stats = numpy.nanpercentile(mags[good_values], [16, 50, 84])
        #     median = _stats[1]
        #     sigma = 0.5 * (_stats[2] - _stats[0])
        #
        #     # good data needs to be within 3 sigma of median
        #     good_values = (mags < (median + outlier_threshold * sigma)) & (mags > (median - outlier_threshold * sigma))

        out_name = "avg_" + mag_name
        avgmag = numpy.average(mags[good_values])
        output_df.loc[i, out_name] = avgmag
        output_df.loc[i, mag_name + "_ngood"] = numpy.sum(good_values)
        output_df.loc[i, mag_name + "_nbad"] = numpy.sum(~good_values)

        single_epoch_lightcurve[mag_name] = (
        mags, good_values, mjds, median, sigma, lightcurve_this_epoch['magzero_std'])

    epoch_lightcurves.append(single_epoch_lightcurve)
    print("made single epoch lightcurve")



# output_df.info()
# display(output_df)
display(output_df)
output_fn = fits_fn[:-5] + ".csv"
print(output_fn)
output_df.to_csv(output_fn)



fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('mjd_mean   vs   avg_mag_auto')
ax.scatter(output_df['mjd_mean'], output_df['avg_mag_auto'])
ax.plot(output_df['mjd_mean'], output_df['avg_mag_auto'])



fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('mjd_mean   vs   avg_phot_aper2')
ax.scatter(output_df['mjd_mean'], output_df['avg_phot_aper2'])
ax.plot(output_df['mjd_mean'], output_df['avg_phot_aper2'])