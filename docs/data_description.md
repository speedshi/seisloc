# Data description

Here we describe the formats of the files used by `migrationloc`.

## Input data

The following files are used as input to `migrationloc`.  Their names are hard-coded:

* `migpara.dat`
* `soupos.dat`
* `travelp.dat`
* `travels.dat`
* `waveform.dat`

All the input files are stored in the folder `./data`. Double precision is normally used.

1. `migpara.dat`, format: text file. This file contains key parameters (such as the number of stations, time window length and characteristic function type) adopted in `migrationloc`. Detailed information can be found in the [provided example file](../migpara.dat), or in the [section below](#migration-parameter-file-migparadat)

2. `soupos.dat`, format: binary file. This file contains the position information of all the imaging points. The data unit: meter, size: `nsr*3*precision` byte.
   It contains the X, Y and Z coordinates of all the potential source points (imaging points). First is the X coordinates of all the potential source points, then is the Y coordinates of all the potential source points, last is the Z coordinates of all the potential source points (format detail shown in  Figure 1). The storage order of all the source points much keep consistent with the traveltime table file. The number of source points must be same as described in the `migpara.dat` file.

![Data format of `soupos.dat`](images/1.png)  
*Figure 1. Data format of `soupos.dat`*

3. `travelp.dat`, format: binary file. This file is the traveltime table of direct P-waves. The data unit: second, size: `nsr*nre*precision` byte.
   It contains the travetime information of the direct P-waves. First is the traveltimes of different potential source points for the first station (receiver), then for the second station, and so on (format detail shown in Figure 2). The order for different potential source points must keep consistent with the source position file `soupos.dat`. The number of source points and stations must be same as described in the `migpara.dat`` file.

![Data format of `travelp.dat`](images/2.png)  
*Figure 2. Data format of `travelp.dat`*

4. `travels.dat`, format: binary file. This file is the traveltime table of direct S-waves. The data unit: second, size: `nsr*nre*precision`` byte.
   It contains the travetime information of the direct S-waves. First is the traveltimes of different potential source points for the first station (receiver), then for the second station, and so on (format is the same as `travelp.dat`). The order for different potential source points must keep consistent with the source position file `soupos.dat`. The number of source points and stations must be same as described in the `migpara.dat` file.

5. `waveform.dat`, format: binary file. It contains seismic waveform data. The data unit:  unitless, size: `nre*nt*precision` byte.
   It contains the recorded waveform data of different stations. First is the recorded waveforms of different stations for the first recording time, then for the second recording time, and so on (format detail shown in Figure 3). The order for different stations must keep consistent with the traveltime table file `travelp.dat` and `travels.dat`. The number of stations must be same as described in the `migpara.dat` file. The number of sampling times (recording times) must be the same as described in the `migpara.dat` file.

![Data format of `waveformat.dat`](images/3.png)  
*Figure 3. Data format of `waveform.dat`*


## Migration parameter file `migpara.dat`

Below is an example parameter file read by `migrationloc`.

```
0              | migtp (integer)
2              | phasetp (integer)
0              | cfuntp (integer)
15             | nre (integer)
10000          | nsr (integer)
waveform.dat   | dfname (character)
0.001          | dt (real)
3600           | tdatal (real)
1              | tpwind (real)
1              | tswind (real)
0.1            | dt0 (real)
0.4            | vthrd (real)
2              | mcmdim (integer)
200            | spaclim (real)
2              | timelim (real)
1              | nssot (integer)
```

Each parameter must appear on its own on each line.  Comments after the parameter are allowed,
but comment lines are not.

The migration parameters are:

* `migtp` (integer): specify the migration method:
    - `0` for MCM
    - `1` for conventional migration.
* `phasetp` (integer): specify the phases used for migration:
    - `0` for P-phase only
    - `1` for S-phase only
    - `2` for P+S-phases
* `cfuntp` (integer): specify the datatype of characteristic function used for migration:
    - `0` for original data
    - `1` for envelope
    - `2` for absolute value
    - `3` for non-negative value
    - `4` for squared value
    - *For other CFs such as STA/LTA or kurtosis, please pre-process the data.*
* `nre` (integer): number of stations or receivers.
* `nsr` (integer): number of imaging points or potential source locations.
* `dfname` (character): the filename of the input seismic data.  Cannot exceed 20 characters.
* `dt` (real): time sampling interval of the input seismic data in second (s).
* `tdatal` (real): time length of the whole seismic data in second (s).
* `tpwind` (real): P-phase time window length in second (s) used for migration.
* `tswind` (real): S-phase time window length in second (s) used for migration.
* `dt0` (real): time sampling interval of searching origin times in second (s).
* `vthrd` (real): threshold value for identifying seismic event in the migration volume.
    - For MCM, should between `0` and `1`.
    - If `vthrd <= 0`, adaptive thresholding is used: `mean + 3*std`.
    - If `vthrd >= 1`, adaptive thresholding is as: `mean + vthrd*std`.
* `mcmdim` (integer): the dimension of MCM, should `>= 2`. Need this parameter when migtp is 0. For large arrays, do not exceed 5.
* `spaclim` (real): the space limit in searching for potential seismic events, in meter (m).
* `timelim` (real): the time limit in searching for potential seismic events, in second (s).
* `nssot` (integer): the maximum number of potential seismic events can be accept for a single origin time. If it is less then 0, then use default value: 100.
