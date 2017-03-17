# gnps-calour
[Calour](https://github.com/amnona/Calour) interface to GNPS metabolomics database for MS1 tables

## Installation:
1. Install Calour:

Follow instructions [Here](https://github.com/amnona/Calour)

2. Install gnps-calour:

```
pip install git+git://github.com/amnona/gnps-calour
```

## Usage:
Easiest way to use calour is via the [EZCalour](https://github.com/amnona/EZCalour) GUI.

Alternatively, Calour can be used via Jupyter notebook/ipython.

GNPS-Calour is automatically used for metabolomics tables if a corresponding gnps annotations file is available and specified (see below)

Linking between the MS1 table and the GNPS spectra is done via an MZ/RT window (default is MZ+-0.1, RT+-30).

If availble, you will be able to see the GNPS feature annotations as the feature name in the heatmap, and double clicking on a gnps description in the database list window will link to the corresponding GNPS webpage for the metabolite.

**How to obtain the GNPS annotations file**
In order to get GNPS annotations in the heatmap, you need to first upload your study to [GNPS]((gnps.ucsd.edu) and

Once ready, download from [gnps](gnps.ucsd.edu) / full job description / Export / Download clustered data,

Select Download/Download tab delimited all results download

unzip the zip file

use the file "cluterinfosummarygroup_attributes_withids_withcomponentID/XXXX.tsv" as the GNPS file in the load menu

