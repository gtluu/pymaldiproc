# pyMALDIproc

## About

pyMALDIproc is a library designed for mass spectrometrists that aggregates methods from different libraries to allow users to easily process and work with MALDI-TOF dried droplet data, similar to the [MALDIquant](https://cran.r-project.org/web/packages/MALDIquant/index.html) R package.

## Installation

#### Install Anaconda on Linux

1. Download and install Anaconda for [Linux](https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh). 
Follow the prompts to complete installation. Anaconda3-2021.11 for Linux is used as an example here.
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash /path/to/Anaconda3-2021.11-Linux-x86_64.sh
```
2. Add ```anaconda3/bin``` to PATH.
```
export PATH=$PATH:/path/to/anaconda3/bin
```

#### Install Anaconda on Windows.

1. Download and install Anaconda for [Windows](https://repo.anaconda.com/archive/Anaconda3-2021.11-Windows-x86_64.exe). 
Follow the prompts to complete installation.
2. Run ```Anaconda Prompt (R-MINI~1)``` as Administrator.

#### Set Up ```conda env```

3. Create a conda instance. You must be using Python 3.7.
```
conda create -n pymaldiproc python=3.7
```
4. Activate conda environment.
```
conda activate pymaldiproc
```

#### Install pyMALDIproc

5. Download pyMALDIproc by cloning the Github repo (you will need to have [Git](https://git-scm.com/downloads) and 
ensure that the option to enable symbolic links was checked during installation). It may be necessary to explicitly
allow for the use of symbolic links by adding the ```-c core.symlinks=true``` parameter on Windows.
```
git clone https://www.github.com/gtluu/pyMALDIproc
or
git clone -c core.symlinks=true https://www.github.com/gtluu/pyMALDIproc
```
6. Install dependencies.
```
# pyMALDIproc dependencies
pip install -r /path/to/pyMALDIproc/requirements.txt
```

## Usage

#### Using pyMALDIproc

```python
from pymaldiproc.classes import MALDISpectrum
from pymaldiproc.data_import import import_mzml
from pymaldiproc.preprocessing import *

# Read in MALDI spectra from mzML file.
data = import_mzml('profile.mzML')  # data is a list of spectra

# Data Preprocessing
preprocessed_data = trim_spectra(data, 2000, 20000)
preprocessed_data = transform_intensity(preprocessed_data)
preprocessed_data = smooth_baseline(preprocessed_data)
preprocessed_data = remove_baseline(preprocessed_data)
preprocessed_data = normalize_intensity(preprocessed_data)
preprocessed_data = bin_spectra(preprocessed_data, 10000, 2000, 20000)
preprocessed_data = align_spectra(preprocessed_data)

# Peak Picking
peaks = peak_picking(preprocessed_data)

# Get feature matrix with peak lists for all spectra in the dataset.
feature_matrix = get_feature_matrix(peaks)
# Export feature matrix
export_feature_list(feature_matrix, 'feature_matrix.csv')
```

#### Using the Example Dashboard

Run the following command with the ```pymaldiproc``` venv activated.
```
python /path/to/pymaldiproc/example_dashboard.py
```
The dashboard will be available at ```https://localhost:8050```.
