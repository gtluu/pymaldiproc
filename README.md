# pyMALDIproc

## About

pyMALDIproc is a library designed for mass spectrometrists that aggregates methods from different libraries to allow 
users to easily process and work with MALDI-TOF dried droplet data, similar to the [MALDIquant](https://cran.r-project.org/web/packages/MALDIquant/index.html) R package.

## Installation

#### Install Anaconda and Git on Windows.

1. Download and install Anaconda for [Windows](https://repo.anaconda.com/archive/Anaconda3-2021.11-Windows-x86_64.exe). 
Follow the prompts to complete installation.
2. Download and install [Git](https://git-scm.com/downloads) and ensure that the option to enable symbolic links was 
checked during installation if not already installed.
3. Run ```Anaconda Prompt```.

#### Install Anaconda and Git on Linux.

1. Download and install Anaconda for [Linux](https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Linux-x86_64.sh). 
Follow the prompts to complete installation.
```
wget https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Linux-x86_64.sh
bash Anaconda3-2023.07-2-Linux-x86_64.sh
```
2. Download and install ```git``` if not already installed. On Ubuntu 22.04.3 LTS:
```
sudo apt-get install git
```

#### Set Up ```conda env```

4. Create a conda instance.
```
conda create -n pymaldiproc python=3.11
```
5. Activate conda environment.
```
conda activate pymaldiproc
```

#### Install pyMALDIproc

6. Install dependencies.
```
pip install -r https://raw.githubusercontent.com/gtluu/pymaldiproc/main/requirements.txt
```
7. Install pyMALDIproc.
```
pip install git+https://github.com/gtluu/pymaldiproc
```

## Usage

#### Using pyMALDIproc

```python
from pymaldiproc.classes import OpenMALDISpectrum
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
preprocessed_data = bin_spectra(preprocessed_data, 10000, 2000, 20000)

# Get feature matrix with peak lists for all spectra in the dataset.
feature_matrix = get_feature_matrix(peaks)
# Export feature matrix
export_feature_matrix(feature_matrix, 'feature_matrix.csv')
```

#### Using the pyMALDIproc Dashboard

Run the following command with the ```pymaldiproc``` venv activated.
```
python [path to pymaldiproc]/dashboard.py
```
The dashboard will be available at ```https://localhost:8050```.
