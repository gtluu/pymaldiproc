Installation
============
To run pyMALDIproc and pyMALDIviz, you can set up a ``conda`` environment to do so. Please note that while pyMALDIproc
and pyMALDIviz will run on any OS for open data formats, analysis of raw Bruker data formats is only supported on
Windows and Linux; macOS is not supported for raw bruker data. A single requirements file is provided to install all
dependencies for both pyMALDIproc and pyMALDIviz.

Installing on Windows
---------------------
1. Download and install `Anaconda for Windows <https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Windows-x86_64.exe>`_ if not already installed. Follow the prompts to complete installation.

2. Download and install `Git for Windows <https://github.com/git-for-windows/git/releases/download/v2.42.0.windows.2/Git-2.42.0.2-64-bit.exe>`_ if not already installed.

3. Run ``Anaconda Prompt``.

4. Create a conda instance.

   .. code-block::

        conda create -n pymaldiproc python=3.11

5. Activate conda environment.

   .. code-block::

        conda activate pymaldiproc

6. Install dependencies.

   .. code-block::

        pip install -r https://raw.githubusercontent.com/gtluu/pymaldiproc/main/requirements.txt

7. Install pyMALDIproc/pyMALDIviz.

   .. code-block::

        pip install git+https://github.com/gtluu/pyMALDIproc.git

8. pyMALDIproc and pyMALDIviz are now ready to use. See the :doc:`pyMALDIproc Usage <pymaldiproc_usage>` and :doc:`pyMALDIviz Usage <pymaldiviz_usage>` pages for more information.

Installing on Linux
-------------------
Please note that while these instructions should apply to most Linux distros, pyMALDIproc is tested on Ubuntu 23.04.3
LTS. We recommend using this distro if you encounter compatibility issues in others.

1. If not already installed, download and install `Anaconda for Linux <https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Linux-x86_64.sh>`_. Anaconda3-2023.07-2 for Linux is used as an example here.

   * Alternatively, the script can be downloaded in the ``Terminal`` using the following command.

   .. code-block::

        wget https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Linux-x86_64.sh

2. If not already installed, install ``git``. On Ubuntu 23.04.3 LTS, this can be done using the following command.

   .. code-block::

        sudo apt-get install git

3. Install Anaconda for Linux via the bash script that was downloaded. After installation, restart the terminal (or open a new terminal window).

   .. code-block::

        bash [path to]/Anaconda3-2023.07-2-Linux-x86_64.sh

4. In the terminal, create a conda virtual environment.

   .. code-block::

        conda create -n pymaldiproc python=3.11

5. Activate conda environment.

   .. code-block::

        conda activate pymaldiproc

6. Install dependencies.

   .. code-block::

        pip install -r https://raw.githubusercontent.com/gtluu/pymaldiproc/main/requirements.txt

7. Install pyMALDIproc/pyMALDIviz.

   .. code-block::

        pip install git+https://github.com/gtluu/pyMALDIproc.git

8. pyMALDIproc and pyMALDIviz are now ready to use. See the :doc:`pyMALDIproc Usage <pymaldiproc_usage>` and :doc:`pyMALDIviz Usage <pymaldiviz_usage>` pages for more information.
