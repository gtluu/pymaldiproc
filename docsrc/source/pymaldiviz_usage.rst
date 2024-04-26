pyMALDIviz Usage
=================
pyMALDIviz provides a simple Dash dashboard that can be used to visualize MALDI data.

pyMALDIviz can be launched by downloading the compiled package from _____ and starting ``pyMALDIviz.exe``.
Alternatively, after installing pyMALDIproc and pyMALDIviz as detailed in the
:doc:`installation instructions <installation>`, pyMALDIviz can be launched from the command line via the
``pyMALDIviz`` command with the appropriate ``conda`` environment activated.

Once loaded, you will see the main pyMALDIviz window. Buttons to load in data, perform preprocessing, and edit
preprocessing parameters are present.

.. image:: imgs/pymaldiviz_01.png
   :alt: pyMALDIviz Main Window

Data can be loaded by clicking on the ``Load *.mzML File`` or ``Load Bruker *.d File`` buttons, which will open one of
the file dialogues shown below to open an *.mzML file or *.d directory, respectively.

.. image:: imgs/pymaldiviz_02.png
   :alt: mzML File Selection Dialogue Window

.. image:: imgs/pymaldiviz_03.png
   :alt: Bruker .d Directory Selection Dialogue Window

The selected MALDI data will be loaded into pyMALDIviz, and a dropdown menu will appear where the loaded spectra can be
selected.

.. image:: imgs/pymaldiviz_04.png
   :alt: Spectrum Selection Dropdown Menu

Selecting a spectrum from the dropdown can will load an interactive spectrum plot. The zoom level can be adjusted by
clicking and dragging or reset using the controls on the upper right hand corner of the plot. The controls also provide
various options, including saving the current plot to a *.png file.

.. image:: imgs/pymaldiviz_05.png
   :alt: Interactive Spectrum Plot

Preprocessing can be applied to spectra by clicking on the any of the preprocessing buttons above the spectrum plot.
The preprocessing parameters can be edited via the ``Edit Preprocessing Parameters`` button, which will open another
window with preprocessing parameters for each step.

.. image:: imgs/pymaldiviz_06.png
   :alt: Preprocessing Parameters Part 1

.. image:: imgs/pymaldiviz_07.png
   :alt: Preprocessing Parameters Part 2

Any edited parameters can be saved by clicking the ``Save`` button, and a confirmation message will appear once the
parameters have been saved.

.. image:: imgs/pymaldiviz_08.png
   :alt: Edit Preprocessing Parameters Save Confirmation Message

Once preprocessing and peak picking has been applied, a spectrum will similar to the spectrum shown below. In this
example, spectrum trimming, baseline smoothing, baseline removal, intensity normalization, and peak picking with
deisotoping has been applied to the spectrum.

.. image:: imgs/pymaldiviz_09.png
   :alt: Preprocessed and Peak Picked Spectrum

If at any point, the applied preprocessing proves to be unsatisfactory, preprocessing can be cleared using the
``Undo Preprocessing Button``, which will reset the spectrum plot to its default state in which the raw data is
displayed.

.. image:: imgs/pymaldiviz_05.png
   :alt: Reset Spectrum Plot with Undo Preprocessing Button

In addition to being visualized, the peak list seen in the spectrum plot can also be saved to a CSV file for use with
downstream analysis platforms/pipelines by clicking the ``Export Peak List from Labeled Peaks`` button, which will
open a dialogue to save the CSV file.

.. image:: imgs/pymaldiviz_10.png
   :alt: Export Peak List File Dialogue
