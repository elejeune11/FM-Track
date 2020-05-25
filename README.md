# FM-Track

<img src="https://github.com/elejeune11/FM-Track/blob/master/misc/sampleimage.png" width="400" height="400" align="right">

FM-Track is a feature-based fiducial marker tracking software for applications in cell mechanics.

Research methods in mechanobiology that involve tracking the deformation of fiducial markers in the vicinity of a cell are increasing in popularity. Here we present a software called FM-Track, a software to facilitate feature-based particle tracking tailored to applications in cell mechanics. FM-Track contains functions for pre-processing images, running fiducial marker tracking, and simple post-processing and visualization. We expect that FM-Track will help researchers in mechanobiology and related fields by providing straightforward and extensible software written in the Python language.

For a lengthier description, please refer to the [technical overview](technicaloverview.pdf).

## Installation

To install FM-Track, it is recommended to use Pip within a conda environment. Anaconda allows you to create [environments](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) to make software dependencies and installations easier to manage. By installing FM-Track into a conda environment, Anaconda automatically configures the proper dependencies for you and installs them all into the correct location. Once installed, FM-Track can be run using [Python import statements](https://docs.python.org/3/reference/import.html). This allows you to import FM-Track functions from the Anaconda installation folder so they can be used in an external file.

Currently, FM-Track is only compatible with Python < 3.7. To create a compatible environment, use the following commands:

```
conda create -n fmtrack python=3.6 anaconda
conda activate fmtrack
conda install -c conda-forge latexcodec
conda install pip
```

To install FM-Track, download the zip file for the repository, unzip the folder, and navigate into the folder inside of a bash shell (such as terminal) using the following command:

```
cd <path-to-folder>
```

To complete the installation, run the following command:

```
sudo pip install .
```

Once running this command, you may delete the downloaded folder if you wish. The [source code folder](fmtrack) has now been copied into the proper Anaconda environment location, so you no longer need to keep a local version.

## Usage

Below is a short example script intended for use on the [data](examples/data) folder. It (1) uses FM-Track to construct the cellular boundaries from images, (2) uses FM-Track to locate the bead positions from images, and (3) runs the actual tracking algorithm.

```
import fmtrack

# set microscope dimensions
X_DIM = 149.95; Y_DIM = 149.95; Z_DIM = 140

# (1) compute cellular boundary from images
print('Importing cell data')
cell_init = fmtrack.FMMesh()
cell_init.get_cell_surface('./data/CytoD/Cell/Gel 2 CytoD%s.tif',0, X_DIM, Y_DIM, Z_DIM, 1.0)
cell_final = fmtrack.FMMesh()
cell_final.get_cell_surface('./data/Normal/Cell/Gel 2 Normal%s.tif',0, X_DIM, Y_DIM, Z_DIM, 1.0)

# (2) find bead positions from images
print('Importing bead data')
beads_init = fmtrack.FMBeads()
beads_init.get_bead_centers('./data/CytoD/Beads/Gel 2 CytoD%s.tif', 1, X_DIM, Y_DIM, Z_DIM)
beads_final = fmtrack.FMBeads()
beads_final.get_bead_centers('./data/Normal/Beads/Gel 2 Normal%s.tif', 1, X_DIM, Y_DIM, Z_DIM)

# (3) run tracking algorithm
tracker = fmtrack.FMTracker(cell_init=cell_init, cell_final=cell_final, beads_init=beads_init, beads_final=beads_final)
tracker.run_tracking()
tracker.save_all('./data/Track_CytoD_to_Normal')
```

To summarize the [technical overview](technicaloverview.pdf), the most basic version of FM-Track accepts a set of bead locations in an initial state and a set in a final state, then it determines which beads from each set can be confidently matched. To do this, FM-Track uses a number of objects to make it more easily usable.

* `FMBeads()`: stores a single data set of beads.
* `FMMesh()`: stores the data for a single cell.
* `FMTracker()`: accepts two sets of bead data, one initial and one final, and matches the beads. Optionally, it can additionally accept two cells.
* `FMPlot()`: stores the parameters for particular types of plots you might want to create (2D and 3D)
* `InputInfo()`: enables the user to run the tracking algorithm on multiple datasets in sequence. This class uses the native filesystem, which is specified in the [test.py](examples/test.py) script.

To understand FM-Track more thoroughly, please refer to the [examples](examples) folder.

## Built With

* [py-earth](https://github.com/scikit-learn-contrib/py-earth) - Used to compute multivariate adaptive regression splines
* [scikit-learn](https://scikit-learn.org/stable/) - Used to generate Gaussian Process Regression models
* [pyvista](http://www.pyvista.org/) - Used to generate 3D deformation field graphs

## Authors

* Emma Lejeune
* Alex Khang
* Jake Sansom

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
