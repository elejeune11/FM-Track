# FM-Track

FM-Track is a feature-based fiducial marker tracking software for applications in cell mechanics.

Research methods in mechanobiology that involve tracking the deformation of fiducial markers in the vicinity of a cell are increasing in popularity. Here we present a software called FM-Track, a software to facilitate feature-based particle tracking tailored to applications in cell mechanics. FM-Track contains functions for pre-processing images, running fiducial marker tracking, and simple post-processing and visualization. We expect that FM-Track will help researchers in mechanobiology and related fields by providing straightforward and extensible software written in the Python language.

For a lengthier description, please refer to the [technical overview](technicaloverview.pdf).

## Installation

To install FM-Track, it is recommended to use Pip within a conda environment. Currently, FM-Track is only compatible with Python <= 3.6. To create a compatible environment, use the following commands:

```
conda create -n fmtrack python=3.6 anaconda
conda activate fmtrack
conda install pip
```

To install FM-Track, download the zip file for the repository, unzip the folder, and navigate into the folder inside of a bash shell using the following command:

```
cd <path-to-folder>
```

To complete the installation, run the following command:

```
sudo python setup.py install
```

## Testing and Usage

The easiest way to learn how to use FM-Track is by looking at the [test.py](examples/test.py) script. This script runs FM-Track on a set of sample data so you can see how it works on test data before trying it on your own. To properly use this file, copy the [data](examples/data) folder to your Desktop (or somewhere easily accessible), then change the root_directory variable at the top of test.py to store the path of this data folder. Next, run test.py using an IDE or in terminal by running:

```
conda activate fmtrack
python <path-to-test.py>
```

FM-Track uses the input_info class to store all your filepaths. The test.py file explains this in greater detail. Additionally, input_info stores the values of other tunable parameters. The test.py file explains this as well. The parameters that can be modified are as follows:

* Color channels of beads and cell
* Dimensions of FOV
* Cellular threshold for analyzing boundaries of cellular material
* Number of feature vectors used
* Number of nearest beads analyzed
* Tracking type (translation correction or not)
* Cell buffers (for both the regular and translation correcting steps)
* Data plotting options

For a greater description of all of these, please refer to the [technical overview](technicaloverview.pdf).

## Built With

* [py-earth](https://github.com/scikit-learn-contrib/py-earth) - Used to compute multivariate adaptive regression splines
* [scikit-learn](https://scikit-learn.org/stable/) - Used to generate Gaussian Process Regression models

## Versioning

We use [SemVer](http://semver.org/) for versioning.

## Authors

* Emma Lejeune
* Alex Khang

## Contributors

* Jake Sansom

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
