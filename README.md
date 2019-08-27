# FM-Track

FM-Track is a feature-based fiducial marker tracking software for applications in cell mechanics.

Research methods in mechanobiology that involve tracking the deformation of fiducial markers in the vicinity of a cell are increasing in popularity. Here we present a software called FM-Track, a software to facilitate feature-based particle tracking tailored to applications in cell mechanics. FM-Track contains functions for pre-processing images, running fiducial marker tracking, and simple post-processing and visualization. We expect that FM-Track will help researchers in mechanobiology and related fields by providing straightforward and extensible software written in the Python language.

For a lengthier description, please refer to the [technical overview](technicaloverview.pdf).

## Installation

To install FM-Track, it is recommended to use Pip within a conda environment. Currently, FM-Track is only compatible with Python <= 3.6. To create a compatible environment, use the following commands:

```
conda create -n <environmentname> python=3.6 anaconda
conda activate <environmentname>
conda install pip
```

To install FM-Track, simply use the following command:

```
pip install fm-track
```

You can also download the zip file for the repository, unzip the folder, navigate into the folder inside of a bash shell, and run the following command:

```
sudo python setup.py install
```

## Usage

FM-Track uses the input_info class to store filepaths for easy IO. To initialize an input_info object, pass your root file path into the initialization function like so:

```Python
info = input_info('/Users/<username>/Desktop/data')
```

In this example, we will be interacting with data inside of the data folder, located on the Desktop. All of the input and output files will be located in this data directory.

Input_info stores the names of your subdirectories and other tunable parameters. By default, input_info objects are initialized to work with with the sample data provided in the software package (see the Testing section of this document for an example). To change the names of these properties, simply modify the properties of the object like such:

```Python
info.num_feat = 6
```

To view all of the properties you can change, simply look inside the [input_info.py](input_info.py) file.

## Testing

To test FM-Track, it is easiest to run the sofware on the data folder included within the fmtrack directory. You can either download this data from GitHub as a .zip file or copy the folder from your directory. Copy this data directory to an easily accessible testing location. Try running the following code after modifying the \<filepath\> to match the path of your data folder.

```Python
from fmtrack.input_info import input_info
from fmtrack.run_tracking_all_steps import run_tracking_all_steps

info = input_info('<filepath>')
run_tracking_all_steps(True,True,True,info)
```

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
