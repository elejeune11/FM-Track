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
