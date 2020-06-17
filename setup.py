from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='fm-track',
      version='2.0.0-dev',
      description='Feature-based fiducial marker tracking software for applications in cell mechanics',
      long_description=readme(),
      long_description_content_type='text/markdown',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research'
      ],
      keywords='fiducial marker cell mechanics tracking feature',
      url='https://github.com/elejeune11/FM-Track',
      author='Emma Marie Lejeune',
      author_email='emma.lejeune.11@gmail.com',
      license='MIT',
      packages=['fmtrack'],
      python_requires='>=3.0, <3.7',
      package_data={
        'fmtrack': ['data/*.tif'],
        'fmtrack': ['el_papers.mplstyle']
      },
      install_requires=[
          'markdown',
          'numpy',
          'matplotlib',
          'scipy',
          'sklearn-contrib-py-earth',
          'scikit-image',
          'scikit-learn',
          'pyvista',
          'tqdm',
          'meshio',
          'trimesh',
          'latex'
      ],
      setup_requires=['numpy'],
      include_package_data=True,
      zip_safe=False)
