import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hydroRaVENS",
    version="1.1.0",
    author="Andrew D. Wickert",
    author_email="awickert@umn.edu",
    description="Rain and Variable Evapotranspiration, Nieve, and Streamflow",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MNiMORPH/hydroRaVENS/",
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Intended Audience :: Science/Research",
    ],
    entry_points = {
        'console_scripts': [
            'hydroravens = hydroravens.cli:RunHydroRaVENS',
        ],
    },
    keywords='hydrology linear reservoir rainfall runoff snow streamflow discharge',
    project_urls={
        'Model page': 'https://csdms.colorado.edu/wiki/Model:hydroRaVENS',
    },
)
