
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4383138.svg)](https://doi.org/10.5281/zenodo.4383138)

![genGEO logo](media/logo.jpg)

genGEO is a coupled reservoir-electricity-cost geothermal simulator. It is an open-source, object-oriented python code which may be extended to include any geothermal or power generation scenario.

#### Example
An example file is provided with [exampleCPGSystemStudy.py](exampleCPGSystemStudy.py). This example implements a CPG system and solves for a number of reservoir depths and permeabilities. It contains a large number of scenarios and therefore runs for a long time. Hence, it should only be used as an example to derive other studies from.

Install
---
genGEO uses [conda](https://docs.conda.io/projects/conda/en/latest/index.html) for easy package management and comes with an yml file to create an environment with all required packages.

For conda users please run the following for installation:
- `conda env create -f genGEO.yml`

and for updates:
- `conda env update -f genGEO.yml`

For non conda users, please install the following packages:
- python=3.7
- pandas
- scipy
- coolprop
- xlrd

Testing
---
To test if genGEO works fine on your system, please run the tests with the following command:

`cd <path_to_genGEO> ; ./runTests.sh`

For the full test which also tests the optimization mode please run:

`cd <path_to_genGEO> ; ./runTests.sh full`


Contact
---
For contact please visit: [https://geg.ethz.ch/people/](https://geg.ethz.ch/people/)
