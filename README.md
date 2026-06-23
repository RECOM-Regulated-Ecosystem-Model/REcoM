# Regulated Ecosystem Model (REcoM)

The Regulated Ecosystem Model (REcoM) is an ocean biogeochemistry model that represents the plankton ecosystem with medium complexity, including variable stoichiometry and multiple plankton groups.

## Table of contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Quick start](#quick-start)
* [Usage](#usage)
* [Known issues and limitations](#known-issues-and-limitations)
* [Getting help](#getting-help)
* [Contributing](#contributing)
* [License](#license)
* [Acknowledgments](#acknowledgments)


## Introduction

The Regulated Ecosystem Model (REcoM) is a marine biogeochemistry model that simulates the coupled cycles of carbon, nitrogen, silicic acid, iron and oxygen in the ocean. It represents the plankton ecosystem with medium complexity, including variable stoichiometry and multiple plankton groups. The user can choose between configurations with varying complexities. REcoM uses mocsy to solve the carbonate chemistry and compute air-sea CO2 and O2 fluxes. It is used for studies of plankton ecosystems, primary production, the global carbon budget, and carbon-climate feedbacks spanning time-scales from recent decades, future centuries and paleo timescales. It is currently used in its version coupled to FESOM and AWI-ESM.

## Installation

## Quick start

## Usage

### Basic operation

### More options

## Known issues and limitations

## Getting help

## Contributing

In order to contribute to REcoM, please start by forking the repository and creating a branch off of `master`. After you have developed your work, please create a pull request with a clear description of the changes you have made. The pull request will then undergo automated CI testing and review by the maintainers. The CI tests check two main aspects. The first is whether the code still builds correctly with your introduced changes. The second is a number of code formatting and quality aspects, evaluated by Codee and Fortitude through pre-commit hooks.

We highly recommend that you install [pre-commit](https://pre-commit.com/) in your local environment and that you enable it for your clone of the repository. Then before every commit, the pre-commit hooks will automatically check the formatting and quality of your code, and will prevent you from committing if any issues are found. This way, you can ensure that your code adheres to the standards of the project and that it will pass the CI tests when you create a pull request.

For convenience you can also use [act](https://github.com/nektos/act) to run the CI tests locally before pushing your changes to GitHub. This will allow you to catch any issues early on and avoid delays in the review process.

You can learn more about the described tools here:

* [Pre-commit](https://pre-commit.com/)
* [Codee](https://codee.com)
* [Fortitude](https://fortitude.readthedocs.io/en/stable/)
* [ACT](https://github.com/nektos/act)

## License

## Acknowledgments

This README was based on the [readmine](https://mhucka.github.io/readmine/) template.
