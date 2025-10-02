# GRTeclyn

GRTeclyn (previously referred to as GRAMReX) is a new numerical relativity code developed by the [GRTL Collaboration](www.grtlcollaboration.org) that is currently still under development.  It is a port of the [GRChombo code](https://github.com/GRChombo/GRChombo) (based on the Chombo libraries) to the [AMReX](https://amrex-codes.github.io/) library in order to leverage AMReX's support for GPUs and ongoing active development.

The AMReX documentation can be found [here](https://amrex-codes.github.io/amrex/docs_html).

The name follows a similar pattern to GRChombo, namely "GR" + "\<Tool in a foreign language\>". In this case, "Teclyn" is a Welsh word for "Tool".

## Development status

The current timeline for releasing GRTeclyn is December 2025 - the goal is that at this point it should contain well tested methods for vacuum binaries with gravitational wave extraction, and capabilities to add scalar fields. Further features will be ported later.

A detailed list of the status of particular features can be found in the [documentation](https://grtlcollaboration.github.io/GRTeclyn/#summary-of-features).

## Documentation

Documentation can be found [here](https://grtlcollaboration.github.io/GRTeclyn/) (under construction). Note that the GitHub wiki is no longer in use.

The documentation contains useful information on obtaining and building the code, prerequisities and running the binary black hole example.

## License

GRTeclyn is licensed under the BSD 3-Clause License. Please see [LICENSE](LICENSE) for details.