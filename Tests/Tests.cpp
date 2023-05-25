#include <iostream>
// Catch2 header
#include "catch_amalgamated.hpp"

#include "AMReX_REAL.H"

int main(int argc, char *argv[])
{
    constexpr int cout_precision          = 17;
    Catch::StringMaker<double>::precision = cout_precision;
    Catch::StringMaker<float>::precision  = cout_precision;

    int result = Catch::Session().run(argc, argv);

    return result;
}