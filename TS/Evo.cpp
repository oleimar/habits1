#include "EvoSim.hpp"
#include <iostream>

// The EvoLearn program runs evolutionary simulations
// Copyright (C) 2025  Olof Leimar
// See Readme.md in top repo directory for copyright notice

int main(int argc, char* argv[])
{
    // Open input file and read indata
    EvoInpData eid(argv[1]);
    if (!eid.OK) {
        std::cout << "Input failed!" << "\n";
        return -1;
    }
    // Run the iteration
    Evo evo(eid);
    evo.Run();
    return 0;
}
