/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include <cstdlib>
#include <iostream>
#include <vector>

inline void run_terminal_command(const std::string& command) {
    int result = system(command.c_str());
    if (result != 0) {
        std::cerr << "Command failed with exit code: " << result << std::endl;
    }
}

inline void run_terminal_python(){
    const std::vector<std::string> commands_ = {"python3 ../../src/python_scripts/lattice_generation_RGB.py"};
    for (const auto &command : commands_){
        run_terminal_command(command);
    }
}

