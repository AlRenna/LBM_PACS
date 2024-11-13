/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */


#include <cstdlib>
#include <iostream>

inline void run_terminal_python(std::string script_path) {
    std::string command = {"python3 " + script_path};
    int result = system(command.c_str());
    if (result != 0) {
        std::cerr << "Pyhton terminal command failed with exit code: " << result << std::endl;
    }
}

