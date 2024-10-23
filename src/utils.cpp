#include <cstdlib>
#include <iostream>
#include <vector>

void run_command(const std::string& command) {
    int result = system(command.c_str());
    if (result != 0) {
        std::cerr << "Command failed with exit code: " << result << std::endl;
    }
}

void run_python(const std::vector<std::string> &commands){
    const std::vector<std::string> commands_ = {"source env/bin/activate","python3 generate_lattice_RGB.py", "deactivate" };
    for (const auto &command : commands_){
        run_command(command);
    }
}
