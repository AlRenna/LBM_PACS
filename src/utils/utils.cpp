/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "src/utils/utils.hpp"

std::vector<double> 
lid_driven(double val, unsigned int nx, unsigned int ny, const std::string& filename_nodes)
{
    std::vector<double> ux_in(nx*ny, 0.0);

    std::ifstream file(filename_nodes);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file");
    }

    std::string line;
    // Skip the first line
    std::getline(file, line);

    while (std::getline(file, line)) {
        if (line.empty()) {
        break;
        }
        std::stringstream ss(line);
        std::string token;

        // Read coordinates
        std::getline(ss, token, ',');
        unsigned int x = std::stoi(token);
        std::getline(ss, token, ',');
        unsigned int y = std::stoi(token);

        // Read NodeType
        std::getline(ss, token, ',');
        unsigned int type = std::stoi(token);

        if(type == 1){
            ux_in[y*nx+x] = val;
        }
        else if (type == 4){
            for (int i = y*nx; i < y*nx+x; ++i) {
                ux_in[i] = 0.0;
            }
            break;
        }
    }
    file.close();
    return ux_in;
}