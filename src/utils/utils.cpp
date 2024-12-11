/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "src/utils/utils.hpp"
#include <cmath>

void 
writeResults(std::ofstream &file_u, std::ofstream &file_ux, std::ofstream &file_uy, std::ofstream &file_rho, 
             const std::vector<double>& ux_out, const std::vector<double>& uy_out, const std::vector<double>& rho_out, 
             unsigned int nx, unsigned int ny) {
  // Save ux_out
  for(unsigned int index = 0; index < nx * ny; ++index){
    file_u << std::sqrt(ux_out[index] * ux_out[index] + uy_out[index] * uy_out[index]) << " ";
    file_ux << ux_out[index] << " ";
    file_uy << uy_out[index] << " ";
    file_rho << rho_out[index] << " ";
  }
  file_u << "\n";
  file_ux << "\n";
  file_uy << "\n";
  file_rho << "\n";
}


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

        if(type == 2){
            ux_in[y*nx+x] = val;
        }
        else if (type == 5){
            for (int i = y*nx; i < y*nx+x; ++i) {
                ux_in[i] = 0.0;
            }
            break;
        }
    }
    file.close();
    return ux_in;
}

std::vector<double>
uniform_left_inlet(double val, unsigned int nx, unsigned int ny, const std::string& filename_nodes)
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

        if(type == 3){
            ux_in[y*nx+x] = val;
        }
    }
    file.close();
    return ux_in;
}

