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
             unsigned int nx, unsigned int ny, double Cx, double Ct, double Crho) {
  // Save ux_out
    double Cvel = Cx / Ct;

    for(unsigned int index = 0; index < nx * ny; ++index){
        file_u << std::sqrt(ux_out[index] * ux_out[index] + uy_out[index] * uy_out[index]) * Cvel<< " ";
        file_ux << ux_out[index] * Cvel << " ";
        file_uy << uy_out[index] * Cvel << " ";
        file_rho << rho_out[index] * Crho << " ";
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
left_inlet(double val, unsigned int nx, unsigned int ny, const std::string& filename_nodes)
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

Filter::Filter() : frequency(1.0) {
    filter_function = unitary_filter;
}

Filter::Filter(const std::string& filter_name, double arg) : frequency(arg) {
    set_filter_function(filter_name, arg);
}

void Filter::set_filter_function(const std::string& filter_name, double arg) {
    if (filter_name == "sigmoid") {
        filter_function = sigmoid_filter;
    } else if (filter_name == "sinusoidal") {
        this->frequency = arg;
        filter_function = [this](double value) { return sinusoidal_filter(value, this->frequency); };
    } else {
        filter_function = unitary_filter;
    }
}

double Filter::apply(double value) const {
    return filter_function(value);
}


/// @name User defined filter functions, must input values between [0,1] and output values between [0,1]
///@{
double Filter::unitary_filter(double value) {
    return 1.0;
}

double Filter::sigmoid_filter(double value) {
    return 1 / (1 + std::exp(-25 * (value - 0.2)));
}

double Filter::sinusoidal_filter(double value, double frequency) {
    return 0.5 * (1 + std::sin(2 * M_PI * frequency * value ));
}
/// }@
