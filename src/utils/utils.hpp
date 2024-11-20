/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <vector>
#include <string>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>


/**
 * @brief Operator to multiply a scalar by a vector.
 * 
 * @param scalar 
 * @param vec 
 * @return std::vector<double> 
 */
inline std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    std::vector<double> result(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
}

/**
 * @brief Operator to multiply a vector by a scalar.
 * 
 * @param vec 
 * @param scalar 
 * @return std::vector<double> 
 */
inline std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    return scalar * vec;
}

/**
 * @brief Function to generate the initial velocity field for the lid driven cavity.
 * 
 * @param val u_lid value
 * @param nx 
 * @param ny 
 * @param filename_nodes file containing the nodes (output of the python script)
 */
std::vector<double> lid_driven(double val, unsigned int nx, unsigned int ny, const std::string& filename_nodes);

std::vector<double> uniform_left_inlet(double val, unsigned int nx, unsigned int ny, const std::string& filename_nodes);




#endif // __UTILS_HPP__