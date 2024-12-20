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
#include <cmath>


/// @name Utility functions
/// @{

/**
 * @brief Function to write the results of the simulation to a file.
 * 
 * @param file_u Output file for velocity magnitude
 * @param file_ux Output file for x-component of velocity
 * @param file_uy Output file for y-component of velocity
 * @param file_rho Output file for density
 * @param ux_out Vector containing x-component of velocity
 * @param uy_out Vector containing y-component of velocity
 * @param rho_out Vector containing density
 * @param nx Number of nodes in x-direction
 * @param ny Number of nodes in y-direction
 */
void writeResults(std::ofstream &file_u, std::ofstream &file_ux, std::ofstream &file_uy, std::ofstream &file_rho, 
                  const std::vector<double>& ux_out, const std::vector<double>& uy_out, const std::vector<double>& rho_out, 
                  unsigned int nx, unsigned int ny, double Cx, double Ct, double Crho);


/**
 * @brief Function to convert a double array to a std::vector<double>.
 * 
 * @param array Pointer to the double array
 * @param size Size of the array
 */
inline std::vector<double> arrayToVector(const double* array, std::size_t size) {
    return std::vector<double>(array, array + size);
}

/**
 * @brief Function to convert a std::vector<double> to a double array.
 * 
 * @param vec The input vector
 */
inline double* vectorToArray(const std::vector<double>& vec) {
    double* array = new double[vec.size()];
    std::copy(vec.begin(), vec.end(), array);
    return array;
}

/**
 * @brief Function to convert a std::vector<int> to an int array.
 * 
 * @param vec The input vector
 */
inline int* vectorToArray(const std::vector<int>& vec) {
    int* array = new int[vec.size()];
    std::copy(vec.begin(), vec.end(), array);
    return array;
}

/**
 * @brief Function to convert a std::vector<std::vector<double>> to a double array.
 * 
 * @param vec The input 2D vector
 */
inline double* vector2DToArray(const std::vector<std::vector<double>>& vec) {
    std::size_t total_size = 0;
    for (const auto& inner_vec : vec) {
        total_size += inner_vec.size();
    }
    
    double* array = new double[total_size];
    std::size_t index = 0;
    for (const auto& inner_vec : vec) {
        std::copy(inner_vec.begin(), inner_vec.end(), array + index);
        index += inner_vec.size();
    }
    
    return array;
}

/// @}

/// @name Operators
/// @{

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

/// @}

/// @name Tests utilities functions
/// @{

/**
 * @brief Function to generate the initial velocity field for the lid driven cavity.
 * 
 * @param val u_lid value
 * @param nx 
 * @param ny 
 * @param filename_nodes file containing the nodes (output of the python script)
 */
std::vector<double> lid_driven(double val, unsigned int nx, unsigned int ny, const std::string& filename_nodes);

/**
 * @brief Function to generate the initial velocity field for the uniform left inlet case.
 * 
 * @param val u_in value
 * @param nx 
 * @param ny 
 * @param filename_nodes 
 */
std::vector<double> uniform_left_inlet(double val, unsigned int nx, unsigned int ny, const std::string& filename_nodes);

/// @}


#endif // __UTILS_HPP__