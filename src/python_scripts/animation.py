import os
import json 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# TODO: figure out how add a filter for wall and obstacle cells

def read_file_to_array(filename, nx, ny):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    t = len(lines)
    data_array = np.zeros((nx, ny, t))
    
    for i, line in enumerate(lines):
        values = list(map(float, line.split()))
        data_array[:, :, i] = np.array(values).reshape((nx, ny))
    
    return data_array


def animate_array(array, title):
    output_folder = "output_animations"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output_file = os.path.join(output_folder, title + '.mp4')
    max_val = np.max(array)
    min_val = min(0,np.min(array))
    
    def update(frame):
        plt.clf()
        plt.imshow(array[:, :, frame], origin='upper', cmap='RdBu_r', interpolation='spline16', vmin=min_val, vmax=max_val)
        plt.colorbar()
        plt.title(title)
    
    
    ani = FuncAnimation(plt.figure(), update, frames=array.shape[2], interval=100)
    
    # Save the animation to a file
    ani.save(output_file, writer='ffmpeg')
    print(f"Animation saved in {output_file}")

def read_params():
    # Read the JSON file
    with open('params.json', 'r') as file:
        params = json.load(file)
    
    # Extract the variables
    image_path = params['image_path']
    num_points_x = params["lattice"]['new_nx']
    num_points_y = params["lattice"]['new_ny']
    
    return num_points_x, num_points_y

def main(directory, nx, ny):

    velocity = read_file_to_array(directory + '/velocity_out.txt', nx, ny)
    rho = read_file_to_array(directory + '/rho_out.txt', nx, ny)

    animate_array(velocity, 'Velocity Animation')
    animate_array(rho, 'Rho Animation')

if __name__ == "__main__":
    directory = "output_results" # input("Enter the directory containing the CSV files: ")
    nx, ny = read_params() 
    
    main(directory, nx, ny)