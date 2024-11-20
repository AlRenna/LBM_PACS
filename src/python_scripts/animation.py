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

def read_lattice_map(filepath, nx, ny):
    lattice_map = np.zeros((nx, ny), dtype=bool)
    df = pd.read_csv(filepath)
    for _, row in df.iterrows():
        if row['type'] == 4 or row['type'] == 0:
            # TODO: when passing different val
            lattice_map[int(row['coord_x']), int(row['coord_y'])] = True 
    return lattice_map

def animate_array(array, dt, title, lattice_map):
    output_folder = "output_animations"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output_file = os.path.join(output_folder, title + '.mp4')
    max_val = np.max(array)
    min_val = min(0, np.min(array))
    
    def update(frame):
        plt.clf()
        frame_data = np.where(lattice_map, array[:, :, frame], np.nan)
        plt.imshow(frame_data, origin='upper', cmap='RdBu_r', interpolation='spline16', vmin=min_val, vmax=max_val)
        plt.colorbar()
        plt.title(f"{title} - Time: {frame * dt :.2f}")
    
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
    num_points_x = params["generated_variables"]['new_nx']
    num_points_y = params["generated_variables"]['new_ny']
    dt = params["generated_variables"]['dt']
    save_iter = params["time"]['save_iter']
    dt = dt * save_iter
    
    return num_points_x, num_points_y, dt

def main(directory, nx, ny, dt):
    lattice_map = read_lattice_map('lattice.csv', nx, ny)
    velocity = read_file_to_array(os.path.join(directory, 'velocity_out.txt'), nx, ny)
    rho = read_file_to_array(os.path.join(directory, 'rho_out.txt'), nx, ny)

    animate_array(velocity, dt, 'Velocity Animation', lattice_map)
    animate_array(rho, dt, 'Rho Animation', lattice_map)

if __name__ == "__main__":
    directory = "output_results" # input("Enter the directory containing the CSV files: ")
    nx, ny, dt, = read_params() 
    # print(nx, ny, dt)
    main(directory, nx, ny, dt)