import os
import json 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def read_txt_files(directory):
    print("Sorting files ...........")
    ux_files = sorted([f for f in os.listdir(directory) if f.startswith('ux_out_') and f.endswith('.txt')])
    uy_files = sorted([f for f in os.listdir(directory) if f.startswith('uy_out_') and f.endswith('.txt')])
    rho_files = sorted([f for f in os.listdir(directory) if f.startswith('rho_out_') and f.endswith('.txt')])
    print("Loading files ...........")
    print(ux_files)
    ux_data = [np.loadtxt(os.path.join(directory, f), delimiter=',') for f in ux_files]
    uy_data = [np.loadtxt(os.path.join(directory, f), delimiter=',') for f in uy_files]
    rho_data = [np.loadtxt(os.path.join(directory, f), delimiter=',') for f in rho_files]
    
    return ux_data, uy_data, rho_data

# def create_3d_arrays(ux_data, uy_data, rho_data, nx, ny):
#     T = len(ux_data)
#     velocity = np.zeros((nx, ny, T, 2))  # 4D array to store velocity components
#     rho = np.zeros((nx, ny, T))
    
#     for t in range(T):
#         velocity[:, :, t, 0] = ux_data[t].values.reshape((nx, ny))
#         velocity[:, :, t, 1] = uy_data[t].values.reshape((nx, ny))
#         rho[:, :, t] = rho_data[t].values.reshape((nx, ny))
    
#     return velocity, rho

def create_3d_arrays(ux_data, uy_data, rho_data, nx, ny):
    T = len(ux_data)
    velocity = np.zeros((nx, ny, T, 2))  # 4D array to store velocity components
    velocity_mag = np.zeros((nx, ny, T))  
    rho = np.zeros((nx, ny, T))
    
    for t in range(T):
        velocity[:, :, t, 0] = ux_data[t].reshape((nx, ny))
        velocity[:, :, t, 1] = uy_data[t].reshape((nx, ny))
        velocity_mag[:, :, t] = np.sqrt(ux_data[t].reshape((nx, ny))**2 +  uy_data[t].reshape((nx, ny))**2)
        rho[:, :, t] = rho_data[t].reshape((nx, ny))
    
    return velocity, velocity_mag, rho

# def animate_array(array, title):
#     fig, ax = plt.subplots()
#     cax = ax.imshow(array[:, :, 0], cmap='viridis')
#     fig.colorbar(cax)
#     ax.set_title(title)
    
#     def update(frame):
#         cax.set_array(array[:, :, frame])
#         return cax,
    
#     ani = FuncAnimation(fig, update, frames=array.shape[2], blit=True)
#     plt.show()

def animate_array(array, title):
    output_folder = "output_animations"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output_file = os.path.join(output_folder, title + '.mp4')
    fig, ax = plt.subplots()
    cax = ax.imshow(array[:, :, 0], cmap='viridis')
    fig.colorbar(cax)
    ax.set_title(title)
    
    def update(frame):
        cax.set_array(array[:, :, frame])
        return cax,
    
    ani = FuncAnimation(fig, update, frames=array.shape[2], blit=True)
    
    # Save the animation to a file
    ani.save(output_file, writer='ffmpeg')
    print(f"Animation saved in {output_file}")

def read_params():
    # Read the JSON file
    with open('params.json', 'r') as file:
        params = json.load(file)
    
    # Extract the variables
    image_path = params['image_path']
    num_points_x = params['new_nx']
    num_points_y = params['new_ny']
    
    return num_points_x, num_points_y

def main(directory, nx, ny):
    ux_data, uy_data, rho_data = read_txt_files(directory)
    velocity, velocity_mag, rho = create_3d_arrays(ux_data, uy_data, rho_data, nx, ny)
    
    animate_array(velocity[:, :, :, 0], 'Ux Animation')
    animate_array(velocity[:, :, :, 1], 'Uy Animation')
    animate_array(velocity_mag, 'Velocity Animation')
    animate_array(rho, 'Rho Animation')

if __name__ == "__main__":
    directory = "output_results" # input("Enter the directory containing the CSV files: ")
    nx, ny = read_params() 
    
    main(directory, nx, ny)