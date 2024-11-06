import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def read_csv_files(directory):
    ux_files = sorted([f for f in os.listdir(directory) if f.startswith('ux_out_') and f.endswith('.csv')])
    uy_files = sorted([f for f in os.listdir(directory) if f.startswith('uy_out_') and f.endswith('.csv')])
    rho_files = sorted([f for f in os.listdir(directory) if f.startswith('rho_out_') and f.endswith('.csv')])
    
    ux_data = [pd.read_csv(os.path.join(directory, f), header=None) for f in ux_files]
    uy_data = [pd.read_csv(os.path.join(directory, f), header=None) for f in uy_files]
    rho_data = [pd.read_csv(os.path.join(directory, f), header=None) for f in rho_files]
    
    return ux_data, uy_data, rho_data

def create_3d_arrays(ux_data, uy_data, rho_data, nx, ny):
    T = len(ux_data)
    velocity = np.zeros((nx, ny, T, 2))  # 4D array to store velocity components
    rho = np.zeros((nx, ny, T))
    
    for t in range(T):
        velocity[:, :, t, 0] = ux_data[t].values.reshape((nx, ny))
        velocity[:, :, t, 1] = uy_data[t].values.reshape((nx, ny))
        rho[:, :, t] = rho_data[t].values.reshape((nx, ny))
    
    return velocity, rho

def animate_array(array, title):
    fig, ax = plt.subplots()
    cax = ax.imshow(array[:, :, 0], cmap='viridis')
    fig.colorbar(cax)
    ax.set_title(title)
    
    def update(frame):
        cax.set_array(array[:, :, frame])
        return cax,
    
    ani = FuncAnimation(fig, update, frames=array.shape[2], blit=True)
    plt.show()

def main(directory, nx, ny):
    ux_data, uy_data, rho_data = read_csv_files(directory)
    velocity, rho = create_3d_arrays(ux_data, uy_data, rho_data, nx, ny)
    
    animate_array(velocity[:, :, :, 0], 'Ux Animation')
    animate_array(velocity[:, :, :, 1], 'Uy Animation')
    animate_array(rho, 'Rho Animation')

if __name__ == "__main__":
    directory = input("Enter the directory containing the CSV files: ")
    nx = int(input("Enter the value of nx: "))
    ny = int(input("Enter the value of ny: "))
    main(directory, nx, ny)