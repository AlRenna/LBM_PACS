import os
import json 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def read_file_to_array(filename, nx, ny):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    t = len(lines)
    data_array = np.zeros((ny, nx, t))
    
    for i, line in enumerate(lines):
        values = list(map(float, line.split()))
        data_array[:, :, i] = np.array(values).reshape((ny, nx))
    return data_array

def read_lattice_map(filepath, nx, ny):
    lattice_map = np.zeros((ny, nx), dtype=bool)
    df = pd.read_csv(filepath)
    for _, row in df.iterrows():
        if row['type'] == 4 or row['type'] == 0:
            lattice_map[int(row['coord_y']), int(row['coord_x'])] = True 
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


def animate_velocity_field(v_array, ux_array, uy_array, dt, title, lattice_map):
    output_folder = "output_animations"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # u_array = np.flip(u_array, axis=0)
    # v_array = np.flip(v_array, axis=0)
    
    output_file = os.path.join(output_folder, title + '.mp4')
    X, Y = np.meshgrid(np.arange(ux_array.shape[1]), np.arange(ux_array.shape[0]))


    def update(frame):
        plt.clf()
        u_frame = np.where(lattice_map, ux_array[:, :, frame], np.nan)
        v_frame = np.where(lattice_map, uy_array[:, :, frame], np.nan)
        vel_mag = np.where(lattice_map, v_array[:, :, frame], np.nan)
        u_frame = np.nan_to_num(u_frame, nan=0.0) + 1e-10  
        v_frame = np.nan_to_num(v_frame, nan=0.0) + 1e-10  
        plt.imshow(vel_mag, origin='upper', cmap='RdBu_r', interpolation='spline16', vmin=0, vmax=1)
        plt.quiver(X, Y, u_frame, v_frame, vel_mag, cmap='viridis')
        plt.title(f"{title} - Time: {frame * dt :.2f}")
    
    ani = FuncAnimation(plt.figure(), update, frames=ux_array.shape[2], interval=100)
    
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
    T_final = params["time"]['T_final']
    max_iter = params["generated_variables"]['iterations']
    dt = dt * save_iter
    
    return num_points_x, num_points_y, dt, T_final, max_iter

def plot_lift_drag(filename, T_final, max_iter, dt):
    output_folder = "output_animations"
    time = np.array(range(0, max_iter))
    time = time * dt
    output_file = os.path.join(output_folder, 'lift_drag.png')
    
    if not os.path.exists(filename):
        return
    
    with open(filename, 'r') as file:
        lines = file.readlines()
    

    lift = list(map(float, lines[1].split()))
    lift = np.array(lift)
    drag = list(map(float, lines[3].split()))
    drag = np.array(drag)
    
    plt.clf()
    plt.plot(time , lift, label='Lift')
    plt.plot(time, drag, label='Drag')
    plt.xlabel('Time')
    plt.ylabel('Force')
    plt.title('Lift and Drag vs Time')
    plt.legend()
    # Save the plot to a file
    plt.savefig(os.path.join(output_folder, 'lift_drag.png'), format='png')
    print(f"Plot saved in {output_file}")

def main(directory, nx, ny, dt, T_final, max_iter):
    lattice_map = read_lattice_map('lattice.csv', nx, ny)
    velocity = read_file_to_array(os.path.join(directory, 'velocity_out.txt'), nx, ny)
    rho = read_file_to_array(os.path.join(directory, 'rho_out.txt'), nx, ny)
    u_array = read_file_to_array(os.path.join(directory, 'ux_out.txt'), nx, ny)
    v_array = read_file_to_array(os.path.join(directory, 'uy_out.txt'), nx, ny)

    animate_velocity_field(velocity, u_array, v_array, dt, 'Vectorial Animation', lattice_map)
    animate_array(velocity, dt, 'Velocity Animation', lattice_map)
    animate_array(rho, dt, 'Rho Animation', lattice_map)
    plot_lift_drag(os.path.join(directory, 'lift_&_drag.txt'), T_final, max_iter, dt)


if __name__ == "__main__":
    directory = "output_results" # input("Enter the directory containing the CSV files: ")
    nx, ny, dt, T_final, max_iter = read_params() 
    # print(nx, ny, dt)
    main(directory, nx, ny, dt, T_final, max_iter)