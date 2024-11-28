import cv2
import numpy as np
import csv
import json 

def read_params():
    # Read the JSON file
    with open('params.json', 'r') as file:
        params = json.load(file)
    
    # Extract the variables
    image_path = params['image_path']
    num_points_x = params["lattice"]['nx']
    num_points_y = params["lattice"]['ny']
    
    return image_path, num_points_x, num_points_y

def adapt_nx_ny(image_path, nx, ny):
    image = cv2.imread(image_path)
    if image is None:
        raise ValueError("Image not found or unable to load.")
    
    # Get image dimensions
    height, width, _ = image.shape
    height = height - 1
    width = width - 1

    # Calculate the spacing between points
    x_spacing = width // (nx - 1)
    y_spacing = height // (ny - 1)

    nx = width // x_spacing + 1
    ny = height // y_spacing + 1
    
    print(f"NEW nx: {nx}, ny: {ny}")

    # Read the JSON file
    with open('params.json', 'r') as file:
        params = json.load(file)
    
    # Update the JSON file with new values
    params["generated_variables"]['new_nx'] = nx
    params["generated_variables"]['new_ny'] = ny

    # Calculate the time interval
    length = params["lattice"]['Length']
    dt  = np.sqrt(3) *  length/ np.sqrt(nx * nx + ny * ny)
    iterations = int(np.ceil(params["time"]['T_final'] / dt))
    params["generated_variables"]['dt'] = dt
    params["generated_variables"]['iterations'] = iterations
    
    # Write the updated JSON file
    with open('params.json', 'w') as file:
        json.dump(params, file, indent=4)
    
    return nx, ny

def preprocess_image(image_path, cutoff_value=128):

    # Split the original image into its RGB components and the grayscale image
    original_image = cv2.imread(image_path)
    if original_image is None:
        raise ValueError("Image not found or unable to load.")
    b, g, r = cv2.split(original_image)
    # gray = cv2.cvtColor(original_image, cv2.COLOR_BGR2GRAY)

    _, r = cv2.threshold(r, cutoff_value, 255, cv2.THRESH_BINARY)
    _, g = cv2.threshold(g, cutoff_value, 255, cv2.THRESH_BINARY)
    _, b = cv2.threshold(b, cutoff_value, 255, cv2.THRESH_BINARY)
    # _, gray = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY)
    gray = (r + g + b) 

    # cv2.imwrite('r.png', r)
    # cv2.imwrite('g.png', g)
    # cv2.imwrite('b.png', b)
    # cv2.imwrite('gray.png', gray)

    return r, g, b, gray

def classify_points(image_path, num_points_x, num_points_y):
    
    # Preprocess the image to r,g,b adn grayscale each corresponding to wall, inlet, outlet and fluid pixels
    r,g,b,gray = preprocess_image(image_path)
    
    # Get image dimensions
    height, width = r.shape
    height = height - 1
    width = width - 1

    # Calculate the spacing between points 
    x_spacing = width // (num_points_x - 1)
    y_spacing = height // (num_points_y - 1)

    # Initialize lists to hold fluid, solid, inlet and outlet points
    fluid_points = []
    obstacle_points = []
    solid_points = []
    inlet_points = []
    outlet_points = []
    external_points = []

    # Classify points
    for i in range(num_points_y):
        for j in range(num_points_x):
            x_px = int(j * x_spacing) 
            y_px = int(i * y_spacing) 
            if (r[y_px, x_px] == 0 and g[y_px, x_px] == 0 and b[y_px, x_px] == 255):  # Wall pixel
                solid_points.append((j, i))
                external_points.append((j, i))
            if (r[y_px, x_px] == 0 and g[y_px, x_px] == 255 and b[y_px, x_px] == 0):  # inlet pixel
                inlet_points.append((j, i))
                external_points.append((j, i))
            if (r[y_px, x_px] == 255 and g[y_px, x_px] == 0 and b[y_px, x_px] == 0):  # outlet pixel
                outlet_points.append((j, i))
                external_points.append((j, i))
            if (r[y_px, x_px] == 255 and g[y_px, x_px] == 255 and b[y_px, x_px] == 255): # obstacle pixel 
                obstacle_points.append((j, i))
                external_points.append((j, i))
            if (r[y_px, x_px] == 0 and g[y_px, x_px] == 0 and b[y_px, x_px] == 0):  # Fluid pixel
                fluid_points.append((j, i))

    return fluid_points, obstacle_points, solid_points, inlet_points, outlet_points, external_points, gray

def relative_cutoff_distance(image, coord1, coord2):
    x1 = coord1[0]
    y1 = coord1[1]
    x2 = coord2[0] 
    y2 = coord2[1] 

    if image[y1, x1] == image[y2, x2]:
        raise ValueError("Both coordinates are in the same region (both black or both white).")

    # Calculate the Euclidean distance between the two points
    total_distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    # Traverse the line between the two points to find the cutoff point
    num_steps = int(80*total_distance)
    for step in range(num_steps + 1):
        t = step / num_steps
        x = int(x1 + t * (x2 - x1))
        y = int(y1 + t * (y2 - y1))
        if image[y, x] != image[y1, x1]:
            # cutoff_distance = np.sqrt((x - x1) ** 2 + (y - y1) ** 2)
            # print(f"Distance between points: ({x1}, {y1}) and ({x2}, {y2}) at ({x}, {y}): t = {t} \n\n")
            return  t #cutoff_distance / total_distance
        

    raise ValueError("No color change detected between the two coordinates.")

def remove_boundary_points_from_fluid(fluid_points_set, boundary_points_distances):
    for boundary_point in boundary_points_distances:
        point = (boundary_point[0], boundary_point[1])
        if point in fluid_points_set:
            fluid_points_set.remove(point)
    return list(fluid_points_set)

def identify_boundary_points_and_distances(fluid_points, external_points, num_points_x, num_points_y, image):
    boundary_points_distances = []

    # Define the 8 possible directions (E, N, W, S, NE, NW, SW, SE)
    directions = [(1,0), (0,-1), (-1,0), (0,1), (1,-1), (-1,-1), (-1,1), (1,1)]

    # Get image dimensions
    height, width = image.shape
    height = height - 1
    width = width - 1

    # Calculate the spacing between points (remove //)
    x_spacing = width // (num_points_x - 1)
    y_spacing = height // (num_points_y - 1)

    # Convert fluid_points to a set for faster lookup
    fluid_points_set = set(fluid_points)

    for ext_point in external_points:
        x, y = ext_point
        for (dx, dy) in directions :
            distance = 0.0

            # Check if the adjacent point is fluid
            x_adj, y_adj = x + dx, y + dy
            if (x_adj, y_adj) in fluid_points_set:
                #Convert to pixel coordinates
                x_px = int(x * x_spacing) 
                y_px = int(y * y_spacing) 
                x_adj_px = int(x_adj * x_spacing) 
                y_adj_px = int(y_adj * y_spacing) 

                #fluid_points_set.remove((adj_x, adj_y))
                distance = relative_cutoff_distance(image, (x_adj_px,y_adj_px), (x_px,y_px))
                boundary_points_distances.append((x_adj, y_adj, -dx, -dy, distance))

    # Remove boundary points from fluid points
    fluid_points = remove_boundary_points_from_fluid(fluid_points_set, boundary_points_distances)

    return fluid_points, boundary_points_distances

def create_csv_with_point_types_and_distances(fluid_points, obstacle_points, solid_points, inlet_points, outlet_points, boundary_points_distances, num_points_x, num_points_y, output_csv_path):
    # Create a dictionary to store the type and distances for each point
    point_data = {}

    # Mark fluid points
    for x, y in fluid_points:
        point_data[(x, y)] = [0] + [0.0] * 8

    # Mark obstacle points
    for x, y in obstacle_points:
        point_data[(x, y)] = [1] + [0.0] * 8

    # Mark solid points
    for x, y in solid_points:
        point_data[(x, y)] = [2] + [0.0] * 8
    
    # Mark inlet points
    for x, y in inlet_points:
        point_data[(x, y)] = [3] + [0.0] * 8

    # Mark outlet points
    for x, y in outlet_points:
        point_data[(x, y)] = [4] + [0.0] * 8

    # Mark boundary points and store distances
    directions = [(1,0), (0,-1), (-1,0), (0,1), (1,-1), (-1,-1), (-1,1), (1,1)]

    for x, y, dx, dy, distance in boundary_points_distances:
        if (x, y) not in point_data:
            point_data[(x, y)] = [5] + [0.0] * 8
        direction_index = directions.index((dx, dy))
        point_data[(x, y)][1 + direction_index] = distance

    # Write the data to a CSV file
    with open(output_csv_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['coord_x', 'coord_y', 'type', 'E', 'N', 'W', 'S', 'NE', 'NW', 'SW', 'SE'])
        for i in range(num_points_y):
            for j in range(num_points_x):
                if (j, i) in point_data:
                    writer.writerow([j, i] + point_data[(j, i)])
                else:
                    writer.writerow([j, i, 1] + [0.0] * 8)

def draw_lattice(image_path, nx, ny, output_image_path):
    # Load the image
    image = cv2.imread(image_path)
    if image is None:
        raise ValueError("Image not found or unable to load.")
    
    # Get image dimensions
    height, width, _ = image.shape
    height = height - 1
    width = width - 1

    # Calculate the spacing between points (remove //)
    x_spacing = width // (nx - 1)
    y_spacing = height // (ny - 1)

    # Define the size of the purple squares
    square_size = 1

    # Draw purple squares on the grid
    for i in range(ny):
        for j in range(nx):
            x = int(j * x_spacing)
            y = int(i * y_spacing)
            top_left = (x - square_size // 2, y - square_size // 2)
            bottom_right = (x + square_size // 2, y + square_size // 2)
            cv2.rectangle(image, top_left, bottom_right, (255, 0, 255), -1)  # Purple color in BGR

    # Save the modified image
    cv2.imwrite(output_image_path, image)




def main():
    image_path, num_points_x, num_points_y = read_params()
    # Adapt nx and ny to the image size
    num_points_x, num_points_y = adapt_nx_ny(image_path, num_points_x, num_points_y)

    # Classify points to fluid, obstacle, solid, inlet, outlet and external (which is the sum of the previous three)
    fluid_points, obstacle_points, solid_points, inlet_points, outlet_points, external_points, original_image = classify_points(image_path, num_points_x, num_points_y)
    
    # Identify boundary points and calculate distances wrt external points
    fluid_points, boundary_points_distances = identify_boundary_points_and_distances(fluid_points, external_points, num_points_x, num_points_y, original_image)
    
    # Create a CSV file with the lattice information
    create_csv_with_point_types_and_distances(fluid_points, obstacle_points, solid_points, inlet_points, outlet_points, boundary_points_distances, num_points_x, num_points_y, 'lattice.csv')
    
    # Create an image overlaid with the lattice
    draw_lattice(image_path, num_points_x, num_points_y, 'Lattice_nodes.png')
    
    return num_points_x, num_points_y

if __name__ == "__main__":
    #image_path = 'lid_driven.png'
    #num_points_x = 50
    #num_points_y = 50
    nx, ny = main()
    print(f"CSV file created successfully with nx: {nx}, ny: {ny}")


