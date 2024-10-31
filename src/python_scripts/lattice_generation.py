import cv2
import numpy as np
import csv

def preprocess_image_black_white(image_path, output_path, cutoff_value=128):
    # Load the image
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise ValueError("Image not found or unable to load.")
    # Convert the image to grayscale if it is not already
    if len(image.shape) == 3:
        image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    # Apply thresholding to convert grey pixels to black or white
    _, thresholded_image = cv2.threshold(image, cutoff_value, 255, cv2.THRESH_BINARY)

    # Save the thresholded image
    #cv2.imwrite(output_path, thresholded_image)

    return thresholded_image

def relative_cutoff_distance(image, coord1, coord2, x_spacing, y_spacing):
    x1 = coord1[0] * x_spacing
    y1 = coord1[1] * y_spacing
    x2 = coord2[0] * x_spacing
    y2 = coord2[1] * y_spacing

    if image[y1, x1] == image[y2, x2]:
        raise ValueError("Both coordinates are in the same region (both black or both white).")

    # Calculate the Euclidean distance between the two points
    total_distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    # Traverse the line between the two points to find the cutoff point
    num_steps = int(total_distance)
    for step in range(num_steps + 1):
        t = step / num_steps
        x = int(x1 + t * (x2 - x1))
        y = int(y1 + t * (y2 - y1))
        if image[y, x] != image[y1, x1]:
            cutoff_distance = np.sqrt((x - x1) ** 2 + (y - y1) ** 2)
            return cutoff_distance / total_distance

    raise ValueError("No color change detected between the two coordinates.")

def remove_boundary_points_from_internal(internal_points_set, boundary_points_distances):
    for boundary_point in boundary_points_distances:
        point = (boundary_point[0], boundary_point[1])
        if point in internal_points_set:
            internal_points_set.remove(point)
    return list(internal_points_set)

def identify_boundary_points_and_distances(internal_points, external_points, x_spacing, y_spacing, image):
    boundary_points_distances = []
    

    # Define the 8 possible directions (E, N, W, S, NE, NW, SW, SE)
    directions = [(1,0), (0,-1), (-1,0), (0,1), (1,-1), (-1,-1), (-1,1), (1,1)]

    # Convert internal_points to a set for faster lookup
    internal_points_set = set(internal_points)

    for ext_point in external_points:
        x, y = ext_point
        for (dx, dy) in directions :
            distance = 0.0
            adj_x, adj_y = x + dx, y + dy
            if (adj_x, adj_y) in internal_points_set:
                #internal_points_set.remove((adj_x, adj_y))
                distance = relative_cutoff_distance(image, (adj_x,adj_y), (x,y), x_spacing, y_spacing,)
                boundary_points_distances.append((adj_x, adj_y, -dx, -dy, distance))

    # Remove boundary points from internal points
    internal_points = remove_boundary_points_from_internal(internal_points_set, boundary_points_distances)

    return internal_points, boundary_points_distances

def classify_points(image_path, num_points_x, num_points_y):
    # Load the image
    # Preprocess the image to black and white
    
    preprocessed_image_path = 'preprocessed_image.png'
    preprocess_image_black_white(image_path, preprocessed_image_path)
    
    # Load the preprocessed image
    image = cv2.imread(preprocessed_image_path, cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise ValueError("Preprocessed image not found or unable to load.")

    # Get image dimensions
    height, width = image.shape

    # Calculate the spacing between points
    x_spacing = width // (num_points_x - 1)
    y_spacing = height // (num_points_y - 1)

    # Initialize lists to hold internal and external points
    internal_points = []
    external_points = []

    # Classify points
    for i in range(num_points_y):
        for j in range(num_points_x):
            x = j * x_spacing
            y = i * y_spacing
            if image[y, x] == 255:  # White pixel
                internal_points.append((j, i))
            else:  # Black pixel
                external_points.append((j, i))

    return internal_points, external_points, x_spacing, y_spacing, image

def create_csv_with_point_types_and_distances(internal_points, external_points, boundary_points_distances, num_points_x, num_points_y, output_csv_path):
    # Create a dictionary to store the type and distances for each point
    point_data = {}

    # Mark internal points
    for x, y in internal_points:
        point_data[(x, y)] = [0] + [0.0] * 8

    # Mark external points
    for x, y in external_points:
        point_data[(x, y)] = [1] + [0.0] * 8

    # Mark boundary points and store distances
    directions = [(1,0), (0,-1), (-1,0), (0,1), (1,-1), (-1,-1), (-1,1), (1,1)]

    for x, y, dx, dy, distance in boundary_points_distances:
        if (x, y) not in point_data:
            point_data[(x, y)] = [2] + [0.0] * 8
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

# Example usage
if __name__ == "__main__":
    image_path = 'image.png'
    num_points_x = 100
    num_points_y = 100
    internal_points, external_points, x_spacing, y_spacing, image = classify_points(image_path, num_points_x, num_points_y)
    internal_points, boundary_points_distances = identify_boundary_points_and_distances(internal_points, external_points, x_spacing, y_spacing, image)
    create_csv_with_point_types_and_distances(internal_points, external_points, boundary_points_distances, num_points_x, num_points_y, 'output.csv')
    print("CSV file created successfully.")

