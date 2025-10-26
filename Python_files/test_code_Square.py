from fem_student_Square import Node, Element, Region
from mesh_2D_Square import Mesh2D

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import time

def create_refined_mesh_with_boundary_conditions():
    # Creation an object "Region" containing the physical properties (thermal conductivity and source)
    wall_region = Region(" wall ", 1.0, 100)

    # the list of initial mesh nodes with explicit Dirichlet boundary conditions
    n0 = Node(0., 0., 10.0)
    n0.set_dirichlet(10.0)
    n1 = Node(1., 0., 20.0)
    n1.set_dirichlet(20.0)
    n2 = Node(1., 1., 20.0)
    n2.set_dirichlet(20.0)
    n3 = Node(0., 1., 10.0)
    n3.set_dirichlet(10.0)
    nodes = list()
    nodes.append(n0)
    nodes.append(n1)
    nodes.append(n2)
    nodes.append(n3)

    # the list of initial mesh elements
    elements = list()
    elements.append(Element([nodes[0], nodes[1], nodes[2]], wall_region))
    elements.append(Element([nodes[0], nodes[2], nodes[3]], wall_region))

    # Creation of a mesh object from the nodes and the elements lists
    mesh = Mesh2D(nodes, elements)

    # mesh refine using Dichotomy algorithm
    dicot = 2  # Mesh refinement level
    mesh.refine(dicot)

    return mesh

def solve_heat_distribution(mesh):
    # Separate known and unknown nodes
    free_nodes = [node for node in mesh.nodes if not node.is_dirichlet]
    dirichlet_nodes = [node for node in mesh.nodes if node.is_dirichlet]

    # Construct global stiffness matrix
    M = np.zeros((len(mesh.nodes), len(mesh.nodes)))
    # New: Force vector to capture source term contributions
    F = np.zeros(len(mesh.nodes))

    for element in mesh.elements:
        # Stiffness matrix contribution
        small_matrix = [
            [element.m00(), element.m01(), element.m02()],
            [element.m01(), element.m11(), element.m12()],
            [element.m02(), element.m12(), element.m22()]
        ]

        # Source term contribution
        # Distribute the source term equally among element nodes
        source_value = element.region.source
        for i in range(3):
            node_id = element.nodes[i].ID
            # Assuming uniform source distribution in the element
            # The factor 1/3 distributes the source term equally to each node of the triangle
            F[node_id] += source_value * element.S * (1/3)

        # Assemble stiffness matrix
        for i in range(3):
            id_i = element.nodes[i].ID
            for j in range(3):
                id_j = element.nodes[j].ID
                M[id_i, id_j] += small_matrix[i][j]

    # Separate known and unknown nodes
    temp = len(free_nodes)
    NewM = np.zeros((temp, temp))
    NewA = np.zeros((temp, len(dirichlet_nodes)))
    NewF = np.zeros(temp)

    # Stiffness matrix for free nodes
    for i in range(temp):
        for j in range(temp):
            NewM[i][j] = M[free_nodes[i].ID, free_nodes[j].ID]

    # Boundary condition matrix
    for i in range(temp):
        for j in range(len(dirichlet_nodes)):
            NewA[i][j] = M[free_nodes[i].ID, dirichlet_nodes[j].ID]

    # Force vector for free nodes
    # Subtract the contribution from Dirichlet nodes
    for i in range(temp):
        NewF[i] = F[free_nodes[i].ID]
        for j in range(len(dirichlet_nodes)):
            NewF[i] -= NewA[i][j] * dirichlet_nodes[j].value

    # Prepare boundary temperatures
    MatrixT = np.zeros((len(dirichlet_nodes), 1))
    for i, node in enumerate(dirichlet_nodes):
        MatrixT[i][0] = node.value

    # Compute solution
    T = np.linalg.solve(NewM, NewF)

    # Assign temperatures to nodes
    node_values = np.zeros(len(mesh.nodes))
    for node in dirichlet_nodes:
        node_values[node.ID] = node.value

    for i, node in enumerate(free_nodes):
        # Extract the scalar value from the 1D array
        node_values[node.ID] = T[i][0] if T[i].ndim > 0 else T[i]

    return mesh, node_values


def plot_temperature_distribution(mesh, node_values):
    # Prepare data for plotting
    x = []
    y = []
    values = []

    for node in mesh.nodes:
        x.append(node.x)
        y.append(node.y)
        values.append(node_values[node.ID])

    # Create 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot surface with color representing temperature
    scatter = ax.plot_trisurf(x, y, values, linewidth=0.2, antialiased=True, cmap=cm.coolwarm)

    # Add colorbar
    plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=10, label='Temperature')

    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Temperature')
    ax.set_title('Temperature Distribution')

    plt.tight_layout()
    plt.show()

    # Print temperature details for verification
    print("\nTemperature Details:")
    print("Min value:", min(values))
    print("Max value:", max(values))
    print("Unique values:", set(values))


def print_mesh_details(mesh, node_values):
    """
    Prints the details of the mesh, including element IDs, node coordinates, and node temperature values.
    """
    print("\nMesh Details:")
    k = 0
    for element in mesh.elements:
        print(f"Element {k}:")
        for node in element.nodes:
            print(f"  Node ID: {node.ID}, Coordinates: ({node.x}, {node.y}), Temperature: {node_values[node.ID]}")
        k += 1

def calculate_thermal_resistance(mesh, node_values):
    """
    Calculate thermal resistance based on the given formula
    R_th = (T_1 - T_0)^2 / (∬_Ω λ(grad(T)).(grad(T)) dΩ)

    :param mesh: Mesh2D object
    :param node_values: Array of node temperature values
    :return: Thermal resistance value
    """
    # Find max and min temperature nodes
    max_temp_node = max(mesh.nodes, key=lambda node: node_values[node.ID])
    min_temp_node = min(mesh.nodes, key=lambda node: node_values[node.ID])

    # Temperature difference
    delta_T = max_temp_node.value - min_temp_node.value

    # Compute the volume integral term
    volume_integral = 0.0

    # Iterate through all elements
    for element in mesh.elements:
        # Compute gradient for each shape function
        grad_phi1x = element.grad_phi1x()
        grad_phi1y = element.grad_phi1y()
        grad_phi2x = element.grad_phi2x()
        grad_phi2y = element.grad_phi2y()
        grad_phi3x = element.grad_phi3x()
        grad_phi3y = element.grad_phi3y()

        # Compute temperature gradients at nodes
        T1 = node_values[element.nodes[0].ID]
        T2 = node_values[element.nodes[1].ID]
        T3 = node_values[element.nodes[2].ID]

        # Gradient of temperature in the element
        grad_T_x = T1 * grad_phi1x + T2 * grad_phi2x + T3 * grad_phi3x
        grad_T_y = T1 * grad_phi1y + T2 * grad_phi2y + T3 * grad_phi3y

        # Conductivity of the element
        lamda = element.lamda

        # Add to volume integral
        # The element area is used as the differential volume
        volume_integral += lamda * (grad_T_x**2 + grad_T_y**2) * element.S

    # Compute thermal resistance
    # Add a small epsilon to avoid division by zero
    thermal_resistance = (delta_T**2) / (volume_integral + 1e-10)

    return {
        'thermal_resistance': thermal_resistance,
        'max_temp': max_temp_node.value,
        'min_temp': min_temp_node.value,
        'delta_T': delta_T
    }

def analyze_mesh_refinement(source_term=100, max_refinement=5):
    """
    Analyze thermal resistance and computation time
    for different mesh refinement levels

    :param source_term: Source term value
    :param max_refinement: Maximum number of refinement iterations
    :return: Tuple of lists (refinement levels, thermal resistances, computation times)
    """
    refinement_levels = []
    thermal_resistances = []
    computation_times = []

    for dicot in range(max_refinement + 1):
        # Create mesh with boundary conditions
        wall_region = Region(" wall ", 1.0, source_term)

        # Initial mesh creation
        n0 = Node(0., 0., 10.0)
        n0.set_dirichlet(10.0)
        n1 = Node(1., 0., 20.0)
        n1.set_dirichlet(20.0)
        n2 = Node(1., 1., 20.0)
        n2.set_dirichlet(20.0)
        n3 = Node(0., 1., 10.0)
        n3.set_dirichlet(10.0)
        nodes = [n0, n1, n2, n3]

        elements = [
            Element([nodes[0], nodes[1], nodes[2]], wall_region),
            Element([nodes[0], nodes[2], nodes[3]], wall_region)
        ]

        # Create mesh
        mesh = Mesh2D(nodes, elements)

        # Refine mesh
        mesh.refine(dicot)

        # Solve heat distribution and measure time
        start_time = time.time()
        mesh, node_values = solve_heat_distribution(mesh)
        solve_time = time.time() - start_time

        # Calculate thermal resistance
        resistance_data = calculate_thermal_resistance(mesh, node_values)

        # Store results
        refinement_levels.append(dicot)
        thermal_resistances.append(resistance_data['thermal_resistance'])
        computation_times.append(solve_time)

        # Optional: print details for each refinement level
        print(f"\nRefinement Level {dicot}:")
        print(f"Nodes: {len(mesh.nodes)}")
        print(f"Elements: {len(mesh.elements)}")
        print(f"Thermal Resistance: {resistance_data['thermal_resistance']}")
        print(f"Computation Time: {solve_time:.6f} seconds")
        print(f"Temperature Range: {resistance_data['min_temp']:.2f} - {resistance_data['max_temp']:.2f}")

    # Plot results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Thermal Resistance vs Refinement Level
    ax1.plot(refinement_levels, thermal_resistances, marker='o')
    ax1.set_xlabel('Refinement Level')
    ax1.set_ylabel('Thermal Resistance')
    ax1.set_title('Thermal Resistance vs Mesh Refinement')

    # Computation Time vs Refinement Level
    ax2.plot(refinement_levels, computation_times, marker='o', color='red')
    ax2.set_xlabel('Refinement Level')
    ax2.set_ylabel('Computation Time (seconds)')
    ax2.set_title('Computation Time vs Mesh Refinement')

    plt.tight_layout()
    plt.show()

    return refinement_levels, thermal_resistances, computation_times



def main():
    # Create mesh with boundary conditions
    mesh = create_refined_mesh_with_boundary_conditions()

    # Solve heat distribution
    mesh, node_values = solve_heat_distribution(mesh)

    # Print mesh details (integrating the last part of the code)
    print_mesh_details(mesh, node_values)

    # Plot temperature distribution
    plot_temperature_distribution(mesh, node_values)

    analyze_mesh_refinement()

if __name__ == '__main__':
    main()

#
