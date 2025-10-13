# FEM Thermal Analysis for 2D Square Conduction

This project applies the **Finite Element Method (FEM)** to solve **steady-state heat conduction** problems in a 2D square domain using **Python**.  
It models the temperature distribution across a square plate by incorporating boundary conditions and various levels of mesh refinement for improved accuracy.

---

## 🧩 Components

### **1. 2D Heat Transfer in a Square Domain**
- Simulates heat conduction in a square region using **triangular elements**.  
- Implements **Dirichlet boundary conditions** (fixed temperatures) on the square edges.  
- Solves for **temperature distribution** under both uniform and non-uniform heat sources.

---

### **2. Mesh and Elements**
- The domain is discretized using **triangular elements**, with each node representing a temperature point.  
- **Mesh refinement** is applied using a **dichotomy algorithm**, subdividing triangles to increase accuracy.  
- The refined mesh results in a more precise temperature profile but increases computation time.

---

### **3. Thermal Resistance Analysis**
- **Thermal resistance** is computed from the temperature distribution.  
- Results show how increasing mesh density improves solution accuracy while raising computational cost.

---

### **4. Python Implementation**
- Built with **NumPy** and **SciPy** for numerical computation.  
- **Matplotlib** is used for temperature distribution visualization.  
- The algorithm assembles the **global stiffness matrix** and **force vector**, applying boundary conditions and heat source terms.

### **5. Customizable Parameters**
The simulation offers several parameters that can be modified to explore different thermal behaviors:

- **Heat Source:** Can be included or excluded to observe the effect of internal heat generation.  
- **Mesh Refinement:** The level of mesh refinement can be adjusted for higher accuracy.  
- **Convection:** Optional convection boundary conditions simulate heat exchange with the surrounding environment.  
- **Thermal Conductivity & Source Term:** Both thermal conductivity and internal heat source can be varied to analyze their impact on the temperature distribution.
---

## 📁 Python Files

| File | Description |
|------|--------------|
| `FEM_ST~1.PY` | Implements the FEM formulation for the 2D heat conduction problem. It defines nodes, elements, and material properties. |
| `MESH_2~1.PY` | Handles 2D mesh generation and refinement (triangle subdivision). |
| `TEST_C~1.PY` | Executes the main simulation, applies mesh refinement, and visualizes the temperature distribution. |

## 📊 Results Visualization

### Figure 1. Temperature Distribution without heat source
![Temperature Distribution](https://github.com/HenokMD/FEM-Thermal-Heat-Transfer/blob/main/Temperature%20distribution%20in%20the%20square%20area%2C%20without%20a%20heat%20source.png)

### Figure 2. Temperature Distribution with a heat source
![Temperature Distribution](https://github.com/HenokMD/FEM-Thermal-Heat-Transfer/blob/main/Temperature%20distribution%20in%20the%20square%20area%20with%20a%20heat%20source.png)

### Figure 3. Thermal Analysis of Square Area
![Thermal Analysis](https://github.com/HenokMD/FEM-Thermal-Heat-Transfer/blob/main/Thermal%20Analysis%20of%20Square%20Area.png)


---

## ⚙️ Requirements

Install dependencies using:
```bash
pip install numpy scipy matplotlib
