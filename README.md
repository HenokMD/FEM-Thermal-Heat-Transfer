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

---

## 📁 Python Files

| File | Description |
|------|--------------|
| `FEM_ST~1.PY` | Implements the FEM formulation for the 2D heat conduction problem — defines nodes, elements, and material properties. |
| `MESH_2~1.PY` | Handles 2D mesh generation and refinement (triangle subdivision). |
| `TEST_C~1.PY` | Executes the main simulation, applies mesh refinement, and visualizes the temperature distribution. |

---

## ⚙️ Requirements

Install dependencies using:
```bash
pip install numpy scipy matplotlib
