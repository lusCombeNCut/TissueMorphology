# Lumen Pressure Force — Derivation

This document derives the force on each cell from the hydrostatic lumen pressure,
starting from the work–energy relation and ending at the expressions implemented in
`LumenPressureForce.hpp`.

---

## 1. Starting point: pressure–volume work

A fluid-filled lumen at uniform hydrostatic pressure $P$ does work on its boundary
when the enclosed volume (or area, in 2-D) changes:

$$
\delta W = P \, \delta V
$$

In the overdamped (zero-inertia) regime used by Chaste the "energy" view and
the "virtual-work" view are interchangeable: the generalised force on degree
of freedom $\mathbf{r}_i$ is

$$
\mathbf{F}_i
= P \, \frac{\partial V}{\partial \mathbf{r}_i}
= P \, \nabla_i V
$$

where $V$ is the volume (3-D) or area (2-D) of the lumen cavity.

> **Convention.** $P > 0$ means the lumen pushes the boundary *outward*.
> Because $\nabla_i V$ points in the direction that increases $V$
> when $\mathbf{r}_i$ moves outward, the sign is automatically correct.

---

## 2. Two-dimensional case — the shoelace formula

### 2.1 Polygon area from cell positions

In 2-D the lumen cross-section is the polygon whose vertices are the cell
centres $\mathbf{r}_0, \mathbf{r}_1, \dots, \mathbf{r}_{N-1}$, ordered
counter-clockwise by angle around the centroid.

The signed area of this polygon is given by the **shoelace formula**:

$$
A = \frac{1}{2} \sum_{i=0}^{N-1}
  \bigl( x_i \, y_{i+1} - x_{i+1} \, y_i \bigr)
$$

where indices are taken modulo $N$, i.e. $\mathbf{r}_N \equiv \mathbf{r}_0$.

Writing out the terms that involve vertex $i$:

$$
A = \frac{1}{2} \Bigl[
  \;\cdots
  + \bigl( x_{i-1}\, y_i - x_i \, y_{i-1} \bigr)
  + \bigl( x_i \, y_{i+1} - x_{i+1} \, y_i \bigr)
  + \cdots
\Bigr]
$$

Only the $(i{-}1, i)$ and $(i, i{+}1)$ edge terms contain $x_i$ or $y_i$.

### 2.2 Differentiating with respect to vertex $i$

$$
\frac{\partial A}{\partial x_i}
= \frac{1}{2}
  \bigl( y_{i+1} - y_{i-1} \bigr)
$$

$$
\frac{\partial A}{\partial y_i}
= \frac{1}{2}
  \bigl( x_{i-1} - x_{i+1} \bigr)
$$

Collecting into a gradient vector:

$$
\nabla_i A
= \begin{pmatrix}
    \tfrac{1}{2}(y_{i+1} - y_{i-1}) \\[4pt]
    \tfrac{1}{2}(x_{i-1} - x_{i+1})
  \end{pmatrix}
$$

### 2.3 Force on cell $i$

$$
\boxed{
\mathbf{F}_i = P \, \nabla_i A
= \frac{P}{2}
  \begin{pmatrix}
    y_{i+1} - y_{i-1} \\
    x_{i-1} - x_{i+1}
  \end{pmatrix}
}
$$

### 2.4 Geometric interpretation

The gradient $\nabla_i A$ can be decomposed into contributions from the two
edges meeting at vertex $i$:

- Edge $(i{-}1, i)$ has outward normal $\mathbf{n}^- = \tfrac{1}{2}(y_i - y_{i-1},\; x_{i-1} - x_i)^\top$
- Edge $(i, i{+}1)$ has outward normal $\mathbf{n}^+ = \tfrac{1}{2}(y_{i+1} - y_i,\; x_i - x_{i+1})^\top$

Summing:

$$
\mathbf{n}^- + \mathbf{n}^+
= \frac{1}{2}
  \begin{pmatrix}
    y_{i+1} - y_{i-1} \\
    x_{i-1} - x_{i+1}
  \end{pmatrix}
= \nabla_i A
$$

So the shoelace gradient on vertex $i$ is exactly the sum of the
(un-normalised, length-weighted) outward normals of the two half-edges
adjacent to that vertex. This is the standard pressure–boundary integral
result: hydrostatic pressure acts normal to the surface, with magnitude
proportional to the edge length.

---

## 3. Three-dimensional case — spherical approximation

### 3.1 Volume gradient for a closed surface

For a closed surface $S$ enclosing volume $V$, the divergence theorem gives

$$
V = \frac{1}{3} \oint_S \mathbf{r} \cdot \hat{\mathbf{n}} \, dS
$$

If vertex $i$ at position $\mathbf{r}_i$ "owns" a surface patch of area
$\Delta S_i$ with outward unit normal $\hat{\mathbf{n}}_i$, then

$$
\nabla_i V \approx \Delta S_i \; \hat{\mathbf{n}}_i
$$

### 3.2 Area element on a sphere

For $N$ cells distributed roughly uniformly on a sphere of radius $r_i$
(the distance from the centroid to cell $i$):

$$
\Delta S_i = \frac{4 \pi r_i^2}{N}
$$

The outward unit normal at cell $i$ is simply the radial direction:

$$
\hat{\mathbf{n}}_i = \frac{\mathbf{r}_i - \mathbf{c}}{|\mathbf{r}_i - \mathbf{c}|}
$$

where $\mathbf{c}$ is the lumen centroid.

### 3.3 Force on cell $i$

$$
\boxed{
\mathbf{F}_i
= P \, \nabla_i V
\approx P \; \frac{4 \pi r_i^2}{N} \; \hat{\mathbf{n}}_i
}
$$

This is the expression used in the 3-D branch of `LumenPressureForce`.

> **Note.** This is an approximation valid for roughly spherical cell
> arrangements. For vertex-based 3-D models where a proper surface
> triangulation is available, `LumenPressureSubForce` (in OrganoidChaste)
> computes the exact surface gradient from the mesh faces instead.

---

## 4. Area / volume gradient per model type

The four model types differ in what the "degree of freedom" $\mathbf{r}_i$ is and how
the gradient $\nabla_i A$ (or $\nabla_i V$) is computed.

---

### 4.1 Node 2D (`LumenPressureForce<2>`, node-based population)

| | |
|---|---|
| **DOF** | Cell centre $\mathbf{r}_i = (x_i, y_i)$ — the single point representing cell $i$ |
| **Lumen boundary** | Polygon formed by sorting all cell centres by angle around the centroid |
| **Gradient** | Exact shoelace derivative (Section 2) |

$$
\nabla_i A =
\begin{pmatrix}
  \tfrac{1}{2}(y_{i+1} - y_{i-1}) \\[3pt]
  \tfrac{1}{2}(x_{i-1} - x_{i+1})
\end{pmatrix}
$$

Force applied **directly to the node**.

---

### 4.2 Node 3D (`LumenPressureForce<3>`, node-based population)

| | |
|---|---|
| **DOF** | Cell centre $\mathbf{r}_i$ on the surface of a ~spherical organoid |
| **Lumen boundary** | Inferred implicitly from the set of cell positions |
| **Gradient** | Spherical approximation (Section 3) |

$$
\nabla_i V \approx \frac{4\pi r_i^2}{N}\,\hat{\mathbf{n}}_i,
\qquad
\hat{\mathbf{n}}_i = \frac{\mathbf{r}_i - \mathbf{c}}{|\mathbf{r}_i - \mathbf{c}|}
$$

Force applied **directly to the node**.

This is an approximation; it becomes exact only in the limit of a uniform spherical
distribution of cells.

---

### 4.3 Vertex 2D (`LumenPressureForce<2>`, vertex-based population)

In a vertex model each "cell" is a polygonal element with several mesh nodes.
There is no single cell-centre degree of freedom — the actual DOFs are the
mesh **node** positions.

**Approach used here:** treat the element centroid as the representative position and
apply the shoelace gradient at that level, then distribute the resulting force equally
to all nodes of the element.

**Step 1 — centroid-level gradient.** Form the polygon of element centroids sorted by
angle (identical to the Node 2D case) and compute

$$
\mathbf{F}_i^{\text{cell}} = P \, \nabla_i A
= \frac{P}{2}
\begin{pmatrix}
  y_{i+1} - y_{i-1} \\
  x_{i-1} - x_{i+1}
\end{pmatrix}
$$

**Step 2 — node distribution.** If element $i$ has $n_i$ nodes, the force on each
node $k$ of that element is

$$
\mathbf{F}_{i,k}^{\text{node}} = \frac{\mathbf{F}_i^{\text{cell}}}{n_i}
$$

This equal-split is a simplification: the thermodynamically exact approach would
differentiate the area directly with respect to each individual mesh node position.
However, since all $n_i$ nodes are approximately equidistant from the centroid and
the lumen, the equal split is a good approximation and avoids coupling the lumen
geometry calculation to the vertex mesh topology.

---

### 4.4 Vertex 3D (`LumenPressureSubForce<3>`, OrganoidChaste monolayer mesh)

This model uses a proper 3-D vertex mesh where each cell element has explicitly
labelled **apical**, **basal**, and **lateral** faces. Only the apical face borders
the lumen.

**Exact gradient from face triangulation.**
The `CalculateLumenVolGradient` method in `MonolayerVertexMesh` computes
$\nabla_i V$ from the divergence theorem applied to the apical surface.

For a closed surface the divergence theorem gives

$$
V = \frac{1}{3} \oint_S \mathbf{r} \cdot \hat{\mathbf{n}} \, dS
$$

Differentiating with respect to the position of node $i$ (which lies on one or more
apical face triangles) yields contributions from the two triangles adjacent to node $i$
in the fan triangulation of its apical face:

$$
\nabla_i V
= \sum_{\substack{\text{apical triangles} \\ \text{containing node } i}}
  \frac{1}{2}\bigl(\mathbf{r}_a \times \mathbf{r}_b\bigr)
$$

where $\mathbf{r}_a$ and $\mathbf{r}_b$ are the positions of the other two vertices of
each triangle (the cross product gives the outward face-area vector of that triangle).
The $\tfrac{1}{2}$ comes from differentiating a triangle area $\tfrac{1}{2}|\mathbf{a} \times \mathbf{b}|$
with respect to a vertex.

The code distributes contribution symmetrically across the $n_f$ nodes of the face using
centre-of-face triangulation, so that all nodes receive an equal share of the face-area
vector regardless of face shape.

Force applied **directly to the mesh node** (no secondary distribution step needed,
because the DOFs are the actual mesh nodes).

---

### 4.5 Comparison summary

| Model | DOF | $\nabla_i A$ or $\nabla_i V$ | Exact? |
|-------|-----|-------------------------------|--------|
| Node 2D | cell centre | shoelace on sorted centres | **Yes** |
| Node 3D | cell centre | $\tfrac{4\pi r^2}{N}\hat{\mathbf{n}}$ | Approx (spherical) |
| Vertex 2D | centroid → equally to nodes | shoelace on centroids then $\div n_i$ | Approx (centroid proxy) |
| Vertex 3D (OrganoidChaste) | mesh node | cross-product sum over apical face triangles | **Yes** |

---

## 5. Summary of notation

| Symbol | Meaning |
|--------|---------|
| $P$ | Constant hydrostatic lumen pressure (parameter `mPressure`) |
| $V$ / $A$ | Lumen volume (3-D) or area (2-D) |
| $\mathbf{r}_i$ | Position of cell $i$ |
| $\nabla_i$ | Gradient with respect to $\mathbf{r}_i$ |
| $N$ | Number of cells on the lumen boundary |
| $\mathbf{c}$ | Lumen centroid (auto-tracked from cell positions) |
| $r_i$ | Distance from cell $i$ to centroid |
| $\hat{\mathbf{n}}_i$ | Outward unit normal at cell $i$ |

---

## 6. Code reference

The derivation above maps directly to `LumenPressureForce<DIM>::AddForceContribution`:

- **2-D:** cells sorted by angle → shoelace gradient → `force = mPressure * grad`
- **3-D:** radial direction + spherical area element → `force = mPressure * area_element * n_hat`

For vertex-based 2-D populations the per-cell force is distributed equally
to all nodes of the vertex element.
