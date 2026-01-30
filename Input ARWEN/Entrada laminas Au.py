import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Polygon, box, MultiPolygon
from shapely.ops import unary_union
import triangle as tr
import itertools
import os

# === FUNCIONES AUXILIARES ===
def get_rotated_rectangle(x, y, l, e, angle):
    angle_rad = np.radians(angle)
    A = np.array([x, y])
    D = A + [l * np.cos(angle_rad), l * np.sin(angle_rad)]
    C = D + [-e * np.sin(angle_rad), e * np.cos(angle_rad)]
    B = A + [-e * np.sin(angle_rad), e * np.cos(angle_rad)]
    return Polygon([A, D, C, B])

def polygon_to_triangle_data(polygon):
    coords = list(polygon.exterior.coords)[:-1]
    vertices = np.array(coords)
    segments = [(i, (i+1) % len(vertices)) for i in range(len(vertices))]

    holes_coords = []
    for interior in polygon.interiors:
        hole_coords = list(interior.coords)[:-1]
        holes_coords.append(np.mean(hole_coords, axis=0))
        hole_indices = range(len(vertices), len(vertices) + len(hole_coords))
        segments.extend([(hole_indices[i], hole_indices[(i+1) % len(hole_coords)]) for i in range(len(hole_coords))])
        vertices = np.vstack([vertices, hole_coords])

    return {
        'vertices': vertices,
        'segments': segments,
        'holes': holes_coords
    }

def is_convex(polygon):
    coords = list(polygon.exterior.coords)[:-1]
    if len(coords) < 3:
        return False
    sign = None
    for i in range(len(coords)):
        dx1 = coords[(i+1) % len(coords)][0] - coords[i][0]
        dy1 = coords[(i+1) % len(coords)][1] - coords[i][1]
        dx2 = coords[(i+2) % len(coords)][0] - coords[(i+1) % len(coords)][0]
        dy2 = coords[(i+2) % len(coords)][1] - coords[(i+1) % len(coords)][1]
        cross = dx1 * dy2 - dy1 * dx2
        new_sign = cross >= 0
        if sign is None:
            sign = new_sign
        elif sign != new_sign:
            return False
    return True

def group_triangles_into_convex_regions(triangles):
    remaining = triangles.copy()
    changed = True
    while changed:
        changed = False
        new_polys = []
        used = set()
        for i, j in itertools.combinations(range(len(remaining)), 2):
            if i in used or j in used:
                continue
            tri1 = remaining[i]
            tri2 = remaining[j]
            shared = tri1.intersection(tri2)
            if shared.length > 0:
                union = tri1.union(tri2)
                if isinstance(union, Polygon) and is_convex(union):
                    new_polys.append(union)
                    used.update([i, j])
                    changed = True
                    break
        for idx, tri in enumerate(remaining):
            if idx not in used:
                new_polys.append(tri)
        remaining = new_polys
    return remaining

# === PARÁMETROS Y GENERACIÓN DE GEOMETRÍA ===

T = 300  # K
ux = 0.0
uy = 0.0

# Material A (CH)
A_width = 600  # um
A_height = 400  # um (duplicado)
rho_A = 0.7  # kg/m^3
material_A_poly = box(0, 0, A_width, A_height)

# Material B (Cu)
num_rectangles = 31
l = 35 # um
e = 8  # um
rho_B = 8970  # kg/m^3

material_B_polys = []
placed = 0
max_attempts = 5000

# Material C (He)
C_width = A_width
C_height = 100  # um
rho_C = 0.1785  # kg/m^3

# Lámina inferior de He
material_C_poly = box(0, 0, C_width, -C_height)
# Lámina superior de He (nueva)
material_D_poly = box(0, A_height, C_width, A_height + C_height)

# Generación de Au en CH
while placed < num_rectangles and max_attempts > 0:
    xr = np.random.uniform(0, A_width)
    yr = np.random.uniform(0, A_height)
    angle = np.random.uniform(0, 180)
    x = xr - l / 2 * np.cos(np.radians(angle))
    y = yr - l / 2 * np.sin(np.radians(angle))
    poly = get_rotated_rectangle(x, y, l, e, angle)
    if material_A_poly.contains(poly):
        if all(not poly.intersects(b) for b in material_B_polys):
            material_B_polys.append(poly)
            placed += 1
    max_attempts -= 1

# Crear regiones A - B trianguladas
material_B_union = unary_union(material_B_polys)
A_minus_B = material_A_poly.difference(material_B_union)

triangulated_polys = []
if isinstance(A_minus_B, Polygon):
    data = polygon_to_triangle_data(A_minus_B)
    t = tr.triangulate(data, 'p')
    for tri_indices in t['triangles']:
        pts = t['vertices'][tri_indices]
        triangulated_polys.append(Polygon(pts))
elif isinstance(A_minus_B, MultiPolygon):
    for poly in A_minus_B.geoms:
        data = polygon_to_triangle_data(poly)
        t = tr.triangulate(data, 'p')
        for tri_indices in t['triangles']:
            pts = t['vertices'][tri_indices]
            triangulated_polys.append(Polygon(pts))

convex_parts = group_triangles_into_convex_regions(triangulated_polys)

# === CUENTAS DE REGIONES ===
num_A = len(convex_parts)
num_B = len(material_B_polys)
num_C = 1  # He inferior
num_D = 1  # He superior
total_regions = num_A + num_B + num_C + num_D

# === PROPIEDADES ===
densities = ([rho_A] * num_A) + ([rho_B] * num_B) + ([rho_C] * (num_C + num_D))
temperatures = [T] * total_regions
uxs = [ux] * total_regions
uys = [uy] * total_regions
materials = ([1] * num_A) + ([2] * num_B) + ([3] * num_C) + ([3] * num_D)  # He = material 3

# === DENSIDAD MEDIA DE CH + Au CON ROTACIÓN ===
import math

S_foam = A_width * A_height
S_B = num_B * l * e
S_A = S_foam - S_B

rho_promedio = (S_A * rho_A + S_B * rho_B) / S_foam if S_foam > 0 else 0

print(f"Densidad media CH+Au (revolución x=0): {rho_promedio:.4f} kg/m³")

# === VISUALIZACIÓN ===
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(0, A_width)
ax.set_ylim(-C_height, A_height + C_height)
ax.set_aspect('equal')

for tri in convex_parts:
    coords = np.array(tri.exterior.coords)
    ax.add_patch(patches.Polygon(coords, closed=True, facecolor='lightblue', edgecolor='black', alpha=0.8))

for poly_B in material_B_polys:
    coords = np.array(poly_B.exterior.coords)
    ax.add_patch(patches.Polygon(coords, closed=True, facecolor='gold', edgecolor='black', alpha=0.9))

coords_C = np.array(material_C_poly.exterior.coords)
coords_D = np.array(material_D_poly.exterior.coords)
ax.add_patch(patches.Polygon(coords_C, closed=True, facecolor='lightgreen', edgecolor='black', alpha=0.7))
ax.add_patch(patches.Polygon(coords_D, closed=True, facecolor='lightgreen', edgecolor='black', alpha=0.7))

plt.title("CH con Au + Láminas de He arriba y abajo")
plt.xlabel("X (μm)")
plt.ylabel("Y (μm)")
plt.tight_layout()
plt.show()

output_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "Subregiones_input_rod_ARWEN_optimo.txt"
)

# === GUARDAR EN ARCHIVOS ===
with open(output_path, "w") as f:
    f.write(f'# REGIONES EN LAS QUE SE DIVIDE EL MATERIAL CH: {len(convex_parts)}\n')
    f.write(f'# REGIONES EN LAS QUE SE DIVIDE EL MATERIAL Cu: {len(material_B_polys)}\n')
    f.write(f'# REGIONES EN LAS QUE SE DIVIDE EL MATERIAL He: 2\n')  #

    f.write(f'# Número de rectángulos de Cu N = {num_rectangles}\n')
    f.write(f'# Longitud de los rectángulos l = {l} um\n')
    f.write(f'# Espesor de los rectángulos e = {e} um\n')

    f.write(f'# Densidad CH: {rho_A} kg/m3\n')
    f.write(f'# Densidad Cu: {rho_B} kg/m3\n')
    f.write(f'# Densidad media de la espuma: {rho_promedio} kg/m3\n')


    # Ahora imprimir en formato solicitado
    f.write(f"reg.regions = {total_regions}\n")
    f.write("reg.den     = " + " ".join(f"{d:.1f}" for d in densities) + " # kg/m^3")
    f.write('\n')
    f.write("reg.tem     = " + " ".join(f"{t:.1f}" for t in temperatures))
    f.write('\n')
    f.write("reg.ux      = " + " ".join(f"{v:.1f}" for v in uxs))
    f.write('\n')
    f.write("reg.uy      = " + " ".join(f"{v:.1f}" for v in uys))
    f.write('\n')
    f.write("reg.mat     = " + " ".join(str(m) for m in materials))
    f.write('\n')

    f.write('\nreg.np = ')
    vertex_counts = []
    # Regiones de material A
    for poly in convex_parts:
        num_vertices = len(list(poly.exterior.coords))
        vertex_counts.append(str(num_vertices))

    # Regiones de material B
    for poly in material_B_polys:
        num_vertices = len(list(poly.exterior.coords))
        vertex_counts.append(str(num_vertices))

    # Región de material C (solo 1 polígono)
    num_vertices_C = len(list(material_C_poly.exterior.coords))
    vertex_counts.append(str(num_vertices_C))

    # Región de material D (solo 1 polígono)
    num_vertices_D = len(list(material_D_poly.exterior.coords))
    vertex_counts.append(str(num_vertices_D))

    # Escribir en una sola línea
    f.write(" ".join(vertex_counts) + "\n")

    f.write('\nreg.points = \n')
    
    # Coordenadas material A
    for poly in convex_parts:
        coords = list(poly.exterior.coords)
        for x, y in coords:
            x_m = x * 1e-6
            y_m = y * 1e-6
            f.write(f"{x_m:.15e} {y_m:.15e}\n")
        f.write("\n")

    # Coordenadas material B
    for poly in material_B_polys:
        coords = list(poly.exterior.coords)
        for x, y in coords:
            x_m = x * 1e-6
            y_m = y * 1e-6
            f.write(f"{x_m:.15e} {y_m:.15e}\n")
        f.write("\n")

    # Coordenadas material C
    coords_C = list(material_C_poly.exterior.coords)
    for x, y in coords_C:
        x_m = x * 1e-6
        y_m = y * 1e-6
        f.write(f"{x_m:.15e} {y_m:.15e}\n")
    f.write("\n")

    # Coordenadas material D
    coords_C = list(material_D_poly.exterior.coords)
    for x, y in coords_D:
        x_m = x * 1e-6
        y_m = y * 1e-6
        f.write(f"{x_m:.15e} {y_m:.15e}\n")
    f.write("\n")

