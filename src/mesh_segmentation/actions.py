import random
import bpy

def assignMaterials(mesh, k, idx):
    """Assigns a random colored material for each found segment"""

    # clear all existing materials
    mesh.materials.clear()

    for i in range(k):
        material = bpy.data.materials.new(''.join(['mat', mesh.name, str(i)]))
        material.diffuse_color = (random.random(), random.random(),
                                  random.random(), 1.0)
        mesh.materials.append(material)

    for i, id in enumerate(idx):
        mesh.polygons[i].material_index = id
