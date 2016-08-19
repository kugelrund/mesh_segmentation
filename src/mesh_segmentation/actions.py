import random
import bpy

def assignMaterials(mesh, k, idx):
    """Assigns a random colored material for each found segment"""

    # clear all existing materials
    while mesh.materials:
        mesh.materials.pop(0, update_data=True)

    for i in range(k):
        material = bpy.data.materials.new(''.join(['mat', mesh.name, str(i)]))
        material.diffuse_color = (random.random(),
                                  random.random(),
                                  random.random())
        mesh.materials.append(material)

    for i, id in enumerate(idx):
        mesh.polygons[i].material_index = id
