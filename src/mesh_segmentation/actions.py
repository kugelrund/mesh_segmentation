import random
import bpy

def applyMaterials(mesh, k, idx):
    """Applys a random material for each found segment"""
    
    # clear all existing materials
    while mesh.materials:
        mesh.materials.pop(0, update_data=True)
      
    for i in range(k):
        material = bpy.data.materials.new(''.join(['mat', mesh.name, str(i)]))
        material.diffuse_color = (random.random(), random.random(), random.random())
        mesh.materials.append(material)
    
    for i in range(len(idx)):
        mesh.polygons[i].material_index = idx[i]        
    