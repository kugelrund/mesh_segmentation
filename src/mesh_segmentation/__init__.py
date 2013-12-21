bl_info = {
    "name": "Mesh Segmentation",
    "description": "Segments an object and applies an action on each segment",
    "category": "Object"}
    
if "bpy" in locals():
    import imp
    imp.reload(segmentation)
    imp.reload(actions)
else:
    from . import segmentation
    from . import actions

import bpy

#### developing purpose ###
import imp
imp.reload(segmentation)
imp.reload(actions)
###########################


class MeshSegmentation(bpy.types.Operator):
    """Segment a mesh"""
    
    bl_idname = "object.mesh_segmentation"
    bl_label = "Segment Mesh"
    bl_options = {'REGISTER', 'UNDO'}

    # properties set by the user
    exe = bpy.props.BoolProperty(name="Execute",
        description="If it shall actually run.",
        default=False)
    k = bpy.props.IntProperty(name="Clusters",
        description="Amount of clusters",
        default=2)
    delta = bpy.props.FloatProperty(name="Delta",
        description="Set close to zero for more importance on the angular distance, set close to one for more importance on the geodesic distance.",
        default=0.02,
        min=0,
        max=1,
        subtype='FACTOR')
    eta = bpy.props.FloatProperty(name="Weight of concativity",
        description="",
        default=0.15,
        min=0,
        max=1,
        subtype='FACTOR')

    def execute(self, context):
        """Executes the segmentation"""
        if self.exe:
            segmentation.segment_mesh(context.active_object.data, self.k, (self.delta, self.eta), actions.applyMaterials)
            self.exe = False
            
        return {'FINISHED'}
    
    
def register():
    """Registers the MeshSegmentation operator in blender"""
    bpy.utils.register_class(MeshSegmentation)

def unregister():
    """Unregisters the MeshSegmentation operator from blender"""
    bpy.utils.unregister_class(MeshSegmentation)

# developing purpose for registering when run from text editor    
if __name__ == "__main__":
    register()
