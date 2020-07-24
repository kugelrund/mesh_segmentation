bl_info = {
    "name": "Mesh Segmentation",
    "description": "Segments an object and applies an action on each segment",
    "blender": (2, 80, 0),
    "category": "Mesh"}

import bpy
from mesh_segmentation import segmentation
from mesh_segmentation import actions
# developing purpose for reloading modules if already imported
import imp
imp.reload(segmentation)
imp.reload(actions)


class MeshSegmentation(bpy.types.Operator):
    """Segment a mesh"""

    bl_idname = "mesh.mesh_segmentation"
    bl_label = "Segment Mesh"
    bl_options = {'REGISTER', 'UNDO'}

    # parameters
    action: bpy.props.EnumProperty(name ="Action",
                                   items = [('assignMaterials',
                                             "Assign materials",
                                             "Assigns a different material for "
                                             "each found segment")],
                                   description = "What to do with the "
                                                 "segmentation",
                                   default = 'assignMaterials')
    k: bpy.props.IntProperty(name = "Clusters",
                             description = "Amount of clusters",
                             min = 2,
                             default = 2)
    delta: bpy.props.FloatProperty(name = "Delta",
                                   description = "Set close to zero for more "
                                                 "importance on the angular "
                                                 "distance, set close to one "
                                                 "for more importance on the "
                                                 "geodesic distance.",
                                   default = 0.03,
                                   min = 0,
                                   max = 1,
                                   subtype = 'FACTOR')
    eta: bpy.props.FloatProperty(name = "Weight of convexity",
                                 description = "Set close to zero for more "
                                               "importance on concave angles, "
                                               "set close to one to treat "
                                               "concave and convex angles "
                                               "equally.",
                                 default = 0.15,
                                 min = 1e-10,
                                 max = 1,
                                 subtype = 'FACTOR')
    ev_method: bpy.props.EnumProperty(name = "EV method",
        items = [('sparse', "Sparse", "Sparse method for eigenvector "
                                      "computation (scipy.sparse.linalg.eigsh)"),
                 ('dense', "Dense", "Dense method for eigenvector computation "
                                    "(scipy.linalg.eigh)")],
        description = "Method to use for eigenvector computation. 'Sparse' "
                      "should usually preferred, as it tends to be much faster "
                      "without sacrificing quality. 'Dense' can be tried as "
                      "fallback. Default",
        default = 'sparse')
    kmeans_init: bpy.props.EnumProperty(name = "k-means initialization",
        items = [('liu_zhang', "Liu & Zhang", "Initialization by Liu & Zhang"),
                 ('kmeans++', "k-means++", "Initialization from k-means++")],
        description = "Method to use for initializing centroids for k-means.",
        default = 'liu_zhang')

    def execute(self, context):
        """Executes the segmentation"""
        if bpy.ops.mesh.separate(type='LOOSE') != {'CANCELLED'}:
            self.report({'ERROR'}, "Separated not connected parts, choose "
                                   "one of them for segmentation!")
            return {'CANCELLED'}
        else:
            segmentation.segment_mesh(mesh = context.active_object.data,
                                      k = self.k,
                                      coefficients = (self.delta, self.eta),
                                      action = getattr(actions, self.action),
                                      ev_method = self.ev_method,
                                      kmeans_init = self.kmeans_init)
            return {'FINISHED'}

    def invoke(self, context, event):
        if context.active_object.type == 'MESH':
            return context.window_manager.invoke_props_dialog(self)
        else:
            self.report({'ERROR'}, "Selected object is not a mesh!")
            return {'CANCELLED'}


def register():
    """Registers the addon in blender"""
    bpy.utils.register_class(MeshSegmentation)


def unregister():
    """Unregisters the addon from blender"""
    bpy.utils.unregister_class(MeshSegmentation)


# developing purpose for registering when run from blender texteditor
if __name__ == "__main__":
    register()
