bl_info = {
    "name": "2D Mesh Visualizer",
    "description": "Visualize 2D mesh data",
    "author": "TheBusyTypist",
    "location": "View3D > Add > Mesh",
    "category": "Add Mesh"
}

import bpy
import bmesh

class Visualizer2D(bpy.types.Operator):
    """Visualize 2D mesh data"""
    bl_idname = "object.contour_2d"
    bl_label = "Visualize 2D mesh data"
    bl_options = {"REGISTER", "UNDO"}

    Source = bpy.props.StringProperty(
        name="File path",
        subtype="FILE_PATH")

    def execute(self, context):
        if self.Source != "":
            with open(self.Source, "r") as file:
                lines = file.readlines()
                vcnt, ecnt = map(int, lines[0].split())

                vertices = []
                for i in range(vcnt):
                    x, y = map(float, lines[i + 1].split())
                    vertices.append((x, y, 0.0))

                edges = []
                for i in range(ecnt):
                    a, b = map(int, lines[1 + ecnt + i].split())
                    edges.append((a, b))

                mesh = bpy.data.meshes.new("Mesh")
                mesh.from_pydata(vertices, edges, [])
                obj = bpy.data.objects.new("MeshObj", mesh)
                context.scene.objects.link(obj)
        
        return {"FINISHED"}

def register():
    bpy.utils.register_class(Visualizer2D)

def unregister():
    bpy.utils.unregister_class(Visualizer2D)

if __name__ == "__main__":
    register()