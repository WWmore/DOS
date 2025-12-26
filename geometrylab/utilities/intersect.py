import numpy as np

# -----------------------------------------------------------------------------

from geometrylab.utilities.arrayutilities import orthogonal_vector

from geometrylab.geometry.frame import Frame
# -----------------------------------------------------------------------------


def points_plane_projection(points, plane_point=None,
                            plane_normal=np.array([0, 0, 1]),
                            light_direction=np.array([1, 1, 1]),
                            offset=0):
    """
    see:
    https://www.scratchapixel.com/lessons/3d-basic-rendering/
    ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
    """
    if plane_point is None:
        plane_point = np.array([0, 0, 0]) - plane_normal * offset
    try:
        e1 = orthogonal_vector(plane_normal[0])
    except:
        e1 = orthogonal_vector(plane_normal)
    e2 = np.cross(plane_normal, e1)
    A = np.tile(plane_point, (len(points), 1))
    O = points
    E_1 = np.tile(e1, (len(points), 1))
    E_2 = np.tile(e2, (len(points), 1))
    D = np.array(light_direction)
    T = O - A
    P = np.cross(D, E_2, axisb=1)
    Q = np.cross(T, E_1, axisb=1)
    f = (np.einsum('ij,ij -> i', P, E_1) + 1e-20) ** -1
    # t = f * np.einsum('ij,ij -> i', Q, E_2)
    u = f * np.einsum('ij,ij -> i', P, T)
    v = f * np.einsum('ij,j -> i', Q, D)
    i = (np.einsum('i,ij -> ij', u, E_1) + np.einsum('i,ij -> ij', v, E_2)) + A
    return i

def mesh_plane_projection(mesh, frame=None, light_direction=[1,-1,1], offset=None):
    '''
    see:
    https://www.scratchapixel.com/lessons/3d-basic-rendering/
    ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
    '''
    if frame is None:
        frame = Frame()
        if offset is None:
            frame.origin[:,2] = np.min(mesh.vertices[:,2])
        else:
            frame.origin =- frame.e3*offset
    A = np.tile(frame.origin, (mesh.V,1))
    O = mesh.vertices
    E_1 = np.tile(frame.e1, (mesh.V,1))
    E_2 = np.tile(frame.e2, (mesh.V,1))
    D = np.array(light_direction)
    T = O - A
    P = np.cross(D, E_2, axisb=1)
    Q = np.cross(T, E_1, axisb=1)
    f = (np.einsum('ij,ij -> i', P, E_1) + 1e-20)**-1
    #t = f * np.einsum('ij,ij -> i', Q, E_2)
    u = f * np.einsum('ij,ij -> i', P, T)
    v = f * np.einsum('ij,j -> i', Q, D)
    i = (np.einsum('i,ij -> ij', u, E_1) + np.einsum('i,ij -> ij', v, E_2)) + A
    return i


def mesh_plane_shadow(mesh, frame=None, light_direction=[1,-1,1], offset=None): ##Hui add
    P = mesh_plane_projection(mesh, frame=frame,
                              light_direction=light_direction, offset=offset)
    M = mesh.copy_mesh()
    M.vertices = P
    return M