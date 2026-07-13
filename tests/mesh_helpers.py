"""Small MeshLib/NumPy fixtures shared by geometry tests."""

from __future__ import annotations

import meshlib.mrmeshnumpy as mrmeshnumpy
import meshlib.mrmeshpy as mrmeshpy
import numpy as np

MeshArrays = tuple[np.ndarray, np.ndarray]


def arrays(mesh: mrmeshpy.Mesh) -> MeshArrays:
    return (
        np.asarray(mrmeshnumpy.getNumpyVerts(mesh), dtype=np.float64).copy(),
        np.asarray(mrmeshnumpy.getNumpyFaces(mesh.topology), dtype=np.int32).copy(),
    )


def box(extents=(1.0, 1.0, 1.0)) -> MeshArrays:
    size = mrmeshpy.Vector3f(*map(float, extents))
    base = mrmeshpy.Vector3f(*(-0.5 * np.asarray(extents, dtype=float)))
    return arrays(mrmeshpy.makeCube(size, base))


def sphere(subdivisions: int = 2, radius: float = 1.0) -> MeshArrays:
    resolution = 8 * (2 ** subdivisions)
    return arrays(mrmeshpy.makeUVSphere(float(radius), resolution, resolution // 2))


def torus(
    major_radius: float,
    minor_radius: float,
    major_sections: int,
    minor_sections: int,
) -> MeshArrays:
    return arrays(
        mrmeshpy.makeTorus(
            float(major_radius),
            float(minor_radius),
            int(major_sections),
            int(minor_sections),
        )
    )


def transformed(
    mesh: MeshArrays,
    *,
    translation=(0.0, 0.0, 0.0),
    rotation_z_deg: float = 0.0,
) -> MeshArrays:
    vertices, faces = mesh
    angle = np.deg2rad(rotation_z_deg)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=float,
    )
    moved = vertices @ rotation.T + np.asarray(translation, dtype=float)
    return moved, faces.copy()


def concatenate(*meshes: MeshArrays) -> MeshArrays:
    vertices = []
    faces = []
    offset = 0
    for mesh_vertices, mesh_faces in meshes:
        vertices.append(mesh_vertices)
        faces.append(mesh_faces + offset)
        offset += len(mesh_vertices)
    return np.vstack(vertices), np.vstack(faces).astype(np.int32)


def write(path, mesh: MeshArrays) -> str:
    vertices, faces = mesh
    value = mrmeshnumpy.meshFromFacesVerts(
        np.ascontiguousarray(faces, dtype=np.int32),
        np.ascontiguousarray(vertices, dtype=np.float64),
    )
    mrmeshpy.saveMesh(value, str(path))
    return str(path)


def volume(mesh: MeshArrays) -> float:
    vertices, faces = mesh
    value = mrmeshnumpy.meshFromFacesVerts(faces.astype(np.int32), vertices)
    return abs(float(value.volume()))


def mean_adjacency_angle(mesh: MeshArrays) -> float:
    vertices, faces = mesh
    triangles = vertices[faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    edges: dict[tuple[int, int], list[int]] = {}
    for face_id, face in enumerate(faces):
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edges.setdefault(tuple(sorted((int(a), int(b)))), []).append(face_id)
    pairs = [ids for ids in edges.values() if len(ids) == 2]
    dots = np.asarray([np.clip(np.dot(normals[a], normals[b]), -1.0, 1.0) for a, b in pairs])
    return float(np.degrees(np.arccos(dots)).mean())
