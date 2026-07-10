"""The repair command is part of every installation."""

import trimesh

from dicom_surface import repair, validate


def test_repair_closes_an_open_mesh(tmp_path):
    box = trimesh.creation.box()
    box.update_faces([False, False] + [True] * (len(box.faces) - 2))

    source = str(tmp_path / "open.stl")
    output = str(tmp_path / "repaired.stl")
    box.export(source)

    before = validate.validate(source, self_intersections=False)
    assert not before["watertight"]

    stats = repair.repair(source, output)
    after = validate.validate(output, self_intersections=False)

    assert stats["holes_in"] > 0
    assert stats["holes_out"] == 0
    assert stats["holes_filled"] > 0
    assert after["watertight"]
