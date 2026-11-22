import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np


class ObstacleLayerGraphTests(unittest.TestCase):
    def test_obstacle_patch_layers_expand_by_face_adjacency(self):
        from pressure_SM_mesh_graph._3D.mesh_graph import protected_cells_from_layers

        owner = np.arange(7, dtype=np.int64)
        neighbour = np.arange(1, 8, dtype=np.int64)
        protected = protected_cells_from_layers(8, owner, neighbour, np.array([0]), 3)
        np.testing.assert_array_equal(protected, [True, True, True, False, False, False, False, False])

    def test_protected_cells_remain_singletons_at_each_level(self):
        from pressure_SM_mesh_graph._3D.mesh_graph import coarsen_graph, protected_cells_from_layers

        centres = np.column_stack((np.arange(8.), np.zeros((8, 2))))
        owner = np.arange(7, dtype=np.int64)
        neighbour = np.arange(1, 8, dtype=np.int64)
        protected = protected_cells_from_layers(8, owner, neighbour, np.array([0]), 3)
        graph = coarsen_graph(centres, np.ones(8), owner, neighbour, levels=2,
                              protected_cells=protected, protect_after=0)

        aggregate_sizes = np.bincount(graph['mapping'], minlength=graph['coarse_count'])
        np.testing.assert_array_equal(aggregate_sizes[graph['mapping'][protected]], 1)
        self.assertLess(graph['coarse_count'], 8)

    def test_obstacle_protection_configuration_is_validated(self):
        from pressure_SM_mesh_graph._3D.graph_io import validate_config

        config = {'add_gradDpPrev_input': True,
                  'graph_obstacle_patch': 'obstacle',
                  'graph_obstacle_protection_layers': 15}
        validate_config(config)

        for invalid in (-1, 1.5, True):
            with self.subTest(invalid=invalid):
                with self.assertRaisesRegex(ValueError, 'nonnegative integer'):
                    validate_config({**config, 'graph_obstacle_protection_layers': invalid})

        with self.assertRaisesRegex(ValueError, 'cannot be combined'):
            validate_config({**config, 'graph_adaptive_coarsening': True,
                             'graph_adaptive_extra_levels': 2,
                             'graph_adaptive_regions': [(0, 1, 0, 1, 0, 1)]})

    def test_graph_cache_tracks_obstacle_protection(self):
        from pressure_SM_mesh_graph._3D.graph_io import load_mesh_graph

        centres = np.column_stack((np.arange(8.), np.zeros((8, 2))))
        mesh_data = {
            'owner': np.r_[np.arange(7), 0],
            'neighbour': np.arange(1, 8),
            'boundary': [('obstacle', {'startFace': 7, 'nFaces': 1, 'type': 'wall'})],
        }

        class FakeFoamFile:
            def __init__(self, path):
                self.name = Path(path).name

            def __getitem__(self, key):
                return mesh_data[self.name]

        with tempfile.TemporaryDirectory() as directory, patch(
                'pressure_SM_mesh_graph._3D.graph_io.FoamFile', FakeFoamFile):
            unprotected, unprotected_signature = load_mesh_graph(centres, np.ones(8), directory, directory, 2)
            protected, protected_signature = load_mesh_graph(centres, np.ones(8), directory, directory, 2,
                                                               obstacle_layers=3)
            self.assertNotEqual(protected_signature, unprotected_signature)
            self.assertEqual(int(protected['protected_cell_count']), 3)
            self.assertLess(unprotected['coarse_count'], protected['coarse_count'])
            cached, cached_signature = load_mesh_graph(centres, np.ones(8), directory, directory, 2,
                                                        obstacle_layers=3)
            self.assertEqual(cached_signature, protected_signature)
            np.testing.assert_array_equal(cached['mapping'], protected['mapping'])
            with self.assertRaisesRegex(ValueError, 'was not found'):
                load_mesh_graph(centres, np.ones(8), directory, directory, 2,
                                obstacle_patch='missing', obstacle_layers=3)


if __name__ == '__main__':
    unittest.main()
