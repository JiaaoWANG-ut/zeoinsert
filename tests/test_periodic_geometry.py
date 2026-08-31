import itertools
import unittest

import numpy as np

from periodic_geometry import minimum_image


def cell_from_parameters(a, b, c, alpha, beta, gamma):
    alpha, beta, gamma = np.deg2rad([alpha, beta, gamma])
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * np.cos(gamma), b * np.sin(gamma), 0.0])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return np.vstack([va, vb, [cx, cy, cz]])


def brute_force(vec, cell, radius=6):
    frac = vec @ np.linalg.inv(cell)
    shifts = np.asarray(
        list(itertools.product(range(-radius, radius + 1), repeat=3)),
        dtype=float,
    )
    candidates = (frac[None, :] - shifts) @ cell
    return candidates[np.argmin(np.einsum("ni,ni->n", candidates, candidates))]


class MinimumImageTests(unittest.TestCase):
    def test_orthogonal_cell(self):
        cell = np.diag([10.0, 12.0, 14.0])
        vec = np.array([6.0, -7.0, 8.0])
        np.testing.assert_allclose(
            minimum_image(vec, cell, np.linalg.inv(cell)),
            [-4.0, 5.0, -6.0],
        )

    def test_vectorized_shape_is_preserved(self):
        cell = np.diag([10.0, 10.0, 10.0])
        vecs = np.array([[[6.0, 0.0, 0.0], [0.0, -6.0, 0.0]]])
        result = minimum_image(vecs, cell)
        self.assertEqual(result.shape, vecs.shape)
        np.testing.assert_allclose(result, [[[-4.0, 0.0, 0.0], [0.0, 4.0, 0.0]]])

    def test_project_cells_match_brute_force(self):
        cells = [
            cell_from_parameters(14.50863, 14.50863, 14.50863, 60, 60, 60),
            cell_from_parameters(18.126, 18.126, 7.567, 90, 90, 120),
            cell_from_parameters(29.701, 29.701, 6.9204, 90, 90, 120),
        ]
        rng = np.random.default_rng(20260831)
        for cell in cells:
            for vec in rng.normal(size=(40, 3)) @ cell:
                actual = minimum_image(vec, cell)
                expected = brute_force(vec, cell)
                np.testing.assert_allclose(
                    np.linalg.norm(actual), np.linalg.norm(expected), atol=1e-10
                )

    def test_skew_cell_uses_exact_fallback(self):
        cell = np.array([
            [-0.08649378, -0.72901298, 0.26733954],
            [0.13792046, -1.46350011, -0.22432642],
            [-1.40883874, -2.21768125, -2.56051537],
        ])
        vec = np.array([3.49286597, -1.1532251, 7.17281376])
        actual = minimum_image(vec, cell)
        expected = brute_force(vec, cell)
        np.testing.assert_allclose(actual, expected, atol=1e-8)

    def test_invalid_cell_is_rejected(self):
        with self.assertRaises(ValueError):
            minimum_image([1.0, 2.0, 3.0], np.zeros((3, 3)))


if __name__ == "__main__":
    unittest.main()
