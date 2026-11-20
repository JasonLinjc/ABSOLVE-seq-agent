import unittest

import numpy as np

from power_analysis import estimate, panel_wald


class ConditionalPowerTest(unittest.TestCase):
    def test_three_donor_wald_statistics_match_brglm2(self):
        # Independently fitted with glm(..., method=brglmFit, type="AS_mean").
        design = np.array([[1, a == 1, a == 2, d == 1, d == 2]
                           for d in range(3) for a in range(3)], dtype=float)
        depths = np.array([90, 85, 79, 73, 99, 121, 55, 84, 93])
        outcomes = np.array([[0, 1, 2, 1, 3, 7, 0, 2, 2], [0] * 9,
                             [65, 85, 56, 72, 97, 103, 47, 71, 77]])
        expected = [[1.38798686704, 1.97722956898], [-.12888437242, -.18602602161],
                    [3.34468685888, -1.53876290691]]
        np.testing.assert_allclose(panel_wald(outcomes, depths, design), expected, atol=2e-7)

    def test_pair_power_responds_to_depth_and_handles_unavailable_effects(self):
        def sim(n, background=0.001):
            return estimate("t", "Dedud", np.array([n, n]), background, ("HiFi",), None, 1000, 13)
        deep = sim(100000)
        self.assertEqual(deep, sim(100000))
        sparse = sim(1)
        self.assertEqual((sparse[0]["power"], sparse[0]["mc_low95"]), (0.0, 0.0))
        certain = sim(10000000)
        self.assertEqual((certain[1]["power"], certain[1]["mc_high95"]), (1.0, 1.0))
        self.assertGreater(deep[1]["power"], sim(100)[1]["power"])
        self.assertGreater(deep[1]["power"], deep[0]["power"])
        self.assertTrue(all(r["power"] is None and r["simulations"] == 0 for r in sim(0)))
        self.assertEqual(sim(100, 0.9995)[1]["status"], "effect exceeds probability range")
