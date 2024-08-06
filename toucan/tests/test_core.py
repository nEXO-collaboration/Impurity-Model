import unittest
import json
import tempfile
import os
from toucan.setup import OutgassingSetup
from toucan.calculations import get_diff_temp, get_initial_impurities
from toucan.time_evolution import TimeEvolution
from toucan.utils import get_time_stamps


class TestOutgassingSetup(unittest.TestCase):
    def setUp(self):
        # Create a temporary JSON file with test data
        self.temp_dir = tempfile.mkdtemp()
        self.json_file = os.path.join(self.temp_dir, "test_library.json")

        test_data = {
            "Material": {
                "Teflon": {
                    "Oxygen": {
                        "Diffusion Constant": 31.4e-8,
                        "Solubility": 0.22,
                        "Activation Energy": 0.17,
                    }
                }
            },
            "System": {
                "EXO-200": {
                    "Xenon Mass": 200000,
                    "Teflon": {
                        "EXO-Teflon": {
                            "Volume": 0.693,
                            "Area": 9200.0,
                            "Thickness": 0.15,
                        }
                    },
                }
            },
            "Gas": {"Oxygen": {"Abundance in Air": 0.21, "Molar Mass": 32}},
        }

        with open(self.json_file, "w") as f:
            json.dump(test_data, f)

    def test_outgassing_setup_initialization(self):
        setup = OutgassingSetup(
            "EXO-200", "Teflon", "Oxygen", "EXO-Teflon", data_file=self.json_file
        )

        self.assertEqual(setup.name, "EXO-200")
        self.assertEqual(setup.material, "Teflon")
        self.assertEqual(setup.solute, "Oxygen")
        self.assertEqual(setup.diffusion, 31.4e-8)
        self.assertEqual(setup.solubility, 0.22)
        self.assertEqual(setup.activation_energy, 0.17)
        self.assertEqual(setup.xe_mass, 200000)
        self.assertEqual(setup.volume, 0.693)
        self.assertEqual(setup.area, 9200.0)
        self.assertEqual(setup.thickness, 0.15)
        self.assertEqual(setup.abundance, 0.21)
        self.assertEqual(setup.molar_mass, 32)

    def test_get_diff_temp(self):
        setup = OutgassingSetup(
            "EXO-200", "Teflon", "Oxygen", "EXO-Teflon", data_file=self.json_file
        )
        setup.temperatures = [293, 300]
        get_diff_temp(setup)
        self.assertEqual(len(setup.diffusion_constants), 2)
        self.assertAlmostEqual(setup.diffusion_constants[0], 31.4e-8, places=10)
        self.assertGreater(setup.diffusion_constants[1], setup.diffusion_constants[0])

    def test_get_initial_impurities(self):
        setup = OutgassingSetup(
            "EXO-200", "Teflon", "Oxygen", "EXO-Teflon", data_file=self.json_file
        )
        get_initial_impurities(setup, units="ppm")
        self.assertIsNotNone(setup.initial_impurities)
        self.assertGreater(setup.initial_impurities, 0)

    def test_get_impurities_vs_time(self):
        setup = OutgassingSetup(
            "EXO-200", "Teflon", "Oxygen", "EXO-Teflon", data_file=self.json_file
        )
        setup.temperatures = [293]
        get_diff_temp(setup)
        get_initial_impurities(setup, units="#")
        time_stamps = get_time_stamps(points=[0, 100], spacing=10, time_scale="Seconds")
        impurities = TimeEvolution.get_impurities_vs_time(setup, time=time_stamps)
        self.assertEqual(len(impurities), 1)  # One temperature
        self.assertEqual(len(impurities[0]), 11)  # 11 time points
        self.assertLess(
            impurities[0][-1], impurities[0][0]
        )  # Impurities should decrease over time

    def test_get_flow_rate_vs_time(self):
        setup = OutgassingSetup(
            "EXO-200", "Teflon", "Oxygen", "EXO-Teflon", data_file=self.json_file
        )
        setup.temperatures = [293]
        get_diff_temp(setup)
        get_initial_impurities(setup, units="#")
        time_stamps = get_time_stamps(points=[0, 100], spacing=10, time_scale="Seconds")
        impurities = TimeEvolution.get_impurities_vs_time(setup, time=time_stamps)
        flow_rates = TimeEvolution.get_flow_rate_vs_time(
            setup, time=time_stamps, impurities=impurities
        )
        self.assertEqual(len(flow_rates), 1)  # One temperature
        self.assertEqual(len(flow_rates[0]), 11)  # 11 time points
        self.assertLess(
            flow_rates[0][-1], flow_rates[0][0]
        )  # Flow rate should decrease over time

    def test_get_electron_lifetime_vs_time(self):
        setup = OutgassingSetup(
            "EXO-200", "Teflon", "Oxygen", "EXO-Teflon", data_file=self.json_file
        )
        time_stamps = get_time_stamps(points=[0, 100], spacing=10, time_scale="Seconds")
        lifetimes, params = TimeEvolution.get_electron_lifetime_vs_time(
            setup,
            time=time_stamps,
            initial_impurities=1,
            circulation_rate=1,
            purification_efficiency=1,
            out_diffusion=0.1,
        )
        self.assertEqual(len(lifetimes), 11)  # 11 time points
        self.assertGreater(
            lifetimes[-1], lifetimes[0]
        )  # Lifetime should increase over time
        self.assertIsInstance(params, dict)

    def tearDown(self):
        # Remove the temporary directory and its contents
        for file in os.listdir(self.temp_dir):
            os.remove(os.path.join(self.temp_dir, file))
        os.rmdir(self.temp_dir)


if __name__ == "__main__":
    unittest.main()
