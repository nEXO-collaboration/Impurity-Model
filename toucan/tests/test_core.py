import unittest
import json
import tempfile
import os
from toucan.core import Outgassing_setup

class TestOutgassingSetup(unittest.TestCase):
    def setUp(self):
        # Create a temporary JSON file with test data
        self.temp_dir = tempfile.mkdtemp()
        self.json_file = os.path.join(self.temp_dir, 'test_library.json')
        
        test_data = {
            "Material": {
                "Teflon": {
                    "Oxygen": {
                        "Diffusion Constant": 31.4e-8,
                        "Solubility": 0.22,
                        "Activation Energy": 0.17
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
                            "Thickness": 0.15
                        }
                    }
                }
            },
            "Gas": {
                "Oxygen": {
                    "Abundance in Air": 0.21,
                    "Molar Mass": 32
                }
            }
        }
        
        with open(self.json_file, 'w') as f:
            json.dump(test_data, f)

    def test_outgassing_setup_initialization(self):
        setup = Outgassing_setup("EXO-200", "Teflon", "Oxygen", data_file=self.json_file)
        
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

    def tearDown(self):
        # Remove the temporary directory and its contents
        for file in os.listdir(self.temp_dir):
            os.remove(os.path.join(self.temp_dir, file))
        os.rmdir(self.temp_dir)

if __name__ == '__main__':
    unittest.main()