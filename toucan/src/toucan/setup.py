import json
from typing import List, Union, Optional


class OutgassingSetup:
    def __init__(
        self,
        setup: str,
        material: Optional[str] = None,
        solute: Optional[str] = None,
        version: Optional[str] = None,
        data_file: str = "library.json",
    ):
        # Load data from the JSON file
        try:
            with open(data_file, "r") as file:
                data = json.load(file)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"The data file '{data_file}' was not found. Please ensure it's in the current working directory or provide the correct path."
            )

        self.name: str = setup
        self.material: Optional[str] = material
        self.solute: Optional[str] = solute
        self.version: Optional[str] = version

        self.temperatures: List[Union[int, float]] = []
        self.diffusion_constants: List[float] = []

        # Retrieve material properties safely from JSON data
        if material and solute:
            material_props = data["Material"].get(material, {}).get(solute, {})
            self.diffusion: Optional[float] = material_props.get("Diffusion Constant")
            self.solubility: Optional[float] = material_props.get("Solubility")
            self.activation_energy: Optional[float] = material_props.get(
                "Activation Energy"
            )

        # Retrieve system properties safely from JSON data
        if setup and material and version:
            system_props = (
                data["System"].get(setup, {}).get(material, {}).get(version, {})
            )
            self.volume: float = system_props.get("Volume")
            self.area: float = system_props.get("Area")
            self.thickness: float = system_props.get("Thickness")

        # Retrieve gas properties safely from JSON data
        if solute:
            gas_props = data["Gas"].get(solute, {})
            self.abundance: Optional[float] = gas_props.get("Abundance in Air")
            self.molar_mass: Optional[float] = gas_props.get("Molar Mass")

        # Retrieve Xenon Mass and Field Factor from JSON data
        self.xe_mass: Optional[float] = data["System"].get(setup, {}).get("Xenon Mass")
        self.field_factor: Optional[float] = (
            data["System"].get(setup, {}).get("Field Factor")
        )
        self.comment: Optional[str] = data["System"].get(setup, {}).get("Comment")

    def __str__(self) -> str:
        attributes = vars(self)
        non_empty_attributes = {
            item: attributes[item]
            for item in attributes
            if attributes[item] not in [None, [], ""]
        }
        return "\n".join(
            f"{key}: {value}" for key, value in non_empty_attributes.items()
        )
