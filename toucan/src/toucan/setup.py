import json
from typing import List, Union, Optional
from dataclasses import dataclass, field


@dataclass
class OutgassingSetup:
    """
    Represents the setup for outgassing experiments.

    Attributes:
        name: Name of the setup.
        material: Material used in the setup.
        solute: Solute used in the experiment.
        version: Version of the setup.
        data_file: Path to the JSON file containing setup data.
        temperatures: List of temperatures for the experiment.
        diffusion_constants: List of diffusion constants.
        diffusion: Diffusion constant of the material.
        solubility: Solubility of the solute in the material.
        activation_energy: Activation energy for the diffusion process.
        volume: Volume of the system.
        area: Surface area of the material.
        thickness: Thickness of the material.
        abundance: Abundance of the solute in air.
        molar_mass: Molar mass of the solute.
        xe_mass: Mass of xenon in the system.
        field_factor: Field factor for the system.
        comment: Additional comments about the setup.
    """

    name: str
    material: Optional[str] = None
    solute: Optional[str] = None
    version: Optional[str] = None
    data_file: str = "library.json"
    temperatures: List[Union[int, float]] = field(default_factory=list)
    diffusion_constants: List[float] = field(default_factory=list)
    diffusion: Optional[float] = None
    solubility: Optional[float] = None
    activation_energy: Optional[float] = None
    volume: Optional[float] = None
    area: Optional[float] = None
    thickness: Optional[float] = None
    abundance: Optional[float] = None
    molar_mass: Optional[float] = None
    xe_mass: Optional[float] = None
    field_factor: Optional[float] = None
    comment: Optional[str] = None

    def __post_init__(self):
        """
        Initializes the OutgassingSetup object by loading data from the JSON file.
        """
        self._load_data()

    def _load_data(self):
        """
        Loads data from the JSON file and populates the object attributes.
        """
        try:
            with open(self.data_file, "r") as file:
                data = json.load(file)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"The data file '{self.data_file}' was not found. "
                "Please ensure it's in the current working directory or provide the correct path."
            )

        self._load_material_properties(data)
        self._load_system_properties(data)
        self._load_gas_properties(data)
        self._load_xenon_properties(data)

    def _load_material_properties(self, data: dict):
        """
        Loads material-specific properties from the data dictionary.
        """
        if self.material and self.solute:
            material_props = (
                data.get("Material", {}).get(self.material, {}).get(self.solute, {})
            )
            self.diffusion = material_props.get("Diffusion Constant")
            self.solubility = material_props.get("Solubility")
            self.activation_energy = material_props.get("Activation Energy")

    def _load_system_properties(self, data: dict):
        """
        Loads system-specific properties from the data dictionary.
        """
        if self.name and self.material and self.version:
            system_props = (
                data.get("System", {})
                .get(self.name, {})
                .get(self.material, {})
                .get(self.version, {})
            )
            self.volume = system_props.get("Volume")
            self.area = system_props.get("Area")
            self.thickness = system_props.get("Thickness")

    def _load_gas_properties(self, data: dict):
        """
        Loads gas-specific properties from the data dictionary.
        """
        if self.solute:
            gas_props = data.get("Gas", {}).get(self.solute, {})
            self.abundance = gas_props.get("Abundance in Air")
            self.molar_mass = gas_props.get("Molar Mass")

    def _load_xenon_properties(self, data: dict):
        """
        Loads xenon-specific properties from the data dictionary.
        """
        system_props = data.get("System", {}).get(self.name, {})
        self.xe_mass = system_props.get("Xenon Mass")
        self.field_factor = system_props.get("Field Factor")
        self.comment = system_props.get("Comment")

    def __str__(self) -> str:
        """
        Returns a string representation of the OutgassingSetup object.
        """
        return "\n".join(
            f"{key}: {value}"
            for key, value in vars(self).items()
            if value not in [None, [], ""] and not key.startswith("_")
        )

    def to_dict(self) -> dict:
        """
        Converts the OutgassingSetup object to a dictionary.

        Returns:
            dict: A dictionary representation of the OutgassingSetup object.
        """
        return {
            key: value
            for key, value in vars(self).items()
            if value not in [None, [], ""] and not key.startswith("_")
        }

    def save_to_json(self, filename: str):
        """
        Saves the OutgassingSetup object to a JSON file.

        Args:
            filename (str): The name of the file to save the data to.
        """
        with open(filename, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load_from_json(cls, filename: str) -> "OutgassingSetup":
        """
        Creates an OutgassingSetup object from a JSON file.

        Args:
            filename (str): The name of the file to load the data from.

        Returns:
            OutgassingSetup: An OutgassingSetup object initialized with the data from the JSON file.
        """
        with open(filename, "r") as f:
            data = json.load(f)
        return cls(**data)
