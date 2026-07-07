# Watershed segmentation for the phased charge density

## Installation
Clone the project to a local folder:
```
git clone https://github.com/pavel-kalouguine/PhasedChargeWatershed.git
```

Launch julia in the project folder and add [PKPackages](https://github.com/pavel-kalouguine/PKPackages) registry in the package manager (press `]` in Julia REPL to switch to Pkg REPL):
```
pkg> registry add https://github.com/pavel-kalouguine/PKPackages
```

Activate the local environment and instantiate the project:
```
(@v1.11) pkg> activate .
(PhasedChargeWatershed) pkg> instantiate
```