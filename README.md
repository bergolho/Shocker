# Shocker

Program to generate patient-specific Purkinje Networks within an endocardium surface given by triangles.

## Dependencies

- C++ compiler
- CMake
- VTK library

## How to install

To install CMake, run the following command in a Fedora machine:

```sh
$ sudo dnf install cmake
```

The VTK library can be installed alongside the Paraview program which will enable the visualization of the generated Purkinje networks. 

To install Paraview, run the following command in a Fedora machine:

```sh
$ sudo dnf install paraview
```

## Compile and run

To compile the source code you can use the available script.

```sh
$ ./recompile_project.sh
```

To run a simulation you need to pass a configuration file as input. There are some samples in the ```inputs``` folder.

Generate a Purkinje network in the LV surface.
```sh
./bin/Shocker inputs/finsberg_LV.ini
```

Generate a Purkinje network in the RV surface.
```sh
./bin/Shocker inputs/finsberg_RV.ini
```

The output files are stored in the ```outputs``` folder. Remember to create this folder before running a simulation. 

The program generates a Purkinje network that connects the active Purkinje-Ventricular-Junctions (PVJs) with their respective locations and Local Activation Times (LATs).

![Screenshot](figs/shocker_snapshot_1.png)

![Screenshot](figs/shocker_snapshot_2.png)
