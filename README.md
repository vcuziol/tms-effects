# TMS Effects model

Python code by V. V. Cuziol to accompany the Computational Model of the Effects of Transcranial Magnetic Stimulation on Cortical Networks *[paper](https://link.springer.com/chapter/10.1007/978-3-030-70601-2_338)* and *[master's thesis](https://www.teses.usp.br/teses/disponiveis/59/59143/tde-22062020-195016/publico/VitorCuziol_dissertacao_corrigida.pdf)*.

This code was published only to be read/consulted and for research purposes; and the instructions below are only meant as a guide and are not guaranteed to work.

## Guide (unfinished; not tested)
1) NEURON 7.5 is required. Compile the mechanisms of the Aberra et al. 2019 models by running the command 'nrnivmodl' (in Linux) at the root folder of the mechanisms (AberraEtAl2018_edit/mechanisms), which contains .mod files. To run the 'nrnivmodl' command, you must have NEURON installed.
2) One or more folders resulting from SimNIBS simulations must be placed inside the 'msh_inputs' folder.
3) Inside the 'motor_cortex_indices' folder, you must create two text files with polygon indices indicating the polygons of the motor cortex region in both white matter (WM) and gray matter (GM) meshes. It is recommended that you do this by importing the "obj" files of both GM and WM models into Blender, select the desired region of the motor cortex in each model by using the Ctrl+'+' command, and run the 'print_selected_faces' script in Blender's Python console. The number of simulated neurons is given by the number of polygons in the motor cortex region selected in the gray matter (GM) mesh. It is better that the WM motor cortex selection covers a wider area than the GM motor cortex, because this helps the calculations of the cortical thickness approximation and subsequent positioning of the neurons according to the thickness.
4) The first file to run is 'preprocessing.py', which reads both the SimNIBS .msh file (or files) and the GM/WM motor cortex indices and generates the corresponding meshes to be used in the simulations.
5) If you wish to run simulations for more than one coil orientation, then you must run 'generate_and_run_simnibs_files.py' with a .mat file (generated from SimNIBS) as a basis file, so that the other files with different coil orientations configurations can be generated.
